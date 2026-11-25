#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "Data_Type.h"
#include "HighDim_Alg.h"
#include "Alg_Solve.h"
#include "Alg_BandSolve.h"

// Verification of the sparse direct solver LinBandSys (include/Alg_BandSolve.h)
//
// Usage: ./TestBandSolve [-v]
//   -v  print every check, not only the failed ones
//
// Part 1: storage of the matrix (add, get, reset, norms)
// Part 2: solutions of matrices of the kinds met in TDSolveSacado and of
//         matrices that need row exchanges, against the solution by the
//         inverse of the full matrix (HMAT::inv, the former direct solver)
// Part 3: the renumbering: bandwidths that are found
// Part 4: reuse of a system with the same and with a new sparsity pattern
// Part 5: singular matrices
// Part 6: exchange with LinMat/LinVec and LinSolve
//
// The program returns 0 if all checks pass and 1 otherwise.

// --------------------------- Bookkeeping ---------------------------

static int  num_pass = 0;
static int  num_fail = 0;
static bool verbose  = false;

void check(const bool ok, const std::string& what) {
  if ( ok ) {
    num_pass++;
    if ( verbose ) std::cout << "  pass: " << what << std::endl;
  }
  else {
    num_fail++;
    std::cout << "  FAIL: " << what << std::endl;
  }
  return;
}

void check_small(const DType val, const DType tol, const std::string& what) {
  std::ostringstream os;
  os.precision(3);
  os << what << " (" << val << " <= " << tol << ")";
  check( val <= tol, os.str() );
  return;
}

void check_equal(const int a, const int b, const std::string& what) {
  std::ostringstream os;
  os << what << " (" << a << " vs " << b << ")";
  check( a == b, os.str() );
  return;
}

// Reproducible pseudo-random numbers in [0,1)
static unsigned long rand_state = 12345;
DType rand_unit(void) {
  rand_state = ( rand_state * 1103515245UL + 12345UL ) % 2147483648UL;
  return rand_state / 2147483648.0;
}

// ----------------------- A matrix to be tested ----------------------

// A list of entries (i,j,val) and a right hand side
struct TestMat {
  int n;
  std::vector<int>   row, col;
  std::vector<DType> val;
  std::vector<DType> rhs;

  explicit TestMat(const int n_in) : n(n_in), rhs(n_in, 0.0) {
    for ( int i = 0; i < n; i++ ) rhs[i] = 2.0*rand_unit() - 1.0;
  }
  void add(const int i, const int j, const DType v) {
    row.push_back( i ); col.push_back( j ); val.push_back( v );
  }
  void fill(LinBandSys& sys) const {
    for ( std::size_t k = 0; k < val.size(); k++ ) sys.add( row[k], col[k], val[k] );
    for ( int i = 0; i < n; i++ ) sys.rhs(i) = rhs[i];
  }
};

// Solve by LinBandSys and compare with the solution by the inverse of the
// full matrix. Returns the system through sys_ot for further checks.
void check_solve(const std::string& name, const TestMat& mat, LinBandSys& sys, const DType tol = 1.0e-10) {
  const int n = mat.n;
  sys.reset();
  mat.fill( sys );

  const bool regular = sys.solve();
  check( regular, name + ": solved" );

  // Reference solution
  HMAT mat_ref( n, 0.0 );
  HVEC rhs_ref( n, 0.0 );
  for ( std::size_t k = 0; k < mat.val.size(); k++ ) mat_ref[mat.row[k]][mat.col[k]] += mat.val[k];
  for ( int i = 0; i < n; i++ ) rhs_ref[i] = mat.rhs[i];
  HVEC sol_ref = mat_ref.inv() * rhs_ref;

  DType diff = 0.0, scal = 0.0, rhs_nrm = 0.0;
  for ( int i = 0; i < n; i++ ) {
    diff = std::max( diff, std::fabs( sys.sol(i) - sol_ref[i] ) );
    scal = std::max( scal, std::fabs( sol_ref[i] ) );
    rhs_nrm += mat.rhs[i] * mat.rhs[i];
  }
  rhs_nrm = std::sqrt( rhs_nrm );

  check_small( diff/std::max( scal, 1.0 ), tol, name + ": difference to the solution by the inverse" );
  check_small( sys.getResNorm()/std::max( rhs_nrm, 1.0 ), tol, name + ": residual" );
  check_small( std::fabs( sys.getRhsNorm() - rhs_nrm ), 1.0e-14*std::max( rhs_nrm, 1.0 ), name + ": norm of the right hand side" );
  return;
}

// Tridiagonal matrix as for one chemical: mass matrix + diffusion
TestMat mat_tridiag(const int n) {
  TestMat mat(n);
  for ( int i = 0; i < n; i++ ) {
    mat.add( i, i, 1.0 );
    for ( int j = std::max(0,i-1); j <= std::min(n-1,i+1); j++ )
      mat.add( i, j, ( i == j ) ? 0.8 + 0.4*rand_unit() : -0.4*rand_unit() );
  }
  return mat;
}

// Several species, each with num unknowns, numbered species by species
// (unknown i of species s is s*num+i), as in TDSolveSacado:
//   sten >= 0: species are coupled to each other in the same interval and
//              to themselves in the intervals i-sten to i+sten
//   sten <  0: species are coupled to all species in all intervals (full)
TestMat mat_species(const int num_spec, const int num, const int sten) {
  TestMat mat(num_spec*num);
  for ( int s = 0; s < num_spec; s++ ) {
    for ( int i = 0; i < num; i++ ) {
      mat.add( s*num+i, s*num+i, 1.0 );
      for ( int t = 0; t < num_spec; t++ )
        mat.add( s*num+i, t*num+i, 0.2*(2.0*rand_unit()-1.0) );
      if ( sten >= 0 ) {
        for ( int j = std::max(0,i-sten); j <= std::min(num-1,i+sten); j++ )
          mat.add( s*num+i, s*num+j, 0.2*(2.0*rand_unit()-1.0) );
      }
      else {
        for ( int t = 0; t < num_spec; t++ )
          for ( int j = 0; j < num; j++ )
            mat.add( s*num+i, t*num+j, 0.5*(2.0*rand_unit()-1.0)/(num_spec*num) );
      }
    }
  }
  return mat;
}

// A banded matrix whose unknowns are numbered at random. The diagonal is
// not dominant, so that rows have to be exchanged.
TestMat mat_shuffled_band(const int n, const int kl, const int ku, std::vector<int>* perm_ot = NULL) {
  std::vector<int> perm(n);
  for ( int i = 0; i < n; i++ ) perm[i] = i;
  for ( int i = n-1; i > 0; i-- ) std::swap( perm[i], perm[(int)(rand_unit()*(i+1))] );
  TestMat mat(n);
  for ( int i = 0; i < n; i++ )
    for ( int j = std::max(0,i-kl); j <= std::min(n-1,i+ku); j++ )
      mat.add( perm[i], perm[j], 2.0*rand_unit() - 1.0 + ( i == j ? 0.5 : 0.0 ) );
  if ( perm_ot != NULL ) *perm_ot = perm;
  return mat;
}

// ------------------------------ Part 1 ------------------------------

void test_storage(void) {
  std::cout << "Part 1: storage" << std::endl;
  LinBandSys sys(4);
  check_equal( sys.getSize(), 4, "size" );
  check_equal( sys.getNumEntries(), 0, "no entry at the beginning" );
  check( sys.get(1,2) == 0.0, "an entry that is not registered is zero" );

  sys.add( 1, 2, 3.0 );
  sys.add( 1, 0, -1.0 );
  sys.add( 1, 2, 0.5 );
  sys.add( 3, 3, 2.0 );
  check_equal( sys.getNumEntries(), 3, "three entries registered" );
  check( sys.get(1,2) == 3.5, "add accumulates" );
  check( sys.get(1,0) == -1.0 && sys.get(3,3) == 2.0, "entries added in any order" );
  check( sys.get(2,1) == 0.0 && sys.get(1,1) == 0.0, "other entries are zero" );

  sys.rhs(0) = 3.0;
  sys.rhs(2) = -4.0;
  check( sys.getRhsNorm() == 5.0, "norm of the right hand side" );
  check( sys.rhs(0) == 3.0 && sys.rhs(1) == 0.0, "right hand side" );

  sys.reset();
  check_equal( sys.getNumEntries(), 3, "reset keeps the entries registered" );
  check( sys.get(1,2) == 0.0 && sys.get(3,3) == 0.0, "reset sets the entries to zero" );
  check( sys.getRhsNorm() == 0.0 && sys.getSolInfNorm() == 0.0, "reset sets right hand side and solution to zero" );
  return;
}

// ------------------------------ Part 2 ------------------------------

void test_solve(void) {
  std::cout << "Part 2: solutions" << std::endl;

  // A single unknown, and a 2 x 2 matrix with zero diagonal
  {
    TestMat mat(1);
    mat.add( 0, 0, -4.0 );
    LinBandSys sys(1);
    check_solve( "1 x 1", mat, sys );
  }
  {
    TestMat mat(2);
    mat.add( 0, 1, 2.0 );
    mat.add( 1, 0, -3.0 );
    LinBandSys sys(2);
    check_solve( "2 x 2 with zero diagonal", mat, sys );
  }

  // One chemical
  const int sizes[5] = { 2, 3, 10, 57, 400 };
  for ( int k = 0; k < 5; k++ ) {
    TestMat mat = mat_tridiag( sizes[k] );
    LinBandSys sys( sizes[k] );
    check_solve( "tridiagonal, n = " + std::to_string(sizes[k]), mat, sys );
    check( sys.isBanded() || sizes[k] < 10, "tridiagonal, n = " + std::to_string(sizes[k]) + ": banded solver" );
  }

  // Cell species
  for ( int num_spec = 1; num_spec <= 5; num_spec += 2 ) {
    for ( int sten = 0; sten <= 3; sten += ( sten == 0 ? 1 : 2 ) ) {
      TestMat mat = mat_species( num_spec, 40, sten );
      LinBandSys sys( mat.n );
      check_solve( std::to_string(num_spec) + " species x 40 intervals, stencil " + std::to_string(sten), mat, sys );
    }
  }
  {
    TestMat mat = mat_species( 3, 30, -1 );
    LinBandSys sys( mat.n );
    check_solve( "3 species x 30 intervals, full", mat, sys );
    check( !sys.isBanded(), "3 species x 30 intervals, full: full solver" );
  }

  // Row exchanges
  for ( int k = 0; k < 6; k++ ) {
    const int n  = 20 + 37*k;
    const int kl = 1 + k % 3;
    const int ku = 1 + (k+1) % 4;
    TestMat mat = mat_shuffled_band( n, kl, ku );
    LinBandSys sys( n );
    check_solve( "band (" + std::to_string(kl) + "," + std::to_string(ku) + ") numbered at random, n = " + std::to_string(n), mat, sys );
  }
  {
    // Zero diagonal: every elimination step needs an exchange
    const int n = 40;
    TestMat mat(n);
    for ( int i = 0; i < n; i++ ) {
      if ( i > 0 )   mat.add( i, i-1, 1.0 + 0.1*rand_unit() );
      if ( i < n-1 ) mat.add( i, i+1, 1.0 + 0.1*rand_unit() );
      mat.add( i, i, 0.0 );
    }
    LinBandSys sys(n);
    check_solve( "tridiagonal with zero diagonal", mat, sys );
  }

  // Independent blocks, and a matrix with one full row and column
  {
    const int n = 30;
    TestMat mat(n);
    for ( int i = 0; i < n; i++ ) mat.add( i, i, 1.0 + rand_unit() );
    for ( int i = 0; i+1 < n; i += 3 ) { mat.add( i, i+1, 0.3 ); mat.add( i+1, i, -0.2 ); }
    LinBandSys sys(n);
    check_solve( "independent blocks", mat, sys );
    check_equal( std::max( sys.getLowerBand(), sys.getUpperBand() ), 1, "independent blocks: bandwidth" );
  }
  {
    const int n = 50;
    TestMat mat(n);
    for ( int i = 0; i < n; i++ ) {
      mat.add( i, i, 2.0 + rand_unit() );
      mat.add( 0, i, 0.1*rand_unit() );
      mat.add( i, 0, 0.1*rand_unit() );
    }
    LinBandSys sys(n);
    check_solve( "one full row and column", mat, sys );
  }

  // Sparse at random
  for ( int k = 0; k < 4; k++ ) {
    const int n = 25 + 40*k;
    TestMat mat(n);
    for ( int i = 0; i < n; i++ ) {
      mat.add( i, i, 3.0 + rand_unit() );
      for ( int m = 0; m < 3; m++ ) mat.add( i, (int)(rand_unit()*n), 2.0*rand_unit() - 1.0 );
    }
    LinBandSys sys(n);
    check_solve( "sparse at random, n = " + std::to_string(n), mat, sys );
  }
  return;
}

// ------------------------------ Part 3 ------------------------------

void test_numbering(void) {
  std::cout << "Part 3: renumbering" << std::endl;

  {
    TestMat mat = mat_tridiag( 200 );
    LinBandSys sys( 200 );
    check_solve( "tridiagonal", mat, sys );
    check_equal( sys.getLowerBand(), 1, "tridiagonal: sub-diagonals" );
    check_equal( sys.getUpperBand(), 1, "tridiagonal: super-diagonals" );
  }

  // The bandwidth of coupled species does not depend on the number of intervals
  for ( int num = 25; num <= 400; num *= 4 ) {
    TestMat mat = mat_species( 3, num, 1 );
    LinBandSys sys( mat.n );
    check_solve( "3 species x " + std::to_string(num) + " intervals", mat, sys );
    check( sys.isBanded(), "3 species x " + std::to_string(num) + " intervals: banded solver" );
    check( std::max( sys.getLowerBand(), sys.getUpperBand() ) <= 5,
        "3 species x " + std::to_string(num) + " intervals: at most 5 diagonals on each side (" +
        std::to_string( sys.getLowerBand() ) + "," + std::to_string( sys.getUpperBand() ) + ")" );
  }
  {
    TestMat mat = mat_species( 10, 60, 1 );
    LinBandSys sys( mat.n );
    check_solve( "10 species x 60 intervals", mat, sys );
    check( sys.isBanded() && std::max( sys.getLowerBand(), sys.getUpperBand() ) <= 20,
        "10 species x 60 intervals: at most 20 diagonals on each side (" +
        std::to_string( sys.getLowerBand() ) + "," + std::to_string( sys.getUpperBand() ) + ")" );
  }

  // A band that is hidden by the numbering is found again
  {
    TestMat mat = mat_shuffled_band( 300, 2, 2 );
    LinBandSys sys( 300 );
    check_solve( "band (2,2) numbered at random", mat, sys );
    check( sys.isBanded() && std::max( sys.getLowerBand(), sys.getUpperBand() ) <= 4,
        "band (2,2) numbered at random: at most 4 diagonals on each side (" +
        std::to_string( sys.getLowerBand() ) + "," + std::to_string( sys.getUpperBand() ) + ")" );
  }
  return;
}

// ------------------------------ Part 4 ------------------------------

void test_reuse(void) {
  std::cout << "Part 4: reuse" << std::endl;
  const int n = 60;
  LinBandSys sys(n);

  // Same sparsity, new values, several times
  for ( int k = 0; k < 3; k++ ) {
    TestMat mat = mat_tridiag( n );
    check_solve( "tridiagonal, system " + std::to_string(k+1), mat, sys );
  }
  check_equal( sys.getLowerBand(), 1, "still tridiagonal" );

  // New entries far from the diagonal
  {
    TestMat mat = mat_tridiag( n );
    mat.add( 3, n-2, 0.7 );
    mat.add( n-5, 1, -0.4 );
    check_solve( "tridiagonal with two more entries", mat, sys );
  }
  // The former entries are still registered, with value zero
  {
    TestMat mat = mat_tridiag( n );
    check_solve( "tridiagonal again", mat, sys );
    check( sys.get(3,n-2) == 0.0, "entry that is not set any more is zero" );
  }
  // Now a full matrix in the same system
  {
    TestMat mat(n);
    for ( int i = 0; i < n; i++ )
      for ( int j = 0; j < n; j++ )
        mat.add( i, j, ( i == j ? 4.0 : 0.0 ) + 2.0*rand_unit() - 1.0 );
    check_solve( "full matrix in the same system", mat, sys );
    check( !sys.isBanded(), "full matrix in the same system: full solver" );
  }
  return;
}

// ------------------------------ Part 5 ------------------------------

void test_singular(void) {
  std::cout << "Part 5: singular matrices (one warning is expected for each of the three)" << std::endl;
  {
    // A zero row, banded
    const int n = 20;
    TestMat mat = mat_tridiag( n );
    for ( std::size_t k = 0; k < mat.val.size(); k++ ) if ( mat.row[k] == 7 ) mat.val[k] = 0.0;
    LinBandSys sys(n);
    mat.fill( sys );
    check( !sys.solve(), "zero row in a banded matrix is detected" );
    check( sys.getSolInfNorm() == 0.0, "zero row in a banded matrix: zero solution" );
    // The system can still be used
    TestMat mat_ok = mat_tridiag( n );
    check_solve( "regular matrix afterwards", mat_ok, sys );
  }
  {
    // Two equal rows, full
    const int n = 12;
    TestMat mat(n);
    std::vector<DType> row0(n);
    for ( int j = 0; j < n; j++ ) row0[j] = 1.0 + j;
    for ( int i = 0; i < n; i++ )
      for ( int j = 0; j < n; j++ )
        mat.add( i, j, ( i < 2 ) ? row0[j] : ( i == j ? 3.0 : 0.0 ) + rand_unit() );
    LinBandSys sys(n);
    mat.fill( sys );
    check( !sys.solve(), "two equal rows in a full matrix are detected" );
    check( sys.getSolInfNorm() == 0.0, "two equal rows in a full matrix: zero solution" );
  }
  {
    // No entry at all
    LinBandSys sys(5);
    sys.rhs(2) = 1.0;
    check( !sys.solve(), "empty matrix is detected" );
  }
  return;
}

// ------------------------------ Part 6 ------------------------------

void test_exchange(void) {
  std::cout << "Part 6: exchange with LinMat, LinVec and LinSolve" << std::endl;
  const int n = 30;
  TestMat mat = mat_species( 3, 10, 1 );
  LinBandSys sys(n);
  mat.fill( sys );

  LinMat* lmat = new LinSqSparseMat;
  lmat->setBlckDim(1);
  lmat->setEqualSize(n);
  LinMat* lprc = new LinSqSparseMat;
  lprc->setBlckDim(1);
  lprc->setEqualSize(n);
  LinVec* lsol = new LinVec(1,n);
  LinVec* lrhs = new LinVec(1,n);
  LinVec* lres = new LinVec(1,n);

  sys.copyTo( *lmat, *lrhs );
  bool same = true;
  for ( int i = 0; i < n; i++ ) {
    for ( int j = 0; j < n; j++ )
      if ( lmat->getGlobVal(i,j) != sys.get(i,j) ) same = false;
    if ( lrhs->getGlobVal(i) != sys.rhs(i) ) same = false;
  }
  check( same, "copyTo: same matrix and right hand side" );

  // The direct solver of LinSolve (inverse of the full matrix)
  LinSolve linSolve;
  linSolve.setBlckDim(1);
  linSolve.setSize(n);
  linSolve.linkMat(lmat);
  linSolve.linkPrc(lprc);
  linSolve.linkSol(lsol);
  linSolve.linkRhs(lrhs);
  linSolve.linkRes(lres);
  linSolve.setPrintStatus(false);
  linSolve.setSolvePrec(LSYS_PREC_OFF);
  linSolve.setSolveDir();
  linSolve.solve();

  LinBandSys sys_copy(n);
  mat.fill( sys_copy );
  sys_copy.copySolFrom( *lsol );
  sys.solve();
  DType diff = 0.0;
  for ( int i = 0; i < n; i++ ) diff = std::max( diff, std::fabs( sys.sol(i) - sys_copy.sol(i) ) );
  check_small( diff, 1.0e-12, "copySolFrom: solution of LinSolve agrees with the banded solver" );
  check_small( sys_copy.getResNorm(), 1.0e-12, "copySolFrom: residual of the solution of LinSolve" );

  delete lmat;
  delete lprc;
  delete lsol;
  delete lrhs;
  delete lres;
  return;
}

int main(int argc, char** argv) {
  if ( argc > 1 && std::strcmp( argv[1], "-v" ) == 0 ) verbose = true;

  test_storage();
  test_solve();
  test_numbering();
  test_reuse();
  test_singular();
  test_exchange();

  std::cout << std::endl << num_pass << " checks passed, " << num_fail << " failed." << std::endl;
  return ( num_fail == 0 ) ? 0 : 1;
}
