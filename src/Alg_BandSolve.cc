#include <cmath>
#include <cassert>
#include <cstdio>
#include "Alg_BandSolve.h"
#include "Alg_Solve.h"

// LAPACK: solution of a general (full) linear system by LU factorization
// with partial pivoting
extern "C" void dgesv_(const int* N, const int* NRHS, double* A, const int* LDA, int* IPIV, double* B, const int* LDB, int* INFO);
extern "C" void sgesv_(const int* N, const int* NRHS, float*  A, const int* LDA, int* IPIV, float*  B, const int* LDB, int* INFO);

namespace {

  inline void lapack_gesv(const int* N, const int* NRHS, double* A, const int* LDA, int* IPIV, double* B, const int* LDB, int* INFO) {
    dgesv_( N, NRHS, A, LDA, IPIV, B, LDB, INFO );
    return;
  }

  inline void lapack_gesv(const int* N, const int* NRHS, float* A, const int* LDA, int* IPIV, float* B, const int* LDB, int* INFO) {
    sgesv_( N, NRHS, A, LDA, IPIV, B, LDB, INFO );
    return;
  }

  // Breadth-first search from root through the unknowns that are not yet
  // numbered (numbered[i] == false). The unknowns are appended to queue in
  // the order they are reached; the neighbours of an unknown are visited in
  // the order given by nbr, i.e. by increasing degree.
  // Returns the number of levels; last_beg is the position in queue of the
  // first unknown of the last level.
  int search_levels( const int root, const std::vector< std::vector<int> >& nbr,
      const std::vector<bool>& numbered, std::vector<int>& mark, const int stamp,
      std::vector<int>& queue, int& last_beg ) {
    queue.clear();
    queue.push_back( root );
    mark[root] = stamp;
    int num_level = 0;
    int level_beg = 0;
    last_beg = 0;
    while ( level_beg < (int)queue.size() ) {
      const int level_end = (int)queue.size();
      last_beg = level_beg;
      num_level++;
      for ( int q = level_beg; q < level_end; q++ ) {
        const std::vector<int>& loc_nbr = nbr[queue[q]];
        for ( std::size_t k = 0; k < loc_nbr.size(); k++ ) {
          const int j = loc_nbr[k];
          if ( numbered[j] || mark[j] == stamp ) continue;
          mark[j] = stamp;
          queue.push_back( j );
        }
      }
      level_beg = level_end;
    }
    return num_level;
  }

}

LinBandSys::LinBandSys(const int size_in) {
  assert( size_in > 0 );
  size = size_in;
  row_cols.resize( size );
  row_vals.resize( size );
  vec_rhs.assign( size, 0.0 );
  vec_sol.assign( size, 0.0 );
  new_idx.resize( size );
  for ( int i = 0; i < size; i++ ) new_idx[i] = i;
  pattern_changed = true;
  use_band = false;
  kl = 0;
  ku = 0;
  warned_singular = false;
}

void LinBandSys::reset(void) {
  for ( int i = 0; i < size; i++ ) {
    std::vector<DType>& vals = row_vals[i];
    for ( std::size_t k = 0; k < vals.size(); k++ ) vals[k] = 0.0;
    vec_rhs[i] = 0.0;
    vec_sol[i] = 0.0;
  }
  return;
}

DType LinBandSys::get(const int i, const int j) const {
  const std::vector<int>& cols = row_cols[i];
  const std::vector<int>::const_iterator it = std::lower_bound( cols.begin(), cols.end(), j );
  if ( it == cols.end() || *it != j ) return 0.0;
  return row_vals[i][it-cols.begin()];
}

DType LinBandSys::getRhsNorm(void) const {
  DType val = 0.0;
  for ( int i = 0; i < size; i++ ) val += vec_rhs[i] * vec_rhs[i];
  return std::sqrt( val );
}

DType LinBandSys::getSolInfNorm(void) const {
  DType val = 0.0;
  for ( int i = 0; i < size; i++ ) val = std::max( val, (DType)std::fabs( vec_sol[i] ) );
  return val;
}

DType LinBandSys::getResNorm(void) const {
  DType val = 0.0;
  for ( int i = 0; i < size; i++ ) {
    const std::vector<int>&   cols = row_cols[i];
    const std::vector<DType>& vals = row_vals[i];
    DType res = 0.0;
    for ( std::size_t k = 0; k < cols.size(); k++ ) res += vals[k] * vec_sol[cols[k]];
    res = vec_rhs[i] - res;
    val += res * res;
  }
  return std::sqrt( val );
}

int LinBandSys::getNumEntries(void) const {
  int num = 0;
  for ( int i = 0; i < size; i++ ) num += (int)row_cols[i].size();
  return num;
}

bool LinBandSys::solve(void) {
  if ( pattern_changed ) setup_numbering();

  bool regular = use_band ? solve_band() : solve_full();

  if ( !regular ) {
    // HMAT::inv returns a zero matrix in this case, i.e. also a zero solution
    if ( !warned_singular ) {
      fprintf( stderr, "Warning: LinBandSys::solve - the matrix is singular, the solution is set to zero.\n" );
      warned_singular = true;
    }
    for ( int i = 0; i < size; i++ ) vec_sol[i] = 0.0;
  }
  return regular;
}

// Renumber the unknowns by the reverse Cuthill-McKee algorithm to reduce the
// bandwidth, and decide whether the banded or the full solver is used
void LinBandSys::setup_numbering(void) {
  const int n = size;
  pattern_changed = false;

  for ( int i = 0; i < n; i++ ) new_idx[i] = i;

  // A matrix that is (almost) full cannot be turned into a banded one
  if ( 2.0*getNumEntries() > 1.0*n*n ) {
    calc_bandwidth();
    use_band = false;
    return;
  }

  // Graph of the matrix: i and j are neighbours if mat[i][j] or mat[j][i]
  // is registered
  std::vector< std::vector<int> > nbr( n );
  for ( int i = 0; i < n; i++ ) {
    const std::vector<int>& cols = row_cols[i];
    for ( std::size_t k = 0; k < cols.size(); k++ ) {
      const int j = cols[k];
      if ( j == i ) continue;
      nbr[i].push_back( j );
      nbr[j].push_back( i );
    }
  }
  std::vector<int> degree( n );
  for ( int i = 0; i < n; i++ ) {
    std::sort( nbr[i].begin(), nbr[i].end() );
    nbr[i].erase( std::unique( nbr[i].begin(), nbr[i].end() ), nbr[i].end() );
    degree[i] = (int)nbr[i].size();
  }
  // Neighbours are visited by increasing degree (and index, to be definite)
  for ( int i = 0; i < n; i++ ) {
    std::sort( nbr[i].begin(), nbr[i].end(),
        [&degree](const int a, const int b) {
          return ( degree[a] != degree[b] ) ? ( degree[a] < degree[b] ) : ( a < b );
        } );
  }

  std::vector<bool> numbered( n, false );
  std::vector<int>  mark( n, 0 );
  std::vector<int>  order, queue;
  order.reserve( n );
  queue.reserve( n );
  int stamp = 0;

  // The unknowns that have not been numbered yet, by increasing degree
  std::vector<int> by_degree( n );
  for ( int i = 0; i < n; i++ ) by_degree[i] = i;
  std::sort( by_degree.begin(), by_degree.end(),
      [&degree](const int a, const int b) {
        return ( degree[a] != degree[b] ) ? ( degree[a] < degree[b] ) : ( a < b );
      } );

  // One connected part of the graph at a time
  for ( int s = 0; s < n; s++ ) {
    int root = by_degree[s];
    if ( numbered[root] ) continue;

    // Look for a starting point at the "end" of the graph: move to an
    // unknown of smallest degree in the last level as long as this gives
    // more levels (levels are then narrower, hence a smaller bandwidth)
    int last_beg;
    int num_level = search_levels( root, nbr, numbered, mark, ++stamp, queue, last_beg );
    for ( int trial = 0; trial < 8; trial++ ) {
      int cand = queue[last_beg];
      for ( std::size_t q = last_beg; q < queue.size(); q++ )
        if ( degree[queue[q]] < degree[cand] ) cand = queue[q];
      if ( cand == root ) break;
      int cand_last_beg;
      std::vector<int> cand_queue;
      const int cand_level = search_levels( cand, nbr, numbered, mark, ++stamp, cand_queue, cand_last_beg );
      if ( cand_level <= num_level ) break;
      root      = cand;
      num_level = cand_level;
      last_beg  = cand_last_beg;
      queue.swap( cand_queue );
    }

    // Cuthill-McKee numbering of this part: the order of the search
    for ( std::size_t q = 0; q < queue.size(); q++ ) {
      numbered[queue[q]] = true;
      order.push_back( queue[q] );
    }
  }
  assert( (int)order.size() == n );

  // Reverse the order
  for ( int k = 0; k < n; k++ ) new_idx[order[k]] = n-1-k;

  calc_bandwidth();

  // Banded elimination costs about size*kl*(kl+ku) operations and a full
  // one size^3/3; the full one is done by LAPACK
  use_band = ( 4*(kl+ku) < n );

  if ( !use_band )
    for ( int i = 0; i < n; i++ ) new_idx[i] = i;

  return;
}

void LinBandSys::calc_bandwidth(void) {
  kl = 0;
  ku = 0;
  for ( int i = 0; i < size; i++ ) {
    const std::vector<int>& cols = row_cols[i];
    for ( std::size_t k = 0; k < cols.size(); k++ ) {
      const int d = new_idx[i] - new_idx[cols[k]];
      if (  d > kl ) kl =  d;
      if ( -d > ku ) ku = -d;
    }
  }
  return;
}

// Gaussian elimination with partial pivoting of the banded matrix.
// Row r of the renumbered matrix has its entries in the columns r-kl to
// r+ku. Exchanging rows fills in up to kl more columns to the right, so
// that the columns r-kl to r+ku+kl are stored for each row:
//   a[r*wid + c-r+kl] = mat[r][c]
bool LinBandSys::solve_band(void) {
  const int n   = size;
  const int wid = 2*kl + ku + 1;

  work_mat.assign( (std::size_t)n*wid, 0.0 );
  work_vec.resize( n );
  DType* a = work_mat.data();
  DType* b = work_vec.data();

  for ( int i = 0; i < n; i++ ) {
    const std::vector<int>&   cols = row_cols[i];
    const std::vector<DType>& vals = row_vals[i];
    const int r = new_idx[i];
    for ( std::size_t k = 0; k < cols.size(); k++ )
      a[(std::size_t)r*wid + new_idx[cols[k]]-r+kl] = vals[k];
    b[r] = vec_rhs[i];
  }

  // Elimination
  for ( int k = 0; k < n; k++ ) {
    const int r_max = std::min( k+kl, n-1 );    // last row with an entry in column k
    const int c_max = std::min( k+ku+kl, n-1 ); // last column with an entry in these rows

    // Pivot: the largest entry of column k
    int   p     = k;
    DType a_max = std::fabs( a[(std::size_t)k*wid + kl] );
    for ( int r = k+1; r <= r_max; r++ ) {
      const DType a_loc = std::fabs( a[(std::size_t)r*wid + k-r+kl] );
      if ( a_loc > a_max ) {
        a_max = a_loc;
        p     = r;
      }
    }
    if ( a_max == 0.0 ) return false; // Singular

    DType* row_k = a + (std::size_t)k*wid + kl - k; // row_k[c] = mat[k][c]
    if ( p != k ) {
      DType* row_p = a + (std::size_t)p*wid + kl - p;
      for ( int c = k; c <= c_max; c++ ) std::swap( row_k[c], row_p[c] );
      std::swap( b[k], b[p] );
    }

    const DType pivot = row_k[k];
    for ( int r = k+1; r <= r_max; r++ ) {
      DType* row_r = a + (std::size_t)r*wid + kl - r;
      if ( row_r[k] == 0.0 ) continue;
      const DType mult = row_r[k] / pivot;
      for ( int c = k+1; c <= c_max; c++ ) row_r[c] -= mult * row_k[c];
      b[r] -= mult * b[k];
    }
  }

  // Back substitution
  for ( int k = n-1; k >= 0; k-- ) {
    const int c_max = std::min( k+ku+kl, n-1 );
    const DType* row_k = a + (std::size_t)k*wid + kl - k;
    DType val = b[k];
    for ( int c = k+1; c <= c_max; c++ ) val -= row_k[c] * b[c];
    b[k] = val / row_k[k];
  }

  for ( int i = 0; i < n; i++ ) vec_sol[i] = b[new_idx[i]];

  return true;
}

// LU factorization with partial pivoting of the full matrix by LAPACK
bool LinBandSys::solve_full(void) {
  const int n = size;

  work_mat.assign( (std::size_t)n*n, 0.0 );
  work_vec.resize( n );
  work_piv.resize( n );

  // LAPACK stores a matrix column by column
  for ( int i = 0; i < n; i++ ) {
    const std::vector<int>&   cols = row_cols[i];
    const std::vector<DType>& vals = row_vals[i];
    for ( std::size_t k = 0; k < cols.size(); k++ )
      work_mat[(std::size_t)cols[k]*n + i] = vals[k];
    work_vec[i] = vec_rhs[i];
  }

  int nrhs = 1;
  int info = 0;
  lapack_gesv( &n, &nrhs, work_mat.data(), &n, work_piv.data(), work_vec.data(), &n, &info );
  if ( info != 0 ) return false;

  for ( int i = 0; i < n; i++ ) vec_sol[i] = work_vec[i];

  return true;
}

void LinBandSys::copyTo(LinMat& mat_ot, LinVec& rhs_ot) const {
  assert( mat_ot.getBlckDim() == 1 && rhs_ot.getBlckDim() == 1 );
  assert( mat_ot.getSize() == size && rhs_ot.getSize() == size );
  Entry elem(1,0.0);
  VComp comp(1,0.0);
  for ( int i = 0; i < size; i++ ) {
    const std::vector<int>&   cols = row_cols[i];
    const std::vector<DType>& vals = row_vals[i];
    for ( std::size_t k = 0; k < cols.size(); k++ ) {
      elem = vals[k];
      mat_ot.set( i, cols[k], elem );
    }
    comp = vec_rhs[i];
    rhs_ot.set( i, comp );
  }
  return;
}

void LinBandSys::copySolFrom(const LinVec& sol_in) {
  assert( sol_in.getBlckDim() == 1 && sol_in.getSize() == size );
  for ( int i = 0; i < size; i++ ) vec_sol[i] = sol_in.getGlobVal( i );
  return;
}
