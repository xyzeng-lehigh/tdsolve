#ifndef ALG_BANDSOLVE_H
#define ALG_BANDSOLVE_H

// Direct solution of a sparse linear system with scalar entries
//
// LinBandSys holds a square matrix, a right hand side and a solution vector.
// The matrix is stored row by row and only the entries that have been added
// are kept, so that assembling it costs neither a search in a tree nor any
// memory allocation once all entries have been visited for the first time
// (LinSparseMat of Alg_Solve.h keeps each scalar in an Entry, i.e. a
// separately allocated block matrix, inside a std::map).
//
// solve() is a direct solver that takes advantage of the sparsity:
//   1. The unknowns are renumbered by the reverse Cuthill-McKee algorithm,
//      which gathers the non-zero entries in a band around the diagonal.
//      E.g. the Jacobian of several species that are coupled cell by cell
//      and to the neighbouring cells becomes a band of a few diagonals,
//      independent of the number of cells. The numbering only depends on
//      the sparsity pattern and is reused until a new entry is added.
//   2. The banded matrix is factorized by Gaussian elimination with partial
//      (row) pivoting. The cost is proportional to size * bandwidth^2,
//      instead of size^3 for a full matrix.
//   3. If the band is not narrow enough to pay off (e.g. the Jacobian of a
//      velocity that is an integral over all cells), the full matrix is
//      factorized by LAPACK instead (xGESV).
//
// Typical usage (see TDSolveSacado.cc):
//   LinBandSys sys(n);
//   sys.reset();                 // start a new system, same sparsity
//   sys.add(i, j, val);          // mat[i][j] += val
//   sys.rhs(i) = val;
//   sys.solve();
//   x_i = sys.sol(i);
//
// The iterative solvers and preconditioners of LinSolve (Alg_Solve.h) work
// on LinMat/LinVec; copyTo() and copySolFrom() transfer the system there
// and the solution back.

#include <vector>
#include <algorithm>
#include "Data_Type.h"

class LinMat;
class LinVec;

class LinBandSys {
  public:
    explicit LinBandSys(const int size_in);
    ~LinBandSys(void) {}

    int getSize(void) const { return size; }

    // Set the matrix, the right hand side and the solution to zero.
    // The entries of the matrix stay registered, i.e. the sparsity pattern
    // (and with it the numbering used by the solver) is kept.
    void reset(void);

    // mat[i][j] += val; the entry is registered if it is a new one
    inline void add(const int i, const int j, const DType val);

    // mat[i][j]; zero if the entry is not registered
    DType get(const int i, const int j) const;

    // Right hand side and solution
    DType& rhs(const int i)       { return vec_rhs[i]; }
    DType  rhs(const int i) const { return vec_rhs[i]; }
    DType  sol(const int i) const { return vec_sol[i]; }

    // || rhs ||_2, || sol ||_inf and || rhs - mat * sol ||_2
    DType getRhsNorm(void) const;
    DType getSolInfNorm(void) const;
    DType getResNorm(void) const;

    // Solve mat * sol = rhs. Returns false if the matrix is singular, in
    // which case the solution is set to zero (a warning is printed the first
    // time this happens).
    bool solve(void);

    // Information on the system and on the last call of solve()
    int  getNumEntries(void) const;                // registered entries
    int  getLowerBand(void) const { return kl; }   // sub-diagonals after renumbering
    int  getUpperBand(void) const { return ku; }   // super-diagonals after renumbering
    bool isBanded(void) const { return use_band; } // banded (true) or full (false) solver

    // Exchange with the general storage of Alg_Solve.h (block dimension 1):
    // copyTo sets mat_ot[i][j] for all registered entries and rhs_ot,
    // copySolFrom fetches the solution computed by a LinSolve.
    void copyTo(LinMat& mat_ot, LinVec& rhs_ot) const;
    void copySolFrom(const LinVec& sol_in);

  private:
    int size;

    // Row i of the matrix: column indices in increasing order and values
    std::vector< std::vector<int> >   row_cols;
    std::vector< std::vector<DType> > row_vals;

    std::vector<DType> vec_rhs;
    std::vector<DType> vec_sol;

    // Renumbering; valid if !pattern_changed
    bool pattern_changed;
    bool use_band;
    bool warned_singular;
    int  kl, ku;
    std::vector<int> new_idx; // new_idx[i]: number of unknown i in the banded system

    // Work space of the solvers
    std::vector<DType> work_mat;
    std::vector<DType> work_vec;
    std::vector<int>   work_piv;

    void setup_numbering(void);
    void calc_bandwidth(void);
    bool solve_band(void);
    bool solve_full(void);
};

inline void LinBandSys::add(const int i, const int j, const DType val) {
  std::vector<int>& cols = row_cols[i];
  const std::vector<int>::iterator it = std::lower_bound( cols.begin(), cols.end(), j );
  const std::size_t k = it - cols.begin();
  if ( it == cols.end() || *it != j ) {
    cols.insert( it, j );
    row_vals[i].insert( row_vals[i].begin()+k, 0.0 );
    pattern_changed = true;
  }
  row_vals[i][k] += val;
  return;
}

#endif
