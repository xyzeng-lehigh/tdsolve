#ifndef ALG_SOLVE_H
#define ALG_SOLVE_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include "Data_Type.h"
#include "Comput_Choices.h"
#include "LowDim_Alg.h"
#include "HighDim_Alg.h"
#include "Mat_Entry.h"
#include "LeastSquares_Solve.h"
using namespace std;

// Index pair object
class sm_idx {
    public:
        int i;
        int j;

        sm_idx(void) {
            sm_idx(0,0);
        }

        sm_idx(int i_in, int j_in) {
            i = i_in;
            j = j_in;
        }

        bool operator==(const sm_idx &other) const {
            return ((i==other.i)&&(j==other.j));
        }

        bool operator<(const sm_idx &other) const {
            return ((i<other.i)||
                    ((i==other.i)&&(j<other.j)));
        }
};

/*
 * Column vector storage
 */
class LinVec {
    protected:
        int     size;   // Number of components
        int     dim;    // Size of each component

        VComp** vec;

    public:
        LinVec(void);
        LinVec(const int dim_in, const int size_in) {
            size = -1;
            dim  = -1;
            vec  = NULL;
            setBlckDim(dim_in);
            setSize(size_in);
        }
        ~LinVec(void);

        int getSize(void)    const {return size;}
        int getBlckDim(void) const {return dim;}
        // Always set dim before setSize!!
        void setBlckDim(int dim_in) {
            assert(dim_in>0);
            dim = dim_in;
            return;
        }
        void setSize(int);
        VComp  get(int i) const {return *vec[i];}
        void   set(int, VComp&);
        void reset(int, DType val=0.0);
        void reset(DType val);
        DType& getGlob(int n) {
            int i_blck = n/dim;
            int i_comp = n-i_blck*dim;
            return (*vec[i_blck]).get(i_comp);
        }
        const DType& getGlob(int n) const {
            int i_blck = n/dim;
            int i_comp = n-i_blck*dim;
            return (*vec[i_blck]).get(i_comp);
        }
        DType getGlobVal(int n) const {
            int i_blck = n/dim;
            int i_comp = n-i_blck*dim;
            return (*vec[i_blck]).get(i_comp);
        }

        // 2-norm
        DType getNorm(void);
        // this = vec_in
        void copy(LinVec& vec_in);
        // this = a*this
        void scale(const DType &);
        // this = a*this + b*vec
        void scale_inc(const DType&,const DType&,LinVec& vec_in);
        // this = (I-2vec vec^T)*this
        // vec_proj needs to be normalized
        void proj(LinVec& vec_proj);

};

// Inner product
DType calc_inner_prod(const LinVec& vec_1, const LinVec& vec_2);

/*
 * Matrix storage
 */

class LinMat {
    protected:
        int     size_r; // Number of rows of the block matrix
        int     size_c; // Number of columes of the block matrix
        int     size;   // Only used for square matrices
        int     dim;    // Size of each block (a square block)
        std::map<sm_idx,Entry> dummy_mat;

        // mat_ot is different from either this or mat_in; most straightforward implementation
        virtual void xMat_new(LinMat*& mat_in, LinMat*& mat_ot);
        // mat_ot is the same as mat_in; this is a square matrix but different from the other two
        virtual void xMat_pst(LinMat*& mat_io);
        // mat_ot is the same as this, but different from mat_in that is square
        virtual void xMat_pre(LinMat*& mat_in);
        // mat = mat X mat ...
        // Currently, a temp matrix is created and 
        // use the subroutine xMat_new
        // Thus this is pure virtual
        virtual void xMat_sqr(void) = 0;

        // vec_ot is different from vec_in
        virtual void xVec_new(LinVec*& vec_in, LinVec*& vec_ot);
        virtual void xVec_t_new(LinVec*& vec_in, LinVec*& vec_ot);
        // vec_ot is the same as vec_in
        virtual void xVec_rpl(LinVec*& vec_io);
        virtual void xVec_t_rpl(LinVec*& vec_io);

    public:
        enum    matType {LINMAT_FULL=0, LINMAT_SPARSE, LINMAT_TRID} mtype;
        enum    dimType {LINMAT_GEN=0, LINMAT_SQ=1} dtype;
        bool isFull(void) const {return mtype==LINMAT_FULL;}
        bool isSparse(void) const {return mtype==LINMAT_SPARSE;}
        bool isTriD(void) const {return mtype==LINMAT_TRID;}
        bool isGen(void) const {return dtype==LINMAT_GEN;}
        bool isSquare(void) const {return dtype==LINMAT_SQ;}

        LinMat(void) {size_r=-1; size_c=-1; dim=-1; dummy_mat.clear();}
        virtual ~LinMat(void) {}

        int  getSizeCol(void) const {return size_r;}
        int  getSizeRow(void) const {return size_c;}
        int  getBlckDim(void) const {return dim;}
        // Always set dim before setSize!!
        void setBlckDim(int dim_in) {
            assert(dim_in>0);
            dim = dim_in; 
            return;
        }

        virtual DType getGlobVal(int i, int j) {
            int i_blck = i/dim;
            int j_blck = j/dim;
            int i_comp = i-i_blck*dim;
            int j_comp = j-j_blck*dim;
            return (get(i_blck,j_blck).getVal())[i_comp][j_comp];
        }

        matType getType(void) {return mtype;}
        virtual void setSize(int,int);
        virtual void setEqualSize(int i) {
            size = i;
            this->setSize(i,i);
            return;
        }
        int getSize(void) const {return size;}
        // get entry of mat[i][j]
        virtual Entry get(int,int) = 0;
        // set entry of mat[i][j]
        virtual void   set(int,int,Entry&) = 0;
        // reset entry of mat[i][j] to a certain value
        virtual void reset(int,int,DType val_diag=0.0) = 0;
        // reset entire row mat[i] to the same value
        // For sparse matrix, this will not erase entries
        virtual void reset(int,DType val_diag=0.0) = 0;
        // reset main diagonal to multiple of identity 
        virtual void reset(DType val_diag=0.0) = 0;
        // mat.xMat(mat_in,mat_ot);
        // mat_ot = mat X mat_in
        // The dimension rule applies
        // This support various storage options
        virtual void xMat(LinMat*& mat_in, LinMat*& mat_ot);
        // mat.xVec(vec_in,vec_ot);
        // vec_ot = mat X vec_in
        // The dimension rule applies
        // This support various storage options
        virtual void xVec(LinVec*& vec_in, LinVec*& vec_ot);
        // mat.xVec_t(vec_in,vec_ot);
        // vec_ot = transpose(mat) X vec_in
        // The dimension rule applies
        // This support various storage options
        virtual void xVec_t(LinVec*& vec_in, LinVec*& vec_ot);

        // Only for LinSparseMat
        virtual int getCol(std::map<sm_idx,Entry>&,int) {
            return 0;}
        virtual int getRow(std::map<sm_idx,Entry>&,int) {
            return 0;}
        virtual std::map<sm_idx,Entry>::iterator it_beg(void) {
            return dummy_mat.begin();}
        virtual std::map<sm_idx,Entry>::iterator it_end(void) {
            return dummy_mat.end();}

        // Debug purpose
        void print(void) {
            for ( int i_r = 0; i_r < size_r; i_r++ ) 
                for ( int i_comp = 0; i_comp < dim; i_comp++ ) {
                    int i = i_r * dim + i_comp;
                    for ( int j_c = 0; j_c < size_c; j_c++ ) 
                        for ( int j_comp = 0; j_comp < dim; j_comp++ ) {
                            int j = j_c * dim + j_comp;
                            cout << getGlobVal(i,j) << " ";
                        }
                    cout << endl;
                }
        }

        void print_skip(int skip_i, int skip_j) {
            for ( int s_i = 0; s_i < skip_i; s_i++ )
                for ( int i_r = s_i; i_r < size_r; i_r += skip_i )
                    for ( int i_comp = 0; i_comp < dim; i_comp++ ) {
                        int i = i_r * dim + i_comp;
                        for ( int s_j = 0; s_j < skip_j; s_j++ )
                            for ( int j_c = s_j; j_c < size_c; j_c += skip_j )
                                for ( int j_comp = 0; j_comp < dim; j_comp++ ) {
                                    int j = j_c * dim + j_comp;
                                    cout << getGlobVal(i,j) << " ";
                                }
                        cout << endl;
                    }
        }

        void print_reverse(int skip_i, int skip_j) {
            for ( int s_i = skip_i-1; s_i > -1; s_i-- )
                for ( int i_r = s_i; i_r < size_r; i_r += skip_i )
                    for ( int i_comp = 0; i_comp < dim; i_comp++ ) {
                        int i = i_r * dim + i_comp;
                        for ( int s_j = skip_j-1; s_j > -1; s_j-- )
                            for ( int j_c = s_j; j_c < size_c; j_c += skip_j )
                                for ( int j_comp = 0; j_comp < dim; j_comp++ ) {
                                    int j = j_c * dim + j_comp;
                                    cout << getGlobVal(i,j) << " ";
                                }
                        cout << endl;
                    }
        }
};

// Square matrices
class LinSqMat : public LinMat {
    protected:
        void xMat_sqr(void) {
            fprintf(stderr,"Error: xMat_sqr should not be called for a LinSqMat instance\n");
            exit(-1);
        }
    public:
        LinSqMat(void) {dtype = LINMAT_SQ;}
        ~LinSqMat(void) {};

        int getSize(void) const {return size;}

        virtual void setEqualSize(int);

        Entry get(int i_r, int i_c) {
            fprintf(stderr,"Error: get(int,int) should not be called for a LinSqMat instance\n");
            exit(-1);
        }
        void set(int i_r, int i_c, Entry& ent_in) {
            fprintf(stderr,"Error: set(int,int,Entry&) should not be called for a LinSqMat instance\n");
            exit(-1);
        }
        void reset(int i_r, int i_c, DType val_in=0.0) {
            fprintf(stderr,"Error: reset(int,int,DType&) should not be called for a LinSqMat instance\n");
            exit(-1);
        }
        void reset(int i_r, DType val_in=0.0) {
            fprintf(stderr,"Error: reset(int,DType&) should not be called for a LinSqMat instance\n");
            exit(-1);
        }
        void reset(DType val_in=0.0) {
            fprintf(stderr,"Error: reset() should not be called for a LinSqMat instance\n");
            exit(-1);
        }
};

// Full matrix
class LinFullMat : public LinMat {
    protected:
        Entry**  mat;

        void xMat_sqr(void);
    public:
        LinFullMat(void);
        ~LinFullMat(void);

        void  setSize(int,int);
        Entry get(int,int);
        void  set(int,int,Entry&);
        void  reset(int,int,DType val_in=0.0);
        void  reset(int,DType val_in=0.0);
        void  reset(DType val_in=0.0);
};

// Tri-diagonal matrix
class LinTriMat : public LinSqMat {
    protected:
        Entry** lower;
        Entry** upper;
        Entry** mdiag;

        void xMat_new(LinMat*&,LinMat*&);
        void xMat_pst(LinMat*&);
        void xMat_pre(LinMat*&);
        void xMat_sqr(void);

        void xVec_new(LinVec*&,LinVec*&);
        void xVec_t_new(LinVec*&,LinVec*&);
        //void xVec_rpl(LinVec*&);
        //void xVec_t_rpl(LinVec*&);
    public:
        LinTriMat(void);
        ~LinTriMat(void);

        void  setEqualSize(int);
        Entry get(int,int);
        void  set(int,int,Entry&);
        void  reset(int,int,DType val_in=0.0);
        void  reset(int,DType val_in=0.0);
        void  reset(DType val_in=0.0);
};

// Sparse matrix
class LinSparseMat : public LinMat {
    protected:
        // This map is used to store the entries
        // during matrix construction
        // It may not be appropriate for some direct
        // solve, in which case the data structure
        // need to be transformed
        std::map<sm_idx,Entry> mat;

        void xMat_new(LinMat*&,LinMat*&);
        void xMat_pst(LinMat*&);
        void xMat_pre(LinMat*&);
        void xMat_sqr(void);

        void xVec_new(LinVec*&,LinVec*&);
        void xVec_t_new(LinVec*&,LinVec*&);
        //void xVec_rpl(LinVec*&);
        //void xVec_t_rpl(LinVec*&);
    public:
        LinSparseMat(void);
        ~LinSparseMat(void) {};

        void  setSize(int,int);
        // If the computed entry is zero, 
        // entry.sig = true
        Entry get(int,int);
        void  set(int,int,Entry&);
        void  reset(int,int,DType val_in=0.0);
        void  reset(int,DType val_in=0.0);
        void  reset(DType val_in=0.0);

        // Get a column, number of non-zero entries returned
        int getCol(std::map<sm_idx,Entry>&,int);
        // Get a row, number of non-zero entries returned
        int getRow(std::map<sm_idx,Entry>&,int);

        // Get the initial/final iterator
        std::map<sm_idx,Entry>::iterator it_beg(void) {
            return mat.begin();}
        std::map<sm_idx,Entry>::iterator it_end(void) {
            return mat.end();}
};

// Square full matrix
class LinSqFullMat : public LinFullMat {
    protected:
    public:
        LinSqFullMat(void) {dtype = LINMAT_SQ;}
        ~LinSqFullMat(void) {};

        void setEqualSize(int i) {
            size = i;
            LinFullMat::setSize(i,i);
            return;
        }
};

// Square sparse matrix
class LinSqSparseMat : public LinSparseMat {
    protected:
    public:
        LinSqSparseMat(void) {dtype = LINMAT_SQ;}
        ~LinSqSparseMat(void) {};

        void setEqualSize(int i) {
            size = i;
            LinSparseMat::setSize(i,i);
            return;
        }
};

/*
 * Linear solve of PAx=Pb (P is preconditioner) 
 *  A : mat -- size X size; 
 *  x : sol -- size X 1; 
 *  b : rhs -- size X 1
 *  P : prc -- size X size
 * Memory strategies
 *  (F) full   storage -- mat is size X size unit of memory
 *  (S) sparse storage -- mat only have nonzeros (dynamic size)
 * Solver strategies
 *  (D) decomposition-based -- equivalent to finding the inverse of PA
 *  (I) iterative-based     -- avoiding any inverse
 * Combinations
 *  (F-D), (F-I), (S-I); 
 *  (S-D) may be available (like tri-diagonal matrices)
 * Available solvers (to double check)
 *  CG (PCG) -- (S-I);
 *  GS       -- (F-D);
 */
class LinSolve {
    protected:
        int     size;
        int     dim;     // Dimension of each block
        bool    preCond; // if to apply left preconditioner
        LinMat*    mat;     // size X size
        LinVec*    sol;     // size X 1
        LinMat*    prc;     // size X size (left)
        LinVec*    rhs;     // size X 1
        LinVec*    res;     // size X 1
        // mat and prc are full or sparse

        bool    print_stat;

        // Control parameters
        int     max_iter;
        int     prec_ilu_fill;
        DType   target_res;
        DType   rhs_nrm_ini;
        std::vector<DType> abs_res;     // ||res||
        std::vector<DType> rel_res_ini; // ||res||/||res_0||
        std::vector<DType> rel_res_rhs; // ||res||/||rhs||
        DType   coef_abs, coef_ini, coef_rhs;
        // res =  coef_abs * abs_res + coef_ini * rel_res_ini
        //      + coef_rhs * rel_res_rhs
        std::vector<DType> res_hist;

        ChoiceLSysType   lsysType;
        ChoiceLSysIter   lsysIter;
        ChoiceLSysPrec   lsysPrec;
        ChoiceLSysKrylov lsysKrylov;

        typedef void (LinSolve::*funIter)(DType);
        typedef void (LinSolve::*funKrylov)(void);
        typedef void (LinSolve::*funPrec)(int);

        // Copmute the residual
        // res = prc * ( rhs - mat * sol )
        DType compute_res(void);
        DType compute_tot_res(DType);
        // Update various measures of residuals
        bool update_res(DType res_in = -1.0);
        // Check if to continue by iteration number
        bool check_max_iter(int);

        // Need to be performed when both
        // mat, rhs and prc are setup correctly if prc is
        // linked from outside, or
        // compute prc first (the space is still needed)
        // and then apply to mat, rhs
        funPrec apply_prec;
        void apply_prec_link(int dummy=0);
        void apply_prec_eye(int dummy=0);
        void apply_prec_ilu_0(int dummy_level=0);
        void apply_prec_ilu_p(int fill_level);
        void apply_prec_bilu_0(int dummy_level=0);
        void apply_prec_bilu_p(int fill_level);
        // ilu_nzpat[i][k] = j if (i,j) is non-zero
        // each ilu_nzpat[i] is in increasing order
        // (i,i) is always listed inside
        std::vector<int>* ilu_nzpat;
        // for k=0; k<ilu_idx[i], ilu_nzpat[i][k]<i
        // for k=ilu_idx[i], ilu_nzpat[i][k] = i
        // for k=ilu_idx[i]+1, k<ilu_nzpat[i].size()
        //     ilu_nzpat[i][k]>i
        int* ilu_idx;

        void calc_nz_pattern(LinMat*&,std::vector<int>*&,int*&,DType);
        void calc_nz_bpattern(LinMat*&,std::vector<int>*&,int*&,DType);

        // Start the olver
        bool solve_start(void);

        // Solver 0: direct method (stable GS)
        void solve_dir(bool,int);

        // Solver 1 : stationary iterative solver
        void solve_iter(bool,int);

        funIter solve_iter_single;

        void solve_iter_jac(DType dummy);

        void solve_iter_gs(DType dummy);
        void solve_iter_bgs(DType dummy);

        DType sor_par;
        void solve_iter_sor(DType);
        void solve_iter_bsor(DType);

        // Solver 2 : Krylov space solver
        void solve_krylov(bool,int);
        LstSqr lsSolve;

        funKrylov solve_krylov_pre;
        funKrylov solve_krylov_single;

        // Maximum allowed dimension of Krylov space
        int  kspace_dim;
        // If the residual vector is computed intermediately
        bool kspace_res;
        // Signal of exiting
        bool kspace_signal;

        std::vector<LinVec*> vec_v;  // DQGMRES, LANCZOS, QMR
        std::vector<LinVec*> vec_w;  // GMRES, LANCZOS, QMR
        std::vector<LinVec*> vec_p;  // DQGMRES, QMR
        void clear_krylov_spaces(std::vector<LinVec*>&);

        // GMRES: 
        // kspace_dim=0:
        //   GMRES (Full GMRES) 
        // kspace_dim>0: 
        //   DQGMRES (Direct quasi-GMRES)
        // Both with Householder orthogonalization 
        // for enhancing robustness
        LinVec*  gmres_g;
        LinVec*  gmres_z;
        LinVec*  gmres_v;
        UTriMAT* gmres_h;
        LinVec*  gmres_sol0;
        std::vector<DType> gmres_cos;
        std::vector<DType> gmres_sin;
        void solve_krylov_gmres_pre(void);
        void solve_krylov_dqgmres_pre(void);
        DType gmres_gamma;
        void solve_krylov_gmres(void);
        void solve_krylov_dqgmres(void);

        // CG: BICGSTAB (Biconjugate CG Stabilized)
        LinVec* vec_cg_r;    // size X 1
        LinVec* vec_cg_p;    // size X 1
        void solve_krylov_cg_pre(void);
        void solve_krylov_cg(void);

        // LANCZOS: restarted two-sided Lanczos
        std::vector<DType> lanczos_alpha;
        std::vector<DType> lanczos_delta;
        std::vector<DType> lanczos_beta;
        DType   lanczos_beta0;
        LinVec* lanczos_sol0;
        void solve_krylov_lanczos_pre(void);
        void solve_krylov_lanczos(void);
        // QMR: restarted quasi-minimal residual
        DType qmr_gamma;
        DType qmr_cos[2];
        DType qmr_sin[2];
        DType qmr_beta;
        DType qmr_delta;
        int   qmr_iter;
        void solve_krylov_qmr_pre(void);
        void solve_krylov_qmr(void);

    public:
        LinSolve(void);
        virtual ~LinSolve(void);

        HVEC convert(LinVec*&);
        HMAT convert(LinMat*&);
        void convert(LinVec*&, HVEC, int, int);
        void convert(LinMat*&, HMAT, int, int, DType conv_tol=TOL_ZERO);

        void setSORPar(DType par_in) {
            sor_par = par_in;
            return;
        }

        void setSolveDir() {
            lsysType = LSYS_TYPE_DIR;
            return;
        }

        void setSolveIter(ChoiceLSysIter lsysIter_in) {
            lsysType = LSYS_TYPE_ITER;
            lsysIter = lsysIter_in;
            switch (lsysIter) {
                case LSYS_ITER_JAC:
                    solve_iter_single = &LinSolve::solve_iter_jac;
                    break;
                case LSYS_ITER_GS:
                    solve_iter_single = &LinSolve::solve_iter_gs;
                    break;
                case LSYS_ITER_BGS:
                    solve_iter_single = &LinSolve::solve_iter_bgs;
                    break;
                case LSYS_ITER_SOR:
                    solve_iter_single = &LinSolve::solve_iter_sor;
                    break;
                case LSYS_ITER_BSOR:
                    solve_iter_single = &LinSolve::solve_iter_bsor;
                    break;
                default:
                    fprintf(stderr,"Warning: LinSolve the iterative solver is not defined, will use Gauss-Seidel\n");
                    solve_iter_single = &LinSolve::solve_iter_gs;
                    break;
            }
            return;
        }

        void setSolvePrec(ChoiceLSysPrec lsysPrec_in, int prec_ilu_fill_in=0) {
            lsysPrec        = lsysPrec_in;
            prec_ilu_fill   = prec_ilu_fill_in;
            preCond         = true;
            switch (lsysPrec) {
                case LSYS_PREC_OFF:
                    apply_prec = NULL;
                    preCond    = false;
                    break;
                case LSYS_PREC_LINK:
                    apply_prec = &LinSolve::apply_prec_link;
                    break;
                case LSYS_PREC_EYE:     // Identity preconditioner (test only)
                    apply_prec = &LinSolve::apply_prec_eye; 
                    break;
                case LSYS_PREC_ILU:
                    if (prec_ilu_fill==0)
                        apply_prec = &LinSolve::apply_prec_ilu_0;
                    else
                        apply_prec = &LinSolve::apply_prec_ilu_p;
                    break;
                case LSYS_PREC_BILU:
                    if (prec_ilu_fill==0)
                        apply_prec = &LinSolve::apply_prec_bilu_0;
                    else
                        apply_prec = &LinSolve::apply_prec_bilu_p;
                    break;
            }
            return;
        }

        void setSolveKrylov(ChoiceLSysKrylov lsysKrylov_in, int kspace_dim_in) {
            lsysType   = LSYS_TYPE_KRYLOV;
            lsysKrylov = lsysKrylov_in;
            kspace_dim = max(0,kspace_dim_in);
            kspace_res = true;
            switch (lsysKrylov) {
                case LSYS_KRYLOV_CG:
                    solve_krylov_pre    = &LinSolve::solve_krylov_cg_pre;
                    solve_krylov_single = &LinSolve::solve_krylov_cg;
                    break;
                case LSYS_KRYLOV_GMRES:
                    if (kspace_dim==0) {
                        solve_krylov_pre    = &LinSolve::solve_krylov_gmres_pre;
                        solve_krylov_single = &LinSolve::solve_krylov_gmres;
                    }
                    else {
                        assert(kspace_dim>1);
                        kspace_res          = false;
                        solve_krylov_pre    = &LinSolve::solve_krylov_dqgmres_pre;
                        solve_krylov_single = &LinSolve::solve_krylov_dqgmres;
                    }
                    break;
                case LSYS_KRYLOV_LANCZOS:
                    solve_krylov_pre    = &LinSolve::solve_krylov_lanczos_pre;
                    solve_krylov_single = &LinSolve::solve_krylov_lanczos;
                    break;
                case LSYS_KRYLOV_QMR:
                    solve_krylov_pre    = &LinSolve::solve_krylov_qmr_pre;
                    solve_krylov_single = &LinSolve::solve_krylov_qmr;
                    break;
                default:
                    fprintf(stderr,"Warning: LinSolve the Krylov space solver is not defined, will use GMRES\n");
                    solve_krylov_pre    = &LinSolve::solve_krylov_gmres_pre;
                    solve_krylov_single = &LinSolve::solve_krylov_gmres;
                    break;
            }
            return;
        }

        std::vector<DType>& getAbsResHist(void) {return abs_res;}
        std::vector<DType>& getRelResIniHist(void) {return rel_res_ini;}
        std::vector<DType>& getRelResRhsHist(void) {return rel_res_rhs;}
        std::vector<DType>& getResHist(void) {return res_hist;}

        void setTargetRes(DType res_in) {
            target_res = (res_in>0.0) ? res_in : TOL_ZERO;
            return;
        }
        void setMaxIter(int max_iter_in) {
            max_iter = (max_iter_in>=0) ? max_iter_in : -1;
            return;
        }

        void setResCoefAbs(DType coef_in) {
            coef_abs = coef_in; return;}
        void setResCoefIni(DType coef_in) {
            coef_ini = coef_in; return;}
        void setResCoefRhs(DType coef_in) {
            coef_rhs = coef_in; return;}

        void setSize(int size_in) {
            assert(size_in>0);
            size = size_in;
            return;
        }
        void setBlckDim(int dim_in) {
            assert(dim_in>0);
            dim = dim_in;
            return;
        }
        void setPrintStatus(bool prt_in) {
            print_stat = prt_in;
            return;
        }
        virtual void linkMat(LinMat*&); 
        virtual void linkPrc(LinMat*&);
        void linkSol(LinVec*&);
        void linkRhs(LinVec*&);
        void linkRes(LinVec*&);

        // return whether target residual is met
        bool solve(void);
};

#endif
