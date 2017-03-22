#ifndef LEASTSQUARES_SOLVE_H
#define LEASTSQUARES_SOLVE_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <vector>
#include "Data_Type.h"
#include "LowDim_Alg.h"
#include "HighDim_Alg.h"
using namespace std;

/*
 * General least squares solve
 * To solve
 *      min ||Ax-b||_(2,W) such that Cx=d
 * Entities
 *      A : m x n (known)
 *      x : n x 1 (unknown)
 *      b : m x 1 (known)
 *      C : k x n (known)
 *      d : k x 1 (known)
 * Weighted norm (stored as a vector)
 *      W : m x m (known)
 * Lagrangian multiplier
 *      l : k x 1 (unknown)
 * Method: to find (x,l) such that
 *      (Ax-b)^TW^TW(Ax-b) + 2l^T(Cx-d)
 *
 * Case 1: k=0 (no constraint)
 *      H = A^TW^2A
 *      r = A^TW^2b
 *      x = H.inv() * r
 *
 * Case 2: k>0
 * i.e. the saddle point problem
 *      [ H    C^T ][ x ]   [ r ]
 *      [ C     0  ][ l ] = [ d ]
 * or equivalently
 *      l = (CH.inv()C^T).inv() * [ CH.inv()r - d ]
 *      x = H.inv() * [ r - C^Tl ]
 *
 * Remark: exact (LU) solve is used to compute
 * the inverse of matrices; thus n and k should be small
 */

class LstSqr {
    protected:
        int Nvars;      // Number of unkonws (n)
        int Neqns;      // Number of equations (m)
        int Ncstr;      // Number of constraints (k)

        HVEC*  sol;     // Nvars X 1 (x)
        HVEC*  rhs;     // Neqns X 1 (b)
        HVEC** mat;     // Neqns X Nvars (A)
        HVEC*  wgt;     // Neqns X 1 (W.diag)

        HMAT*  coef_mat;    // Nvars X Nvars (H)
        HMAT*  coef_inv;    // Nvars X Nvars (H.inv)
        HVEC*  coef_rhs;    // Nvars X 1 (r)

        bool   use_cstr;    // Add constraints ?
        HVEC** c_mat;       // Ncstr X Nvars (C)
        HVEC*  c_rhs;       // Ncstr X 1 (d)
        HVEC*  c_lag;       // Ncstr X 1 (l)

        HMAT* cstr_mat; // Ncstr X Ncstr (CH.inv C^T)
        HMAT* cstr_inv; // Ncstr X Ncstr (above inv)
        HVEC* cstr_rhs; // Ncstr X 1 (CH.inv r - d)

        void clear(void);      // Clear arrays

        void cstr_solve(void); // Given
                               //      c_mat
                               //      c_rhs
                               // and
                               //      coef_inv
                               //      coef_rhs
                               // Compute
                               //      cstr_mat
                               //      cstr_inv
                               //      cstr_rhs
                               //      c_lag

    public:
        LstSqr(void);
        ~LstSqr(void);

        void setNvars(int Nvars_in) {
            assert(Nvars_in>0);
            Nvars = Nvars_in;
            return;
        }
        void setNeqns(int Neqns_in) {
            assert(Neqns_in>0);
            Neqns = Neqns_in;
            return;
        }
        void setNcstr(int Ncstr_in) {
            if (Ncstr_in>0) {
                Ncstr = Ncstr_in;
                use_cstr = true;
            }
            else {
                Ncstr = 0;
                use_cstr = false;
            }
            return;
        }

        // Reset everything to initial status
        void reset(void);
        // Actual solve step, typical usage:
        // reset()
        // setNvars()
        // setNeqns()
        // setNcstr()
        // setup()
        // set weights/eqns/constraints, etc
        // solve()
        void solve(void);

        // Must be performed after setting N***'s
        // and before loading any entries
        // This step allocate spaces
        void setup(void);

        void setEqnsWeight(int,DType);

        void setEqnsCoef(int,int,DType);
        void setEqnsCoef(int,HVEC*);
        void setEqnsRhs(int,DType);
        void setEqns(int,HVEC*,DType);

        void setCstrCoef(int,int,DType);
        void setCstrCoef(int,HVEC*);
        void setCstrVal(int,DType);
        void setCstr(int,HVEC*,DType);

        DType getSol(int i) const {
            return (*sol)[i];
        }
        HVEC getSol() const {
            return *sol;
        }
};

/*
 * Local polynomial coefficients
 */

// Base class

class PolyCoef {
    protected:
        int     p_order;    // Order of polynomials
        DType** p_coef;     // Polynomial coefficients
        int*    n_coef;     // Number of coefficients for each order
        int     Ncoef;      // Number of coefficients
    public:
        PolyCoef(void);
        ~PolyCoef(void);

        void clear(void);
        virtual void setup(int);

        int     getNcoef(void) {return Ncoef;}
        int     getOrder(void) {return p_order;}

        void    setLocCoef(DType,int,int);
        void    setGlobCoef(DType,int);
};

// 1D case

class PolyCoef1D : public PolyCoef {
    public:
        PolyCoef1D(void) {}
        ~PolyCoef1D(void) {}

        void    setup(int);
        DType   getCoef(int);
        int     getNumCoef(int);

        int     getLocCoefIdx(int);
        int     getGlobCoefIdx(int);
        void    globIdx2Der(int,int&);

        void    setCoef(DType,int);
};

// 2D case

class PolyCoef2D : public PolyCoef {
    public:
        PolyCoef2D(void) {}
        ~PolyCoef2D(void) {}

        void    setup(int);
        DType   getCoef(int,int);
        int     getNumCoef(int,int);

        int     getLocCoefIdx(int,int);
        int     getGlobCoefIdx(int,int);
        void    globIdx2Der(int,int&,int&);

        void    setCoef(DType,int,int);
};

// 3D case

class PolyCoef3D : public PolyCoef {
    public:
        PolyCoef3D(void) {}
        ~PolyCoef3D(void) {}

        void    setup(int);
        DType   getCoef(int,int,int);

        int     getLocCoefIdx(int,int,int);
        int     getGlobCoefIdx(int,int,int);
        void    globIdx2Der(int,int&,int&,int&);

        void    setCoef(DType,int,int,int);
};

/*
 * Local least squares solver for polynomial coefficients
 */

// Base class

class LocLstSqr {
    protected:
        int     p_order;    // Order of polynomials 
        int     Npars;      // Number of unknowns
        bool    use_val;    // Nodal val @ center ?
        // coef_mat * coef_par = coef_rhs
        // coef_par : unknowns
        HVEC*   coef_rhs;   // Npars X 1   
        HVEC*   coef_par;   // Npars X 1
        HMAT*   coef_mat;   // Npars X Npars
        HMAT*   coef_inv;   // Npars X Npars -- inv of coef_mat

        // The linear constraints 
        // (implemented by Lagrangian multipliers)
        // Added as nodal constraint
        int     Ncstr;      // Number of constraints
        bool    use_cstr;   // Add constraints ?
        std::vector<DType> cstr_vars;
                            // Ncstr X Npars vectors
        std::vector<DType> cstr_vals;
                            // Ncstr vectors
        HMAT*   cstr_mat;   // Ncstr X Ncstr
        HMAT*   cstr_inv;   // Ncstr X Ncstr
        HVEC*   cstr_lag;   // Ncstr X 1 (Lagrangian multiplier)
        HVEC*   cstr_rhs;   // Ncstr X 1

        void add2Arrays(DType,DType*,DType);
        void add2Mat(DType*,DType);
        void add2Rhs(DType,DType*,DType);
        void add2CstrLists(DType,DType*);

        void cstr_solve(void);  // Given
                                //      cstr_vars
                                //      cstr_vals
                                // and
                                //      coef_inv
                                //      coef_rhs
                                // Compute
                                //      cstr_mat
                                //      cstr_inv
                                //      cstr_rhs
                                //      cstr_lag
    public:
        LocLstSqr(void);

        void setOrder(int order_in, bool val_in) {
            assert(order_in>=0);
            p_order = order_in;
            use_val = val_in;
            return;
        }

        void clear(void);
        void resetMats(void);
        void resetRhs(void);
        virtual void solve(void) = 0;

        void cstr_reset(void);
        void cstr_start(void);  // The same as above
        void cstr_end(void);    // Check consistency
};

// 1D case

class LocLstSqr1D : public LocLstSqr {
    protected:
        PolyCoef1D* polyCoef;
        DType       x_c;
        DType       val_c;
    public:
        void setup(int order, PolyCoef1D*, bool val_in);
        void setCenterX(DType x_c_in) {x_c=x_c_in;}
        void setCenterV(DType);

        void addNodalMat(DType,
                DType weight=1.0);
        void addCellMat(DType,DType,
                DType weight=1.0);

        void addNodalVal(DType,DType,
                DType weight=1.0);
        void addCellVal(DType,DType,DType, 
                DType weight=1.0);

        void addNodalArrays(DType,DType,
                DType weight=1.0);
        void addCellArrays(DType,DType,DType, 
                DType weight=1.0);

        // Add nodal/cell value constraint
        void addNodalCstr(DType,DType);
        void addCellCstr(DType,DType,DType);

        // Add derivative constraints
        // To be added

        void solve(void);

};

// 2D case

class LocLstSqr2D : public LocLstSqr {
    protected:
        PolyCoef2D* polyCoef;
        VEC<2>      x_c;
        DType       val_c;
        int        *o_x;
        int        *o_y;
        
        void calcOrd(void);
    public:
        void setup(int order_in, PolyCoef2D*, bool val_in);
        void setCenterX(VEC<2> x_c_in) {x_c=x_c_in;}
        void setCenterV(DType);

        // a point in 2D
        void addNodalMat(VEC<2>,
                DType weight=1.0);
        // a rectangle in 2D
        // x_{lower left}, x_{upper right}
        void addCellMat(VEC<2>,VEC<2>,
                DType weight=1.0);
        // a triangle in 2D
        // counter clock wise
        void addCellMat(VEC<2>,VEC<2>,VEC<2>,
                DType weight=1.0);

        // a point in 2D
        void addNodalVal(DType,VEC<2>,
                DType weight=1.0);
        // a rectangle in 2D
        // x_{lower left}, x_{upper right}
        void addCellVal(DType,VEC<2>,VEC<2>,
                DType weight=1.0);
        // a triangle in 2D
        // counter clock wise
        void addCellVal(DType,VEC<2>,VEC<2>,VEC<2>,
                DType weight=1.0);

        // a point in 2D
        void addNodalArrays(DType,VEC<2>,
                DType weight=1.0);
        // a rectangle in 2D
        // x_{lower left}, x_{upper right}
        void addCellArrays(DType,VEC<2>,VEC<2>,
                DType weight=1.0);
        // a triangle in 2D
        // counter clock wise
        void addCellArrays(DType,VEC<2>,VEC<2>,VEC<2>,
                DType weight=1.0);

        // Add nodal/cell value constraint
        void addNodalCstr(DType,VEC<2>);
        void addCellCstr(DType,VEC<2>,VEC<2>);
        void addCellCstr(DType,VEC<2>,VEC<2>,VEC<2>);

        // Add derivative constraints
        // To be added

        void solve(void);
};

#endif
