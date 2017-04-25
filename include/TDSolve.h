#ifndef TDSOLVE_H
#define TDSOLVE_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"

class  HVEC;
struct TDCoef;
struct TDMesh;
struct TDSetup;
struct TDLin;
class TDInput;
class LinMat;
class LinVec;
class LinSolve;

class TDSolve {
  public:
    TDSolve(TDInput* tdInput_ptr);
    ~TDSolve(void);

    virtual void setup(void) = 0;
    virtual void set_ic(void) = 0;
    virtual void calc_dt(const int nt) = 0;
    virtual void output(const int nt) = 0;
    void diagnostic(const int nt);

    virtual void update_sol(const int fr, const int to) = 0;
    void average_sol(const DType, const int, 
                     const DType, const int, const int);
    virtual void copy_data(const int fr, const int to) = 0;

    bool nextIter(void) {
      iter++;
      return !reach_term_now;
    }
  protected:
    HVEC** data_g;
    HVEC** data_h;
    HVEC** data_n;
    HVEC** data_a;
    HVEC** data_v;
    DType  data_r[3];

    HVEC** data_dv; // (dv/dx)/R+(2v)/(xR)
    HVEC** data_na; // (1/x^2/R^2) d(x^2 N dA/dx )/dx

    TDInput* tdInput;
    TDCoef*  tdCoef;
    TDMesh*  tdMesh;
    TDSetup* tdSetup;
    TDLin*   tdLin;

    int    Nx;
    int    Na_m;
    int    Na;
    DType  dx;
    DType* x_f;
    DType* x_c;
    DType* x_a_f;
    DType* x_a_c;

    DType  dtn, tsol;
    bool   reach_term_next;
    bool   reach_term_now;

    // Linear solver for anything else other than A
    LinMat*   mat;
    LinMat*   prc;
    LinVec*   sol;
    LinVec*   rhs;
    LinVec*   res;
    LinSolve* linSolve;

    // Linear solver for A
    LinMat*   mat_a;
    LinMat*   prc_a;
    LinVec*   sol_a;
    LinVec*   rhs_a;
    LinVec*   res_a;
    LinSolve* linSolve_a;

    int iter;

    void setupLinSolve( LinMat*&, LinMat*&, LinVec*&, LinVec*&, LinVec*&, LinSolve*&, int );
    virtual void update_vel(const int nt) = 0;
    void update_rad(const int fr, const int to);
};

#endif
