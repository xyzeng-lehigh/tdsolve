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
    HVEC** data_g;// G = number density of glimoa cells
    HVEC** data_h;// H = number density of necrotic cells
    HVEC** data_n;// N --> M = number density of immune cells infiltrated into the tumor 
    HVEC** data_a;// A = concentration of chemotatic attractants by glimoa cells
    HVEC** data_v;// V = veocity field within the tumor boundary
    HVEC** data_p1;//concentration of osteopontin secreted by immunue cells
    HVEC** data_p2;// concentration of osteopontin secretedd by glimoa cells
    DType  data_r[3];

    HVEC** data_dv; // (dv/dx)/R+(2v)/(xR)
    HVEC** data_na; // (1/x^2/R^2) d(x^2 N dA/dx )/dx

    TDInput* tdInput;
    TDCoef*  tdCoef;
    TDMesh*  tdMesh;
    TDSetup* tdSetup;
    TDLin*   tdLin;

    int    Nx;// will be used for Np1 and Np2
    int    Na_m;
    int    Na;
    

    DType  dx;
    DType* x_f;
    DType* x_c;
    DType* x_a_f;
    DType* x_a_c;

    DType  dtn, tsol, t_output, dtn_old;
    bool   reach_term_next;
    bool   reach_output_now;
    bool   reach_output_next;
    bool   reach_term_now;
    bool   use_tfreq;

    // Linear solver for anything else other than A, P1, P2
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

    // Linear solver for P1 OPN model                                                                                                                        
    LinMat*   mat_p1;
    LinMat*   prc_p1;
    LinVec*   sol_p1;
    LinVec*   rhs_p1;
    LinVec*   res_p1;
    LinSolve* linSolve_p1;

   // Linear solver for P2 OPN model

    LinMat*   mat_p2;
    LinMat*   prc_p2;
    LinVec*   sol_p2;
    LinVec*   rhs_p2;
    LinVec*   res_p2;
    LinSolve* linSolve_p2;


    int iter;

    void setupLinSolve( LinMat*&, LinMat*&, LinVec*&, LinVec*&, LinVec*&, LinSolve*&, int );
    virtual void update_vel(const int nt) = 0;
    void update_rad(const int fr, const int to);
};

#endif
