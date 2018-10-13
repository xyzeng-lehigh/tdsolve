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
    virtual ~TDSolve(void);

    virtual void setup(void) = 0;
    virtual void set_ic(void) = 0;
    virtual void calc_dt(const int nt) = 0;
    virtual void output(const int nt) = 0;
    void diagnostic(const int nt);

    // The solver centers at tsol
    // In case of higher-order time integration, one may need to update the solutions from time collocation
    //   tsol + shift_dt_fr * dtn
    // to
    //   tsol + shift_dt_to * dtn
    // Return value is true if the update is accepted
    // o.w., false is returned, and typically dt is reduced
    virtual bool update_sol(const int fr, const int to, const DType shift_dt_fr=0.0, const DType shift_dt_to=1.0) = 0;
    // tsol + shift_dt_to * dtn is the targeted time collocation for the averaged solution
    void average_sol(const DType coef_1, const int nt_1, 
                     const DType coef_2, const int nt_2, const int to, const DType shift_dt_to=1.0);
    virtual void copy_data(const int fr, const int to) = 0;

    bool nextIter(void) {
      iter++;
      return !reach_term_now;
    }

    // offset_dt is the difference between the target tsol and the current tsol
    virtual void update_vel(const int fr, const int to, const DType offset_dt=0.0, const bool last_iter=false) = 0;
    void update_rad(const int fr, const int to, const int v_nt=-1);

    void scale_dtn(const DType scale) { dtn *= scale; }
    void reduce_dtn(const DType scale);

    // Functions called by API
    bool check_Rterm(const int nt);
    DType get_tsol(void) {return tsol;}
    virtual DType get_nperc(const int nt) {return 0.0;};

  protected:
    HVEC** data_g;  // G = number density of glimoa cells
    HVEC** data_h;  // H = number density of necrotic cells
    HVEC** data_n;  // N --> M = number density of immune cells infiltrated into the tumor 
    HVEC** data_a;  // A = concentration of chemotatic attractants by glimoa cells
    HVEC** data_v;  // V = veocity field within the tumor boundary
    DType  data_r[3];
    DType  rprime;  // R'

    // OPN model
    HVEC** data_p1; // concentration of osteopontin secreted by immunue cells
    HVEC** data_p2; // concentration of osteopontin secreted by glimoa cells

    HVEC** data_dv; // (dv/dx)/R+(2v)/(xR)
    HVEC** data_na; // (1/x^2/R^2) d(x^2 N dA/dx )/dx

    TDInput* tdInput;
    TDCoef*  tdCoef;
    TDMesh*  tdMesh;
    TDSetup* tdSetup;
    TDLin*   tdLin;

    int    Nx; // will be used for Np1 and Np2
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
    LinMat*   mat_p;
    LinMat*   prc_p;
    LinVec*   sol_p;
    LinVec*   rhs_p;
    LinVec*   res_p;
    LinSolve* linSolve_p;

    int iter;

    void setupLinSolve( LinMat*&, LinMat*&, LinVec*&, LinVec*&, LinVec*&, LinSolve*&, int );

    void upscale_var(HVEC&, HVEC&, DType*&, int, DType, int scal_pow=2);
    void downscale_var(HVEC&, HVEC&, DType*&, int, DType, int scal_pow=2);

};

#endif
