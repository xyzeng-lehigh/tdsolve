#ifndef TDSOLVE_H
#define TDSOLVE_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"
#include "RK_Coef.h"

class  HVEC;
class  HMAT;
struct TDCoef;
struct TDMesh;
struct TDSetup;
struct TDLin;
class TDInput;
class LinMat;
class LinVec;
class LinSolve;
struct RKButcher;

class TDSolve {
  public:
    TDSolve(TDInput* tdInput_ptr);
    virtual ~TDSolve(void);

    virtual void setup(void) = 0;
    virtual void set_ic(void) = 0;
    virtual void calc_dt(const int nt) = 0;
    virtual void output(const int nt) = 0;
    void diagnostic(const int nt);
    virtual void alloc(void) = 0;

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
    virtual void copy_aux(const int fr, const int to) {};

    virtual void compute_rhs(const DType alpha, const int nt) {};
    virtual void update_sol(const int rkns, RKButcher& rkButcher, const int fr, const int to) {};
    virtual void update_vel(const DType alpha, const int nt) {};

    virtual void setup_rk_lin( RKButcher& rkButcher ) {};
    virtual void update_rk_sol( RKButcher& rkButcher, const bool diagnos=false ) {};

    bool nextIter(void) {
      iter++;
      return !reach_term_now;
    }

    // offset_dt is the difference between the target tsol and the current tsol
    virtual void update_vel(const int fr, const int to, const DType offset_dt=0.0, const bool last_iter=false) = 0;
    void update_rad(const int fr, const int to, const int v_nt=-1);
    void update_rad_radau(RKButcher&, const bool diagnos=false);
    void eval_rad_radau_fval(HVEC&, const RKButcher&, const DType*, const HVEC&);
    void eval_rad_radau_hess(HMAT&, const RKButcher&, const DType*, const HVEC&);
    void check_conv_radau(const RKButcher&, int&, bool&, bool&, const bool diagnos=false);
    void check_conv_radius(const RKButcher&, int&, bool&, bool&, const bool diagnos=false);

    void scale_dtn(const DType scale) { dtn *= scale; }
    void reduce_dtn(const DType scale);

    // Functions called by API
    bool check_Rterm(const int nt);
    DType get_tsol(void) {return tsol;}
    virtual DType get_nperc(const int nt) {return 0.0;};

    void downscale_sol(const int s_nt, const int r_nt);
    void upscale_sol(const int s_nt, const int r_nt);
    void print_warning(const char*); 

    TDInput*& getTDInput(void) {
      return tdInput;
    }

  protected:
    HVEC** data_g;  // G = number density of glimoa cells
    HVEC** data_h;  // H = number density of necrotic cells
    HVEC** data_n;  // N --> M = number density of immune cells infiltrated into the tumor 
    HVEC** data_a;  // A = concentration of chemotatic attractants by glimoa cells
    HVEC** data_v;  // V = veocity field within the tumor boundary
    DType* data_r;  // R
    DType* data_rp; // R'

    HVEC** data_cell; // Temporary array for cellular values
    HVEC** data_node; // Temporary array for nodal values

    DType* diff_g;  // Difference in G
    DType* diff_h;  // Difference in H
    DType* diff_n;  // Difference in N
    DType* diff_v;  // Difference in V
    DType* diff_r;  // Difference in R

    // OPN model
    HVEC** data_p1; // concentration of osteopontin secreted by immunue cells
    HVEC** data_p2; // concentration of osteopontin secreted by glimoa cells

    HVEC** data_dv; // (dv/dx)/R+(2v)/(xR)
    HVEC** data_na; // (1/x^2/R^2) d(x^2 N dA/dx )/dx

    bool   use_resd;

    HVEC** resd_g;
    HVEC** resd_h; 
    HVEC** resd_n;
    HVEC** resd_a;

    // OPN model
    HVEC** resd_p1;
    HVEC** resd_p2;

    TDInput* tdInput;
    TDCoef*  tdCoef;
    TDMesh*  tdMesh;
    TDSetup* tdSetup;
    TDLin*   tdLin;

    int    Ndata; // Number of copies of data

    int    Nx; // will be used for Np1 and Np2
    int    Na_m;
    int    Na;
    DType  dx;
    DType* x_f;
    DType* x_c;
    DType* x_a_f;
    DType* x_a_c;

    DType  dtn, tsol, t_output, dtn_old;
    DType  cur_cfl;
    bool   reach_term_next;
    bool   reach_output_now;
    bool   reach_output_next;
    bool   reach_term_now;
    bool   use_tfreq;

    // Used by check_conv_radau
    DType init_err;

    // Linear solver for anything else other than A, P1, P2
    LinMat*   mat;
    LinMat*   prc;
    LinVec*   sol;
    LinVec*   rhs;
    LinVec*   res;
    LinSolve* linSolve;

    // Linear solver for Radau RK for cell species
    LinMat*   matBlk;
    LinMat*   matFlt;
    LinMat*   matSlm;
    LinMat*   matRK;
    LinMat*   prcRK;
    LinVec*   solRK;
    LinVec*   rhsRK;
    LinVec*   resRK;
    LinSolve* linSolveRK;

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
