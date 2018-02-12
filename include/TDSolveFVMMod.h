#ifndef TDSOLVE_FVM_MOD_H
#define TDSOLVE_FVM_MOD_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"
#include "TDSolve.h"
#include "TDSolveFVM.h"

class TDSolveFVMMod : public TDSolveFVM {
  public:
    void set_ic(void);
    void update_sol(const int fr, const int to);

    TDSolveFVMMod(TDInput* tdInput_ptr);
    ~TDSolveFVMMod(void);

  protected:
    // lim_g, lim_h, lim_n are not actually used,
    // but they may be implemented in the slopes for
    // comparison
    HVEC** lim_g;
    HVEC** lim_h;
    HVEC** lim_n; // or M
    HVEC** lim_s; // Sync

    HVEC** flux_g_v;
    HVEC** flux_g_r;
    HVEC** flux_h_v;
    HVEC** flux_h_r;
    HVEC** flux_n_v; // or M
    HVEC** flux_n_r; // or M
    HVEC** flux_n_u; // or M

    HVEC** data_u;   // U = velocity due to chemical gradient

    void update_vel(const int nt);
    void print_custom_vel_u( FILE*, const int nt );

  private:
    // Original model
    void compute_vel_u(const int nt);
    void compute_vel_u_chem(const int nt);

    void compute_cells_fluxes_u(const int nt);
    void compute_cells_fluxes_u_upw(const int nt);
    void compute_cells_fluxes_u_muscl(const int nt);
    void compute_cells_fluxes_v(const int nt);
    void compute_cells_fluxes_v_upw(const int nt);
    void compute_cells_fluxes_v_muscl(const int nt);
    void compute_cells_fluxes_r(const int nt);
    void update_cells_common(const int fr, const int to);
    void update_vel_int(const int nt);

    // OPN model
    void update_vel_int_opn(const int nt);

    // The minmod limiter
    DType LimFun_Minmod(DType dval_l, DType dval_r) {
      if ( dval_l * dval_r <= 0.0 )
        return 0.0;
      else
        return min(dval_l/dval_r, 1.0);
    }

    // The PPM limiter
    DType LimFun_PPM_BC(DType dval_c, DType dval_r); 
    void  LimFun_PPM_INT(DType alpha1, DType alpha2, DType& lim_l, DType& lim_r);

};

#endif

