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

    void update_vel(const int nt);

    void update_cells_upw(const int fr, const int to);
    void update_cells_muscl(const int fr, const int to);

  private:
    // Original model
    void compute_cells_fluxes_u(const int nt);
    void compute_cells_fluxes_u_upw(const int nt);
    void compute_cells_fluxes_u_muscl(const int nt);
    void compute_cells_fluxes_v(const int nt);
    void compute_cells_fluxes_v_upw(const int nt);
    void compute_cells_fluxes_v_muscl(const int nt);
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

};

#endif

