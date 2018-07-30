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
    void calc_dt(const int nt);
    bool update_sol(const int fr, const int to, const DType shift_dt_fr=0.0, const DType shift_dt_to=1.0);
    void update_vel(const int fr, const int to, const DType shift_dt=0.0, const bool last_iter=false);


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
    HVEC** lim_muscl; // Sync limiter for MUSCL

    HVEC** flux_g_v;
    HVEC** flux_g_r;
    HVEC** flux_h_v;
    HVEC** flux_h_r;
    HVEC** flux_n_v; // or M
    HVEC** flux_n_r; // or M
    HVEC** flux_n_u; // or M

    HVEC** data_u;   // U = velocity due to chemical gradient
    HVEC** vel_tmp;  // Temporary storage for velocities

    void calc_total_diff_model(const int nt);
    void print_custom_vel_u( FILE*, const int nt );

  private:
    // Original model
    void compute_vel_u(const int nt, const DType tsol_loc);
    void compute_vel_u_chem(const int nt);

    void compute_cells_fluxes_u(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_u_upw(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_u_muscl(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_v(const int nt);
    void compute_cells_fluxes_v_upw(const int nt);
    void compute_cells_fluxes_v_muscl(const int nt);
    void compute_flux_muscl_lims(const int nt);
    void compute_flux_r_ppm_lims(const int nt);
    void compute_cells_fluxes_r(const int nt);
    void update_cells_common(const int fr, const int to);
    void update_vel_int(const int nt);

    // implicit methods
    void update_cells_orig_imp_upw(const int fr, const int to, const DType tsol_loc);
    void update_cells_orig_imp_muscl(const int fr, const int to, const DType tsol_loc);
    void update_mat_upw_v(const int nt, HVEC**& vel, const int v_nt, const bool rbc=false);
    void update_mat_orig_muscl_v(const int nt, HVEC**& vel, const int v_nt, HVEC**& lim, const bool rbc=false);
    void update_mat_muscl_v(const int nt, HVEC**& vel, const int v_nt, const bool rbc=false);
    void update_mat_ppm_r(const int nt);
    void update_cells_imp_upw(const int fr, const int to, const DType tsol_loc);
    void update_cells_imp_muscl(const int fr, const int to, const DType tsol_loc);

    // Classical FVM
    void compute_cells_fluxes_orig(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_orig_upw(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_orig_muscl(const int nt, const DType tsol_loc);
    void update_vel_int_orig(const int nt);

    // Update chems
    void update_chem_imp( const int fr, const int to, const DType shift_dt_fr=0.0, const DType shift_dt_to=1.0 );

    void compute_flux_chem_muscl_lims(const int nt, HVEC**& data, HVEC**& lim, const int N_loc, const DType* xc_loc, const DType* xf_loc, const int rbc=0, const DType bc_val=0.0);
    void update_mat_chem_muscl_r(const int nt, LinMat*& mat_loc, HVEC**& lim, const int N_loc, const DType* xc_loc, const DType* xf_loc, const int rbc=0);
    void update_mat_chem_diff(const int nt, LinMat*& mat_loc, const int N_loc, const DType* xc_loc, const DType* xf_loc, const DType diff_coef, const int rbc = 0);
    void update_chem_muscl_be(const int fr, const int to);
    void update_chem_opn_muscl_be(const DType tsol_loc, const int fr, const int to, 
                                const DType coef_D, const DType coef_alpha, 
                                const DType coef_beta, const DType coef_gamma, 
                                HVEC**& data_source, HVEC**& data_P, HVEC**& flux_P, HVEC**& slope_P, 
                                const bool treatment=false, const bool treatBC=false);

    // OPN model
    void update_vel_int_opn(const int nt);

    // Other functions
    DType calc_n_bc( const DType tsol_loc );

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

