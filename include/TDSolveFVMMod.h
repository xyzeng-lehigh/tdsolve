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
    void alloc(void);
    bool update_sol(const int fr, const int to, const DType shift_dt_fr=0.0, const DType shift_dt_to=1.0);
    void update_vel(const int fr, const int to, const DType shift_dt=0.0, const bool last_iter=false);

    void compute_rhs(const DType alpha, const int nt);
    void compute_rhs_common(const int nt);
    void update_sol(const int rkns, RKButcher& rkButcher, const int fr, const int to);
    void update_vel(const DType alpha, const int nt);

    void setup_rk_lin( RKButcher& rkButcher );
    bool update_rk_sol( RKButcher& rkButcher, const bool diagnos=false );

    void copy_aux(const int fr, const int to);

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
    void compute_vel_u(const int nt, const DType tsol_loc, const int u_nt=0);
    void compute_vel_u_chem(const int nt, const int u_nt=0);

    void compute_cells_fluxes_u(const int nt, const DType tsol_loc, const int u_nt=0);
    void compute_cells_fluxes_u_upw(const int nt, const DType tsol_loc, const int u_nt=0);
    void compute_cells_fluxes_u_muscl(const int nt, const DType tsol_loc, const int u_nt=0);
    void compute_cells_fluxes_v(const int nt, const int aux_nt=0);
    void compute_cells_fluxes_v_upw(const int nt);
    void compute_cells_fluxes_v_muscl(const int nt, const int slp_nt=0);
    void compute_flux_muscl_lims(const int nt);
    void compute_flux_r_ppm_lims(const int nt);
    void compute_cells_fluxes_r(const int nt, const int rp_nt=0);
    void update_cells_common(const int fr, const int to, bool only_res = false);
    void update_vel_int(const int nt, const int u_nt=0);

    // implicit methods
    void update_cells_orig_imp_upw(const int fr, const int to, const DType tsol_loc);
    void update_cells_orig_imp_muscl(const int fr, const int to, const DType tsol_loc);
    void update_mat_upw_v(const DType multip, const int nt, HVEC**& vel, const int v_nt, const bool rbc=false, const int rp_nt=0);
    void update_mat_orig_muscl_v(const int nt, HVEC**& vel, const int v_nt, HVEC**& lim, const bool rbc=false);
    void update_mat_muscl_v(const DType multip, const int nt, HVEC**& vel, const int v_nt, const bool rbc=false, const int rp_nt=0);
    void update_mat_ppm_r(const DType multip, const int nt, const int rp_nt=0);
    void update_cells_imp_upw(const int fr, const int to, const DType tsol_loc);
    void update_cells_imp_muscl(const int fr, const int to, const DType tsol_loc);

    // Classical FVM
    void compute_cells_fluxes_orig(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_orig_upw(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_orig_upw_auto_g(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_orig_upw_auto_h(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_orig_upw_auto_n(const int nt, const DType tsol_loc);
    void compute_cells_fluxes_orig_muscl(const int nt, const DType tsol_loc);
    void update_vel_int_orig(const int nt);
    void update_vel_int_opn_orig(const int nt);

    void update_cells_common_auto(const int fr, const int to);
    void update_cells_common_auto_g(const int fr, const int to);
    void update_cells_common_auto_h(const int fr, const int to);
    void update_cells_common_auto_n(const int fr, const int to);

    // Update chems
    void update_chem_imp( const int fr, const int to, const DType shift_dt_fr=0.0, const DType shift_dt_to=1.0 );

    void compute_flux_chem_muscl_lims(const int nt, HVEC**& data, HVEC**& lim, const int N_loc, const DType* xc_loc, const DType* xf_loc, const int rbc=0, const DType bc_val=0.0);
    void update_mat_chem_muscl_r(const int nt, LinMat*& mat_loc, HVEC**& lim, const int N_loc, const DType* xc_loc, const DType* xf_loc, const int rbc=0, const int rp_nt=0);
    void update_mat_chem_diff(const int nt, LinMat*& mat_loc, const int N_loc, const DType* xc_loc, const DType* xf_loc, const DType diff_coef, const int rbc = 0);
    void update_chem_muscl_be(const int fr, const int to, const int rp_nt=0);
    void update_chem_opn_muscl_be(const DType tsol_loc, const int fr, const int to, 
                                const DType coef_D, const DType coef_alpha, 
                                const DType coef_beta, const DType coef_gamma, 
                                HVEC**& data_source, HVEC**& data_P, HVEC**& flux_P, HVEC**& slope_P, 
                                const bool treatment=false, const bool treatBC=false, const int rp_nt=0);

    // Right hand sides
    void update_chem_rhs( const int fr, const int to, const DType shift_dtr=0.0 );
    void update_chem_upw_rhs( const int fr, const int to );
    void update_chem_muscl_rhs( const int fr, const int to );

    // OPN model
    void update_vel_int_opn(const int nt, const int u_nt=0);

    // Other functions
    DType calc_n_bc( const DType tsol_loc );

    // The minmod limiter
    DType LimFun_Minmod(DType dval_l, DType dval_r, bool& no_lim, bool& zero_lim); 

    // The PPM limiter
    DType LimFun_PPM_BC(DType dval_c, DType dval_r); 
    void  LimFun_PPM_INT(DType alpha1, DType alpha2, DType& lim_l, DType& lim_r);

    // RK Radau methods
    void update_sol_radau_exp(const int rkns, RKButcher& rkButcher, const int fr, const int to);
    void update_sol_radau_imp_upw(const int rkns, RKButcher& rkButcher, const int fr, const int to);
    void update_sol_radau_imp_muscl(const int rkns, RKButcher& rkButcher, const int fr, const int to);

    void copy_rk_data_to_sol( LinVec*& vec, RKButcher& rkButcher );
    void add_sol_to_rk_data( const DType scale, LinVec*& vec, RKButcher& rkButcher );

    void eval_rk_radau_fval( LinVec*& vec, RKButcher& rkButcher );
    void compute_rhs_primitive(const DType alpha, const int nt);
    void compute_cells_fluxes_u_primitive(const int nt, const DType tsol_loc, const int u_nt=0);
    void compute_cells_fluxes_v_primitive(const int nt, const int slp_nt=0);
    void compute_cells_fluxes_r_primitive(const int nt, const int rp_nt=0);
    void compute_rhs_common_primitive(const int nt);

    void eval_rk_radau_fjac( LinMat*& lmat, RKButcher& rkButcher );

    void update_vel_primitive( RKButcher& rkButcher );
    void update_vel_int_primitive(const int nt, const int u_nt=0);
    void update_vel_int_opn_primitive(const int nt, const int u_nt=0);

    void assemble_rk_mat_dcell( LinMat*& lmat, RKButcher& rkButcher, const int nt );
    void assemble_matblock_cells_g( const int nt, LinMat*& matblk );
    void assemble_matblock_cells_h( const int nt, LinMat*& matblk );
    void assemble_matblock_cells_n( const int nt, LinMat*& matblk );

    void update_matblock_advective( const int offset_r, const int offset_c, LinMat*& matblk, HVEC**& vel, const int nt, const bool rbc=false );
    void update_matblock_geometric( const int offset_r, const int offset_c, LinMat*& matblk, const int nt );

    void assemble_rk_mat_dvel( LinMat*& lmat, RKButcher& rkButcher, const int nt );

    void assemble_matblock_dv( const int nt, LinMat*& matslm, const int offset, HVEC**& data, HVEC**& slope );

    void assemble_matblock_dvel_dg( const int nt, LinMat*& matflt );
    void assemble_matblock_dvel_dh( const int nt, LinMat*& matflt );
    void assemble_matblock_dvel_dn( const int nt, LinMat*& matflt );

};

#endif

