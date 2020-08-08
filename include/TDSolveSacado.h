#ifndef TDSOLVE_SACADO_H
#define TDSOLVE_SACADO_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"
#include "TDSolve.h"
#include "Sacado.hpp"

// Essentially FVMMod, but with Sacado arrays
class TDSolveSacado : public TDSolve {
  public:
    TDSolveSacado(TDInput* tdInput_ptr);
    ~TDSolveSacado(void);

    void setup(void);
    void set_ic(void);

    void calc_dt(const int nt);
    void output(const int nt);
    void diagnostic(const int nt);

    void alloc(void);
    bool update_sol(const int fr, const int to, const DType shift_dt_fr=0.0, const DType shift_dt_to=1.0);
    void copy_data(const int fr, const int to);
    // Needs to update whenever a or n is modified
    // Before calling, the S_ independent arrays should be associated with fr
    // After calling, (*data_v)[to] will be filled
    void update_vel(const int fr, const int to, const DType shift_dt=0.0, const bool last_iter=false);

    // Functions called by APE
    DType get_gperc(const int nt);
    DType get_nperc(const int nt);

  protected:
    struct ChemPar {
      DType decay;
      DType diffusion;
      // Michaelis Mentun coefficients
      DType mm_ceil, mm_half;
      std::string driver;
      bool dirichlet_bc;
      // Testing options
      int test_opt;
    };
    struct LSPack {
      LinMat*   mat;
      LinMat*   prc;
      LinVec*   sol;
      LinVec*   rhs;
      LinVec*   res;
      LinSolve* linSolve;
    };
    // Variable names
    std::vector<std::string> cellnames;   // All cell species name
    std::vector<std::string> cellnames_u; // Cell species subject to infiltration
    std::vector<std::string> velnames;    // All velocity names
    std::vector<std::string> chemnames;   // All chemoattractant/OPN names
    std::vector<std::string> ale_chemnames;
    std::map<std::string, int> varsizes;
    std::map<std::string, DType*> varcoord;
    std::map<std::string, DType*> varflxcoord; // Only exist for cells/chems
    std::map<std::string, ChemPar> chempars;
    std::map<std::string, LSPack> varlspacks;
    std::map<int, LSPack> dimlspacks; // Avoid duplicating linear solvers of the same size
    LSPack cell_lspack; // A linear solver for all cells
    std::map<std::string, int> H_cell_offset;
    // A map from cell names to non-geometric velocities
    std::map<std::string,std::set<std::string>> cellvels;

    void add_vel_to_cellvels(std::string, std::string);

    // Independent variables (cell numbers and chemicals)
    typedef std::vector<Sacado::Fad::DFad<DType>> SacadoIndVarArray;
    std::map<std::string, SacadoIndVarArray> S_data_var;
    std::map<std::string, int> S_data_offset;
    int S_data_num;
    // HVEC version
    std::map<std::string, HVEC**> H_data_var;
    std::map<std::string, HVEC**> H_data_res;

    // Infiltration velocity
    HVEC** data_u;
    std::map<std::string, HVEC**> H_data_vel;

    // Dependent variables
    typedef Sacado::Fad::DFad<DType> SacadoDepVar;
    typedef std::vector<SacadoDepVar> SacadoDepVarArray;
    typedef std::map<std::string, SacadoDepVarArray> SacadoDepVarMap;
    typedef std::map<std::string, SacadoDepVarMap> SacadoDepVarMapMap;
    // R, R', etc (they may not need to be allocated as Sacado array)
    SacadoDepVar S_data_r;
    SacadoDepVar S_data_rp;
    // Velocities
    SacadoDepVarMap S_data_vel;
    // Cell numbers
    HVEC** data_inc; // Incompressibility violation
    SacadoDepVarMap    S_flux_cell;
    SacadoDepVarMap    S_flux_cell_r;
    SacadoDepVarMapMap S_flux_cell_v;
    SacadoDepVarMap    S_srce_cell;
    SacadoDepVarMap    S_rec_ppm_cell;
    SacadoDepVarMap    S_lim_ppm_cell[2];
    SacadoDepVarArray  S_lim_ppm_sync[2];
    //std::vector<Sacado::Fad::DFad<DType>> S_lim_cell_sync; // MUSCL
    //SacadoDepVarMap S_slope_cell; // MUSCL

    // Chemoattractant/OPN 
    SacadoDepVarMap S_srce_chem;
    SacadoDepVarMap S_flux_chem;
    SacadoDepVarMap S_flux_chem_r;
    SacadoDepVarMap S_flux_chem_d;
    //std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_slope_chem; // MUSCL or PPM
    //std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_slope_p;
    
    // Source term function handles
    typedef void (TDSolveSacado::*fCalcSrce)(SacadoDepVarArray&,const DType);
    std::map<std::string, fCalcSrce> F_calcSrceMap;

    // Additional setup functions
    void setup_chem_par(std::string chemname, const DType decay, const DType diffusion,
        std::string cellname, const DType mm_ceil, const DType mm_half, const bool dbc, const int test_opt=0);

    // Initial condition for three cell numbers
    void setic_unif_cell(std::string name, DType val, int nt=0);
    void setic_sine_cell(std::string name, DType c0, DType c1, DType cs, int nt=0);
    void setic_step_cell(std::string name, DType x_sep, DType val_l, DType val_r, int nt=0);
    void setic_expd_chem(DType multip, DType decay, int nt=0);
    void setic_polyexp_chem(DType poly_0, DType poly_2, int nt=0);

    // Alloc/deallocate additional arrays
    void alloc_hvec_res(void);

    // Alloc Sacado arrays
    void alloc_sacado_ind(void); // Independent variables
    void alloc_sacado_ind_array(std::string, const int);
    void alloc_sacado_dep(void); // Dependent variables
    void alloc_sacado_dep_array(SacadoDepVarArray&, const int);

    // Alloc functions
    void alloc_functions(void);

    // Other Sacado utilities
    void set_sacado_ind(const int);
    void set_sacado_ind(std::vector<std::string>& names, const int);
    void set_sacado_ind_array(std::string, const int);
    void set_sacado_r(const int);
    void set_sacado_v(const int);
    void set_from_sacado_ind(const int);
    void set_from_sacado_ind_array(std::string, HVEC&);
    void set_from_sacado_r(const int);
    void set_from_sacado_vel(std::string, const int);
    void reset_sacado_ind_array(void);
    void reset_sacado_ind_array(SacadoIndVarArray& s_var, const int s_offset, DType val=0.0);
    void reset_sacado_dep_array(void);
    void reset_sacado_dep_array(SacadoDepVarMapMap&, DType val=0.0);
    void reset_sacado_dep_array(SacadoDepVarMap&, DType val=0.0);
    void reset_sacado_dep_array(SacadoDepVarArray&, DType val=0.0);
    void reset_sacado_dep_array(SacadoDepVarArray&, DType val, int num);
    void sum_in_sacado_dep_array(SacadoDepVarArray&, DType, SacadoDepVarArray&);
    void print_sacado_ind_array(SacadoIndVarArray&);
    void print_sacado_dep_array(SacadoDepVarArray&);

  private:
    void setup_new_lin_solve(const int num);

    void calc_total_inc_viol(const int nt);
    void compute_vel_u(const DType tsol_loc);
    void compute_vel_u_pres(const DType tsol_loc);
    void compute_vel_u_chem(void);

    void compute_vel_v(void);

    void compute_cells_varflux(std::string velname, std::string cellname, const DType tsol_loc, const bool rbc=false);
    void compute_cells_varflux_upw(SacadoDepVarArray& s_vel, SacadoIndVarArray& s_var, SacadoDepVarArray& s_flux, const DType* rdata=NULL);

    void compute_cells_vflux(void);

    void LimFun_PPM_BC(DType dval_c, DType dval_r, SacadoDepVar& s_lim); 
    void LimFun_PPM_BC(SacadoDepVar& s_dval_c, SacadoDepVar& s_dval_r, SacadoDepVar& s_lim); 
    void LimFun_PPM_INT(DType alpha1, DType alpha2, SacadoDepVar& s_lim_l, SacadoDepVar& s_lim_r);
    void LimFun_PPM_INT(SacadoDepVar& s_alpha1, SacadoDepVar& s_alpha2, SacadoDepVar& s_lim_l, SacadoDepVar& s_lim_r);

    void compute_cells_rflux(void);
    void calc_rflux_cell_var_upw(std::string cellname);
    void calc_rflux_cell_var_wupw(std::string cellname);
    inline DType calc_wupw( const DType l_x, const DType l_dx );
    void compute_cells_rflux_ppm(void);
    void calc_ppm_lim_var(SacadoIndVarArray& s_var, SacadoDepVarArray& s_rec, SacadoDepVarArray& s_lim_l, SacadoDepVarArray& s_lim_r, DType* l_x_c, const int nt=-1, const bool tosync=false, std::vector<bool>* set_zero=NULL, SacadoDepVarArray* s_alpha=NULL);
    void calc_ppm_lim_sync(SacadoDepVarArray* s_sync, SacadoDepVarMap* s_lims, std::vector<bool>& set_zero, SacadoDepVarArray* s_alpha);
    void calc_ppm_flux_var(SacadoIndVarArray& s_var, SacadoDepVarArray& s_rec, SacadoDepVarArray& s_flx, SacadoDepVarArray* s_lim, DType* l_x_c, const int nt=-1);

    void compute_cells_flux(void);

    void compute_chems_dflux(void);
    void calc_dflux_chem_var(std::string chemname);

    void compute_chems_rflux(const int nt);
    void calc_rflux_chem_var_upw(std::string chemname, const int nt);
    void calc_rflux_chem_var_ppm(std::string chemname, const int nt);
    void calc_rflux_chem_var_wupw(std::string chemname, const int nt);

    void compute_chems_flux(void);

    void reset_residual(void);
    void compute_cells_residual(const int fr, const int to);
    void assemble_inertia_var(HVEC& h_res, DType r_fr, HVEC& h_var, DType r_to, SacadoIndVarArray s_var, DType* l_x_c);
    void assemble_flux_var(HVEC& h_res, SacadoDepVarArray& s_flx);
    void assemble_srce_var(HVEC& h_res, SacadoDepVarArray& s_src);
    void compute_chems_residual(const int fr, const int to);

    void update_cells_by_res(const int nt);
    void update_sol_by_res_var_exp(HVEC& h_var, HVEC& h_res, DType h_r, DType* l_x_c, const int num);

    void update_cells_by_res_imp(const int nt);
    void update_cell_linsys_by_res_var(HVEC& h_res, DType h_r, DType* l_x_c, SacadoDepVarArray& s_flx, SacadoDepVarArray& s_src, const int num, const int s_offset, const int l_offset, LSPack& lsPack);
    void update_hsol_by_linsol(HVEC& h_var, const int l_offset, const int num, LSPack& lsPack);

    void update_chems_by_res_imp(const int nt);
    void update_chem_by_res_var_imp(HVEC& h_var, HVEC& h_res, DType h_r, DType* l_x_c, SacadoDepVarArray& s_flx, SacadoDepVarArray& s_src, const int num, const int offset, LSPack& lsPack);

    void compute_rad(const int fr, const int to);

    void update_source_geom(int nt);
    void calc_srce_geom(SacadoIndVarArray& s_var, SacadoDepVarArray& s_srce, DType* l_x_c);
    void calc_srce_geom(SacadoIndVarArray& s_var, SacadoDepVarArray& s_srce, DType* l_x_c, const int nt);

    void compute_source_cell_prob(const DType tsol_loc);
    void calc_srce_cell_base(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_g_base(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_base(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_base(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_g_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_g_opn(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_opn(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_opn(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_g_model(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_model(SacadoDepVarArray& s_srce, const DType tsol_loc);

    inline SacadoDepVar icalc_scal_logist_mm(DType,Sacado::Fad::DFad<DType>,DType,DType,DType,Sacado::Fad::DFad<DType>);

    void compute_source_chem_prob(const int nt, const DType tsol_loc);
    void calc_srce_chem_var(std::string chemname, const int nt, const DType tsol_loc);

    // Other functions
    DType calc_n_bc( const DType tsol_loc );
    DType get_cell_perc(std::string cellname, const int nt);

    HVEC** srce_max; // Incompressibility violation
    void reset_srce_max(void);
    void adjust_srce_max(const int i, DType lip_loc); 
    void adjust_srce_max_geom(void);

  private: // Testing: positivity-preserving
    void adjust_pos_ppm_lim( SacadoDepVarArray* s_alpha );

  private: // Testing: fixed range for A
    DType* x_ale_c[2];
    DType* x_ale_f[2];
    DType* dx_ale[2];
    void update_ale_grids(const int nt);
    // Solve (g^num-1)/(g-1) = val
    bool calc_ale_growth(DType& g, const int num, const DType val);
    void print_data_array(std::string,DType*,const int);
    bool check_chem_skip(std::string);

    void compute_ale_source_chem_prob(const int nt, const DType tsol_loc);
    void calc_ale_srce_chem_var(std::string chemname, const int nt, const DType tsol_loc);

    void compute_ale_chems_dflux(const int nt);
    void calc_ale_dflux_chem_var(std::string chemname,const int nt);

    void compute_ale_chems_rflux(const int fr, const int to);
    void calc_ale_rflux_chem_var(std::string chemname, const int fr, const int to);

    void compute_ale_chems_residual(const int fr, const int to);
    void assemble_ale_residual_var(HVEC& h_res, const int fr, const int to, HVEC& h_var, SacadoIndVarArray& s_var, SacadoDepVarArray& s_flx_r, SacadoDepVarArray& s_flx_d, SacadoDepVarArray& s_src);

    void update_ale_chems_by_res_imp(const int nt);
    void update_ale_chem_by_res_var_imp(HVEC& h_var, HVEC& h_res, SacadoDepVarArray& s_flx_r, SacadoDepVarArray& s_flx_d, SacadoDepVarArray& s_src, const int offset, LSPack& lsPack, const int nt);
};

inline DType TDSolveSacado::calc_wupw( const DType l_x, const DType l_dx ) {
  return l_x * ( l_x * l_x - 0.25 * l_dx * l_dx );
}

inline TDSolveSacado::SacadoDepVar TDSolveSacado::icalc_scal_logist_mm(DType scal, Sacado::Fad::DFad<DType> s_logist, DType mm_const, DType mm_ceil, DType mm_half, Sacado::Fad::DFad<DType> s_mm) {
  SacadoDepVar s_val;
  if ( s_logist >= 1.0 )
    s_val = 0.0;
  else if ( s_logist <= 0.0 )
    s_val = scal * (mm_const + mm_ceil * s_mm/(mm_half+s_mm));
  else
    s_val = scal * (1.0-s_logist) * (mm_const + mm_ceil * s_mm/(mm_half+s_mm));
  return s_val;
}

#endif
