#ifndef TDSOLVE_SACADO_H
#define TDSOLVE_SACADO_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <string>
#include "Data_Type.h"
#include "HighDim_Alg.h"
#include "TDSolve.h"

class LinBandSys;

// Automatic differentiation: the built-in AutoDiff::Fad::DFad by default.
// Configuring with -DUSE_TRILINOS=ON switches back to the Sacado package of
// Trilinos, which provides the same class as Sacado::Fad::DFad; this is only
// meant for cross-checking the two.
#ifdef TDS_USE_TRILINOS
#include "Sacado.hpp"
namespace AutoDiff = Sacado;
#else
#include "AutoDiff.h"
#endif

// Essentially FVMMod, but with Sacado arrays
class TDSolveSacado : public TDSolve {
  public:
    TDSolveSacado(TDInput* tdInput_ptr);
    ~TDSolveSacado(void);

    std::string prob_name;

    typedef bool (TDSolveSacado::*fSetupFun)(void);
    std::map<int, fSetupFun> setupVarMap;
    std::map<int, fSetupFun> setupICMap;
    std::map<int, fSetupFun> setupSrceMap;
    std::map<int, fSetupFun> allocVelUVarMap;
    std::map<int, fSetupFun> computeVelUFlxMap;

    typedef bool (TDSolveSacado::*fSetupTFun)(const DType);
    std::map<int, fSetupTFun> setupVelUMap;
    std::map<int, fSetupTFun> setupExtraVFlxMap;
    std::map<int, fSetupTFun> setupAuxMap;

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

    // Functions called by API
    DType get_gperc(const int nt);
    DType get_nperc(const int nt);
    DType get_argmax(const std::string name, const int nt);
    DType get_total(const std::string name, const int nt);

  protected:
    bool enable_time_delay;
    DType time_delay_max;
    // Solutions at previous times: prev_cell[name][n] is the variable name
    // at the time prev_time[n]
    std::vector<DType>                       prev_time; // Increasing order
    std::map<std::string,std::vector<HVEC*>> prev_cell;
    // Linear interpolation in time between two consecutive stored solutions
    struct PrevInterp {
      int   i_left, i_rght;
      DType coef_left, coef_rght;
    };
    void clear_prev_cell(int i);
    void push_to_prev_cells(DType cur_time, int nt);
    bool get_prev_interp(DType pst_time, PrevInterp& interp);
    bool get_prev_cell(DType pst_time, std::map<std::string,HVEC*>& pst_cell);
    bool get_prev_cell(DType pst_time, int i, std::map<std::string,DType>& pst_data);
    // Variable stored in data (= prev_cell[name]) in the i-th interval
    inline DType calc_prev_data(const std::vector<HVEC*>& data, const int i, const PrevInterp& interp) {
      return interp.coef_left * (*data[interp.i_left])[i] + interp.coef_rght * (*data[interp.i_rght])[i];
    }
    
    // Specific toggle options for oncolyticOld model
    bool in_oncolyticOld_stage[3];

  protected:
    // Setting up variables
    void add_cellular_var(std::vector<std::string>& names, std::string name, const int size, DType* x_var, DType* x_flx);
    void add_nodal_var(std::vector<std::string>& names, std::string name, const int size, DType* x_var);

    bool setup_var_model(void);
    bool setup_var_pdgf(void);
    bool setup_var_opn(void);
    bool setup_var_oncolyticOld(void);
    bool setup_var_potassium(void);
    bool setup_var_oncolytic(void);
    bool setup_var_onc_model(void);

    bool setup_ic_model(void);
    bool setup_ic_pdgf(void);
    bool setup_ic_opn(void);
    bool setup_ic_oncolyticOld(void);
    bool setup_ic_potassium(void);
    bool setup_ic_oncolytic(void);
    bool setup_ic_onc_model(void);

    bool setup_srce_model(void);
    bool setup_srce_pdgf(void);
    bool setup_srce_opn(void);
    bool setup_srce_oncolyticOld(void);
    bool setup_srce_potassium(void);
    bool setup_srce_oncolytic(void);
    bool setup_srce_onc_model(void);

    bool alloc_uvel_map_all(void);
    bool alloc_uvel_map_oncolytic(void);

    bool add_vupdate_uflx(void);
    bool add_vupdate_uflx_oncolytic(void);
    bool add_vupdate_uflx_gen(std::string velname, std::string cellname);

    bool compute_extra_vflx_infiltration(const DType);
    bool compute_extra_vflx_oncolytic(const DType);

    bool setup_aux_oncolyticOld(const DType);
    bool setup_aux_oncolytic(const DType);

  protected:
    // Dependent variables
    typedef AutoDiff::Fad::DFad<DType> SacadoDepVar;
    typedef std::vector<SacadoDepVar> SacadoDepVarArray;
    typedef std::map<std::string, SacadoDepVarArray> SacadoDepVarMap;
    typedef std::map<std::string, SacadoDepVarMap> SacadoDepVarMapMap;

    // Source term function handles
    typedef void (TDSolveSacado::*fCalcSrce)(SacadoDepVarArray&,const DType);
    typedef void (TDSolveSacado::*fCalcTSrce)(std::string,const int,const DType);
    std::map<std::string, fCalcSrce> F_calcSrceMap;

    // Inhomogeneous bc function handle
    typedef DType (TDSolveSacado::*fCalcBC)(const DType);
    struct ChemPar {
      DType decay;
      DType diffusion;
      fCalcTSrce srceFun;
      // Michaelis Mentun coefficients
      DType mm_ceil, mm_half;
      std::string driver;
      bool dirichlet_bc;
      fCalcBC F_chemBC;
      // Testing options
      int test_opt;
    };
    // A linear system is always assembled in sys (see Alg_BandSolve.h). It
    // is solved there if a direct solver is requested (lsysType = DIR with
    // lsysDirect = LU and without preconditioner); otherwise it is copied
    // to mat and rhs, solved by linSolve, and the solution is copied back.
    // mat, ..., linSolve are shared by all variables of the same size, while
    // each use has its own sys so that sparsity patterns do not mix.
    struct LSPack {
      LinMat*   mat;
      LinMat*   prc;
      LinVec*   sol;
      LinVec*   rhs;
      LinVec*   res;
      LinSolve* linSolve;
      LinBandSys* sys;
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
    typedef std::vector<AutoDiff::Fad::DFad<DType>> SacadoIndVarArray;
    std::map<std::string, SacadoIndVarArray> S_data_var;
    std::map<std::string, int> S_data_offset;
    int S_data_num;
    // Derivative components of the independent variables
    //   S_seed_wid = 0: The i-th unknown of a variable is the independent
    //     variable number S_data_offset[name]+i, and S_data_num is the total
    //     number of unknowns. A derivative costs S_data_num operations.
    //   S_seed_wid > 0: The unknowns i, i+S_seed_wid, i+2*S_seed_wid, ...
    //     of a variable share the derivative component
    //     S_data_offset[name] + i % S_seed_wid, so that S_data_num is only
    //     S_seed_wid times the number of variables, whatever Nx is. Then
    //     f.dx( seed_idx(S_data_offset[name],j) ) is the derivative of f
    //     with respect to the j-th unknown of the variable name PROVIDED
    //     THAT f depends on no other unknown j+m*S_seed_wid of that variable.
    //     This holds for all derivatives used for the Jacobians if
    //       1. source terms in an interval only depend on the unknowns in
    //          the same interval, and
    //       2. the fluxes of the interval i only depend on the unknowns in
    //          the intervals i-sten to i+sten, with 2*sten+1 <= S_seed_wid.
    //     A setup_var_* function may set S_seed_wid if its problem meets
    //     these conditions (see setup_var_potassium); the input toggle
    //     Jacobian_compress = 0 keeps it zero.
    int S_seed_wid;
    // Number of derivative components used by a variable with num unknowns
    inline int num_seed(const int num) const {
      return ( S_seed_wid > 0 && num > S_seed_wid ) ? S_seed_wid : num;
    }
    // Derivative component of the i-th unknown of a variable
    inline int seed_idx(const int s_offset, const int i) const {
      return s_offset + ( S_seed_wid > 0 ? i % S_seed_wid : i );
    }
    // HVEC version
    std::map<std::string, HVEC**> H_data_var;
    std::map<std::string, HVEC**> H_data_res;

    // Infiltration velocity
    HVEC** data_u;
    std::map<std::string, HVEC**> H_data_vel;

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
    //std::vector<AutoDiff::Fad::DFad<DType>> S_lim_cell_sync; // MUSCL
    //SacadoDepVarMap S_slope_cell; // MUSCL

    // Chemoattractant/OPN 
    SacadoDepVarMap S_srce_chem;
    SacadoDepVarMap S_flux_chem;
    SacadoDepVarMap S_flux_chem_r;
    SacadoDepVarMap S_flux_chem_d;
    //std::map<std::string, std::vector<AutoDiff::Fad::DFad<DType>>> S_slope_chem; // MUSCL or PPM
    //std::map<std::string, std::vector<AutoDiff::Fad::DFad<DType>>> S_slope_p;
    
    // Additional setup functions
    void setup_chem_par(std::string chemname, const DType decay, const DType diffusion,
        std::string cellname, const DType mm_ceil, const DType mm_half, const bool dbc, const int test_opt=0);
    bool setup_chem_par_gen(std::string chemname, const DType decay, 
        const DType diffusion, const bool dbc, fCalcTSrce srceFun );

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
    void set_from_sacado_vel(const int);
    void set_from_sacado_vel_ind(std::string, const int);
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
    void solve_linsys(LSPack& lsPack);

    void calc_total_inc_viol(const int nt);
    void compute_vel_u(const DType tsol_loc);
    bool compute_vel_u_pres(const DType tsol_loc);
    bool compute_vel_u_pdgf(const DType tsol_loc);
    bool compute_vel_u_oncolytic(const DType tsol_loc); 

    // mode = 'Y': Stencil of driving variable is large enough
    //        'D': Dirichlet data
    //        'N': Neumann data
    void compute_vel_grad(const char mode, std::string velname, std::vector<DType> coef, std::vector<std::string> trgtnames, const DType bcVal = 0.0);

    void compute_vel_v(void);

    void compute_cells_varflux(std::string velname, std::string cellname, const DType tsol_loc, fCalcBC bcFun = NULL);
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

    void compute_chems_dflux(const int nt, const DType tsol_loc);
    void calc_dflux_chem_var(std::string chemname, const int nt, const DType tsol_loc);

    void compute_chems_rflux(const int nt, const DType tsol_loc);
    void calc_rflux_chem_var_upw(std::string chemname, const int nt);
    void calc_rflux_chem_var_ppm(std::string chemname, const int nt);
    void calc_rflux_chem_var_wupw(std::string chemname, const int nt, const DType tsol_loc);

    void compute_chems_flux(void);

    void reset_residual(void);
    void compute_cells_residual(const int fr, const int to);
    void assemble_inertia_var(HVEC& h_res, DType r_fr, HVEC& h_var, DType r_to, SacadoIndVarArray& s_var, DType* l_x_c);
    void assemble_flux_var(HVEC& h_res, SacadoDepVarArray& s_flx);
    void assemble_srce_var(HVEC& h_res, SacadoDepVarArray& s_src);
    void compute_chems_residual(const int fr, const int to);

    void update_cells_by_res(const int nt);
    void update_sol_by_res_var_exp(HVEC& h_var, HVEC& h_res, DType h_r, DType* l_x_c, const int num);

    void update_cells_by_res_imp(const int nt, DType& cur_res);
    void update_cell_linsys_by_res_var(HVEC& h_res, DType h_r, DType* l_x_c, SacadoDepVarArray& s_flx, SacadoDepVarArray& s_src, const int num, const int s_offset, const int l_offset, LSPack& lsPack);
    void update_hsol_by_linsol(HVEC& h_var, const int l_offset, const int num, LSPack& lsPack, DType& cur_res);

    void update_chems_by_res_imp(const int nt, DType& loc_res);
    void update_chem_by_res_var_imp(HVEC& h_var, HVEC& h_res, DType h_r, DType* l_x_c, SacadoDepVarArray& s_flx, SacadoDepVarArray& s_src, const int num, const int offset, LSPack& lsPack, DType& loc_res);

    void compute_rad(const int fr, const int to);

    void update_source_geom(int nt);
    void calc_srce_geom(SacadoIndVarArray& s_var, SacadoDepVarArray& s_srce, DType* l_x_c);
    void calc_srce_geom(SacadoIndVarArray& s_var, SacadoDepVarArray& s_srce, DType* l_x_c, const int nt);

    void compute_source_cell_prob(const DType tsol_loc);
    void calc_srce_cell_base(SacadoDepVarArray& s_srce, const DType tsol_loc);

    void calc_srce_g_base(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_base(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_base(SacadoDepVarArray& s_srce, const DType tsol_loc);

    void calc_srce_g_model(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_model(SacadoDepVarArray& s_srce, const DType tsol_loc);

    void calc_srce_g_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc);

    void calc_srce_g_opn(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_opn(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_opn(SacadoDepVarArray& s_srce, const DType tsol_loc);

    void calc_srce_b_oncolyticOld(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_oncolyticOld(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_oncolyticOld(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_g_oncolyticOld(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_mv_oncolyticOld(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_mb_oncolyticOld(SacadoDepVarArray& s_srce, const DType tsol_loc);

    void calc_srce_b1_potassium(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_b2_potassium(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_b3_potassium(SacadoDepVarArray& s_srce, const DType tsol_loc);

    inline SacadoDepVar icalc_scal_logist_mm(DType,const AutoDiff::Fad::DFad<DType>&,DType,DType,DType,const AutoDiff::Fad::DFad<DType>&);
    inline SacadoDepVar icalc_mm(DType mm_ceil, DType mm_half, const AutoDiff::Fad::DFad<DType>& s_var) {
      SacadoDepVar result = 0.0;
      if ( std::fabs( s_var ) > 1.0 )
        result = mm_ceil * ( 1.0 / ( (mm_half/s_var) + 1.0 ) );
      else if ( s_var != 0.0 )
        result = mm_ceil * ( s_var / (mm_half + s_var) );
      return result;
    }

    int idx_bdry_oncolytic;
    void calc_srce_b_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_q_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    inline SacadoDepVar icalc_q_growth_oncolytic(const AutoDiff::Fad::DFad<DType>& s_q) { SacadoDepVar val = s_q * 0.0; return val;}
    void calc_srce_n_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_m1_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_m1e_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_m2_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_m2e_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_m3_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_m3e_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_s_oncolytic(SacadoDepVarArray& s_srce, const DType tsol_loc);

    void calc_srce_m_onc_model(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_s_onc_model(SacadoDepVarArray& s_srce, const DType tsol_loc);

    void compute_source_chem_prob(const int nt, const DType tsol_loc);
    void calc_srce_chem_var(std::string chemname, const int nt, const DType tsol_loc);

    void calc_chem_srce_vi_oncolyticOld(std::string chemname, const int nt, const DType tsol_loc);
    void calc_chem_srce_a_oncolyticOld(std::string chemname, const int nt, const DType tsol_loc);

    void calc_chem_srce_p_potassium(std::string chemname, const int nt, const DType tsol_loc);
    void calc_chem_srce_v_potassium(std::string chemname, const int nt, const DType tsol_loc);

    void calc_chem_srce_vi_oncolytic(std::string chemname, const int nt, const DType tsol_loc );

    inline SacadoDepVar icalc_scal_block_potassium(DType,DType,const AutoDiff::Fad::DFad<DType>&,DType heaviside_scale);
    //inline SacadoDepVar icalc_scal_block_potassium(DType,DType,AutoDiff::Fad::DFad<DType>,DType heaviside_scale = 1.0);
    // The same function for values only (no derivatives)
    inline DType icalc_scal_block_potassium(DType,DType,DType,DType heaviside_scale);
    // The same function, storing the result in s_val. In a loop this avoids
    // the allocation of a new variable for every result.
    inline void iset_scal_block_potassium(SacadoDepVar& s_val,DType,DType,const AutoDiff::Fad::DFad<DType>&,DType heaviside_scale);

    // Parameters of the potassium model, read once from
    // tdCoef->probVarParVals["potassium"] by setup_par_potassium
    struct PotassiumPar {
      bool  has_rho;  // global/ie_ratio is given
      DType rho;      // global/ie_ratio (1 if not given)
      DType b1_decay; // B1/grow_neg_virus
      DType b1_block; // B1/grow_block (1 if not given)
      DType b2_decay; // B2/grow_neg
      DType b2_block; // B2/grow_block (0 if not given)
      DType b2_thres; // B2/threshold (P/threshold if not given)
      DType p_slope;  // P/slope
      DType p_thres;  // P/threshold
      DType p_burst;  // P/burst
      DType v_burst;  // V/burst
      DType v_block;  // V/burst_block (0 if not given)
    };
    PotassiumPar potPar;
    void setup_par_potassium(void);

    // Other functions
    DType calc_n_bc( const DType tsol_loc );
    DType get_cell_perc(std::string cellname, const int nt);
    DType calc_treatment_P1(const DType tsol_loc);
    DType calc_m1_bc_oncolytic( const DType tsol_loc ); 
    DType calc_m2_bc_oncolytic( const DType tsol_loc ); 
    DType calc_m3_bc_oncolytic( const DType tsol_loc ); 
    SacadoDepVar calc_m1_exter_oncolytic( const int icell, SacadoDepVar r_loc, const DType dr, const DType tsol_loc );
    SacadoDepVar calc_m2_exter_oncolytic( const int icell, SacadoDepVar r_loc, const DType dr, const DType tsol_loc );
    SacadoDepVar calc_m3_exter_oncolytic( const int icell, SacadoDepVar r_loc, const DType dr, const DType tsol_loc );
    SacadoDepVar calc_m_exter_onc_model( SacadoDepVar r_loc, const DType dr, const DType tsol_loc );

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
  // XYZ - NEED TO MODIFY IT FOR GENERAL tdSetup->dim 
  return l_x * ( l_x * l_x - 0.25 * l_dx * l_dx );
}

inline TDSolveSacado::SacadoDepVar TDSolveSacado::icalc_scal_logist_mm(DType scal, const AutoDiff::Fad::DFad<DType>& s_logist, DType mm_const, DType mm_ceil, DType mm_half, const AutoDiff::Fad::DFad<DType>& s_mm) {
  SacadoDepVar s_val;
  if ( s_logist >= 1.0 )
    s_val = 0.0;
  else if ( s_logist <= 0.0 )
    s_val = scal * (mm_const + mm_ceil * s_mm/(mm_half+s_mm));
  else
    s_val = scal * (1.0-s_logist) * (mm_const + mm_ceil * s_mm/(mm_half+s_mm));
  return s_val;
}

inline TDSolveSacado::SacadoDepVar TDSolveSacado::icalc_scal_block_potassium(DType k, DType Ps, const AutoDiff::Fad::DFad<DType>& s_P, DType heaviside_scale) {
  SacadoDepVar s_val;
  iset_scal_block_potassium( s_val, k, Ps, s_P, heaviside_scale );
  return s_val;
}

inline void TDSolveSacado::iset_scal_block_potassium(SacadoDepVar& s_val, DType k, DType Ps, const AutoDiff::Fad::DFad<DType>& s_P, DType heaviside_scale) {
  if ( k > 0.0 )
    s_val = 1.0 - heaviside_scale/(1.0+std::exp(-2.0*k*(s_P-Ps)));
  else if ( k == 0.0 ) {
    s_val = 1.0;
  }
  else {
    if ( s_P < Ps )
      s_val = 1.0;
    else
      s_val = 1.0 - heaviside_scale;
  }
  return;
}

inline DType TDSolveSacado::icalc_scal_block_potassium(DType k, DType Ps, DType P, DType heaviside_scale) {
  DType val;
  if ( k > 0.0 )
    val = 1.0 - heaviside_scale/(1.0+std::exp(-2.0*k*(P-Ps)));
  else if ( k == 0.0 ) {
    val = 1.0;
  }
  else {
    if ( P < Ps )
      val = 1.0;
    else
      val = 1.0 - heaviside_scale;
  }
  return val;
}

#endif
