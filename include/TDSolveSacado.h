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

  protected:
    struct ChemPar {
      DType decay;
      DType diffusion;
      // Michaelis Mentun coefficients
      DType mm_ceil, mm_half;
      std::string driver;
      bool dirichlet_bc;
    };
    // Variable names
    std::vector<std::string> cellnames;   // All cell species name
    std::vector<std::string> cellnames_u; // Cell species subject to infiltration
    std::vector<std::string> velnames;    // All velocity names
    std::vector<std::string> chemnames;   // All chemoattractant/OPN names
    std::map<std::string, int> varsizes;
    std::map<std::string, DType*> varcoord;
    std::map<std::string, DType*> varflxcoord; // Only exist for cells/chems
    std::map<std::string, ChemPar> chempars;

    // Independent variables (cell numbers and chemicals)
    typedef std::vector<Sacado::Fad::DFad<DType>> SacadoIndVarArray;
    std::map<std::string, SacadoIndVarArray> S_data_var;
    std::map<std::string, int> S_data_offset;
    int S_data_num;
    // HVEC version
    std::map<std::string, HVEC**> H_data_var;
    std::map<std::string, HVEC**> H_data_res;

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
    void setup_chem_par( std::string chemname, const DType decay, const DType diffusion,
        std::string cellname, const DType mm_ceil, const DType mm_half, const bool dbc );

    // Initial condition for three cell numbers
    void setic_unif_cell(DType val_g, DType val_h, DType val_n, int nt=0);
    void setic_expd_chem(DType multip, DType decay, int Na, int nt=0);

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
    void set_sacado_ind_array(std::string, const int);
    void set_sacado_r(const int);
    void set_sacado_v(const int);
    void set_from_sacado_ind(const int);
    void set_from_sacado_ind_array(std::string, HVEC&);
    void set_from_sacado_r(const int);
    void set_from_sacado_v(const int);
    void reset_sacado_ind_array(void);
    void reset_sacado_ind_array(SacadoIndVarArray& s_var, const int s_offset, DType val=0.0);
    void reset_sacado_dep_array(void);
    void reset_sacado_dep_array(SacadoDepVarMapMap&, DType val=0.0);
    void reset_sacado_dep_array(SacadoDepVarMap&, DType val=0.0);
    void reset_sacado_dep_array(SacadoDepVarArray&, DType val=0.0);
    void sum_in_sacado_dep_array(SacadoDepVarArray&, DType, SacadoDepVarArray&);
    void print_sacado_ind_array(SacadoIndVarArray&);
    void print_sacado_dep_array(SacadoDepVarArray&);

  private:
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
    void calc_ppm_lim_var(SacadoIndVarArray& s_var, SacadoDepVarArray& s_rec, SacadoDepVarArray& s_lim_l, SacadoDepVarArray& s_lim_r, DType* l_x_c, const bool tosync=false, std::vector<bool>* set_zero=NULL, SacadoDepVarArray* s_alpha=NULL);
    void calc_ppm_lim_sync(SacadoDepVarArray* s_sync, SacadoDepVarMap* s_lims, std::vector<bool>& set_zero, SacadoDepVarArray* s_alpha);
    void calc_ppm_flux_var(SacadoIndVarArray& s_var, SacadoDepVarArray& s_rec, SacadoDepVarArray& s_flx, SacadoDepVarArray* s_lim, DType* l_x_c);

    void compute_cells_flux(void);

    void compute_chems_dflux(void);
    void calc_dflux_chem_var(std::string chemname);

    void compute_chems_rflux(void);
    void calc_rflux_chem_var_upw(std::string chemname);
    void calc_rflux_chem_var_ppm(std::string chemname);

    void compute_chems_flux(void);

    void reset_residual(void);
    void compute_cells_residual(const int fr, const int to);
    void assemble_inertia_var(HVEC& h_res, DType r_fr, HVEC& h_var, SacadoDepVar& s_r_to, SacadoIndVarArray s_var, DType* l_x_c);
    void assemble_flux_var(HVEC& h_res, SacadoDepVarArray& s_flx);
    void assemble_srce_var(HVEC& h_res, SacadoDepVarArray& s_src);
    void compute_chems_residual(const int fr, const int to);

    void update_sol_by_res(const int nt);
    void update_sol_by_res_var(HVEC& h_var, HVEC& h_res, DType h_r, DType* l_x_c, const int num);

    void compute_rad(const int fr, const int to);

    void update_source_geom(void);
    void calc_srce_geom(SacadoIndVarArray& s_var, SacadoDepVarArray& s_srce, DType* l_x_c);

    void compute_source_cell_prob(const DType tsol_loc);
    void calc_srce_g_base(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_base(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_base(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_g_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_pdgf(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_g_opn(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_h_opn(SacadoDepVarArray& s_srce, const DType tsol_loc);
    void calc_srce_n_opn(SacadoDepVarArray& s_srce, const DType tsol_loc);

    void compute_source_chem_prob(const DType tsol_loc);
    void calc_srce_chem_var(std::string chemname, const DType tsol_loc);

    // Other functions
    DType calc_n_bc( const DType tsol_loc );

};

#endif
