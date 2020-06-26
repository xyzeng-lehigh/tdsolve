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

    void alloc(void);
    bool update_sol(const int fr, const int to, const DType shift_dt_fr=0.0, const DType shift_dt_to=1.0);
    void copy_data(const int fr, const int to);
    // Needs to update whenever a or n is modified
    // Before calling, the S_ independent arrays should be associated with fr
    // After calling, (*data_v)[to] will be filled
    void update_vel(const int fr, const int to, const DType shift_dt=0.0, const bool last_iter=false);

  protected:
    // Independent variables
    typedef std::vector<Sacado::Fad::DFad<DType>> SacadoIndVarArray;
    std::map<std::string, SacadoIndVarArray> S_data_var;
    std::map<std::string, int> S_data_offset;
    std::map<std::string, int> S_data_size;
    int S_data_num;

    // Dependent variables
    typedef Sacado::Fad::DFad<DType> SacadoDepVar;
    typedef std::vector<SacadoDepVar> SacadoDepVarArray;
    typedef std::map<std::string, SacadoDepVarArray> SacadoDepVarMap;
    // Velocities 
    SacadoDepVarMap S_data_vel;
    // Cell numbers
    HVEC** data_inc;  // Incompressibility violation
    SacadoDepVarMap S_flux_cell;
    SacadoDepVarMap S_flux_cell_v;
    SacadoDepVarMap S_flux_cell_r;
    SacadoDepVarMap S_flux_cell_u;
    SacadoDepVarMap S_lim_ppm_cell[2];
    SacadoDepVarArray S_lim_ppm_sync[2];
    //std::vector<Sacado::Fad::DFad<DType>> S_lim_cell_sync; // MUSCL
    //SacadoDepVarMap S_slope_cell; // MUSCL

    // Chemoattractant/OPN 
    SacadoDepVarMap S_flux_chem;
    //std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_slope_chem; // MUSCL or PPM
    //std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_slope_p;

    // Initial condition for three cell numbers
    void setic_unif_cell(DType val_g, DType val_h, DType val_n, int nt=0);
    void setic_expd_chem(DType multip, DType decay, int Na, int nt=0);

    // Setup Sacado arrays
    void alloc_sacado_ind(void); // Independent variables
    void alloc_sacado_ind_array(std::string, const int);
    void alloc_sacado_dep(void); // Dependent variables
    void alloc_sacado_dep_array(SacadoDepVarArray&, const int);

    // Other Sacado utilities
    void set_sacado_ind(const int);
    void set_sacado_ind_array(std::string, HVEC&);
    void set_sacado_ind_array(std::string, DType);

  private:
    void calc_total_inc_viol(const int nt);
    void compute_vel_u(const DType tsol_loc);
    void compute_vel_u_pres(const DType tsol_loc);
    void compute_vel_u_chem(void);

};

#endif
