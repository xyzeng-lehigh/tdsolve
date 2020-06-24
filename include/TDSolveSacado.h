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
    void update_vel(const int fr, const int to, const DType shift_dt=0.0, const bool last_iter=false);

  protected:
    // Independent variables
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_data_var;
    std::map<std::string, int> S_data_offset;

    // Infiltration velocity
    std::vector<Sacado::Fad::DFad<DType>> S_data_u;
    // Cell numbers
    HVEC** data_inc;  // Incompressibility violation
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_flux_cell;
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_flux_cell_v;
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_flux_cell_r;
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_flux_cell_u;
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_lim_cell;
    std::vector<Sacado::Fad::DFad<DType>> S_lim_cell_sync;
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_slope_cell;

    // Chemotaxis attractants
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_slope_chem;
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_flux_chem;

    // OPN
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_flux_p;
    std::map<std::string, std::vector<Sacado::Fad::DFad<DType>>> S_slope_p;

    // Initial condition for three cell numbers
    void setic_unif_cell(DType val_g, DType val_h, DType val_n, int nt=0);
    void setic_expd_chem(DType multip, DType decay, int Na, int nt=0);

    // Setup Sacado arrays
    void alloc_sacado_ind(void); // Independent variables
    void alloc_sacado_ind_array(std::string, const int);
    void alloc_sacado_dep(void); // Dependent variables


  private:
    void calc_total_inc_viol(const int nt);
    void compute_vel_u(const int nt, const DType tsol_loc);
    void compute_vel_u_pres(const int nt, const DType tsol_loc);
    void compute_vel_u_chem(const int nt );

};

#endif
