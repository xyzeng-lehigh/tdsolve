#ifndef TDSOLVE_FVM_H
#define TDSOLVE_FVM_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"
#include "TDSolve.h"

class TDSolveFVM : public TDSolve {
  public:
    void setup(void);
    void set_ic(void);
    void calc_dt(const int nt);
    void output(const int nt);
    void update_sol(const int fr, const int to);
    void copy_data(const int fr, const int to);

    TDSolveFVM(TDInput* tdInput_ptr);
    ~TDSolveFVM(void);
  protected:
    HVEC** flux_g;
    HVEC** flux_h;
    HVEC** flux_n;

    // Needs to update whenever a or n is modified
    void update_vel(const int nt);
  private:
    void update_vel_int(const int nt);
    void update_vel_int_mod(const int nt);
    void update_vel_int_test(const int nt);

    int* N_dir;
    void update_cells_upw(const int fr, const int to);
    void update_cells_upw_test(const int fr, const int to);
    void update_cells_upw_be(const int fr, const int to);

    void update_chem_upw_be(const int fr, const int to);

    void upscale_var(HVEC&, HVEC&, DType*&, int, DType);
    void downscale_var(HVEC&, HVEC&, DType*&, int, DType);

    // Try the non-dimensionalized version
    void set_ic_dim(void);
    void update_sol_dim(const int fr, const int to);
    void update_chem_upw_be_dim(const int fr, const int to);
    void update_vel_dim(const int nt);
};

#endif

