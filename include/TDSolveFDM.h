#ifndef TDSOLVE_FDM_H
#define TDSOLVE_FDM_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"
#include "TDSolve.h"

class TDSolveFDM : public TDSolve {
  public:
    void setup(void);
    void set_ic(void);
    void calc_dt(const int nt);
    void output(const int nt);
    bool update_sol(const int fr, const int to, const DType shift_dt_fr=0.0, const DType shift_dt_to=1.0);
    void copy_data(const int fr, const int to);
    // Needs to update whenever a or n is modified
    void update_vel(const int fr, const int to, const DType shift_dt=0.0, const bool last_iter=false);


    TDSolveFDM(TDInput* tdInput_ptr);
    ~TDSolveFDM(void);
  
  private:
    void update_cells_upw(const int fr, const int to);
    void update_cells_upw_be(const int fr, const int to);
    void update_cells_pade_be(const int fr, const int to);

    void update_chem_upw_be(const int fr, const int to);
    void update_chem_pade_be(const int fr, const int to);
    
    void update_vel_int(const int nt);
    void update_vel_upw(const int nt);
    void update_vel_bdf2(const int nt);
    void update_vel_cfd(const int nt) {};
};

#endif

