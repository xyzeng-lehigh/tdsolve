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
    void update_sol(const int fr, const int to, const int shift_dt=0);
    void copy_data(const int fr, const int to);

    TDSolveFVM(TDInput* tdInput_ptr);
    ~TDSolveFVM(void);

  protected:
    HVEC** slope_g;
    HVEC** slope_h;
    HVEC** slope_n; // or M
    HVEC** flux_g;
    HVEC** flux_h;
    HVEC** flux_n; // or M
    HVEC** flux_a;
    HVEC** slope_a;

    //OPN model
    HVEC** flux_p1;
    HVEC** slope_p1;
    HVEC** flux_p2;
    HVEC** slope_p2;

    // Needs to update whenever a or n is modified
    void update_vel(const int nt);

    void upscale_var(HVEC&, HVEC&, DType*&, int, DType, int scal_pow=2);
    void downscale_var(HVEC&, HVEC&, DType*&, int, DType, int scal_pow=2);

    // Various output function
    void output_pdgf(const int nt);
    void output_opn(const int nt);
    virtual void print_custom_vel_u( FILE*, const int nt ) {};

    void update_chem_upw_be(const int fr, const int to);
    void update_chem_muscl_imex(const int fr, const int to);
    // OPN model -- to be reused between P1 and P2
    void update_chem_opn_upw_be(const DType tsol_loc, const int fr, const int to, 
                                const DType coef_D, const DType coef_alpha, 
                                const DType coef_beta, const DType coef_gamma, 
                                HVEC**& data_source, HVEC**& data_P,
                                const bool treatment=false, const bool treatBC=false);
    void update_chem_opn_muscl_imex(const DType tsol_loc, const int fr, const int to, 
                                const DType coef_D, const DType coef_alpha, 
                                const DType coef_beta, const DType coef_gamma, 
                                HVEC**& data_source, HVEC**& data_P, HVEC**& flux_P, HVEC**& slope_P, 
                                const bool treatment=false, const bool treatBC=false);

    virtual void update_cells_upw(const int fr, const int to);
    virtual void update_cells_muscl(const int fr, const int to);
    // OPN model
    virtual void update_cells_opn_add(const int fr, const int to); 

    // Treatment for OPN model
    DType calc_treatment(const DType tsol_loc);

  private:
    int*   N_dir;

    void update_vel_int(const int nt);
    void update_vel_int_mod(const int nt);
    //OPN model
    void update_vel_prbid1(const int nt);

    void update_cells_upw_test(const int fr, const int to);
    void update_cells_upw_be(const int fr, const int to);
    void update_cells_muscl_be(const int fr, const int to);

};

#endif

