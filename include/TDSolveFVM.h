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
    void alloc(void);
    bool update_sol(const int fr, const int to, const DType shift_dt_fr=0.0, const DType shift_dt_to=1.0);
    void copy_data(const int fr, const int to);
    // Needs to update whenever a or n is modified
    void update_vel(const int fr, const int to, const DType shift_dt=0.0, const bool last_iter=false);

    // Functions called by API
    DType get_gperc(const int nt);
    DType get_nperc(const int nt);

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
    virtual void update_cells_opn_add(const int fr, const int to, const bool only_res = false ); 

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

