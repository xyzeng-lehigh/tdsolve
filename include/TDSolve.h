#ifndef TDSOLVE_H
#define TDSOLVE_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"

class  HVEC;
struct TDCoef;
struct TDMesh;
struct TDSetup;
struct TDLin;
class TDInput;
class LinMat;
class LinVec;
class LinSolve;

class TDSolve {
  private:
    HVEC** data_g;
    HVEC** data_h;
    HVEC** data_n;
    HVEC** data_a;
    HVEC** data_v;
    DType  data_r[3];

    HVEC** data_na;

    TDInput* tdInput;
    TDCoef*  tdCoef;
    TDMesh*  tdMesh;
    TDSetup* tdSetup;
    TDLin*   tdLin;

    int    Nx;
    DType  dx;
    DType* x_f;
    DType* x_c;

    DType  dtn, tsol;
    bool   reach_term;
    // Needs to update whenever a or n is modified
    void update_na(const int nt); 
    void update_vel(const int nt);

    void update_rad(const int fr, const int to);
    void update_cells_upw(const int fr, const int to);
    void update_chem_upw_be(const int fr, const int to);
    void update_vel_cd(const int nt);

    // Linear solver
    LinMat*   mat;
    LinMat*   prc;
    LinVec*   sol;
    LinVec*   rhs;
    LinVec*   res;
    LinSolve* linSolve;

    int iter;

  public:
    TDSolve(TDInput*);
    ~TDSolve(void);

    void setup(void);
    void set_ic(void);
    void update_sol(const int fr, const int to);
    void average_sol(const DType, const int, 
                     const DType, const int, const int);
    void copy_data(const int fr, const int to);
    void calc_dt(const int nt);
    void diagnostic(const int nt);
    void output(const int nt);


    bool nextIter(void) {
      iter++;
      return !reach_term;
    }
};

#endif

