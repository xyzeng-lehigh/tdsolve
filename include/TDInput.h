#ifndef TDINPUT_H
#define TDINPUT_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"
#include "Comput_Choices.h"

struct TDCoef {
  DType lambda, mu, delta, alpha, theta, sigma;
  DType D, m, beta, gamma;
  DType n_bc;
};

struct TDMesh {
  DType  x_l, x_r;
  DType  dx;
  int    Nx;
  DType* x_f;
  DType* x_c;

  int    Na_m;
  int    Na;
  DType* x_a_f;
  DType* x_a_c;
};

struct TDSetup {
  bool  bmark;
  DType cfl;
  int   ic;
  int   scell, schem, vel, rbc;
  DType Tinit, Tterm;
  int   freq;
  int   maxit;
  DType dTinit, dTgrow;
};

struct TDLin {
  ChoiceLSysType   Type;
  bool             Print;
  ChoiceLSysIter   Iter;
  ChoiceLSysKrylov Krylov;
  int              KrylovCap;
  ChoiceLSysPrec   Prec;
  int              PrecPar;
  DType            TargetRes;
  int              MaxIter;
};

class TDInput {
  public:
    char  inputFile[255];
    char  delims[40];
    char  delims_nospace[40];
    char  outputFile[255];
    char  radiusFile[255];

    TDCoef  tdCoef;
    TDMesh  tdMesh;
    TDSetup tdSetup;
    TDLin   tdLin;

  public:
    TDInput(void);
    ~TDInput(void);
    
    void  setInputFile(const char*);
    void  parse(void);
};

#endif
