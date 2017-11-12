#ifndef TDINPUT_H
#define TDINPUT_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include "Data_Type.h"
#include "Comput_Choices.h"

// Coefficients table
// TDCoef       PDGF (default)      OPN
// lambda       lambda              lambda 
// mu           mu                  mu
// delta_h      delta               delta_1
// alpha        alpha               rho
// theta        theta               theta
// sigma        sigma
// delta_n      rho                 delta_2
// D            D                   D_3
// alpha_a      m                   alpha_3
// beta_a       beta                beta_3
// gamma_a      gamma               gamma_3
// D_1                              D_1
// alpha_1                          alpha_1
// beta_1                           beta_1
// gamma_1                          gamma_1
// D_2                              D_2
// alpha_2                          alpha_2
// beta_2                           beta_2
// gamma_2                          gamma_2
// a                                a
// a_prime                          a_prime
// b                                b
// b_prime                          b_prime

struct TDCoef {
  DType lambda, mu, delta_h, alpha, theta, sigma;
  DType delta_n;
  DType mu_slp;
  DType D, alpha_a, beta_a, gamma_a;
  DType n_bc;
  DType R_0;

  // Additional parameters of the OPN model
  DType a, b, D_1, D_2; 
  DType a_prime, b_prime;
  DType alpha_1, alpha_2, beta_1, beta_2, gamma_1, gamma_2;
};

struct TDMesh {
  DType  x_l, x_r;
  DType  dx;
  int    Nx; // Np1 & Np2 is same as Nx in the OPN model
  DType* x_f;
  DType* x_c;

  int    Na_m;
  int    Na;
  DType* x_a_f;
  DType* x_a_c;
};

struct TDSetup {
  int   frame;
  DType cfl;
  int   ic;
  int   scell, schem, vel, rbc, lbc;
  int   tmeth;
  DType Tinit, Tterm, Tfreq;
  int   freq;
  int   maxit;
  DType dTinit, dTgrow;

  bool  bc_lim;
  bool  nonDim;

  int   prob_id;
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
