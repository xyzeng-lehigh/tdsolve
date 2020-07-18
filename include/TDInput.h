#ifndef TDINPUT_H
#define TDINPUT_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include "Data_Type.h"
#include "Comput_Choices.h"

// Coefficients table
// TDCoef       PDGF (default)      OPN
// lambda       lambda              lambda 
// mu           mu                  mu
// delta_h      delta               delta_1
// alpha        alpha               alpha
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
// q                                q
// mu_coef                          mu_coef
// p                                p

struct TDCoef {
  DType lambda, mu, delta_h, alpha, theta, sigma;
  DType delta_n;
  DType mu_slp;
  DType D, alpha_a, beta_a, gamma_a;
  DType n_bc;
  DType R_0;
  DType q;
  DType mu_coef;
  DType p;

  // Some coefficients for user defined ics
  DType ic_sine_freq;

  // Additional parameters of the OPN model
  DType a, b, D_1, D_2; 
  DType a_prime, b_prime;
  DType alpha_1, alpha_2, beta_1, beta_2, gamma_1, gamma_2;
  DType omega;
  DType treatP1_T, treatP1_h, treatP1_bc;
  std::vector<DType> treatP1_t0;
  bool  treatP1, treatP1BC;
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
  int   assemb;
  bool  chmppm;
  int   ppm_reset_num;

  // PMC for iterates (velocity, radius, etc.)
  bool  vCoup;
  int   vIter; 
  DType vTol;
  int   nrmPow;
  bool  vOnly;

  // Newton iterations
  DType ntTol;
  int   ntIter;

  DType Tinit, Tterm, Tfreq;
  DType Rterm, Rdead;
  int   freq;
  int   maxit;
  DType dTinit, dTgrow, dTmax;
  int   dr_opt, r_opt;
  bool  jac_lim;
  bool  jac_vel;
  bool  slim_tot;
  bool  diagn_vcoup;
  bool  silence_screen;
  bool  silence_output;
  bool  auto_red_dt;

  bool  use_pmc;

  int   logistic;
  int   prob_id;

  // Test cases
  int   infil_u;
  DType init_u_par;
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

    std::map<std::string, DType*> tdParameters;
    std::map<std::string, int*>   tdOptions;
    std::map<std::string, bool*>  tdToggles;

    bool set_parameter(std::string, DType);
    bool set_option(std::string, int);
    bool set_toggle(std::string, bool);
    bool get_parameter(std::string, DType&);
    bool get_option(std::string, int&);
    bool get_toggle(std::string, bool&);
    DType get_parameter(std::string);
    int   get_option(std::string);
    bool  get_toggle(std::string);

    void  display_parameter(void);
    void  display_option(void);
    void  display_toggle(void);

  private:
    void setup_parameters(void);
    void setup_options(void);
    void setup_toggles(void);

  public:
    TDInput(void);
    TDInput(const TDInput&);
    ~TDInput(void);
    
    void  setInputFile(const char*);
    void  parse(void);
};

#endif
