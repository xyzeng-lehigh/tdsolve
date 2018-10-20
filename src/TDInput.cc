#include <fstream>
#include <climits>
#include <algorithm>
#include "TDInput.h"
#include "HighDim_Alg.h"
#include "Alloc_Utils.h"
using namespace std;

TDInput::TDInput(void) {
  tdCoef.lambda      = 1.0;
  tdCoef.mu          = 1.0;
  tdCoef.delta_h     = 1.0;
  tdCoef.alpha       = 1.0;
  tdCoef.theta       = 1.0;
  tdCoef.sigma       = 1.0;
  tdCoef.delta_n     = 0.0;
  tdCoef.D           = 1.0;
  tdCoef.alpha_a     = 1.0;
  tdCoef.beta_a      = 1.0;
  tdCoef.gamma_a     = 1.0;
  tdCoef.n_bc        = 5.0e-3;
  tdCoef.R_0         = 1.0;
  tdCoef.mu_slp      = 0.0;
  tdCoef.q           = 0.0;
  tdCoef.p           = 0.0;
  tdCoef.mu_coef     = 1.0;
  
  //OPN model
  tdCoef.a           = 100.0;
  tdCoef.b           = 100.0;
  tdCoef.D_1         = 0.01;
  tdCoef.D_2         = 0.01;
  tdCoef.a_prime     = 1.0;
  tdCoef.b_prime     = 1.0;
  tdCoef.alpha_1     = 10.0;
  tdCoef.alpha_2     = 10.0;
  tdCoef.beta_1      = 1000.0;
  tdCoef.beta_2      = 1000.0;
  tdCoef.gamma_1     = 1.0;
  tdCoef.gamma_2     = 1.0;
  tdCoef.omega       = 0.0;
  tdCoef.treatP1     = false;
  tdCoef.treatP1BC   = false;
  tdCoef.treatP1_T   = 0.0;
  tdCoef.treatP1_h   = 1.0;
  tdCoef.treatP1_bc  = 0.0;
  tdCoef.treatP1_t0.clear();

  tdMesh.x_l     = 0.0;
  tdMesh.x_r     = 1.0;
  tdMesh.Nx      = 100;
  tdMesh.Na_m    = 2;
  tdMesh.Na      = 200;
  tdMesh.dx      = 0.01;
  tdMesh.x_f     = NULL;
  tdMesh.x_c     = NULL;
  tdMesh.x_a_f   = NULL;
  tdMesh.x_a_c   = NULL;


  tdSetup.frame  = 0;
  tdSetup.cfl    = 0.9;
  tdSetup.ic     = 0;
  tdSetup.schem  = 0;
  tdSetup.scell  = 0;
  tdSetup.vel    = 0;
  tdSetup.lbc    = 0;
  tdSetup.rbc    = 0;
  tdSetup.tmeth  = 1;     // 0 - forward Euler
                          // 1 - RK2
  tdSetup.vCoup  = false; // The next two options are relevant
                          // if vCoup == true
  tdSetup.vIter  = 10;    // -1  -- iterate until converge
                          // >=0 -- up to tIter iterations
  tdSetup.vTol   = 1e-6;  // Iteration tolerance for tIter=-1
  tdSetup.Tinit  = 0.0;
  tdSetup.Tterm  = 1.0;
  tdSetup.Rterm  = 6.0;
  tdSetup.Rdead  = 0.0;
  tdSetup.Tfreq  = -1.0;
  tdSetup.freq   = 1;
  tdSetup.maxit  = INT_MAX;
  tdSetup.dTinit = 1e-10;
  tdSetup.dTgrow = 1.2;
  tdSetup.dr_opt = 0; // 0, default, dr = V
                      // 1, dx independent, dR = V(1+0.5*dt*V/R)
                      // 2, dx dependent, dR = (1/1-dx^2/4)V
  tdSetup.r_opt  = 0; // 0, default, R_{n+1} = R_n + dt*V
                      // 1, consistent, R_{n+1}^2 = R_n^2 + 2*dt*R_n*dR
  tdSetup.slim_tot = false;
  tdSetup.logistic = 0;
  tdSetup.prob_id  = 0;  // 0, default, single chemotaxis
                         // 1, OPN
  tdSetup.init_u_par  = 1.0; 
  tdSetup.diagn_vcoup = true;
  tdSetup.silence     = false; // Turn off all outputs (on-screen or file)

  tdLin.Type      = LSYS_TYPE_DIR;
  tdLin.Print     = true;
  tdLin.Iter      = LSYS_ITER_JAC;
  tdLin.Krylov    = LSYS_KRYLOV_GMRES;
  tdLin.KrylovCap = 50;
  tdLin.Prec      = LSYS_PREC_OFF;
  tdLin.PrecPar   = 0;
  tdLin.TargetRes = 1.0e-8;
  tdLin.MaxIter   = 100;

  strcpy(delims," =:\t");
  strcpy(delims_nospace,"=:");
  strcpy(outputFile,"sol.dat");
  strcpy(radiusFile,"rad.dat");

  setup_parameters();
  setup_options();
  setup_toggles();
}

TDInput::TDInput(const TDInput& tdInput) {
  setup_parameters();
  setup_options();
  setup_toggles();
  strcpy(inputFile,      tdInput.inputFile);
  strcpy(delims,         tdInput.delims);
  strcpy(delims_nospace, tdInput.delims_nospace);
  strcpy(outputFile,     tdInput.outputFile);
  strcpy(radiusFile,     tdInput.radiusFile);

  tdCoef  = tdInput.tdCoef;
  tdSetup = tdInput.tdSetup;
  tdLin   = tdInput.tdLin;

  tdMesh.x_l  = tdInput.tdMesh.x_l;
  tdMesh.x_r  = tdInput.tdMesh.x_r;
  tdMesh.dx   = tdInput.tdMesh.dx;
  tdMesh.Nx   = tdInput.tdMesh.Nx;
  tdMesh.Na_m = tdInput.tdMesh.Na_m;
  tdMesh.Na   = tdInput.tdMesh.Na;

  int Nx = tdMesh.Nx;
  int Na = tdMesh.Na;

  tdMesh.x_f = new DType[Nx+1];
  for ( int i = 0; i < Nx+1; i++ )
    tdMesh.x_f[i] = tdInput.tdMesh.x_f[i];

  tdMesh.x_c = new DType[Nx];
  for ( int i = 0; i < Nx; i++ )
    tdMesh.x_c[i] = tdInput.tdMesh.x_c[i];

  tdMesh.x_a_f = new DType[Na+1];
  for ( int i = 0; i < Na+1; i++ )
    tdMesh.x_a_f[i] = tdInput.tdMesh.x_a_f[i];

  tdMesh.x_a_c = new DType[Na];
  for ( int i = 0; i < Na; i++ )
    tdMesh.x_a_c[i] = tdInput.tdMesh.x_a_c[i];
}

TDInput::~TDInput(void) {
  erase_array( tdMesh.x_f );
  erase_array( tdMesh.x_c );
  erase_array( tdMesh.x_a_c );
  erase_array( tdMesh.x_a_f );
}

void TDInput::setInputFile(const char* file_in) {
  strcpy(inputFile, file_in);
  if (strlen(inputFile)<=0) {
    fprintf(stderr,"Error: Input file not provided.\n");
    exit(1);
  }
  return;
}

void TDInput::parse_old(void) {
  ifstream fin;
  string file_line_s;
  fin.open(inputFile);
  char* name;
  char* value;

  while (getline(fin,file_line_s)) {
    // This ignore a short line
    if (file_line_s.length()<2)
      continue;
    // This indicates the comment lines
    if ((file_line_s[0]=='/')&&(file_line_s[1]=='/'))
      continue;
    // This ignore a non-comment line with no separation tokens
    if (file_line_s.find_first_of((const char*)delims_nospace)==string::npos)
      continue;
    name  = strtok((char*)file_line_s.c_str(),delims);
    value = strtok(NULL,delims);
    if ((name==NULL)||(value==NULL)) {
      fprintf(stderr,"Warning: the following line in the input file %s is not recogonized:\n",inputFile);
      fprintf(stderr,"         %s\n",file_line_s.c_str());
      continue;
    }

    if (strcmp(name,"lambda")==0)
      tdCoef.lambda = atof(value);
    else if (strcmp(name,"mu")==0)
      tdCoef.mu = atof(value);
    else if ( (strcmp(name,"delta")==0) || 
              (strcmp(name,"delta_h")==0) ||
              (strcmp(name,"delta_1")==0) )
      tdCoef.delta_h = atof(value);
    else if (strcmp(name,"alpha")==0)
      tdCoef.alpha = atof(value);
    else if (strcmp(name,"theta")==0)
      tdCoef.theta = atof(value);
    else if (strcmp(name,"sigma")==0)
      tdCoef.sigma = atof(value);
    else if ( (strcmp(name,"delta_n")==0) ||
              (strcmp(name,"delta_2")==0) ||
              (strcmp(name,"rho")==0) )
      tdCoef.delta_n = atof(value);
    else if ( (strcmp(name,"D")==0) ||
              (strcmp(name,"D_3")==0) )
      tdCoef.D = atof(value);
    else if ( (strcmp(name,"m")==0) ||
              (strcmp(name,"alpha_a")==0) ||
              (strcmp(name,"alpha_3")==0) )
      tdCoef.alpha_a = atof(value);
    else if ( (strcmp(name,"beta")==0) ||
              (strcmp(name,"beta_a")==0) ||
              (strcmp(name,"beta_3")==0) )
      tdCoef.beta_a = atof(value);
    else if ( (strcmp(name,"gamma")==0) ||
              (strcmp(name,"gamma_a")==0) ||
              (strcmp(name,"gamma_3")==0) )
      tdCoef.gamma_a = atof(value);
    else if (strcmp(name,"n_bc")==0)
      tdCoef.n_bc= atof(value);
    else if (strcmp(name,"r_0")==0)
      tdCoef.R_0 = atof(value);
    else if (strcmp(name,"mu_slope")==0)
      tdCoef.mu_slp  = atof(value);
    else if (strcmp(name,"a")==0)
      tdCoef.a = atof(value);
    else if (strcmp(name,"b")==0)
      tdCoef.b = atof(value);
    else if (strcmp(name,"a_prime")==0)
      tdCoef.a_prime = atof(value);
    else if (strcmp(name,"b_prime")==0)
      tdCoef.b_prime = atof(value);
    else if (strcmp(name,"D_1")==0)
      tdCoef.D_1 = atof(value);
    else if (strcmp(name,"D_2")==0)
      tdCoef.D_2 = atof(value);
    else if (strcmp(name,"alpha_1")==0)
      tdCoef.alpha_1 = atof(value);
    else if (strcmp(name,"alpha_2")==0)
      tdCoef.alpha_2 = atof(value);
    else if (strcmp(name,"beta_1")==0)
      tdCoef.beta_1 = atof(value);
    else if (strcmp(name,"beta_2")==0)
      tdCoef.beta_2 = atof(value);
    else if (strcmp(name,"gamma_1")==0)
      tdCoef.gamma_1 = atof(value);
    else if (strcmp(name,"gamma_2")==0)
      tdCoef.gamma_2 = atof(value);
    else if (strcmp(name,"q")==0)
      tdCoef.q = atof(value);
    else if (strcmp(name,"p")==0)
      tdCoef.p = atof(value);
    else if (strcmp(name,"mu_coef")==0)
      tdCoef.mu_coef = atof(value);
    else if (strcmp(name,"omega")==0)
      tdCoef.omega = atof(value);
    else if (strcmp(name,"use_treatment_P1")==0)
      tdCoef.treatP1 = (atoi(value)!=0);
    else if (strcmp(name,"use_treatmentBC_P1")==0)
      tdCoef.treatP1BC = (atoi(value)!=0);
    else if (strcmp(name,"treatment_P1_t0")==0)
      tdCoef.treatP1_t0.push_back( atof(value) );
    else if (strcmp(name,"treatment_P1_T")==0)
      tdCoef.treatP1_T = atof(value);
    else if (strcmp(name,"treatment_P1_h")==0)
      tdCoef.treatP1_h = atof(value);
    else if (strcmp(name,"treatment_P1_bc")==0)
      tdCoef.treatP1_bc = atof(value);
    else if (strcmp(name,"x_l")==0)
      tdMesh.x_l = atof(value);
    else if (strcmp(name,"x_r")==0)
      tdMesh.x_r = atof(value);
    else if (strcmp(name,"Nx")==0)
      tdMesh.Nx = atoi(value);
    else if (strcmp(name,"A_Multiplier")==0)
      tdMesh.Na_m = atoi(value);
    else if (strcmp(name,"framework")==0)
      tdSetup.frame = atoi(value);
    else if (strcmp(name,"cfl")==0)
      tdSetup.cfl = atof(value);
    else if (strcmp(name,"ic")==0)
      tdSetup.ic = atoi(value);
    else if (strcmp(name,"space_cell")==0)
      tdSetup.scell = atoi(value);
    else if (strcmp(name,"space_chem")==0)
      tdSetup.schem = atoi(value);
    else if (strcmp(name,"vel")==0)
      tdSetup.vel = atoi(value);
    else if (strcmp(name,"lbc")==0)
      tdSetup.lbc = atoi(value); // 0 -- hom. Neumann; 1 -- ODE bc.
    else if (strcmp(name,"rbc")==0)
      tdSetup.rbc = atoi(value); // 0 -- hom. Neumann; 1 -- ODE bc.
    else if (strcmp(name,"time_method")==0)
      tdSetup.tmeth = atoi(value);
    else if (strcmp(name,"vel_couple")==0)
      tdSetup.vCoup = (atoi(value)!=0);
    else if (strcmp(name,"vel_iter")==0)
      tdSetup.vIter = atoi(value);
    else if (strcmp(name,"vel_iter_tol")==0)
      tdSetup.vTol = atof(value);
    else if (strcmp(name,"Tinit")==0)
      tdSetup.Tinit = atof(value);
    else if (strcmp(name,"Tterm")==0)
      tdSetup.Tterm = atof(value);
    else if (strcmp(name,"Rterm")==0)
      tdSetup.Rterm = atof(value);
    else if (strcmp(name,"Rdead")==0)
      tdSetup.Rdead = atof(value);
    else if (strcmp(name,"Tfreq")==0)
      tdSetup.Tfreq = atof(value);
    else if (strcmp(name,"frequency")==0)
      tdSetup.freq = atoi(value);
    else if (strcmp(name,"max_iter")==0)
      tdSetup.maxit = atoi(value);
    else if (strcmp(name,"dt_init")==0)
      tdSetup.dTinit = atof(value);
    else if (strcmp(name,"dt_growth")==0)
      tdSetup.dTgrow = atof(value);
    else if (strcmp(name,"dr_option")==0)
      tdSetup.dr_opt = atoi(value);
    else if (strcmp(name,"r_option")==0)
      tdSetup.r_opt = atoi(value);
    else if (strcmp(name,"lim_sync_tot")==0)
      tdSetup.slim_tot = (atoi(value)!=0);
    else if (strcmp(name,"logistic")==0)
      tdSetup.logistic = atoi(value);
    else if (strcmp(name,"prob_id")==0) 
      tdSetup.prob_id = atoi(value);
    else if (strcmp(name,"init_u_par")==0)
      tdSetup.init_u_par = atof(value);
    else if (strcmp(name,"diagnostic_vcoup")==0)
      tdSetup.diagn_vcoup = (atoi(value)!=0);
    else if (strcmp(name,"silence")==0)
      tdSetup.silence = (atoi(value)!=0);
    else if (strcmp(name,"lsysType")==0)
      if (strcmp(value,"DIR")==0)
        tdLin.Type = LSYS_TYPE_DIR;
      else if (strcmp(value,"ITER")==0)
        tdLin.Type = LSYS_TYPE_ITER;
      else if (strcmp(value,"KRYLOV")==0)
        tdLin.Type = LSYS_TYPE_KRYLOV;
      else
        tdLin.Type = LSYS_TYPE_DIR;
    else if (strcmp(name,"lsysIter")==0)
      if (strcmp(value,"Jacobi")==0)
        tdLin.Iter = LSYS_ITER_JAC;
      else if (strcmp(value,"GaussSeidel")==0)
        tdLin.Iter = LSYS_ITER_GS;
      else if (strcmp(value,"BackGaussSeidel")==0)
        tdLin.Iter = LSYS_ITER_BGS;
      else if (strcmp(value,"SOR")==0)
        tdLin.Iter = LSYS_ITER_SOR;
      else if (strcmp(value,"BSOR")==0)
        tdLin.Iter = LSYS_ITER_BSOR;
      else
        tdLin.Iter = LSYS_ITER_JAC;
    else if (strcmp(name,"lsysKrylov")==0)
      if (strcmp(value,"CG")==0)
        tdLin.Krylov = LSYS_KRYLOV_CG;
      else if (strcmp(value,"GMRES")==0)
        tdLin.Krylov = LSYS_KRYLOV_GMRES;
      else if (strcmp(value,"LANCZOS")==0)
        tdLin.Krylov = LSYS_KRYLOV_LANCZOS;
      else if (strcmp(value,"QMR")==0)
        tdLin.Krylov = LSYS_KRYLOV_QMR;
      else
        tdLin.Krylov = LSYS_KRYLOV_GMRES;
    else if (strcmp(name,"lsysPrec")==0)
      if (strcmp(value,"OFF")==0)
        tdLin.Prec = LSYS_PREC_OFF;
      else if (strcmp(value,"External")==0)
        tdLin.Prec = LSYS_PREC_LINK;
      else if (strcmp(value,"Eye")==0)
        tdLin.Prec = LSYS_PREC_EYE;
      else if (strcmp(value,"ILU")==0)
        tdLin.Prec = LSYS_PREC_ILU;
      else if (strcmp(value,"BILU")==0)
        tdLin.Prec = LSYS_PREC_BILU;
      else
        tdLin.Prec = LSYS_PREC_OFF;
    else if (strcmp(name,"lsysPrint")==0)
      tdLin.Print = (atoi(value)!=0);
    else if (strcmp(name,"lsysTargetRes")==0)
      tdLin.TargetRes = atof(value);
    else if (strcmp(name,"lsysPrecPar")==0)
      tdLin.PrecPar = atoi(value);
    else if (strcmp(name,"lsysMaxIter")==0)
      tdLin.MaxIter = atoi(value);
    else if (strcmp(name,"lsysKrylovCap")==0)
      tdLin.KrylovCap = atoi(value);
    else if (strcmp(name,"outputDataFile")==0)
      strcpy(outputFile,(const char*)value);
    else if (strcmp(name,"radiusDataFile")==0)
      strcpy(radiusFile,(const char*)value);
    else if (strcmp(name,"recompute_vel")==0 ||
             strcmp(name,"bc_lim") == 0 ) {
      // Obsolete options, just to get rid of the warning below
    }
    else {
      fprintf(stderr,"Warning: the following line in the input file %s is not recogonized:\n",inputFile);
      fprintf(stderr,"         %s\n",file_line_s.c_str());
      continue;
    }
  }

  if ( tdCoef.treatP1_t0.size() > 0 )
    std::sort( tdCoef.treatP1_t0.begin(), tdCoef.treatP1_t0.end() );

  tdMesh.dx = (tdMesh.x_r-tdMesh.x_l)/(1.0*tdMesh.Nx);
  alloc_array_full( tdMesh.x_f, tdMesh.Nx+1 );
  alloc_array_full( tdMesh.x_c, tdMesh.Nx   );

  tdMesh.Na = tdMesh.Nx * tdMesh.Na_m;
  alloc_array_full( tdMesh.x_a_f, tdMesh.Na + 1 );
  alloc_array_full( tdMesh.x_a_c, tdMesh.Na );

  for ( int i = 0; i < tdMesh.Nx+1; i++ )
    tdMesh.x_f[i] = tdMesh.x_l + i * tdMesh.dx;
  for ( int i = 0; i < tdMesh.Nx; i++ )
    tdMesh.x_c[i] = tdMesh.x_l + (i+0.5) * tdMesh.dx;

  for ( int i = 0; i < tdMesh.Na+1; i++ )
    tdMesh.x_a_f[i] = tdMesh.x_l + i * tdMesh.dx;
  for ( int i = 0; i < tdMesh.Na; i++ )
    tdMesh.x_a_c[i] = tdMesh.x_l + (i+0.5) * tdMesh.dx;

  fin.close();
  return;
}

void TDInput::setup_parameters(void) {
  tdParameters.clear();
  tdParameters["lambda"]          = &tdCoef.lambda;
  tdParameters["mu"]              = &tdCoef.mu;
  tdParameters["delta"]           = &tdCoef.delta_h;
  tdParameters["delta_h"]         = &tdCoef.delta_h;
  tdParameters["delta_1"]         = &tdCoef.delta_h;
  tdParameters["alpha"]           = &tdCoef.alpha;
  tdParameters["theta"]           = &tdCoef.theta;
  tdParameters["sigma"]           = &tdCoef.sigma;
  tdParameters["delta_n"]         = &tdCoef.delta_n;
  tdParameters["delta_2"]         = &tdCoef.delta_n;
  tdParameters["rho"]             = &tdCoef.delta_n;
  tdParameters["D"]               = &tdCoef.D;
  tdParameters["D_3"]             = &tdCoef.D;
  tdParameters["m"]               = &tdCoef.alpha_a;
  tdParameters["alpha_a"]         = &tdCoef.alpha_a;
  tdParameters["alpha_3"]         = &tdCoef.alpha_a;
  tdParameters["beta"]            = &tdCoef.beta_a;
  tdParameters["beta_a"]          = &tdCoef.beta_a;
  tdParameters["beta_3"]          = &tdCoef.beta_a;
  tdParameters["gamma"]           = &tdCoef.gamma_a;
  tdParameters["gamma_a"]         = &tdCoef.gamma_a;
  tdParameters["gamma_3"]         = &tdCoef.gamma_a;
  tdParameters["n_bc"]            = &tdCoef.n_bc;
  tdParameters["r_0"]             = &tdCoef.R_0;
  tdParameters["mu_slope"]        = &tdCoef.mu_slp;
  tdParameters["a"]               = &tdCoef.a;
  tdParameters["b"]               = &tdCoef.b;
  tdParameters["a_prime"]         = &tdCoef.a_prime;
  tdParameters["b_prime"]         = &tdCoef.b_prime;
  tdParameters["D_1"]             = &tdCoef.D_1;
  tdParameters["D_2"]             = &tdCoef.D_2;
  tdParameters["alpha_1"]         = &tdCoef.alpha_1;
  tdParameters["alpha_2"]         = &tdCoef.alpha_2;
  tdParameters["beta_1"]          = &tdCoef.beta_1;
  tdParameters["beta_2"]          = &tdCoef.beta_2;
  tdParameters["gamma_1"]         = &tdCoef.gamma_1;
  tdParameters["gamma_2"]         = &tdCoef.gamma_2;
  tdParameters["q"]               = &tdCoef.q;
  tdParameters["p"]               = &tdCoef.p;
  tdParameters["mu_coef"]         = &tdCoef.mu_coef;
  tdParameters["omega"]           = &tdCoef.omega;
  tdParameters["treatment_P1_T"]  = &tdCoef.treatP1_T;
  tdParameters["treatment_P1_h"]  = &tdCoef.treatP1_h;
  tdParameters["treatment_P1_bc"] = &tdCoef.treatP1_bc;
  tdParameters["x_l"]             = &tdMesh.x_l;
  tdParameters["x_r"]             = &tdMesh.x_r;
  tdParameters["cfl"]             = &tdSetup.cfl;
  tdParameters["vel_iter_tol"]    = &tdSetup.vTol;
  tdParameters["Tinit"]           = &tdSetup.Tinit;
  tdParameters["Tterm"]           = &tdSetup.Tterm;
  tdParameters["Rterm"]           = &tdSetup.Rterm;
  tdParameters["Rdead"]           = &tdSetup.Rdead;
  tdParameters["Tfreq"]           = &tdSetup.Tfreq;
  tdParameters["dt_init"]         = &tdSetup.dTinit;
  tdParameters["dt_growth"]       = &tdSetup.dTgrow;
  tdParameters["init_u_par"]      = &tdSetup.init_u_par;
  tdParameters["lsysTargetRes"]   = &tdLin.TargetRes;

  return;
}

void TDInput::setup_options(void) {
  tdOptions.clear();
  tdOptions["Nx"]            = &tdMesh.Nx;
  tdOptions["A_Multiplier"]  = &tdMesh.Na_m;
  tdOptions["framework"]     = &tdSetup.frame;
  tdOptions["ic"]            = &tdSetup.ic;
  tdOptions["space_cell"]    = &tdSetup.scell;
  tdOptions["space_chem"]    = &tdSetup.schem;
  tdOptions["vel"]           = &tdSetup.vel;
  tdOptions["lbc"]           = &tdSetup.lbc;
  tdOptions["rbc"]           = &tdSetup.rbc;
  tdOptions["time_method"]   = &tdSetup.tmeth;
  tdOptions["vel_iter"]      = &tdSetup.vIter;
  tdOptions["frequency"]     = &tdSetup.freq;
  tdOptions["max_iter"]      = &tdSetup.maxit;
  tdOptions["dr_option"]     = &tdSetup.dr_opt;
  tdOptions["r_option"]      = &tdSetup.r_opt;
  tdOptions["logistic"]      = &tdSetup.logistic;
  tdOptions["prob_id"]       = &tdSetup.prob_id;
  tdOptions["lsysPrecPar"]   = &tdLin.PrecPar;
  tdOptions["lsysMaxIter"]   = &tdLin.MaxIter;
  tdOptions["lsysKrylovCap"] = &tdLin.KrylovCap;

  return;
}

void TDInput::setup_toggles(void) {
  tdToggles.clear();
  tdToggles["use_treatment_P1"]   = &tdCoef.treatP1;
  tdToggles["use_treatmentBC_P1"] = &tdCoef.treatP1BC;
  tdToggles["vel_couple"]         = &tdSetup.vCoup;
  tdToggles["lim_sync_tot"]       = &tdSetup.slim_tot;
  tdToggles["diagnostic_vcoup"]   = &tdSetup.diagn_vcoup;
  tdToggles["silence"]            = &tdSetup.silence;
  tdToggles["lsysPrint"]          = &tdLin.Print;
  return;
}

void TDInput::parse(void) {
  ifstream fin;
  string file_line_s;
  fin.open(inputFile);
  char* name;
  char* value;

  while (getline(fin,file_line_s)) {
    // This ignore a short line
    if (file_line_s.length()<2)
      continue;
    // This indicates the comment lines
    if ((file_line_s[0]=='/')&&(file_line_s[1]=='/'))
      continue;
    // This ignore a non-comment line with no separation tokens
    if (file_line_s.find_first_of((const char*)delims_nospace)==string::npos)
      continue;
    name  = strtok((char*)file_line_s.c_str(),delims);
    value = strtok(NULL,delims);
    if ((name==NULL)||(value==NULL)) {
      fprintf(stderr,"Warning: the following line in the input file %s is not recogonized:\n",inputFile);
      fprintf(stderr,"         %s\n",file_line_s.c_str());
      continue;
    }

    std::map<std::string, DType*>::iterator it_par = tdParameters.find(name);
    std::map<std::string, int*>::iterator   it_opt = tdOptions.find(name);
    std::map<std::string, bool*>::iterator  it_tog = tdToggles.find(name);

    if ( it_par != tdParameters.end() ) 
      *(it_par->second) = atof(value);
    else if ( it_opt != tdOptions.end() )
      *(it_opt->second) = atoi(value);
    else if ( it_tog != tdToggles.end() )
      *(it_tog->second) = (atoi(value)!=0);
    else if (strcmp(name,"treatment_P1_t0")==0)
      tdCoef.treatP1_t0.push_back( atof(value) );
    else if (strcmp(name,"lsysType")==0)
      if (strcmp(value,"DIR")==0)
        tdLin.Type = LSYS_TYPE_DIR;
      else if (strcmp(value,"ITER")==0)
        tdLin.Type = LSYS_TYPE_ITER;
      else if (strcmp(value,"KRYLOV")==0)
        tdLin.Type = LSYS_TYPE_KRYLOV;
      else
        tdLin.Type = LSYS_TYPE_DIR;
    else if (strcmp(name,"lsysIter")==0)
      if (strcmp(value,"Jacobi")==0)
        tdLin.Iter = LSYS_ITER_JAC;
      else if (strcmp(value,"GaussSeidel")==0)
        tdLin.Iter = LSYS_ITER_GS;
      else if (strcmp(value,"BackGaussSeidel")==0)
        tdLin.Iter = LSYS_ITER_BGS;
      else if (strcmp(value,"SOR")==0)
        tdLin.Iter = LSYS_ITER_SOR;
      else if (strcmp(value,"BSOR")==0)
        tdLin.Iter = LSYS_ITER_BSOR;
      else
        tdLin.Iter = LSYS_ITER_JAC;
    else if (strcmp(name,"lsysKrylov")==0)
      if (strcmp(value,"CG")==0)
        tdLin.Krylov = LSYS_KRYLOV_CG;
      else if (strcmp(value,"GMRES")==0)
        tdLin.Krylov = LSYS_KRYLOV_GMRES;
      else if (strcmp(value,"LANCZOS")==0)
        tdLin.Krylov = LSYS_KRYLOV_LANCZOS;
      else if (strcmp(value,"QMR")==0)
        tdLin.Krylov = LSYS_KRYLOV_QMR;
      else
        tdLin.Krylov = LSYS_KRYLOV_GMRES;
    else if (strcmp(name,"lsysPrec")==0)
      if (strcmp(value,"OFF")==0)
        tdLin.Prec = LSYS_PREC_OFF;
      else if (strcmp(value,"External")==0)
        tdLin.Prec = LSYS_PREC_LINK;
      else if (strcmp(value,"Eye")==0)
        tdLin.Prec = LSYS_PREC_EYE;
      else if (strcmp(value,"ILU")==0)
        tdLin.Prec = LSYS_PREC_ILU;
      else if (strcmp(value,"BILU")==0)
        tdLin.Prec = LSYS_PREC_BILU;
      else
        tdLin.Prec = LSYS_PREC_OFF;
    else if (strcmp(name,"outputDataFile")==0)
      strcpy(outputFile,(const char*)value);
    else if (strcmp(name,"radiusDataFile")==0)
      strcpy(radiusFile,(const char*)value);
    else if (strcmp(name,"recompute_vel")==0 ||
             strcmp(name,"bc_lim") == 0 ) {
      // Obsolete options, just to get rid of the warning below
    }
    else {
      fprintf(stderr,"Warning: the following line in the input file %s is not recogonized:\n",inputFile);
      fprintf(stderr,"         %s\n",file_line_s.c_str());
      continue;
    }
  }

  if ( tdCoef.treatP1_t0.size() > 0 )
    std::sort( tdCoef.treatP1_t0.begin(), tdCoef.treatP1_t0.end() );

  tdMesh.dx = (tdMesh.x_r-tdMesh.x_l)/(1.0*tdMesh.Nx);
  alloc_array_full( tdMesh.x_f, tdMesh.Nx+1 );
  alloc_array_full( tdMesh.x_c, tdMesh.Nx   );

  tdMesh.Na = tdMesh.Nx * tdMesh.Na_m;
  alloc_array_full( tdMesh.x_a_f, tdMesh.Na + 1 );
  alloc_array_full( tdMesh.x_a_c, tdMesh.Na );

  for ( int i = 0; i < tdMesh.Nx+1; i++ )
    tdMesh.x_f[i] = tdMesh.x_l + i * tdMesh.dx;
  for ( int i = 0; i < tdMesh.Nx; i++ )
    tdMesh.x_c[i] = tdMesh.x_l + (i+0.5) * tdMesh.dx;

  for ( int i = 0; i < tdMesh.Na+1; i++ )
    tdMesh.x_a_f[i] = tdMesh.x_l + i * tdMesh.dx;
  for ( int i = 0; i < tdMesh.Na; i++ )
    tdMesh.x_a_c[i] = tdMesh.x_l + (i+0.5) * tdMesh.dx;

  fin.close();
  return;
}

bool TDInput::set_parameter(std::string name, DType value) {
  std::map<std::string, DType*>::iterator it = tdParameters.find(name);
  if ( it == tdParameters.end() ) {
    return false;
  }
  else {
    *(it->second) = value;
    return true;
  }
}

bool TDInput::get_parameter(std::string name, DType& value) {
  std::map<std::string, DType*>::iterator it = tdParameters.find(name);
  if ( it == tdParameters.end() )
    return false;
  else {
    value = *(it->second);
    return true;
  }
}

DType TDInput::get_parameter(std::string name) {
  std::map<std::string, DType*>::iterator it = tdParameters.find(name);
  if ( it == tdParameters.end() )
    return 0.0;
  else
    return *(it->second);
}

bool TDInput::set_option(std::string name, int value) {
  std::map<std::string, int*>::iterator it = tdOptions.find(name);
  if ( it == tdOptions.end() )
    return false;
  else {
    *(it->second) = value;
    return true;
  }
}

bool TDInput::get_option(std::string name, int& value) {
  std::map<std::string, int*>::iterator it = tdOptions.find(name);
  if ( it == tdOptions.end() )
    return false;
  else {
    value = *(it->second);
    return true;
  }
}

int TDInput::get_option(std::string name) {
  std::map<std::string, int*>::iterator it = tdOptions.find(name);
  if ( it == tdOptions.end() )
    return 0;
  else
    return *(it->second);
}

bool TDInput::set_toggle(std::string name, bool value) {
  std::map<std::string, bool*>::iterator it = tdToggles.find(name);
  if ( it == tdToggles.end() )
    return false;
  else {
    *(it->second) = value;
    return true;
  }
}

bool TDInput::get_toggle(std::string name, bool& value) {
  std::map<std::string, bool*>::iterator it = tdToggles.find(name);
  if ( it == tdToggles.end() )
    return false;
  else {
    value = *(it->second);
    return true;
  }
}

bool TDInput::get_toggle(std::string name) {
  std::map<std::string, bool*>::iterator it = tdToggles.find(name);
  if ( it == tdToggles.end() )
    return true;
  else
    return *(it->second);
}

void TDInput::display_parameter(void) {
  cout << "TDInput parameters: " << std::endl;
  for ( std::map<std::string, DType*>::iterator it = tdParameters.begin();
        it != tdParameters.end(); it++ )
    cout << "\t" << it->first << ": " << *(it->second) << std::endl;
  return;
}

void TDInput::display_option(void) {
  cout << "TDInput options: " << std::endl;
  for ( std::map<std::string, int*>::iterator it = tdOptions.begin();
        it != tdOptions.end(); it++ )
    cout << "\t" << it->first << ": " << *(it->second) << std::endl;
  return;
}

void TDInput::display_toggle(void) {
  cout << "TDInput toggles: " << std::endl;
  for ( std::map<std::string, bool*>::iterator it = tdToggles.begin();
        it != tdToggles.end(); it++ )
    cout << "\t" << it->first << ": " << ((*(it->second)) ? "True" : "False") << std::endl;
  return;
}

