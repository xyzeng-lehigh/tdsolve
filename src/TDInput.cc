#include <fstream>
#include <climits>
#include "TDInput.h"
#include "HighDim_Alg.h"
#include "Alloc_Utils.h"
using namespace std;

TDInput::TDInput(void) {
  tdCoef.lambda = 1.0;
  tdCoef.mu     = 1.0;
  tdCoef.delta  = 1.0;
  tdCoef.alpha  = 1.0;
  tdCoef.theta  = 1.0;
  tdCoef.sigma  = 1.0;
  tdCoef.D      = 1.0;
  tdCoef.m      = 1.0;
  tdCoef.beta   = 1.0;
  tdCoef.gamma  = 1.0;
  tdCoef.n_bc   = 5.0e-3;

  tdMesh.x_l    = 0.0;
  tdMesh.x_r    = 1.0;
  tdMesh.Nx     = 100;
  tdMesh.Na_m   = 2;
  tdMesh.Na     = 200;
  tdMesh.dx     = 0.01;
  tdMesh.x_f    = NULL;
  tdMesh.x_c    = NULL;
  tdMesh.x_a_f  = NULL;
  tdMesh.x_a_c  = NULL;

  tdSetup.bmark = false;
  tdSetup.cfl   = 0.9;
  tdSetup.ic    = 0;
  tdSetup.schem = 0;
  tdSetup.scell = 0;
  tdSetup.vel   = 0;
  tdSetup.rbc   = 0;
  tdSetup.Tinit = 0.0;
  tdSetup.Tterm = 1.0;
  tdSetup.freq  = 1;
  tdSetup.maxit = INT_MAX;
  tdSetup.dTinit = 1e-10;
  tdSetup.dTgrow = 1.2;

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
}

TDInput::~TDInput(void) {
  erase_array( tdMesh.x_f );
  erase_array( tdMesh.x_c );
}

void TDInput::setInputFile(const char* file_in) {
  strcpy(inputFile, file_in);
  if (strlen(inputFile)<=0) {
    fprintf(stderr,"Error: Input file not provided.\n");
    exit(1);
  }
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

    if (strcmp(name,"lambda")==0)
      tdCoef.lambda = atof(value);
    else if (strcmp(name,"mu")==0)
      tdCoef.mu = atof(value);
    else if (strcmp(name,"delta")==0)
      tdCoef.mu = atof(value);
    else if (strcmp(name,"alpha")==0)
      tdCoef.alpha = atof(value);
    else if (strcmp(name,"theta")==0)
      tdCoef.theta = atof(value);
    else if (strcmp(name,"sigma")==0)
      tdCoef.sigma= atof(value);
    else if (strcmp(name,"D")==0)
      tdCoef.D = atof(value);
    else if (strcmp(name,"m")==0)
      tdCoef.m = atof(value);
    else if (strcmp(name,"beta")==0)
      tdCoef.beta = atof(value);
    else if (strcmp(name,"gamma")==0)
      tdCoef.gamma = atof(value);
    else if (strcmp(name,"n_bc")==0)
      tdCoef.n_bc= atof(value);
    else if (strcmp(name,"x_l")==0)
      tdMesh.x_l = atof(value);
    else if (strcmp(name,"x_r")==0)
      tdMesh.x_r = atof(value);
    else if (strcmp(name,"Nx")==0)
      tdMesh.Nx = atoi(value);
    else if (strcmp(name,"A_Multiplier")==0)
      tdMesh.Na_m = atoi(value);
    else if (strcmp(name,"benchmark")==0)
      tdSetup.bmark = (atoi(value) == 1);
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
    else if (strcmp(name,"rbc")==0)
      tdSetup.rbc = atoi(value); // 0 -- hom. Neumann; 1 -- ODE bc.
    else if (strcmp(name,"Tinit")==0)
      tdSetup.Tinit = atof(value);
    else if (strcmp(name,"Tterm")==0)
      tdSetup.Tterm = atof(value);
    else if (strcmp(name,"frequency")==0)
      tdSetup.freq = atoi(value);
    else if (strcmp(name,"max_iter")==0)
      tdSetup.maxit = atoi(value);
    else if (strcmp(name,"dt_init")==0)
      tdSetup.dTinit = atof(value);
    else if (strcmp(name,"dt_growth")==0)
      tdSetup.dTgrow = atof(value);
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
    else {
      fprintf(stderr,"Warning: the following line in the input file %s is not recogonized:\n",inputFile);
      fprintf(stderr,"         %s\n",file_line_s.c_str());
      continue;
    }
  }

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

