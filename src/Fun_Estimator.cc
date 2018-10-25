#include "Fun_Estimator.h"
#include <cmath>
#include "TDInput.h"
#include "HighDim_Alg.h"
#include <iomanip>
#include <fstream>

FunEstimator::FunEstimator(char* fitMeth) {
  data.clear();
  ikeys.clear();
  irefs.clear();
  okeys.clear();
  orefs.clear();
  ofuns.clear();
  basis.clear();
  fcoef.clear();

  if ( std::strcmp(fitMeth,"LS0D")==0 ||
       std::strcmp(fitMeth,"LS0T")==0 ) {
    poly_tensor = (fitMeth[3] == 'T');
    ls_degree = 0;
    estimate = &FunEstimator::est_ls0;}
  else if ( std::strcmp(fitMeth,"LS1D")==0 ||
            std::strcmp(fitMeth,"LS1T")==0 ) {
    poly_tensor = (fitMeth[3] == 'T');
    ls_degree = 1;
    estimate = &FunEstimator::est_ls1;}
  else if ( std::strcmp(fitMeth,"LS2D")==0 || 
            std::strcmp(fitMeth,"LS2T")==0 ) {
    poly_tensor = (fitMeth[3] == 'T');
    ls_degree = 2;
    estimate = &FunEstimator::est_ls2;}
  else if ( std::strcmp(fitMeth,"LS3D")==0 ||
            std::strcmp(fitMeth,"LS3T")==0 ) {
    poly_tensor = (fitMeth[3] == 'T');
    ls_degree = 3;
    estimate = &FunEstimator::est_ls3;}
  else if ( std::strcmp(fitMeth,"LS4D")==0 ||
            std::strcmp(fitMeth,"LS4T")==0 ) {
    poly_tensor = (fitMeth[3] == 'T');
    ls_degree = 4;
    estimate = &FunEstimator::est_ls4;}
  else if ( std::strcmp(fitMeth,"LS5D")==0 ||
            std::strcmp(fitMeth,"LS5T")==0 ) {
    poly_tensor = (fitMeth[3] == 'T');
    ls_degree = 5;
    estimate = &FunEstimator::est_ls5;}
  else { // Default
    poly_tensor = false;
    ls_degree = 2;
    estimate = &FunEstimator::est_ls2;}
}

void FunEstimator::readData(const char* dataFile) {

  std::ifstream fin;
  std::string file_line_s;
  fin.open(dataFile);
  char* word;
  char  delims[40];
  std::strcpy(delims, " \t");

  while (getline(fin,file_line_s)) {
    std::vector<DType> loc_input;
    std::vector<DType> loc_output;
    loc_input.clear();
    loc_output.clear();

    bool line_start = false;
    for ( std::set<std::string>::iterator it = ikeys.begin();
          it != ikeys.end(); it++ ) {
      DType var_i;
      if ( !line_start ) {
        word = strtok((char*)file_line_s.c_str(),delims);
        line_start = true;
      }
      else
        word = strtok(NULL,delims);
      var_i = std::atof(word);
      loc_input.push_back( var_i/irefs[*it] );
    }

    getline(fin,file_line_s);
    line_start = false;
    for ( std::set<std::string>::iterator it = okeys.begin();
          it != okeys.end(); it++ ) {
      DType var_o;
      if ( !line_start ) {
        word = strtok((char*)file_line_s.c_str(),delims);
        if ( word[0] == 'U' ) {
          word = strtok(NULL,delims);
          break;
        }
        else
          var_o = std::atof(word);
        line_start = true;
      }
      else {
        word = strtok(NULL,delims);
        var_o = std::atof(word);
      }
      loc_output.push_back( var_o/orefs[*it] );
    }
    if ( loc_output.size() != 0 )
      data.push_back( VarPair( loc_input, loc_output ) );
  }

  fin.close();

  return;
}

void FunEstimator::preprocess(void) {
  for ( std::vector<VarPair>::iterator it = data.begin();
        it != data.end(); it++ ) {
    std::vector<DType>& oVars = it->second;
    unsigned cur_i = 0;
    for ( std::set<std::string>::iterator it_key = okeys.begin();
          it_key != okeys.end(); it_key++, cur_i++ ) {
      switch ( ofuns[*it_key] ) {
        case FUN_POLY:
          // Do nothing
          break;
        case FUN_EXP_POLY:
          oVars[cur_i] = std::log( oVars[cur_i] );
          break;
        default:
          // Do nothing
          break;
      }
    }
  }
  return;
}

void FunEstimator::getReferenceInput(const bool use_ref, TDInput& tdInput) {
  irefs.clear();
  for ( std::set<std::string>::iterator it = ikeys.begin(); it != ikeys.end(); it++ ) {
    DType vref;
    if ( !(tdInput.get_parameter(*it, vref)) )
      fprintf(stderr, "Error: No reference variable %s found... Exit.\n", (*it).c_str());
    irefs[*it] = use_ref ? vref : 1.0;
    //cout << " (" << *it << ") reference value: " << vref << std::endl;
  }
  return;
}

void FunEstimator::getReferenceOutput(std::map<std::string,DType>& ref_vals) {
  for ( std::set<std::string>::iterator it = okeys.begin(); it != okeys.end(); it++ ) {
    std::map<std::string,DType>::iterator it_ref = ref_vals.find( *it );
    if ( it_ref == ref_vals.end() )
      orefs[*it] = 1.0;
    else
      orefs[*it] = it_ref->second;
    //cout << " (" << *it << ") reference value: " << orefs[*it] << std::endl;
  }
  return;
}

void FunEstimator::getFittingMethod(std::map<std::string,std::string>& ref_funs) {
  for ( std::set<std::string>::iterator it = okeys.begin(); it != okeys.end(); it++ ) {
    std::map<std::string,std::string>::iterator it_ref = ref_funs.find( *it );
    if ( it_ref == ref_funs.end() )
      ofuns[*it] = FUN_POLY;
    else if ( it_ref->second == "Poly" )
      ofuns[*it] = FUN_POLY;
    else if ( it_ref->second == "ExpPoly" )
      ofuns[*it] = FUN_EXP_POLY;
    else
      ofuns[*it] = FUN_POLY;
    //cout << " (" << *it << ") fitting function: " << ofuns[*it] << std::endl;
  }
  return;
}

void FunEstimator::est_ls(const unsigned p) {
  setup_ls_basis(p);

  unsigned Nfunc = basis.size();

  HMAT mat(Nfunc);
  HVEC** rhs = NULL;
  alloc_array_full(rhs, Nfunc, okeys.size());
  reset_array(rhs, okeys.size(), 0.0);

  mat = 0.0;

  for ( std::vector<VarPair>::iterator it = data.begin();
        it != data.end(); it++ ) {
    std::vector<DType> pvals = eval_poly_basis(p, it->first);
    for ( unsigned i = 0; i < Nfunc; i++ )
      for ( unsigned j = 0; j < Nfunc; j++ )
        mat[i][j] += pvals[i] * pvals[j];
    std::vector<DType>& ovals = it->second;
    for ( unsigned i = 0; i < ovals.size(); i++ ) {
      DType loc_oval = ovals[i];
      for ( unsigned j = 0; j < Nfunc; j++ )
        (*rhs[i])[j] += loc_oval * pvals[j];
    }
  }

  char UPLO = 'L';
  int  N    = Nfunc;
  int  NRHS = okeys.size();
  int  LDA  = Nfunc;
  int  LDB  = Nfunc;
  int  INFO;
#ifndef FLOAT_ALG
  double* A = new double[Nfunc*Nfunc];
  double* B = new double[Nfunc*okeys.size()];
  for ( unsigned i = 0; i < Nfunc; i++ )
    for ( unsigned j = 0; j < Nfunc; j++ )
      A[i*Nfunc+j] = mat[i][j];
  for ( unsigned i = 0; i < okeys.size(); i++ )
    for ( unsigned j = 0; j < Nfunc; j++ )
      B[i*Nfunc+j] = (*rhs[i])[j];
  dposv_(&UPLO, &N, &NRHS, &(A[0]), &LDA, &(B[0]), &LDB, &INFO);
  for ( unsigned i = 0; i < okeys.size(); i++ )
    for ( unsigned j = 0; j < Nfunc; j++ )
      (*rhs[i])[j] = B[i*Nfunc+j];
  delete[] A;
  delete[] B;
#else
  float*  A = new float[Nfunc*Nfunc];
  float*  B = new float[Nfunc*okeys.size()];
  for ( unsigned i = 0; i < Nfunc; i++ )
    for ( unsigned j = 0; j < Nfunc; j++ )
      A[i*Nfunc+j] = mat[i][j];
  for ( unsigned i = 0; i < okeys.size(); i++ )
    for ( unsigned j = 0; j < Nfunc; j++ )
      B[i*Nfunc+j] = (*rhs[i])[j];
  sposv_(&UPLO, &N, &NRHS, &(A[0]), &LDA, &(B[0]), &LDB, &INFO);
  for ( unsigned i = 0; i < okeys.size(); i++ )
    for ( unsigned j = 0; j < Nfunc; j++ )
      (*rhs[i])[j] = B[i*Nfunc+j];
  delete[] A;
  delete[] B;
#endif

  int cur_i = 0;
  for ( std::set<std::string>::iterator it = okeys.begin();
        it != okeys.end(); it++, cur_i++ ) {
    std::vector<DType> coef;
    coef.resize(Nfunc,0.0);
    for ( unsigned i = 0; i < Nfunc; i++ )
      coef[i] = (*rhs[cur_i])[i];
    fcoef[*it] = coef;
  }

  erase_array(rhs, okeys.size());

  return;
}

void FunEstimator::setup_ls_basis(const unsigned p) {
  basis.clear();
  unsigned Nfunc = 1;
  unsigned Nvar  = ikeys.size();
  for (unsigned i = 0; i < Nvar; i++) Nfunc *= (p+1);
  for (unsigned n = 0; n < Nfunc; n++) {
    PolyTerm pterm;
    pterm.resize(Nvar,0);
    unsigned loc_n = n;
    unsigned tot_deg = 0;
    for (unsigned pos = 0; pos < Nvar; pos++) {
      pterm[pos] = loc_n % (p+1);
      loc_n = (loc_n - pterm[pos])/(p+1);
      tot_deg += pterm[pos];
    }
    if ( !poly_tensor && tot_deg > p )
      Nfunc -= 1;
    else
      basis.push_back(pterm);
  }
  return;
}

std::vector<DType> FunEstimator::eval_poly_basis(const unsigned p, std::vector<DType>& ivars) {
  std::vector<DType> pvals;
  pvals.resize(basis.size(), 0.0);

  DType** pbasis = NULL;
  alloc_array_full( pbasis, ikeys.size(), p+1 );
  for (unsigned i = 0; i < ikeys.size(); i++ ) {
    pbasis[i][0] = 1.0;
    for (unsigned q = 0; q < p; q++)
      pbasis[i][q+1] = pbasis[i][q] * ivars[i];
  }

  for (unsigned n = 0; n < basis.size(); n++ ) {
    PolyTerm pterm = basis[n];
    DType fval = 1.0;
    for (unsigned i = 0; i < ikeys.size(); i++)
      fval *= pbasis[i][pterm[i]];
    pvals[n] = fval;
  }

  erase_array( pbasis, ikeys.size() );
  return pvals;
}

void FunEstimator::write(const char* oFile) {
  FILE* fp = fopen(oFile, "w");
  write_header( fp );
  for ( std::set<std::string>::iterator it = okeys.begin();
        it != okeys.end(); it++ )
    write_variable( fp, *it );

  fclose( fp );
  return;
}

void FunEstimator::write_header( FILE* fp ) {
  fprintf(fp, "Input variables: ");
  for ( std::set<std::string>::iterator it = ikeys.begin();
        it != ikeys.end(); it++ )
    fprintf(fp, "%s ", (*it).c_str());
  fprintf(fp, "\n");
  fprintf(fp, "Output variables: ");
  for ( std::set<std::string>::iterator it = okeys.begin();
        it != okeys.end(); it++ )
    fprintf(fp, "%s ", (*it).c_str());
  fprintf(fp, "\n\n");
  return;
}

void FunEstimator::write_variable( FILE* fp, const std::string key ) {
  fprintf( fp, "%s: ", key.c_str() );
  switch ( ofuns[key] ) {
    case FUN_POLY:
      fprintf( fp, "multivariate polynomial with degree %u in %u variables.\n", 
               ls_degree, ikeys.size() );
      break;
    case FUN_EXP_POLY:
      fprintf( fp, "exponential of multivariate polynomial with degree %u in %u variables.\n", 
               ls_degree, ikeys.size() );
      break;
    default:
      fprintf( fp, "non-recognized fitting type.\n" );
      return;
  }
  switch ( ofuns[key] ) {
    case FUN_POLY:
    case FUN_EXP_POLY:
      {
        unsigned cur_i = 0;
        std::vector<DType>& pcoef = fcoef[key];
        for ( std::vector<PolyTerm>::iterator it = basis.begin();
              it != basis.end(); it++, cur_i++ ) {
          fprintf( fp, "\t[" );
          PolyTerm pterm = *it;
          for ( unsigned i = 0; i < pterm.size(); i++ ) {
            if ( i != pterm.size()-1 )
              fprintf( fp, "%u,", pterm[i] );
            else
              fprintf( fp, "%u]", pterm[i] );
          }
          fprintf( fp, ": %lf\n", pcoef[cur_i] );
        }
        fprintf( fp, "\n" );
      }
      break;
    default:
      break;
  }
  return;
}
