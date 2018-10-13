#include "Fun_Estimator.h"
#include <cmath>
#include "TDInput.h"

FunEstimator::FunEstimator(char* fitMeth, char* fitType) {
  Ninput  = 0;
  Noutput = 0;
  data.clear();
  refs.clear();
  o_refs.clear();

  if ( std::strcpy(fitMeth,"LS0")==0 )
    estimate = &FunEstimator::est_ls0;
  else if ( std::strcpy(fitMeth,"LS1")==0 )
    estimate = &FunEstimator::est_ls1;
  else if ( std::strcpy(fitMeth,"LS2")==0 )
    estimate = &FunEstimator::est_ls2;
  else if ( std::strcpy(fitMeth,"LS3")==0 )
    estimate = &FunEstimator::est_ls3;
  else if ( std::strcpy(fitMeth,"LS4")==0 )
    estimate = &FunEstimator::est_ls4;
  else if ( std::strcpy(fitMeth,"LS5")==0 )
    estimate = &FunEstimator::est_ls5;
  else // Default
    estimate = &FunEstimator::est_ls0;

  if ( std::strcpy(fitType,"Poly") ) {
    preprocess = &FunEstimator::prec_none;
    write      = &FunEstimator::write_poly;
  }
  else if ( std::strcpy(fitType,"ExpPOly") ) {
    preprocess = &FunEstimator::prec_log;
    write      = &FunEstimator::write_exp_poly;
  }
  else { // Default
    preprocess = &FunEstimator::prec_none;
    write      = &FunEstimator::write_poly;
  }
}

void FunEstimator::readData(const char* dataFile) {

  FILE* iFile = fopen(dataFile, "r");

  while ( feof(iFile) ) {
    std::vector<DType> loc_input;
    std::vector<DType> loc_output;
    loc_input.clear();
    loc_output.clear();
    for ( int i = 0; i < Ninput; i++ ) {
      DType var_i;
      fscanf( iFile, "%lf\t", &var_i );
      // Normalizing the inputs
      loc_input.push_back( var_i/refs[i] );
    }
    fscanf( iFile, "\n" );
    for ( int i = 0; i < Noutput; i++ ) {
      DType var_o;
      fscanf( iFile, "%lf\t", &var_o );
      if ( o_refs[i] == 0.0 )
        loc_output.push_back( var_o );
      else // Normalize the output
        loc_output.push_back( var_o/o_refs[i] );
    }
    data.push_back( VarPair( loc_input, loc_output ) );
  }

  fclose(iFile);

  return;
}

void FunEstimator::prec_log(void) {

  for ( std::vector<VarPair>::iterator it = data.begin(); 
        it != data.end(); it++ ) {
    std::vector<DType>& oVars = it->second;
    for ( std::vector<DType>::iterator varIt = oVars.begin();
          varIt != oVars.end(); varIt++ )
      *varIt = std::log(*varIt);
  }

  return;
}

void FunEstimator::getReferenceInput(TDInput& tdInput, std::set<std::string>& keys) {
  refs.clear();
  for ( std::set<std::string>::iterator it = keys.begin(); it != keys.end(); it++ ) {
    DType vref;
    if ( !(tdInput.get_parameter(*it, vref)) )
      fprintf(stderr, "Error: No reference variable %s found... Exit.\n", (*it).c_str());
    refs.push_back(vref);
  }
  return;
}

void FunEstimator::getReferenceOutput(std::map<std::string,DType>& ref_vals, std::set<std::string>& okeys) {
  for ( std::set<std::string>::iterator it = okeys.begin(); it != okeys.end(); it++ ) {
    std::map<std::string,DType>::iterator it_ref = ref_vals.find( *it );
    if ( it_ref == ref_vals.end() )
      o_refs.push_back( 0.0 );
    else
      o_refs.push_back( it_ref->second );
  }
  return;
}

void FunEstimator::est_ls(const unsigned n) {
}

void FunEstimator::write_poly(const char* oFile) {
}

void FunEstimator::write_exp_poly(const char* oFile) {
}

