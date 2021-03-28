#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "Sacado.hpp"
#define DType double

int main(int argc, char** argv) {
  int Ndup = atoi(argv[1]);
  int Nvar = atoi(argv[2]);

  // Independent variables
  std::vector<Sacado::Fad::DFad<DType>> S_vars;
  for ( int i = 0; i < Nvar; i++ )
    S_vars.push_back( Sacado::Fad::DFad<DType>(Nvar,i,0.0) );
  std::vector<std::vector<Sacado::Fad::DFad<DType>>> S_AllVars;
  for ( int n = 0; n < Ndup; n++ )
    S_AllVars.push_back( S_vars );
  std::vector<Sacado::Fad::DFad<DType>> S_vars_red;

  // Assign value to the Ndup copies of S_vars
  for ( int n = 0; n < Ndup; n++ )
    for ( int i = 0; i < Nvar; i++ ) {
      Sacado::Fad::DFad<DType>& s_var = S_AllVars[n][i];
      s_var = Sacado::Fad::DFad<DType>(Nvar,i,(double)i);
    }

  std::cout << "Independent variables: " << std::endl;
  for ( int n = 0; n < Ndup; n++ ) {
    std::cout << "\tGroup " << n << ":" << std::endl;
    for ( int i = 0; i < Nvar; i++ ) {
      std::cout << "\t\t" << S_AllVars[n][i] << std::endl;
    }
  }

  // Dependent variables
  std::vector<Sacado::Fad::DFad<DType>> S_vals;
  for ( int n = 0; n < Ndup; n++ ) {
    Sacado::Fad::DFad<DType> s_val(0.0);
    for ( int i = 0; i < Nvar; i++ )
      s_val += std::pow(S_AllVars[n][i], 1+n);
    S_vals.push_back( s_val );
  }

  std::cout << "Dependent variables: " << std::endl;
  for ( int n = 0; n < Ndup; n++ )
    std::cout << "\t" << S_vals[n] << std::endl;

  // Reduced size
  if ( Nvar > 1 ) {
    Nvar += 1;
    // Independent variables 
    std::vector<Sacado::Fad::DFad<DType>> S_vars;
    for ( int i = 0; i < Nvar; i++ )
      S_vars.push_back( Sacado::Fad::DFad<DType>(Nvar,i,0.0) );
    std::vector<std::vector<Sacado::Fad::DFad<DType>>> S_AllVars;
    for ( int n = 0; n < Ndup; n++ )
      S_AllVars.push_back( S_vars );
    std::vector<Sacado::Fad::DFad<DType>> S_vars_red;

    // Assign value to the Ndup copies of S_vars
    for ( int n = 0; n < Ndup; n++ )
      for ( int i = 0; i < Nvar; i++ ) {
        Sacado::Fad::DFad<DType>& s_var = S_AllVars[n][i];
        s_var = Sacado::Fad::DFad<DType>(Nvar,i,(double)i);
      }

    std::cout << "Independent variables: " << std::endl;
    for ( int n = 0; n < Ndup; n++ ) {
      std::cout << "\tGroup " << n << ":" << std::endl;
      for ( int i = 0; i < Nvar; i++ ) {
        std::cout << "\t\t" << S_AllVars[n][i] << std::endl;
      }
    }

    // Dependent variables
    std::vector<Sacado::Fad::DFad<DType>> S_vals;
    for ( int n = 0; n < Ndup; n++ ) {
      Sacado::Fad::DFad<DType> s_val(0.0);
      for ( int i = 0; i < Nvar; i++ )
        s_val += std::pow(S_AllVars[n][i], 1+n);
      S_vals.push_back( s_val );
    }

    std::cout << "Dependent variables: " << std::endl;
    for ( int n = 0; n < Ndup; n++ )
      std::cout << "\t" << S_vals[n] << std::endl;
  }

  return 0;
}
