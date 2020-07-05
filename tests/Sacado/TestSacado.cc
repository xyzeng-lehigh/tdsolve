#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "Sacado.hpp"
#define DType double

int main(int argc, char** argv) {
  int num = atoi(argv[1]);
  
  // Independent variables
  std::vector<Sacado::Fad::DFad<DType>> S_vars;

  std::cout << "Stage 1: Initializing independent variables" << std::endl;

  for ( int n = 0; n < num; n++ ) {
    Sacado::Fad::DFad<DType> value(num,n,(double)n+1.0);
    S_vars.push_back(value);
  }

  for ( int n = 0; n < num; n++ ) {
    std::cout << "Ind variable " << n << ": " << S_vars[n] << std::endl;
  }

  // Dependent variables
  Sacado::Fad::DFad<DType> S_fin(1.0);
  std::cout << "Initializing a dependent variable f with value 1.0: " << S_fin << std::endl;
  std::cout << "Calculating f as 1 + x.*x:" << std::endl;
  for ( int n = 0; n < num; n++ ) {
    S_fin += S_vars[n] * S_vars[n];
    std::cout << "f after summing in " << n+1 << " variables: " << S_fin << std::endl;
  }

  std::cout << "Stage 2: Re-assign the independent variables" << std::endl;
  
  for ( int n = 0; n < num; n++ ) {
    //S_vars[n] = -(double)n-1.0;
    S_vars[n] = Sacado::Fad::DFad<DType>(num,n,-1.0-(double)n);
    std::cout << "After re-assigning " << n+1 << " independent variables:" << std::endl;
    for ( int m = 0; m < num; m++ ) {
      std::cout << " The ind. variable " << m << " is: " << S_vars[m] << std::endl;
      std::cout << " The dep. variable f is: " << S_fin << std::endl;
    }
  }

  S_fin = 1.0;
  std::cout << "Re-compute the dependent function." << S_fin << std::endl;
  for ( int n = 0; n < num; n++ ) {
    S_fin += S_vars[n] * S_vars[n];
    std::cout << "f after summing in " << n+1 << " variables: " << S_fin << std::endl;
  }

  return 0;
}
