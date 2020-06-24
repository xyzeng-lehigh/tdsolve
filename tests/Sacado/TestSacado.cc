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

  for ( int n = 0; n < num; n++ ) {
    Sacado::Fad::DFad<DType> value(num,n,(double)n+1.0);
    S_vars.push_back(value);
  }


  // Dependent variables
  Sacado::Fad::DFad<DType> S_fin(0.0);
  for ( int n = 0; n < num; n++ ) {
    S_fin += S_vars[n] * S_vars[n];
  }

  std::cout << "Variable:   " << std::endl;
  for ( int n = 0; n < num; n++ ) 
    std::cout << "\t" << n << "\t = " << S_vars[n].val() << std::endl;
  std::cout << "Function:   " << S_fin.val() << std::endl;
  std::cout << "Derivatives: " << std::endl;
  for ( int n = 0; n < num; n++ )
    std::cout << "\t" << n << "\t = " << S_fin.dx(n) << std::endl;

  return 0;
}
