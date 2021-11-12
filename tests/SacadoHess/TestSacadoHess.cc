#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "Sacado.hpp"
#define DType double
#define NTEST 3
#define NVAR  4

// Test 1:
//   f(x,y,z,w) = 0.5 * ( x^2 + 2 * y^2 + 3 * z^2 + 4 * w^2 );
//   (x,y,z,w) = (0.0, 1.0, 2.0, 3.0);
// Test 2:
//   g(x,y,z,w) = exp( x * y * z * w );
//   (x,y,z,w) = (0.0, 2.0, 4.0, 6.0);
// Test 3:
//   h(x,y,z,w) = log( x + y ) + z / w^1.2;
//   (x,y,z,w) = (0.0, 1.0, 2.0, 3.0);
//   Hess(h) = [-1 -1 0 0; -1 -1 0 0; 0 0 0 -0.1070; 0 0 -0.1070 0.1570]
// Both evaluated at:

int main(int argc, char** argv) {

  // Independent variables
  std::vector<Sacado::Fad::DFad<Sacado::Fad::DFad<DType>>> S_vars;
  for ( int n = 0; n < NVAR; n++ ) {
    S_vars.push_back( Sacado::Fad::DFad<Sacado::Fad::DFad<DType>>( NVAR, n, Sacado::Fad::DFad<DType>( NVAR, n, double(n) ) ) );
  }
  std::vector<std::vector<Sacado::Fad::DFad<Sacado::Fad::DFad<DType>>>> S_AllVars;
  for ( int t = 0; t < NTEST; t++ ) {
    S_AllVars.push_back( S_vars );
    if ( t < 2 ) {
      for ( int n = 0; n < NVAR; n++ )
        S_AllVars[t][n] *= (double)(t+1);
    }
  }

  // Dependent variables
  std::vector<Sacado::Fad::DFad<Sacado::Fad::DFad<DType>>> S_vals;
  for ( int t = 0; t < NTEST; t++ ) {
    Sacado::Fad::DFad<Sacado::Fad::DFad<DType>> s_var = 0.0;
    S_vals.push_back( s_var );
  }

  for ( int n = 0; n < NVAR; n++ ) 
    S_vals[0] += 0.5 * ( (double)n + 1.0 ) * S_AllVars[0][n] * S_AllVars[0][n];

  S_vals[1] = 1.0;
  for ( int n = 0; n < NVAR; n++ ) 
    S_vals[1] *= S_AllVars[1][n];
  S_vals[1] = std::exp( S_vals[1] ); 

  S_vals[2] = std::log( S_AllVars[2][0] + S_AllVars[2][1] ) + S_AllVars[2][2] / std::pow( S_AllVars[2][3], 1.2 );

  // Display the results
  for ( int t = 0; t < NTEST; t++ ) {
    std::cout << "Test " << t << ": " << std::endl;
    for ( int n1 = 0; n1 < NVAR; n1++ ) {
      for ( int n2 = 0; n2 < NVAR; n2++ ) {
        std::cout << "\t" << S_vals[t].dx(n1).dx(n2);
      }
      std::cout << std::endl;
    }
  }

  return 0;
}
