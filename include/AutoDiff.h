#ifndef AUTODIFF_H
#define AUTODIFF_H

// Forward-mode automatic differentiation
//
// A self-contained replacement of the Sacado package of Trilinos, for the
// part of it that TDSolve relies on: the forward-mode variable with a
// dynamically sized derivative array, Sacado::Fad::DFad<T>. The usage is the
// same, with the namespace Sacado replaced by AutoDiff:
//
//   #include "AutoDiff.h"
//
//   typedef AutoDiff::Fad::DFad<DType> ADVar;
//
//   // Three independent variables x0, x1, x2 with values 1, 2, 3
//   std::vector<ADVar> x;
//   for ( int i = 0; i < 3; i++ )
//     x.push_back( ADVar(3, i, 1.0+i) );
//
//   // A dependent variable
//   ADVar f = x[0] * std::exp(x[1]) + std::pow(x[2], 2);
//
//   f.val();   // value of f
//   f.dx(1);   // df/dx1
//
// The files are organized as:
//   AutoDiff_Expr.h   expression templates (how an expression is evaluated)
//   AutoDiff_DFad.h   the DFad variable (storage, assignment, +=, ...)
//   AutoDiff_Ops.h    operators and math functions (differentiation rules)
//   AutoDiff.h        this file; the only one that needs to be included
//   src/AutoDiff.cc   error handling
//
// See tests/AutoDiff for examples and verification.

#include "AutoDiff_Expr.h"
#include "AutoDiff_DFad.h"
#include "AutoDiff_Ops.h"

// Make the math functions callable with the std:: prefix, e.g. std::exp(x)
// and std::pow(x,2) where x is a DFad or an expression, exactly as for a
// DType. (Sacado does the same.)
namespace std {
  using AutoDiff::Fad::exp;
  using AutoDiff::Fad::log;
  using AutoDiff::Fad::log10;
  using AutoDiff::Fad::sqrt;
  using AutoDiff::Fad::cbrt;
  using AutoDiff::Fad::sin;
  using AutoDiff::Fad::cos;
  using AutoDiff::Fad::tan;
  using AutoDiff::Fad::asin;
  using AutoDiff::Fad::acos;
  using AutoDiff::Fad::atan;
  using AutoDiff::Fad::sinh;
  using AutoDiff::Fad::cosh;
  using AutoDiff::Fad::tanh;
  using AutoDiff::Fad::asinh;
  using AutoDiff::Fad::acosh;
  using AutoDiff::Fad::atanh;
  using AutoDiff::Fad::abs;
  using AutoDiff::Fad::fabs;
  using AutoDiff::Fad::pow;
  using AutoDiff::Fad::atan2;
  using AutoDiff::Fad::max;
  using AutoDiff::Fad::min;
}

#endif
