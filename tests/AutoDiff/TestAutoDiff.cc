#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <type_traits>
#include "AutoDiff.h"
#define DType double

// Verification of the built-in automatic differentiation (include/AutoDiff.h)
//
// Usage: ./TestAutoDiff [-v]
//   -v  print every check, not only the failed ones
//
// Part 1: storage and bookkeeping of AutoDiff::Fad::DFad
// Part 2: every operator and math function, against finite differences and,
//         where a simple closed form exists, against the exact derivative
// Part 3: expressions of the kind used in TDSolveSacado
// Part 4: +=, -=, *=, /= and assignments where the left hand side also
//         appears on the right hand side
// Part 5: second derivatives with DFad< DFad<DType> >
//         (the three cases of tests/SacadoHess)
// Part 6: error handling
//
// The program returns 0 if all checks pass and 1 otherwise.

typedef AutoDiff::Fad::DFad<DType> ADVar;
typedef AutoDiff::Fad::DFad<ADVar> ADVar2;

// --------------------------- Bookkeeping ---------------------------

static int  num_pass = 0;
static int  num_fail = 0;
static bool verbose  = false;

void check(const bool ok, const std::string& what) {
  if ( ok ) {
    num_pass++;
    if ( verbose ) std::cout << "  pass: " << what << std::endl;
  }
  else {
    num_fail++;
    std::cout << "  FAIL: " << what << std::endl;
  }
  return;
}

bool is_close(const DType a, const DType b, const DType tol) {
  return std::fabs(a-b) <= tol * ( 1.0 + std::fabs(a) + std::fabs(b) );
}

void check_close(const DType a, const DType b, const DType tol, const std::string& what) {
  std::ostringstream os;
  os.precision(17);
  os << what << " (" << a << " vs " << b << ")";
  check( is_close(a, b, tol), os.str() );
  return;
}

// A function of several variables that can be evaluated with DType as well
// as with ADVar: FUN( x[0]*std::exp(x[1]) ) is a function of the vector x.
// TYPE_OF(x) is DType or ADVar accordingly
#define TYPE_OF(X) typename std::decay<decltype(X)>::type::value_type
#define FUN(EXPR) [=](const auto& x) -> TYPE_OF(x) { return EXPR; }

// Check the value and all first derivatives of f at x0:
//   value       against the evaluation with DType
//   derivatives against central finite differences
//   derivatives against grad, if given
template <typename Fun>
void check_fun(const std::string& name, Fun f, const std::vector<DType>& x0,
               const std::vector<DType>& grad = std::vector<DType>()) {
  const int num = x0.size();
  std::vector<DType> x_d(x0);
  std::vector<ADVar> x_a;
  for ( int i = 0; i < num; i++ )
    x_a.push_back( ADVar(num, i, x0[i]) );

  const DType f_d = f(x_d);
  const ADVar f_a = f(x_a);
  check_close( f_a.val(), f_d, 1.0e-14, name + ": value" );
  check( f_a.size() == num, name + ": number of derivative components" );

  for ( int i = 0; i < num; i++ ) {
    const DType h = 1.0e-6 * ( 1.0 + std::fabs(x0[i]) );
    x_d[i] = x0[i] + h;
    const DType f_p = f(x_d);
    x_d[i] = x0[i] - h;
    const DType f_m = f(x_d);
    x_d[i] = x0[i];
    check_close( f_a.dx(i), (f_p-f_m)/(2.0*h), 1.0e-7,
        name + ": d/dx" + std::to_string(i) + " against finite difference" );
    if ( !grad.empty() )
      check_close( f_a.dx(i), grad[i], 1.0e-14,
          name + ": d/dx" + std::to_string(i) + " against exact derivative" );
  }
  return;
}

// f and g must give the same value and the same derivatives
void check_same(const ADVar& f, const ADVar& g, const std::string& name) {
  check_close( f.val(), g.val(), 1.0e-14, name + ": value" );
  check( f.size() == g.size(), name + ": number of derivative components" );
  if ( f.size() == g.size() )
    for ( int i = 0; i < f.size(); i++ )
      check_close( f.dx(i), g.dx(i), 1.0e-14, name + ": d/dx" + std::to_string(i) );
  return;
}

// ------------------ Part 1: storage and bookkeeping -----------------

void test_storage(void) {
  std::cout << "Part 1: storage and bookkeeping" << std::endl;
  const int num = 4;

  // Constants
  ADVar c_def;
  ADVar c_val(2.5);
  ADVar c_int = 3;
  check( c_def.val() == 0.0 && c_def.size() == 0, "default constructor: zero constant" );
  check( c_val.val() == 2.5 && c_val.size() == 0, "constant from a DType" );
  check( c_int.val() == 3.0 && c_int.size() == 0, "constant from an int" );
  check( c_val.dx(0) == 0.0 && c_val.dx(num-1) == 0.0, "derivative of a constant is zero" );
  check( c_val.isPassive() && !c_val.hasFastAccess(), "a constant is passive" );

  // Independent variables
  std::vector<ADVar> x;
  for ( int i = 0; i < num; i++ )
    x.push_back( ADVar(num, i, 1.0+i) );
  bool ok = true;
  for ( int i = 0; i < num; i++ ) {
    ok = ok && x[i].val() == 1.0+i && x[i].size() == num;
    for ( int j = 0; j < num; j++ )
      ok = ok && x[i].dx(j) == ( i == j ? 1.0 : 0.0 ) && x[i].fastAccessDx(j) == x[i].dx(j);
  }
  check( ok, "independent variables: dx(i) is one for itself and zero otherwise" );
  check( !x[0].isPassive() && x[0].hasFastAccess(), "an independent variable is not passive" );

  // Value with zero derivatives
  ADVar z(num, 1.5);
  ok = z.val() == 1.5 && z.size() == num;
  for ( int j = 0; j < num; j++ ) ok = ok && z.dx(j) == 0.0;
  check( ok, "constructor with a size and a value: all derivatives zero" );
  ADVar z0(0, 1.5);
  check( z0.val() == 1.5 && z0.size() == 0, "constructor with size zero: a constant" );

  // Copy, assignment, move
  ADVar f = x[0] * x[1] + x[2];
  ADVar f_copy(f);
  ADVar f_assign;
  f_assign = f;
  check_same( f_copy,   f, "copy constructor" );
  check_same( f_assign, f, "copy assignment" );
  ADVar f_move( std::move(f_copy) );
  check_same( f_move, f, "move constructor" );
  f_assign = 1.0;
  f_assign = std::move(f_move);
  check_same( f_assign, f, "move assignment" );
  f_assign = f_assign;
  check_same( f_assign, f, "self assignment" );

  // Assigning a constant drops the derivatives; the variable can be reused
  f = 2.0;
  check( f.val() == 2.0 && f.size() == 0 && f.dx(1) == 0.0, "assigning a DType makes a constant" );
  f = 7;
  check( f.val() == 7.0 && f.size() == 0, "assigning an int makes a constant" );
  f = x[3];
  check( f.val() == 4.0 && f.size() == num && f.dx(3) == 1.0 && f.dx(0) == 0.0, "reuse after assigning a constant" );
  f = c_val;
  check( f.val() == 2.5 && f.size() == 0, "assigning a constant variable makes a constant" );

  // Re-assigning independent variables (as done at every Newton iteration)
  for ( int i = 0; i < num; i++ )
    x[i] = ADVar(num, i, -1.0-i);
  ok = true;
  for ( int i = 0; i < num; i++ )
    for ( int j = 0; j < num; j++ )
      ok = ok && x[i].val() == -1.0-i && x[i].dx(j) == ( i == j ? 1.0 : 0.0 );
  check( ok, "re-assigning independent variables" );

  // Changing a value keeps the derivatives
  f = x[0] * x[1];
  f.val() = 10.0;
  check( f.val() == 10.0 && f.dx(0) == x[1].val() && f.dx(1) == x[0].val(), "val() can be modified" );
  f.fastAccessDx(2) = 5.0;
  check( f.dx(2) == 5.0, "fastAccessDx() can be modified" );

  // diff, zero, resize, expand
  ADVar g(3.0);
  g.diff(2, num);
  ok = g.val() == 3.0 && g.size() == num;
  for ( int j = 0; j < num; j++ ) ok = ok && g.dx(j) == ( j == 2 ? 1.0 : 0.0 );
  check( ok, "diff() makes an independent variable" );
  g.zero();
  ok = g.size() == num;
  for ( int j = 0; j < num; j++ ) ok = ok && g.dx(j) == 0.0;
  check( ok, "zero() sets all derivatives to zero" );
  g.diff(1, num);
  g.expand(num+2);
  ok = g.size() == num+2;
  for ( int j = 0; j < num+2; j++ ) ok = ok && g.dx(j) == ( j == 1 ? 1.0 : 0.0 );
  check( ok, "expand() keeps the derivatives and adds zeros" );
  g.resize(2);
  check( g.size() == 2, "resize() changes the number of derivatives" );
  g.resizeAndZero(num);
  check( g.size() == num && g.dx(2) == 0.0 && g.dx(3) == 0.0, "resizeAndZero() zeros the added derivatives" );
  g.resize(0);
  check( g.size() == 0 && g.dx(0) == 0.0 && g.val() == 3.0, "resize(0) makes a constant" );

  // Mixing constants and variables
  ADVar m = c_val * x[0] + c_int;
  check( m.size() == num && m.val() == 2.5*x[0].val() + 3.0 && m.dx(0) == 2.5 && m.dx(1) == 0.0,
      "constant variables mixed with independent variables" );
  ADVar cc = c_val * c_int - 1.0;
  check( cc.size() == 0 && cc.val() == 6.5, "an expression of constants is a constant" );
  ADVar zx = z + x[1];
  check( zx.size() == num && zx.dx(1) == 1.0 && zx.dx(0) == 0.0, "variable with zero derivatives plus a variable" );

  // Comparisons act on values
  ADVar a(num, 0, 1.0), b(num, 1, 2.0);
  check(  (a <  b) &&  (a <= b) && !(a >  b) && !(a >= b) && !(a == b) &&  (a != b), "comparison of two variables" );
  check(  (a <  1.5) && (1.5 >  a) && (a == 1.0) && (1.0 == a) && (a != 2) && (a >= 1) && (a <= 1.0), "comparison with a DType or an int" );
  check(  (a*b > a) && (a + b == 3.0) && (2.0 <= a*b) && !(a - b >= 0.0), "comparison of expressions" );
  check(  !(a - 1.0) && !!(a), "operator! tests the value against zero" );

  // Output
  std::ostringstream os;
  os << a << "|" << c_val << "|" << a + b;
  check( os.str() == "1 [ 1 0 0 0 ]|2.5 [ ]|3 [ 1 1 0 0 ]", "output as: value [ derivatives ]" );

  // Single precision
  AutoDiff::Fad::DFad<float> s(2, 0, 1.5f);
  AutoDiff::Fad::DFad<float> t = s * 2.0 + std::sin(s) / 3;
  check( t.size() == 2 && std::fabs(t.dx(0) - (2.0f + std::cos(1.5f)/3.0f)) < 1.0e-6f && t.dx(1) == 0.0f, "DFad<float>" );

  return;
}

// --------------- Part 2: operators and math functions ---------------

void test_functions(void) {
  std::cout << "Part 2: operators and math functions" << std::endl;
  const DType a = 0.7, b = -1.3, c = 2.1;

  // Arithmetic between variables
  check_fun( "x+y", FUN( x[0]+x[1] ), {a,b}, {1.0,1.0} );
  check_fun( "x-y", FUN( x[0]-x[1] ), {a,b}, {1.0,-1.0} );
  check_fun( "x*y", FUN( x[0]*x[1] ), {a,b}, {b,a} );
  check_fun( "x/y", FUN( x[0]/x[1] ), {a,b}, {1.0/b,-a/(b*b)} );
  check_fun( "-x",  FUN( -x[0] ), {a}, {-1.0} );
  check_fun( "+x",  FUN( +x[0] ), {a}, {1.0} );

  // Arithmetic with a constant on either side
  check_fun( "x+c", FUN( x[0]+2.1 ), {a}, {1.0} );
  check_fun( "c+x", FUN( 2.1+x[0] ), {a}, {1.0} );
  check_fun( "x-c", FUN( x[0]-2.1 ), {a}, {1.0} );
  check_fun( "c-x", FUN( 2.1-x[0] ), {a}, {-1.0} );
  check_fun( "x*c", FUN( x[0]*2.1 ), {a}, {c} );
  check_fun( "c*x", FUN( 2.1*x[0] ), {a}, {c} );
  check_fun( "x/c", FUN( x[0]/2.1 ), {a}, {1.0/c} );
  check_fun( "c/x", FUN( 2.1/x[0] ), {a}, {-c/(a*a)} );
  check_fun( "integer constants", FUN( 2*x[0] + x[1]/2 - 1 + 3/x[0] ), {a,b}, {2.0-3.0/(a*a),0.5} );

  // Compound expressions
  check_fun( "x*y*z", FUN( x[0]*x[1]*x[2] ), {a,b,c}, {b*c,a*c,a*b} );
  check_fun( "rational", FUN( (x[0]+x[1])*(x[2]-x[0])/(x[1]*x[2]+3.0) ), {a,b,c} );
  check_fun( "repeated variable", FUN( x[0]*x[0]*x[0] - 2.0*x[0]*x[0] + x[0] ), {a}, {3.0*a*a-4.0*a+1.0} );
  check_fun( "nested", FUN( -(x[0]-x[1]*(x[2]+x[0]/(x[1]-4.0))) ), {a,b,c} );

  // Functions of one variable
  check_fun( "exp",   FUN( std::exp(x[0]) ),   {a}, {std::exp(a)} );
  check_fun( "log",   FUN( std::log(x[0]) ),   {a}, {1.0/a} );
  check_fun( "log10", FUN( std::log10(x[0]) ), {a}, {1.0/(a*std::log(10.0))} );
  check_fun( "sqrt",  FUN( std::sqrt(x[0]) ),  {a}, {0.5/std::sqrt(a)} );
  check_fun( "cbrt",  FUN( std::cbrt(x[0]) ),  {a}, {1.0/(3.0*std::cbrt(a*a))} );
  check_fun( "sin",   FUN( std::sin(x[0]) ),   {a}, {std::cos(a)} );
  check_fun( "cos",   FUN( std::cos(x[0]) ),   {a}, {-std::sin(a)} );
  check_fun( "tan",   FUN( std::tan(x[0]) ),   {a}, {1.0/(std::cos(a)*std::cos(a))} );
  check_fun( "asin",  FUN( std::asin(x[0]) ),  {a}, {1.0/std::sqrt(1.0-a*a)} );
  check_fun( "acos",  FUN( std::acos(x[0]) ),  {a}, {-1.0/std::sqrt(1.0-a*a)} );
  check_fun( "atan",  FUN( std::atan(x[0]) ),  {b}, {1.0/(1.0+b*b)} );
  check_fun( "sinh",  FUN( std::sinh(x[0]) ),  {b}, {std::cosh(b)} );
  check_fun( "cosh",  FUN( std::cosh(x[0]) ),  {b}, {std::sinh(b)} );
  check_fun( "tanh",  FUN( std::tanh(x[0]) ),  {b}, {1.0/(std::cosh(b)*std::cosh(b))} );
  check_fun( "asinh", FUN( std::asinh(x[0]) ), {b}, {1.0/std::sqrt(1.0+b*b)} );
  check_fun( "acosh", FUN( std::acosh(x[0]) ), {c}, {1.0/std::sqrt(c*c-1.0)} );
  check_fun( "atanh", FUN( std::atanh(x[0]) ), {a}, {1.0/(1.0-a*a)} );
  check_fun( "fabs, positive", FUN( std::fabs(x[0]) ), {a}, {1.0} );
  check_fun( "fabs, negative", FUN( std::fabs(x[0]) ), {b}, {-1.0} );
  check_fun( "abs, positive",  FUN( std::abs(x[0]) ),  {a}, {1.0} );
  check_fun( "abs, negative",  FUN( std::abs(x[0]) ),  {b}, {-1.0} );
  check_fun( "function of an expression", FUN( std::exp(std::sin(x[0]*x[1]) - std::sqrt(x[2])) ), {a,b,c} );
  check_fun( "function of functions", FUN( std::log(1.0 + std::fabs(x[0]*x[1])) * std::tanh(x[2]-x[0]) ), {a,b,c} );

  // pow
  check_fun( "pow(x,2)",    FUN( std::pow(x[0],2) ),    {b}, {2.0*b} );
  check_fun( "pow(x,3)",    FUN( std::pow(x[0],3) ),    {b}, {3.0*b*b} );
  check_fun( "pow(x,1)",    FUN( std::pow(x[0],1) ),    {b}, {1.0} );
  check_fun( "pow(x,0)",    FUN( std::pow(x[0],0) ),    {b}, {0.0} );
  check_fun( "pow(x,-2)",   FUN( std::pow(x[0],-2) ),   {b}, {-2.0/(b*b*b)} );
  check_fun( "pow(x,1.2)",  FUN( std::pow(x[0],1.2) ),  {a}, {1.2*std::pow(a,0.2)} );
  check_fun( "pow(x,-0.5)", FUN( std::pow(x[0],-0.5) ), {a}, {-0.5*std::pow(a,-1.5)} );
  check_fun( "pow(0,2)",    FUN( std::pow(x[0],2) ),    {0.0}, {0.0} );
  check_fun( "pow(0,1)",    FUN( std::pow(x[0],1) ),    {0.0}, {1.0} );
  check_fun( "pow(c,x)",    FUN( std::pow(2.1,x[0]) ),  {b}, {std::log(c)*std::pow(c,b)} );
  check_fun( "pow(x,y)",    FUN( std::pow(x[0],x[1]) ), {a,b}, {b*std::pow(a,b-1.0),std::log(a)*std::pow(a,b)} );
  check_fun( "pow of expressions", FUN( std::pow(x[0]*x[0]+1.0,x[1]*x[2]) ), {a,b,c} );

  // atan2
  check_fun( "atan2(y,x)", FUN( std::atan2(x[0],x[1]) ), {a,b}, {b/(a*a+b*b),-a/(a*a+b*b)} );
  check_fun( "atan2(y,c)", FUN( std::atan2(x[0],2.1) ),  {a}, {c/(a*a+c*c)} );
  check_fun( "atan2(c,x)", FUN( std::atan2(2.1,x[0]) ),  {a}, {-c/(a*a+c*c)} );

  // max and min
  check_fun( "max(x,y), x<y", FUN( std::max(x[0],x[1]) ), {b,a}, {0.0,1.0} );
  check_fun( "max(x,y), x>y", FUN( std::max(x[0],x[1]) ), {a,b}, {1.0,0.0} );
  check_fun( "min(x,y), x<y", FUN( std::min(x[0],x[1]) ), {b,a}, {1.0,0.0} );
  check_fun( "min(x,y), x>y", FUN( std::min(x[0],x[1]) ), {a,b}, {0.0,1.0} );
  check_fun( "max(x,c), x<c", FUN( std::max(x[0],2.1) ),  {a}, {0.0} );
  check_fun( "max(c,x), x>c", FUN( std::max(-2.1,x[0]) ), {a}, {1.0} );
  check_fun( "min(x,c), x<c", FUN( std::min(x[0],2.1) ),  {a}, {1.0} );
  check_fun( "min(c,x), x>c", FUN( std::min(-2.1,x[0]) ), {a}, {0.0} );
  check_fun( "max of expressions", FUN( std::max(x[0]*x[1],x[1]+x[2]) ), {a,b,c} );
  check_fun( "min of expressions", FUN( std::min(x[0]*x[1],x[1]+x[2]) ), {a,b,c} );
  check_fun( "min(c,expression)",  FUN( std::min(1.0,5.0/(std::fabs(x[1])/std::fabs(x[0])-1.0)) ), {a,c} );

  // Without the std:: prefix
  check_fun( "unqualified calls", FUN( exp(x[0]) + sqrt(x[0]) * pow(x[1],2) - fabs(x[1]) ), {a,b},
      {std::exp(a)+0.5*b*b/std::sqrt(a),2.0*b*std::sqrt(a)+1.0} );

  return;
}

// ----------- Part 3: expressions of the kind used in TDSolve ---------

// Michaelis-Menten term, taking and returning variables by value
template <typename T>
T calc_mm(const DType mm_ceil, const DType mm_half, T s_var) {
  T result = 0.0;
  if ( std::fabs( s_var ) > 1.0 )
    result = mm_ceil * ( 1.0 / ( (mm_half/s_var) + 1.0 ) );
  else if ( s_var != 0.0 )
    result = mm_ceil * ( s_var / (mm_half + s_var) );
  return result;
}

// Smoothed switch
template <typename T>
T calc_block(const DType k, const DType Ps, T s_P, const DType scale) {
  T s_val;
  if ( k > 0.0 )
    s_val = 1.0 - scale/(1.0+std::exp(-2.0*k*(s_P-Ps)));
  else
    s_val = 1.0;
  return s_val;
}

// Source term in the transformed coordinate, where the radius is stored in
// a variable without derivatives
template <typename T>
T calc_srce(const std::vector<T>& x, const DType l_x, const DType radius, const int dim) {
  T s_r = radius;
  T s_srce = 0.0;
  s_srce += std::pow( l_x * s_r, dim-1 ) *
    ( 0.8 * x[0] * ( 1.0 - x[0] - x[1] ) - calc_mm( 1.5, 0.4, x[2] ) * x[1]
      + calc_block( 3.0, 0.2, x[2], 0.9 ) * calc_mm<T>( 0.7, 2.0, x[0]*s_r ) );
  return s_srce;
}

// Upwind flux with a PPM-type limited reconstruction
template <typename T>
T calc_flux(const std::vector<T>& x, const DType l_x, const DType radius, const DType radius_p) {
  T s_r = radius, s_rp = radius_p;
  T s_diff_l = x[1] - x[0];
  T s_diff_r = x[2] - x[1];
  T s_alpha1 = 2.0*std::fabs(s_diff_l)/std::fabs(s_diff_r);
  T s_alpha2 = 2.0*std::fabs(s_diff_r)/std::fabs(s_diff_l);
  T s_lim = std::min( 1.0, std::min( s_alpha1, s_alpha2 ) );
  T loc_vel = - l_x * s_rp * s_r + x[3];
  T s_flx;
  if ( loc_vel >= 0.0 )
    s_flx = loc_vel * ( x[1] + 0.5 * s_lim * s_diff_r );
  else
    s_flx = loc_vel * ( x[2] - 0.5 * s_lim * s_diff_r );
  return s_flx;
}

void test_solver_expressions(void) {
  std::cout << "Part 3: expressions of the kind used in TDSolve" << std::endl;
  const std::vector<DType> x0 = {0.3, 0.45, 1.7, 0.6};

  check_fun( "Michaelis-Menten, small argument", FUN( calc_mm(1.5, 0.4, x[0]) ), {0.3}, {1.5*0.4/(0.7*0.7)} );
  check_fun( "Michaelis-Menten, large argument", FUN( calc_mm(1.5, 0.4, x[0]) ), {1.7}, {1.5*0.4/(2.1*2.1)} );
  check_fun( "Michaelis-Menten of an expression", FUN( calc_mm<TYPE_OF(x)>(1.5, 0.4, x[0]*x[1]) ), {0.3,0.45} );
  check_fun( "smoothed switch", FUN( calc_block(3.0, 0.2, x[0], 0.9) ), {0.3} );
  for ( int dim = 1; dim <= 3; dim++ )
    check_fun( "source term, dim = " + std::to_string(dim), FUN( calc_srce(x, 0.35, 2.2, dim) ), x0 );
  check_fun( "limited flux, positive velocity", FUN( calc_flux(x, 0.35, 2.2, 0.4) ), x0 );
  check_fun( "limited flux, negative velocity", FUN( calc_flux(x, 0.35, 2.2, 4.0) ), x0 );
  check_fun( "limited flux, limiter not active", FUN( calc_flux(x, 0.35, 2.2, 0.4) ), {0.3, 0.45, 0.61, 0.6} );

  // Summing over cells and assembling a Jacobian entry, as for the total
  // number of cells: d(sum_i w_i x_i^2)/dx_j = 2 w_j x_j
  const int num = 50;
  std::vector<ADVar> s_var;
  for ( int i = 0; i < num; i++ )
    s_var.push_back( ADVar(num, i, 0.1*i - 2.0) );
  ADVar tot = 0.0;
  for ( int i = 0; i < num; i++ )
    tot += ( 1.0 + 0.01*i ) * s_var[i] * s_var[i];
  bool ok = tot.size() == num;
  for ( int j = 0; j < num && ok; j++ )
    ok = is_close( tot.dx(j), 2.0 * ( 1.0 + 0.01*j ) * ( 0.1*j - 2.0 ), 1.0e-14 );
  check( ok, "sum over 50 cells" );

  // Only neighboring cells contribute to a flux difference
  std::vector<ADVar> s_flx(num+1);
  for ( int i = 1; i < num; i++ )
    s_flx[i] = 0.5 * ( s_var[i-1] + s_var[i] ) * ( s_var[i] - s_var[i-1] );
  ok = true;
  for ( int i = 1; i < num-1; i++ )
    for ( int j = 0; j < num; j++ ) {
      const DType elem = s_flx[i+1].dx(j) - s_flx[i].dx(j);
      DType exact = 0.0;
      if ( j == i-1 ) exact =  s_var[i-1].val();
      if ( j == i   ) exact = -2.0 * s_var[i].val();
      if ( j == i+1 ) exact =  s_var[i+1].val();
      ok = ok && is_close( elem, exact, 1.0e-14 );
    }
  check( ok, "tridiagonal Jacobian of a flux difference" );

  return;
}

// --------- Part 4: compound assignment and self reference -----------

void test_assignment(void) {
  std::cout << "Part 4: compound assignment and self reference" << std::endl;
  const int num = 3;
  std::vector<ADVar> x;
  for ( int i = 0; i < num; i++ )
    x.push_back( ADVar(num, i, 0.4 + 0.7*i) );
  const ADVar c = 1.6;
  ADVar f, g;

  // With a variable or an expression
  f = x[0];          f += x[1]*x[2];   g = x[0] + x[1]*x[2];     check_same( f, g, "f += expression" );
  f = x[0];          f -= x[1]*x[2];   g = x[0] - x[1]*x[2];     check_same( f, g, "f -= expression" );
  f = x[0];          f *= x[1]+x[2];   g = x[0] * (x[1]+x[2]);   check_same( f, g, "f *= expression" );
  f = x[0];          f /= x[1]+x[2];   g = x[0] / (x[1]+x[2]);   check_same( f, g, "f /= expression" );
  f = x[0]*x[1];     f += x[2];        g = x[0]*x[1] + x[2];     check_same( f, g, "f += variable" );
  f = x[0]*x[1];     f *= x[2];        g = x[0]*x[1]*x[2];       check_same( f, g, "f *= variable" );

  // With a DType
  f = x[0]*x[1];     f += 2.5;         g = x[0]*x[1] + 2.5;      check_same( f, g, "f += DType" );
  f = x[0]*x[1];     f -= 2.5;         g = x[0]*x[1] - 2.5;      check_same( f, g, "f -= DType" );
  f = x[0]*x[1];     f *= 2.5;         g = x[0]*x[1] * 2.5;      check_same( f, g, "f *= DType" );
  f = x[0]*x[1];     f /= 2.5;         g = x[0]*x[1] / 2.5;      check_same( f, g, "f /= DType" );

  // With a constant variable
  f = x[0]*x[1];     f += c;           g = x[0]*x[1] + 1.6;      check_same( f, g, "f += constant variable" );
  f = x[0]*x[1];     f -= c;           g = x[0]*x[1] - 1.6;      check_same( f, g, "f -= constant variable" );
  f = x[0]*x[1];     f *= c;           g = x[0]*x[1] * 1.6;      check_same( f, g, "f *= constant variable" );
  f = x[0]*x[1];     f /= c;           g = x[0]*x[1] / 1.6;      check_same( f, g, "f /= constant variable" );

  // Starting from a constant
  f = 1.6;           f += x[1]*x[2];   g = 1.6 + x[1]*x[2];      check_same( f, g, "constant += expression" );
  f = 1.6;           f -= x[1]*x[2];   g = 1.6 - x[1]*x[2];      check_same( f, g, "constant -= expression" );
  f = 1.6;           f *= x[1]*x[2];   g = 1.6 * x[1]*x[2];      check_same( f, g, "constant *= expression" );
  f = 1.6;           f /= x[1]*x[2];   g = 1.6 / (x[1]*x[2]);    check_same( f, g, "constant /= expression" );

  // The left hand side appears on the right hand side
  const ADVar f0 = x[0]*x[1] + 0.5;
  f = f0;            f = f + x[2];                g = f0 + x[2];              check_same( f, g, "f = f + x" );
  f = f0;            f = f * x[2] + f;            g = f0 * x[2] + f0;         check_same( f, g, "f = f * x + f" );
  f = f0;            f = x[2] - f;                g = x[2] - f0;              check_same( f, g, "f = x - f" );
  f = f0;            f = 2.0 * f;                 g = 2.0 * f0;               check_same( f, g, "f = 2 f" );
  f = f0;            f = std::exp(f) / f;         g = std::exp(f0) / f0;      check_same( f, g, "f = exp(f) / f" );
  f = f0;            f += f;                      g = f0 + f0;                check_same( f, g, "f += f" );
  f = f0;            f -= 0.5 * f;                g = f0 - 0.5 * f0;          check_same( f, g, "f -= f / 2" );
  f = f0;            f *= f;                      g = f0 * f0;                check_same( f, g, "f *= f" );
  f = f0;            f /= f + 1.0;                g = f0 / (f0 + 1.0);        check_same( f, g, "f /= f + 1" );
  f = 1.6;           f = f + x[1]*x[2];           g = 1.6 + x[1]*x[2];        check_same( f, g, "constant: f = f + expression" );
  f = 1.6;           f = f * x[1];                g = 1.6 * x[1];             check_same( f, g, "constant: f = f * x" );
  f = 1.6;           f = x[1] / f;                g = x[1] / 1.6;             check_same( f, g, "constant: f = x / f" );
  f = 1.6;           f = std::max(f, x[2]);       g = std::max(1.6, x[2]);    check_same( f, g, "constant: f = max(f,x)" );

  // Accumulation into an array that starts as constants
  std::vector<ADVar> s_tot(num);
  for ( int i = 0; i < num; i++ ) s_tot[i] = 0.0;
  for ( int rep = 0; rep < 2; rep++ )
    for ( int i = 0; i < num; i++ )
      s_tot[i] = s_tot[i] + 0.25 * x[i] * x[(i+1)%num];
  for ( int i = 0; i < num; i++ ) {
    g = 0.5 * x[i] * x[(i+1)%num];
    check_same( s_tot[i], g, "accumulation into array entry " + std::to_string(i) );
  }

  return;
}

// ------------------- Part 5: second derivatives ---------------------

void test_hessian(void) {
  std::cout << "Part 5: second derivatives" << std::endl;
  const int num = 4;
  const DType pts[3][4] = { {0.0, 1.0, 2.0, 3.0}, {0.0, 2.0, 4.0, 6.0}, {0.0, 1.0, 2.0, 3.0} };
  DType hess[3][4][4];
  DType grad[3][4];
  for ( int t = 0; t < 3; t++ )
    for ( int i = 0; i < num; i++ ) {
      grad[t][i] = 0.0;
      for ( int j = 0; j < num; j++ ) hess[t][i][j] = 0.0;
    }

  // Test 1: f(x,y,z,w) = 0.5 * ( x^2 + 2 * y^2 + 3 * z^2 + 4 * w^2 )
  for ( int i = 0; i < num; i++ ) {
    grad[0][i]    = (i+1.0) * pts[0][i];
    hess[0][i][i] = i+1.0;
  }
  // Test 2: g(x,y,z,w) = exp( x * y * z * w ) at x = 0
  grad[1][0] = 48.0;
  hess[1][0][0] = 2304.0;
  hess[1][0][1] = hess[1][1][0] = 24.0;
  hess[1][0][2] = hess[1][2][0] = 12.0;
  hess[1][0][3] = hess[1][3][0] = 8.0;
  // Test 3: h(x,y,z,w) = log( x + y ) + z / w^1.2
  grad[2][0] = grad[2][1] = 1.0;
  grad[2][2] = std::pow(3.0,-1.2);
  grad[2][3] = -2.4 * std::pow(3.0,-2.2);
  hess[2][0][0] = hess[2][0][1] = hess[2][1][0] = hess[2][1][1] = -1.0;
  hess[2][2][3] = hess[2][3][2] = -1.2 * std::pow(3.0,-2.2);
  hess[2][3][3] = 2.0 * 1.2 * 2.2 * std::pow(3.0,-3.2);

  const char* names[3] = { "quadratic", "exp(x*y*z*w)", "log(x+y)+z/w^1.2" };
  for ( int t = 0; t < 3; t++ ) {
    std::vector<ADVar2> s_var;
    for ( int i = 0; i < num; i++ )
      s_var.push_back( ADVar2( num, i, ADVar( num, i, pts[t][i] ) ) );

    ADVar2 s_val = 0.0;
    if ( t == 0 ) {
      for ( int i = 0; i < num; i++ )
        s_val += 0.5 * ( (double)i + 1.0 ) * s_var[i] * s_var[i];
    }
    else if ( t == 1 ) {
      s_val = 1.0;
      for ( int i = 0; i < num; i++ )
        s_val *= s_var[i];
      s_val = std::exp( s_val );
    }
    else {
      s_val = std::log( s_var[0] + s_var[1] ) + s_var[2] / std::pow( s_var[3], 1.2 );
    }

    bool ok_grad = true, ok_hess = true;
    for ( int i = 0; i < num; i++ ) {
      ok_grad = ok_grad && is_close( s_val.dx(i).val(), grad[t][i], 1.0e-13 )
                        && is_close( s_val.val().dx(i), grad[t][i], 1.0e-13 );
      for ( int j = 0; j < num; j++ )
        ok_hess = ok_hess && is_close( s_val.dx(i).dx(j), hess[t][i][j], 1.0e-13 );
    }
    check( ok_grad, std::string(names[t]) + ": gradient" );
    check( ok_hess, std::string(names[t]) + ": Hessian" );
    if ( verbose ) {
      for ( int i = 0; i < num; i++ ) {
        for ( int j = 0; j < num; j++ )
          std::cout << "\t" << s_val.dx(i).dx(j);
        std::cout << std::endl;
      }
    }
  }

  // d^2/dx^2 of sin(x) * x^3 by differentiating twice
  const DType p = 0.8;
  ADVar2 s_x( 1, 0, ADVar( 1, 0, p ) );
  ADVar2 s_f = std::sin( s_x ) * std::pow( s_x, 3 );
  check_close( s_f.dx(0).dx(0), ( 6.0*p - p*p*p ) * std::sin(p) + 6.0*p*p * std::cos(p), 1.0e-13, "(sin(x) x^3)''" );

  // Scaling a variable by a DType: y = 2 x, (y^3)'' = 48 x
  ADVar2 s_y( 1, 0, ADVar( 1, 0, p ) );
  s_y *= 2.0;
  ADVar2 s_g = s_y * s_y * s_y;
  check_close( s_g.dx(0).dx(0), 48.0*p, 1.0e-13, "((2x)^3)''" );

  return;
}

// ---------------------- Part 6: error handling ----------------------

void test_errors(void) {
  std::cout << "Part 6: error handling" << std::endl;
  ADVar a(3, 0, 1.0), b(4, 1, 2.0), f;
  bool caught;

  caught = false;
  try { f = a * b; } catch ( std::logic_error& e ) { caught = true; }
  check( caught, "variables with different numbers of derivatives cannot be combined" );

  caught = false;
  try { f = a; f += b; } catch ( std::logic_error& e ) { caught = true; }
  check( caught, "variables with different numbers of derivatives cannot be accumulated" );

  // A failed assignment leaves the left hand side untouched, also if it
  // appears on the right hand side
  caught = false;
  f = a;
  try { f = f + b; } catch ( std::logic_error& e ) { caught = true; }
  check( caught, "variables with different numbers of derivatives cannot be combined with the left hand side" );
  check( f.size() == 3 && f.val() == 1.0 && f.dx(0) == 1.0 && f.dx(1) == 0.0, "a failed assignment does not modify the left hand side" );

  caught = false;
  try { ADVar bad(3, 3, 1.0); } catch ( std::logic_error& e ) { caught = true; }
  check( caught, "the index of an independent variable must be smaller than the size" );

  caught = false;
  try { f.diff(-1, 3); } catch ( std::logic_error& e ) { caught = true; }
  check( caught, "the index of an independent variable must not be negative" );

  return;
}

int main(int argc, char** argv) {
  for ( int i = 1; i < argc; i++ )
    if ( strcmp(argv[i], "-v") == 0 ) verbose = true;

  test_storage();
  test_functions();
  test_solver_expressions();
  test_assignment();
  test_hessian();
  test_errors();

  std::cout << std::endl << num_pass << " checks passed, " << num_fail << " failed." << std::endl;
  return num_fail == 0 ? 0 : 1;
}
