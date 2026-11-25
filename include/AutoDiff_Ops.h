#ifndef AUTODIFF_OPS_H
#define AUTODIFF_OPS_H

// Forward-mode automatic differentiation: operators and math functions
//
// This file has two parts.
//
// Part 1 - the "Op" classes. Each one implements the value and the derivative
// of one operation; they are the only place where differentiation rules are
// written down. There are three flavors (see AutoDiff_Expr.h):
//   f(e)                       e.g. ExpOp       built from x = e.val()
//   e (op) c  and  c (op) e    e.g. MulConstOp  built from x = e.val() and c
//   e1 (op) e2                 e.g. MulOp       built from the two values
// Whatever depends on values only is computed in the constructor, so that
// dx() is as cheap as possible.
//
// Part 2 - the overloaded operators and functions that build expressions out
// of the Op classes:
//   arithmetic   + - * /  (and unary + -)
//   comparison   == != < > <= >=   (compare values, return bool)
//   functions    exp log log10 sqrt cbrt sin cos tan asin acos atan
//                sinh cosh tanh asinh acosh atanh abs fabs
//                pow atan2 max min
//
// To add a new function f: write an Op class with val() and dx(d) in part 1,
// add one AUTODIFF_UNARY_FUN line in part 2, and one "using" line in
// AutoDiff.h if it should also be callable as std::f.
//
// The derivative formulas are written in the same form as those of
// Sacado::Fad::DFad (Trilinos), so that the two produce identical results
// to the last bit (provided the compiler does not contract a*b+c into a
// fused multiply-add; see manual/dev/autodiff.txt). Keep this in mind
// before "simplifying" a formula, e.g. d/x_ into d*(1/x_).

#include <cmath>
#include <cstdlib>
#include <type_traits>
#include "AutoDiff_Expr.h"
#include "AutoDiff_DFad.h"

namespace AutoDiff {

  namespace Fad {

    // =================================================================
    // Part 1a: f(e)
    //   x: value of e;  d: derivative of e
    // =================================================================

    // +e
    template <typename T> class PlusOp {
      public:
        explicit PlusOp(const T& x) : f_(x) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d; }
      private:
        T f_;
    };

    // -e
    template <typename T> class MinusOp {
      public:
        explicit MinusOp(const T& x) : f_(-x) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return -d; }
      private:
        T f_;
    };

    // exp(e)' = exp(x) d
    template <typename T> class ExpOp {
      public:
        explicit ExpOp(const T& x) { using std::exp; f_ = exp(x); }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return f_*d; }
      private:
        T f_;
    };

    // log(e)' = d / x
    template <typename T> class LogOp {
      public:
        explicit LogOp(const T& x) : x_(x) { using std::log; f_ = log(x); }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/x_; }
      private:
        T x_, f_;
    };

    // log10(e)' = d / (log(10) x)
    template <typename T> class Log10Op {
      public:
        explicit Log10Op(const T& x) {
          using std::log10; using std::log;
          f_ = log10(x);
          a_ = log(T(10))*x;
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/a_; }
      private:
        T f_, a_;
    };

    // sqrt(e)' = d / (2 sqrt(x))
    template <typename T> class SqrtOp {
      public:
        explicit SqrtOp(const T& x) {
          using std::sqrt;
          f_ = sqrt(x);
          a_ = T(2)*f_;
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/a_; }
      private:
        T f_, a_;
    };

    // cbrt(e)' = d / (3 cbrt(x^2))
    template <typename T> class CbrtOp {
      public:
        explicit CbrtOp(const T& x) {
          using std::cbrt;
          f_ = cbrt(x);
          a_ = T(3)*cbrt(x*x);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/a_; }
      private:
        T f_, a_;
    };

    // sin(e)' = cos(x) d
    template <typename T> class SinOp {
      public:
        explicit SinOp(const T& x) {
          using std::sin; using std::cos;
          f_ = sin(x);
          a_ = cos(x);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d*a_; }
      private:
        T f_, a_;
    };

    // cos(e)' = -sin(x) d
    template <typename T> class CosOp {
      public:
        explicit CosOp(const T& x) {
          using std::sin; using std::cos;
          f_ = cos(x);
          a_ = sin(x);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return -d*a_; }
      private:
        T f_, a_;
    };

    // tan(e)' = (1 + tan(x)^2) d
    template <typename T> class TanOp {
      public:
        explicit TanOp(const T& x) {
          using std::tan;
          f_ = tan(x);
          a_ = T(1) + f_*f_;
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d*a_; }
      private:
        T f_, a_;
    };

    // asin(e)' = d / sqrt(1 - x^2)
    template <typename T> class ASinOp {
      public:
        explicit ASinOp(const T& x) {
          using std::asin; using std::sqrt;
          f_ = asin(x);
          a_ = sqrt(T(1) - x*x);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/a_; }
      private:
        T f_, a_;
    };

    // acos(e)' = -d / sqrt(1 - x^2)
    template <typename T> class ACosOp {
      public:
        explicit ACosOp(const T& x) {
          using std::acos; using std::sqrt;
          f_ = acos(x);
          a_ = sqrt(T(1) - x*x);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return -d/a_; }
      private:
        T f_, a_;
    };

    // atan(e)' = d / (1 + x^2)
    template <typename T> class ATanOp {
      public:
        explicit ATanOp(const T& x) {
          using std::atan;
          f_ = atan(x);
          a_ = T(1) + x*x;
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/a_; }
      private:
        T f_, a_;
    };

    // sinh(e)' = cosh(x) d
    template <typename T> class SinhOp {
      public:
        explicit SinhOp(const T& x) {
          using std::sinh; using std::cosh;
          f_ = sinh(x);
          a_ = cosh(x);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d*a_; }
      private:
        T f_, a_;
    };

    // cosh(e)' = sinh(x) d
    template <typename T> class CoshOp {
      public:
        explicit CoshOp(const T& x) {
          using std::sinh; using std::cosh;
          f_ = cosh(x);
          a_ = sinh(x);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d*a_; }
      private:
        T f_, a_;
    };

    // tanh(e)' = (1 - tanh(x)^2) d
    template <typename T> class TanhOp {
      public:
        explicit TanhOp(const T& x) {
          using std::tanh;
          f_ = tanh(x);
          a_ = T(1) - f_*f_;
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d*a_; }
      private:
        T f_, a_;
    };

    // asinh(e)' = d / sqrt(1 + x^2)
    template <typename T> class ASinhOp {
      public:
        explicit ASinhOp(const T& x) {
          using std::asinh; using std::sqrt;
          f_ = asinh(x);
          a_ = sqrt(T(1) + x*x);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/a_; }
      private:
        T f_, a_;
    };

    // acosh(e)' = d / sqrt((x-1)(x+1))
    template <typename T> class ACoshOp {
      public:
        explicit ACoshOp(const T& x) {
          using std::acosh; using std::sqrt;
          f_ = acosh(x);
          a_ = sqrt((x - T(1))*(x + T(1)));
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/a_; }
      private:
        T f_, a_;
    };

    // atanh(e)' = d / (1 - x^2)
    template <typename T> class ATanhOp {
      public:
        explicit ATanhOp(const T& x) {
          using std::atanh;
          f_ = atanh(x);
          a_ = T(1) - x*x;
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/a_; }
      private:
        T f_, a_;
    };

    // abs(e)' = d if x >= 0, -d otherwise
    template <typename T> class AbsOp {
      public:
        explicit AbsOp(const T& x) : pos_(x >= T(0)) { using std::abs; f_ = abs(x); }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return pos_ ? d : T(-d); }
      private:
        bool pos_;
        T f_;
    };

    // fabs(e)' = d if x >= 0, -d otherwise
    template <typename T> class FAbsOp {
      public:
        explicit FAbsOp(const T& x) : pos_(x >= T(0)) { using std::fabs; f_ = fabs(x); }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return pos_ ? d : T(-d); }
      private:
        bool pos_;
        T f_;
    };

    // =================================================================
    // Part 1b: e (op) c  and  c (op) e, where c is a constant
    //   x: value of e;  d: derivative of e
    // =================================================================

    // e + c
    template <typename T> class AddConstOp {
      public:
        AddConstOp(const T& x, const T& c) : f_(x + c) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d; }
      private:
        T f_;
    };

    // c + e
    template <typename T> class ConstAddOp {
      public:
        ConstAddOp(const T& x, const T& c) : f_(c + x) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d; }
      private:
        T f_;
    };

    // e - c
    template <typename T> class SubConstOp {
      public:
        SubConstOp(const T& x, const T& c) : f_(x - c) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d; }
      private:
        T f_;
    };

    // c - e
    template <typename T> class ConstSubOp {
      public:
        ConstSubOp(const T& x, const T& c) : f_(c - x) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return -d; }
      private:
        T f_;
    };

    // e * c
    template <typename T> class MulConstOp {
      public:
        MulConstOp(const T& x, const T& c) : c_(c), f_(x*c) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d*c_; }
      private:
        T c_, f_;
    };

    // c * e
    template <typename T> class ConstMulOp {
      public:
        ConstMulOp(const T& x, const T& c) : c_(c), f_(c*x) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return c_*d; }
      private:
        T c_, f_;
    };

    // e / c
    template <typename T> class DivConstOp {
      public:
        DivConstOp(const T& x, const T& c) : c_(c), f_(x/c) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return d/c_; }
      private:
        T c_, f_;
    };

    // (c / e)' = -c d / x^2
    template <typename T> class ConstDivOp {
      public:
        ConstDivOp(const T& x, const T& c) : c_(c), f_(c/x), a_(x*x) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return -d*c_/a_; }
      private:
        T c_, f_, a_;
    };

    // pow(e,c)' = c x^(c-1) d, evaluated as c d / x * x^c.
    // Special cases: c = 1 gives d; otherwise x = 0 gives 0.
    template <typename T> class PowConstOp {
      public:
        PowConstOp(const T& x, const T& c) : x_(x), c_(c) {
          using std::pow;
          f_ = pow(x, c);
          if ( c == T(1) )
            kind_ = IDENTITY;
          else if ( x == T(0) )
            kind_ = ZERO;
          else
            kind_ = GENERAL;
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const {
          if ( kind_ == GENERAL )
            return c_*d/x_*f_;
          else if ( kind_ == IDENTITY )
            return d;
          else
            return T(0);
        }
      private:
        enum Kind { GENERAL = 0, IDENTITY, ZERO };
        T x_, c_, f_;
        Kind kind_;
    };

    // pow(c,e)' = log(c) c^x d. Special case: c = 0 gives 0.
    template <typename T> class ConstPowOp {
      public:
        ConstPowOp(const T& x, const T& c) : zero_(c == T(0)), a_(T(0)) {
          using std::pow; using std::log;
          f_ = pow(c, x);
          if ( !zero_ ) a_ = log(c);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return zero_ ? T(0) : T(d*a_*f_); }
      private:
        bool zero_;
        T f_, a_;
    };

    // atan2(e,c)' = c d / (x^2 + c^2)
    template <typename T> class Atan2ConstOp {
      public:
        Atan2ConstOp(const T& x, const T& c) : c_(c), a_(x*x + c*c) {
          using std::atan2;
          f_ = atan2(x, c);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return (c_*d)/a_; }
      private:
        T c_, f_, a_;
    };

    // atan2(c,e)' = -c d / (c^2 + x^2)
    template <typename T> class ConstAtan2Op {
      public:
        ConstAtan2Op(const T& x, const T& c) : c_(c), a_(c*c + x*x) {
          using std::atan2;
          f_ = atan2(c, x);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return (-c_*d)/a_; }
      private:
        T c_, f_, a_;
    };

    // max(e,c): e if x >= c, c otherwise
    template <typename T> class MaxConstOp {
      public:
        MaxConstOp(const T& x, const T& c) : first_(x >= c), f_(first_ ? x : c) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return first_ ? d : T(0); }
      private:
        bool first_;
        T f_;
    };

    // max(c,e): c if c >= x, e otherwise
    template <typename T> class ConstMaxOp {
      public:
        ConstMaxOp(const T& x, const T& c) : first_(c >= x), f_(first_ ? c : x) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return first_ ? T(0) : d; }
      private:
        bool first_;
        T f_;
    };

    // min(e,c): e if x <= c, c otherwise
    template <typename T> class MinConstOp {
      public:
        MinConstOp(const T& x, const T& c) : first_(x <= c), f_(first_ ? x : c) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return first_ ? d : T(0); }
      private:
        bool first_;
        T f_;
    };

    // min(c,e): c if c <= x, e otherwise
    template <typename T> class ConstMinOp {
      public:
        ConstMinOp(const T& x, const T& c) : first_(c <= x), f_(first_ ? c : x) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d) const { return first_ ? T(0) : d; }
      private:
        bool first_;
        T f_;
    };

    // =================================================================
    // Part 1c: e1 (op) e2
    //   x, y: values of e1, e2;  d1, d2: derivatives of e1, e2
    //   dx(d1,d2): both have derivatives
    //   dx1(d1):   only e1 has derivatives (e2 is a constant)
    //   dx2(d2):   only e2 has derivatives (e1 is a constant)
    // =================================================================

    // e1 + e2
    template <typename T> class AddOp {
      public:
        AddOp(const T& x, const T& y) : f_(x + y) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d1, const T& d2) const { return d1 + d2; }
        AUTODIFF_INLINE T dx1(const T& d1) const { return d1; }
        AUTODIFF_INLINE T dx2(const T& d2) const { return d2; }
      private:
        T f_;
    };

    // e1 - e2
    template <typename T> class SubOp {
      public:
        SubOp(const T& x, const T& y) : f_(x - y) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d1, const T& d2) const { return d1 - d2; }
        AUTODIFF_INLINE T dx1(const T& d1) const { return d1; }
        AUTODIFF_INLINE T dx2(const T& d2) const { return -d2; }
      private:
        T f_;
    };

    // (e1 * e2)' = x d2 + d1 y
    template <typename T> class MulOp {
      public:
        MulOp(const T& x, const T& y) : x_(x), y_(y), f_(x*y) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d1, const T& d2) const { return x_*d2 + d1*y_; }
        AUTODIFF_INLINE T dx1(const T& d1) const { return d1*y_; }
        AUTODIFF_INLINE T dx2(const T& d2) const { return x_*d2; }
      private:
        T x_, y_, f_;
    };

    // (e1 / e2)' = (d1 y - d2 x) / y^2
    template <typename T> class DivOp {
      public:
        DivOp(const T& x, const T& y) : x_(x), y_(y), f_(x/y), a_(y*y) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d1, const T& d2) const { return (d1*y_ - d2*x_)/a_; }
        AUTODIFF_INLINE T dx1(const T& d1) const { return d1/y_; }
        AUTODIFF_INLINE T dx2(const T& d2) const { return -d2*x_/a_; }
      private:
        T x_, y_, f_, a_;
    };

    // pow(e1,e2)' = (d2 log(x) + y d1 / x) x^y. Special case: x = 0 gives 0.
    // If only e1 has derivatives, the rule of pow(e,c) is used.
    template <typename T> class PowOp {
      public:
        PowOp(const T& x, const T& y) :
          x_(x), y_(y), zero_(x == T(0)), one_(y == T(1)), has_log_(false), a_(T(0)) {
          using std::pow;
          f_ = pow(x, y);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d1, const T& d2) const {
          return zero_ ? T(0) : T((d2*log_x() + y_*d1/x_)*f_);
        }
        AUTODIFF_INLINE T dx1(const T& d1) const {
          return one_ ? d1 : ( zero_ ? T(0) : T(y_*d1/x_*f_) );
        }
        AUTODIFF_INLINE T dx2(const T& d2) const {
          return zero_ ? T(0) : T(d2*log_x()*f_);
        }
      private:
        // log(x) is only needed (and only valid) if e2 has derivatives
        const T& log_x(void) const {
          if ( !has_log_ ) {
            using std::log;
            a_ = log(x_);
            has_log_ = true;
          }
          return a_;
        }
        T x_, y_, f_;
        bool zero_, one_;
        mutable bool has_log_;
        mutable T a_;
    };

    // atan2(e1,e2)' = (y d1 - x d2) / (x^2 + y^2)
    template <typename T> class Atan2Op {
      public:
        Atan2Op(const T& x, const T& y) : x_(x), y_(y), a_(x*x + y*y) {
          using std::atan2;
          f_ = atan2(x, y);
        }
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d1, const T& d2) const { return (y_*d1 - x_*d2)/a_; }
        AUTODIFF_INLINE T dx1(const T& d1) const { return y_*d1/a_; }
        AUTODIFF_INLINE T dx2(const T& d2) const { return -x_*d2/a_; }
      private:
        T x_, y_, f_, a_;
    };

    // max(e1,e2): e1 if x >= y, e2 otherwise
    template <typename T> class MaxOp {
      public:
        MaxOp(const T& x, const T& y) : first_(x >= y), f_(first_ ? x : y) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d1, const T& d2) const { return first_ ? d1 : d2; }
        AUTODIFF_INLINE T dx1(const T& d1) const { return first_ ? d1 : T(0); }
        AUTODIFF_INLINE T dx2(const T& d2) const { return first_ ? T(0) : d2; }
      private:
        bool first_;
        T f_;
    };

    // min(e1,e2): e1 if x <= y, e2 otherwise
    template <typename T> class MinOp {
      public:
        MinOp(const T& x, const T& y) : first_(x <= y), f_(first_ ? x : y) {}
        AUTODIFF_INLINE const T& val(void) const { return f_; }
        AUTODIFF_INLINE T dx(const T& d1, const T& d2) const { return first_ ? d1 : d2; }
        AUTODIFF_INLINE T dx1(const T& d1) const { return first_ ? d1 : T(0); }
        AUTODIFF_INLINE T dx2(const T& d2) const { return first_ ? T(0) : d2; }
      private:
        bool first_;
        T f_;
    };

    // =================================================================
    // Part 2: operators and functions
    // =================================================================

    // FUN(e)
#define AUTODIFF_UNARY_FUN(FUN,OP)                                           \
    template <typename E>                                                    \
    inline UnaryExpr< OP<typename E::value_type>, E >                        \
    FUN (const Expr<E>& e) {                                                 \
      return UnaryExpr< OP<typename E::value_type>, E >(e.derived());        \
    }

    // FUN(e1,e2), FUN(e,c) and FUN(c,e)
    // In the last two, c is anything convertible to the value type of e.
#define AUTODIFF_BINARY_FUN(FUN,OP,OP_EXPR_CONST,OP_CONST_EXPR)              \
    template <typename E1, typename E2>                                      \
    inline typename std::enable_if<                                          \
      std::is_same<typename E1::value_type,typename E2::value_type>::value,  \
      BinaryExpr< OP<typename E1::value_type>, E1, E2 > >::type              \
    FUN (const Expr<E1>& e1, const Expr<E2>& e2) {                           \
      return BinaryExpr< OP<typename E1::value_type>, E1, E2 >(              \
          e1.derived(), e2.derived());                                       \
    }                                                                        \
    template <typename E>                                                    \
    inline UnaryExpr< OP_EXPR_CONST<typename E::value_type>, E >             \
    FUN (const Expr<E>& e, const typename E::value_type& c) {                \
      return UnaryExpr< OP_EXPR_CONST<typename E::value_type>, E >(          \
          e.derived(), c);                                                   \
    }                                                                        \
    template <typename E>                                                    \
    inline UnaryExpr< OP_CONST_EXPR<typename E::value_type>, E >             \
    FUN (const typename E::value_type& c, const Expr<E>& e) {                \
      return UnaryExpr< OP_CONST_EXPR<typename E::value_type>, E >(          \
          e.derived(), c);                                                   \
    }

    // Comparisons e1 (op) e2, e (op) c and c (op) e: only values are compared
#define AUTODIFF_RELATIONAL_OP(OP)                                           \
    template <typename E1, typename E2>                                      \
    inline typename std::enable_if<                                          \
      std::is_same<typename E1::value_type,typename E2::value_type>::value,  \
      bool >::type                                                           \
    operator OP (const Expr<E1>& e1, const Expr<E2>& e2) {                   \
      return e1.derived().val() OP e2.derived().val();                       \
    }                                                                        \
    template <typename E>                                                    \
    inline bool                                                              \
    operator OP (const Expr<E>& e, const typename E::value_type& c) {        \
      return e.derived().val() OP c;                                         \
    }                                                                        \
    template <typename E>                                                    \
    inline bool                                                              \
    operator OP (const typename E::value_type& c, const Expr<E>& e) {        \
      return c OP e.derived().val();                                         \
    }

    AUTODIFF_UNARY_FUN(operator+, PlusOp)
    AUTODIFF_UNARY_FUN(operator-, MinusOp)
    AUTODIFF_UNARY_FUN(exp,   ExpOp)
    AUTODIFF_UNARY_FUN(log,   LogOp)
    AUTODIFF_UNARY_FUN(log10, Log10Op)
    AUTODIFF_UNARY_FUN(sqrt,  SqrtOp)
    AUTODIFF_UNARY_FUN(cbrt,  CbrtOp)
    AUTODIFF_UNARY_FUN(sin,   SinOp)
    AUTODIFF_UNARY_FUN(cos,   CosOp)
    AUTODIFF_UNARY_FUN(tan,   TanOp)
    AUTODIFF_UNARY_FUN(asin,  ASinOp)
    AUTODIFF_UNARY_FUN(acos,  ACosOp)
    AUTODIFF_UNARY_FUN(atan,  ATanOp)
    AUTODIFF_UNARY_FUN(sinh,  SinhOp)
    AUTODIFF_UNARY_FUN(cosh,  CoshOp)
    AUTODIFF_UNARY_FUN(tanh,  TanhOp)
    AUTODIFF_UNARY_FUN(asinh, ASinhOp)
    AUTODIFF_UNARY_FUN(acosh, ACoshOp)
    AUTODIFF_UNARY_FUN(atanh, ATanhOp)
    AUTODIFF_UNARY_FUN(abs,   AbsOp)
    AUTODIFF_UNARY_FUN(fabs,  FAbsOp)

    AUTODIFF_BINARY_FUN(operator+, AddOp,   AddConstOp,   ConstAddOp)
    AUTODIFF_BINARY_FUN(operator-, SubOp,   SubConstOp,   ConstSubOp)
    AUTODIFF_BINARY_FUN(operator*, MulOp,   MulConstOp,   ConstMulOp)
    AUTODIFF_BINARY_FUN(operator/, DivOp,   DivConstOp,   ConstDivOp)
    AUTODIFF_BINARY_FUN(pow,       PowOp,   PowConstOp,   ConstPowOp)
    AUTODIFF_BINARY_FUN(atan2,     Atan2Op, Atan2ConstOp, ConstAtan2Op)
    AUTODIFF_BINARY_FUN(max,       MaxOp,   MaxConstOp,   ConstMaxOp)
    AUTODIFF_BINARY_FUN(min,       MinOp,   MinConstOp,   ConstMinOp)

    AUTODIFF_RELATIONAL_OP(==)
    AUTODIFF_RELATIONAL_OP(!=)
    AUTODIFF_RELATIONAL_OP(<)
    AUTODIFF_RELATIONAL_OP(>)
    AUTODIFF_RELATIONAL_OP(<=)
    AUTODIFF_RELATIONAL_OP(>=)

#undef AUTODIFF_UNARY_FUN
#undef AUTODIFF_BINARY_FUN
#undef AUTODIFF_RELATIONAL_OP

    // !e is true if the value of e is zero
    template <typename E>
    inline bool operator! (const Expr<E>& e) {
      return !e.derived().val();
    }

  } // namespace Fad

} // namespace AutoDiff

#endif
