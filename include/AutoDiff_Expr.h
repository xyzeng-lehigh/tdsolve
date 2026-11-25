#ifndef AUTODIFF_EXPR_H
#define AUTODIFF_EXPR_H

// Forward-mode automatic differentiation: expression templates
//
// Every arithmetic operation or math function applied to AutoDiff::Fad::DFad
// objects returns a light-weight "expression" object instead of a DFad. The
// expression is only turned into numbers when it is assigned to (or used to
// construct) a DFad, at which point one loop over the derivative components
// evaluates the whole right hand side. This avoids allocating a temporary
// derivative array for every intermediate result.
//
// There are three kinds of expressions:
//   DFad<T>               - the leaf: a value and an array of derivatives
//                           (AutoDiff_DFad.h)
//   UnaryExpr<Op,E>       - one sub-expression: f(e), e+c, c*e, pow(e,c), ...
//   BinaryExpr<Op,E1,E2>  - two sub-expressions: e1+e2, e1*e2, pow(e1,e2), ...
// The "Op" classes hold the actual formulas and live in AutoDiff_Ops.h.
//
// All expressions provide the same interface:
//   val()            value of the expression
//   size()           number of derivative components (0: a constant)
//   prepare()        must be called once before dx()/fastAccessDx()
//   hasFastAccess()  true if every leaf has derivative components
//   dx(i)            i-th derivative component, valid for any leaves
//   fastAccessDx(i)  i-th derivative component, valid if hasFastAccess()
//
// The value of an expression (and everything in the derivative formula that
// only depends on values, e.g. exp(x) in d/dx exp(x)) is computed once when
// the expression is built, so that the loop over the derivative components
// only does additions, multiplications and divisions.
//
// WARNING: expressions keep references to the DFad objects they are built
// from. Never store one in an "auto" variable; always assign to a DFad.

#include <type_traits>

// The functions that are called for every derivative component (val, dx,
// fastAccessDx) must be inlined for the derivative loop to be efficient.
// Compilers do not always do so on their own inside large functions.
#if defined(__GNUC__) || defined(__clang__)
#define AUTODIFF_INLINE inline __attribute__((always_inline))
#else
#define AUTODIFF_INLINE inline
#endif

namespace AutoDiff {

  namespace Impl {
    // Error handling, kept out of line (src/AutoDiff.cc)
    [[noreturn]] void throw_size_mismatch(const int sz1, const int sz2);
    [[noreturn]] void throw_bad_index(const int i, const int sz);
  }

  namespace Fad {

    template <typename T> class DFad;

    // Base class of all expressions (curiously recurring template pattern)
    template <typename Derived>
    class Expr {
      public:
        const Derived& derived(void) const { return static_cast<const Derived&>(*this); }
    };

    // How an expression stores a sub-expression:
    //   leaves (DFad) by reference - no copy of the derivative array
    //   other expressions by value - they are small and often temporaries
    template <typename E> struct ExprHold { typedef const E type; };
    template <typename T> struct ExprHold< DFad<T> > { typedef const DFad<T>& type; };

    // Innermost scalar type, e.g. double for DFad< DFad<double> >
    template <typename T> struct ScalarType { typedef T type; };
    template <typename T> struct ScalarType< DFad<T> > { typedef typename ScalarType<T>::type type; };

    // An expression of one sub-expression.
    // Op is constructed from the value of the sub-expression (and from the
    // constant c for operations like e+c), and provides
    //   val()   the value
    //   dx(d)   the derivative, given the derivative d of the sub-expression
    template <typename Op, typename E>
    class UnaryExpr : public Expr< UnaryExpr<Op,E> > {
      public:
        typedef typename E::value_type  value_type;
        typedef typename E::scalar_type scalar_type;

        explicit UnaryExpr(const E& arg) :
          arg_(arg), op_(arg.val()) {}
        UnaryExpr(const E& arg, const value_type& c) :
          arg_(arg), op_(arg.val(), c) {}

        AUTODIFF_INLINE const value_type& val(void) const { return op_.val(); }
        int  size(void) const { return arg_.size(); }
        int  prepare(void) const { return arg_.prepare(); }
        bool hasFastAccess(void) const { return arg_.hasFastAccess(); }

        AUTODIFF_INLINE value_type dx(const int i) const { return op_.dx(arg_.dx(i)); }
        AUTODIFF_INLINE value_type fastAccessDx(const int i) const { return op_.dx(arg_.fastAccessDx(i)); }

      private:
        typename ExprHold<E>::type arg_;
        const Op op_;
    };

    // An expression of two sub-expressions.
    // Op is constructed from the two values, and provides
    //   val()        the value
    //   dx(d1,d2)    the derivative if both sub-expressions have derivatives
    //   dx1(d1)      the derivative if only the first one has derivatives
    //   dx2(d2)      the derivative if only the second one has derivatives
    // The last two matter because DFad objects that hold constants (size 0)
    // are freely mixed with those that hold derivatives.
    template <typename Op, typename E1, typename E2>
    class BinaryExpr : public Expr< BinaryExpr<Op,E1,E2> > {
      public:
        typedef typename E1::value_type  value_type;
        typedef typename E1::scalar_type scalar_type;

        BinaryExpr(const E1& arg1, const E2& arg2) :
          arg1_(arg1), arg2_(arg2), op_(arg1.val(), arg2.val()), mode_(ACTIVE_2) {}

        AUTODIFF_INLINE const value_type& val(void) const { return op_.val(); }

        int size(void) const {
          const int sz1 = arg1_.size(), sz2 = arg2_.size();
          return sz1 > sz2 ? sz1 : sz2;
        }

        // Decide which of the two sub-expressions carry derivatives. This is
        // done right before the derivative loop rather than at construction,
        // because the left hand side of an assignment may also appear on the
        // right hand side and may be resized by the assignment.
        int prepare(void) const {
          const int sz1 = arg1_.prepare(), sz2 = arg2_.prepare();
          if ( sz1 > 0 && sz2 > 0 ) {
            if ( sz1 != sz2 ) Impl::throw_size_mismatch(sz1, sz2);
            mode_ = ACTIVE_BOTH;
          }
          else if ( sz1 > 0 )
            mode_ = ACTIVE_1;
          else
            mode_ = ACTIVE_2;
          return sz1 > sz2 ? sz1 : sz2;
        }

        bool hasFastAccess(void) const {
          return arg1_.hasFastAccess() && arg2_.hasFastAccess();
        }

        AUTODIFF_INLINE value_type dx(const int i) const {
          if ( mode_ == ACTIVE_BOTH )
            return op_.dx(arg1_.dx(i), arg2_.dx(i));
          else if ( mode_ == ACTIVE_1 )
            return op_.dx1(arg1_.dx(i));
          else
            return op_.dx2(arg2_.dx(i));
        }

        AUTODIFF_INLINE value_type fastAccessDx(const int i) const {
          return op_.dx(arg1_.fastAccessDx(i), arg2_.fastAccessDx(i));
        }

      private:
        enum Mode { ACTIVE_BOTH = 0, ACTIVE_1, ACTIVE_2 };

        typename ExprHold<E1>::type arg1_;
        typename ExprHold<E2>::type arg2_;
        const Op op_;
        mutable Mode mode_;
    };

  } // namespace Fad

} // namespace AutoDiff

#endif
