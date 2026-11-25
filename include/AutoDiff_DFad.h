#ifndef AUTODIFF_DFAD_H
#define AUTODIFF_DFAD_H

// Forward-mode automatic differentiation: the DFad variable
//
// AutoDiff::Fad::DFad<T> holds a value of type T and a dynamically sized
// array of derivatives with respect to a set of independent variables.
// It is a drop-in replacement of Sacado::Fad::DFad<T> from Trilinos for the
// functionality used in TDSolve:
//
//   DFad<DType> c(3.0);        // constant: no derivative components
//   DFad<DType> x(n, i, 3.0);  // i-th of n independent variables with value 3
//   DFad<DType> f = x * x + c; // dependent variable
//   f.val();                   // value
//   f.dx(j);                   // derivative w.r.t. the j-th independent variable
//   f.size();                  // number of derivative components (0 or n)
//
// A DFad of size 0 is a constant; its dx(j) is zero for all j. Constants and
// variables of size n may be mixed freely in an expression, but all DFad's
// with derivatives must have the same size.
//
// DFad< DFad<T> > can be used to compute second derivatives.

#include <ostream>
#include <type_traits>
#include <utility>
#include "AutoDiff_Expr.h"

namespace AutoDiff {

  namespace Fad {

    template <typename T>
    class DFad : public Expr< DFad<T> > {
      public:
        typedef T value_type;
        typedef typename ScalarType<T>::type scalar_type;

        // ----------------------- Constructors ------------------------

        // Constant with value zero
        DFad(void) : val_(T(0.)), sz_(0), len_(0), dx_(NULL) {}

        // Constant with value x
        DFad(const T& x) : val_(x), sz_(0), len_(0), dx_(NULL) {}

        // Constant with value x, where x is anything convertible to T
        // (e.g. an int or a double for DFad<DType> and DFad< DFad<DType> >)
        template <typename S,
                  typename std::enable_if< std::is_convertible<S,T>::value, int >::type = 0>
        DFad(const S& x) : val_(x), sz_(0), len_(0), dx_(NULL) {}

        // Value x and sz derivative components, all zero
        DFad(const int sz, const T& x) : val_(x), sz_(sz > 0 ? sz : 0), len_(sz_), dx_(NULL) {
          if ( sz_ > 0 ) dx_ = new T[sz_]();
        }

        // The i-th of sz independent variables, with value x:
        // sz derivative components, all zero except the i-th which is one
        DFad(const int sz, const int i, const T& x) : val_(x), sz_(sz), len_(sz), dx_(NULL) {
          if ( i < 0 || i >= sz ) Impl::throw_bad_index(i, sz);
          dx_ = new T[sz_]();
          dx_[i] = T(1.);
        }

        // Copy and move
        DFad(const DFad& x) : val_(x.val_), sz_(x.sz_), len_(x.sz_), dx_(NULL) {
          if ( sz_ > 0 ) {
            dx_ = new T[sz_];
            for ( int i = 0; i < sz_; i++ ) dx_[i] = x.dx_[i];
          }
        }

        DFad(DFad&& x) noexcept :
          val_(std::move(x.val_)), sz_(x.sz_), len_(x.len_), dx_(x.dx_) {
          x.sz_ = 0; x.len_ = 0; x.dx_ = NULL;
        }

        // Evaluate an expression, e.g. DFad<DType> f = x * x + c;
        template <typename E,
                  typename std::enable_if< std::is_same<typename E::value_type,T>::value, int >::type = 0>
        DFad(const Expr<E>& expr) : val_(T(0.)), sz_(0), len_(0), dx_(NULL) {
          const E& x = expr.derived();
          const int xsz = x.prepare();
          if ( xsz > 0 ) {
            dx_ = new T[xsz];
            sz_ = len_ = xsz;
            if ( x.hasFastAccess() )
              for ( int i = 0; i < xsz; i++ ) dx_[i] = x.fastAccessDx(i);
            else
              for ( int i = 0; i < xsz; i++ ) dx_[i] = x.dx(i);
          }
          val_ = x.val();
        }

        ~DFad(void) { delete[] dx_; }

        // --------------------- Value and derivatives ------------------

        AUTODIFF_INLINE const T& val(void) const { return val_; }
        AUTODIFF_INLINE T&       val(void)       { return val_; }

        // Number of derivative components
        int size(void) const { return sz_; }

        // True if there are derivative components, i.e. fastAccessDx is valid
        bool hasFastAccess(void) const { return sz_ != 0; }

        // True if this is a constant
        bool isPassive(void) const { return sz_ == 0; }

        // The i-th derivative component; zero if this is a constant
        AUTODIFF_INLINE T dx(const int i) const { return sz_ ? dx_[i] : T(0.); }

        // The i-th derivative component without checking that it exists
        AUTODIFF_INLINE const T& fastAccessDx(const int i) const { return dx_[i]; }
        AUTODIFF_INLINE T&       fastAccessDx(const int i)       { return dx_[i]; }

        // The array of derivative components (NULL if never allocated)
        const T* dx(void) const { return dx_; }

        // -------------------- Managing the derivatives ----------------

        // Make this the ith of n independent variables, keeping the value
        void diff(const int ith, const int n) {
          if ( ith < 0 || ith >= n ) Impl::throw_bad_index(ith, n);
          if ( sz_ != n ) resize(n);
          zero();
          dx_[ith] = T(1.);
          return;
        }

        // Change the number of derivative components. Existing components are
        // not preserved and new ones are not initialized.
        void resize(const int sz) {
          if ( sz > len_ ) {
            delete[] dx_;
            dx_ = new T[sz];
            len_ = sz;
          }
          sz_ = sz;
          return;
        }

        // Same as resize, but the components that are added are set to zero
        void resizeAndZero(const int sz) {
          if ( sz > len_ ) {
            delete[] dx_;
            dx_ = new T[sz]();
            len_ = sz;
          }
          else if ( sz > sz_ ) {
            for ( int i = sz_; i < sz; i++ ) dx_[i] = T(0.);
          }
          sz_ = sz;
          return;
        }

        // Increase the number of derivative components, preserving the
        // existing ones and setting the new ones to zero
        void expand(const int sz) {
          if ( sz > len_ ) {
            T* dx_new = new T[sz]();
            for ( int i = 0; i < sz_; i++ ) dx_new[i] = dx_[i];
            delete[] dx_;
            dx_ = dx_new;
            len_ = sz;
          }
          else if ( sz > sz_ ) {
            for ( int i = sz_; i < sz; i++ ) dx_[i] = T(0.);
          }
          if ( sz > sz_ ) sz_ = sz;
          return;
        }

        // Set all derivative components to zero
        void zero(void) {
          for ( int i = 0; i < sz_; i++ ) dx_[i] = T(0.);
          return;
        }

        // Expression interface (see AutoDiff_Expr.h)
        int prepare(void) const { return sz_; }

        // ------------------------- Assignment --------------------------

        // Assign a constant: the derivatives are dropped
        DFad& operator= (const T& v) {
          val_ = v;
          sz_  = 0;
          return *this;
        }

        template <typename S>
        typename std::enable_if< std::is_convertible<S,T>::value, DFad& >::type
        operator= (const S& v) {
          val_ = T(v);
          sz_  = 0;
          return *this;
        }

        DFad& operator= (const DFad& x) {
          if ( this != &x ) {
            if ( sz_ != x.sz_ ) resize(x.sz_);
            for ( int i = 0; i < sz_; i++ ) dx_[i] = x.dx_[i];
            val_ = x.val_;
          }
          return *this;
        }

        DFad& operator= (DFad&& x) noexcept {
          if ( this != &x ) {
            delete[] dx_;
            val_ = std::move(x.val_);
            sz_  = x.sz_;
            len_ = x.len_;
            dx_  = x.dx_;
            x.sz_ = 0; x.len_ = 0; x.dx_ = NULL;
          }
          return *this;
        }

        // Evaluate an expression. The expression may contain this DFad
        // itself, e.g. f = f + x * x;
        template <typename E>
        typename std::enable_if< std::is_same<typename E::value_type,T>::value, DFad& >::type
        operator= (const Expr<E>& expr) {
          const E& x = expr.derived();
          const int xsz = x.prepare();
          if ( xsz != sz_ ) {
            // If this DFad appears in the expression it is a constant there
            // (otherwise the sizes would agree); it becomes a variable with
            // zero derivatives, and the expression has to be prepared again.
            resizeAndZero(xsz);
            x.prepare();
          }
          if ( xsz > 0 ) {
            if ( x.hasFastAccess() )
              for ( int i = 0; i < xsz; i++ ) dx_[i] = x.fastAccessDx(i);
            else
              for ( int i = 0; i < xsz; i++ ) dx_[i] = x.dx(i);
          }
          val_ = x.val();
          return *this;
        }

        // ------------------ Arithmetic with a constant -----------------

        DFad& operator+= (const T& v) {
          val_ += v;
          return *this;
        }

        DFad& operator-= (const T& v) {
          val_ -= v;
          return *this;
        }

        DFad& operator*= (const T& v) {
          for ( int i = 0; i < sz_; i++ ) dx_[i] *= v;
          val_ *= v;
          return *this;
        }

        DFad& operator/= (const T& v) {
          for ( int i = 0; i < sz_; i++ ) dx_[i] /= v;
          val_ /= v;
          return *this;
        }

        // ---------------- Arithmetic with an expression ----------------

        template <typename E>
        typename std::enable_if< std::is_same<typename E::value_type,T>::value, DFad& >::type
        operator+= (const Expr<E>& expr) {
          const E& x = expr.derived();
          const int xsz = x.prepare(), sz = sz_;
          if ( xsz > 0 ) {
            if ( sz > 0 ) {
              if ( xsz != sz ) Impl::throw_size_mismatch(sz, xsz);
              if ( x.hasFastAccess() )
                for ( int i = 0; i < sz; i++ ) dx_[i] += x.fastAccessDx(i);
              else
                for ( int i = 0; i < sz; i++ ) dx_[i] += x.dx(i);
            }
            else {
              resizeAndZero(xsz);
              x.prepare();
              if ( x.hasFastAccess() )
                for ( int i = 0; i < xsz; i++ ) dx_[i] = x.fastAccessDx(i);
              else
                for ( int i = 0; i < xsz; i++ ) dx_[i] = x.dx(i);
            }
          }
          val_ += x.val();
          return *this;
        }

        template <typename E>
        typename std::enable_if< std::is_same<typename E::value_type,T>::value, DFad& >::type
        operator-= (const Expr<E>& expr) {
          const E& x = expr.derived();
          const int xsz = x.prepare(), sz = sz_;
          if ( xsz > 0 ) {
            if ( sz > 0 ) {
              if ( xsz != sz ) Impl::throw_size_mismatch(sz, xsz);
              if ( x.hasFastAccess() )
                for ( int i = 0; i < sz; i++ ) dx_[i] -= x.fastAccessDx(i);
              else
                for ( int i = 0; i < sz; i++ ) dx_[i] -= x.dx(i);
            }
            else {
              resizeAndZero(xsz);
              x.prepare();
              if ( x.hasFastAccess() )
                for ( int i = 0; i < xsz; i++ ) dx_[i] = -x.fastAccessDx(i);
              else
                for ( int i = 0; i < xsz; i++ ) dx_[i] = -x.dx(i);
            }
          }
          val_ -= x.val();
          return *this;
        }

        template <typename E>
        typename std::enable_if< std::is_same<typename E::value_type,T>::value, DFad& >::type
        operator*= (const Expr<E>& expr) {
          const E& x = expr.derived();
          const int xsz = x.prepare(), sz = sz_;
          const T xval = x.val();
          const T v    = val_;
          if ( xsz > 0 ) {
            if ( sz > 0 ) {
              if ( xsz != sz ) Impl::throw_size_mismatch(sz, xsz);
              if ( x.hasFastAccess() )
                for ( int i = 0; i < sz; i++ ) dx_[i] = v*x.fastAccessDx(i) + dx_[i]*xval;
              else
                for ( int i = 0; i < sz; i++ ) dx_[i] = v*x.dx(i) + dx_[i]*xval;
            }
            else {
              resizeAndZero(xsz);
              x.prepare();
              if ( x.hasFastAccess() )
                for ( int i = 0; i < xsz; i++ ) dx_[i] = v*x.fastAccessDx(i);
              else
                for ( int i = 0; i < xsz; i++ ) dx_[i] = v*x.dx(i);
            }
          }
          else {
            for ( int i = 0; i < sz; i++ ) dx_[i] *= xval;
          }
          val_ *= xval;
          return *this;
        }

        template <typename E>
        typename std::enable_if< std::is_same<typename E::value_type,T>::value, DFad& >::type
        operator/= (const Expr<E>& expr) {
          const E& x = expr.derived();
          const int xsz = x.prepare(), sz = sz_;
          const T xval = x.val();
          const T v    = val_;
          if ( xsz > 0 ) {
            const T xval2 = xval*xval;
            if ( sz > 0 ) {
              if ( xsz != sz ) Impl::throw_size_mismatch(sz, xsz);
              if ( x.hasFastAccess() )
                for ( int i = 0; i < sz; i++ ) dx_[i] = ( dx_[i]*xval - v*x.fastAccessDx(i) ) / xval2;
              else
                for ( int i = 0; i < sz; i++ ) dx_[i] = ( dx_[i]*xval - v*x.dx(i) ) / xval2;
            }
            else {
              resizeAndZero(xsz);
              x.prepare();
              if ( x.hasFastAccess() )
                for ( int i = 0; i < xsz; i++ ) dx_[i] = -v*x.fastAccessDx(i) / xval2;
              else
                for ( int i = 0; i < xsz; i++ ) dx_[i] = -v*x.dx(i) / xval2;
            }
          }
          else {
            for ( int i = 0; i < sz; i++ ) dx_[i] /= xval;
          }
          val_ /= xval;
          return *this;
        }

      private:
        T   val_; // Value
        int sz_;  // Number of derivative components
        int len_; // Allocated length of dx_ (>= sz_)
        T*  dx_;  // Derivative components
    };

    // Output as: value [ dx(0) dx(1) ... ]
    template <typename E>
    std::ostream& operator<< (std::ostream& os, const Expr<E>& expr) {
      const E& x = expr.derived();
      const int sz = x.prepare();
      os << x.val() << " [";
      for ( int i = 0; i < sz; i++ )
        os << " " << x.dx(i);
      os << " ]";
      return os;
    }

  } // namespace Fad

} // namespace AutoDiff

#endif
