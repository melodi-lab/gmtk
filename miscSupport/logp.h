//
//
// Log probability class for single or double precision numbers.
// Some inspiration for this code was taken from HTK's HMath.c
//
// $Header$
//
// Written by: Jeff A. Bilmes
//             bilmes@ee.washington.edu


#ifndef LOGP_H
#define LOGP_H

#include <math.h>
#include <assert.h>

#include "error.h"

// use a lookup table? If so, define the following,
// or otherwise, undefine to use log(1+exp(x)) directly.
// #define _TABLE_

// What we call ~log(0). 
#define LZERO  (-1.0E10)

// Things that result in log values < LSMALL are set to LZERO.
// I.e., this is the smallest log value that we represent
// after the result of an addition, etc. 
#define LSMALL (-0.5E10)

// Lowest exp() arg  = log(MINLARG)
// I.e., exp(z) where z < MINEARG results in 0
#define MINEARG (-708.3)

// Largest exp() arg, i.e., if v is much
// greater than this, then we'll get an infinity.
#define MAXEARG (700)

// lowest log() arg  = exp(MINEARG)
// I..e, log(z) where z < MINLARG results in -Infinity
#define MINLARG 2.45E-308  

// NOTE: On some machines it may be necessary to reduce the
//       values of MINEARG and MINLARG depending on the
//       implementation of log() and exp()


// The types logp<float> and logp<double> are used for representing
// real numbers on a log scale.  LZERO is used for log(0) 
// in log arithmetic, any log real value <= LSMALL is 
// considered to be zero.

extern double logp_minLogExp;
extern double log_FLT_MIN;
extern double log_DBL_MIN;

#ifdef _TABLE_
  const int table_size = 50000;
#endif

//
// FT should either be float or double.
template <
  class FT,          // the storage type
  class iFT = double // the internal temporary type
>
class logp {

private:
  FT v;
#ifdef _TABLE_
  static FT table[table_size];
  static double inc;
#endif

public:

#ifdef _TABLE_
  static bool table_init();
#endif

  logp(const float p) { 
    assert ( p >= 0.0 );

    if (p<0) // check regardless
      coredump("logp: negagive constructor argument.");
    v = (p<MINLARG?LZERO:log(p));
    // Note: LZERO is the smallest log that we represent
    // and LZERO < log(MINLARG), but we can't take
    // the log of anything smaller than MINLARG. Therefore,
    // for conversion purposes, we check p against
    // MINLARG rather than LZERO or LSMALL
  }
  // Same as above, but construct with a double, so we do the log() 
  // before converting to a float rather than converting to a float 
  // and then doing  the log().
  logp(const double p) { 
    assert ( p >= 0.0 );
    if (p<0)
      coredump("logp: negagive constructor argument.");
    v = (p<MINLARG?LZERO:log(p));
  }

  // No arguments causes zero value.
  inline logp() : v(LZERO) {} 

  // A constructor that does no initialization.
  inline logp(void *dummy_not_used) {}

  // Construct from another logp object.
  inline logp(const logp<FT,iFT>& lp) { v = lp.v; }

  // A constructor to initialize from a log probability value
  // stored in a floating point value.
  inline logp(void *dummy_not_used,const float logpval) {
    if (logpval < LSMALL)
      v = LZERO;
    else 
      v = logpval;
  }
  // Have an extra double version since we can then
  // check if (logpval < LSMALL) before conv. to float.
  inline logp(void *dummy_not_used,const double logpval) {
    if (logpval < LSMALL)
      v = LZERO;
    else 
      v = logpval;
  }

  // these are here for quick use so as to avoid compiler temporaries
  // and/or conversion with comparing with float values
  inline bool zero() const { return (v == LZERO); }
  inline bool essentially_zero() const { return (v < LSMALL); }
  inline bool not_essentially_zero() const { return (v >= LSMALL); }
  inline bool almost_zero() const { return (v == LSMALL); }
  inline void set_to_zero() { v = LZERO; }
  inline void set_to_almost_zero() { v = LSMALL; }
  inline void set_to_one() { v = 0; }

  // log addition: Return sum x + y on log scale, 
  // sum < LSMALL is floored to LZERO
  friend inline 
    logp<FT,iFT> operator+(logp<FT,iFT> x, logp<FT,iFT> y) {
    if (x.v<y.v) {
      FT temp = x.v; x.v = y.v; y.v = temp;
    }
    // now y<=x
    iFT diff = y.v-x.v;
    if (diff<logp_minLogExp) {
      // logp_minLogExp == -log(-LZERO)
      // So, if y-x = log(Py/Px) < -log(-LZERO)
      //              log(Px/Py) > log(-LZERO)
      //                  Px/Py > -LZERO
      // LZERO is kind of like log(0) = -inf, so this is saying;
      //                  Px > Inf*Py
      // If this condition is true, than x is *much* larger
      // than y, so we don't even bother adding y.
      if (x.v<LSMALL) {
	// Since x is so small, we return x's value but
	// changed to what we call zero.
	x.v = LZERO;
      }
      // Just return x, disregarding y.
      return x;
    } else {
      // note:
      // x = log(px), y = log(py) where px >= py
      // We want z = log(px+py) = 
      //             log(exp(log(px))+exp(log(py))) =
      //    (1)      log(exp(x)+exp(y)) =
      //             log[exp{x}+exp{x+y-x}] =
      //             log[exp{x} x (1 + exp(y-x))] =
      //             x + log[1 + exp(y-x)] 
      // Why do we do this rather than equation (1)?
      //  a) It is two rather than 3 calls to transcendentals.
      //  b) Is ammenable to integer representation of log(probs)
      //  c) log(1+e^n) can be pre-computed in a table, so
      //     to compute it, we just do a table lookup with y-x.
      //  d) It is numerically better behaved than (1) above, e.g.:
      //      I: if x ~== y and are very small, (i.e., less
      //         than MINEARG), result is still good (unlike in (1))
      //      II: If py happens to be *much* smaller than px, then
      //          we don't loose precision by needlessly evaluating
      //          x = log(exp(x) + ~0)
      //          Even more importantly, if x or y < MINEARG, (which
      //          could easily happen), then z = log(exp(x) + exp(y))
      //          would at best produce garbage.
      logp<FT,iFT> z((void*)0);
      // could use table lookup here, or if we could
      // write a function to do log(1+exp(x)) directly w/o two calls.
#ifdef _TABLE_
      z.v = x.v + table[int(diff*inc)]; 
#else
      z.v = x.v+log(1.0+exp(diff)); 
#endif

      return z;
    }
  }

  inline logp<FT,iFT>& operator +=(logp<FT,iFT> z)
  { *this = (*this + z); return *this; }

  // log subtraction: Return diff x - y on log scale, 
  // diff < LSMALL is floored to LZERO
  friend inline 
    logp<FT,iFT> operator-(logp<FT,iFT> x, logp<FT,iFT> y) {
    if (x.v<y.v)    
      coredump("logp::operator-: result would be negative probability.");
    // now y<=x
    iFT diff = y.v-x.v;
    if (diff<logp_minLogExp) {
      // then x is *much* greater than y, so we don't
      // even bother subtracting off y.
      if (x.v < LSMALL)
	x.v = LZERO;
      return x;
    } else {
      // We want z = log(px-py) = 
      //             x + log[1 - exp(y-x)] 
      iFT tmp = 1.0 - exp(diff);
      logp<FT,iFT> z((void*)0);
      if (tmp<MINLARG) 
	// then y ~== x, return zero.
	z.v = LZERO;
      else 
	z.v = x.v+log(tmp);  // again, could use table lookup here
      return z;
    }
  }

  inline logp<FT,iFT>& operator -=(logp<FT,iFT> z)
  { *this = (*this - z); return *this; }
  
  friend inline 
    logp<FT,iFT> operator*(logp<FT,iFT> x, logp<FT,iFT> y)  
    { // if (x.zero()) return x; 
      // if (y.zero()) return y; 
      // The above probably isn't needed since it will
      // eventually be floored to LZERO on a subsequent add.
      x.v += y.v; return x; }
  inline logp<FT,iFT>& operator *=(logp<FT,iFT> z)
  { *this = (*this * z); return *this; }

  friend inline 
    logp<FT,iFT> max(logp<FT,iFT> x, logp<FT,iFT> y)  
  { 
    if (x.v > y.v)
      return x;
    else 
      return y;
  }
  inline void assign_if_greater(logp<FT,iFT> z)
  { 
    if (z.v > v)
      v = z.v;
  }

  inline logp<FT,iFT> inverse() const {
    void *dummy=NULL;
    return logp<FT,iFT>(dummy,-v);
  }

  inline logp<FT,iFT> pow(double pwr) const {
    void *dummy=NULL;
    return logp<FT,iFT>(dummy,pwr*v);
  }




  friend inline
    logp<FT,iFT> operator/(logp<FT,iFT> x, logp<FT,iFT> y) 
    { 
      if (y.zero()) coredump("logp: divide by zero.");
      if (x.zero()) return x;
      x.v -= y.v; // might need to worry about underflow here.
      return x; }
  inline logp<FT,iFT>& operator /=(logp<FT,iFT> z)
    { *this = (*this / z); return *this; }
  

  // unlog support: Convert log(x) to double, result is
  // floored to 0.0 if x < MINEARG
  inline iFT unlog() const { return (v<MINEARG) ? 0.0 : exp(v); }
  // return true if we can safely call unlog() on this object.
  bool unlogable() { return (v < MAXEARG); }
  // return true if we can safely call unlog() on the inverse of
  // this object. This can also be used as a weaker check for zero
  // than the zero() function (i.e., if zero(), then 
  // inverseUnlogable() will be true)
  bool inverseUnlogable() { return (-v < MAXEARG); }

  inline FT val() const { return v; }

  // this routine breaks all protection of 'v'
  inline FT& valref() { return v; }

  // Set our internal value from normal probability p
  logp<FT,iFT>& setFromP (FT p) {   
    if (p<0)
      coredump("logp: negagive constructor argument.");
    v = (p<MINLARG?LZERO:log(p));
    return *this;
  }
  // Set our internal value from a log'ed probability lp
  logp<FT,iFT>& setFromLogP (FT lp) {   
    if (lp < LSMALL)
      v = LZERO;
    else 
      v = lp;
    return *this;
  }

  // floor self to zero.
  void floor() {
    if (v < LSMALL)
      v = LZERO;
  }

  logp<FT,iFT>& operator = (float p) { 
    if (p<0)
      coredump("logp: negagive constructor argument.");
    v = (p<MINLARG?LZERO:log(p));
    return *this;
  }
  logp<FT,iFT>& operator = (double p) { 
    if (p<0)
      coredump("logp: negagive constructor argument.");
    v = (p<MINLARG?LZERO:log(p));
    return *this;
  }

  inline bool operator== (logp<FT,iFT>y) const { return v == y.v; }
  inline bool operator!= (logp<FT,iFT>y) const { return v != y.v; }
  inline bool operator<  (logp<FT,iFT>y) const { return v < y.v; }
  inline bool operator<= (logp<FT,iFT>y) const { return v <= y.v; }
  inline bool operator>  (logp<FT,iFT>y) const { return v > y.v; }
  inline bool operator>= (logp<FT,iFT>y) const { return v >= y.v; }
  

};

/////////////////////////////////////////////////////
// The default log probability type used by users 
// of this class.
typedef logp<double,double> logpr;
/////////////////////////////////////////////////////

#endif
