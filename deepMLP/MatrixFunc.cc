
#if defined(HAVE_CONFIG_H)
#  include <config.h>
#endif

#if defined(HAVE_MKL)
#  include "mkl.h"
#  include "mkl_lapacke.h"
#  include "mkl_spblas.h"
#  include "mkl_trans.h"
#elif defined(HAVE_BLAS)
extern "C" {            /* Assume C declarations for C++ */
#  include <cblas.h>
}
#endif

#include "MatrixFunc.h"

double Max(const Vector & vec) {
  double max = -INFTY;
  vec.Visit([&](double x) { if (x > max) max = x; });
  return max;
}

double Sum(const Vector & vec) {
  double sum = 0;
  vec.Visit([&](double x) { sum += x; });
  return sum;
}

double NormL1(const Vector & vec) {
  return cblas_dasum(vec.Len(), vec.Start(), vec.Inc());
}

double LogSumFastDestroy(MutableVector & vec) {
  assert (vec.Inc() == 1);

  double max = Max(vec);

  double *s = vec.Start();
  double *tP = vec.Start();
  vec.Visit([&](double x) {
    double val = x - max;
    if (val > -30) *tP++ = val;
  });

  int nT = (int)(tP - vec.Start());
#if HAVE_MKL
  vdExp(nT, s, s);
#else
  MutableVector temp(s, nT);
  temp.Apply([](double x)->double {return exp(x);});
#endif
  double sumExp = cblas_dasum(nT, s, 1);

  return max + log(sumExp);
}

void LogSumWithNeg(MutableVector & vec, MutableVector & temp) {
  int len = vec.Len();
  assert (temp.Len() == len);
  assert (vec.Inc() == 1 && temp.Inc() == 1);
#if HAVE_MKL	
  vdAbs(len, vec.Start(), vec.Start());
#else
  vec.Apply([](double x)->double {return abs(x);});
#endif
#if HAVE_CBLAS_DAXPBY
  cblas_daxpby(len, -2.0, vec.Start(), 1, 0, temp.Start(), 1);
#else
  cblas_dscal(len, 0, temp.Start(), 1);
  cblas_daxpy(len, -2.0, vec.Start(), 1, temp.Start(), 1);
#endif
#if HAVE_MKL
  temp.ApplyVML(vdExp, temp);
  temp.ApplyVML(vdLog1p, temp); // note: might be faster to increment and log!
  vdAdd(len, temp.Start(), vec.Start(), vec.Start());
#else
  temp.Apply([](double x)->double {return exp(x);});
  temp.Apply([](double x)->double {return log1p(x);});
  cblas_daxpy(len, 1.0, temp.Start(), 1, vec.Start(), 1);
#endif
}

void Trunc(MutableVector & vec, double maxVal) {
  vec.Apply([maxVal] (double x) { return (x < -maxVal) ? -maxVal : (x > maxVal) ? maxVal : x; });
}

void Shrink(MutableVector & vec, double d) {
  vec.Apply([d](double x) { return (x>d) ? x-d : (x<-d) ? x+d : 0; });
}

double Norm(const Vector & vec) {
  return cblas_dnrm2(vec.Len(), vec.Start(), vec.Inc());
}

double Decay(MutableVector vec, double stepSize, double decayRate, int decayType) {
  if (decayType == 2) {
    double val = 0.5 * decayRate * vec * vec;
    vec *= (1.0 - stepSize * decayRate);
    return val;
  } else {
    double val = decayRate * NormL1(vec);
    Shrink(vec, stepSize * decayRate);
    return val;
  }
}

void Print(const Vector & vec) {
  vec.Visit([] (double x) { cout << x << endl; });
}

int ArgMax(const Vector & vec) {
  double max = -INFTY;
  int argMax = -1;
  int i = 0;
  vec.Visit([&] (double x) { if (x > max) { max = x; argMax = i; } ++i; });
  return argMax;
}


