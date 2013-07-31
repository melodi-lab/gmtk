#pragma once

#if defined(HAVE_CONFIG_H)
#  include <config.h>
#endif

#if !defined(HAVE_BLAS)
#  include <math.h>

void
cblas_dcopy(int n, double const* x, int incx, double *y, int incy) {
  double const *end = x + n * incx;
  do {
    *y = *x;
    x+=incx; y+=incy;
  } while (x != end);
}

void
cblas_dscal(int n double alpha, double *x, int incx) {
  double *end = x + n * incx;
  do {
    *x = alpha * *x;
    x+=incx;
  } while ( x != end);
}

void
cblas_daxpy(int n, double alpha, double const *x, int incx, double *y, int incy) {
  double const *end = x + n * incx;
  do {  
    *y = alpha * *x + *y;
    x+=incx; y+=incy;
  } while (x != end);
}

double
cblas_ddot(int n, double const *x, int incx, double const *y, int incy) {
  double ddot = 0.0;
  double const *end = x + n * incx; 
  do {
    ddot += *x * *y;
    x+=incx; y+=incy;
  } while (x != end);
  return ddot;
}

double
cblas_dasum(int n, double const *x, int incx) {
  double const *end = x + n * incx;
  double sum = 0.0;
  do {
    sum += fabs(*x);
    x += incx;
  } while (x != end);
  return sum;
}

double
cblas_dnrm2(int n, double const *x, int incx) {
  double const *end = x + n * incx;
  double nrm = 0.0;
  do {
    nrm += *x * *x;
    x += incx;
  } while (x != end);
  return nrm;
}
#endif
