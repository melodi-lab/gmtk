#pragma once

#if defined(HAVE_CONFIG_H)
#  include <config.h>
#endif

#if !defined(HAVE_BLAS) && !defined(HAVE_MKL)
#  include <math.h>

void cblas_dcopy(int n, double const* x, int incx, double *y, int incy);

void cblas_dscal(int n, double alpha, double *x, int incx);

void cblas_daxpy(int n, double alpha, double const *x, int incx, double *y, int incy);

double cblas_ddot(int n, double const *x, int incx, double const *y, int incy);

double cblas_dasum(int n, double const *x, int incx);

double cblas_dnrm2(int n, double const *x, int incx);

#endif
