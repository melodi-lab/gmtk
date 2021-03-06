//
//  Copyright (C) 2001 Jeff Bilmes
//  Licensed under the Open Software License version 3.0
//  See COPYING or http://opensource.org/licenses/OSL-3.0
//
//
// Linear equation solver routines.
// 
//          Jeff Bilmes
//          bilmes@ee.washington.edu


#if HAVE_CONFIG_H
#include <config.h>
#endif
#if defined(HAVE_SYSEXITS_H)
#  include <sysexits.h>
#endif
#ifndef EX_TEMPFAIL
// BSDish exit status encouraging the user to try again
#  define EX_TEMPFAIL 75
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define LINEQSOLVE_USE_MALLOC

#ifndef LINEQSOLVE_USE_MALLOC
#include <alloca.h>
#endif

#include "hgstamp.h"
#include "general.h"
VCID(HGID)
#include "error.h"
#include "lineqsolve.h"


#define NRANSI
#define SINGLE_TINY 1.0e-20;
#define DOUBLE_TINY 1.0e-150;

///////////////////////////////////////////////////////////////
// Single precision versions
///////////////////////////////////////////////////////////////

void ludcmp(float *a, // nXn matrix
	    int n, 
	    int *indx, // row permutation by partial piviting.
	    float *d)
{
  int i,imax=0,j,k;
  float big,dum,sum,temp;
  float *vv;

  if ((vv = (float*) malloc((size_t) (n*sizeof(float)))) == NULL) {
    fprintf(stderr, "ludcmp: Can't allocate memory\n");
    exit(EX_TEMPFAIL);
  }
  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i*n+j])) > big) 
	big=temp;
    if (big == 0.0) {
      error("ERROR: LU Decomposition routine given a singular matrix");
    }
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i*n+j];
      for (k=0;k<i;k++)
	sum -= a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i*n+j];
      for (k=0;k<j;k++)
	sum -= a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
	dum=a[imax*n+k];
	a[imax*n+k]=a[j*n+k];
	a[j*n+k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j*n+j] == 0.0) 
      a[j*n+j]=SINGLE_TINY;
    if (j != n) {
      dum=1.0/(a[j*n+j]);
      for (i=j+1;i<n;i++) 
	a[i*n+j] *= dum;
    }
  }
  free((void*)vv);
}


void lubksb(float *a, 
	    int n, 
	    int *indx, 
	    float *b)
{
  int i,ii=-1,ip,j;
  float sum;

  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0) {
      for (j=ii;j<=(i-1);j++) 
	sum -= a[i*n+j]*b[j];
    } else if (sum) 
      ii=i;
    b[i]=sum;
  }
  for (i=(n-1);i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) 
      sum -= a[i*n+j]*b[j];
    b[i]=sum/a[i*n+i];
  }
}


/*
 *
 * Solves AX=B
 * where A is nXn
 *       X is nXr
 *       B is nXr
 * 
 * The resulting X is placed in B.
 */
void
lineqsolve(const int n, const int nrhs,
	   float *a,  /* the nXn A matrix */
	   float *b)  /* the rXn B^T matrix  and
			  the rXN resulting X^T. */
{

  float d;
  int i;
  float *bp;

#ifdef LINEQSOLVE_USE_MALLOC
  static int *indx = NULL;
  static int indx_len = 0;
  if (indx_len < n) {
    if (indx_len > 0) {
      free((void*)indx);
    }
    indx = (int*)malloc((size_t)(n*sizeof(int)));
    indx_len = n;
  }
#else
  int *indx = (int*)alloca((size_t)(n*sizeof(int)));
#endif
  if (indx == NULL) {
    fprintf(stderr, "lineqsolve: Can't allocate memory\n");
    exit(EX_TEMPFAIL);
  }
  ludcmp(a,n,indx,&d);

  bp = b;
  for (i=0;i<nrhs;i++) {
    lubksb(a,n,indx,bp);
    bp += n;
  }
}



///////////////////////////////////////////////////////////////
// Double precision versions
///////////////////////////////////////////////////////////////

void ludcmp(double *a, // nXn matrix
	    int n, 
	    int *indx, // row permutation by partial piviting.
	    double *d)
{
  unsigned int i,imax=0,j,k, nn = (unsigned) n;
  double big,dum,sum,temp;
  double *vv;

  if ((vv= (double*) malloc((size_t) (n*sizeof(double)))) == NULL) {
    fprintf(stderr, "ludcmp: Can't allocate memory\n");
    exit(EX_TEMPFAIL);
  }
  *d=1.0;
  for (i=0;i<nn;i++) {
    big=0.0;
    for (j=0;j<nn;j++)
      if ((temp=fabs(a[i*nn+j])) > big) 
	big=temp;
    if (big == 0.0) {
      error("ERROR: LU Decomposition routine given a singular matrix");
    }
    vv[i]=1.0/big;
  }
  for (j=0;j<nn;j++) {
    for (i=0;i<j;i++) {
      sum=a[i*nn+j];
      for (k=0;k<i;k++)
	sum -= a[i*nn+k]*a[k*nn+j];
      a[i*nn+j]=sum;
    }
    big=0.0;
    for (i=j;i<nn;i++) {
      sum=a[i*nn+j];
      for (k=0;k<j;k++)
	sum -= a[i*nn+k]*a[k*nn+j];
      a[i*nn+j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<nn;k++) {
	dum=a[imax*nn+k];
	a[imax*nn+k]=a[j*nn+k];
	a[j*nn+k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j*nn+j] == 0.0) 
      a[j*nn+j]=DOUBLE_TINY;
    if (j != nn) {
      dum=1.0/(a[j*nn+j]);
      for (i=j+1;i<nn;i++) 
	a[i*nn+j] *= dum;
    }
  }
  free((void*)vv);
}


void lubksb(double *a, 
	    int n, 
	    int *indx, 
	    double *b)
{
  int i,ii=-1,ip,j;
  double sum;

  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0) {
      for (j=ii;j<=(i-1);j++) 
	sum -= a[i*n+j]*b[j];
    } else if (sum) 
      ii=i;
    b[i]=sum;
  }
  for (i=(n-1);i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) 
      sum -= a[i*n+j]*b[j];
    b[i]=sum/a[i*n+i];
  }
}


/*
 *
 * Solves AX=B
 * where A is nXn
 *       X is nXr
 *       B is nXr
 * 
 * The resulting X is placed in B.
 */
void
lineqsolve(const int n, const int nrhs,
	   double *a,  /* the nXn A matrix */
	   double *b)  /* the rXn B^T matrix  and
			  the rXN resulting X^T. */
{

  double d;
  int i;
  double *bp;

#ifdef LINEQSOLVE_USE_MALLOC
  static int *indx = NULL;
  static int indx_len = 0;
  if (indx_len < n) {
    if (indx_len > 0) {
      free((void*)indx);
    }
    indx = (int*)malloc((size_t)(n*sizeof(int)));
    indx_len = n;
  }
#else
  int *indx = (int*)alloca((size_t)(n*sizeof(int)));
#endif
  if (indx == NULL) {
    fprintf(stderr, "lineqsolve: Can't allocate memory\n");
    exit(EX_TEMPFAIL);
  }
  ludcmp(a,n,indx,&d);

  bp = b;
  for (i=0;i<nrhs;i++) {
    lubksb(a,n,indx,bp);
    bp += n;
  }
}

//////////////////////////////////////////////
// Driver program
//////////////////////////////////////////////



#ifdef MAIN

int 
main(int argc, char *argv[])
{
  int n = 5;
  int r = 5;
  float *a;
  float *b;
  int i,j;

  if (argc > 1) {
    n = atoi(argv[1]);
  }
  if (argc > 2) {
    r = atoi(argv[2]);
  }

  a = (float*)malloc(n*n*sizeof(float));
  b = (float*)malloc(r*n*sizeof(float));
  for (i=0;i<n*n;i++) {
    a[i] = drand48();
  }
  for (i=0;i<n*r;i++) {
    b[i] = drand48();
  }

  printf("A = [\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      printf("%.10f ",a[i*n+j]);
    }
    printf("\n");
  }
  printf("];\n");

  printf("B = [\n");
  for (i=0;i<n;i++) {
    for (j=0;j<r;j++) {
      printf("%.10f ",b[j*n+i]);
    }
    printf("\n");
  }
  printf("];\n");

  lineqsolve(n,r,a,b);

  printf("X = [\n");
  for (i=0;i<n;i++) {
    for (j=0;j<r;j++) {
      printf("%.10f ",b[j*n+i]);
    }
    printf("\n");
  }
  printf("];\n");

  printf("A*X-B\n");

  free((void*)a);
  free((void*)b);

  return 0;
}


#endif
