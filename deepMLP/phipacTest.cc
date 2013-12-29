
/*
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#if defined(HAVE_CONFIG_H)
#  include <config.h>
#endif

#include <stdio.h>

#if defined(USE_PHIPAC)
extern "C" {
void
phipac_dgemm(char* transA, char* transB,
      int* M, int* N, int* K,
      double* alpha,
      double* A, int* Astride,
      double* B, int* Bstride,
      double* beta,
      double* C, int* Cstride);
}
#endif

void
rowMajPrint(double *A, int M, int N, int stride) {
  for (unsigned i=0; i < M; i+=1) {
    for (unsigned j=0; j < N; j+=1) {
      printf("%f ", A[ i * stride + j ]);
    }
    printf("\n");
  }
}


void
colMajPrint(double *A, int M, int N, int stride) {
  for (unsigned i=0; i < M; i+=1) {
    for (unsigned j=0; j < N; j+=1) {
      printf("%f ", A[ j * stride + i ]);
    }
    printf("\n");
  }
}


void
linPrint(double *A, int len) {
  for (unsigned i=0; i < len; i+=1) {
    printf("%f ", A[i]);
  }
  printf("\n");
}


int
main(int argc, char *argv[]) {

  char op = 'n';
  int M=3, K=3, N=3;
  double alpha = 1.0, beta = 1.0;

  double X[9] = {1,1,1,0,0,0,0,0,0};
  double Y[9] = {1,2,3,4,5,6,7,8,9};
  double Z[9] = {0,0,0,0,0,0,0,0,0};

  int Xstride = 3, Ystride = 3, Zstride = 3;

#if defined(USE_PHIPAC)
  phipac_dgemm(&op, &op, 
	       &M, &N, &K, 
	       &alpha, 
	       X, &Xstride, 
	       Y, &Ystride, 
	       &beta, 
	       Z, &Zstride);
#endif

  printf("Row Major Order:\n\n");
  rowMajPrint(X, 3,3,3);
  printf("\n    *\n\n");
  rowMajPrint(Y, 3,3,3);
  printf("\n    =\n\n");
  rowMajPrint(Z, 3,3,3);
  printf("\n\n");
  
  printf("==============================\n\n");
  printf("Column Major Order:\n\n");
  colMajPrint(X, 3,3,3);
  printf("\n    *\n\n");
  colMajPrint(Y, 3,3,3);
  printf("\n    =\n\n");
  colMajPrint(Z, 3,3,3);
  printf("\n\n");
  return 0;
}
