
#include <stdio.h>

#include "Matrix.h"

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

int
main(int argc, char *argv[]) {
  double A[10] = {0,1,2,3,4,5,6,7,8,9};
  double B[10] = {4,5,3,6,8,2,1,2,2,2};

  double C[25];

  char op = 'n';
  int M = 5, K = 2, N = 5;
  int Astride = 5, Bstride = 2, Cstride = 5;
  double alpha = 1.0, beta = 1.0;

  phipac_dgemm(&op, &op, 
	       &M, &N, &K, 
	       &alpha, 
	       A, &Astride, 
	       B, &Bstride, 
	       &beta, 
	       C, &Cstride);


  MutableMatrix mm(C, 5, 5, 5, false);
  printf("Dgemm %d x %d:\n", mm.NumR(), mm.NumC());
  mm.Print(1);
  printf("\n\n");

  return 0;
}
