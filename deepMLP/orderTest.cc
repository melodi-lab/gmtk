
#include "Matrix.h"
#include "stdio.h"

int
main(int argc, char *argv[]) {
  double data10[10] = {0,1,2,3,4,5,6,7,8,9};
  Matrix m10(data10, 5, 2, 5, false);
  m10.Print(1);
  printf("\n\n");

  double mx_data[10] = {4,5,3,6,8,2,1,2,2,2};
  MutableMatrix mx(mx_data, 2, 5, 2, false);
  mx.Print(1);
  printf("\n\n");

  double mm_d[25];
  MutableMatrix mm(mm_d, 5, 5, 5, false);
  mm = m10 * mx;
  printf("Dgemm: %d x %d:\n", mm.NumR(), mm.NumC());
  mm.Print(1);
  printf("\n\n");

  double mmT_d[4] = {0, 0, 0, 0};
  MutableMatrix mmT(mmT_d, 2, 2, 2, false);
  mmT = m10.Trans() * mx.Trans();
  printf("Dgemm T %d x %d:\n", mmT.NumR(), mmT.NumC());
  mmT.Print(1);
  printf("\n\n");

  return 0;
}
