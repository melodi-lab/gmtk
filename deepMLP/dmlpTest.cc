
#include "Matrix.h"
#include "stdio.h"

int
main(int argc, char *argv[]) {
  double data10[10] = {0,1,2,3,4,5,6,7,8,9};
  Matrix m10(data10, 5, 2, 5, false);

  printf("m10 %d x %d:\n", m10.NumR(), m10.NumC());
  m10.Print();
  printf("\n\n");

  Matrix m10T = m10.Trans();
  printf("m10T %d x %d:\n", m10T.NumR(), m10T.NumC());
  m10T.Print();
  printf("\n\n");

  AllocatingMatrix m10c;
  m10c.Resize(5,2);
  m10c.CopyFrom(m10, 2.0);
  printf("my_Domatcopy %d x %d:\n", m10c.NumR(), m10c.NumC());
  m10c.Print();
  printf("\n\n");

  AllocatingMatrix m10TT;
  m10TT.Resize(2,5);
  m10TT.CopyFrom(m10T);
  printf("my_Domatcopy T %d x %d:\n", m10TT.NumR(), m10TT.NumC());
  m10TT.Print();
  printf("\n\n");

  double m10s_d[10];
  MutableMatrix m10s(m10s_d, 5, 2, 5, false);
  m10s = 1.0 * m10 + 2.0 * m10c;
  printf("my_Domatadd T %d x %d:\n", m10s.NumR(), m10s.NumC());
  m10s.Print();
  printf("\n\n");

  double mm_d[25];
  MutableMatrix mm(mm_d, 5, 5, 5, false);
  mm = m10 * m10T;
  printf("Dgemm: %d x %d:\n", mm.NumR(), mm.NumC());
  mm.Print();
  printf("\n\n");

  mm = m10c * m10T;
  printf("Dgemm: %d x %d:\n", mm.NumR(), mm.NumC());
  mm.Print();
  printf("\n\n");

  double mx_data[10] = {4,5,3,6,8,2,1,2,2,2};
  MutableMatrix mx(mx_data, 2, 5, 2, false);
  mm = m10 * mx;
  printf("Dgemm: %d x %d:\n", mm.NumR(), mm.NumC());
  mm.Print();
  printf("\n\n");

  double mmT_d[4] = {0, 0, 0, 0};
  MutableMatrix mmT(mmT_d, 2, 2, 2, false);
  mmT = m10.Trans() * mx.Trans();
  printf("Dgemm T %d x %d:\n", mmT.NumR(), mmT.NumC());
  mmT.Print();
  printf("\n\n");

  return 0;
}
