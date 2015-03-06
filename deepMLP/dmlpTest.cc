
/*
 * Written by Richard Rogers rprogers@uw.edu
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "Matrix.h"
#include "stdio.h"

int
main(int argc, char *argv[]) {
  double data10[14] = {0,1,2,3,4, 0, 0, 5,6,7,8,9, 0, 0};
  Matrix m10(data10, 5, 2, 7, false);

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
  printf("my_Domatadd %d x %d:\n", m10s.NumR(), m10s.NumC());
  m10s.Print();
  printf("\n\n");



  double m10T_add_data[20]={0,10,0,0, 2,12,0,0, 4,14,0,0, 6,16,0,0, 8,18,0,0};
  MutableMatrix m10T_add(m10T_add_data, 2,5,4, false);
  //  m10T_add = m10c.Trans();
  printf("m10T_add %d x %d:\n", m10T_add.NumR(), m10T_add.NumC());
  m10T_add.Print();

  // double add_data[20] = {0,5,0,0, 1,6,0,0, 2,7,0,0, 3,8,0,0, 4,9,0,0};
  //MutableMatrix add(add_data, 2,5,4, false);
  double add_data[14] = {0,1,2,3,4, 20,30, 5,6,7,8,9, 40,50};
  MutableMatrix add(add_data, 5,2,7, false);
  printf("add %d x %d:\n", add.NumR(), add.NumC());
  add.Print();

  //  add = add.Trans();
  printf("add T %d x %d:\n", add.NumC(), add.NumR());
  add.Trans().Print();

  double sum_data[20];
  MutableMatrix sum(sum_data, 2,5,2, false);
  sum = /*1.0 * */ m10T_add + 2.0 * add.Trans();
  printf("m10T_add + 2 * add:\n");
  sum.Print();
  printf("\n\n\n");



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

  double mx_data[20] = {4,5, 1,1, 3,6, 2,2, 8,2, 3,3, 1,2, 4,4, 2,2, 5,5};
  MutableMatrix mx(mx_data, 2, 5, 4, false);
  printf("mx 2.4 x 5:\n");
  mx.Print();
  printf("\n\n");

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
