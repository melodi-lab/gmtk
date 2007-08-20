//
// A very simple array class.
// This class is not meant for protection, it is only
// for convenience.
//
// $Header$
//
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <assert.h>

#include "general.h"
VCID("$Header$")
#include "error.h"

#include "cArray.h"


#ifdef MAIN

volatile int i;
void func11() { i = 0x343434;} 

class foo {



};


int
main()
{
  cArray<int> iar;
  cArray<float> far;

  iar.resize(10);
  far.resize(10);
  // for (int i=0;i<10;i++) {
  int i=0;
  do {
    func11();
    iar[i] = i;
    far[i] = i;
  } while (++i < 10);
  for (int i=0;i<10;i++) {
    printf("iar[%d]=%d,far[%d]=%f\n",i,iar[i],i,far[i]);
  }
  printf("--\n");
  iar.resizeAndCopy(5);
  far.resizeAndCopy(5);
  for (int i=0;i<5;i++) {
    printf("iar[%d]=%d,far[%d]=%f\n",i,iar[i],i,far[i]);
  }

  printf("--\n");
  iar.resizeAndCopy(15);
  far.resizeAndCopy(15);
  for (int i=0;i<15;i++) {
    printf("iar[%d]=%d,far[%d]=%f\n",i,iar[i],i,far[i]);
  }


  cArray<cArray<float> > ffar;
  ffar.resize(10);
  for (int i=0;i<10;i++) {
    ffar[i].resize(10);
  }
  for (int i=0;i<10;i++) {
    for (int j=0;j<10;j++) {
      ffar[i][j] = i*10+j;
    }
  }
  for (int i=0;i<10;i++) {
    for (int j=0;j<10;j++) {
      printf("ffar[%d][%d] = %f\n",i,j,ffar[i][j]);
    }
  }

}

#endif
