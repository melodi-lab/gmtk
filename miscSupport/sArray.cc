//
// A very simple array class.
// This class is not meant for protection, it is only
// for convenience.
//
// $Header$
//
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <values.h>
#include <math.h>
#include <ieeefp.h>
#include <float.h>
#include <assert.h>

#include "general.h"
VCID("$Header$");
#include "error.h"

#include "sArray.h"



#ifdef MAIN

volatile int i;
void func11() { i = 0x343434;} 


int
main()
{
  sArray<int> iar;
  sArray<float> far;

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
  iar.resizeAndCopy(5);
  far.resizeAndCopy(5);
  for (int i=0;i<5;i++) {
    printf("iar[%d]=%d,far[%d]=%f\n",i,iar[i],i,far[i]);
  }

  sArray<sArray<float> > ffar;
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
