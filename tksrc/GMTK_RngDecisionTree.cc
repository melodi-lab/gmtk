/*-
 * GMTK_RngDecisionTree.cc
 *     General Class to map from vectors of integers to some basic type.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"

#include "GMTK_RngDecisionTree.h"


VCID("$Header$");

////////////////////////////////////////////////////////////////////
//        Static Data
////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////







#ifdef MAIN

////////////////////////////////////////////////////////////////////
//        test code Support
////////////////////////////////////////////////////////////////////

int
main()
{
  // read in a file
  iDataStreamFile is ("foo.dt",false);
  RngDecisionTree<int> dt;
  dt.read(is);

  printf("Found decision tree\n");
  oDataStreamFile os("-",false);
  dt.write(os);

  sArray<int> vec;
  vec.resize(dt.numFeatures());
  iDataStreamFile stin ("-",false);

  while (1) {
    printf("Enter a length %d intvec:",dt.numFeatures());
    stin.readArray(vec.ptr,dt.numFeatures());
    printf("Result of querying with vector: ");
    for (int i=0;i<dt.numFeatures();i++) {
      printf(" %d",vec[i]);
    }
    printf("\n is %d\n",dt.query(vec));
  }
}



#endif
