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

/*
 * A decision tree looks something like:
 * <numFtrs>
 * <ftr> <num_splits> r1 r2 ... rs 
 * <ftr> <num_splits> r1 r2 ... rs # split 1
 * <ftr> <num_splits> r1 r2 ... rs # split 2
 * -1 value
 * <ftr> <num_splits> r1 r2 ... rs # split s
 */

char *dtStr =
"% this is a decision tree file\n"
"%\n"
"dt_name 3  % number of features\n"
"0 10 0:5 6:9 10 11 12 13 14 15 50: default\n"
"  1 2 0:10 default\n"
"    2 2 0:10 default\n"
"      -1 expand\n"
"      -1 (p0+1)\n"
"    2 2 0:5 default\n"
"      -1 ( c0 + 1 )\n"
"      -1 (m0 +1)\n"
"  1 2 0:10 default\n"
"    2 2 0:10 default\n"
"      -1 (p0+p1+ 5)\n"
"      -1 6\n"
"    2 2 0:5 default\n"
"      -1 7\n"
"      -1 8\n"
"  -1 10 % when feature[0] = 10, map to 10 regardless of all else\n"
"  -1 11 % when feature[0] = 11, map to 11 regardless of all else\n"
"  -1 12 % when feature[0] = 12, map to 12 regardless of all else\n"
"  -1 13 % when feature[0] = 13, map to 13 regardless of all else\n"
"  -1 14 % when feature[0] = 14, map to 14 regardless of all else\n"
"  -1 15 % when feature[0] = 15, map to 15 regardless of all else\n"
"  -1 16 % when feature[0] >= 50, map to 16 regardless of all else\n"
"  1 2 0:10 default\n"
"    2 2 0:10 default\n"
"      -1 9\n"
"      -1 10\n"
"    2 2 0:5 default\n"
"      -1 11\n"
"      -1 12\n";


int
main()
{

  // first write out the file
  {
    oDataStreamFile dtfile ("/tmp/foo.dt",false);
    dtfile.write(dtStr);
  }

  // read it in again
  iDataStreamFile is ("/tmp/foo.dt",false);
  RngDecisionTree<int> dt;
  dt.read(is);

  printf("Found decision tree\n");
  oDataStreamFile os("-",false);
  dt.write(os);

  vector<int> vec;
  vector<int> card;
  vec.resize(dt.numFeatures());
  card.resize(dt.numFeatures());
  iDataStreamFile stin ("-",false,false);

  // first test iterating through all leaf values.
  for (RngDecisionTree<int>::iterator it = dt.begin();
       it != dt.end(); it++) {
    printf("leaf value = %d\n",it.value());
  }

  while (1) {
    printf("Enter a length %d intvec:",dt.numFeatures());
    fflush(stdout);
    stin.read(vec,dt.numFeatures());
    printf("Enter a length %d set of cardinalities:",dt.numFeatures());
    fflush(stdout);
    stin.read(card,dt.numFeatures());

    printf("Result of querying with vector and cards: ");
    fflush(stdout);
    for (unsigned i=0;i<dt.numFeatures();i++) {
      printf(" %d:%d",vec[i],card[i]);
    }
    printf("\n is %d\n",dt.query(vec,card));
  }

}



#endif
