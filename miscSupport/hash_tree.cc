/*-
 * GMTK_HashTree.cc
 *     Hash Tree
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
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
#include <ctype.h>

#include <iterator>
#include <map>
#include <set>
#include <algorithm>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"

#include "hash_tree.h"

VCID("$Header$");

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Member functions
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


#ifdef MAIN

///////////////////////////////////////////
// main driver debugger for hash table.

#include <string>

int
main()

{
  const unsigned depth = 4;

  sArray<unsigned> start_sizes(depth);
  for (unsigned i=0;i<depth;i++) {
    start_sizes[i] = 20;
  }

  hash_tree<unsigned,double> t(depth,start_sizes.ptr);

  sArray<unsigned> v(depth);
  v[0] = 0;
  v[1] = 1;
  v[2] = 2;
  v[3] = 3;

  printf("before any insert, find = %d\n",(bool)t.find(v));

  // insert a bunch of entries
  unsigned i0,i1,i2,i3;
  for (i0=0;i0<8;i0++) {
    v[0] = i0;
    double val0 = i0;
    for (i1=0;i1<8;i1++) {
      double val1 = 10*val0 + i1;
      v[1] = i1;
      for (i2=0;i2<8;i2++) { 
	double val2 = 10*val1 + i2;
	v[2] = i2;
	for (i3=0;i3<8;i3++) {
	  double val3 = 10*val2 + i3;
	  v[3] = i3;
	  t.insert(v,val3);
	}
      }
    }
  } 

  printf("after insert of %d entries, find = %d\n",
	 t.numEntries(),(bool)t.find(v));
  
  printf("finding a bunch of random elements\n");
  for (unsigned count=0;count<100;count++) {
    for (unsigned vi=0;vi<depth;vi++) {
      v[vi] = rand() % 10;
    }
    printf("find of: ");
    for (unsigned vi=0;vi<depth;vi++)
      printf("%d ",v[vi]);
    printf(": = %d\n",(bool)t.find(v));
  }

  printf("iterating over all elements\n");
  fflush(stdout);
  hash_tree<unsigned,double>::iterator_vector it;
  for (t.begin(it,v.ptr);it!=t.end();it++) {
    for (unsigned i=0;i<depth;i++)
      printf("%d ",v[i]);
    printf(":value = %f\n",(*it));
  }

  printf("iterating over subset of elements elements\n");
  fflush(stdout);
  unsigned start[2] = { 4, 4 };
  for (t.begin(start,2,it,v.ptr);it!=t.end();it++) {
    for (unsigned i=0;i<depth;i++)
      printf("%d ",v[i]);
    printf(":value = %f\n",(*it));
  }


  printf("iterating over all elements using vecp\n");
  fflush(stdout);
  hash_tree<unsigned,double>::iterator_vectorp itp;
  sArray<unsigned*> vp(depth);
  for (unsigned i=0;i<depth;i++) {
    vp[i] = &v[i];
  }
  for (t.begin(itp,vp.ptr);itp!=t.end();itp++) {
    for (unsigned i=0;i<depth;i++)
      printf("%d ",v[i]);
    printf(":value = %f\n",(*itp));
  }

  printf("iterating over subset of elements elements vecp\n");
  fflush(stdout);
  for (t.begin(start,2,itp,vp.ptr);itp!=t.end();itp++) {
    for (unsigned i=0;i<depth;i++)
      printf("%d ",v[i]);
    printf(":value = %f\n",(*itp));
  }


}


#endif /* defined MAIN */
