/*-
 * vhash_set.cc
 *     vhash_set driver class class
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

#include "config.h"
#include "general.h"
VCID("$Header$")

#ifdef MAIN
#define COLLECT_COLLISION_STATISTICS
#include <string>

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#if HAVE_LIMITS_H
#  include <limits.h>
#endif

#include "vhash_set.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Main variables
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

#ifdef MAIN

void
myprintv(int* v,int len)
{
  for (int i=0;i<len;i++) {
    printf("%d ",v[i]);
  }
}


int
main(int argc,char*argv[])
{

  // first insert some numbers.
  int count= 100000;
  int vsize = 3;
  int maxCard = 100000;

  srand(0);
  
  if (argc > 1)
    count = atoi(argv[1]);
  if (argc > 2)
    vsize = atoi(argv[2]);
  if (argc > 3)
    maxCard = atoi(argv[3]);
  if (argc > 4)
    hash_abstract::loadFactor = atof(argv[4]);

  assert ( hash_abstract::loadFactor > 0 && hash_abstract::loadFactor < 1.0 );


  if (maxCard == 0)
    maxCard = INT_MAX;

  printf("Using %d entries, each of length %d, maxCard = %d\n",count,vsize,maxCard);

  vhash_set< UInt32 > ht(vsize,1000);

  // printf("creating %d random entries of size %d\n",count,vsize);
  UInt32* data = new UInt32[vsize*count];
  UInt32* vp = data;
  for (int i=0;i<vsize*count;i++) {
    *vp++ = rand() % maxCard;
  }


  printf("inserting %d random entries of size %d\n",count,vsize);
  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);
  vp = data;
  // UInt32 * res;
  for (int i=0;i<count;i++) {
    bool foundp;
    // *vp = (UInt32)res;
    ht.insert(vp,foundp);
    // printf("inserted ");
    // myprintv(vi,vsize);
    // printf(", foundp = %d\n",foundp);
    vp+= vsize;
  }
  getrusage(RUSAGE_SELF,&rue);
  double userTime,sysTime;
  reportTiming(rus,rue,userTime,sysTime,NULL);


  printf("%d inserts, size %d, lf %f, %d total cols, %d max cols, avg  %e, entries  %d, userTime  %f\n",
	 ht.numInserts,vsize,hash_abstract::loadFactor,
	 ht.numCollisions,
	 ht.maxCollisions,
	 (double)ht.numCollisions/(double)ht.numInserts,
	 ht.totalNumberEntries(),
	 userTime);

#if 0
  
  printf("max %d cols, %d total collisions with %d inserts, avg = %e, num entries = %d\n",
	 ht.maxCollisions,ht.numCollisions,ht.numInserts,
	 (double)ht.numCollisions/(double)ht.numInserts,
	 ht.totalNumberEntries());


  // now go through each entry in the table to test if it is there.
  printf("ensuring inserted entries are all contained\n");
  for (unsigned pos=0;pos<ht.table.size();pos++) {
    if (ht.table[pos].key != NULL) {
      // found an entry, make sure it is there.
      // printf("found at pos %d: ",pos);
      // myprintv(ht.table[pos].key,vsize);
      // printf("\n");
      if (!ht.find(ht.table[pos].key))
	printf("entry should be in hash table but isn't\n");
    }
  }

  // next go through each entry in table to make sure it
  // returns same pointer to itself.
  printf("ensuring pointers returned are proper\n");
  for (unsigned pos=0;pos<ht.table.size();pos++) {
    if (ht.table[pos].key != NULL) {
      // found an entry, make sure it is there.
      if (ht.table[pos].key != ht.insert(ht.table[pos].key))
	printf("entry should have same pointers");
    }
  }



  // finally, try some new random entries and count
  // the number of hits and mises.
  printf("checking if any random entries are contained\n");
  int nhits = 0;
  int *vi = new int[vsize];;
  for (int i=0;i<count;i++) {
    for (int j=0;j<vsize;j++) {
      vi[j] = rand() % maxCard;
    }
    if (ht.find(vi))
      nhits++;
  }
  printf("Found %d out of %d with random vectors\n",nhits,count);


  // next, try hashing some C++ strings.
  count = 100;
  printf("string hash table: inserting %d random strings of various sizes\n",count);
  unsigned size = (rand() % 30) + 1;
  vhash_set< char > sht(size);
  char *str = new char[size];;
  for (int i=0;i<count;i++) {
    for (unsigned j=0;j<size;j++) {
      str[j] = ' ' + (rand() % ('~' - ' '));
    }
    sht.insert(str);
  }
#endif

}

#endif
