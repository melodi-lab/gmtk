/*-
 * hash_abstract.cc
 *     abstract hash class
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


#include "general.h"
VCID("$Header$");

#ifdef MAIN
#define COLLECT_COLLISION_STATISTICS
#endif

#include "hash_map_list.h"

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

int
main(int argc,char*argv[])
{
  // first insert some numbers.
  unsigned count= 100000;
  unsigned maxCard = 100000;

  if (argc > 1)
    count = atoi(argv[1]);
  if (argc > 2)
    maxCard = atoi(argv[2]);

  hash_map_list < unsigned, double > hash_table(10);
  {
    printf("iterator over empty hash table\n");
    hash_map_list < unsigned, double >::iterator it;
    for (it=hash_table.begin();it != hash_table.end(); it++) {
      printf("hash table has %f\n",(*it));
    }
    for (unsigned i=0;i<101;i++) {
      double d = i;
      hash_table.insert(i,d);
    }
    printf("iterator after inserting elements\n");
    unsigned count=0;
    for (it=hash_table.begin();it != hash_table.end(); it++) {
      unsigned count2 = 2*count;
      printf("hash table iterator says %f, find(%d),find(%d) says %d,%d\n",(*it),count,2*count,(bool)hash_table.find(count),(bool)hash_table.find(count2));
      count++;
    } 
  }

  printf("clearing hash table\n");
  hash_table.clear();
  hash_table.clearStats();
  printf("inserting %d random entries\n",count);
  for (unsigned i=0;i<count;i++) {
    unsigned key = rand() % maxCard;
    double val = (double)key;
    hash_table.insert(key,val);
  }
  
  printf("max %d cols, %d total collisions with %d inserts, avg = %e, num entries = %d\n",
	 hash_table.maxCollisions,
	 hash_table.numCollisions,
	 hash_table.numInserts,
	 (double)hash_table.numCollisions/(double)hash_table.numInserts,
	 hash_table.totalNumberEntries());

  printf("Checking that all inserted elements are found\n");
  hash_map_list < unsigned, double >::iterator it;
  unsigned num_ells_it = 0;
  for (it=hash_table.begin();it != hash_table.end(); it++) {
    if (!hash_table.find(it.key()))
      printf("Didn't find key %d in hash table\n",it.key());
    else 
      num_ells_it ++;
  }
  printf("checked that %d elements from iterator are in hash table\n",num_ells_it); 

  // now print out a bunch of values
  for (unsigned i=0;i<100;i++) {
    unsigned key = rand() % maxCard;
    printf("hash_table[key=%d] gives find()=%d, val %f\n",
	   key,(bool)hash_table.find(key),hash_table[key]);
  }


}

#endif
