/*-
 * vshash_map2.cc
 *     vshash_map2 driver class class
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
VCID("$Header$")

#ifdef MAIN
#define COLLECT_COLLISION_STATISTICS
#include <string>
#endif

#include "vshash_map.h"

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

void myprintv(int* v,int len) {
	for (int i=0;i<len;i++) {
		printf("%d ",v[i]);
	}
}


int main(int argc,char*argv[]) {
	// first insert some numbers.
	int count= 100000;
	int vsize = 3;
	int maxCard = 100000;

	if ( argc > 1 )
		count = atoi(argv[1]);
	if ( argc > 2 )
		vsize = atoi(argv[2]);
	if ( argc > 3 )
		maxCard = atoi(argv[3]);

	printf("Using %d entries, each of length %d\n", count, vsize);

	vshash_map<int, double> ht(vsize, 3);

	printf("inserting %d random entries of size %d\n", count, vsize);
	int* vi = new int[vsize];
	for ( int i = 0; i < count; i++ ) {
		double d = 0;
		for ( int j = 0; j < vsize; j++ ) {
			vi[j] = rand() % maxCard;
			d = maxCard*d + vi[j];
		}
		bool foundp;
		ht.insert(vi, d, foundp);
	}

	printf("max %d cols, %d total collisions with %d inserts, avg = %e, num entries = %d\n",
		ht.maxCollisions,ht.numCollisions,ht.numInserts,
		(double)ht.numCollisions/(double)ht.numInserts,
		ht. totalNumberEntries());


	// now go through each entry in the table to test if it is there.
	printf("ensuring inserted entries are all contained\n");
	for ( unsigned pos = 0; pos < ht.table.size(); pos++ ) {
		if ( ht.table[pos].key != NULL ) {
			// found an entry, make sure it is there.
			double *dp = ht.find(ht.table[pos].key);
			if ( dp == NULL )
				printf("ERROR: entry should be in hash table but isn't\n");
		}
	}

	// next go through each entry in table to make sure it
	// returns same pointer to itself.
	printf("ensuring data values returned are proper\n");
	for ( unsigned pos = 0; pos < ht.table.size(); pos++ ) {
		if ( ht.table[pos].key != NULL ) {
			// found an entry, make sure it is there.
			double d = 0;
			for ( int j = 0; j < vsize; j++ ) {
				d = maxCard*d + ht.table[pos].key[j];
			}
			double* d2 = ht.find(ht.table[pos].key);
			if ( d2 == NULL )
				printf("ERROR: didnt find entry when should\n");
			if ( d != *d2 )
				printf("ERROR: entry should have same data values, d=%f,d2=%f\n",d,*d2);
		}
	}

	// finally, try some new random entries and count
	// the number of hits and mises.
	printf("checking if any random entries are contained\n");
	fflush(stdout);
	int nhits = 0;
	for ( int i = 0; i < count; i++ ) {
		double d;
		for ( int j = 0; j < vsize; j++ ) {
			vi[j] = rand() % maxCard;
			d = maxCard*d + vi[j];
		}
		if ( ht.find(vi) )
			nhits++;
	}
	printf("Found %d out of %d with random vectors\n",nhits,count);

	delete [] vi;
	vi = NULL;
}

#endif
