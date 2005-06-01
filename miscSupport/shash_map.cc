/*-
 * shash_map.cc
 *     shash_map driver class class
 *
 * Written by Gang Ji <gang@ee.washington.edu>
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

#include "shash_map.h"

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

int main(int argc,char*argv[]) {
	// first insert some numbers.
	int count= 100000;
	int maxCard = 100000;

	if ( argc > 1 )
		count = atoi(argv[1]);

	if ( argc > 2 )
		maxCard = atoi(argv[3]);

	printf("Using %d entries\n", count);

	shash_map< int, double > ht(3);

	printf("inserting %d random entries\n", count);
	for ( int i = 0; i < count; i++ ) {
		int key = rand() % maxCard;
		double data = key;
		bool foundp;

		ht.insert(key, data, foundp);
		// printf("inserted %d, foundp = %d\n", key, foundp);
	}

	printf("max %d cols, %d total collisions with %d inserts, avg = %e, num entries = %d\n",
		ht.maxCollisions, ht.numCollisions, ht.numInserts,
		(double)ht.numCollisions / (double)ht.numInserts,
		ht.totalNumberEntries());


	// now go through each entry in the table to test if it is there.
	printf("ensuring inserted entries are all contained\n");
	for ( unsigned pos = 0; pos < ht.table.size(); pos++ ) {
		if ( ht.table[pos].active ) {
			// found an entry, make sure it is there.
			// printf("found at pos %d: ",pos);
			// myprintv(ht.table[pos].key,vsize);
			// printf("\n");
			double *dp = ht.find(ht.table[pos].key);
			if ( dp == NULL )
				printf("ERROR: entry should be in hash table but isn't\n");
		}
	}

	// next go through each entry in table to make sure it
	// returns same pointer to itself.
	printf("ensuring data values returned are proper\n");
	for ( unsigned pos = 0; pos < ht.table.size(); pos++) {
		if ( ht.table[pos].active ) {
			// found an entry, make sure it is there.
			double d = ht.table[pos].key;
			// double d2 = ht[ht.table[pos].key];
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
		int d = rand() % maxCard;
		if ( ht.find(d) )
			nhits++;
	}
	printf("Found %d out of %d with random vectors\n",nhits,count);
}

#endif
