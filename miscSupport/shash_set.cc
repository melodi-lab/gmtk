/*-
 * vhash_set.cc
 *     shash_set driver class class
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * Modified by Gang Ji <gang@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */


#if HAVE_CONFIG_H
#include <config.h>
#endif
#include "hgstamp.h"
#include "general.h"
VCID(HGID)

#ifdef MAIN
#define COLLECT_COLLISION_STATISTICS
#include <string>
#endif

#include "shash_set.h"

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
	maxCard = atoi(argv[2]);

	printf("Using %d entries\n", count);

	shash_set< int > ht(3);

	printf("inserting %d random entries\n", count);
	for ( int i = 0; i < count; i++ ) {
		int vi = rand() % maxCard;
		bool foundp;
		ht.insert(vi, foundp);
	}

	printf("max %d cols, %d total collisions with %d inserts, avg = %e, num entries = %d\n",
		ht.maxCollisions,ht.numCollisions,ht.numInserts,
		(double)ht.numCollisions/(double)ht.numInserts,
		ht.totalNumberEntries());


	// now go through each entry in the table to test if it is there.
	printf("ensuring inserted entries are all contained\n");
	for ( unsigned pos = 0; pos < ht.table.size(); pos++ ) {
		if ( ht.table[pos].active ) {
			if ( ! ht.find(ht.table[pos].key) )
				printf("entry should be in hash table but isn't\n");
		}
	}

	// next go through each entry in table to make sure it
	// returns same pointer to itself.
	printf("ensuring pointers returned are proper\n");
	for ( unsigned pos = 0; pos < ht.table.size(); pos++ ) {
		if ( ht.table[pos].active ) {
			// found an entry, make sure it is there.
			if ( ht.table[pos].key != *ht.insert(ht.table[pos].key) )
				printf("entry should have same pointers");
		}
	}



	// finally, try some new random entries and count
	// the number of hits and mises.
	printf("checking if any random entries are contained\n");
	int nhits = 0;
	int vi;
	for ( int i = 0; i < count; i++ ) {
		vi = rand() % maxCard;
		if ( ht.find(vi) )
			nhits++;
	}
	printf("Found %d out of %d with random vectors\n", nhits, count);
}

#endif
