/*-
 * strhash_map.cc
 *     strhash_map driver class class
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
VCID("$Header$");

#ifdef MAIN
#define COLLECT_COLLISION_STATISTICS
#include <string>
#endif

#include "strhash_map.h"

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
	int count= 10000;

	printf("Using %d entries\n", count);

	strhash_map<std::string> ht(3);

	printf("inserting %d random entries\n", count);
	// from 32 to 126
	char key[20];
	for ( int i = 0; i < count; i++ ) {
		int keylen = rand() % 10 + 8;
		for ( int j = 0; j < keylen; j++ ) {
			key[j] = rand() % ('z' - '/') + '/';
		}
		key[keylen] = '\0';
		std::string data = key;

		ht.insert(key, data);
		//printf("inserted %s\n", key);
	}

	//ht.print();
	printf("total size: %d, table size: %d\n", ht.size(), ht.tableSize());
}

#endif

