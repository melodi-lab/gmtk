/*-
 * strhash_map.cc
 *     strhash_map driver class class
 *
 * Written by Gang Ji <gang@ee.washington.edu>
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
VCID(HGID);

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

