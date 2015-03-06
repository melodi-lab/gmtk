/*-
 * hash_abstract.cc
 *     abstract hash class
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
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
#endif

#include "hash_abstract.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


// NOTE: this variable makes the code non-reentrant for threads
bool hash_abstract::global_foundp = false;

const unsigned
hash_abstract::HashTableDefaultApproxStartingSize = 20000;

float
hash_abstract::loadFactor = 0.5;


#if defined(HASH_PRIME_SIZE)

/*
 * Ten random 10 digit primes
 *
 * unsigned long tenDigitPrimes[] = {
 * 5915587277u, 
 * 1500450271u, 
 * 3267000013u, 
 * 5754853343u, 
 * 4093082899u, 
 * 9576890767u, 
 * 3628273133u, 
 * 2860486313u, 
 * 5463458053u, 
 * 3367900313u };
 * 
 */

/*
other primes

7, 13, 31, 61, 127, 251, 509, 1021, 2017, 4093,
  5987, 9551, 15683, 19609, 31397,
  65521, 131071, 262139, 524287, 1048573, 2097143,
  4194301, 8388593, 16777213, 33554393, 67108859,
  134217689, 268435399, 536870909, 1073741789
*/

/*

  1,
  2,
  3,
  7,
  13,
  23,
  59,
  113,
  241,
  503,
  1019,
  2039,
  4091,
  8179,
  11587,
  16369,
  23143,
  32749,
  46349,
  65521,
  92683,
  131063,
  185363,
  262139,
  330287,
  416147,
  524269,
  660557,
  832253,
  1048571,
  1321109,
  1664501,
  2097143,
  2642201,
  3328979,
  4194287,
  5284393,
  6657919,
  8388593,
  10568797,
  13315831,
  16777199,
  33554393,
  67108859,
  134217689,
  268435399,
  536870879,
  1073741789,
  2147483629
  };
*/

//
// the possible sizes for the has table, a list of prime numbers each
// entry roughly doubling in size of the previous previous entry (but
// as the sizes get larger, we increase by smaller amounts, see
// comments below). Comment out the entries at the beginning of the
// table to get the desired starting size of the hash table.
// 

const unsigned
hash_abstract::HashTable_PrimesArray[] = {
3,
5,
11,
23,
53,
107,
211,
421,
743,
1531,
3709,
6317,
12437,
24547,
47779,
75403,
134053,
262139,
524287,
1048573, // after a million, we start increasing by roughly powers of sqrt(2)
1664501,
2097143,
3328979,
4194301,
6657919,
8388593, // after a 8 million, we start increasing by roughly powers of 2^(1/3)
10568981,
13316071,
16777213,
21137981,
26632181,
33554393,
42275897,
53264317,
67108859,  // after 64 million, we start increasing by roughly powers of 2^(1/4)
79806341,
94906297,
112863197,
134217689,
159612647,
189812477,
225726349,
268435399,
319225297,
379624981,
451452737,
536870909,
638450719,
759250133,
902905657,
1073741789, 
1276901389,
1518500213,
1805811251,
2147483647,
/* 4294967291L */
((unsigned long) 2147483647) + ((unsigned long) 2147483642)
};

/////////////////////////////////////////////////////////////
// the size of the above table, determined automatically.
const unsigned 
hash_abstract::HashTable_SizePrimesArray
   = sizeof(HashTable_PrimesArray)/sizeof(unsigned);

#endif // defined(HASH_PRIME_SIZE)
