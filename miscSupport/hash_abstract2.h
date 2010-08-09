/*
 * abstract_hash.h
 *   Generic stuff for hash tables 2.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * $Header$
 *
 */

#ifndef HASH_ABSTRACT2_H
#define HASH_ABSTRACT2_H

#include "debug.h"
#include "error.h"
#include "general.h"

class hash_abstract2 {

protected:

  // NOTE: this variable makes the code non-reentrant for threads
  static bool global_foundp;
  static const unsigned HashTableDefaultApproxStartingSize;

  //////////////////////////////////////////////////////////////////////////////
  // total number of entries currently inserted in the hash table
  // (distinct from the number that is allocated which is maintained locally)
  unsigned numberUniqueEntriesInserted;

  // Number of entries that when reached we should do a resize (kept
  // here so that we don't need to recompute using a divide, etc.).
  unsigned numEntriesToCauseResize;


#if defined(HASH_PRIME_SIZE)
  static const unsigned HashTable_SizePrimesArray;
  static const unsigned HashTable_PrimesArray[];
  ////////////////////////////////
  // the index into the above prime array of the current size.
  unsigned primesArrayIndex; 
  ////////////////////////////////
  // the initial starting index into the above prime 
  // array based on the contructors approximate starting
  // size.
  unsigned initialPrimesArrayIndex; 


  //////////////////////////////////
  // find the starting index in the array of primes starting
  // at the approximateStartingSize (i.e., finds the next prime
  // greater).
  void findPrimesArrayIndex(unsigned approximateStartingSize) {
    // do a linear search for now, but could do bin-search.
    for (primesArrayIndex=0;
	 primesArrayIndex<HashTable_SizePrimesArray;
	 primesArrayIndex++) {
      if (HashTable_PrimesArray[primesArrayIndex] >= approximateStartingSize)
	break;
    }
    if (primesArrayIndex == HashTable_SizePrimesArray)
      error("ERROR: hash_map_list can't create hash table of such a large approximate size %u\n",approximateStartingSize);
  }

#endif


protected:



/*
--------------------------------------------------------------------
mix, by Bob Jenkins, December 1996, Public Domain.
hash(), hash2(), hash3, and mix() are externally useful functions.
Routines to test the hash are included if SELF_TEST is defined.
You can use this free for any purpose.  It has no warranty.
--------------------------------------------------------------------
*/
#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}
  ////////////////////////////////////////////////////////////////////////////////////
  // Jenkens's hash function.
  ////////////////////////////////////////////////////////////////////////////////////
  /*
    --------------------------------------------------------------------
    by Bob Jenkins, December 1996, Public Domain.
    This works on all machines.  hash2() is identical to hash() on 
    little-endian machines, except that the length has to be measured
    in UInt32s instead of bytes.  It is much faster than hash().  It 
    requires
    -- that the key be an array of UInt32's, and
    -- that all your machines have the same endianness, and
    -- that the length be the number of UInt32's in the key
    --------------------------------------------------------------------
  */
  inline UInt32 hash_jenkens( register UInt32 *k,  /* the key */
			      register unsigned ksize, /* the length of the key, in UInt32s */
			      register UInt32 initval /* the previous hash, or an arbitrary value */
			      )
  {
    register UInt32 a,b,c,len;

    /* Set up the internal state */
    len = ksize;
    a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
    c = initval;           /* the previous hash value */

    /*---------------------------------------- handle most of the key */
    while (len >= 3)
      {
	a += k[0];
	b += k[1];
	c += k[2];
	mix(a,b,c);
	k += 3; len -= 3;
      }

    /*-------------------------------------- handle the last 2 UInt32's */
    c += ksize;
    switch(len)              /* all the case statements fall through */
      {
	/* c is reserved for the ksize */
      case 2 : b+=k[1];
      case 1 : a+=k[0];
	/* case 0: nothing left to add */
      }
    mix(a,b,c);
    /*-------------------------------------------- report the result */
    return c;
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // FNV style hash function.
  ////////////////////////////////////////////////////////////////////////////////////
  inline UInt32 hash_fnv1( register UInt32 * k,
			   register unsigned ksize,
			   register UInt32 hval )
  {
    register UInt32* ke= k+ksize;
    do {
      hval ^= *k++;
      // hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24);
      hval *= 16777619; // 32 bit prime
    } while (k!=ke);
    return hval;
  }
  inline UInt32 hash_fnv2 ( register UInt32 * k,
			    register unsigned ksize,
			    register UInt32  hval )
  {
    register UInt32*ke=k+ksize;
    do {
      // hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24);
      hval *= 16777619; // 32 bit prime
      hval ^= *k++;
    } while (k!=ke);
    return hval;
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // GMTK-A hash function.
  ////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // since this is a double hash table, we define two address
  // functions, h1() and h2().  h1() gives the starting position of in
  // the array of the key, and h2() gives the increment when we have a
  // collision.
  inline unsigned hash_gmtk_a1(const UInt32* key, 
			       const unsigned ksize)
  {
    const UInt32* keyp = key + ksize;

    unsigned long a = ksize;
    // unsigned long a = 0;
    do {
      --keyp;
      // a =65599*a + (*keyp) + 1;
      // a = 402223*a + (*keyp) + 1;
      // a = 611953*a + (*keyp) + 1;
      a = 3367900313ul*a + (*keyp) + 1;
    } while (keyp != key);
    return a;
  }
  ///////////////////////////////////////////////////////////////////////
  // h2() is the increment of of key when we have a collision.  "The
  // value of h2(key) must be relatively prime to the hash-table m for
  // the entire hash table to be searched" (from Corman, Leiserson,
  // Rivest). Therefore, we have result of the h2() function satisfy
  //         1) it must be greater than zero 
  //         2) it must be strictly less than table.size()
  // 
  inline unsigned hash_gmtk_a2(const UInt32* key, 
			       const unsigned ksize, 
			       const unsigned start = 0)
  {

    unsigned long a=start;
    const UInt32* keyp = key + ksize;
    do {
      --keyp;
      // a =65599*a + (*keyp) + 1;
      // a = 402223*a + (*keyp) + 1;
      // a = 611953*a + (*keyp) + 1;
      // a = 1500450271ul*a + (*keyp) + 1;
      a = 3267000013ul*a + (*keyp) + 1;
    } while (keyp != key);
    return a;
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // GMTK-B hash function.
  ////////////////////////////////////////////////////////////////////////////////////
  inline unsigned hash_gmtk_b1(UInt32* key, 
			       const unsigned ksize)
  {
    register unsigned long a = 0;
    const UInt32* keyp = key + ksize;
    do {
      a += (a <<3) + (*key++);
    } while (keyp != key);
    return a;
  }
  inline unsigned hash_gmtk_b2(UInt32* key, 
			       const unsigned ksize, 
			       const unsigned start = 0)
  {
    register unsigned long a=start;
    const UInt32* keyp = key + ksize;
    do {
      a += (a <<3) + (*key++);
    } while (keyp != key);
    return a;
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // GMTK-C hash function.
  ////////////////////////////////////////////////////////////////////////////////////
  inline unsigned hash_gmtk_c1(UInt32* key, const unsigned ksize)
  {
    register unsigned long a = ksize;
    const UInt32* keyp = key + ksize;
    do {
      a += (a <<3) + (a>>(8*sizeof(unsigned)-3)) + (*key++);
    } while (keyp != key);
    return a;
  }
  inline unsigned hash_gmtk_c2(UInt32* key, 
			       const unsigned ksize, 
			       const unsigned start = 0)
  {
    register unsigned long a=start;
    const UInt32* keyp = key + ksize;
    do {
      a += (a <<3) + (a>>(8*sizeof(unsigned)-3)) + (*key++);
    } while (keyp != key);
    return a;
  }


  ////////////////////////////////////////////////////////////////////////////////////
  // GMTK-D hash function.
  ////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  inline unsigned hash_gmtk_d1(Uint32* key, const unsigned ksize)
  {
    const UInt32* keyp = key + ksize;
    unsigned long a = ksize;
    do {
      a = 3367900313ul*a + (*key++);
    } while (keyp != key);
    return a;
  }
  ///////////////////////////////////////////////////////////////////////
  inline unsigned hash_gmtk_d2(UInt32* key, const unsigned ksize, const unsigned start = 0)
  {

    unsigned long a=start;
    const UInt32* keyp = key + ksize;
    do {
      a = 3267000013ul*a + (*key++);
    } while (keyp != key);
    return a;
  }



  ////////////////////////////////////////////////////////////////////////////////////
  // GMTK-E hash function.
  ////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  inline unsigned hash_gmtk_e1(UInt32* key, const unsigned ksize)
  {
    const UInt32* keyp = key + ksize;
    unsigned long a = ksize;
    do {
      a = 3367900313ul*a + (*key++);
    } while (keyp != key);
    return a;
  }
  ///////////////////////////////////////////////////////////////////////
  inline unsigned hash_gmtk_e2(UInt32* key, const unsigned ksize, const unsigned start = 0)
  {

    unsigned long a=start;
    const UInt32* keyp = key + ksize;
    do {
      a += (a <<3) + (a>>(8*sizeof(unsigned)-3)) + (*key++);
    } while (keyp != key);
    return a;
  }


#ifndef HASH_PRIME_SIZE
   // mask to 'and' with to get valid hash key.
   unsigned sizeMask;
#endif


public:

  hash_abstract2() : numberUniqueEntriesInserted(0) {}

  ////////////////////////////////////////////////////////////////
  // return the total number of entries that have been inserted into
  // the hash table (this is different than the total allocated size
  // of any table).
  unsigned totalNumberEntries() { return numberUniqueEntriesInserted; }

  // how full does the hash table need to be before a resize occurs.
  static float loadFactor;


 #ifdef COLLECT_COLLISION_STATISTICS
   unsigned maxCollisions;
   unsigned numCollisions;
   unsigned numInserts;
 #endif


   //////////////////////////////////////////////////////////////////
   // return the entry of key in table a_table
  template <class TableType>
  inline unsigned entryOf(const UInt32* key,
		   const unsigned ksize,
		   TableType& a_table)
  {

 #if defined(HASH_FNV)
     unsigned aa = hash_fnv1((UInt32*)key,ksize,0);
 #elif defined(HASH_JENKENS)
     unsigned aa = hash_jenkens((UInt32*)key,ksize,0);
 #elif defined(HASH_GMTK_A)
     unsigned aa = hash_gmtk_a1((UInt32*)key,ksize);
 #elif defined(HASH_GMTK_B)
     unsigned aa = hash_gmtk_b1((UInt32*)key,ksize);
 #elif defined(HASH_GMTK_C)
     unsigned aa = hash_gmtk_c1((UInt32*)key,ksize);
 #elif defined(HASH_GMTK_D)
     unsigned aa = hash_gmtk_d1((UInt32*)key,ksize);
 #elif defined(HASH_GMTK_E)
     unsigned aa = hash_gmtk_e1((UInt32*)key,ksize);
 #else /* use several */
     unsigned aa = hash_fnv1((UInt32*)key,ksize,hash_jenkens((UInt32*)key,ksize,hash_gmtk_a1((UInt32*)key,ksize)));
 #endif

#if defined(HASH_PRIME_SIZE)
     unsigned a = aa % a_table.size();
#elif defined(HASH_LOC_FOLD)
     // do a fast 16-bit fold for all occasions. Ideally, we'd fold according to current table size.
     unsigned a = (aa ^ (aa >> 16)) & sizeMask;
#else
     unsigned a = aa & sizeMask;
#endif


 #ifdef COLLECT_COLLISION_STATISTICS
     unsigned collisions=0;
 #endif


     // printf("entryOf: size=%d, ksize=%d, sizeMask = 0x%X, key= ",size,ksize, sizeMask);
     // for (unsigned j=0;j<ksize;j++)
     // printf("%X ",key[j]);
     // printf("a1=%d, empty = %d,",a,empty(a_table[a]));
     // if (!empty(a_table[a]))
     // printf(" ,ke=%d",keyEqual(a_table[a].key,key));
     // printf("\n");
     
     const UInt32* const key_endp = key+ksize;
     if (a_table.ptr[a].empty() || a_table.ptr[a].keyEqual(key,key_endp)) {
       return a;
     }

 #if defined(HASH_FNV)
     unsigned inc = (hash_fnv2((UInt32*)key,ksize,aa));
 #elif defined(HASH_JENKENS)
     unsigned inc = (hash_jenkens((UInt32*)key,ksize,aa));
 #elif defined(HASH_GMTK_A)
     unsigned inc = (hash_gmtk_a2((UInt32*)key,ksize,aa));
 #elif defined(HASH_GMTK_B)
     unsigned inc = (hash_gmtk_b2((UInt32*)key,ksize,aa));
 #elif defined(HASH_GMTK_C)
     unsigned inc = (hash_gmtk_c2((UInt32*)key,ksize,aa));
 #elif defined(HASH_GMTK_D)
     unsigned inc = (hash_gmtk_d2((UInt32*)key,ksize,aa));
 #elif defined(HASH_GMTK_E)
     unsigned inc = (hash_gmtk_e2((UInt32*)key,ksize,aa));
 #else
     unsigned inc = (hash_fnv2((UInt32*)key,ksize,hash_jenkens((UInt32*)key,ksize,hash_gmtk_a2((UInt32*)key,ksize,aa))));
 #endif

     // The value of inc must be relatively prime to the hash-table m
     // for the entire hash table to be searched (see Corman,
     // Leiserson, Rivest).
#if defined(HASH_PRIME_SIZE)
     // Therefore, we must have that inc satisfy
     //    1) it must be greater than zero 
     //    2) it must be strictly less than table.size()
     inc = 
       ((inc % (a_table.size()-1)) // this gives [0 : (a_table.size()-2)]
	+
	1                      // this gives [1 : (a_table.size()-1)] 
	);
#else
     // For a power of 2 sized table, we make sure that the increment
     // is non-zero and odd.
     inc |= 0x1;
#endif

     do {

#ifdef COLLECT_COLLISION_STATISTICS
      collisions++;
#endif

#if defined(HASH_PRIME_SIZE)
      a = (a+inc) % a_table.size();
#else
      a = (a+inc) & sizeMask;
#endif

      // printf("entryOf: C. inc=%u,now at entry %d\n",inc,a);
    } while ( !a_table.ptr[a].empty()
	      &&
	      (!a_table.ptr[a].keyEqual(key,key_endp) )
	      );

#ifdef COLLECT_COLLISION_STATISTICS
    if (collisions > maxCollisions)
      maxCollisions = collisions;
    numCollisions += collisions;
#endif

    return a;
  }

};


#endif // defined HASH_ABSTRACT2

