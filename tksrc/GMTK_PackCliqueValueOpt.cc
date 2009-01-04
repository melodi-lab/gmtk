/*-
 * GMTK_PackCliqueValueOpt
 *     Pack clique value, optimized routine
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

#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_MaxClique.h"
#include "GMTK_PackCliqueValue.h"

VCID("$Header$")



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Support for distance routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*
 * returns the number of bits set in the unsigned
 *
 */
inline unsigned numBitsDifferent(unsigned long u1, unsigned long u2) 
{
  register unsigned count=0;
  register unsigned i = 8*sizeof(unsigned long);
  do {
    count += ( (u1&0x1) != (u2&0x1) );
    u1 >>=1; u2>>=1;
  } while (--i);
  return count;
}

inline unsigned hamdist(unsigned x, unsigned y)
{
  unsigned dist = 0, val = x ^ y;
  while(val)
  {
    ++dist; 
    val &= val - 1;
  }
  return dist;
}


unsigned 
PackCliqueValue::hamming_bit_distance(const unsigned *const packed_vec1,
				      const unsigned *const packed_vec2)
{
  register unsigned dist = 0;
  // No need to unpack since we are just doing bit distance.
  const unsigned *const packed_vec1_endp = packed_vec1 + numUnsignedInPackedVector;
  const unsigned* packed_vec1_p = packed_vec1;
  const unsigned* packed_vec2_p = packed_vec2;
  do {
    dist += hamdist(*packed_vec1_p,*packed_vec2_p);
    packed_vec1_p++;
    packed_vec2_p++;
  } while (packed_vec1_p != packed_vec1_endp);
  return dist;
}


unsigned 
PackCliqueValue::hamming_entry_distance(const unsigned *const packed_vec1,
					const unsigned *const packed_vec2)
{
  register ValLocator* vl_p = valLocators.ptr;
  register const ValLocator *vl_nwb_endp = member_vl_nwb_endp;
  register const ValLocator *vl_endp = member_vl_endp;
  register unsigned dist = 0;

  // do the unpacking, first the ones
  // that do not span a word boundary
  {
    do {
      register unsigned res1 =
	(packed_vec1[vl_p->start] & vl_p->startMask)
	>> vl_p->startRightShift;
      register unsigned res2 =
	(packed_vec2[vl_p->start] & vl_p->startMask)
	>> vl_p->startRightShift;
      vl_p++;

      dist += (res1 != res2);

    } while (vl_p != vl_nwb_endp);
  }
  // next the ones that do span a word boundary
  {
    while (vl_p != vl_endp) {
      register unsigned res1 =
	(packed_vec1[vl_p->start] & vl_p->startMask)
	>> vl_p->startRightShift;
      res1 |=
	((packed_vec1[vl_p->start+1] & vl_p->nextMask) <<
	 vl_p->nextLeftShift);

      register unsigned res2 =
	(packed_vec2[vl_p->start] & vl_p->startMask)
	>> vl_p->startRightShift;
      res2 |=
	((packed_vec2[vl_p->start+1] & vl_p->nextMask) <<
	 vl_p->nextLeftShift);

      dist += (res1 != res2);

      vl_p++;
    }
  }  
  return dist;
}

unsigned
PackCliqueValue::hamming_weighted_entry_distance(const unsigned *const packed_vec1,
						 const unsigned *const packed_vec2)
{

  register ValLocator* vl_p = valLocators.ptr;
  register const ValLocator *vl_nwb_endp = member_vl_nwb_endp;
  register const ValLocator *vl_endp = member_vl_endp;
  register unsigned dist = 0;
  register unsigned* valBits_p = valBits.ptr;

  // do the unpacking, first the ones
  // that do not span a word boundary
  {
    do {
      // we don't need to right shift since we are only testing equality below
      register unsigned res1 =
	(packed_vec1[vl_p->start] & vl_p->startMask);
      register unsigned res2 =
	(packed_vec2[vl_p->start] & vl_p->startMask);
      vl_p++;

      // dist += (res1 != res2)*valBits[i];
      if (res1 != res2) {
	// BP will predict that forward branch will not be taken, or condition is true.
	// It is more common that res1 will not equal res2.
	dist += (*valBits_p);
      } 

      valBits_p++;
    } while (vl_p != vl_nwb_endp);
  }
  // next the ones that do span a word boundary
  {
    while (vl_p != vl_endp) {
      // we again don't need to right shift since we are only testing
      // equality below.  (i.e., the assumption here is that the stuff
      // at the bottom of one word and the top of another word, when
      // or'd together, will not overlap).
      register unsigned res1 =
	(packed_vec1[vl_p->start] & vl_p->startMask);
      res1 |=
	(packed_vec1[vl_p->start+1] & vl_p->nextMask);

      register unsigned res2 =
	 (packed_vec2[vl_p->start] & vl_p->startMask);
      res2 |=
	(packed_vec2[vl_p->start+1] & vl_p->nextMask);

      // dist += (res1 != res2)*valBits[i];
      if (res1 != res2) {
	// BP will predict that forward branch will not be taken, or condition is true.
	// It is more common that res1 will not equal res2.
	dist += (*valBits_p);
      } 

      vl_p++;
      valBits_p++;
    }
  }
  return dist;
}



