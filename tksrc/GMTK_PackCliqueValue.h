/*
 * GMTK_PackCliqueValue.h
 *   PackCliqueValue:
 *
 *     Class to support packing/unpacking of a vector of unsigned into
 *     a packed machine word representation.  It uses an array of
 *     pointers to functions to do the packing under the assumption
 *     that the branch-always-taken condition for a routine call will
 *     be cheaper than using the if-then branch but which will be less
 *     predictable for the branch prediction logic.
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

#ifndef GMTK_PACK_CLIQUE_VALUE_H
#define GMTK_PACK_CLIQUE_VALUE_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_RandomVariable.h"
#include "GMTK_MaxClique.h"

#include "debug.h"

// class mention for forward references.
class Partition;
class GMTemplate;

class PackCliqueValue {

  friend class MaxClique;


  unsigned numUnsignedInPackedVector;

  const unsigned unpackedVectorLength;

  // Structure to hold information in going
  // from a packed word to an unpacked word.
  struct ValLocator {
    // the starting word in a packed vector
    unsigned start;
    // the mask in this first word
    unsigned startMask;
    // the amount to right-shift the packed word
    // to get an unpacked word
    unsigned startRightShift;

    // if we overlap to next machine word, the mask
    // in that next word, these bits are to be or'd in
    // the higher order bits of the unpacked word.
    unsigned nextMask;
    // the amount to left shift these masked bits to
    // form the high order bits of the unpacked word.
    unsigned nextLeftShift;
  };


  sArray< ValLocator> valLocators;

  // array where lower locations hold indices for
  // which there is no word-boundary overlap and
  // upper locations hold indices for which there
  // are word-boundary overlaps.
  sArray< unsigned> iterations;

  // index in 'iterations' where word-boundary overlaps occur.
  unsigned wordBoundaryOverlapLocation;


public:

  PackCliqueValue(MaxClique& maxClique);

  PackCliqueValue(const unsigned len, const unsigned *const cards); 

  ~PackCliqueValue() {}


  // return the number of unsigneds that
  // hold a packed clique value
  unsigned packedLen() { return numUnsignedInPackedVector; }
  unsigned unPackedLen() { return unpackedVectorLength; }

  // pack()
  // this routine assumes that both
  //         unpackedVectorLength > 0 
  // and that
  //      packedVectorLength > 0.
  void pack(const unsigned *const unpacked_vec,
	    unsigned *const packed_vec) {
    // zero out packed vector
    unsigned *packed_vecp = packed_vec;
    const unsigned *const packed_vec_endp = 
      packed_vec + numUnsignedInPackedVector;
    do {
      *packed_vecp++ = 0;
    } while (packed_vecp != packed_vec_endp);

    // do the packing, first the ones that do not
    // span a word boundaries
    const unsigned *it = iterations.ptr;
    ValLocator* vl_p = valLocators.ptr;
    {
      const unsigned *it_endp = iterations.ptr+wordBoundaryOverlapLocation;
      // assume that at least one case did not span a word boundary
      do {
	const unsigned val = unpacked_vec[*it];
	packed_vec[vl_p->start] |= (val << vl_p->startRightShift);
	vl_p++;
      } while (++it != it_endp);
    }
    // next the ones that span word boundaries
    {
      const unsigned *it_endp = iterations.ptr+unpackedVectorLength;
      // assume that at least one case did not span a word boundary
      while (it != it_endp) {
	const unsigned val = unpacked_vec[*it];
	packed_vec[vl_p->start] |= (val << vl_p->startRightShift);
	packed_vec[vl_p->start+1] |= 
	  ((val&(vl_p->nextMask<<vl_p->nextLeftShift)) >> vl_p->nextLeftShift);
	vl_p++;
	it++;
      }
    }
  }

  // same as above, but that packs from an array of pointers to ints
  void pack(const unsigned **const unpacked_vec,
	    unsigned *const packed_vec);


  void unpack(const unsigned *const packed_vec,
	      unsigned *const unpacked_vec) {
    unsigned *it = iterations.ptr;
    ValLocator* vl_p = valLocators.ptr;
    // do the unpacking, first the ones
    // that do not span a word boundary
    {
      const unsigned *const it_endp = iterations.ptr+wordBoundaryOverlapLocation;
      do {
	unpacked_vec[*it] = 
	  (packed_vec[vl_p->start] & vl_p->startMask)
	    >> vl_p->startRightShift;
	vl_p++;
	it++;
      } while (it != it_endp);
    }
    // next the ones that do span a word boundary
    {
      const unsigned *const it_endp = iterations.ptr+unpackedVectorLength;
      while (it != it_endp) {
	register unsigned res =
	  (packed_vec[vl_p->start] & vl_p->startMask)
	    >> vl_p->startRightShift;
	res |=
	  ((packed_vec[vl_p->start+1] & vl_p->nextMask) <<
	   vl_p->nextLeftShift);
	unpacked_vec[*it] = res;	
	vl_p++;
	it++;
      }
    }
  }

  // same as above, but version that unpacks to array of pointers to
  // ints.
  void unpack(const unsigned **const packed_vec,
	      unsigned *const unpacked_vec);





};


#endif

