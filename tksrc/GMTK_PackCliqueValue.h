/*
 * GMTK_PackCliqueValue.h
 *   PackCliqueValue:
 *
 *     Class to support packing/unpacking of a vector of unsigned into
 *     a packed machine word representation.  It uses an array, the
 *     first portion for those elements that fit within a word
 *     boundary and the second portion for those that cross word
 *     boundaries.  This will obtain better branch prediction behavior
 *     (since otherwise, it'll be harder for the processor to
 *     accurately predict branch outcomes).
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


#include "debug.h"

// class mention for forward references.
class Partition;
class GMTemplate;
class RV;

class PackCliqueValue {

  friend class MaxClique;

  unsigned numUnsignedInPackedVector;

  unsigned totalNumBits;

  // @@@ should be const
  // const unsigned unpackedVectorLength;
  unsigned unpackedVectorLength;

  // Structure to hold information in going
  // from a packed word to an unpacked word.
  struct ValLocator {
    // the starting word in a packed vector
    unsigned start;
    // the amount to right-shift the packed word
    // to get an unpacked word
    unsigned startRightShift;
    // the mask in this first word
    unsigned startMask;

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

  // initialize
  void init(const unsigned *const cards);

public:

  PackCliqueValue(vector<RV*>& nodes);

  PackCliqueValue(const unsigned len, const unsigned *const cards); 

  // create an empty one for re-construction later
  PackCliqueValue() : numUnsignedInPackedVector(0),totalNumBits(0),unpackedVectorLength(0) {}

  ~PackCliqueValue() {}


  // Return the number of unsigned words that are required to hold a
  // packed clique value.
  unsigned packedLen() { return numUnsignedInPackedVector; }
  // Return the number of words needed to hold an unpacked clique
  // value (i.e., the number of hidden variables in the clique).
  unsigned unPackedLen() { return unpackedVectorLength; }
  // Return the number of bits needed to hold a packed clique value.
  unsigned packedLenBits() { return totalNumBits; }
  // Return the number of bytes needed to hold a packed clique value.
  unsigned packedLenBytes() { return (totalNumBits+7)/8; }


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
  void pack(const unsigned *const *const unpacked_vec,
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
	const unsigned val = (*unpacked_vec[*it]);
	packed_vec[vl_p->start] |= (val << vl_p->startRightShift);
	vl_p++;
      } while (++it != it_endp);
    }
    // next the ones that span word boundaries
    {
      const unsigned *it_endp = iterations.ptr+unpackedVectorLength;
      // assume that at least one case did not span a word boundary
      while (it != it_endp) {
	const unsigned val = (*unpacked_vec[*it]);
	packed_vec[vl_p->start] |= (val << vl_p->startRightShift);
	packed_vec[vl_p->start+1] |= 
	  ((val&(vl_p->nextMask<<vl_p->nextLeftShift)) >> vl_p->nextLeftShift);
	vl_p++;
	it++;
      }
    }
  }

  // version to do type change.
  // inline void pack(const int *const *const unpacked_vec,
  // int *const packed_vec) {
  // pack((const unsigned *const *const) unpacked_vec,
  //	 (unsigned*const)packed_vec);
  // }

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
  void unpack(const unsigned *const packed_vec,
	      unsigned **const unpacked_vec) {
    unsigned *it = iterations.ptr;
    ValLocator* vl_p = valLocators.ptr;
    // do the unpacking, first the ones
    // that do not span a word boundary
    {
      const unsigned *const it_endp = iterations.ptr+wordBoundaryOverlapLocation;
      do {
	(*unpacked_vec[*it]) = 
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
	(*unpacked_vec[*it]) = res;	
	vl_p++;
	it++;
      }
    }
  }


};


#endif

