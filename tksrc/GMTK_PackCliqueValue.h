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

  typedef void (PackCliqueValue::*PackFunction)(const unsigned val,
						    unsigned* const,
						    ValLocator&);
  typedef unsigned (PackCliqueValue::*UnPackFunction)(const unsigned* const,
						      ValLocator&);


  void packToSingleWord(const unsigned val,
			unsigned* const pack_vec, 
			ValLocator& vl) 
  {
    pack_vec[vl.start] |= (val << vl.startRightShift);
  }

  void packToDoubleWord(const unsigned val,
			unsigned *const pack_vec, 
			ValLocator& vl)
  {
    // high order bits will fall off top of word so no need to 'and' here.
    pack_vec[vl.start]   |= (val << vl.startRightShift);
    pack_vec[vl.start+1] |= 
      ((val&(vl.nextMask<<vl.nextLeftShift)) >> vl.nextLeftShift);
  }


  unsigned unPackFromSingleWord(const unsigned* const pack_vec, 
				ValLocator& vl) 
  {
    return ( (pack_vec[vl.start] & vl.startMask) >> vl.startRightShift );
  }


  unsigned unPackFromDoubleWord(const unsigned* const pack_vec,
				ValLocator& vl) 
  {
    register unsigned res = 
      ( (pack_vec[vl.start] & vl.startMask) >> vl.startRightShift );
    res |= ( (pack_vec[vl.start+1] & vl.nextMask) << vl.nextLeftShift );
    return res;
  }

  sArray< ValLocator> valLocators;

  // array where lower locations hold indices for
  // which there is no word-boundary overlap and
  // upper locations hold indices for which there
  // are word-boundary overlaps.
  sArray< unsigned> iterations;

  // index in 'iterations' where word-boundary overlaps occur.
  unsigned wordBoundayOverlapLocation;


  sArray< PackFunction> packFunctions;
  sArray< UnPackFunction> unPackFunctions;

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
    // do the packing.
    const unsigned *vecp = unpacked_vec;
    const unsigned *vec_endp = unpacked_vec+unpackedVectorLength;
    ValLocator* vl_p = valLocators.ptr;
    PackFunction* pf_p = packFunctions.ptr;
    do {
      (this->*(*pf_p))(*vecp,packed_vec,*vl_p);
      pf_p++;
      vecp++;
      vl_p++;
    } while (vecp != vec_endp);
  }

  void unpack(const unsigned *const packed_vec,
	      unsigned *unpacked_vec) {
    unsigned *vecp = unpacked_vec;
    const unsigned *const vec_endp = unpacked_vec+unpackedVectorLength;
    ValLocator* vl_p = valLocators.ptr;
    UnPackFunction* upf_p = unPackFunctions.ptr;
    do {
      *vecp++ = (this->*(*upf_p))(packed_vec,*vl_p);
      upf_p++;
      vl_p++;
    } while (vecp != vec_endp);
  }


};


#endif

