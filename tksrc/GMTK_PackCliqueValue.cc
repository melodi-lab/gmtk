/*-
 * GMTK_PackCliqueValue
 *     Pack clique value
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

#include "GMTK_RandomVariable.h"
#include "GMTK_MaxClique.h"
#include "GMTK_PackCliqueValue.h"

VCID("$Header$");

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Support for building a tree from clique graph
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createPartitionJunctionTree()
 *   Create a mini-junction tree from the cliques in the given partition.
 *   This uses Kruskal's greedy (but optimal) algorithm for MST generation.
 *
 * Preconditions:
 *   The partition must be instantiated with cliques 
 *
 * Postconditions:
 *   The cliques in the partition are now such that they
 *   form a junction tree over cliques within that partition.
 *
 * Side Effects:
 *   Modifies all neighbors variables within the cliques within the
 *   partition.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */


PackCliqueValue::PackCliqueValue(const unsigned len, 
				 const unsigned *const cards)
  : unpackedVectorLength(len)
{

  // 1.0/log10(2) pre-defined, as a 64-bit double
  const double inv_log_2_base_e = 1.4426950408889633870E0;
 
  const unsigned numBitsPerUnsigned = sizeof(unsigned)*8;

  unsigned totalNumBits = 0;
  for (unsigned i = 0; i< len; i++) {
    totalNumBits += 
      (unsigned)ceil(log((double)cards[i])*inv_log_2_base_e);
  }

  numUnsignedInPackedVector = (unsigned)
    ceil((double)totalNumBits/(double)numBitsPerUnsigned);

  valLocators.resize(len);
  iterations.resize(len);

  unsigned curUnsignedLocation = 0;
  // number of unused bits in the current packed word
  unsigned curNumberUnusedBits = numBitsPerUnsigned;

  unsigned wordBoundaryNoOverlapLocation = 0;
  wordBoundaryOverlapLocation = len;
  for (unsigned i=0; i<len; i++) {

    unsigned curNumberBits = 
      (unsigned)ceil(log((double)cards[i])*inv_log_2_base_e);

    if (curNumberBits <= curNumberUnusedBits) {
      // use bits only in current word
      valLocators[wordBoundaryNoOverlapLocation].start = curUnsignedLocation;
      valLocators[wordBoundaryNoOverlapLocation].startRightShift = 
	(numBitsPerUnsigned - curNumberUnusedBits);
      valLocators[wordBoundaryNoOverlapLocation].startMask =
	((1 << curNumberBits)-1) << 
	valLocators[wordBoundaryNoOverlapLocation].startRightShift;

      iterations[wordBoundaryNoOverlapLocation] = i;

      curNumberUnusedBits -= curNumberBits;

      if (curNumberUnusedBits == 0) {
	curUnsignedLocation ++;
	curNumberUnusedBits = numBitsPerUnsigned;
      }
      wordBoundaryNoOverlapLocation++;
    } else {
      // use up remaining bits in this word

      wordBoundaryOverlapLocation--;
      
      valLocators[wordBoundaryOverlapLocation].start = curUnsignedLocation;
      valLocators[wordBoundaryOverlapLocation].startRightShift = 
	(numBitsPerUnsigned - curNumberUnusedBits);
      valLocators[wordBoundaryOverlapLocation].startMask =
	((1 << curNumberUnusedBits)-1) << 
	valLocators[wordBoundaryOverlapLocation].startRightShift;

      const unsigned numBitsRemaining = curNumberBits - curNumberUnusedBits;

      valLocators[wordBoundaryOverlapLocation].nextLeftShift = curNumberUnusedBits;
      valLocators[wordBoundaryOverlapLocation].nextMask = 
	((1 << numBitsRemaining)-1);

      iterations[wordBoundaryOverlapLocation] = i;

      curNumberUnusedBits = numBitsPerUnsigned - numBitsRemaining;
      curUnsignedLocation ++;

    }

  }
  assert( wordBoundaryNoOverlapLocation == wordBoundaryOverlapLocation);

} 


#ifdef MAIN


///////////////////////////////////////////
// main driver debugger for hash table.

#include <string>



int main(int argc,char*argv[])
{

  RAND myrnd(false);

  const unsigned numEpochs = 50;
  for (unsigned epoch=0;epoch<numEpochs;epoch++) {

    const unsigned len = myrnd.uniform(21)+1;
    sArray<unsigned> cards(len);
    for (unsigned i=0;i<len;i++) {
      cards[i] = myrnd.uniform(2,500000);
    }
    PackCliqueValue pcl(len,cards.ptr);
    const unsigned numExamples = 100000;

    printf("Epoch %d: Testing %d examples, len = %d, plen = %d, cards:",
	   epoch,numExamples,len,pcl.packedLen());
    for (unsigned i=0;i<len;i++) {
      printf(" %d",cards[i]);
    }
    printf("\n"); fflush(stdout);


    sArray<unsigned> vec(len);
    sArray<unsigned> packed_vec(pcl.packedLen());
    sArray<unsigned> unpacked_vec(len);
    for (unsigned ex=0;ex<numExamples;ex++) {

      for (unsigned i=0;i<len;i++) {
	vec.ptr[i] = myrnd.uniform(cards[i]-1);
      }

      // #define TIME_PACK

#ifdef TIME_PACK
      pcl.pack(vec.ptr,packed_vec.ptr);
      pcl.unpack(packed_vec.ptr,unpacked_vec.ptr);
#else
      for (unsigned i=0;i<len;i++) {
	unpacked_vec.ptr[i] = vec.ptr[i];
      }
#endif

      for (unsigned i=0;i<len;i++) {
	if (vec.ptr[i] != unpacked_vec.ptr[i])
	  error("ERROR: epoch %d, ex %d, location %d, initial packed %d and after packed %d, len=%d,plen=%d\n",epoch,ex,i,vec[i],unpacked_vec[i],len,pcl.packedLen()); 
      }
    }
  }

}


#endif



