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

#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_MaxClique.h"
#include "GMTK_PackCliqueValue.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Support for finding good binary packing arrangement 
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


struct binAssignment {

  //////////////////////////
  // input variables.
  //////////////////////////

  // the location in the original array
  unsigned loc;
  // number of bits in this entry (i.e., the int requires val=ceil(log2(card)) bits)
  unsigned val;


  //////////////////////////
  // output variables.
  //////////////////////////

  // the bin number where this entry is assigned.
  unsigned bin;
  // if this is being split across bins 'bin' and 'bin+1' or not
  bool split;

  binAssignment() { split = false ; }
  bool operator<(const binAssignment& other) const { return val < other.val; } 

  struct descendByVal {
    bool operator()(const binAssignment& x, const binAssignment& y) { return (x.val) > (y.val); }
  };

  struct ascendByBin {
    bool operator()(const binAssignment& x, const binAssignment& y) {
       if (x.bin < y.bin) {
	 return true;
       } else if (x.bin == y.bin) {
	 if (  x.split == true  && x.split == y.split && x.loc != y.loc ) {
	   error("ERROR: bin packing: x.bin==y.bin = %d and both splits are true\n",x.bin);
	 }
	 // if y is split, then it goes at the end of the bin, it is considered greater.
	 return ( y.split ); 
       } else
	 return false;
    }
  };
};


void
printBinAssignmentArray(vector <binAssignment>& ba) 
{
  for (unsigned i=0;i<ba.size();i++) {
    printf("b=%d,s=%d,l=%d,v=%d, ",ba[i].bin,ba[i].split,ba[i].loc,ba[i].val);
  }
  printf("\n");
}



/*-
 *-----------------------------------------------------------------------
 * findNaive()
 *   find an assignmnet of clique values storage locations (enough bits) to 
 *   packed array locations. This uses a naive algorithm that is
 *   always guaranteed to succeed.
 *
 * Preconditions:
 *   bins must be of length so that all clique values can fit in packed
 *   array of this length, and ints must be initialized with the number
 *   of bits needed for each RV val, and the original array pointer.
 * 
 * Postconditions:
 *   the remaining members of elemnets of ints are filled in giving
 *   the locations of where packing should occur.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   number of splits (fewer is better)
 *
 *-----------------------------------------------------------------------
 */
unsigned
findNaive(vector<unsigned> bins,
	  vector<binAssignment>& ints)
{
  unsigned res = 0;
  for (unsigned i=0;i<ints.size();i++) {
    // find first non-empty bin.
    unsigned j =0;
    for (;j<bins.size();j++)
      if (bins[j] > 0)
	break;
    assert (j < bins.size());
    if (bins[j] >= ints[i].val) {
      bins[j] -= ints[i].val;
      ints[i].bin = j;
      ints[i].split = false;
    } else {
      // use up all of first slot.
      unsigned remainder = bins[j];
      bins[j] = 0;
      // assumption here is that an item is never
      // larger than a single bin.
      bins[j+1] -= (ints[i].val-remainder);
      // where does it go
      ints[i].bin = j;
      ints[i].split = true;

      res++;
    }
  }
  return res;
}


/*-
 *-----------------------------------------------------------------------
 * findApproximateBestRetry()
 *   find an assignmnet of clique values storage locations (enough
 *   bits) to packed array locations. This uses a smart algorithm that
 *   first sorts descending by bits and then finds the first entry
 *   where it fits.  This routine is not guaranteed to succeed (and
 *   retries with a permutation of the vals if it doesn't at
 *   first). If it fails, it returns ~0x0u.
 *
 * Preconditions:
 *   bins must be of length so that all clique values can fit in packed
 *   array of this length, and ints must be initialized with the number
 *   of bits needed for each RV val, and the original array pointer.
 * 
 * Postconditions:
 *   the remaining members of elemnets of ints are filled in giving
 *   the locations of where packing should occur.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   number of splits (fewer is better)
 *
 *-----------------------------------------------------------------------
 */
unsigned
findApproximateBestRetry(vector<unsigned> bins,
		    vector<binAssignment>& ints)
{
  sort(ints.begin(),ints.end(),binAssignment::descendByVal());

  unsigned secondTry = 0;
 top:

  unsigned res = 0;
  for (unsigned i=0;i<ints.size();i++) {
    // printf("ints[%d] = %d\n",i,ints[i].val);
    // first try to find a bin where it fits.
    bool found = false;
    for (unsigned j=0;j<bins.size();j++) {
      if (bins[j] >= ints[i].val) {
	bins[j] -= ints[i].val;

	ints[i].bin = j;
	ints[i].split = false;

	found = true;
	break;
      }
    }
    if (!found) {
      // need to do a split.
      for (unsigned j=0;(j+1)<bins.size();j++) {
	if (bins[j]+bins[j+1] >= ints[i].val) {
	  // use up all of first slot.
	  unsigned remainder = bins[j];
	  bins[j] = 0;
	  bins[j+1] -= (ints[i].val-remainder);

	  ints[i].bin = j;
	  ints[i].split = true;

	  found = true;
	  res++;
	  break;
	}
      }
    }
    if (!found) {
      if (secondTry == 10000)
	return ~0x0u;
      secondTry++;
      for (unsigned k=0;k<bins.size();k++)
	bins[k] = 8*sizeof(unsigned);
      random_shuffle(ints.begin(),ints.end());
      goto top;
    }

  }
  return res;
}


/*-
 *-----------------------------------------------------------------------
 * findApproximateBest2Retry()
 *   find an assignmnet of clique values storage locations (enough
 *   bits) to packed array locations. This uses a smart algorithm that
 *   first sorts descending by bits and then finds the best entry
 *   where it fits.  This routine is not guaranteed to succeed (and
 *   retries with a permutation of the vals if it doesn't at
 *   first).. If it fails, it returns ~0x0u.
 *
 * Preconditions:
 *   bins must be of length so that all clique values can fit in packed
 *   array of this length, and ints must be initialized with the number
 *   of bits needed for each RV val, and the original array pointer.
 * 
 * Postconditions:
 *   the remaining members of elemnets of ints are filled in giving
 *   the locations of where packing should occur.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   number of splits (fewer is better)
 *
 *-----------------------------------------------------------------------
 */
unsigned
findApproximateBest2Retry(vector<unsigned> bins,
		    vector<binAssignment>& ints)

{

  sort(ints.begin(),ints.end(),binAssignment::descendByVal());

  unsigned secondTry = 0;

 top:

  unsigned res = 0;
  for (unsigned i=0;i<ints.size();i++) {
    // printf("ints[%d] = %d\n",i,ints[i].val);
    // first try to find a bin where it best fits.
    unsigned bestIndex = bins.size();
    unsigned leastDiff = sizeof(unsigned)*8;
    bool found = false;
    for (unsigned j=0;j<bins.size();j++) {
      if (bins[j] >= ints[i].val) {
	unsigned diff = (bins[j] - ints[i].val);
	if (diff < leastDiff) {
	  leastDiff = diff;
	  bestIndex = j;
	}
      }
    }
    if (bestIndex < bins.size()) {
      // then found a place
      bins[bestIndex] -= ints[i].val;

      ints[i].bin = bestIndex;
      ints[i].split = false;

      found = true;
    } else {
      // need to do a split. Find a location
      // where it best fits in a split.
      leastDiff = 2*sizeof(unsigned)*8;
      for (unsigned j=0;(j+1)<bins.size();j++) {
	if (bins[j]+bins[j+1] >= ints[i].val) {
	  unsigned diff = (bins[j] + bins[j+1] - ints[i].val);
	  if (diff < leastDiff) {
	    leastDiff = diff;
	    bestIndex = j;
	  }
	}
      }
      if (bestIndex < bins.size()) {
	// then found a place
	unsigned remainder = bins[bestIndex];
	bins[bestIndex] = 0;
	bins[bestIndex+1] -= (ints[i].val-remainder);

	ints[i].bin = bestIndex;
	ints[i].split = true;

	found = true;
	res++;
      }
    }
    if (!found) {
      if (secondTry == 10000)
	return ~0x0u;
      secondTry++;
      for (unsigned k=0;k<bins.size();k++)
	bins[k] = 8*sizeof(unsigned);
      random_shuffle(ints.begin(),ints.end());
      goto top;

      // TODO: other heuristics to implement:
      //   1) don't give up so quickly. I.e., rather than immediately going for a split,
      //      if we don't find an entry, we shuffle and keep doing that a number of times
      //      looking for a perfect packing rather than immediately backing off to a split.
      //   2) rather than randomly_shuffling, do a random perturbation of the sorted entries,
      //      i.e., either sorted, or almost sorted is probably ok. We could randoml swap
      //      a few of the entries in the ints list each time. 

    }

  }
  return res;
}

/*-
 *-----------------------------------------------------------------------
 * PackCliqueValue::init()
 *   create an object that can be used to pack the clique values
 *   into a packed array of words.
 *
 * Preconditions:
 *   cards must be an length 'len' array of unsigned integers correspoinding
 *   to the cardinalities of the corresponding random variables
 *   in a clique.
 * 
 * Postconditions:
 *   Object is constructed. 
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
PackCliqueValue::init(const unsigned *const cards, const bool useNaive)
{
  // for easy access
  unsigned len = unpackedVectorLength;

  const unsigned numBitsPerUnsigned = sizeof(unsigned)*8;

  totalNumBits = 0;
  vector<binAssignment> ints(len);

  for (unsigned i = 0; i< len; i++) {
    unsigned bits = bitsRequiredUptoNotIncluding(cards[i]);
    totalNumBits += bits;
    ints[i].loc = i;
    ints[i].val = bits;
  }

  numUnsignedInPackedVector = 
    (totalNumBits+numBitsPerUnsigned-1)/numBitsPerUnsigned;

  vector<unsigned> bins(numUnsignedInPackedVector);

  for (unsigned i=0;i<numUnsignedInPackedVector;i++) {
    bins[i] = sizeof(unsigned)*8;
  }
#if 0
  // unused
  unsigned splits;
#endif

  if (useNaive || numUnsignedInPackedVector == 1) {
    // splits = findNaive(bins,ints);
    (void) findNaive(bins,ints);
  } else {
    vector<binAssignment> ints1 = ints;
    for (unsigned i=0;i<numUnsignedInPackedVector;i++) {
      bins[i] = sizeof(unsigned)*8;
    }
    unsigned tr1 = findApproximateBestRetry(bins,ints1);
    if (tr1 == 0) {
      ints = ints1;
      // splits = tr1;
    } else {
      vector<binAssignment> ints2 = ints;
      for (unsigned i=0;i<numUnsignedInPackedVector;i++) {
	bins[i] = sizeof(unsigned)*8;
      }
      unsigned tr2 = findApproximateBest2Retry(bins,ints2);
      if (tr2 == 0) {
	ints = ints2;
	// splits = tr2;
      } else {
	// try naive case as well.
	vector<binAssignment> ints_n = ints;	
	for (unsigned i=0;i<numUnsignedInPackedVector;i++) {
	  bins[i] = sizeof(unsigned)*8;
	}
	unsigned trn = findNaive(bins,ints_n);

	// printf("tr1=%d,tr2=%d,trn=%d\n",tr1,tr2,trn);
	
	// use one that has minimum.
	if (tr1 <= tr2) {
	  if (trn <= tr1) {
	    // use naive
	    ints = ints_n;
	    // splits = trn;
	  } else {
	    ints = ints1;
	    // splits = tr1;
	  }
	} else { // tr2  < tr1
	  if (trn <= tr2) {
	    // use naive
	    ints = ints_n;
	    // splits = trn;
	  } else {
	    ints = ints2;
	    // splits = tr2;
	  }
	}
      }
    }
  }
  // sort ascending by bin in order to:
  //  1) have close to unit stride during pack
  //  2) the below algorithm settingi shifts/masks will work.
  sort(ints.begin(),ints.end(),binAssignment::ascendByBin());

  // 
  // printBinAssignmentArray(ints);
  // printf("splits = %d\n",splits);

  valLocators.resize(len);
  valBits.resize(len);

  // number of unused bits in the current packed word
  unsigned curNumberUnusedBits = numBitsPerUnsigned;

  unsigned wordBoundaryNoOverlapLocation = 0;

  // index in 'iterations' where word-boundary overlaps occur.
  unsigned wordBoundaryOverlapLocation = len;

  for (unsigned i=0; i<len; i++) {

    const unsigned curNumberBits = ints[i].val;
    const unsigned curUnsignedLocation = ints[i].bin;

    if (ints[i].split == false) {
      assert (curNumberBits <= curNumberUnusedBits);

      valBits[wordBoundaryNoOverlapLocation] = curNumberBits;

      // use bits only in current word
      valLocators[wordBoundaryNoOverlapLocation].start = curUnsignedLocation;
      valLocators[wordBoundaryNoOverlapLocation].startRightShift = 
	(numBitsPerUnsigned - curNumberUnusedBits);
      valLocators[wordBoundaryNoOverlapLocation].startMask =
	((1 << curNumberBits)-1) << 
	valLocators[wordBoundaryNoOverlapLocation].startRightShift;

      valLocators[wordBoundaryNoOverlapLocation].loc = ints[i].loc;

      if (i+1 < len && ints[i].bin != ints[i+1].bin) {
	curNumberUnusedBits = numBitsPerUnsigned;
      } else 
	curNumberUnusedBits -= curNumberBits;

      wordBoundaryNoOverlapLocation++;

    } else {

      //if (curNumberBits <= curNumberUnusedBits) {
      // printf("i=%d, curNumberBits=%d,curNumberUnusedBits=%d\n",i,curNumberBits,curNumberUnusedBits);
      // fflush(stdout);
      // }

      assert (curNumberBits > curNumberUnusedBits);
      // use up remaining bits in this word

      wordBoundaryOverlapLocation--;

      valBits[wordBoundaryOverlapLocation] = curNumberBits;
      
      valLocators[wordBoundaryOverlapLocation].start = curUnsignedLocation;
      valLocators[wordBoundaryOverlapLocation].startRightShift = 
	(numBitsPerUnsigned - curNumberUnusedBits);
      valLocators[wordBoundaryOverlapLocation].startMask =
	((1 << curNumberUnusedBits)-1) << 
	valLocators[wordBoundaryOverlapLocation].startRightShift;
      
      // and use bits from next word

      const unsigned numBitsRemaining = curNumberBits - curNumberUnusedBits;

      valLocators[wordBoundaryOverlapLocation].nextLeftShift = curNumberUnusedBits;
      valLocators[wordBoundaryOverlapLocation].nextMask = 
	((1 << numBitsRemaining)-1);

      valLocators[wordBoundaryOverlapLocation].loc = ints[i].loc;

      curNumberUnusedBits = numBitsPerUnsigned - numBitsRemaining;

    }

  }
  assert( wordBoundaryNoOverlapLocation == wordBoundaryOverlapLocation);

  // lastly, initialize member pointer endpointers so 
  // that the need not be recomputed each time.
  member_vl_nwb_endp = valLocators.ptr+wordBoundaryOverlapLocation;
  member_vl_endp = valLocators.ptr+unpackedVectorLength;

}



/*-
 *-----------------------------------------------------------------------
 * PackCliqueValue::numWordsRequiredFor()
 * return the number of bits required for storing a vector of nodes
 * (or array with cardinalities).
 *
 * Preconditions:
 *   input must be valid
 * 
 * Postconditions:
 *   see above.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   number of bits n
 *
 *-----------------------------------------------------------------------
 */
unsigned PackCliqueValue::numWordsRequiredFor(vector<RV*>& nodes,RV* extra)
{
  if (nodes.size() == 0 && (extra == NULL)) return 0;
  sArray < unsigned > cards(nodes.size() + ((extra != NULL)?1:0));
  vector<RV*>::iterator it;
  unsigned i=0;
  for (it=nodes.begin();
       it != nodes.end();it++) {
    RV*rv = (*it);
    assert ( rv->discrete() );
    DiscRV*drv = RV2DRV(rv);
    cards.ptr[i] = drv->cardinality;
    i++;
  }
  if (extra != NULL) {
    RV*rv = extra;
    assert ( rv->discrete() );
    DiscRV*drv = RV2DRV(rv);
    cards.ptr[i] = drv->cardinality;
  }
  return numWordsRequiredFor(cards.size(),cards.ptr);
}
unsigned PackCliqueValue::numWordsRequiredFor(set<RV*>& nodes,RV* extra)
{
  if (nodes.size() == 0 && (extra == NULL)) return 0;
  sArray < unsigned > cards(nodes.size() + ((extra != NULL)?1:0));
  set<RV*>::iterator it;
  unsigned i=0;
  for (it=nodes.begin();
       it != nodes.end();it++) {
    RV*rv = (*it);
    assert ( rv->discrete() );
    DiscRV*drv = RV2DRV(rv);
    cards.ptr[i] = drv->cardinality;
    i++;
  }
  if (extra != NULL) {
    RV*rv = extra;
    assert ( rv->discrete() );
    DiscRV*drv = RV2DRV(rv);
    cards.ptr[i] = drv->cardinality;
  }
  return numWordsRequiredFor(cards.size(),cards.ptr);
}
unsigned PackCliqueValue::numWordsRequiredFor(const unsigned len, const unsigned *const cards)
{
  unsigned totalNumBits = 0;
  for (unsigned i = 0; i< len; i++) {
    unsigned bits = bitsRequiredUptoNotIncluding(cards[i]);
    totalNumBits += bits;
  }
  const unsigned numBitsPerUnsigned = sizeof(unsigned)*8;
  unsigned numUnsignedInPackedVector = 
    (totalNumBits+numBitsPerUnsigned-1)/numBitsPerUnsigned;
  return numUnsignedInPackedVector;
}




/*-
 *-----------------------------------------------------------------------
 * PackCliqueValue::PackCliqueValue()
 *   Construct an object that can be used to pack the clique values
 *   into a packed array of words. Uses an array of
 *   cardinalities
 *
 * Preconditions:
 *   see init()
 * 
 * Postconditions:
 *   see above.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
PackCliqueValue::PackCliqueValue(const unsigned len, 
				 const unsigned *const cards,
				 bool useNaive)
  : unpackedVectorLength(len)
{
  assert (unpackedVectorLength > 0);
  init(cards,useNaive);
}


/*-
 *-----------------------------------------------------------------------
 * PackCliqueValue::PackCliqueValue()
 *   Construct an object that can be used to pack the clique values
 *   into a packed array of words. Same as above version
 *   but works directly with cliques.
 *
 * Preconditions:
 *   see above constructor + nodes's variable must have
 *   been instantiated.
 * 
 * Postconditions:
 *   see above.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */

// PackCliqueValue::PackCliqueValue(MaxClique& maxClique)
//   : unpackedVectorLength(maxClique.nodes.size())
// {
//   assert (maxClique.nodes.size() != 0);

//   sArray < unsigned > cards(maxClique.nodes.size());

//   set<RV*>::iterator it;
//   unsigned i=0;
//   for (it=maxClique.nodes.begin();
//        it != maxClique.nodes.end();it++) {
//     RV*rv = (*it);
//     cards.ptr[i] = rv->cardinality;
//     i++;
//   }
//   init(cards.ptr);
// }

PackCliqueValue::PackCliqueValue(vector<RV*>& nodes)
  : unpackedVectorLength(nodes.size())
{


  assert (nodes.size() > 0);

  sArray < unsigned > cards(nodes.size());

  vector<RV*>::iterator it;
  unsigned i=0;
  for (it=nodes.begin();
       it != nodes.end();it++) {

    RV*rv = (*it);
    assert ( rv->discrete() );
    DiscRV*drv = RV2DRV(rv);
    cards.ptr[i] = drv->cardinality;
    i++;
  }
  init(cards.ptr);
}




#ifdef MAIN


///////////////////////////////////////////
// main driver debugger for hash table.

#include <string>



int main(int argc,char*argv[])
{

  RAND myrnd(true);
  printf("sizeof(PackCliqueValue::ValLocator) = %lu\n",sizeof(PackCliqueValue::ValLocator));

  const unsigned numEpochs = 50;
  for (unsigned epoch=0;epoch<numEpochs;epoch++) {

    const unsigned len = myrnd.uniform(1,22);
    sArray<unsigned> cards(len);
    for (unsigned i=0;i<len;i++) {
      cards[i] = myrnd.uniform(2,50000);
    }
    PackCliqueValue pcl(len,cards.ptr);
    const unsigned numExamples = 1000000;

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

#define TIME_PACK

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



