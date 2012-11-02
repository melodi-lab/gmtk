
/*-
 * GMTK_JunctionTree.cc
 *     Junction Tree, message passing routines.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003-2009.
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
#include <new>
#include <typeinfo>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_INTTYPES_H
   // The ISO C99 standard specifies that the macros in inttypes.h must
   //  only be defined if explicitly requested. 
#  define __STDC_FORMAT_MACROS 1
#  include <inttypes.h>
#endif
#if HAVE_STDINT_H
#  include <stdint.h>
#endif

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_GMParms.h"
#include "GMTK_Dlinks.h"

#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_FileSource.h"
#  include "GMTK_StreamSource.h"
#endif

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)


/*
 *  init_CC_CE_rvs(it)
 *    given an iterator representing an unrolled segment,
 *    we create the pair of random variables CC and CE that
 *    will be shifted around in time/position in order
 *    to do inference.
 *
 * Preconditions:
 *    1) partitionStructureArray must be set up and initialized before
 *    this routine is called.
 *    2) All random variables should currently be shifted to their 
 *       zero position (otherwise thigns will be totally out of sync).
 *
 */
void
JunctionTree::init_CC_CE_rvs(ptps_iterator& ptps_it)
{
  if (ptps_it.num_c_partitions() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // partitions (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = partitionStructureArray[1].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,cur_CC_rvs);
    tmp1 = partitionStructureArray[3].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,cur_CE_rvs);
  } else {
    cur_CC_rvs.clear();
    cur_CE_rvs.clear();
  } 
  cur_cc_shift = cur_ce_shift = 0;
  ptps_it.go_to_part_no(0);
}

/*
 * shiftCCtoPosition(int pos)
 *  Given an absolute partition number 'pos', we shift the pair of
 *  random variables CC so that the right position (the 2nd
 *  C of CC) is now at position 'pos'. If the CE is shifted,
 *  we first shift CE back to position 0 before shifting CC.
 */
void
JunctionTree::shiftCCtoPosition(int pos)
{
  if (cur_ce_shift != 0) {
    assert ( cur_cc_shift == 0 );
    // we need to get CE back to zero.
    adjustFramesBy (cur_CE_rvs, -cur_ce_shift
		    *gm_template.S*fp.numFramesInC());
    cur_ce_shift = 0;
  }
  int delta = pos - cur_cc_shift;
  adjustFramesBy (cur_CC_rvs, delta
		  *gm_template.S*fp.numFramesInC());
  cur_cc_shift = pos;
}

/*
 * shiftOriginalVarstoPosition(vector<RV*> rvs, int pos, int &prevPos)
 *  Shift the random variables in 'rvs' in time (frames) by
 *  the difference between 'prevPos' and 'pos'. This is used  in 
 *  printSavedViterbiValues to adjust the rvs' frame numbers in
 *  the unpacking buffers.
 */

void
JunctionTree::shiftOriginalVarstoPosition(vector<RV*> rvs, int pos, int &prevPos)
{
  set<RV*> uprvs(rvs.begin(),rvs.end());
  int delta = (pos - prevPos);
  adjustFramesBy(uprvs, delta);
  prevPos = pos;
}

/*
 * shiftCCtoPosition(int pos)
 *  Given an absolute partition number 'pos', we shift the pair of
 *  random variables CE so that the right position (the E
 *  of CE) is now at position 'pos'. If the CC pair is shifted,
 *  we first shift CC back to position 0 before shifting CE.
 */
void
JunctionTree::shiftCEtoPosition(int pos)
{
  if (cur_cc_shift != 0) {
    assert ( cur_ce_shift == 0 );
    // we need to get CC back to zero.
    adjustFramesBy (cur_CC_rvs, -cur_cc_shift
		    *gm_template.S*fp.numFramesInC());
    cur_cc_shift = 0;
  }
  int delta = pos - cur_ce_shift;
  adjustFramesBy (cur_CE_rvs, delta
		  *gm_template.S*fp.numFramesInC());
  cur_ce_shift = pos;
}



/*
 * setCurrentInferenceShiftTo(int pos)
 *
 *  This activates partitions 'pos-1' and 'pos' so the
 *  right of two successive partitions is active at position 'pos'.
 *  This routine could be called:
 *
 *      activateTwoAdjacentPartitionsWithRightPartitionAtPos(pos)
 *      alignRightSideOfPartitionPairToPos(pos)
 * 
 *  Given an absolute partition number 'pos', we shift the pair of
 *  random variables (either CC or CE) so the *right* position (the 2nd
 *  C of CC or the E of CE) is now at position 'pos'. This means that
 *  messages entirely within position 'pos' and entirely within
 *  position 'pos-1' and between positions 'pos-1' and 'pos' will be
 *  correct, but no other messages are guaranteed to be correct.
 *
 * preconditions:
 *    the class member iterator 'inference_it' must be initilzed for the
 *    current and appropriate segment length.
 *
 * side effects
 *   - modifies the random variables corresponding to a partition pair
 *   - modifies the inference_it iterator.
 *
 */
void
JunctionTree::setCurrentInferenceShiftTo(int pos)
{

  if (inference_it.at_entry(pos)) {
    // printf("== Already at shift %d\n",pos);
    return;
  }

  //   printf("========================================\n");
  //   printf("== Setting current inference shift to %d\n",pos);
  //   printf("========================================\n");

  inference_it.go_to_part_no(pos);
  if (inference_it.num_c_partitions() <= 2) {
    // then do nothing, since nothing is never shifted away from zero.
  } else {
    // need to do some work. We need to get the right
    // of the appropriate pair (either the second C of CC in PCCE,
    // or the E of CE in PCCE) to be at position pos.
    if (inference_it.at_p()) {
      // P has an intersection with the first C in PCCE, so we need
      // to get that first C into the right position. We do that
      // by getting both CC's into the right position.
      shiftCCtoPosition(0);
    } else if (inference_it.at_e()) {
      // get the CE into the right position
      shiftCEtoPosition(pos - 3);
    } else {
      // Get CC into the right position, where 'pos' corresponds
      // to the second C in PCCE.
      if (pos <= 2)
	shiftCCtoPosition(0);
      else 
	shiftCCtoPosition(pos-2);
    }
  }
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//   Standard Collect/Distribute Evidence Support
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////



/*
 *
 * Preconditions: 
 *    The linear algorithm deScatterOutofRoot( .. ) 
 *    must have *just* been called and the random variables, and the iterator it has
 *    been set appropriately to the current partition.
 * Relies on:
 *   1) That 'it' is set to the current partition.
 *   2) That all random variables in the current structure (as designated by 'it') are set to 
 *      the current viterbi value.
 *
 */
void
JunctionTree::recordPartitionViterbiValue(ptps_iterator& it)
{
  PartitionStructures& ps = partitionStructureArray[it.ps_i()];
  unsigned partitionLength = ps.packer.packedLen();
  unsigned N_best = 1;
  unsigned num_to_write = N_best * partitionLength;
  if (partitionLength > 0)  {
    // if it is not greater than zero, then the partition has
    // no hidden discrete variables.

    if (binaryViterbiFile) {
      // Do a binary dump of the Viterbi values to a file
      // in chronological order for later printing by a separate
      // program. This is O(1) memory and O(T) disk space.

      if (it.at_p()) {
	if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset, SEEK_SET) == (off_t) -1) {
	  char *err = strerror(errno);
	  error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
	}	
	ps.packer.pack(ps.hrvValuePtrs.ptr,P_partition_values.ptr);
#if 0
printf("P     %03u %08x", (unsigned) ftell(binaryViterbiFile), P_partition_values.ptr[0]);
for (unsigned i=1; i < ps.packer.packedLen(); i+=1)
  printf(" %08x", P_partition_values.ptr[i]);
printf("\n");
#endif
	if (fwrite(P_partition_values.ptr, sizeof(unsigned), num_to_write, binaryViterbiFile) 
	    != num_to_write) 
	{
	  char *err = strerror(errno);
	  error("ERROR: write failed on '%s': %s\n", binaryViterbiFilename, err);
	}
      } else if (it.at_e()) {
	off_t offset = (off_t)  // P size + (T-2) * C size
	  (   (   N_best * partitionStructureArray[0].packer.packedLen()       
	        + N_best * partitionStructureArray[1].packer.packedLen() * it.num_c_partitions()  
	      ) * sizeof(unsigned)   );                         
	if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset + offset, SEEK_SET) == (off_t) -1) {
	  char *err = strerror(errno);
	  error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
	}	
	ps.packer.pack(ps.hrvValuePtrs.ptr,E_partition_values.ptr);
#if 0
printf("E %03llu %03u %08x", offset, (unsigned)ftell(binaryViterbiFile), E_partition_values.ptr[0]);
for (unsigned i=1; i < ps.packer.packedLen(); i+=1)
  printf(" %08x", E_partition_values.ptr[i]);
printf("\n");
#endif
	if (fwrite(E_partition_values.ptr, sizeof(unsigned), num_to_write, binaryViterbiFile) 
	    != num_to_write) 
	{
	  char *err = strerror(errno);
	  error("ERROR: write failed on '%s': %s\n", binaryViterbiFilename, err);
	}
	nextViterbiOffset = gmtk_ftell(binaryViterbiFile); // remember where to start next segment
	if (nextViterbiOffset == (off_t)-1) {
	  char *err = strerror(errno);
	  error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
	}
      } else { // at a C partition
	off_t offset = (off_t)  // P size + (t-1) * C size
          (   (   N_best * partitionStructureArray[0].packer.packedLen()       
		+ N_best * partitionLength * ( it.pt_i() - 1 )
              ) * sizeof(unsigned)   );               
	if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset + offset, SEEK_SET) == (off_t) -1) {
	  char *err = strerror(errno);
	  error("seek failed on '%s': %s\n", binaryViterbiFilename, err);
	}	
	ps.packer.pack(ps.hrvValuePtrs.ptr, C_partition_values.ptr);  // stomp on old values to be O(1) memory
#if 0
printf("C %03llu %03u %08x", offset, (unsigned) ftell(binaryViterbiFile), C_partition_values.ptr[0]);
for (unsigned i=1; i < ps.packer.packedLen(); i+=1)
  printf(" %08x", C_partition_values.ptr[i]);
printf("\n");
#endif
	if (fwrite(C_partition_values.ptr, sizeof(unsigned), num_to_write, binaryViterbiFile) 
	    != num_to_write) 
	{
	  char *err = strerror(errno);
	  error("write failed on '%s': %s\n", binaryViterbiFilename, err);
	}
      }      
    } else {
      // Store the Viterbi values in ?_partition_values in O(T)
      // memory for printing within this program.

      if (it.at_p()) {
	ps.packer.pack(ps.hrvValuePtrs.ptr,P_partition_values.ptr);
      } else if (it.at_e()) {
	ps.packer.pack(ps.hrvValuePtrs.ptr,E_partition_values.ptr);
      } else {
	ps.packer.pack(ps.hrvValuePtrs.ptr,
		       C_partition_values.ptr
		       + 
		       (it.pt_i()-1)*partitionLength);
      }
    }
  }
}



/*
 *
 * This routine saves the viterbi values computed by the most recent
 * linear inference run (assuming its data structures are still valid)
 * to stdout.
 *
 * Preconditions: 
 *
 *    Assumes that distributeEvidence has just been run and all data
 *    structures (such as the compressed viterbi value array) are set
 *    up appropriately. 
 *
 *    Assumes that inference_it is currently set for the current
 *    segment.
 *  
 *    Assumes that the CC and CE partition pair random variables
 *    have been properly set up.
 * 
 *
 */
void
JunctionTree::printSavedPartitionViterbiValues(FILE* f,
					       bool printObserved,
					       regex_t* preg,
					       char* partRangeFilter)
{
  fprintf(f,"Printing random variables from (P,C,E)=(%d,%d,%d) partitions\n",
	  P_partition_values.size(),
	  C_partition_values.size(),
	  E_partition_values.size());

  Range* partRange = NULL;
  if (partRangeFilter != NULL) {
    partRange = new Range(partRangeFilter,0,inference_it.pt_len());
    if (partRange->length() == 0) { 
      warning("WARNING: Part range filter must specify a valid non-zero "
	      "length range within [0:%d]. Range given is %s\n",
	      inference_it.pt_len(), partRangeFilter);
      delete partRange;
      partRange = NULL;
    }
  }
  if (partRange == NULL)
    partRange = new Range("all",0,inference_it.pt_len());

  Range::iterator* partRange_it = new Range::iterator(partRange->begin());

#if 0
  // unused
  int previous_C = -1;
#endif
  while (!partRange_it->at_end()) {

    unsigned part = (*partRange_it);
    setCurrentInferenceShiftTo(part);

    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

    if (inference_it.at_p()) {
      // print P partition
      if (ps.packer.packedLen() > 0) 
	ps.packer.unpack(P_partition_values.ptr,ps.hrvValuePtrs.ptr);
      if (printObserved && ps.allrvs.size() > 0) {
	fprintf(f,"Ptn-%d P': ",part);
	printRVSetAndValues(f,ps.allrvs,true,preg);
      } else if (ps.packer.packedLen() > 0) {
	fprintf(f,"Ptn-%d P': ",part);
	printRVSetAndValues(f,ps.hidRVVector,true,preg);
      }
    } else if (inference_it.at_e()) {
      // print E partition
      if (ps.packer.packedLen() > 0) 
	ps.packer.unpack(E_partition_values.ptr,ps.hrvValuePtrs.ptr);
      if (printObserved && ps.allrvs.size() > 0) {
	fprintf(f,"Ptn-%d E': ",part);
	printRVSetAndValues(f,ps.allrvs,true,preg);
      } else if (ps.packer.packedLen() > 0) {
	fprintf(f,"Ptn-%d E': ",part);
	printRVSetAndValues(f,ps.hidRVVector,true,preg);
      }
    } else {
      assert ( inference_it.at_c() );      
      // print C partition
#if 0
      if ((previous_C == -1)
	  ||
	  ps.packer.compare(C_partition_values.ptr
			    + 
			    (inference_it.pt_i()-1)*ps.packer.packedLen(),
			    C_partition_values.ptr
			    + 
			    (previous_C-1)*ps.packer.packedLen())) 
#endif
	{
	  if (ps.packer.packedLen() > 0)
	    ps.packer.unpack(C_partition_values.ptr
			     + 
			     (inference_it.pt_i()-1)*ps.packer.packedLen(),
			     ps.hrvValuePtrs.ptr);
	  if (printObserved && ps.allrvs.size() > 0) {
	    fprintf(f,"Ptn-%d C': ",part);
	    printRVSetAndValues(f,ps.allrvs,true,preg);
	  } else if (ps.packer.packedLen() > 0) {
	    fprintf(f,"Ptn-%d C': ",part);
	    printRVSetAndValues(f,ps.hidRVVector,true,preg);
	  }
	}
      // previous_C = inference_it.pt_i();
    }
    (*partRange_it)++;
  }
  delete partRange;
}



/*
 * This version of the above reads the Viterbi values from a file.
 *
 */
void
JunctionTree::printSavedPartitionViterbiValues(unsigned numFrames,
					       FILE* vitFile,
					       FILE* f,
					       bool printObserved,
					       regex_t* preg,
					       char* partRangeFilter)
{
  unsigned totalNumberPartitions;
  (void) unroll(numFrames,ZeroTable,&totalNumberPartitions);

  new (&inference_it) ptps_iterator(*this,totalNumberPartitions);
  init_CC_CE_rvs(inference_it);

  Range* partRange = NULL;
  if (partRangeFilter != NULL) {
    partRange = new Range(partRangeFilter,0,inference_it.pt_len());
    if (partRange->length() == 0) { 
      warning("WARNING: Part range filter must specify a valid non-zero length range "
	      "within [0:%d]. Range given is %s\n",inference_it.pt_len(),partRangeFilter);
      delete partRange;
      partRange = NULL;
    }
  }
  if (partRange == NULL)
    partRange = new Range("all",0,inference_it.pt_len());

  Range::iterator* partRange_it = new Range::iterator(partRange->begin());

  while (!partRange_it->at_end()) {
    unsigned part = (*partRange_it);
    setCurrentInferenceShiftTo(part);
    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

#if 0
    // obsoleted by N-best support
    if (inference_it.at_p()) {
      if (fread(P_partition_values.ptr, sizeof(unsigned), ps.packer.packedLen(), binaryViterbiFile) != ps.packer.packedLen())
	error("ERROR: faild to read viterbi values from '%s'\n", binaryViterbiFilename);
    } else if (inference_it.at_e()) {
      if (fread(E_partition_values.ptr, sizeof(unsigned), ps.packer.packedLen(), binaryViterbiFile) != ps.packer.packedLen())
	error("ERROR: faild to read viterbi values from '%s'\n", binaryViterbiFilename);
    } else {
      if (fread(C_partition_values.ptr, sizeof(unsigned), ps.packer.packedLen(), binaryViterbiFile) != ps.packer.packedLen())
	error("ERROR: faild to read viterbi values from '%s'\n", binaryViterbiFilename);
    }
#endif

    unsigned N_best = 1;
    unsigned num_to_read = N_best * ps.packer.packedLen();
    if (inference_it.at_p()) {
      if (fread(P_partition_values.ptr, sizeof(unsigned), num_to_read, binaryViterbiFile) != num_to_read) {
	char *err = strerror(errno);
	error("Error reading '%s': %s\n", binaryViterbiFilename, err);
      }
    } else if (inference_it.at_e()) {
      if (fread(E_partition_values.ptr, sizeof(unsigned), num_to_read, binaryViterbiFile) != num_to_read) {
	char *err = strerror(errno);
	error("Error reading '%s': %s\n", binaryViterbiFilename, err);
      }
    } else {
      if (fread(C_partition_values.ptr, sizeof(unsigned), num_to_read, binaryViterbiFile) != num_to_read) {
	char *err = strerror(errno);
	error("Error reading '%s': %s\n", binaryViterbiFilename, err);
      }
    }

    if (inference_it.at_p()) {
      // print P partition
      if (ps.packer.packedLen() > 0) 
	ps.packer.unpack(P_partition_values.ptr,ps.hrvValuePtrs.ptr);
	fprintf(f,"Ptn-%d P': ",part);
      if (printObserved && ps.allrvs.size() > 0) {
	printRVSetAndValues(f,ps.allrvs,true,preg);
      } else if (ps.packer.packedLen() > 0) {
	printRVSetAndValues(f,ps.hidRVVector,true,preg);
      }
    } else if (inference_it.at_e()) {
      // print E partition
      if (ps.packer.packedLen() > 0) 
	ps.packer.unpack(E_partition_values.ptr,ps.hrvValuePtrs.ptr);
	fprintf(f,"Ptn-%d E': ",part);
      if (printObserved && ps.allrvs.size() > 0) {
	printRVSetAndValues(f,ps.allrvs,true,preg);
      } else if (ps.packer.packedLen() > 0) {
	printRVSetAndValues(f,ps.hidRVVector,true,preg);
      }
    } else {
      assert ( inference_it.at_c() );      
      // print C partition
	{
	  if (ps.packer.packedLen() > 0)
	    ps.packer.unpack(C_partition_values.ptr, ps.hrvValuePtrs.ptr);
	    fprintf(f,"Ptn-%d C': ",part);
	  if (printObserved && ps.allrvs.size() > 0) {
	    printRVSetAndValues(f,ps.allrvs,true,preg);
	  } else if (ps.packer.packedLen() > 0) {
	    printRVSetAndValues(f,ps.hidRVVector,true,preg);
	  }
	}
      // previous_C = inference_it.pt_i();
    }
    (*partRange_it)++;
  }
  delete partRange;
  //clearAfterUnroll();
}





/*
 *
 * This routine saves the viterbi values computed by the most recent
 * linear inference run (assuming its data structures are still valid)
 * to stdout. This routine has a bunch of other options that are useful
 * for user-level printing. For printing by partition, see the partition
 * version of this routine.
 *
 * Preconditions: 
 *
 *    Assumes that distributeEvidence has just been run and all data
 *    structures (such as the compressed viterbi value array) are set
 *    up appropriately. 
 *
 *    Assumes that inference_it is currently set for the current
 *    segment.
 *  
 *    Assumes that the CC and CE partition pair random variables
 *    have been properly set up.
 * 
 *
 */
#if 0
void 
JunctionTree::printSavedViterbiValues(FILE* f,
				      bool printObserved,
				      regex_t *preg,
				      bool reverseOrder,
				      unsigned maxTriggerVars,
				      const char **triggerVars,
				      const char **triggerValSets)

{
  fprintf(f,"Printing random variable Viterbi values from %d frames%s.\n",curNumFrames,
	  (reverseOrder?" in reverse order":""));

  // approach: go through each frame, and search the corresponding partition
  // for the matching random variables 

  // The difficulty here is that a partition (i.e., a C') can have as
  // many diffrent frames as there are numbers of random variables. So
  // the number of frames *could* theoretically grow without bound.
  // Here, we only need to worry about the variables that are actually
  // being selected to be printed via preg and triggers.

  error("ERROR: function not implemented yet\n");

}
#endif


/* FTOC maps frame # f to chunk # given a prologue of P frames and chunks of C frames */
#define FTOC(P,C,f) \
  ( ((f) - (P)) / (C) )

/*
 * Create the printing (P,C,E _rvs) and unpacking (Pprime,Cprime,Eprime _rvs) buffers for printSavedViterbiValues
 */

void JunctionTree::createUnpackingMap(
  vector<RV*> &unrolled_rvs, map<RVInfo::rvParent, unsigned> &unrolled_map,
  vector<RV*> &P_rvs, vector<RV*> &hidP_rvs, vector<RV*> &Pprime_rvs, vector<RV *> &hidPprime_rvs,
  vector<vector<RV*> > &C_rvs, vector<vector<RV*> > &hidC_rvs, vector<vector<RV*> > &Cprime_rvs, vector<vector<RV*> > &hidCprime_rvs,
  vector<RV*> &E_rvs, vector<RV *> &hidE_rvs, vector<vector<RV*> > &Eprime_rvs, vector<vector<RV*> > &hidEprime_rvs,
  sArray<DiscRVType*> &PprimeValuePtrs, 
  vector<sArray<DiscRVType *> > &CprimeValuePtrs, 
  vector<sArray<DiscRVType *> > &EprimeValuePtrs) 
{
  unsigned M = gm_template.M;
  unsigned S = gm_template.S;
  infoMsg(IM::Printing, IM::Info,"creating C' -> C map, (M,S) = (%u,%u)    frames in (P,C,E) = (%d,%d,%d)\n", 
	  M, S, fp.numFramesInP(), fp.numFramesInC(), fp.numFramesInE());

  // the # of modified C' partitions needed  
  unsigned nCprimes = 1 + (M / S) + ( (M % S) ? 1 : 0 );

  // the # of original C partitions needed
  unsigned nCs = nCprimes * S;
  infoMsg(IM::Printing, IM::Moderate, "# C' = %u   # C = %u\n", nCprimes, nCs);

  // unroll to create RV instances shared by the printing & unpacking sets
  fp.unroll(nCs-1, unrolled_rvs, unrolled_map);

  unsigned NP = fp.numFramesInP();
  unsigned NC = fp.numFramesInC();

  P_rvs.clear(); hidP_rvs.clear();
  Pprime_rvs.clear(); hidPprime_rvs.clear();

  C_rvs.resize(nCs); hidC_rvs.resize(nCs);
  vector<set<RV*> > Calready(nCs);
  for (unsigned i=0; i < nCs; i+=1) {
    C_rvs[i].clear();
    hidC_rvs[i].clear();
    Calready[i].clear();
  }
  Cprime_rvs.resize(nCprimes); hidCprime_rvs.resize(nCprimes); 
  for (unsigned i=0; i < nCprimes; i+=1) {
    Cprime_rvs[i].clear();
    hidCprime_rvs[i].clear();
  }
  E_rvs.clear(); hidE_rvs.clear();
  Eprime_rvs.resize(nCprimes); hidEprime_rvs.resize(nCprimes);
  for (unsigned i=0; i < nCprimes; i+=1) {
    Eprime_rvs[i].clear();
    hidEprime_rvs[i].clear();
  }
  CprimeValuePtrs.resize(nCprimes);
  EprimeValuePtrs.resize(nCprimes);

  // the unpacking could end at any of the C' sets, which
  // in turn means that we don't know which of the C sets
  // unpacking E' will complete. So, we setup an E' set for
  // each possible case, and use the sentence length to 
  // determine which to use at unpacking time

  // TODO: Since the mapping is now based on partitionStructureArray,
  //       the number of C's in the unrolled model is known at map
  //       creation time. Thus we should be able to construct only
  //       the E' mapping for only the relevant C'E' transition 

  infoMsg(IM::Printing, IM::Moderate, "\nP':\n");
  vector<RV*> P = partitionStructureArray[0].allrvs_vec;
  for (vector<RV*>::iterator it = P.begin(); it != P.end(); ++it) {
    RV *v = getRV(unrolled_rvs, unrolled_map, *it);
    if ( (*it)->frame() < NP ) { // v is in P and P'
      P_rvs.push_back(v);
      Pprime_rvs.push_back(v);
      if (v->hidden()) {
	hidP_rvs.push_back(v);
	hidPprime_rvs.push_back(v);
      }
      infoMsg(IM::Printing, IM::Moderate, "%s(%u) -> P\n", (*it)->name().c_str(), (*it)->frame());
    } else {                     // v is in P' but not P
      Pprime_rvs.push_back(v);
      unsigned t = FTOC(NP,NC,(*it)->frame()); // the unmodified C index this variable belongs in
      C_rvs[t].push_back(v);
      Calready[t].insert(v);
      if (v->hidden()) {
	hidPprime_rvs.push_back(v);
	hidC_rvs[t].push_back(v);
      }
      infoMsg(IM::Printing, IM::Moderate, "%s(%u) -> C(%u)\n", (*it)->name().c_str(), (*it)->frame(), t);
    }
  }
  PprimeValuePtrs.resize(hidPprime_rvs.size());
  for (unsigned i=0; i < hidPprime_rvs.size(); i+=1) {
    DiscRV *drv = (DiscRV *) hidPprime_rvs[i];
    PprimeValuePtrs[i] = &(drv->val);
  }

  // If the unrolling consists of just P'E' then partitionStructureArray does
  // not contain a C' to base the C' -> C mapping on. But since there are no
  // C's to unpack, we don't need to construct the C' -> C mapping. Ticket #127
  if (partitionStructureArray.size() > 2) {
    // Map the first C' variables to their corresponding C (for printing)
    // and C' (for unpacking) instances
    infoMsg(IM::Printing, IM::Moderate, "\nC':\n");
    vector<RV*> C = partitionStructureArray[1].allrvs_vec;
    for (vector<RV*>::iterator it = C.begin(); it != C.end(); ++it) {
      RV *v = getRV(unrolled_rvs, unrolled_map, *it);
      unsigned t = FTOC(NP,NC,(*it)->frame()); // the unmodified C index this variable belongs in
#if 1
      C_rvs[t].push_back(v);
      Calready[t].insert(v);
      Cprime_rvs[0].push_back(v);
      if (v->hidden()) {
	hidC_rvs[t].push_back(v);
	hidCprime_rvs[0].push_back(v);
      }
#endif
      infoMsg(IM::Printing, IM::Moderate, "C(%u) C'[0] : %s(%u)\n", FTOC(NP,NC,(*it)->frame()), (*it)->name().c_str(), (*it)->frame());
    }

    // Slide each C'[0] variable forward S Cs to compute the
    // remaining C's
    for (unsigned i=1; i < nCprimes; i+=1) {
      infoMsg(IM::Printing, IM::Moderate, "\n");
      for (vector<RV*>::iterator it = Cprime_rvs[i-1].begin(); 
	   it != Cprime_rvs[i-1].end();
	   ++it)
	{
	  RV *v = *it;
	  unsigned f = NP + ( v->frame() - NP + S * NC ) % (NC * nCs); // shift v over S original Cs (mod nCs)
	  RVInfo::rvParent target(v->name(), f);
	  RV *rv = getRV(unrolled_rvs, unrolled_map, target);
	  unsigned t = FTOC(NP,NC,f);

	  // rv may already be in (hid?)C_rvs[t], resulting in excess length when printing C.
	  // We want the ordering of a vector, but the single membership of a set...
	  // Does need to be added to (hid?)Cprime_rvs[i], as that is newly constructed...
	  if (Calready[t].find(rv) == Calready[t].end()) {
	    C_rvs[t].push_back(rv);
	    if (rv->hidden()) {
	      hidC_rvs[t].push_back(rv);
	    }
	    Calready[t].insert(rv);
	  }
	  Cprime_rvs[i].push_back(rv);
	  if (rv->hidden()) {
	    hidCprime_rvs[i].push_back(rv);
	  }
	  infoMsg(IM::Printing, IM::Moderate, "C(%u) C'[%u] : %s(%u)\n", t, i, target.first.c_str(), target.second);
	}
    }
    for (unsigned i=0; i < nCprimes; i+=1) {
      CprimeValuePtrs[i].resize(hidCprime_rvs[i].size());
      for (unsigned j=0; j < hidCprime_rvs[i].size(); j+=1) {
	DiscRV *drv = (DiscRV *) hidCprime_rvs[i][j];
	CprimeValuePtrs[i][j] = &(drv->val);
      }
    }
  }

  infoMsg(IM::Printing, IM::Moderate, "\nE':\n");
  int Eidx  = partitionStructureArray.size() - 1;
  int delta = (nCs - ( (partitionStructureArray.size() - 2) * S + M ) ) * NC;
  unsigned firstEframe = NP + nCs * NC;
  infoMsg(IM::Printing, IM::Moderate, "delta = %d    firstE = %u\n", delta, firstEframe);
  vector<RV*> E = partitionStructureArray[Eidx].allrvs_vec;
  for (vector<RV*>::iterator it = E.begin(); it != E.end(); ++it) {
    RV *v = *it;
    if (v->frame() + delta >= firstEframe) { // v is in E and E'
      RV *rv = getRV(unrolled_rvs, unrolled_map, v, delta);
      E_rvs.push_back(rv);
      if (rv->hidden()) hidE_rvs.push_back(rv);
      for (unsigned i=0; i < nCprimes; i+=1) {
	Eprime_rvs[i].push_back(rv);
	if (rv->hidden()) hidEprime_rvs[i].push_back(rv);
      }
      infoMsg(IM::Printing, IM::Moderate, "%s(%u) -> %s(%u) -> E\n", v->name().c_str(), v->frame()+delta, rv->name().c_str(), rv->frame());
    } else { // v is in E' but not E
      // The C'E' transition could occur after unpacking any of the nCprimes C' sets,
      // so Eprime_rvs[i] handles the C'[i]E' transition
      for (unsigned i=0; i < nCprimes; i+=1) {
	RVInfo::rvParent target(v->name(), NP + (v->frame() + (S * NC * (1 - (Eidx-1))) - NP + S * NC * i) % (NC * nCs));
	RV *rv_shifted = getRV(unrolled_rvs, unrolled_map, target);
	Eprime_rvs[i].push_back(rv_shifted);
	if (rv_shifted->hidden()) hidEprime_rvs[i].push_back(rv_shifted);
      }
      infoMsg(IM::Printing, IM::Moderate, "%s(%u) -> C(%u)\n", v->name().c_str(), v->frame()+delta, FTOC(NP,NC,v->frame()+delta));
    }
  }
  for (unsigned i=0; i < nCprimes; i+=1) {
    EprimeValuePtrs[i].resize(hidEprime_rvs[i].size());
    for (unsigned j=0; j < hidEprime_rvs[i].size(); j+=1) {
      DiscRV *drv = (DiscRV *) hidEprime_rvs[i][j];
      EprimeValuePtrs[i][j] = &(drv->val);
    }
  }

  // print mapping as debug info
  if (IM::messageGlb(IM::Printing, IM::Info)) {
    infoMsg(IM::Printing,IM::Info, "G' -> G mapping:\n\n");

    set<RV*> pset(P_rvs.begin(), P_rvs.end());
    for (vector<RV*>::iterator it = Pprime_rvs.begin(); 
	 it != Pprime_rvs.end();
	 ++it)
      {
	RV *v = *it;
	if (pset.find(v) == pset.end()) {
	  infoMsg(IM::Printing,IM::Info, "P'    C[%u]: %s(%u)\n", FTOC(NP,NC,v->frame()), v->name().c_str(), v->frame());
	} else {
	  infoMsg(IM::Printing,IM::Info, "P'    P   : %s(%u)\n", v->name().c_str(), v->frame());
	}
      }
    for (unsigned i=0; i < nCprimes; i+=1) {
      infoMsg(IM::Printing, IM::Info, "-----------------------------\n");
      for (vector<RV*>::iterator it = Cprime_rvs[i].begin();
	   it != Cprime_rvs[i].end();
	   ++it)
	{
	  RV *v = *it;
	  infoMsg(IM::Printing,IM::Info, "C'[%u] C[%u]: %s(%u)\n", i, FTOC(NP,NC,v->frame()),v->name().c_str(), v->frame());
	}
    }

    for (unsigned i=0; i < nCprimes; i+=1) {
      infoMsg(IM::Printing, IM::Info, "-----------------------------\n");
      for (vector<RV*>::iterator it = Eprime_rvs[i].begin(); 
	   it != Eprime_rvs[i].end();
	   ++it)
	{
	  RV *v = *it;
	  unsigned c = FTOC(NP,NC,v->frame());
	  if (c < nCs) {
	    infoMsg(IM::Printing,IM::Info, "E'[%u] C[%u]: %s(%u)\n", i, c, v->name().c_str(), v->frame());
	  } else {
	    infoMsg(IM::Printing,IM::Info, "E'[%u] E   : %s(%u)\n", i, v->name().c_str(), v->frame());
	  }
	}
    }
    
#if 0
    infoMsg(IM::Printing, IM::Info, "-----------------------------\n");
    for (vector<RV*>::iterator it = E_rvs.begin();
	 it != E_rvs.end();
	 ++it)
      {
	RV *v = *it;
	infoMsg(IM::Printing,IM::Info, "E'    E   : %s(%u)\n", v->name().c_str(), v->frame());
      }
#endif
  }
}




/*
 *
 * This routine saves the viterbi values computed by the most recent
 * linear inference run (assuming its data structures are still valid)
 * to stdout. Unlike printSavedPartitionViterbiValues, this method
 * prints the values ordered by the original P, C, and E, rather than
 * the modified P', C', and E'.
 *
 * Preconditions: 
 *
 *    Assumes that distributeEvidence has just been run and all data
 *    structures (such as the compressed viterbi value array) are set
 *    up appropriately. 
 *
 *    Assumes that inference_it is currently set for the current
 *    segment.
 *  
 *    Assumes that the CC and CE partition pair random variables
 *    have been properly set up.
 * 
 *
 */
void
JunctionTree::printSavedViterbiValues(FILE* f,
				      bool printObserved,
				      regex_t* preg)
{

  vector<RV*> unrolled_rvs;
  map<RVInfo::rvParent, unsigned> unrolled_map;

  vector<RV*> P_rvs;      // original P for printing
  vector<RV*> Pprime_rvs; // modified P' for unpacking
  vector<RV*> hidP_rvs;      // hidden subset of original P for printing
  vector<RV*> hidPprime_rvs; // hidden subset of modified P' for unpacking

  vector<vector<RV*> > C_rvs; // original Cs for printing
  vector<vector<RV*> > Cprime_rvs; // modified C's for unpacking
  vector<vector<RV*> > hidC_rvs; // hidden subset of original Cs for printing
  vector<vector<RV*> > hidCprime_rvs; // hidden subset of modified C's for unpacking

  vector<RV*> E_rvs; // ... printing
  vector<vector<RV*> > Eprime_rvs; // ... unpacking
  vector<RV*> hidE_rvs; // ... printing
  vector<vector<RV*> > hidEprime_rvs; // ... unpacking

  sArray<DiscRVType *>PprimeValuePtrs;
  vector<sArray<DiscRVType *> > CprimeValuePtrs;
  vector<sArray<DiscRVType *> > EprimeValuePtrs;

  createUnpackingMap(unrolled_rvs, unrolled_map, 
		     P_rvs, hidP_rvs, Pprime_rvs, hidPprime_rvs, 
		     C_rvs, hidC_rvs, Cprime_rvs, hidCprime_rvs,
		     E_rvs, hidE_rvs, Eprime_rvs, hidEprime_rvs,
		     PprimeValuePtrs, CprimeValuePtrs, EprimeValuePtrs);

  fprintf(f,"Printing random variables from (P',C',E')=(%d,%d,%d) modified partitions\n",
#if 0
	  (int)P_rvs.size(),
	  (int)C_rvs[0].size(),
	  (int)E_rvs.size());
#else
	  P_partition_values.size(),
	  C_partition_values.size(),
	  E_partition_values.size());
#endif

  Range* partRange = new Range("all",0,inference_it.pt_len());

  Range::iterator* partRange_it = new Range::iterator(partRange->begin());

  vector<int> Cpos(C_rvs.size());
  for (unsigned int i=0; i < Cpos.size(); i+=1) 
    Cpos[i] = fp.numFramesInP() + i * fp.numFramesInC();
  int Epos = fp.numFramesInP() + C_rvs.size() * fp.numFramesInC();


  unsigned primeIndex = 0;
  unsigned originalIndex = 0;
  unsigned Ccount = 1;
#if 0
// unused
  int previous_C = -1;
#endif

  while (!partRange_it->at_end()) {

    unsigned part = (*partRange_it);
    setCurrentInferenceShiftTo(part);

    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
#if 0
    infoMsg(IM::Printing,IM::Info, "ps hid rvs:");
    for (vector<RV*>::iterator it = ps.hidRVVector.begin(); it != ps.hidRVVector.end(); ++it) {
      infoMsg(IM::Printing,IM::Info, " %s(%d)", (*it)->name().c_str(), (*it)->frame());
    }
    infoMsg(IM::Printing,IM::Info, "\n");
#endif
    if (inference_it.at_p()) {
      // print P partition
      if (ps.packer.packedLen() > 0) 
	ps.packer.unpack(P_partition_values.ptr,PprimeValuePtrs.ptr);
      if (hidP_rvs.size() > 0  ||  (printObserved && P_rvs.size() > 0) ) { 
	fprintf(f,"Ptn-0 P: ");
	if (printObserved)
	  printRVSetAndValues(f,P_rvs,true,preg);
	else
	  printRVSetAndValues(f,hidP_rvs,true,preg);
      }
    } else if (inference_it.at_e()) {
      primeIndex = (primeIndex + Eprime_rvs.size() - 1) % Eprime_rvs.size(); // primeIndex -= 1 mod nCprimes
      if (ps.packer.packedLen() > 0) 
	ps.packer.unpack(E_partition_values.ptr,EprimeValuePtrs[primeIndex].ptr);
      // print completed C partitions
      int targetFrame = fp.numFramesInP() + (int)(part-1) * gm_template.S * fp.numFramesInC();
      for (unsigned i=0; i < gm_template.M; i+=1) { // unpacking E' completes the last M Cs
	shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
	if (hidC_rvs[originalIndex].size() > 0 || (printObserved && C_rvs[originalIndex].size() > 0) ) {
	  fprintf(f,"Ptn-%u C: ", Ccount);
	  if (printObserved) 
	    printRVSetAndValues(f,C_rvs[originalIndex],true,preg);
	  else
	    printRVSetAndValues(f,hidC_rvs[originalIndex],true,preg);
	}
	Ccount += 1;
	originalIndex = (originalIndex + 1) % C_rvs.size();
	targetFrame += fp.numFramesInC();
      } 
      // print E partition
      if ( (hidE_rvs.size() > 0)  || (printObserved && E_rvs.size() > 0) ) {
	shiftOriginalVarstoPosition(E_rvs, targetFrame, Epos);
	fprintf(f,"Ptn-%u E: ", Ccount);
	if (printObserved) 
	  printRVSetAndValues(f,E_rvs,true,preg);
	else
	  printRVSetAndValues(f,hidE_rvs,true,preg);
      }
    } else {
      assert ( inference_it.at_c() );
      // print C partition
      {
        if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(C_partition_values.ptr
			   + 
			   (inference_it.pt_i()-1)*ps.packer.packedLen(),
			   CprimeValuePtrs[primeIndex].ptr);
	int targetFrame = fp.numFramesInP() + (int)(part-1) * gm_template.S * fp.numFramesInC();
	for (unsigned i=0; i < gm_template.S; i+=1) { // unpacking a C' completes S Cs
	  shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
	  if ( (hidC_rvs[originalIndex].size() > 0)  || (printObserved && C_rvs[originalIndex].size() > 0) ) {
	    fprintf(f,"Ptn-%u C: ", Ccount);
	    if (printObserved) 
	      printRVSetAndValues(f,C_rvs[originalIndex],true,preg);
	    else
	      printRVSetAndValues(f,hidC_rvs[originalIndex],true,preg);
	  }
	  Ccount += 1;
	  originalIndex = (originalIndex + 1) % C_rvs.size();
	  targetFrame += fp.numFramesInC();
	}
	primeIndex = (primeIndex + 1) % Cprime_rvs.size();
      }
      // previous_C = inference_it.pt_i();
    }
    (*partRange_it)++;
  }

  delete partRange;

}


/*
 * This version of the above reads the Viterbi values from a file.
 *
 */
void
JunctionTree::printSavedViterbiValues(unsigned numFrames,
				      FILE* f,
				      FILE* binVitFile,
				      bool printObserved,
				      regex_t* preg)
{
  unsigned totalNumberPartitions;
  (void) unroll(numFrames,ZeroTable,&totalNumberPartitions);

  new (&inference_it) ptps_iterator(*this,totalNumberPartitions);
  init_CC_CE_rvs(inference_it);

  vector<RV*> unrolled_rvs;
  map<RVInfo::rvParent, unsigned> unrolled_map;

  vector<RV*> P_rvs;      // original P for printing
  vector<RV*> Pprime_rvs; // modified P' for unpacking
  vector<RV*> hidP_rvs;      // hidden subset of original P for printing
  vector<RV*> hidPprime_rvs; // hidden subset of modified P' for unpacking

  vector<vector<RV*> > C_rvs; // original Cs for printing
  vector<vector<RV*> > Cprime_rvs; // modified C's for unpacking
  vector<vector<RV*> > hidC_rvs; // hidden subset of original Cs for printing
  vector<vector<RV*> > hidCprime_rvs; // hidden subset of modified C's for unpacking

  vector<RV*> E_rvs; // ... printing
  vector<vector<RV*> > Eprime_rvs; // ... unpacking
  vector<RV*> hidE_rvs; // ... printing
  vector<vector<RV*> > hidEprime_rvs; // ... unpacking

  sArray<DiscRVType *>PprimeValuePtrs;
  vector<sArray<DiscRVType *> > CprimeValuePtrs;
  vector<sArray<DiscRVType *> > EprimeValuePtrs;

  createUnpackingMap(unrolled_rvs, unrolled_map, 
		     P_rvs, hidP_rvs, Pprime_rvs, hidPprime_rvs, 
		     C_rvs, hidC_rvs, Cprime_rvs, hidCprime_rvs,
		     E_rvs, hidE_rvs, Eprime_rvs, hidEprime_rvs,
		     PprimeValuePtrs, CprimeValuePtrs, EprimeValuePtrs);

  Range* partRange = new Range("all",0,inference_it.pt_len());

  Range::iterator* partRange_it = new Range::iterator(partRange->begin());

  vector<int> Cpos(C_rvs.size());
  for (unsigned int i=0; i < Cpos.size(); i+=1) 
    Cpos[i] = fp.numFramesInP() + i * fp.numFramesInC();
  int Epos = fp.numFramesInP() + C_rvs.size() * fp.numFramesInC();


  unsigned primeIndex = 0;
  unsigned originalIndex = 0;
  unsigned Ccount = 1;

  while (!partRange_it->at_end()) {

    unsigned part = (*partRange_it);
    setCurrentInferenceShiftTo(part);
    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];


    unsigned N_best = 1;
    unsigned num_to_read = N_best * ps.packer.packedLen();
    if (inference_it.at_p()) {
      if (fread(P_partition_values.ptr, sizeof(unsigned), num_to_read, binaryViterbiFile) != num_to_read) {
	char *err = strerror(errno);
	error("Error reading '%s': %s\n", binaryViterbiFilename, err);
      }
    } else if (inference_it.at_e()) {
      if (fread(E_partition_values.ptr, sizeof(unsigned), num_to_read, binaryViterbiFile) != num_to_read) {
	char *err = strerror(errno);
	error("Error reading '%s': %s\n", binaryViterbiFilename, err);
      }
    } else {
      if (fread(C_partition_values.ptr, sizeof(unsigned), num_to_read, binaryViterbiFile) != num_to_read) {
	char *err = strerror(errno);
	error("Error reading '%s': %s\n", binaryViterbiFilename, err);
      }
    }

    if (inference_it.at_p()) {
      // print P partition
      if (ps.packer.packedLen() > 0) 
	ps.packer.unpack(P_partition_values.ptr,PprimeValuePtrs.ptr);
      if (hidP_rvs.size() > 0  ||  (printObserved && P_rvs.size() > 0) ) { 
	fprintf(f,"Ptn-0 P: ");
	if (printObserved)
	  printRVSetAndValues(f,P_rvs,true,preg);
	else
	  printRVSetAndValues(f,hidP_rvs,true,preg);
      }
    } else if (inference_it.at_e()) {
      primeIndex = (primeIndex + Eprime_rvs.size() - 1) % Eprime_rvs.size(); // primeIndex -= 1 mod nCprimes
      if (ps.packer.packedLen() > 0) 
	ps.packer.unpack(E_partition_values.ptr,EprimeValuePtrs[primeIndex].ptr);
      // print completed C partitions
      int targetFrame = fp.numFramesInP() + (int)(part-1) * gm_template.S * fp.numFramesInC();
      for (unsigned i=0; i < gm_template.M; i+=1) { // unpacking E' completes the last M Cs
	shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
	if (hidC_rvs[originalIndex].size() > 0 || (printObserved && C_rvs[originalIndex].size() > 0) ) {
	  fprintf(f,"Ptn-%u C: ", Ccount);
	  if (printObserved) 
	    printRVSetAndValues(f,C_rvs[originalIndex],true,preg);
	  else
	    printRVSetAndValues(f,hidC_rvs[originalIndex],true,preg);
	}
	Ccount += 1;
	originalIndex = (originalIndex + 1) % C_rvs.size();
	targetFrame += fp.numFramesInC();
      } 
      // print E partition
      if ( (hidE_rvs.size() > 0)  || (printObserved && E_rvs.size() > 0) ) {
	shiftOriginalVarstoPosition(E_rvs, targetFrame, Epos);
	fprintf(f,"Ptn-%u E: ", Ccount);
	if (printObserved) 
	  printRVSetAndValues(f,E_rvs,true,preg);
	else
	  printRVSetAndValues(f,hidE_rvs,true,preg);
      }
    } else {
      assert ( inference_it.at_c() );
      // print C partition
      {
        if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(C_partition_values.ptr, CprimeValuePtrs[primeIndex].ptr);
	int targetFrame = fp.numFramesInP() + (int)(part-1) * gm_template.S * fp.numFramesInC();
	for (unsigned i=0; i < gm_template.S; i+=1) { // unpacking a C' completes S Cs
	  shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
	  if ( (hidC_rvs[originalIndex].size() > 0)  || (printObserved && C_rvs[originalIndex].size() > 0) ) {
	    fprintf(f,"Ptn-%u C: ", Ccount);
	    if (printObserved) 
	      printRVSetAndValues(f,C_rvs[originalIndex],true,preg);
	    else
	      printRVSetAndValues(f,hidC_rvs[originalIndex],true,preg);
	  }
	  Ccount += 1;
	  originalIndex = (originalIndex + 1) % C_rvs.size();
	  targetFrame += fp.numFramesInC();
	}
	primeIndex = (primeIndex + 1) % Cprime_rvs.size();
      }
      // previous_C = inference_it.pt_i();
    }
    (*partRange_it)++;
  }

  delete partRange;
  //clearAfterUnroll();
}


void
JunctionTree::readBinaryVitPartition(PartitionStructures& ps, unsigned part) {
  assert(part == inference_it.pt_i());
  unsigned N_best = 1;
  unsigned num_to_read = N_best * ps.packer.packedLen();
  off_t offset;
  if (inference_it.at_p()) {
    if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset, SEEK_SET) == (off_t) -1) {
      char *err = strerror(errno);
      error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
    }       
    if (fread(P_partition_values.ptr, sizeof(unsigned), num_to_read, binaryViterbiFile) != num_to_read) {
      char *err = strerror(errno);
      error("Error reading '%s': %s\n", binaryViterbiFilename, err);
    }
  } else if (inference_it.at_e()) {
    offset = (off_t)  // P' size + (T-2) * C' size
      (   (   N_best * partitionStructureArray[0].packer.packedLen()       
	    + N_best * partitionStructureArray[1].packer.packedLen() * inference_it.num_c_partitions()  
	  ) * sizeof(unsigned)   );                         
    if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset + offset, SEEK_SET) == (off_t) -1) {
      char *err = strerror(errno);
      error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
    }       
    if (fread(E_partition_values.ptr, sizeof(unsigned), num_to_read, binaryViterbiFile) != num_to_read) {
      char *err = strerror(errno);
      error("Error reading '%s': %s\n", binaryViterbiFilename, err);
    }
  } else {
    offset = (off_t)  // P' size + (t-1) * C' size
      (   (   N_best * partitionStructureArray[0].packer.packedLen()
            + N_best * partitionStructureArray[1].packer.packedLen() * ( part - 1 )
	  ) * sizeof(unsigned)   );               
    if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset + offset, SEEK_SET) == (off_t) -1) {
      char *err = strerror(errno);
      error("seek failed on '%s': %s\n", binaryViterbiFilename, err);
    }       
    if (fread(C_partition_values.ptr, sizeof(unsigned), num_to_read, binaryViterbiFile) != num_to_read) {
      char *err = strerror(errno);
      error("Error reading '%s': %s\n", binaryViterbiFilename, err);
    }
  } 
}


void
JunctionTree::printSavedViterbiValues(unsigned numFrames, FILE* f,
				      FILE *binVitFile,
				      bool printObserved,
				      regex_t* preg,
				      char *partRangeFilter)
{

  if (partRangeFilter == NULL) {
    if (binaryViterbiFile) 
      printSavedViterbiValues(numFrames, f, binaryViterbiFile, printObserved, preg);
    else
      printSavedViterbiValues(f, printObserved, preg);
    return;
  }

  if (binaryViterbiFile) {
    unsigned totalNumberPartitions;
    (void) unroll(numFrames,ZeroTable,&totalNumberPartitions);
    
    new (&inference_it) ptps_iterator(*this,totalNumberPartitions);
    init_CC_CE_rvs(inference_it);
  }

  vector<RV*> unrolled_rvs;
  map<RVInfo::rvParent, unsigned> unrolled_map;

  vector<RV*> P_rvs;      // original P for printing
  vector<RV*> Pprime_rvs; // modified P' for unpacking
  vector<RV*> hidP_rvs;      // hidden subset of original P for printing
  vector<RV*> hidPprime_rvs; // hidden subset of modified P' for unpacking

  vector<vector<RV*> > C_rvs; // original Cs for printing
  vector<vector<RV*> > Cprime_rvs; // modified C's for unpacking
  vector<vector<RV*> > hidC_rvs; // hidden subset of original Cs for printing
  vector<vector<RV*> > hidCprime_rvs; // hidden subset of modified C's for unpacking

  vector<RV*> E_rvs; // ... printing
  vector<vector<RV*> > Eprime_rvs; // ... unpacking
  vector<RV*> hidE_rvs; // ... printing
  vector<vector<RV*> > hidEprime_rvs; // ... unpacking

  sArray<DiscRVType *>PprimeValuePtrs;
  vector<sArray<DiscRVType *> > CprimeValuePtrs;
  vector<sArray<DiscRVType *> > EprimeValuePtrs;

  createUnpackingMap(unrolled_rvs, unrolled_map, 
		     P_rvs, hidP_rvs, Pprime_rvs, hidPprime_rvs, 
		     C_rvs, hidC_rvs, Cprime_rvs, hidCprime_rvs,
		     E_rvs, hidE_rvs, Eprime_rvs, hidEprime_rvs,
		     PprimeValuePtrs, CprimeValuePtrs, EprimeValuePtrs);

  fprintf(f,"Printing random variables from (P',C',E')=(%d,%d,%d) modified partitions\n",
	  P_partition_values.size(),
	  C_partition_values.size(),
	  E_partition_values.size());

  unsigned M = gm_template.M;
  unsigned S = gm_template.S;
  unsigned totalOriginalPartitions = 2 + inference_it.num_c_partitions() * S + M;

  infoMsg(IM::Printing,IM::High,"M = %u   S = %u   # orig parts = %u\n",
       M, S, totalOriginalPartitions);

  Range* partRange = NULL;
  partRange = new Range(partRangeFilter,0,totalOriginalPartitions);
  if (partRange->length() == 0) { 
    warning("WARNING: Part range filter must specify a valid non-zero "
	    "length range within [0:%d]. Range given is %s\n",
	    totalOriginalPartitions, partRangeFilter);
    delete partRange;
    partRange = NULL;
  }

  if (partRange == NULL)
    partRange = new Range("all",0,totalOriginalPartitions);

  Range::iterator* partRange_it = new Range::iterator(partRange->begin());

  vector<int> Cpos(C_rvs.size());
  for (unsigned int i=0; i < Cpos.size(); i+=1) 
    Cpos[i] = fp.numFramesInP() + i * fp.numFramesInC();
  int Epos = fp.numFramesInP() + C_rvs.size() * fp.numFramesInC();
  
 
  int primeIndex = 0;     // which of the Cprime_rvs or Eprime_rvs to unpack to
  int originalIndex = 0;  // which of the C_rvs to print from

  while (!partRange_it->at_end()) {

    unsigned part = (*partRange_it);

    /* Before we can print this original partition $C_j$ (j = part), we must 
       unpack the modified partitions 

       $$\left\{ C'_i \left| \, \max\left(-1,\left\lceil\frac{j-s-m+1}{s}\right\rceil\right) \leq i 
         \leq \left\lfloor \frac{j}{s} \right\rfloor \right. \right\}$$

       where $s$ and $m$ are the boundary algorithm parameters, $C'_{-1}=P'$, and $C'_{N_{C'}}=E'$. 
     */
    
    int numerator = ( (int)part - 1 - (int)S - (int)M + 1 );
    int firstPrimePart;
    if (numerator <= -(int)S) {
      firstPrimePart = -1; // max
    } else if (numerator <= 0){
      firstPrimePart = 0;  // ceil
    } else {
      firstPrimePart = numerator / (int)S;
      if (numerator % (int)S)
	firstPrimePart += 1;   // ceil
    }
    firstPrimePart += 1; // account for C_{-1} = P'

    unsigned lastPrimePart =  (part > 0) ?  1 + ((int)part-1) / (int)S : 0;
    if (lastPrimePart >= inference_it.pt_len())  // E original partition # may > # of modified partitions
      lastPrimePart = inference_it.pt_len() - 1;

    infoMsg(IM::Printing,IM::High,"original partition %u requires unpacking modified partitions %u to %u:\n", 
	    part, firstPrimePart, lastPrimePart);

    for (unsigned i = (unsigned) firstPrimePart; i <= lastPrimePart; i += 1) { // unpack C'_{i} set
      infoMsg(IM::Printing,IM::High,"unpack %u'\n", i); 
      setCurrentInferenceShiftTo(i);
      PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

      if (binaryViterbiFile) {
	readBinaryVitPartition(ps, i);
      }
      
      if (inference_it.at_p()) { // P'
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(P_partition_values.ptr,PprimeValuePtrs.ptr);
      } else if (inference_it.at_e()) { // E'
	// -1 to get the preceding C', -1 to account for C'_{-1} = P'
	primeIndex = ((int)i - 2) % (int)Eprime_rvs.size(); 
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(E_partition_values.ptr,EprimeValuePtrs[primeIndex].ptr);
      } else { // C'
	assert ( inference_it.at_c() );
	primeIndex = ((int)i - 1) % (int)Cprime_rvs.size();
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(C_partition_values.ptr  + 
			   ( binaryViterbiFile ? 0 : (inference_it.pt_i()-1)*ps.packer.packedLen() ),
			   CprimeValuePtrs[primeIndex].ptr);
      }
    }

    if (part == 0) { // print P partition
      if (hidP_rvs.size() > 0  ||  (printObserved && P_rvs.size() > 0) ) { 
	fprintf(f,"Ptn-0 P: ");
	if (printObserved)
	  printRVSetAndValues(f,P_rvs,true,preg);
	else
	  printRVSetAndValues(f,hidP_rvs,true,preg);
      }
    } else if (part == totalOriginalPartitions-1) { // print E partition
      if ( (hidE_rvs.size() > 0)  || (printObserved && E_rvs.size() > 0) ) {
	int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
	shiftOriginalVarstoPosition(E_rvs, targetFrame, Epos);
	fprintf(f,"Ptn-%u E: ", totalOriginalPartitions-1);
	if (printObserved) 
	  printRVSetAndValues(f,E_rvs,true,preg);
	else
	  printRVSetAndValues(f,hidE_rvs,true,preg);
      }
    } else {      // print C partition
      int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
      originalIndex = ((int)part - 1) % (int) C_rvs.size();
      shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
      if ( (hidC_rvs[originalIndex].size() > 0)  || (printObserved && C_rvs[originalIndex].size() > 0) ) {
	fprintf(f,"Ptn-%u C: ", part);
	if (printObserved) 
	  printRVSetAndValues(f,C_rvs[originalIndex],true,preg);
	else
	  printRVSetAndValues(f,hidC_rvs[originalIndex],true,preg);
      }
    }
    (*partRange_it)++;
  }
  
  delete partRange;
}


void
JunctionTree::printSavedViterbiFrames(unsigned numFrames, FILE* f,
				      FILE *binVitFile,
				      bool printObserved,
				      regex_t* preg,
				      char *frameRangeFilter)
{
  unsigned numUsableFrames;
  if (binaryViterbiFile) {
    unsigned totalNumberPartitions;
    numUsableFrames = unroll(numFrames,ZeroTable,&totalNumberPartitions);
    
    new (&inference_it) ptps_iterator(*this,totalNumberPartitions);
    init_CC_CE_rvs(inference_it);
  } else {
    numUsableFrames = this->numUsableFrames;
  }

  vector<RV*> unrolled_rvs;
  map<RVInfo::rvParent, unsigned> unrolled_map;

  vector<RV*> P_rvs;      // original P for printing
  vector<RV*> Pprime_rvs; // modified P' for unpacking
  vector<RV*> hidP_rvs;      // hidden subset of original P for printing
  vector<RV*> hidPprime_rvs; // hidden subset of modified P' for unpacking

  vector<vector<RV*> > C_rvs; // original Cs for printing
  vector<vector<RV*> > Cprime_rvs; // modified C's for unpacking
  vector<vector<RV*> > hidC_rvs; // hidden subset of original Cs for printing
  vector<vector<RV*> > hidCprime_rvs; // hidden subset of modified C's for unpacking

  vector<RV*> E_rvs; // ... printing
  vector<vector<RV*> > Eprime_rvs; // ... unpacking
  vector<RV*> hidE_rvs; // ... printing
  vector<vector<RV*> > hidEprime_rvs; // ... unpacking

  sArray<DiscRVType *>PprimeValuePtrs;
  vector<sArray<DiscRVType *> > CprimeValuePtrs;
  vector<sArray<DiscRVType *> > EprimeValuePtrs;

  createUnpackingMap(unrolled_rvs, unrolled_map, 
		     P_rvs, hidP_rvs, Pprime_rvs, hidPprime_rvs, 
		     C_rvs, hidC_rvs, Cprime_rvs, hidCprime_rvs,
		     E_rvs, hidE_rvs, Eprime_rvs, hidEprime_rvs,
		     PprimeValuePtrs, CprimeValuePtrs, EprimeValuePtrs);

  fprintf(f,"Printing random variables from (P',C',E')=(%d,%d,%d) modified partitions\n",
	  P_partition_values.size(),
	  C_partition_values.size(),
	  E_partition_values.size());


  unsigned NP = fp.numFramesInP();
  unsigned NC = fp.numFramesInC();

  unsigned M = gm_template.M;
  unsigned S = gm_template.S;
  unsigned totalOriginalPartitions = 2 + inference_it.num_c_partitions() * S + M;

  infoMsg(IM::Printing,IM::High,"NP = %u   NC = %u   M = %u   S = %u   # orig parts = %u\n",
       NP, NC, M, S, totalOriginalPartitions);

  Range* frameRange = NULL;
  frameRange = new Range(frameRangeFilter,0,numUsableFrames);
  if (frameRange->length() == 0) { 
    warning("WARNING: Frame range filter must specify a valid non-zero "
	    "length range within [0:%d]. Range given is %s\n",
	    numUsableFrames, frameRangeFilter);
    delete frameRange;
    frameRange = NULL;
  }

  if (frameRange == NULL)
    frameRange = new Range("all",0,numUsableFrames);

  Range::iterator* frameRange_it = new Range::iterator(frameRange->begin());

  vector<int> Cpos(C_rvs.size());
  for (unsigned int i=0; i < Cpos.size(); i+=1) 
    Cpos[i] = fp.numFramesInP() + i * fp.numFramesInC();
  int Epos = fp.numFramesInP() + C_rvs.size() * fp.numFramesInC();
  
 
  int primeIndex = 0;     // which of the Cprime_rvs or Eprime_rvs to unpack to
  int originalIndex = 0;  // which of the C_rvs to print from

  int minAvailableFrame = -1; // nothing unpacked yet
  int maxAvailableFrame = -1;

  while (!frameRange_it->at_end()) {
    unsigned ppp = (*frameRange_it);
    infoMsg(IM::Printing,IM::High,"frame %u ", ppp);

    // map frame to original partition
    unsigned part;
    if ( (*frameRange_it) < (int)NP ) {
      part = 0; // P
    } else if ( (*frameRange_it) >= (int)NP + ((int)totalOriginalPartitions-2)*(int)NC ) {
      part = totalOriginalPartitions - 1; // E
    } else {
      part = 1 + ((*frameRange_it) - NP) / NC; // C
    }
    infoMsg(IM::Printing,IM::High,"in original partition %u ", part);

    // already fully unpacked? if so, print it
    if (minAvailableFrame <= (*frameRange_it) && (*frameRange_it) <= maxAvailableFrame) {
      infoMsg(IM::Printing,IM::High,"is available to print:\n");
      if (part == 0) { // print P partition
	if (hidP_rvs.size() > 0  ||  (printObserved && P_rvs.size() > 0) ) { 
	  fprintf(f,"Ptn-0 P: ");
	  if (printObserved)
	    printRVSetAndValues(f,P_rvs,true,preg, frameRange);
	  else
	    printRVSetAndValues(f,hidP_rvs,true,preg, frameRange);
	}
      } else if (part == totalOriginalPartitions-1) { // print E partition
	if ( (hidE_rvs.size() > 0)  || (printObserved && E_rvs.size() > 0) ) {
	  int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
	  shiftOriginalVarstoPosition(E_rvs, targetFrame, Epos);
	  fprintf(f,"Ptn-%u E: ", totalOriginalPartitions-1);
	  if (printObserved) 
	    printRVSetAndValues(f,E_rvs,true,preg, frameRange);
	  else
	    printRVSetAndValues(f,hidE_rvs,true,preg, frameRange);
	}
      } else {      // print C partition
	int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
	originalIndex = ((int)part - 1) % (int) C_rvs.size();
	shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
	if ( (hidC_rvs[originalIndex].size() > 0)  || (printObserved && C_rvs[originalIndex].size() > 0) ) {
	  fprintf(f,"Ptn-%u C: ", part);
	  if (printObserved) 
	    printRVSetAndValues(f,C_rvs[originalIndex],true,preg, frameRange);
	  else
	    printRVSetAndValues(f,hidC_rvs[originalIndex],true,preg, frameRange);
	}
      }
      (*frameRange_it)++;  // move on to next frame
      continue;
    }

    /* Before we can print this original partition $C_j$ (j = part), we must 
       unpack the modified partitions 

       $$\left\{ C'_i \left| \, \max\left(-1,\left\lceil\frac{j-s-m+1}{s}\right\rceil\right) \leq i 
         \leq \left\lfloor \frac{j}{s} \right\rfloor \right. \right\}$$

       where $s$ and $m$ are the boundary algorithm parameters, $C'_{-1}=P'$, and $C'_{N_{C'}}=E'$. 
     */
    
    int numerator = ( (int)part - 1 - (int)S - (int)M + 1 );
    int firstPrimePart;
    if (numerator <= -(int)S) {
      firstPrimePart = -1; // max
    } else if (numerator <= 0){
      firstPrimePart = 0;  // ceil
    } else {
      firstPrimePart = numerator / (int)S;
      if (numerator % (int)S)
	firstPrimePart += 1;   // ceil
    }
    firstPrimePart += 1; // account for C_{-1} = P'

    unsigned lastPrimePart = (part > 0) ?  1 + ((int)part-1) / (int)S : 0;
    if (lastPrimePart >= inference_it.pt_len())  // E original partition # may > # of modified partitions
      lastPrimePart = inference_it.pt_len() - 1;

    infoMsg(IM::Printing,IM::High,"requires unpacking modified partitions %u to %u:\n  unpack:", 
	    firstPrimePart, lastPrimePart);

    for (unsigned i = (unsigned) firstPrimePart; i <= lastPrimePart; i += 1) { // unpack C'_{i} set
      infoMsg(IM::Printing,IM::High,"  %u'", i); 
      setCurrentInferenceShiftTo(i);
      PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

      if (binaryViterbiFile) { // load packed values from disk if not already in memory
	readBinaryVitPartition(ps, i);
      }
      
      // unpack
      if (inference_it.at_p()) { // P'
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(P_partition_values.ptr,PprimeValuePtrs.ptr);
      } else if (inference_it.at_e()) { // E'
	// -1 to get the preceding C', -1 to account for C'_{-1} = P'
	primeIndex = ((int)i - 2) % (int)Eprime_rvs.size(); 
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(E_partition_values.ptr,EprimeValuePtrs[primeIndex].ptr);
      } else { // C'
	assert ( inference_it.at_c() );
	primeIndex = ((int)i - 1) % (int)Cprime_rvs.size();
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(C_partition_values.ptr  + 
			   ( binaryViterbiFile ? 0 : (inference_it.pt_i()-1)*ps.packer.packedLen() ),
			   CprimeValuePtrs[primeIndex].ptr);
      }
    }

    // unpacking firstPrimePart ... lastPrimePart makes these frames available:
    if (firstPrimePart == 0) {
      minAvailableFrame = 0;
    } else {
      minAvailableFrame = NP + (  ( (firstPrimePart-1) * S + M ) * NC  );
    }

    if (lastPrimePart == inference_it.pt_len()-1) {
      maxAvailableFrame = numUsableFrames - 1;
    } else {
      maxAvailableFrame = NP + lastPrimePart * S * NC - 1;
    }
    infoMsg(IM::Printing,IM::High,"  available frames %u to %u\n", minAvailableFrame, maxAvailableFrame);
  }
  
  delete frameRange;
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setRootToMaxCliqueValue()
 *   
 *   This version of probEvidence() merely sets E's root clique to its
 *   current max value, and return the max value.
 *
 * See Also:
 *   probEvidence() above, the other (overloaded) version which is
 *   a const. mem version of collectEvidence().
 *
 * Preconditions:
 *   - collectEvidence() (or the other probEvidence()) must have been called
 *     for this to work.
 *   - The E partition must exist in the partitionStructureArray array
 *     and must have been called for this to work right.
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   The probability of the evidence
 *
 *-----------------------------------------------------------------------
 */
logpr
JunctionTree::setRootToMaxCliqueValue()
{
  ptps_iterator ptps_it(*this);
  ptps_it.set_to_last_entry();
  if (partitionStructureArray[ptps_it.ps_i()].maxCliquesSharedStructure.size() == 0)
    --ptps_it;

  return partitionTableArray[ptps_it.pt_i()].maxCliques[ptps_it.cur_ri()].
    maxProbability(partitionStructureArray[ptps_it.ps_i()].maxCliquesSharedStructure[ptps_it.cur_ri()],
		   true);

}




///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//   Forward Messages
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::ceGatherIntoRoot
 *   
 *   Collect Evidence Gather Into Root: This routine does a collect
 *   evidence pass for this partition, and gathers all messages into
 *   the root within the current partition. It does so using the
 *   message order given in the argument 'message_order', and gathers
 *   into the provided root clique.  
 *
 * See Also:
 *   Dual routine: JunctionTree::deScatterOutofRoot()
 *
 *
 * Preconditions:
 *   It is assumed that either:
 *     1) this is the left-most partition
 *  or 2) that the left interface clique within this partition has had a message
 *        sent to it from the left neighbor partition.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::ceGatherIntoRoot(PartitionStructures& ps,
			       PartitionTables& pt,
			       // index of root clique in the partition
			       const unsigned root,
			       // message order of the JT in this partition
			       vector< pair<unsigned,unsigned> >& message_order,
			       // the name of the partition (for debugging/status msgs)
			       const char*const part_type_name,
			       // number of the partition in unrolled graph 
			       // (for printing/debugging/status msgs only)
			       const unsigned part_num,
			       const bool clearWhenDone,
			       const bool alsoClearOrigins)
{
  // first check that this is not an empty partition.
  if (ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! partitionDebugRange.contains((int)part_num)) {
#if 0
    printf("ceGather [part %u]: lowering inference level to %d\n", 
	   part_num, IM::glbMsgLevel(IM::DefaultModule));
#endif
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }

  bool zeroClique = false;
  try {
    // Now, do partition messages.
    for (unsigned msgNo=0;msgNo < message_order.size(); msgNo ++) {
      const unsigned from = message_order[msgNo].first;
      const unsigned to = message_order[msgNo].second;
      infoMsg(IM::Inference, IM::Med+5,
	      "CE: gathering into %s,part[%d]: clique %d\n",
	      part_type_name,part_num,from);

      // this may now throw an exception on zero clique errors - RR
      pt.maxCliques[from].
	ceGatherFromIncommingSeparators(ps.maxCliquesSharedStructure[from],
					pt.separatorCliques,
					ps.separatorCliquesSharedStructure.ptr);
  
      infoMsg(IM::Inference, IM::Mod,
	      "CE: message %s,part[%d]: clique %d --> clique %d\n",
	      part_type_name,part_num,from,to);
      pt.maxCliques[from].
	ceSendToOutgoingSeparator(ps.maxCliquesSharedStructure[from],
				  pt.separatorCliques,
				  ps.separatorCliquesSharedStructure.ptr);

      // TODO: if we are just computing probE here, we should delete
      // memory in pt.maxCliques[from]. Also, if we're only doing probE,
      // we should not keep the cliques around at all, only the outgoing
      // separator.
      if (clearWhenDone) {
	pt.maxCliques[from].
	  clearCliqueAndIncommingSeparatorMemory(ps.maxCliquesSharedStructure[from],
						 pt.separatorCliques,
						 ps.separatorCliquesSharedStructure.ptr);

	if (alsoClearOrigins) {
	  // then clear out the origin memory used for inference.
	  ps.origin.clearCliqueAndIncommingSeparatorMemoryForClique(from); 
	}
      }
    }
  } catch (ZeroCliqueException &e) {
    zeroClique = true; // abort this partition & segment
  }
  if (!zeroClique) {
    // collect to partition's root clique
    infoMsg(IM::Inference, IM::Med+5,
	    "CE: gathering into partition root %s,part[%d]: clique %d\n",
	    part_type_name,part_num,root);
    try {
      pt.maxCliques[root].
	ceGatherFromIncommingSeparators(ps.maxCliquesSharedStructure[root],
					pt.separatorCliques,
					ps.separatorCliquesSharedStructure.ptr);
    } catch (ZeroCliqueException &e) {
      zeroClique = true; // abort this partition & segment
    }
    if (!zeroClique) {
      if (IM::messageGlb(IM::InferenceMemory, IM::Med+9)) {
	pt.reportMemoryUsageTo(ps,stdout);
      }
    }
  }
  if (! partitionDebugRange.contains((int)part_num)) {
#if 0
    printf("ceGather [part %u]: raising inference level to %d\n", 
	   part_num, inferenceDebugLevel);
#endif
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }
  if (zeroClique) {
    throw ZeroCliqueException(); // continue to abort segment
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::ceSendForwardsCrossPartitions
 *   
 *   Collect Evidence Send To Next Partition: This routine sends a
 *   message from the right interface clique of a left (or previous)
 *   partition to the left interface clique of a right (or next)
 *   partition in the partition series. It is assumed that the right
 *   interface clique has had all its incomming messages sent to it.
 *
 * See Also:
 *   Dual routine: JunctionTree::deSendBackwardsCrossPartitions()
 *
 *
 * Preconditions:
 *   It is assumed that:
 *     1) the right interface of the previous partition must have had
 *        all messages sent to it.
 * 
 * Postconditions:
 *     the left interface of the next partition is now set up.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::ceSendForwardsCrossPartitions(// previous partition
				    PartitionStructures& previous_ps,
				    PartitionTables& previous_pt,
				    // root clique of the previous partition (i.e., the
				    // right interface clique) 
				    const unsigned previous_part_root,
				    // name of previous partition (for debugging/status msgs)
				    const char*const previous_part_type_name,
				    // sequence number (in unrolling) of previous partition
				    // (for debugging/status msgs)
				    const unsigned previous_part_num,
				    // next partition
				    PartitionStructures& next_ps,
				    PartitionTables& next_pt,
				    // leaf clique of next partition (i.e., index number
				    // of the left interface clique of next partition)
				    const unsigned next_part_leaf,
				    // name (debugging/status msgs)
				    const char*const next_part_type_name,
				    // partitiiton number (debugging/status msgs)
				    const unsigned next_part_num)
{
  // check for empty partitions.
  if (previous_ps.maxCliquesSharedStructure.size() == 0 || next_ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! partitionDebugRange.contains((int)next_part_num)) {
#if 0
    printf("ceGather [part %u]: lowering inference level to %d\n", 
	   next_part_num, IM::glbMsgLevel(IM::DefaultModule));
#endif
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }


  infoMsg(IM::Inference, IM::Mod,"CE: message %s,part[%d],clique(%d) --> %s,part[%d],clique(%d)\n",
	  previous_part_type_name,
	  previous_part_num,
	  previous_part_root,
	  next_part_type_name,
	  next_part_num,
	  next_part_leaf);
  previous_pt.maxCliques[previous_part_root].
    ceSendToOutgoingSeparator(previous_ps.maxCliquesSharedStructure[previous_part_root],
			      next_pt.separatorCliques[next_ps.separatorCliquesSharedStructure.size()-1],
			      next_ps.separatorCliquesSharedStructure[next_ps.separatorCliquesSharedStructure.size()-1]);

  if (IM::messageGlb(IM::InferenceMemory, IM::Med+9)) {
    previous_pt.reportMemoryUsageTo(previous_ps,stdout);
  }

  if (! partitionDebugRange.contains((int)next_part_num)) {
#if 0
    printf("ceGather [part %u]: raising inference level to %d\n", 
	   next_part_num, inferenceDebugLevel);
#endif
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectEvidence()
 *   
 *   Collect Evidence: This routine performes a complete collect
 *   evidence pass for this series of partitions that have been
 *   unrolled a certain amount, given by the partitionStructureArray array. It
 *   sets the appropriate names, etc. of the partitions depending on
 *   if this is a left-interface or right-interface form of inference.
 *   the root within the current partition. 
 *  
 *   This routine demonstrates a simple version of linear space
 *   collect evidence inference. It keeps everything in memory at the
 *   same time when doing going forward, so it is suitable for use in
 *   a collect evidence/distribute evidence (forward/backward)
 *   framework.  A companion routine (aptly named
 *   'distributeEvidence') will, assuming collect evidence has been
 *   called, do the distribute evidence stage, thereby leaving all
 *   cliques in all partitions locally and therefore globally
 *   consistent, ready to be used say in EM training.
 *
 * See Also:
 *    1) contant memory (not dept. on time T) version of collect evidence
 *       JunctionTree::probEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *   The parititons must have been created and placed in the array partitionStructureArray.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::collectEvidence()
{
  // This routine handles all of:
  // 
  //    unrolled 0 times: (so there is a single P1, and E1)  
  //    unrolled 1 time: so there is a P1, C1, C3, E1
  //    unrolled 2 or more times: so there is a P1 C1 [C2 ...] C3, E1
  //    etc.

  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this);

  init_CC_CE_rvs(inference_it);

  // we skip the first Co's LI separator if there is no P1
  // partition, since otherwise we'll get zero probability.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.skipLISeparator();
  // gather into the root of the current  partition
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   partitionTableArray[inference_it.pt_i()],
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i());
  // if the LI separator was turned off, we need to turn it back on.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.useLISeparator();

  for (unsigned part=1;part < inference_it.pt_len(); part ++ ) {

    setCurrentInferenceShiftTo(part);

    // inference_it.printState(stdout);

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[inference_it.ps_prev_i()],
			  partitionTableArray[inference_it.pt_prev_i()],
			  inference_it.prev_ri(),
			  inference_it.prev_nm(),
			  inference_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[inference_it.ps_i()],
			  partitionTableArray[inference_it.pt_i()],
			  inference_it.cur_li(),
			  inference_it.cur_nm(),
			  inference_it.pt_i());


    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // it might be that E is the first partition as well, say if this is
    // a static graph, and in this case we need in this case to skip the
    // incomming separator, which doesn't exist.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();
    // next, gather into the root of the final E partition
    ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		     partitionTableArray[inference_it.pt_i()],
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();

    // if the LI separator was turned off, we need to turn it back on.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();

  }
  assert ( inference_it.at_e() );


#if 0
  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this);

  init_CC_CE_rvs(inference_it);
  setCurrentInferenceShiftTo(0);

  // we skip the first Co's LI separator if there is no P1
  // partition, since otherwise we'll get zero probability.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.skipLISeparator();
  // gather into the root of the current  partition
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   partitionTableArray[inference_it.pt_i()],
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i());
  // if the LI separator was turned off, we need to turn it back on.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.useLISeparator();

  for (unsigned part=1;part < inference_it.pt_len(); part ++ ) {
    setCurrentInferenceShiftTo(part);

    // inference_it.printState(stdout);

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[inference_it.ps_prev_i()],
			  partitionTableArray[inference_it.pt_prev_i()],
			  inference_it.prev_ri(),
			  inference_it.prev_nm(),
			  inference_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[inference_it.ps_i()],
			  partitionTableArray[inference_it.pt_i()],
			  inference_it.cur_li(),
			  inference_it.cur_nm(),
			  inference_it.pt_i());


    // it might be that E is the first partition as well, say if this is
    // a static graph, and in this case we need in this case to skip the
    // incomming separator, which doesn't exist.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();
    // next, gather into the root of the final E partition
    ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		     partitionTableArray[inference_it.pt_i()],
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();
  }
  assert ( inference_it.at_e() );

#endif

}



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//   Backwards Messages
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::deScatterOutofRoot
 *   
 *   Distribute Evidence Scatter Outof Root: This routine does a
 *   distribute evidence pass for this partition, and scatters all
 *   messages outof the root within the current partition. It does so
 *   using the message order given in the argument 'message_order',
 *   and scatters out of the provided root clique (which is the right
 *   interface clique of this partition). By "Scatter", I mean it
 *   sends messages from the root clique distributing everything
 *   ultimately to all leaf cliques in this partition.
 *
 * See Also:
 *   Dual routine: JunctionTree::ceGatherIntoRoot()
 *
 *
 * Preconditions:
 *   It is assumed that either:
 *     1) this is the ritht-most partition
 *  or 2) that the right interface clique within this partition has had a message
 *        sent to it from the right neighbor partition.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::deScatterOutofRoot(// the partition
				 PartitionStructures& ps,
				 PartitionTables& pt,
				 // root (right interface clique) of this partition
				 const unsigned root,
				 // message order
				 vector< pair<unsigned,unsigned> >& message_order,
				 // name (debugging/status msgs)
				 const char*const part_type_name,
				 // partition number (debugging/status msgs)
				 const unsigned part_num)
{
  if (ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! partitionDebugRange.contains((int)part_num)) {
#if 0
    printf("deScatter [part %u]: lowering inference level to %d\n", part, IM::glbMsgLevel(IM::DefaultModule));
#endif
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }

  infoMsg(IM::Inference, IM::Med+5,"DE: distributing out of partition root %s,part[%d]: clique %d\n",
	  part_type_name,part_num,root);
  pt.maxCliques[root].
    deScatterToOutgoingSeparators(ps.maxCliquesSharedStructure[root],
				  pt.separatorCliques,
				  ps.separatorCliquesSharedStructure.ptr);
  for (unsigned msgNoP1=message_order.size();msgNoP1 > 0; msgNoP1 --) {
    const unsigned to = message_order[msgNoP1-1].first;
    const unsigned from = message_order[msgNoP1-1].second;
    infoMsg(IM::Inference, IM::Mod,"DE: message %s,part[%d]: clique %d <-- clique %d\n",
	    part_type_name,part_num,to,from);
    pt.maxCliques[to].
      deReceiveFromIncommingSeparator(ps.maxCliquesSharedStructure[to],
				      pt.separatorCliques,
				      ps.separatorCliquesSharedStructure.ptr);

    infoMsg(IM::Inference, IM::Med+5,"DE: distributing out of %s,part[%d]: clique %d\n",
	    part_type_name,part_num,to);
    pt.maxCliques[to].
      deScatterToOutgoingSeparators(ps.maxCliquesSharedStructure[to],
				    pt.separatorCliques,
				    ps.separatorCliquesSharedStructure.ptr);
  }

  if (! partitionDebugRange.contains((int)part_num)) {
#if 0
    printf("deScatter [part %u]: raising inference level to %d\n", part, inferenceDebugLevel);
#endif
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }


}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::deSendBackwardsCrossPartitions()
 *   
 *   Distribute Evidence Receive To Previous Partition: This routine
 *   sends a message from the left interface clique of a right (or
 *   next) partition to the right interface clique of a left (or
 *   previous) partition in the partition series. It is assumed that
 *   the left interface clique has had all its messages sent to it.
 *
 * See Also:
 *   Dual routine: JunctionTree::ceSendForwardsCrossPartitions()
 *
 * Preconditions:
 *   It is assumed that:
 *     1) the left interface of the next partition must have had
 *        all messages sent to it.
 * 
 * Postconditions:
 *     the right interface of the previous partition is now set up.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::deSendBackwardsCrossPartitions(// previous partition
					     PartitionStructures& previous_ps,
					     PartitionTables& previous_pt,
					     // root (right interface clique) of previous partition 
					     const unsigned previous_part_root,
					     // name of prev part (debugging/status msgs)
					     const char*const previous_part_type_name,
					     // number of prev part (debugging/status msgs)
					     const unsigned previous_part_num,

					     // the next partition
					     PartitionStructures& next_ps,
					     PartitionTables& next_pt,
					     // leaf (left interface cliuqe) of next partition
					     const unsigned next_part_leaf,
					     // name of next part (debugging/status msgs)
					     const char*const next_part_type_name,
					     // number of next part (debugging/status msgs)
					     const unsigned next_part_num
					     )
{
  // check for empty partitions.
  if (previous_ps.maxCliquesSharedStructure.size() == 0 || next_ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! partitionDebugRange.contains((int)previous_part_num)) {
#if 0
    printf("deScatter [part %u]: lowering inference level to %d\n", 
	   previous_part_num, IM::glbMsgLevel(IM::DefaultModule));
#endif
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }

  infoMsg(IM::Inference, IM::Mod,"DE: message %s,part[%d],clique(%d) <-- %s,part[%d],clique(%d)\n",
	  previous_part_type_name,previous_part_num,previous_part_root,
	  next_part_type_name,next_part_num,next_part_leaf);
  previous_pt.maxCliques[previous_part_root].
    deReceiveFromIncommingSeparator(previous_ps.maxCliquesSharedStructure[previous_part_root],
				    next_pt.separatorCliques[next_ps.separatorCliquesSharedStructure.size()-1],
				    next_ps.separatorCliquesSharedStructure[next_ps.separatorCliquesSharedStructure.size()-1]);

  if (! partitionDebugRange.contains((int)previous_part_num)) {
#if 0
    printf("deScatter [part %u]: raising inference level to %d\n", 
	   previous_part_num, inferenceDebugLevel);
#endif
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::distributeEvidence()
 *   
 *   Distribute Evidence: This routine performs a complete distribute
 *   evidence pass for this series of partitions that have been
 *   unrolled a certain amount, given by the partitionStructureArray array. It
 *   sets the appropriate names, etc. of the partitions depending on
 *   if this is a left-interface or right-interface form of inference.
 *   the root within the current partition.
 *  
 *   This routine demonstrates a simple version of linear space
 *   distribute evidence inference. It uses data kept in memory when
 *   doing going backwards, and it is assumed that it will be used in
 *   tandem with a previous collect evidence call.  A companion
 *   routine ('collectEvidence') will, do collect evidence appropriately
 *   leaving all data structures set up for distributeEvidence() to be called.
 *
 *   After distributeEvidence() is called, all cliques in all
 *   partitions will be locally (& globally if it is a JT) consistant, and so will
 *   be ready for EM training, etc.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) contant memory (not dept. on time T) version of collect evidence
 *       JunctionTree::probEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *   - collectEvidence() must have been called right before this.
 *   - The parititons must have been created and placed in the array partitionStructureArray.
 *   - inference_it must be initialized to the current segment
 *   - the pair of partition RVS set up correct.
 * 
 *
 * Postconditions:
 *     All cliques are now locally consistant, and will be globally consistant
 *     if we have a junction tree. 
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::distributeEvidence()
{
  for (unsigned part= (inference_it.pt_len()-1) ; part > 0 ; part -- ) {

    setCurrentInferenceShiftTo(part);
    
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();    
    else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();

    deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		       partitionTableArray[inference_it.pt_i()],
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();
    else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();

    if (viterbiScore)
      recordPartitionViterbiValue(inference_it);

    // send backwads message to previous partition
    deSendBackwardsCrossPartitions(partitionStructureArray[inference_it.ps_prev_i()],
				   partitionTableArray[inference_it.pt_prev_i()],
				   inference_it.prev_ri(),
				   inference_it.prev_nm(),
				   inference_it.pt_prev_i(),
				   //
				   partitionStructureArray[inference_it.ps_i()],
				   partitionTableArray[inference_it.pt_i()],
				   inference_it.cur_li(),
				   inference_it.cur_nm(),
				   inference_it.pt_i());
  }

  setCurrentInferenceShiftTo(0);
  // do the final scatter out of root of initial P partition.

  deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		     partitionTableArray[inference_it.pt_i()],
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());

  if (viterbiScore)
    recordPartitionViterbiValue(inference_it);



#if 0

  ptps_iterator ptps_it(*this);
  ptps_it.set_to_last_entry();

  // the set of rvs that need to be adjusted, if any.
  set <RV*> two_part_rvs;
  if (ptps_it.num_c_partitions() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // partitions (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = partitionStructureArray[3].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,two_part_rvs);

    // we start by adjusting C E up to appropriate values.
    adjustFramesBy (two_part_rvs, (ptps_it.num_c_partitions()-2)
		    *gm_template.S*fp.numFramesInC());    

  } 


  while (!ptps_it.at_first_entry()) {

    // ptps_it.printState(stdout);

    if (ptps_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();    
    else if (!ptps_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();

    deScatterOutofRoot(partitionStructureArray[ptps_it.ps_i()],
		       partitionTableArray[ptps_it.pt_i()],
		       ptps_it.cur_ri(),
		       ptps_it.cur_message_order(),
		       ptps_it.cur_nm(),
		       ptps_it.pt_i());
    if (ptps_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();
    else if (!ptps_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();
    if (viterbiScore)
      recordPartitionViterbiValue(ptps_it);

    // send backwads message to previous partition
    deSendBackwardsCrossPartitions(partitionStructureArray[ptps_it.ps_prev_i()],
				   partitionTableArray[ptps_it.pt_prev_i()],
				   ptps_it.prev_ri(),
				   ptps_it.prev_nm(),
				   ptps_it.pt_prev_i(),
				   //
				   partitionStructureArray[ptps_it.ps_i()],
				   partitionTableArray[ptps_it.pt_i()],
				   ptps_it.cur_li(),
				   ptps_it.cur_nm(),
				   ptps_it.pt_i());


    if (ptps_it.num_c_partitions() > 2) {
      // then we need to do some RV shifting.
      if (ptps_it.at_e()) {
	// we need to shift the two C',E' parts back to normal, and then shift
	// the (C,C) part forward appropriately

	// First, shift C',E' rvs back to normal frame
	adjustFramesBy (two_part_rvs, -(ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

	set <RV*> tmp1,tmp2;
	tmp1 = partitionStructureArray[1].returnRVsAndTheirObservedParentsAsSet();
	tmp2 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
	// We're done now with two_part_rvs, so we make a new one. The
	// next routine clears out the original set values.
	unionRVs(tmp1,tmp2,two_part_rvs);
	// now adjust C',C' variables to match.
	adjustFramesBy (two_part_rvs, (ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

      } else if (ptps_it.pt_i() > 2) {
	// shift the two C partitions over by one more
	adjustFramesBy (two_part_rvs, -1*gm_template.S*fp.numFramesInC());
      }
    }

    --ptps_it;
  }

  // do the final scatter out of root of initial P partition.
  deScatterOutofRoot(partitionStructureArray[ptps_it.ps_i()],
		     partitionTableArray[ptps_it.pt_i()],
		     ptps_it.cur_ri(),
		     ptps_it.cur_message_order(),
		     ptps_it.cur_nm(),
		     ptps_it.pt_i());



  if (viterbiScore)
    recordPartitionViterbiValue(ptps_it);

#endif
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::emIncrement()
 *
 *    A version of emIncrement that works with the data structures set up by 
 *    collectEvidence() and distributeEvidence() above. Note, EM training
 *    can also be performed with the island algorithm.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) distributeEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *    collectEvidence() AND distributeEvidence() must have just been
 *    called setting up the data structures. All cliques must exist and
 *    must have their clique tables filled out and ready.
 *    Also, all parametr accumulators must be set up and ready to go.
 *
 * Postconditions:
 *   The accumulators are increment accordingly.
 *
 * Side Effects:
 *   This will update the accumulators of all trainable parameter objects.
 *
 * Results:
 *   None
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::emIncrement(const logpr probE,
			  const bool localCliqueNormalization,
			  const double emTrainingBeam)
{

  // Quite simply, just iterate through all partitions and call emIncrement
  // therein.

  // forward order
  // 
  // for (unsigned part=0;part < inference_it.pt_len(); part ++ ) {
  // 
  // reverse order increment (which is the same as island algorithm order,
  // and this can make slight numeric differences). To avoid
  // questions about this ("why does island EM produce small differences
  // than normal EM???") we use the reverse order here.
  for (unsigned part=inference_it.pt_len();part > 0 ; part --) {
    setCurrentInferenceShiftTo(part-1);
    infoMsg(IM::Training, IM::High-1,
	    "EM: accumulating stats for %s,part[%d]\n",
	    inference_it.cur_nm(),inference_it.pt_i());
    partitionTableArray[inference_it.pt_i()].
      emIncrement(partitionStructureArray[inference_it.ps_i()],
		  probE,
		  localCliqueNormalization,
		  emTrainingBeam);
  }



#if 0
  // Quite simply, just iterate through all partitions and call emIncrement
  // therein.
  ptps_iterator ptps_it(*this);
  ptps_it.set_to_first_entry();

  // the set of rvs that need to be adjusted, if any.
  set <RV*> two_part_rvs;
  if (ptps_it.num_c_partitions() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // partitions (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = partitionStructureArray[1].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,two_part_rvs);
  }

  while (1) {
    infoMsg(IM::Training, IM::High-1,
	    "EM: accumulating stats for %s,part[%d]\n",
	    ptps_it.cur_nm(),ptps_it.pt_i());

    if (ptps_it.num_c_partitions() > 2) {
      // then we need to do some RV shifting.
      if (ptps_it.at_e()) {
	// we need to shift the two C1',C2' parts back to normal, and then shift
	// the (C,E) part forward.

	// First, shift C1',C2' rvs back to normal frame
	adjustFramesBy (two_part_rvs, -(ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

	set <RV*> tmp1,tmp2;
	tmp1 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
	tmp2 = partitionStructureArray[3].returnRVsAndTheirObservedParentsAsSet();
	// We're done now with two_part_rvs, so we make a new one. The
	// next routine clears out the original set values.
	unionRVs(tmp1,tmp2,two_part_rvs);
	// now adjust C2',E' variables to match.
	adjustFramesBy (two_part_rvs, (ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

      } else if (ptps_it.pt_i() > 2) {
	// shift the two C partitions over by one more
	adjustFramesBy (two_part_rvs, 1*gm_template.S*fp.numFramesInC());
      }
    }

    partitionTableArray[ptps_it.pt_i()].
      emIncrement(partitionStructureArray[ptps_it.ps_i()],
		  probE,
		  localCliqueNormalization,
		  emTrainingBeam);
    if (ptps_it.at_last_entry())
      break;
    else
      ++ptps_it;
  }


#endif
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::probEvidence()
 *   
 *   This version of probEvidence() merely sums up the entries in E's
 *   root clique and return the result. This, if collect evidence has
 *   been called all the way to E's root clique, this function will
 *   return the probabilty of the evidence prob(E).
 *
 *   Note that if we are in viterbi mode, then rather than summing
 *   we return the score of the maximum value clique entry.
 *
 * See Also:
 *   probEvidence() above, the other (overloaded) version which is
 *   a const. mem version of collectEvidence().
 *
 * Preconditions:
 *   - collectEvidence() (or the other probEvidence()) must have been called
 *   - The E partition must exist in the partitionStructureArray array
 *     and must have been called for this to work right.
 *   - both the partitionStructureArray and the partitionTableArray must
 *     be in a consistent state (i.e., unroll() must have been called
 *     and one of the above mentioned routines should have been successfully called).
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   The probability of the evidence
 *
 *-----------------------------------------------------------------------
 */
logpr
JunctionTree::probEvidence()
{

  inference_it.set_to_last_entry();
  // first check to see if there are any cliques in the final entry,
  // and if not we use the last C rather than the last E.
  if (partitionStructureArray[inference_it.ps_i()].maxCliquesSharedStructure.size() == 0)
    --inference_it;
  
  // this next routine is not necessary since the max and sum routines
  // do not require access to random variables, all the info for the
  // scores lives in the tables.
  // setCurrentInferenceShiftTo(inference_it.pt_i());

  if (viterbiScore) {
    return partitionTableArray[inference_it.pt_i()].maxCliques[inference_it.cur_ri()].maxProb();
  } else {
    return partitionTableArray[inference_it.pt_i()].maxCliques[inference_it.cur_ri()].sumProbabilities();
  }


#if 0
  ptps_iterator ptps_it(*this);
  ptps_it.set_to_last_entry();

  if (partitionStructureArray[ptps_it.ps_i()].maxCliquesSharedStructure.size() == 0)
    --ptps_it;

  if (viterbiScore) {
    return partitionTableArray[ptps_it.pt_i()].maxCliques[ptps_it.cur_ri()].maxProb();
  } else {
    return partitionTableArray[ptps_it.pt_i()].maxCliques[ptps_it.cur_ri()].sumProbabilities();
  }
#endif

}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::printProbEvidenceAccordingToAllCliques()
 *   
 *   This routine will cycle through all cliques of all partitions and will
 *   print out the sums of the probs. of each cliuqe. Therefore, if 
 *   collect/distribute evidence has been called (and it is working)
 *   this routine should print out exactly the same value for all
 *   cliques (to within numerical precision and logp.h table error roundoff).
 *
 *   This routine really is only used for debugging/status messages.
 *
 * Preconditions:
 *   - collectEvidence() and distributeEvidence() should have been called
 *     and the partitionStructureArray array left set up.
 *
 * Postconditions:
 *   none
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
JunctionTree::printProbEvidenceAccordingToAllCliques()
{
  for (unsigned part=0;part < inference_it.pt_len(); part ++ ) {
    // next call is not needed here.
    // setCurrentInferenceShiftTo(part);
    for (unsigned cliqueNo=0;cliqueNo<partitionStructureArray[inference_it.ps_i()].maxCliquesSharedStructure.size();cliqueNo++) {
      printf("Part no %d: clique no %d: log probE = %f\n",
	     inference_it.ps_i(),cliqueNo,
	     partitionTableArray[inference_it.pt_i()].maxCliques[cliqueNo].sumProbabilities().valref());
    }
  }

#if 0
  ptps_iterator ptps_it(*this);
  ptps_it.set_to_first_entry();

  while (1) {
    for (unsigned cliqueNo=0;cliqueNo<partitionStructureArray[ptps_it.ps_i()].maxCliquesSharedStructure.size();cliqueNo++) {
      printf("Part no %d: clique no %d: log probE = %f\n",
	     ptps_it.ps_i(),cliqueNo,
	     partitionTableArray[ptps_it.pt_i()].maxCliques[cliqueNo].sumProbabilities().valref());
    }
    if (ptps_it.at_last_entry())
      break;
    else
      ++ptps_it;    
  }
#endif

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fixed memory versions of computing the probability of evidence
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::probEvidenceFixedUnroll()
 *
 *    A constant memory (i.e., indep. of T), combination of unroll and
 *    collectEvidence() above. It is constnat memory in that
 *    it keeps no more than two partitions in memory simultaneously.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) distributeEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *    same as unroll() above.
 *
 * Postconditions:
 *   The probability of the evidence is returned to the caller.
 *
 * Side Effects:
 *   This will update the space managers averages and statistics.
 *
 *   TODO: hash table usage, sharingetc. so that hash tables are not re-used 
 *         for the entire length.
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
/*
 * a version of the above routine but that only unrolls by
 * a fixed amount, there is only a constant (rather than a linear cost)
 * for the graph data structures.
 */

logpr 
JunctionTree::probEvidenceFixedUnroll(const unsigned int numFrames,
				      unsigned* numUsableFrames,
				      bool limitTime,
				      unsigned* numPartitionsDone,
				      const bool noE)
{

  FileSource *gomFS;
  // This should be safe since gmtkOnline is the only program
  // that does inference and doesn't use FileSource and gmtkOnline
  // only uses onlineFixedUnroll
  gomFS= static_cast<FileSource *>(globalObservationMatrix);
  assert(typeid(*globalObservationMatrix) == typeid(*gomFS));

  // Unroll, but do not use the long table array (2nd parameter is
  // false) for allocation, but get back the long table length
  // in a local variable for use in our iterator.
  unsigned totalNumberPartitions;
  {
    unsigned tmp = unroll(numFrames,ZeroTable,&totalNumberPartitions);
    gomFS->justifySegment(tmp);
    if (numUsableFrames) 
      *numUsableFrames = tmp;
    // limit scope of tmp.
  }
  if (numPartitionsDone)
    *numPartitionsDone = 0;

  // set up our iterator
  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this,totalNumberPartitions);

  init_CC_CE_rvs(inference_it);

  PartitionTables* prev_part_tab = NULL;
  PartitionTables* cur_part_tab
    = new PartitionTables(inference_it.cur_jt_partition());

  // we skip the first Co's LI separator if there is no P1
  // partition, since otherwise we'll get zero probability.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.skipLISeparator();
  // gather into the root of the current  partition
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   *cur_part_tab,
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i());
  // possibly print the P or C partition information
  if (inference_it.cur_part_clique_print_range() != NULL)
    printAllCliques(partitionStructureArray[inference_it.ps_i()],
		    *cur_part_tab,
		    inference_it.pt_i(),
		    inference_it.cur_nm(),
		    inference_it.cur_part_clique_print_range(),
		    stdout,
		    normalizePrintedCliques);
  // if the LI separator was turned off, we need to turn it back on.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.useLISeparator();

  for (unsigned part=1;part < inference_it.pt_len(); part ++ ) {
    delete prev_part_tab;
    prev_part_tab = cur_part_tab;

    setCurrentInferenceShiftTo(part);
    cur_part_tab
      = new PartitionTables(inference_it.cur_jt_partition());

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[inference_it.ps_prev_i()],
			  *prev_part_tab,
			  inference_it.prev_ri(),
			  inference_it.prev_nm(),
			  inference_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[inference_it.ps_i()],
			  *cur_part_tab,
			  inference_it.cur_li(),
			  inference_it.cur_nm(),
			  inference_it.pt_i());

    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // it might be that E is the first partition as well, say if this is
    // a static graph, and in this case we need in this case to skip the
    // incomming separator, which doesn't exist.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();
    if (!(inference_it.at_e() && noE)) {  
      // we only do this if we're either not at an E, or if we are at
      // an E and noE is false.

      // next, gather into the root of the final E partition
      ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		       *cur_part_tab,
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());
      // possibly print the P or C partition information
      if (inference_it.cur_part_clique_print_range() != NULL)
	printAllCliques(partitionStructureArray[inference_it.ps_i()],
			*cur_part_tab,
			inference_it.pt_i(),
			inference_it.cur_nm(),
			inference_it.cur_part_clique_print_range(),
			stdout,
			normalizePrintedCliques);
    }
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();

    // if the LI separator was turned off, we need to turn it back on.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();


    if (limitTime && probEvidenceTimeExpired)
      goto finished;
  }
  assert ( inference_it.at_e() );

 finished:

  logpr rc;
  if (inference_it.at_e()) {
    // then we finished.
    rc = cur_part_tab->maxCliques[E_root_clique].sumProbabilities();
  }
  if (numPartitionsDone)
    *numPartitionsDone = inference_it.pt_i();

  delete cur_part_tab;

  return rc;


#if 0

  // Unroll, but do not use the long table array (2nd parameter is
  // false) for allocation, but get back the long table length
  // in a local variable for use in our iterator.
  unsigned totalNumberPartitions;
  {
    unsigned tmp = unroll(numFrames,ShortTable,&totalNumberPartitions);
    if (numUsableFrames) 
      *numUsableFrames = tmp;
    // limit scope of tmp.
  }
  if (numPartitionsDone)
    *numPartitionsDone = 0;

  // set up our iterator
  ptps_iterator ptps_it(*this,totalNumberPartitions);
  ptps_it.set_to_first_entry();

  // the set of rvs that need to be adjusted, if any.
  set <RV*> two_part_rvs;
  if (ptps_it.num_c_partitions() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // partitions (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = partitionStructureArray[1].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,two_part_rvs);
  } 

  // adjustment indices to swap the two C tables for
  // reuse. These are either zero (when we have
  // less than 2 C partitions, or swap between 1 and -1.
  // the only effec the table indices.
  int prev_i_adjust = 0;
  int cur_i_adjust = 0;

  while (!ptps_it.at_last_entry()) {

    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (ptps_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // gather into the root of the current  partition
    ceGatherIntoRoot(partitionStructureArray[ptps_it.ps_i()],
		     partitionTableArray[ptps_it.ps_i()+cur_i_adjust],
		     ptps_it.cur_ri(),
		     ptps_it.cur_message_order(),
		     ptps_it.cur_nm(),
		     ptps_it.pt_i());

    // possibly print the P or C partition information
    if (ptps_it.cur_part_clique_print_range() != NULL)
      printAllCliques(partitionStructureArray[ptps_it.ps_i()],
		      partitionTableArray[ptps_it.ps_i()+cur_i_adjust],
		      ptps_it.pt_i(),
		      ptps_it.cur_nm(),
		      ptps_it.cur_part_clique_print_range(),
		      stdout,
		      normalizePrintedCliques);

    if (limitTime && probEvidenceTimeExpired)
      goto finished;

    // if the LI separator was turned off, we need to turn it back on.
    if (ptps_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();

    // advance to next partition
    ++ptps_it;

    if (ptps_it.num_c_partitions() > 2) {
      // then we need to do some RV shifting.
      if (ptps_it.at_e()) {
	// we need to shift the two C1',C2' parts back to normal, and then shift
	// the (C,E) part forward.

	// First, shift C1',C2' rvs back to normal frame
	adjustFramesBy (two_part_rvs, -(ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

	set <RV*> tmp1,tmp2;
	tmp1 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
	tmp2 = partitionStructureArray[3].returnRVsAndTheirObservedParentsAsSet();
	// We're done now with two_part_rvs, so we make a new one. The
	// next routine clears out the original set values.
	unionRVs(tmp1,tmp2,two_part_rvs);
	// now adjust C2',E' variables to match.
	adjustFramesBy (two_part_rvs, (ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

	// previous adjust needs to be what current just was
	prev_i_adjust = cur_i_adjust;
	// and turn off the current one since current is now an E
	cur_i_adjust = 0;
      } else if (ptps_it.pt_i() > 2) {

	assert ( ptps_it.at_c() );

	// shift the two C partitions over by one more
	adjustFramesBy (two_part_rvs, 1*gm_template.S*fp.numFramesInC());

	// swap settings
	if (cur_i_adjust == 0) {
	  cur_i_adjust = -1;
	  prev_i_adjust = 1;
	} else {
	  cur_i_adjust = 0;
	  prev_i_adjust = 0;
	}
      }
    }

    partitionTableArray[ptps_it.ps_i()+cur_i_adjust].init(partitionStructureArray[ptps_it.ps_i()]);

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[ptps_it.ps_prev_i()],
			  partitionTableArray[ptps_it.ps_prev_i()+prev_i_adjust],
			  ptps_it.prev_ri(),
			  ptps_it.prev_nm(),
			  ptps_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[ptps_it.ps_i()],
			  partitionTableArray[ptps_it.ps_i()+cur_i_adjust],
			  ptps_it.cur_li(),
			  ptps_it.cur_nm(),
			  ptps_it.pt_i());


  }

  assert ( ptps_it.at_e() );

  // it might be that E is the first partition as well, say if this is
  // a static graph, and in this case we need in this case to skip the
  // incomming separator, which doesn't exist.
  if (!ptps_it.has_c_partition() && P1.cliques.size() == 0)
    E1.skipLISeparator();
  // next, gather into the root of the final E partition
  if (!noE) 
    ceGatherIntoRoot(partitionStructureArray[ptps_it.ps_i()],
		     partitionTableArray[ptps_it.ps_i()],
		     ptps_it.cur_ri(),
		     ptps_it.cur_message_order(),
		     ptps_it.cur_nm(),
		     ptps_it.pt_i());
  if (!ptps_it.has_c_partition() && P1.cliques.size() == 0)
    E1.useLISeparator();

  // possibly print E partition clique info.
  if (ptps_it.cur_part_clique_print_range() != NULL)
      printAllCliques(partitionStructureArray[ptps_it.ps_i()],
		      partitionTableArray[ptps_it.ps_i()],
		      ptps_it.pt_i(),
		      ptps_it.cur_nm(),
		      ptps_it.cur_part_clique_print_range(),
		      stdout,
		      normalizePrintedCliques);

  // root clique of last partition did not do partition, since it
  // never sent to next separator (since there is none). We explicitly
  // call pruning on the root clique of the last partition.
  // curPart->maxCliques[E_root_clique].ceDoAllPruning();

 finished:

  logpr rc = partitionTableArray[ptps_it.ps_i()].maxCliques[E_root_clique].sumProbabilities();
  if (numPartitionsDone)
    *numPartitionsDone = ptps_it.pt_i();


  // lastly, shift variables back for distribute evidence.
  adjustFramesBy (two_part_rvs, -(ptps_it.num_c_partitions()-2)
		  *gm_template.S*fp.numFramesInC());

  return rc;

#endif

}




// not-quite-right DBN online filtering
logpr 
JunctionTree::onlineFixedUnroll(StreamSource *globalObservationMatrix,
				unsigned *numUsableFrames,
				unsigned *numPartitionsDone,
				const bool noE,
				FILE *f,
				const bool printObserved,
				regex_t *preg,
				char *partRangeFilter)
{
  // Unroll, but do not use the long table array (2nd parameter is
  // false) for allocation, but get back the long table length
  // in a local variable for use in our iterator.
  unsigned totalNumberPartitions;

  unsigned M = gm_template.M;
  unsigned S = gm_template.S;

#if 0
  unsigned numFramesInPprime = fp.numFramesInP() + fp.numFramesInC() * M;
  unsigned numFramesInCprime = fp.numFramesInC() * (S+M);
  unsigned numFramesInEprime = fp.numFramesInE() + fp.numFramesInC() * M;
#endif

#if 0
  unsigned numPreloadFrames = Dlinks::globalMinLag() + 
    numFramesInPprime + 2 * numFramesInCprime + numFramesInEprime;
#else
  
  if ((unsigned) Dlinks::globalMinLag() > globalObservationMatrix->startSkip()) {
    error("ERROR: -startSkip must be at least %d\n", Dlinks::globalMinLag());
  }
  if (fp.numFramesInC() == 0) {
    error("ERROR: gmtkOnline does not support empty chunks\n");
  }

  // Try to read in enough frames for P' C' C' C' E'
  // because we need enough frames in the queue to notice the
  // C' to E' transition before inference needs to start working
  // on E'. Also, if we need to reset the frame # due to overflow,
  // we want to come back to the C'C' state (rather than P'C').
  // So, we "work" on the middle C', with the following C'E' as
  // "runway" to notice the end of stream, and the preceeding C'
  // in case of frame # overflow. When we start doing smoothing,
  // we will also need to include the C's for the "psuedofuture".

  // Now we also need to preload startSkip frames so that they
  // can be skipped :) and the max Dlink lag so that the last
  // partition has sufficient future available.
  
  // The (3S + M) is because P'C' share M|C| frames, then each
  // C' requires S|C| new frames (the last C' and E' share
  // M|C| frames, but they're already counted in the C's). Here
  // |X| is the number of frames in partition X.

  unsigned tau = 0; // # of "future" C's for smoothing
  unsigned numPreloadFrames = 
    globalObservationMatrix->startSkip() + 
    fp.numFramesInP() + 
    ( (3+tau) * S + M ) * fp.numFramesInC() + 
    fp.numFramesInE() +
    Dlinks::globalMaxLag();
  // Assume the above won't over-flow with just 5 partitions
#endif

  unsigned numNewFrames = fp.numFramesInC() * S;

printf("preaload %u frames\n", numPreloadFrames);
  globalObservationMatrix->preloadFrames(numPreloadFrames);

  unsigned truePtLen = 0; // 0 until we know the true number of modified partitions
  unsigned currentMaxFrameNum;
  {
    unsigned T = globalObservationMatrix->numFrames();
    unsigned tmp;

    viterbiScore = false; // avoid allocating space for O(T) viterbi values in unroll()

    if (T > 0) {
      // We already know the length of this segment (it's probably
      // very short, since we only try to pre-load enough frames to
      // process P' C' C' C' E'), so we can pass the true number of
      // frames to unroll...
      
      tmp = unroll(T,ZeroTable,&totalNumberPartitions);
      truePtLen = totalNumberPartitions;
      currentMaxFrameNum  = T;
    } else {
      // We don't know how many frames are in the segment yet,
      // so just assume it's very long. Should be OK, since the
      // ZeroTable option just uses the minimum number of partitions.
      // We'll reset the ptps_iterator's partition length later
      // when we actually know it.
      
#ifndef UINT32_MAX
      // FIXME - this should be in stdint.h as arranged by autoconf
#  define UINT32_MAX		(4294967295U)
#endif

#define MAX_FRAME_NUMBER (1073741824U)

      tmp = unroll(MAX_FRAME_NUMBER,ZeroTable,&totalNumberPartitions);

      currentMaxFrameNum = numPreloadFrames;
    }

printf("onlineFixedUnroll: total # partitions %u\n", totalNumberPartitions);


  viterbiScore = true;  // do compute viterbi values in deScatterOutofRoot() (max-product semiring)
  
  if (numUsableFrames) 
      *numUsableFrames = tmp;
    // limit scope of tmp.
  }
  if (numPartitionsDone)
    *numPartitionsDone = 0;
  
  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this,totalNumberPartitions);

  init_CC_CE_rvs(inference_it);

  PartitionTables* prev_part_tab = NULL;
  PartitionTables* cur_part_tab
    = new PartitionTables(inference_it.cur_jt_partition());

  PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
      

  // we skip the first Co's LI separator if there is no P1
  // partition, since otherwise we'll get zero probability.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.skipLISeparator();
  // gather into the root of the current  partition
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   *cur_part_tab,
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i());

  // Set clique to most probable values given observations up to
  // the current partition
  cur_part_tab->maxCliques[inference_it.cur_ri()].
    maxProbability(ps.maxCliquesSharedStructure[inference_it.cur_ri()], true);

  // Send messages from the root clique to the rest of the cliques
  // in this partition so that they are consistant with the observations
  // in this partition. We originally wanted to send messages only to
  // the cliques actually being printed, but deScatterToOutgoingSeparators()
  // sends messages to all of a clique's outgoing separators (rather
  // than just those on the path to a printing clique) and we decided
  // not to implement a "subset" scatter. We think that in the common
  // cases there won't be much extra work from the full scatter.
  deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		     *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());

  // print filter ("Viterbi") values

  fprintf(f,"Ptn-%d P': ", inference_it.pt_i());
  if (printObserved && ps.allrvs.size() > 0) {
    printRVSetAndValues(f,ps.allrvs,true,preg);
  } else if (ps.packer.packedLen() > 0) {
    printRVSetAndValues(f,ps.hidRVVector,true,preg);
  }

  if (inference_it.cur_part_clique_print_range() != NULL) {
    viterbiScore = false; // sum-product semiring
    
    deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		       *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());
    
    // possibly print the P or C partition information
    printAllCliques(partitionStructureArray[inference_it.ps_i()],
		    *cur_part_tab,
		    inference_it.pt_i(),
		    inference_it.cur_nm(),
		    inference_it.cur_part_clique_print_range(),
		    stdout,
		    normalizePrintedCliques);
    viterbiScore = true;  // do compute viterbi values in deScatterOutofRoot() (max-product semiring)
  }


  // if the LI separator was turned off, we need to turn it back on.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.useLISeparator();

  for (unsigned part=1; part < inference_it.pt_len(); part += 1 ) {
    delete prev_part_tab;
    prev_part_tab = cur_part_tab;

    setCurrentInferenceShiftTo(part);

    cur_part_tab
      = new PartitionTables(inference_it.cur_jt_partition());

    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
      

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[inference_it.ps_prev_i()],
			  *prev_part_tab,
			  inference_it.prev_ri(),
			  inference_it.prev_nm(),
			  inference_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[inference_it.ps_i()],
			  *cur_part_tab,
			  inference_it.cur_li(),
			  inference_it.cur_nm(),
			  inference_it.pt_i());


    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // it might be that E is the first partition as well, say if this is
    // a static graph, and in this case we need in this case to skip the
    // incomming separator, which doesn't exist.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();
    if (!(inference_it.at_e() && noE)) {  
      // we only do this if we're either not at an E, or if we are at
      // an E and noE is false.

      // next, gather into the root of the final E partition
      ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		       *cur_part_tab,
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());

      cur_part_tab->maxCliques[inference_it.cur_ri()].
	maxProbability(ps.maxCliquesSharedStructure[inference_it.cur_ri()], true);

      // Send messages from the root clique to the rest of the cliques
      // in this partition so that they are consistant with the observations
      // in this partition. We originally wanted to send messages only to
      // the cliques actually being printed, but deScatterToOutgoingSeparators()
      // sends messages to all of a clique's outgoing separators (rather
      // than just those on the path to a printing clique) and we decided
      // not to implement a "subset" scatter. We think that in the common
      // cases there won't be much extra work from the full scatter.
      deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
			 *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
			 inference_it.cur_ri(),
			 inference_it.cur_message_order(),
			 inference_it.cur_nm(),
			 inference_it.pt_i());


      // print filter values

      fprintf(f,"Ptn-%d %c': ",inference_it.pt_i(), inference_it.at_e() ? 'E' : 'C');
      if (printObserved && ps.allrvs.size() > 0) {
	printRVSetAndValues(f,ps.allrvs,true,preg);
      } else if (ps.packer.packedLen() > 0) {
	printRVSetAndValues(f,ps.hidRVVector,true,preg);
      }

      // possibly print the P or C partition information
      if (inference_it.cur_part_clique_print_range() != NULL) {
	viterbiScore = false; // sum-product semiring
	
	deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
			   *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
			   inference_it.cur_ri(),
			   inference_it.cur_message_order(),
			   inference_it.cur_nm(),
			   inference_it.pt_i());

	printAllCliques(partitionStructureArray[inference_it.ps_i()],
			*cur_part_tab,
			inference_it.pt_i(),
			inference_it.cur_nm(),
			inference_it.cur_part_clique_print_range(),
			stdout,
			normalizePrintedCliques);
	viterbiScore = true;  // do compute viterbi values in deScatterOutofRoot() (max-product semiring)
      }


    }
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();

    // if the LI separator was turned off, we need to turn it back on.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();

    if (currentMaxFrameNum > MAX_FRAME_NUMBER - numNewFrames) {
      // frame number is about to overflow
      infoMsg(IM::Inference, IM::Info, "resetting frame to %u\n", numNewFrames);
      globalObservationMatrix->resetFrameNumbers(0);
      infoMsg(IM::Inference, IM::Info, "resetting ptps to partition 2\n");
      part = 1; // restart @ C'_1 (1 is C'_0, but about to increment part at top of loop)
      currentMaxFrameNum = numPreloadFrames;
    }

    // read in the same # of frames that we're about to consume to maintain
    // enough queued frames to be sure we don't overshoot the C'->E' transition
    unsigned nQueued = globalObservationMatrix->enqueueFrames(numNewFrames);

    currentMaxFrameNum += nQueued;

    // update the ptps_iterator if we just found out the true length of the segment
    if (truePtLen == 0 && globalObservationMatrix->numFrames() != 0) {
#if 0
printf("learned T=%u at %u %u\n", globalObservationMatrix->numFrames(), part, inference_it.pt_i());
#endif
      unsigned basicTempMaxUnrollAmnt;
      unsigned basicTempMinUnrollAmnt;
      int      modTempMaxUnrollAmnt;
      int      modTempMinUnrollAmnt;
      unsigned numUsableFrm;
      unsigned frmStart;
      if (!gm_template.computeUnrollParameters(globalObservationMatrix->numFrames(),
					       basicTempMaxUnrollAmnt,
					       basicTempMinUnrollAmnt,
					       modTempMaxUnrollAmnt,
					       modTempMinUnrollAmnt,
					       numUsableFrm,
					       frmStart))
	error("Segment of %d frames too short with current GMTK template of length [P=%d,C=%d,E=%d] %d frames, and M=%d,S=%d boundary parameters. Use longer utterances, different template, or decrease M,S if >1.\n",
	      globalObservationMatrix->numFrames(),
	      fp.numFramesInP(),fp.numFramesInC(),fp.numFramesInE(),
	      fp.numFrames(),
	      gm_template.M,gm_template.S);
      
      truePtLen = modTempMaxUnrollAmnt + 3;
      inference_it.set_pt_len(truePtLen);
      if (numUsableFrames) 
	*numUsableFrames = numUsableFrm;

#if 0
    if (nQueued != numFramesInCprime && nQueued != 0) {
      error("Stream segment %u length %u is incompatible with model unrolling",
	    globalObservationMatrix->segmentNumber(), 
	    globalObservationMatrix->numFrames());
    }
#endif
    }
    
    
  }
  assert ( inference_it.at_e() );

 finished:

  logpr rc;
  if (inference_it.at_e()) {
    // then we finished.
    rc = cur_part_tab->maxCliques[E_root_clique].sumProbabilities();
  }
  if (numPartitionsDone)
    *numPartitionsDone = inference_it.pt_i();

  delete cur_part_tab;

  return rc;

}



/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
