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


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

VCID("$Header$")

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
  if (ps.packer.packedLen() > 0)  {
    // if it is not greater than zero, then the partition has
    // no hidden discrete variables.

    if (it.at_p()) {
      ps.packer.pack(ps.hrvValuePtrs.ptr,P_partition_values.ptr);
    } else if (it.at_e()) {
      ps.packer.pack(ps.hrvValuePtrs.ptr,E_partition_values.ptr);
    } else {
      ps.packer.pack(ps.hrvValuePtrs.ptr,
		     C_partition_values.ptr
		     + 
		     (it.pt_i()-1)*ps.packer.packedLen());
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
      warning("WARNING: Part range filter must specify a valid non-zero length range within [0:%d]. Range given is %s\n",inference_it.pt_len(),partRangeFilter);
      delete partRange;
      partRange = NULL;
    }
  }
  if (partRange == NULL)
    partRange = new Range("all",0,inference_it.pt_len());

  Range::iterator* partRange_it = new Range::iterator(partRange->begin());

  int previous_C = -1;
  while (!partRange_it->at_end()) {

    unsigned part = (*partRange_it);
    setCurrentInferenceShiftTo(part);

    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

    if (ps.packer.packedLen() > 0) {
      if (inference_it.at_p()) {
	// print P partition
	fprintf(f,"Ptn-%d P': ",part);
	ps.packer.unpack(P_partition_values.ptr,ps.hrvValuePtrs.ptr);
	if (printObserved) 
	  printRVSetAndValues(f,ps.allrvs,true,preg);
	else
	  printRVSetAndValues(f,ps.hidRVVector,true,preg);
      } else if (inference_it.at_e()) {
	// print E partition
	PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
	ps.packer.unpack(E_partition_values.ptr,ps.hrvValuePtrs.ptr);
	fprintf(f,"Ptn-%d E': ",part);
	if (printObserved) 
	  printRVSetAndValues(f,ps.allrvs,true,preg);
	else
	  printRVSetAndValues(f,ps.hidRVVector,true,preg);
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
	    ps.packer.unpack(C_partition_values.ptr
			     + 
			     (inference_it.pt_i()-1)*ps.packer.packedLen(),
			     ps.hrvValuePtrs.ptr);
	    
	    fprintf(f,"Ptn-%d C': ",part);
	    if (printObserved) 
	      printRVSetAndValues(f,ps.allrvs,true,preg);
	    else
	      printRVSetAndValues(f,ps.hidRVVector,true,preg);
	  }
	previous_C = inference_it.pt_i();
      }
    }

    (*partRange_it)++;
  }

  delete partRange;

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

  // Now, do partition messages.
  for (unsigned msgNo=0;msgNo < message_order.size(); msgNo ++) {
    const unsigned from = message_order[msgNo].first;
    const unsigned to = message_order[msgNo].second;
    infoMsg(IM::Med+5,
	    "CE: gathering into %s,part[%d]: clique %d\n",
	    part_type_name,part_num,from);
    pt.maxCliques[from].
      ceGatherFromIncommingSeparators(ps.maxCliquesSharedStructure[from],
				      pt.separatorCliques,
				      ps.separatorCliquesSharedStructure.ptr);
    infoMsg(IM::Mod,
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
  // collect to partition's root clique
  infoMsg(IM::Med+5,
	  "CE: gathering into partition root %s,part[%d]: clique %d\n",
	  part_type_name,part_num,root);
  pt.maxCliques[root].
    ceGatherFromIncommingSeparators(ps.maxCliquesSharedStructure[root],
				    pt.separatorCliques,
				    ps.separatorCliquesSharedStructure.ptr);

  if (IM::messageGlb(IM::Med+9)) {
    pt.reportMemoryUsageTo(ps,stdout);
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

  infoMsg(IM::Mod,"CE: message %s,part[%d],clique(%d) --> %s,part[%d],clique(%d)\n",
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
  if (IM::messageGlb(IM::Med+9)) {
    previous_pt.reportMemoryUsageTo(previous_ps,stdout);
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

  infoMsg(IM::Med+5,"DE: distributing out of partition root %s,part[%d]: clique %d\n",
	  part_type_name,part_num,root);
  pt.maxCliques[root].
    deScatterToOutgoingSeparators(ps.maxCliquesSharedStructure[root],
				  pt.separatorCliques,
				  ps.separatorCliquesSharedStructure.ptr);
  for (unsigned msgNoP1=message_order.size();msgNoP1 > 0; msgNoP1 --) {
    const unsigned to = message_order[msgNoP1-1].first;
    const unsigned from = message_order[msgNoP1-1].second;
    infoMsg(IM::Mod,"DE: message %s,part[%d]: clique %d <-- clique %d\n",
	    part_type_name,part_num,to,from);
    pt.maxCliques[to].
      deReceiveFromIncommingSeparator(ps.maxCliquesSharedStructure[to],
				      pt.separatorCliques,
				      ps.separatorCliquesSharedStructure.ptr);

    infoMsg(IM::Med+5,"DE: distributing out of %s,part[%d]: clique %d\n",
	    part_type_name,part_num,to);
    pt.maxCliques[to].
      deScatterToOutgoingSeparators(ps.maxCliquesSharedStructure[to],
				    pt.separatorCliques,
				    ps.separatorCliquesSharedStructure.ptr);
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

  infoMsg(IM::Mod,"DE: message %s,part[%d],clique(%d) <-- %s,part[%d],clique(%d)\n",
	  previous_part_type_name,previous_part_num,previous_part_root,
	  next_part_type_name,next_part_num,next_part_leaf);
  previous_pt.maxCliques[previous_part_root].
    deReceiveFromIncommingSeparator(previous_ps.maxCliquesSharedStructure[previous_part_root],
				    next_pt.separatorCliques[next_ps.separatorCliquesSharedStructure.size()-1],
				    next_ps.separatorCliquesSharedStructure[next_ps.separatorCliquesSharedStructure.size()-1]);
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
    infoMsg(IM::High-1,
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
    infoMsg(IM::High-1,
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

  // Unroll, but do not use the long table array (2nd parameter is
  // false) for allocation, but get back the long table length
  // in a local variable for use in our iterator.
  unsigned totalNumberPartitions;
  {
    unsigned tmp = unroll(numFrames,ZeroTable,&totalNumberPartitions);
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
		    false);
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
			false);
    }
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();


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
		      false);

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
		      false);

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



/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
