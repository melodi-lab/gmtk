/*-
 * GMTK_JunctionTree.cc
 *     Junction Tree, message passing routines specific to the island algorithm 
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2009, < fill in later >
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

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Island Algorithm Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



PartitionTables* JunctionTree::
createPartition(const unsigned part)
{
  setCurrentInferenceShiftTo(part);
  infoMsg(IM::Mod,"$$$ createPartition: part = %d, nm = %s\n",
	  part,inference_it.cur_nm());
  if (cur_prob_evidence.not_essentially_zero())
    return new PartitionTables(inference_it.cur_jt_partition());
  else
    return NULL;
}


#if 0
/*
 * This routine saves the viterbi values computed by the most recent island run (assuming island's 
 * data structures are still valid) to stdout.
 *
 * Preconditions: 
 *    Assumes that island has just been run and all data structures are set up appropriately.
 */
void
JunctionTree::saveViterbiValuesIsland(FILE *f)
{

  fprintf(f,"Printing random variables from (P,C,E)=(%d,%d,%d) partitions\n",
	  P_partition_values.size(),
	  C_partition_values.size(),
	  E_partition_values.size());

  inference_it.set_to_first_entry();
  setCurrentInferenceShiftTo(inference_it.pt_i());
  while (!inference_it.at_last_entry()) {
    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
    if (ps.packer.packedLen() > 0) {
      if (inference_it.at_p()) {
	// print P partition
	fprintf(f,"P partition\n");
	ps.packer.unpack(P_partition_values.ptr,ps.hrvValuePtrs.ptr);
      } else {
	assert ( inference_it.at_c() );      
	// print C partition

	ps.packer.unpack(C_partition_values.ptr
			 + 
			 (inference_it.pt_i()-1)*ps.packer.packedLen(),
			 ps.hrvValuePtrs.ptr);
	fprintf(f,"C partition\n");
      }
    }
    printRVSetAndValues(f,ps.hidRVVector,true);
    ++ inference_it;
    setCurrentInferenceShiftTo(inference_it.pt_i());
  }
  PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
  assert ( inference_it.at_e() );
  if (ps.packer.packedLen() > 0) {
    ps.packer.unpack(E_partition_values.ptr,ps.hrvValuePtrs.ptr);
  }
  fprintf(f,"E partition\n");
  printRVSetAndValues(f,ps.hidRVVector,true);
  // print E partition
}
#endif


/*-
 *-----------------------------------------------------------------------
 * JunctionTree:: interface routines for island algorithm.
 *
 *    The various routines are interface routines for island algorithm
 *    to the lower level collect/distribute evidence routines defined
 *    above. These routines each only know the current partition (and
 *    possibly the previous or next partition) number, take the actual
 *    partition information from a pre-allocated array, and call
 *    the appropriate real routine. The routines included here
 *    are:
 *        ceGatherIntoRoot - call CE into the root (RI) of a partition
 *        createPartition - create the partition at location given
 *        ceSendForwardsCrossPartitions - send from RI of left partition to LI of next right partition.
 *        cePruneRootCliqueOfPartition - explicitly prune the root cliuqe of the partition 
 *        deSendBackwardsCrossPartitions - send from LI of right partition to RI of previous left partition
 *        deletePartition - free memory associated with partition, and set array ptr to NULL
 *        deScatterOutofRoot - call DE from root (RI) of a partition
 *        probEvidenceRoot - return the prob of evidence from the root clique of the partition
 *        emIncrementIsland - name says it all
 *        printAllCliques - print out all clique entries according to clique ranges.
 *
 *    Note that these routines only do real work if the current prob evidence is
 *    something other than zero. As if it is zero, that means that
 *    we've already gotten to the end and found it is zero, and we're
 *    in the mode where we're just freeing up memory.
 *
 * Preconditions:
 *   For the non create/destroy routines, each routine assumes that
 *   the array location for the corresponding partition has been
 *   assigned appropriately. The other pre-conditions are the
 *   same as the routines (above) that are called. 
 *   The create routine has the same preconditions.
 *   The destroy routine assumes that the partition at the current location
 *   has been created.
 *
 * Postconditions:
 *   Since these routiens are mainly stubs, see the corresponding function that is called.
 *
 * Side Effects:
 *   These routines will affect internal structures, space managers, hash tables
 *   within the partitions that they refer to.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */


void JunctionTree::
deleteIsland(const unsigned part)
{
  setCurrentInferenceShiftTo(part);
  infoMsg(IM::Mod,"*** deleteIsland: part = %d (%s)\n",
	  part,inference_it.cur_nm());
  map < unsigned, PartitionTables*>::iterator it;
  it = islandsMap.find(part);
  assert (it != islandsMap.end());
  assert ((*it).second != NULL);  
  // delete the partition
  delete (*it).second;
  // and remove it from the map
  islandsMap.erase(it);
}


void JunctionTree::
storeIsland(const unsigned part,
	    PartitionTables *pt)
{
  setCurrentInferenceShiftTo(part);
  infoMsg(IM::Mod,"$$$ storeIsland: part = %d, nm = %s\n",
	  part,inference_it.cur_nm());
  // check that it is not there already.
  map < unsigned, PartitionTables*>::iterator it;
  it = islandsMap.find(part);
  assert (it == islandsMap.end());
  // and add it
  islandsMap[part] = pt;
}


PartitionTables* JunctionTree::
retreiveIsland(const unsigned part)
{
  setCurrentInferenceShiftTo(part);
  infoMsg(IM::Mod,"$$$ retreiveIsland: part = %d, nm = %s\n",
	  part,inference_it.cur_nm());
  map < unsigned, PartitionTables*>::iterator it;
  it = islandsMap.find(part);
  assert (it != islandsMap.end());
  assert ((*it).second != NULL);
  return ((*it).second);
}



void JunctionTree::
ceGatherIntoRoot(const unsigned part,
		 PartitionTables *pt)
{
  setCurrentInferenceShiftTo(part);
  infoMsg(IM::Mod,"==> ceGatherIntoRoot: part = %d, nm = %s\n",
	  part,inference_it.cur_nm());

  // We check here the condition if the partition number is 1 (i.e.,
  // the 2nd partition) and P has no cliques. If this is the case,
  // then we don't want to use the partition 1's left interface
  // separator since it will be empty. Since we don't know
  // at this point if partition one is a C' or an E' we ask to skip
  // both (this will work either in the case when the number of C'
  // partitions is 0 or 1).
  if (inference_it.at_first_c() && P1.cliques.size() == 0) {
    Co.skipLISeparator();
    E1.skipLISeparator();
  }
  if (cur_prob_evidence.not_essentially_zero())
    ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		     *pt,
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());
  // restore the LI interface skip state of the partitions.
  if (inference_it.at_first_c() && P1.cliques.size() == 0) {
    Co.useLISeparator();
    E1.useLISeparator();
  }

}


void JunctionTree::
deScatterOutofRoot(const unsigned part,
		   PartitionTables* pt)
{
  setCurrentInferenceShiftTo(part);
  infoMsg(IM::Mod,"<== deScatterOutofRoot: part = %d (%s)\n",
	  part,inference_it.cur_nm());

  if (inference_it.at_first_c() && P1.cliques.size() == 0) {
    Co.skipLISeparator();
    E1.skipLISeparator();
  }

  if (cur_prob_evidence.not_essentially_zero())
    deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		       *pt,
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());

  if (inference_it.at_first_c() && P1.cliques.size() == 0) {
    Co.useLISeparator();
    E1.useLISeparator();
  }

}


void JunctionTree::
ceSendForwardsCrossPartitions(const unsigned lpart,
			      PartitionTables *lpt,
			      PartitionTables *rpt)
{
  // pre-condiitons. lpart is the left part, and rpart is the
  // partition to the right of us. The assumption here is obviously
  // that we are not at the end of the sequence, so that we can set
  // the inference shift to lpart+1.  We shift the variables so that
  // the right pair is aligned with rpt, and we shift to lpart+1
  // so that both partitions 'lpart' and 'lpart+1' are active, and we
  // need this since we're sending a message from 'lpart' to 'lpart+1'.
  setCurrentInferenceShiftTo(lpart+1);

  infoMsg(IM::Mod,"--> ceSendForwardsCrossPartitions: left part[%d] (%s) --> right part[%d] (%s)\n",
	  inference_it.pt_prev_i(),inference_it.prev_nm(),
	  inference_it.pt_i(),inference_it.cur_nm());

  if (cur_prob_evidence.not_essentially_zero())
    ceSendForwardsCrossPartitions(// left partition
				  partitionStructureArray[inference_it.ps_prev_i()],
				  *lpt,
				  inference_it.prev_ri(),
				  inference_it.prev_nm(),
				  inference_it.pt_prev_i(),
				  // right partition
				  partitionStructureArray[inference_it.ps_i()],
				  *rpt,
				  inference_it.cur_li(),
				  inference_it.cur_nm(),
				  inference_it.pt_i());

}


/*
 * deSendBackwardsCrossPartitions: Note, unlike most of the other
 * routines in this file which leave the random variables shifted so
 * that (left_part-1,left_part) are active, this routine will leave the random
 * variables shifted to (left_part,left_part+1) are active. This is becase
 * left_part is the left-partition, and we need to send a message backwards from 
 * the right partition (which is partition number 'left_part+1') to the left
 * partition (which is partition number 'left_part').
 */
void JunctionTree::
deSendBackwardsCrossPartitions(const unsigned left_part,
			       PartitionTables *lpt,
			       PartitionTables *rpt)
{

  // Some debugging messages. Don't delete in case we want
  // to re-activate them.
  // 
  // printf("<-- deSendBackwardsCrossPartitions: left_part = %d\n",left_part);
  // setCurrentInferenceShiftTo(left_part+1);
  // printf("@@@@@@@@@@@@@ state at left_part+1, left_part = %d @@@@@@@@@@@@\n",left_part);
  // inference_it.go_to_part_no(left_part+1);
  // inference_it.printState(stdout);
  // printf("@@@@@@@@@@@@@ state at left_part, left_part = %d @@@@@@@@@@@@\n",left_part);
  // inference_it.go_to_part_no(left_part);
  // inference_it.printState(stdout);
  // fflush(stdout);

  // we shift the variables so that the right pair is aligned
  // with rpt.
  setCurrentInferenceShiftTo(left_part+1);
  infoMsg(IM::Mod,
         "<-- deSendBackwardsCrossPartitions: left part[%d] (%s) <-- right part[%d] (%s)\n",
	  inference_it.pt_prev_i(),inference_it.prev_nm(),
	  inference_it.pt_i(),inference_it.cur_nm());

  if (cur_prob_evidence.not_essentially_zero())
    deSendBackwardsCrossPartitions(// left partition
				  partitionStructureArray[inference_it.ps_prev_i()],
				  *lpt,
				  inference_it.prev_ri(),
				  inference_it.prev_nm(),
				  inference_it.pt_prev_i(),
				  // right partition
				  partitionStructureArray[inference_it.ps_i()],
				  *rpt,
				  inference_it.cur_li(),
				  inference_it.cur_nm(),
				  inference_it.pt_i());

}



logpr
JunctionTree::probEvidenceRoot(const unsigned part,
			       PartitionTables* pt)
{
  setCurrentInferenceShiftTo(part);
  // return the sum of probs for the root (right interface) clique of the given partition.
  infoMsg(IM::Mod,"^^^ computing evidence for JT root: part = %d (%s)\n",
	  part,inference_it.cur_nm());
  return pt->maxCliques[inference_it.cur_ri()].sumProbabilities();
}


logpr
JunctionTree::setRootToMaxCliqueValue(const unsigned part,
				      PartitionTables* pt)
{
  setCurrentInferenceShiftTo(part);
  // return the sum of probs for the root (right interface) clique of the given partition.
  infoMsg(IM::Mod,"^^^ setting JT root to max clique value: part = %d (%s)\n",
	  part,inference_it.cur_nm());
  return pt->maxCliques[inference_it.cur_ri()].maxProbability(partitionStructureArray[inference_it.ps_i()].maxCliquesSharedStructure[inference_it.cur_ri()]);
}


void
JunctionTree::emIncrementIsland(const unsigned part,
				PartitionTables* pt,
				const logpr cur_prob_evidence,
				const bool localCliqueNormalization)
{
  setCurrentInferenceShiftTo(part);
  // increment for this partition.
  infoMsg(IM::Mod,"^^^ incrementing EM: part = %d (%s)\n",
	  part,inference_it.cur_nm());
  return pt->emIncrement(partitionStructureArray[inference_it.ps_i()],
			 cur_prob_evidence,
			 localCliqueNormalization,
			 curEMTrainingBeam);
}


#if 0

/*
 *
 * Preconditions: 
 *    The island algorithm JunctionTree::deScatterOutofRoot(const unsigned part, PartitionTables* pt)
 *    must have *just* been called and the random variables and inference_it have
 *    been set appropriately to the current partition.
 * Relies on:
 *   1) That inference_it is set to the current partition.
 *   2) That all random variables in the current structure (as designated by inference_it) are set to 
 *      the current viterbi value.
 *
 */
void
JunctionTree::storeViterbiValueIsland()
{
  PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
  if (inference_it.at_p()) {
    ps.packer.pack(ps.hrvValuePtrs.ptr,P_partition_values.ptr);
  } else if (inference_it.at_e()) {
    ps.packer.pack(ps.hrvValuePtrs.ptr,E_partition_values.ptr);
  } else {
    // 
    ps.packer.pack(ps.hrvValuePtrs.ptr,
		   C_partition_values.ptr
		   + 
		   (inference_it.pt_i()-1)*ps.packer.packedLen());
  }
}
#endif


void
JunctionTree::printAllCliques(const unsigned part,
			      PartitionTables* pt,
			      FILE* f,
			      const bool normalize,
			      const bool justPrintEntropy)
{
  setCurrentInferenceShiftTo(part);
  printAllCliques(partitionStructureArray[inference_it.ps_i()],
		  *pt,
		  inference_it.pt_i(),
		  inference_it.cur_nm(),
		  inference_it.cur_part_clique_print_range(),
		  f,
		  normalize,
		  justPrintEntropy);
}


void
JunctionTree::printAllCliqueProbabilties(const unsigned part,
					 PartitionTables* pt)
{
  setCurrentInferenceShiftTo(part);
  PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
  for (unsigned cliqueNo=0;cliqueNo<ps.maxCliquesSharedStructure.size();cliqueNo++) {
    printf("XXX Island: Part no %d: clique no %d: log probE = %f\n",
	   part,
	   cliqueNo,
	   pt->maxCliques[cliqueNo].sumProbabilities().valref());
  }
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectDistributeIslandBase()
 *
 * Preconditions:
 *   1) Guaranteed that partition[start] has
 *   been allocated. If start > 0, then we are also
 *   guaranteed that ceSendNext(start-1,start) has
 *   been called. In either case, we are also guaranteed
 *   taht ceGatherIntoRoot(start) has been called.
 *
 *   2) Also, we are guaranteed that either end == final, or that
 *   partition[end+1] has also already been created and that we've
 *   already sent a forward message to partition[end+1]. Therefore, if
 *   end < final, we are ready to do a deReceiveToPrev(end+1,end) to
 *   the partition at position end.
 *
 *   This function should only be called by collectDistributeIslandRecurse()
 *
 * Postconditions:
 *
 *  All partitions between start and end will exist and will be
 *  complete, meaning they will have had *all* messages sent, and so
 *  they could be used for EM updating. All partitions
 *  from [start+1,end] will be deallocated. 
 *
 *
 * Side Effects:
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::collectDistributeIslandBase(const unsigned start,
					  const unsigned end,
					  const bool runEMalgorithm,
					  const bool runViterbiAlgorithm,
					  const bool localCliqueNormalization)
{

  // First to through the forward part of the linear section from
  // [start,end] inclusive. Specifically: 
  // 
  // 1) partition[start] exists and already has been gathered into, but has not
  //    had a message sent to partition[start+1]. Partition [start] is an
  //    island.
  // 2) We create partitions [start+1,end] inclusive. At iteration
  //    i, we create partition (i+1). We also delete these very
  //    same partitions before we leave.
  //  3) partition 'end' might be the real end of the sequence, so we need
  //    to add a few special checks for that (i.e. don't send a message
  //    beyond the end, do EM and/or viterbi tasks, etc..
  // 4) partition 'end' might not be the real end of the sequence, in such
  //    case, 'end+1' is a partition that is an island, has already been
  //    created (it exists), and it already has had a message sent to 
  //    it from partition 'end', so we don't send to it again. Rather
  //    we can go backwards from it when we start on our backwards pass.
  //
  // The above means that we use the array for the creation and
  // storage of partitions [start+1,end] and the map is used for all
  // other partition locations. We thus create (end - start) total
  // partitions, and use location 0 for the island, so we need:
  assert ( end - start + 1 <= islandPartitionTableArray.size() );

  // we retrive the island and store it in the first array element (this therefore must
  // not be deleted, as we only delete partitions that we create).
  islandPartitionTableArray.ptr[0] = retreiveIsland(start);
  for (unsigned part = start; part <= end; part ++) {

    if (part > start) {
      // First one is created and done already (see pre-conditions),
      // and lives in the map (it is an island), so we don't do it
      // here.
      ceGatherIntoRoot(part,islandPartitionTableArray.ptr[part-start]);
    }

    if (part != end) {
      // then we create the next partition and send starting message
      // to it.
      islandPartitionTableArray.ptr[part+1 - start]
	 = createPartition(part+1);
      ceSendForwardsCrossPartitions(part,
				    islandPartitionTableArray.ptr[part - start],
				    islandPartitionTableArray.ptr[part+1 - start]);
    } else {
      // We are at the end, so we don't create and send to the next
      // partition. This is because either 1)
      // end=(partPArray.size()-1) (we are really at the last
      // partition on the right), so there is no right neighboring
      // partition, or 2) the right neighboring partition has already
      // been sent a message from the previous
      // incarnation/instantiation of the current partition.  
      // 
      // Note in earlier versions of GMTK, the clique gather code
      // didn't do pruning at the end, so we needed to explicity call
      // prunning of the root clique here (since there is no outgong
      // message). Now, we prune more eagerly at the time each clique
      // has been constructed, so there is nothing more to do here.
      // Note, however, that all the partition's root cliques need to
      // be pruned, as otherwise there might be entries in that clique
      // not contained in its previously-created outgoing separator.
    }
  }

  for (unsigned part = end; (1);) {
    if (part == (inference_it.pt_len()-1)) {
      // Do true end case separately. 
      // This is the *true end*, i.e., this is the case that we are at
      // the true final end of the sequence (so that we can finally
      // get p(E), the probability of evidence here for this
      // segment. Note that there is *no* partition to the right of us
      // in this case, so we need to treat it special case.

      cur_prob_evidence = probEvidenceRoot(part,
					   islandPartitionTableArray.ptr[part - start]);

      // We assume part is current since that was done by the
      // 'probEvidenceRoot' call above. If this routine changes
      // we'll need to uncomment the following line.
      // setCurrentInferenceShiftTo(part);
      infoMsg(IM::Low,"XXX Island Finished Inference: part = %d (%s): log probE = %f\n",
	      part,inference_it.cur_nm(),
	      cur_prob_evidence.valref());

      if (runEMalgorithm) {
	// TODO: Only initialize this if we are not doing
	// localCliqueNormalization
	if (cur_prob_evidence.essentially_zero()) {
	  infoMsg(IM::Default,"Island not training segment since probability is essentially zero\n");
	  // Note that we can't freely just jump out as we have to
	  // free up all the memory that we allocated. We thus have to
	  // check a bunch of conditions on the way out and do EM
	  // training only when appropriate, but always delete. But
	  // see below (i.e., we could jump out and have the initial
	  // caller free up all partitions in the map).
	}
      } else if (runViterbiAlgorithm) {
	if (cur_prob_evidence.essentially_zero()) {
	  infoMsg(IM::Default,"Island not decoding segment since probability is essentially zero\n");
	  // note that we can't freely just jump out as we have to
	  // free up all the memory that we allocated. We thus have to
	  // check a bunch of conditions on the way out and do
	  // decoding only when appropriate, but always delete.
	  // TODO: update: since we only have the array and the island
	  // set, we can just jump out and delete all of them, but
	  // the unwinding code is still here for now.
	} else {
	  // TODO: will change for k-best
	  setRootToMaxCliqueValue(part,
				  islandPartitionTableArray.ptr[part - start]);
	}
      }

    } else {
      // Then this is not the true end of the sequence, but this is
      // the end of an island segment (a region between two islands),
      // and there is a partition on the right that we must receive a
      // message from.

      if (part == end) {
	// then we are at the end of a segment, and to
	// the right of us is an island.
	PartitionTables* next_pt = retreiveIsland(part+1);
	deSendBackwardsCrossPartitions(part,
				       islandPartitionTableArray.ptr[part - start],
				       next_pt);
	// we don't delete the island we just retreived, only our
	// caller does that.
      } else {
	// Then  what is on the right is a partiiton that we created here and
	// we can send a normal backwards message.
	deSendBackwardsCrossPartitions(part,
				       islandPartitionTableArray.ptr[part - start],
				       islandPartitionTableArray.ptr[part+1 - start]);
	// We're now done with the partition on the right that we
	// created in this routine level, so we can safely delete it.
	delete islandPartitionTableArray.ptr[part+1 - start];
      }
    }

    // scatter into all cliques within this separator (at the very least)
    deScatterOutofRoot(part,islandPartitionTableArray.ptr[part - start]);

    // Assuming that the complete probability was not zero, we now
    // have a completed partition that we can use for EM, viterbi
    // decoding, scoring, etc.
    if (cur_prob_evidence.not_essentially_zero()) {

      if (IM::messageGlb(IM::Mod)) {
	infoMsg(IM::Mod,"!!! finished partition: part = %d (%s)\n",
		part,inference_it.cur_nm());
	printAllCliqueProbabilties(part,islandPartitionTableArray.ptr[part - start]);
      }

      // and do em updating if appropriate.
      if (runEMalgorithm) {
	emIncrementIsland(part,
			  islandPartitionTableArray.ptr[part - start],
			  cur_prob_evidence,
			  localCliqueNormalization);
      }

      // and save viterbi values if appropriate.
      if (runViterbiAlgorithm) {
	recordPartitionViterbiValue(inference_it);
      }

      printAllCliques(part,
		      islandPartitionTableArray.ptr[part - start],
		      stdout,
		      true);
    }

    if (part == start)
      break;
    part--;
  }
  // note that we do not delete islandPartitionTableArray.ptr[0] as
  // that is an island to be deleted by our caller.
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectDistributeIslandRecurse()
 *
 *  Run a collect distribute evidence stage
 *  between the partitions that are at locations [start,end] INCLUSIVE.
 *  Note if this is the very first call, it should be called
 *  with [start=0,end=(nparts-1)].
 *
 * Preconditions:
 *    For this to work, we must have that:
 *     1) partitions[start] exists (i.e., allocated), and
 *        ceGather(start) has already been called.
 *     2) Either end == partitions[final], 
 *        Or,  end != partitions[final], and
 *             partitions[end+1] has been allocated, and
 *             it is all ready to have deReceiveToPrev(end+1,end)
 *             be called.
 *    This function should only be called by collectDistributeIsland()
 *
 * Postconditions:
 *  All partitions between start and end will have at one time
 *  been in a state where all cliques will have had all
 *  messages, and so could be used for EM updating.
 *
 * Side Effects:
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::collectDistributeIslandRecurse(const unsigned start,
					     const unsigned end,
					     const unsigned base,
					     const unsigned linear_section_threshold,
					     const bool runEMalgorithm,
					     const bool runViterbiAlgorithm,
					     const bool localCliqueNormalization)
{
  // We're doing from [start,end] inclusive, so compute length
  // accordingly
  const unsigned len = end-start + 1;
  if (len <= linear_section_threshold) {
    // do base case.
    collectDistributeIslandBase(start,end,runEMalgorithm,
				runViterbiAlgorithm,
				localCliqueNormalization);
  } else { 
    const unsigned section_size = len/base;

    if (section_size <= 1) {
      infoMsg(IM::Huge,"Island collect/distribute inference, log base (%d) too large for current sect length (%d) & linear sect threshold (%d), it would result in sub sect len (%d). Backing off to linear case.\n",base,len,linear_section_threshold,section_size);
      // First, we might need to reallocate the linear size array since it it is possible that
      // we are jumping down to base before we have reached the linear threshold.
      islandPartitionTableArray.growIfNeeded(end-start+1);
      return collectDistributeIslandBase(start,end,runEMalgorithm,
					 runViterbiAlgorithm,
					 localCliqueNormalization);
    }
    // We are now assured there that section_size is at least two.

    ////////////////////////////////////////////////////////////////////
    // Step 1: we go through and do a forward pass in sections of
    // length section_size (+ a bit of extra to cover any
    // remainder). That is, we have 'base' sections, and within each
    // section we create/delete the partitions necessary to move to
    // the right, but the right edge of each section, we leave an
    // island in place in the map. Note, that we do this step only for
    // 'base-1' sections (and thus create 'base-1' islands), as the
    // right-most section is one that is created via a recursive call
    // (in step 2), and that recursive call will do the appropriate
    // forward pass starting at the right-most island that we create
    // in Step 1. In Step 2, below, we start at the right-most
    // section, and section-by-section move to the left calling this
    // routine recursively (we can move to the left in Step 2 since
    // we've got the islands in place from step 1).

    // compute what is left over at end.
    const unsigned remainder = len - section_size*base;
    unsigned section_start = start;
    // number of sections that get an extra partition.
    unsigned num_to_distribute = remainder;
    // The size of section is stored in this variable.
    unsigned cur_section_size;
    while (1) {
      // Update the current section size based on a uniform
      // allocation of the remaining partitions over the
      // first 'remainder' sections.
      cur_section_size = section_size;
      if (num_to_distribute > 0) {
	// distribute the remainder evenly over the sections.
	cur_section_size ++;
	num_to_distribute--;
      }
      // don't do last section, as is handled by recursive call.
      if (section_start + cur_section_size == (end+1))
	break;

      // the last partition included within section.
      const unsigned section_end = section_start + cur_section_size-1;
      // At this point, we are Guaranteed that:
      //  - partition[section_start] has been allocated and has
      //    already had ceGatherIntoRoot(section) called. I.e.
      //    partition[section_start] is an island and lives in the map.
      // 
      // We now run a constant space collect evidence stage, where we:
      //    a) leave partition[section_start] in place, since it is an
      //       island created from our caller, and its storage is
      //       handled by the map.
      // 
      //    b) create and delete all partitions starting from and
      //       including (section_start+1) up to and including
      //       section_end. We use temporary partition storage for
      //       these sections, since they live only for a short time
      //       (i.e., they are neither stored in in the map nor the
      //       array).
      // 
      //    c) we create one more partition, i.e., the last partition
      //       we create is section[section_end+1] and it gets a
      //       ceSendTo message, and is also gathered (since it is
      //       going to be the start for the recursive call), and we
      //       leave that one allocated in place inside the map. I.e.,
      //       this is a new island created for this section.


      PartitionTables* pt_island = retreiveIsland(section_start);
      PartitionTables* pt_next = createPartition(section_start+1);
      ceSendForwardsCrossPartitions(section_start,pt_island,pt_next);
      ceGatherIntoRoot(section_start+1,pt_next);
      for (unsigned part = section_start+1; part <= section_end; part ++) {
	PartitionTables* pt_cur = pt_next;
	pt_next = createPartition(part+1);
	ceSendForwardsCrossPartitions(part,pt_cur,pt_next);
	delete pt_cur;
	ceGatherIntoRoot(part+1,pt_next);
      }
      // Note: partition[section_end+1] is an island to be stored in
      // the map.
      storeIsland(section_end+1,pt_next);

      // move to point to next partition.
      section_start += cur_section_size;
    }

    ////////////////////////////////////////////////////////////////////
    // Step 2: Now that step 1 is done (see above), we go through each
    // section in right-to-left order, calling this routine
    // recursively. Also, since we are moving from right to left, we
    // no longer need the island partition on the right of each
    // section that we created in step1 so we delete it after the
    // recursive call.

    unsigned num_not_to_distribute = base - remainder;
    section_start = (end+1);
    bool first_iteration = true;
    while (1) {
      cur_section_size = section_size;
      if (num_not_to_distribute > 0) 
	num_not_to_distribute--;
      else {
	// distribute remainder over the first sections
	cur_section_size++;
      }
      section_start -= cur_section_size;
      // the last partition included within section.
      const unsigned section_end = section_start + cur_section_size-1;
      // recurse
      collectDistributeIslandRecurse(section_start,section_end,
				     base,linear_section_threshold,
				     runEMalgorithm,
				     runViterbiAlgorithm,
				     localCliqueNormalization);

      // We need to delete island partition at location
      // section_start+cur_section_size if it is one that we created, since
      // we will not be needing it again (i.e., the partition being deleted
      // is the one that is on the right of each of the secions that we created
      // in step 1 of this routine).
      if (!first_iteration) {
	assert (section_start+cur_section_size != (end+1));
	// this partition section_start+cur_section_size is an island to be
	// removed.
	deleteIsland(section_start+cur_section_size);
      }
      first_iteration = false;
      // if this was the first section, we end now.
      if (section_start == start)
	break;
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectDistributeIsland()
 *
 *  Do a collect evidence/distribute evidence pass using the island
 *  algorithm, a log-space version of collect/distribute evidence.
 *  The original idea was presented in the paper:
 *
 *      J. Binder, K. Murphy, and S. Russel. Space-efficient inference in
 *      dynamic probabilistic networks. In IJCAI, 1997.
 *      http://citeseer.ist.psu.edu/30635.html
 *
 *  By word of mouth, I (bilmes) heard that the idea arose from an
 *  algorithm that is used in genetic sequencing, and the connection
 *  was originally suggested by Paul Horton, but I am not sure of the
 *  original source in genetics (as of 2004. Please let me know if you
 *  are reading this and you know.).
 *
 *  Here, we adopt this idea to GMTK's notion of graph partitions.
 *
 *  Most simply, in the linear case we pay O(T) memory and O(T)
 *  compute. In the island case, we pay O(b*log_b(T)) memory and
 *  O(Tlog_b(T)) compute, where b is the base of the logarithm (note,
 *  all costs are of course multiplied by the average within-partition
 *  cost, we explain things only interms of the time/space complexity
 *  in how it relates to the length of the segment T). This can make
 *  the difference between being able to use collect/distribute
 *  evidence if in going from O(T) to O(logT) we go from *not* being
 *  able to fit in main memory to *being* able to fit in memory, but
 *  the additional time cost of logT can still be quite large. If T =
 *  1024, say, and base = 2, then we pay an extra time cost of a
 *  factor of 10!! (e.g., 1 day to 10 days).
 *
 *  More specifically, for a segment of length T, this algorithm never
 *  keeps simultaneously in memory more than about M partitions, where
 *
 *            M = T/base^k + k*(base-1)
 *
 *  and where k is the *smallest* integer value such that
 *
 *            T/base^k <= lst
 *
 *  and where lst = linear_section_threshold (the threshold at which
 *  we drop down to linear inference).  Therefore, we keep
 *  in memory no more than about
 *
 *            M = lst + k*(base-1)
 *
 *  partitions, where k is defined as above, log_b(T/lst) = k. For
 *  example, if lst == 3, and base == 3, then we keep simultaneously
 *  in memory no more than about 3+(base-1)log3(T/3) partitions,
 *  having logarithmic growth in T.
 *
 *  As mentioned above, this doesn't come for free, however, as we pay
 *  a cost in time complexity to save memory. Specifically, we do
 *  multiple collect evidence stages (unfortunately, the more time
 *  costly of the two, collect vs. distribute). Specifically, we need
 *  to do T' = R(T) collect evidence stages between partitions, where
 * 
 *        R(T) = T + (base-1)*R(T/base), 
 *
 *  This is the recurance relationship corresponding to what is
 *  implemented here in GMTK. To solve it, however, we can simplify by
 *  saying
 *  
 *        R(T) < T + base*R(T/base) 
 *
 *  meaning R(T) < T + T + ... + T, and there are log_{base}(T) terms
 *  in the sum. In actuality, there are k terms in the sum, where k is
 *  defined as above.
 * 
 *  Therefore, we can decrease running time and increase space
 *  requirements either by increasing base (from 2 on up) or
 *  alternatively (or also) by increasing lst (from 2 on up). Note
 *  that we increase mem in increasing b since O(b*log_b(T)) =
 *  O((b/ln(b))*ln(T)), so mem is growing as b/ln(b). Note also that
 *  one should set base and lst such that everything just fits in main
 *  memory (i.e., so we don't start swapping to disk), but it is not
 *  worth it to set these so that it takes any less than main memory,
 *  as that will slow things down further than necessary. Also note
 *  that argmin_{ b \in 2,3,4,... } b/ln(b) = 3, so a base of 3 is
 *  a good starting point. 
 *
 * See Also:
 *    0) collectEvidence()
 *    1) distributeEvidence()
 *    2) probEvidence()
 *
 * Preconditions:
 *   Same as unroll()
 *
 * Postconditions:
 *   Forward/backward has been run.
 *
 * Side Effects:
 *  Updates space manager statitics in cliques. 
 *  TODO: add side effects here as routine evolves.
 *
 * Results:
 *
 *
 *-----------------------------------------------------------------------
 */
logpr
JunctionTree::collectDistributeIsland(// number of frames in this segment.
				      const unsigned int numFrames,
				      // return value, number of frames
				      // that are actually used.
				      unsigned& numUsableFrames,
				      // the base of the logarithm
				      const unsigned base,
				      // the threshold at which we drop
				      // down to the linear collect/distribute
				      // evidence stage.
				      const unsigned linear_section_threshold,
				      const bool runEMalgorithm,
				      const bool runViterbiAlgorithm,
				      const bool localCliqueNormalization)
{

  // cant run both EM and viterbi at the same time.
  assert (!runEMalgorithm || !runViterbiAlgorithm);

  // must have a linear_section_threshold of at least two partitions.
  if (linear_section_threshold < 2)
    error("ERROR: Island algorithm collect/distribute inference. linear section threshold value (%d) is too small.\n",linear_section_threshold);

  // the log base must be a number that actually causes a split.
  if (base <= 1)
    error("ERROR: Island algorithm collect/distribute inference. base of log (%d) is too small.\n",base);

  unsigned totalNumberPartitions;
  numUsableFrames = unroll(numFrames,ZeroTable,&totalNumberPartitions);

  // In the island algorithm, we never hold more than the linear
  // section (stored in islandPartitionTableArray) and the island partitions
  // (stored in the map islandsMap) partition tables at a time.
  // islandPartitionTableArray is constantly being reused. Also, the
  // partitions in islandsMap are deleted (and the entry is deleted
  // from the map) as soon as we know they are no longer needed.
  // Therefore, the only object that is still (as of Jan 2009) linear
  // in the segment length is the observation vector.  
  // 
  // TODO: fix observation code to not load in entire segment at the
  // same time. Ideally, the code would have no linear dependence
  //   Right now it still does w.r.t. the obsevation sequence,
  //   the stored viteri values, and the shared clique value pool (which
  //   is not really linear in the length, but instead holds the clique
  //   values for all partitions even if they don't need to be held at all times.
  // 
  // pre-allocate the array for our needs.
  islandPartitionTableArray.resize(min(linear_section_threshold,totalNumberPartitions));

  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this,totalNumberPartitions);

  init_CC_CE_rvs(inference_it);

  if (inference_it.num_c_partitions() == 0) {
    // In this case, there is:
    //    for LI: a P'=P, E'=[CE]
    //    for RI: a P'=[P C], a E'=E
    // In either case, it might be that P is empty.
    if (P1.cliques.size() == 0)
      E1.skipLISeparator();
  }

  // Start off with unity probability of evidence (we can set it to be
  // anything other than zero actually). This will signal to the
  // island algorithm to keep going and do real work. If we reach the
  // end and if the segment decodes with zero probability, this
  // variable will be set to such, and which will cause the island
  // algorithm to, during its backward phase, just free up the islands
  // of memory that have been allocated rather than doing anything
  // else.
  cur_prob_evidence.set_to_one();

  // The recursion assumes that its first partition is already
  // allocated, so we make sure to do that here.
  PartitionTables* pt = new PartitionTables(inference_it.cur_jt_partition());
  storeIsland(0,pt);
  ceGatherIntoRoot(0,pt);
  collectDistributeIslandRecurse(0,totalNumberPartitions-1,base,linear_section_threshold,
				 runEMalgorithm,
				 runViterbiAlgorithm,
				 localCliqueNormalization);
  deleteIsland(0);

  // TODO: if we get zero probability, right now the code unwinds all
  // the way to delete the islands. Since we have all the islands here
  // in this map data structure, we don't need to do that and can jump
  // right back here to delete the islands via the map.
  assert(islandsMap.empty());

  // turn it back on in all cases.
  E1.useLISeparator();

  return cur_prob_evidence;

}


/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
