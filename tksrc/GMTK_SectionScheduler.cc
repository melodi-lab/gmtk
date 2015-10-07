/*
 * GMTK_SectionScheduler.h
 *   Root class for the inference algorithms at the time series level.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "error.h"

#include "GMTK_SectionScheduler.h"
#include "GMTK_SectionInferenceAlgorithm.h"

// default names of the three sections for printing/debugging messages.
const char* SectionScheduler::P1_n = "P'";
const char* SectionScheduler::Co_n = "C'";
const char* SectionScheduler::E1_n = "E'";




/*
 *  init_CC_CE_rvs(it)
 *    given an iterator representing an unrolled segment,
 *    we create the pair of random variables CC and CE that
 *    will be shifted around in time/position in order
 *    to do inference.
 *
 * Preconditions:
 *    1) section_structure_array must be set up and initialized before
 *    this routine is called.
 *    2) All random variables should currently be shifted to their 
 *       zero position (otherwise thigns will be totally out of sync).
 *
 */
void
SectionScheduler::init_CC_CE_rvs(SectionIterator &it)
{
  if (it.num_c_sections() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // sections (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = section_structure_array[1].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = section_structure_array[2].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,cur_CC_rvs);
    tmp1 = section_structure_array[3].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,cur_CE_rvs);
  } else {
    cur_CC_rvs.clear();
    cur_CE_rvs.clear();
  } 
  cur_cc_shift = cur_ce_shift = 0;
  it.go_to_section_no(0);
}

/*
 * shiftCCtoPosition(int pos)
 *  Given an absolute section number 'pos', we shift the pair of
 *  random variables CC so that the right position (the 2nd
 *  C of CC) is now at position 'pos'. If the CE is shifted,
 *  we first shift CE back to position 0 before shifting CC.
 */
void
SectionScheduler::shiftCCtoPosition(int pos, bool reset_observed)
{
  if (cur_ce_shift != 0) {
    assert ( cur_cc_shift == 0 );
    // we need to get CE back to zero.
    adjustFramesBy (cur_CE_rvs, -cur_ce_shift * gm_template.chunkSkip()*fp.numFramesInC(),
		    reset_observed);  // gmtkOnline can't observe old values - ticket #468
    cur_ce_shift = 0;
  }
  int delta = pos - cur_cc_shift;
  adjustFramesBy (cur_CC_rvs, delta * gm_template.chunkSkip()*fp.numFramesInC());
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
SectionScheduler::shiftOriginalVarstoPosition(vector<RV*> rvs, int pos, int &prevPos)
{
  set<RV*> uprvs(rvs.begin(),rvs.end());
  int delta = (pos - prevPos);
  adjustFramesBy(uprvs, delta);
  prevPos = pos;
}

/*
 * shiftCCtoPosition(int pos)
 *  Given an absolute section number 'pos', we shift the pair of
 *  random variables CE so that the right position (the E
 *  of CE) is now at position 'pos'. If the CC pair is shifted,
 *  we first shift CC back to position 0 before shifting CE.
 */
void
SectionScheduler::shiftCEtoPosition(int pos, bool reset_observed)
{
  if (cur_cc_shift != 0) {
    assert ( cur_ce_shift == 0 );
    // we need to get CC back to zero.
    adjustFramesBy (cur_CC_rvs, -cur_cc_shift * gm_template.chunkSkip()*fp.numFramesInC(),
                    reset_observed);  // gmtkOnline can't observe old values - ticket #468  cur_cc_shift = 0;
  }
  int delta = pos - cur_ce_shift;
  adjustFramesBy (cur_CE_rvs, delta * gm_template.chunkSkip()*fp.numFramesInC());
  cur_ce_shift = pos;
}



/*
 * setCurrentInferenceShiftTo(int pos)
 *
 *  This activates sections 'pos-1' and 'pos' so the
 *  right of two successive sections is active at position 'pos'.
 *  This routine could be called:
 *
 *      activateTwoAdjacentSectionsWithRightSectionAtPos(pos)
 *      alignRightSideOfSectionPairToPos(pos)
 * 
 *  Given an absolute section number 'pos', we shift the pair of
 *  random variables (either CC or CE) so the *right* section (the 2nd
 *  C of CC or the E of CE) is now at position 'pos'. This means that
 *  messages entirely within section 'pos' and entirely within
 *  section 'pos-1' and between sections 'pos-1' and 'pos' will be
 *  correct, but no other messages are guaranteed to be correct.
 *
 * preconditions:
 *    the SectionIterator must be initilzed for the
 *    current and appropriate segment length.
 *
 * side effects
 *   - modifies the random variables corresponding to a section pair
 *   - modifies the section iterator.
 *
 */
void
SectionScheduler::setCurrentInferenceShiftTo(SectionIterator &it, int pos)
{

  if (it.at_entry(pos)) {
    // printf("== Already at shift %d\n",pos);
    return;
  }

  //   printf("========================================\n");
  //   printf("== Setting current inference shift to %d\n",pos);
  //   printf("========================================\n");

  it.go_to_section_no(pos);
  if (it.num_c_sections() <= 2) {
    // then do nothing, since nothing is never shifted away from zero.
  } else {
    // need to do some work. We need to get the right
    // of the appropriate pair (either the second C of CC in PCCE,
    // or the E of CE in PCCE) to be at position pos.
    if (it.at_p()) {
      // P has an intersection with the first C in PCCE, so we need
      // to get that first C into the right position. We do that
      // by getting both CC's into the right position.
      shiftCCtoPosition(0);
    } else if (it.at_e()) {
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





/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::unroll()
 *   sets up the internal data structures and creates a junction tree corresponding
 *   to the graph with C' unrolled k times, k >= 0. (i.e., we have
 *   a graph that has (k+1) copies of C') and place them into this
 *   class's members section_structure_array and section_table_array.
 *
 *   This also set up data structures for inference. Some of what this routine
 *   does is copy data out of STL structures into faster customized structures
 *   (such as hash tables, arrays, etc.).
 *
 *   Sets things up so that inference can take place.
 *   Unrolling is in terms of new C' sections (so that if S=1, k corresponds to frames)
 *   but in general the network will be T(P') + k*T(C') + T(E) frames long, and in general
 *   T(C') >= S. This also depends on if the boundary algorithm was run or not. 
 *
 * Preconditions:
 *   The sections must be validly instantiated with cliques, and
 *   the routine assignRVsToCliques() must have been called.
 *
 * Postconditions:
 *   The sections are such that they now know who they are connected to on
 *   the right and left.
 *   Also, unrolled data structures are set up and ready for inference for this length.
 *
 * Side Effects:
 *   Modifies all neighbors variables within the cliques within the
 *   section.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
unsigned
SectionScheduler::unroll(const unsigned int numFrames,
			 const UnrollTableOptions tableOption,
			 unsigned *totalNumberSections)
{
  // note: the argument name is numFrames, which indicates
  // the number of frames in an observation file segment.


  ////////////////////////////////////////////////////////////////////////
  // various other parameters regarding the current segment. These
  // are all computed by SectionScheduler::unroll(...)
  // 
  // 1) The number of frames given, but in units of the amount by which
  // the basic template would need to be unrolled to match the given
  // number of frames.
  unsigned basicTemplateMaxUnrollAmount;
  // 2) For the number of frames given, the min amount by which we can
  // unroll and still re-use the rvs accross time.
  unsigned basicTemplateMinUnrollAmount;
  // 3) how much the modified [P' C' E'] template would need to be
  // unrolled to match the number of frames given. This
  // value could be negative.
  int modifiedTemplateMaxUnrollAmount;
  // 4) how much the modified [P' C' E'] template would need to be
  // unrolled to match the number of frames given for the structure
  // case (i.e., this would either be [P' E'], [P' C' E'], or [P' C' C' E']
  // so this value is never more than 1 and could be negative.
  int modifiedTemplateMinUnrollAmount;
  // 5) How many usable frames are there in the current segment (for
  //    some boundaries/triangulations and their resulting sections,
  //    it might not be possible to have a segment that is a multiple
  //    of one (1) and so this number is guaranteed to be this.
  unsigned numUsableFrames;
  // 6) The offset by which we should start the observations (should be given
  // to observation matrix). Why this exists, see 5).
  unsigned frameStart;

#if 0
  // FIXME  delete this if unneeded

  // 7) The current total number of real frames.
  unsigned curNumFrames;
#endif
  ////////////////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////////////////
  // collect/distribute/probE support variables used in various routines.
  // current set of unrolled random variables
  vector <RV*> cur_unrolled_rvs;
  // current mapping from 'name+frame' to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > cur_ppf;
  ////////////////////////////////////////////////////////////////////////

  unsigned M = gm_template.maxNumChunksInBoundary();
  unsigned S = gm_template.chunkSkip();

  // first create the unrolled set of random variables corresponding to this JT.
  if (!gm_template.computeUnrollParameters(numFrames,
					   basicTemplateMaxUnrollAmount,
					   basicTemplateMinUnrollAmount,
					   modifiedTemplateMaxUnrollAmount,
					   modifiedTemplateMinUnrollAmount,
					   numUsableFrames,
					   frameStart))
    error("Segment of %d frames too short with current GMTK template of length [P=%d,C=%d,E=%d] %d frames, and M=%d,S=%d boundary parameters. Use segments of at least P+M*C+E = %d frames, different template, or decrease M,S if >1. You can identify which segment(s) in the input file(s) are too short with a command like \"obs-info -p -i1 file1 ... | awk '$2 < %d {print $1}'\"\n",
	  numFrames,
	  fp.numFramesInP(),fp.numFramesInC(),fp.numFramesInE(),
	  fp.numFrames(),
	  M,S, 
	  fp.numFramesInP() + M * fp.numFramesInC() + fp.numFramesInE(),
	  fp.numFramesInP() + M * fp.numFramesInC() + fp.numFramesInE());
  const int numCoSectionTables = modifiedTemplateMaxUnrollAmount+1;
  const int numCoSectionStructures= modifiedTemplateMinUnrollAmount+1;
  const int numSectionStructures= modifiedTemplateMinUnrollAmount+3;

  // we should never have more than 2 C structures since 2 is the max
  // necessary.
  assert ( numCoSectionStructures <= 2 );

  infoMsg(IM::Info,"Number of Current Frames = %d, Number of Currently Usable Frames = %d\n",numFrames,numUsableFrames);
  infoMsg(IM::Tiny,"Number Of Frames = %d, Unrolling Basic Template %d times, Modified Template %d times\n",
	  numFrames,
	  basicTemplateMaxUnrollAmount,
	  modifiedTemplateMaxUnrollAmount);

  // unrolled random variables
  // vector <RV*> unrolled_rvs;
  // mapping from 'name+frame' to integer index into unrolled_rvs.
  // map < RVInfo::rvParent, unsigned > ppf;

  // Only unroll the minimum amount since unrolling is expensive and
  // takes up much memory.
  fp.unroll(basicTemplateMinUnrollAmount,cur_unrolled_rvs,cur_ppf);


#if 0 

  // FIXME  this is just commented out temporarily as refactoring moves code

  // set the observed variables for now, but these may/will be modified later.
  setObservedRVs(cur_unrolled_rvs);

  prepareForNextInferenceRound();
  // this clears the shared caches in the origin cliques
  clearCliqueSepValueCache(perSegmentClearCliqueValueCache);
#endif

  // clear out the old and pre-allocate for new size.
  section_structure_array.clear();
  section_table_array.clear();

  // 
  // For the structure array, we always only need to unroll by the
  // amount corresponding to the amount by which the basic random
  // variables have been unrolled.
  section_structure_array.resize(numSectionStructures);
  unsigned partNo = 0;
  new (&section_structure_array[partNo++]) PartitionStructures(P1
							       ,cur_unrolled_rvs
							       ,cur_ppf
							       ,0*S*fp.numFramesInC(),
							       false // does not have a li separator
							       );
  for (int p=0;p<numCoSectionStructures;p++) {
    new (&section_structure_array[partNo]) PartitionStructures(Co,cur_unrolled_rvs,cur_ppf,p*S*fp.numFramesInC());
    partNo++;
  }
  new (&section_structure_array[partNo++]) 
    PartitionStructures(E1,cur_unrolled_rvs,cur_ppf,
			modifiedTemplateMinUnrollAmount*S*fp.numFramesInC());

   
  // 
  // re-allocate. We add three to modifiedTemplateUnrollAmount
  // since: 
  //  1) when we unroll by n we get n+1 copies (+1) and
  //     modifiedTemplateUnrollAmount is how much to unroll,
  //     but below is how much to allocate.
  //  2) we need a section for P, even if it is empty. (+1) 
  //  3) we need a section for E, even if it is empty. (+1)
  if (tableOption == LongTable) {
    // then we pre-allocate the table array to correspond to the
    // entire length of the segment. I.e., each table is unique.  With
    // this option, it is possible to store all of the tables in
    // memory at the same time (assuming there is sufficient RAM).
    section_table_array.resize(modifiedTemplateMaxUnrollAmount+3);
    unsigned partNo = 0;
    new (&section_table_array[partNo++]) PartitionTables(P1);
    for (int p=0;p<numCoSectionTables;p++) {
      new (&section_table_array[partNo]) PartitionTables(Co);
      partNo++;
    }
    new (&section_table_array[partNo++])
      PartitionTables(E1);
    assert (partNo == section_table_array.size());
  } else  if (tableOption == ShortTable) {
    // Then we pre-allocate the table array to correspond only to the
    // structure array, assuming that the tables are going to be
    // reused in some way. In this option, it is not possible to store
    // all the tables simultaneously in memory (meaning, we
    // are only allocating tables for *at most* a [P C C E] but
    // it might be [P C E] or even [P E] depending on the
    // number of C sections in the structure (which depends on
    // the real unrolling amount).
    section_table_array.resize(numSectionStructures);
    unsigned partNo = 0;
    new (&section_table_array[partNo++]) PartitionTables(P1);
    for (int p=0;p<numCoSectionStructures;p++) {
      new (&section_table_array[partNo]) PartitionTables(Co);
      partNo++;
    }
    new (&section_table_array[partNo++])
      PartitionTables(E1);
    assert (partNo == section_structure_array.size());
  } else if (tableOption == ZeroTable) {
    section_table_array.clear();
  } else {
    // assume NoTouchTable, so leave it alone.
  }
  
 #if 0 
  // FIXME  this is just commented out temporarily as refactoring moves code


 if (viterbiScore == true) {
    // Then we need to allocate the storage to keep the Viterbi options.
    // For now, it is unfortunately as long as the segment iself, even in the
    // island algorithm case, although the values area packed.

    // Once we include N-best list, then this is how it will effect these
    // data structures.
    const unsigned N_best = 1;

    // First do P
    P_partition_values.resize(N_best*section_structure_array[0].packer.packedLen());

    // Next C
    C_partition_values.resize(N_best
			      *section_structure_array[1].packer.packedLen()
			      *  (binaryViterbiFile ? 1 : numCoSectionTables) );

    // Next E
    E_partition_values.resize(N_best
			      *section_structure_array[numSectionStructures-1]
			      .packer.packedLen());

  }
#endif


  if (totalNumberSections != NULL)
    *totalNumberSections = (modifiedTemplateMaxUnrollAmount+3);

  return numUsableFrames;
}

