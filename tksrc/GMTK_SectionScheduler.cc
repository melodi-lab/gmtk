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



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::unroll()
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
  // are all computed by JunctionTree::unroll(...)
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

