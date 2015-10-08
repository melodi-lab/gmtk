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

#include <algorithm>
#include <float.h>

#include "error.h"

#include "GMTK_BoundaryTriangulate.h"

#include "GMTK_SectionScheduler.h"
#include "GMTK_SectionInferenceAlgorithm.h"

// default names of the three sections for printing/debugging messages.
const char* SectionScheduler::P1_n = "P'";
const char* SectionScheduler::Co_n = "C'";
const char* SectionScheduler::E1_n = "E'";


const char* SectionScheduler::junctionTreeMSTpriorityStr = "DSU";


void 
SectionScheduler::prepareForUnrolling(JT_Partition &section) {
}

void 
SectionScheduler::prepareForUnrolling() {
}



  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
void 
SectionScheduler::setUpDataStructures(iDataStreamFile &tri_file,
				      char const *varSectionAssignmentPrior,
				      char const *varCliqueAssignmentPrior,
				      bool checkTriFileCards)
{

  // Utilize both the partition information and elimination order
  // information already computed and contained in the file. This
  // enables the program to use external triangulation programs,
  // where this program ensures that the result is triangulated
  // and where it reports the quality of the triangulation.

  infoMsg(IM::Max,"Reading triangulation file '%s' ...\n", tri_file.fileName());
  if (!fp.readAndVerifyGMId(tri_file, checkTriFileCards)) {
    error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",
	  tri_file.fileName(),fp.fileNameParsing.c_str());
  }
  gm_template.readPartitions(tri_file);
  gm_template.readMaxCliques(tri_file);

  infoMsg(IM::Max,"Triangulating graph...\n");

  // TODO: It looks like this is really just adding the edges to make the
  //       cliques specified in the tri file actual cliques in the "graph"
  //       by adding any missing edges. 
  gm_template.triangulatePartitionsByCliqueCompletion();

  if (1) { 
    // check that graph is indeed triangulated.
    // TODO: perhaps take this check out so that inference code does
    // not need to link to the triangulation code (either that, or put
    // the triangulation check in a different file, so that we only
    // link to tri check code).

    // TODO: Post-refactor, we're not assuming the graph must be triangulated?
    //       Move this into the subset of section inference algorithms that require it.
    BoundaryTriangulate triangulator(fp,
				     gm_template.maxNumChunksInBoundary(),
				     gm_template.chunkSkip(),1.0);
    triangulator.ensurePartitionsAreChordal(gm_template);
  }

  ////////////////////////////////////////////////////////////////////
  // CREATE JUNCTION TREE DATA STRUCTURES
  infoMsg(IM::Default,"Creating Junction Tree\n"); fflush(stdout);

  setUpJTDataStructures(varSectionAssignmentPrior,varCliqueAssignmentPrior);
  prepareForUnrolling();
  infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
  ////////////////////////////////////////////////////////////////////
}



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






// create the three junction trees for the basic sections.
void 
SectionScheduler::createSectionJunctionTrees(const string pStr) {
  infoMsg(IM::Giga,"Creating of P section JT\n");
  createSectionJunctionTree(gm_template.P,pStr);
  infoMsg(IM::Giga,"Creating of C section JT\n");
  createSectionJunctionTree(gm_template.C,pStr);
  infoMsg(IM::Giga,"Creating of E section JT\n");
  createSectionJunctionTree(gm_template.E,pStr);
  infoMsg(IM::Giga,"Done creating P,C,E section JTs\n");
}



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Support for building a tree from clique graph
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////





/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createSectionJunctionTree()
 *   Create a mini-junction tree from the cliques in the given section.
 *   This uses Kruskal's greedy (but optimal) algorithm for MST generation.
 *
 *   TODO: move this routine to a MaxClique class at some point.
 *   TODO: do this during triangulation time.
 *
 * Preconditions:
 *   The section must be instantiated with cliques 
 *
 * Postconditions:
 *   The cliques in the section are now such that they
 *   form a junction tree over cliques within that section.
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
void 
SectionScheduler::createSectionJunctionTree(Section& section, const string junctionTreeMSTpriorityStr) {
  // TODO: junctionTreeMSTpriorityStr is shadowing the SectionScheduler::junctionTreeMSTpriorityStr
  //       member here. In SectionScheduler.h, the parameter is called pStr
  const unsigned numMaxCliques = section.cliques.size();

  infoMsg(IM::Giga,"Starting create JT\n");

  if (numMaxCliques == 0) {
    // Nothing to do.
    // This could happen if the partition is empty which might occur
    // for empty P's and E's. C should never be empty.
    return;
  } else if (numMaxCliques == 1) {
    // then nothing to do
    infoMsg(IM::Giga,"Partition has only one clique\n");
    return;
  } else if (numMaxCliques == 2) {
    // then JT is easy, just connect the two cliques.
    infoMsg(IM::Giga,"Partition has only two cliques\n");
    section.cliques[0].neighbors.push_back(1);
    section.cliques[1].neighbors.push_back(0);
  } else {

    infoMsg(IM::Giga,"Partition has only %d cliques\n",numMaxCliques);

    // Run max spanning tree to construct JT from set of cliques.
    // This is basically Krusgal's algorithm, but without using the
    // fast data structures (it doesn't need to be that fast
    // since it is run one time per partition, for all inference
    // runs of any length.

    // Create a vector of sets, corresponding to the trees associated
    // with each clique. Non-empty set intersection corresponds to
    // tree overlap. Each set contains the clique indices in teh set.
    vector < set<unsigned> >  findSet;
    findSet.resize(numMaxCliques);
    for (unsigned i=0;i<numMaxCliques;i++) {
      set<unsigned> iset;
      iset.insert(i);
      findSet[i] = iset;
    }

    vector< Edge > edges;
    edges.reserve((numMaxCliques*(numMaxCliques-1))/2);

    for (unsigned i=0;i<numMaxCliques;i++) {
      for (unsigned j=i+1;j<numMaxCliques;j++) {
	set<RV*> sep_set;
	set_intersection(section.cliques[i].nodes.begin(),
			 section.cliques[i].nodes.end(),
			 section.cliques[j].nodes.begin(),
			 section.cliques[j].nodes.end(),
			 inserter(sep_set,sep_set.end()));
	Edge e;
	// define the edge
	e.clique1 = i; e.clique2 = j; 
	// !!!************************!!!
	// !!!** MUST DO THIS FIRST **!!!
	// !!!************************!!!
	// First push sep set size. To get a JT, we *MUST*
	// always choose from among the cliques that
	// have the largest intersection size.
	e.weights.push_back((double)sep_set.size());
	// now that the size is there, we have many other options.

	// All gets sorted in decreasing order, so that larger values
	// have priority. Thus, the remaining items we push back in
	// the case of ties.  
	//
	// *** Larger numbers are prefered. ***
	//
	// Options to include by priority of junctionTreeMSTpriorityStr:
	//   * D: number of deterministic nodes in separator 
	//   * E: number of deterministic nodes in union of cliques.
	//   * S: neg. weight of separator
	//   * U: neg. weight of union of cliques
	//   * V: neg. frame number variance in separator
	//   * W: neg. frame number variance in union
	//   * H: number of hidden nodes in separator
	//   * O: number of observed nodes in separator
	//   * L: number of hidden nodes in union
	//   * Q: number of observed nodes in union
	//
	//  Any can be preceeded by a '-' sign to flip effect. E.g., D-E-SU
	//
	// Default case: DSU
	//
	// TODO: Other ideas for this:
	//        a) maximize number of variables in same frame (or near each other) (like variance)
	//        b) minimize number of neighbors in each clique (i.e., 
	//           if cliques already have neighbors, choose the ones with fewer.
	//        c) integrate with RV value assignment to minimize
	//           the number of unassigned clique nodes (since they're
	//           iterated over w/o knowledge of any parents. If this
	//           ends up being a search, make this be offline, in with gmtkTriangulate
	// 

	float mult = 1.0;
	for (unsigned charNo=0;charNo< junctionTreeMSTpriorityStr.size(); charNo++) {
	  const char curCase = toupper(junctionTreeMSTpriorityStr[charNo]);
	  if (curCase == '-') {
	    mult = -1.0;
	    continue;
	  }

	  if (curCase == 'D') {

	    // push back number of deterministic nodes in
	    // the separator
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = sep_set.end();
	    unsigned numDeterministicNodes = 0;
	    for (it = sep_set.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->discrete() && RV2DRV(rv)->deterministic())
		numDeterministicNodes++;
	    }
	    e.weights.push_back(mult*(double)numDeterministicNodes);

	  } else if (curCase == 'E' || curCase == 'L' || curCase == 'Q' || curCase == 'W') {
	    // push back negative weight of two cliques together.
	    set<RV*> clique_union;
	    set_union(section.cliques[i].nodes.begin(),
		      section.cliques[i].nodes.end(),
		      section.cliques[j].nodes.begin(),
		      section.cliques[j].nodes.end(),
		      inserter(clique_union,clique_union.end()));

	    if (curCase == 'E') {
	      // push back number of deterministic nodes in
	      // the union
	      set<RV*>::iterator it;
	      set<RV*>::iterator it_end = clique_union.end();
	      unsigned numDeterministicNodes = 0;
	      for (it = clique_union.begin(); it != it_end; it++) {
		RV* rv = (*it);
		if (rv->discrete() && RV2DRV(rv)->deterministic())
		  numDeterministicNodes++;
	      }
	      e.weights.push_back(mult*(double)numDeterministicNodes);
	    } else if (curCase == 'L') {
	      set<RV*>::iterator it;
	      set<RV*>::iterator it_end = clique_union.end();
	      unsigned numHidden = 0;
	      for (it = clique_union.begin(); it != it_end; it++) {
		RV* rv = (*it);
		if (rv->hidden())
		  numHidden++;
	      }
	      e.weights.push_back(mult*(double)numHidden);
	    } else if (curCase == 'Q') {
	      set<RV*>::iterator it;
	      set<RV*>::iterator it_end = clique_union.end();
	      unsigned numObserved = 0;
	      for (it = clique_union.begin(); it != it_end; it++) {
		RV* rv = (*it);
		if (rv->observed())
		  numObserved++;
	      }
	      e.weights.push_back(mult*(double)numObserved);
	    } else if (curCase == 'W') {
	      // compute frame number variance in union, push back
	      // negative to prefer smalller frame variance (i.e.,
	      // connect things that on average are close to each
	      // other in time).
	      set<RV*>::iterator it;
	      set<RV*>::iterator it_end = clique_union.end();
	      double sum = 0;
	      double sumSq = 0;
	      for (it = clique_union.begin(); it != it_end; it++) {
		RV* rv = (*it);
		sum += rv->frame();
		sumSq += rv->frame()*rv->frame();
	      }
	      double invsize = 1.0/(double)clique_union.size();
	      double variance = invsize*(sumSq - sum*sum*invsize);
	      e.weights.push_back(mult*(double)-variance);
	    }
	  } else if (curCase == 'S') {

	    // push back negative weight of separator, to prefer
	    // least negative (smallest)  weight, since larger numbers
	    // are prefered.
	    e.weights.push_back(-(double)MaxClique::computeWeight(sep_set));

	    // printf("weight of clique %d = %f, %d = %f\n",
	    // i,section.cliques[i].weight(),
	    // j,section.cliques[j].weight());

	  } else if (curCase == 'U') {

	    // push back negative weight of two cliques together.
	    set<RV*> clique_union;
	    set_union(section.cliques[i].nodes.begin(),
		      section.cliques[i].nodes.end(),
		      section.cliques[j].nodes.begin(),
		      section.cliques[j].nodes.end(),
		      inserter(clique_union,clique_union.end()));
	    e.weights.push_back(-(double)MaxClique::computeWeight(clique_union));
	  } else if (curCase == 'V') {
	    // compute frame number variance in separator, push back
	    // negative to prefer smalller frame variance (i.e.,
	    // connect things that on average are close to each
	    // other in time).
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = sep_set.end();
	    double sum = 0;
	    double sumSq = 0;
	    for (it = sep_set.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      sum += rv->frame();
	      sumSq += rv->frame()*rv->frame();
	    }
	    if (sep_set.size() == 0) {
	      e.weights.push_back(mult*(double)-FLT_MAX);
	    } else {
	      double invsize = 1.0/(double)sep_set.size();
	      double variance = invsize*(sumSq - sum*sum*invsize);
	      e.weights.push_back(mult*(double)-variance);
	    }
	  } else if (curCase == 'H') {
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = sep_set.end();
	    unsigned numHidden = 0;
	    for (it = sep_set.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->hidden())
		numHidden++;
	    }
	    e.weights.push_back(mult*(double)numHidden);
	  } else if (curCase == 'O') {
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = sep_set.end();
	    unsigned numObserved = 0;
	    for (it = sep_set.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->observed())
		numObserved++;
	    }
	    e.weights.push_back(mult*(double)numObserved);
	  } else {
	    error("ERROR: Unrecognized junction tree clique sort order letter '%c' in string '%s'\n",curCase,junctionTreeMSTpriorityStr.c_str());
	  }
	  mult = 1.0;
	}

	// add the edge.
	edges.push_back(e);
	if (IM::messageGlb(IM::Giga)) {
	  infoMsg(IM::Giga,"Edge (%d,%d) has sep size %.0f, ",
		  i,j,
		  e.weights[0]);
	  for (unsigned charNo=0;charNo< junctionTreeMSTpriorityStr.size(); charNo++) {
	    const char curCase = toupper(junctionTreeMSTpriorityStr[charNo]);
	    infoMsg(IM::Giga,"%c,weight[%d] = %f, ",curCase,charNo+1,e.weights[charNo+1]);
	  }
	  infoMsg(IM::Giga,"\n");
	}
      }
    }

    // sort in decreasing order by edge weight which in this
    // case is the sep-set size.
    sort(edges.begin(),edges.end(),EdgeCompare());

    unsigned joinsPlusOne = 1;
    for (unsigned i=0;i<edges.size();i++) {
      infoMsg(IM::Giga,"Edge %d has sep size %.0f\n",
	      i,
	      edges[i].weights[0]);

      set<unsigned>& iset1 = findSet[edges[i].clique1];
      set<unsigned>& iset2 = findSet[edges[i].clique2];
      if (iset1 != iset2) {
	// merge the two sets
	set<unsigned> new_set;
	set_union(iset1.begin(),iset1.end(),
		  iset2.begin(),iset2.end(),
		  inserter(new_set,new_set.end()));
	// make sure that all members of the set point to the
	// new set.
	set<unsigned>::iterator ns_iter;
	for (ns_iter = new_set.begin(); ns_iter != new_set.end(); ns_iter ++) {
	  const unsigned clique = *ns_iter;
	  findSet[clique] = new_set;
	}
	infoMsg(IM::Giga,"Joining cliques %d and %d (edge %d) with intersection size %.0f\n",
		edges[i].clique1,edges[i].clique2,i,edges[i].weights[0]);

	if (edges[i].weights[0] == 0.0) {
	  if (IM::messageGlb(IM::High)) {
	    // there is no way to know the difference here if the
	    // graph is non-triangualted or is simply disconnected
	    // (which is ok). A non-triangulated graph might have
	    // resulted from the user editing the trifile, but we
	    // presume that MCS has already checked for this when
	    // reading in the trifiles. We just issue an informative
	    // message just in case.
	    printf("NOTE: junction tree creation joining two cliques (%d and %d) with size 0 set intersection. Either disconnected (which is ok) or non-triangulated (which is bad) graph.\n",
		   edges[i].clique1,edges[i].clique2);
	    // TODO: print out two cliques that are trying to be joined.
	    printf("Clique %d: ",edges[i].clique1);
	    section.cliques[edges[i].clique1].printCliqueNodes(stdout);
	    printf("Clique %d: ",edges[i].clique2);
	    section.cliques[edges[i].clique2].printCliqueNodes(stdout);
	  }
	}

	section.cliques[edges[i].clique1].neighbors.push_back(edges[i].clique2);
	section.cliques[edges[i].clique2].neighbors.push_back(edges[i].clique1);

	if (++joinsPlusOne == numMaxCliques)
	  break;
      }
    }
  }
}

// routine to find the interface cliques of the sections
void 
SectionScheduler::computeSectionInterfaces() {
}

// routine to create the factors in the appropriate sections
void 
SectionScheduler::createFactorCliques() {
}

// routine to find the interface cliques of a section
void 
SectionScheduler::computeSectionInterface(JT_Partition& section1,
					  unsigned int& section1_ric,
					  JT_Partition& section2,
					  unsigned int& section2_lic,
					  bool& icliques_same)
{
}


// root the JT
void 
SectionScheduler::createDirectedGraphOfCliques() {
}

void 
SectionScheduler::createDirectedGraphOfCliques(JT_Partition& section, const unsigned root) {
}



// Assign probability giving random variables to cliques (i.e.,
// these are assigned only to cliques such that the random variables
// and *all* their parents live in the clique, plus some other
// criterion in order to make message passing as efficient as
// possible).
void 
SectionScheduler::assignRVsToCliques(const char* varSectionAssignmentPrior, const char *varCliqueAssignmentPrior) {
}

void 
SectionScheduler::assignRVsToCliques(const char *const sectionName,
				     JT_Partition&section,
				     const unsigned rootClique,
				     const char* varSectionAssignmentPrior,
				     const char *varCliqueAssignmentPrior)
{
}



void 
SectionScheduler::assignFactorsToCliques() {
}

void 
SectionScheduler::assignFactorsToCliques(JT_Partition& section) {
}


// For the three sections, set up the different message passing
// orders that are to be used. This basically just does a tree
// traversal using the previously selected root.
void 
SectionScheduler::setUpMessagePassingOrders() {
}

void 
SectionScheduler::setUpMessagePassingOrder(JT_Partition& section,
					   const unsigned root,
					   vector< pair<unsigned,unsigned> >&order,
					   const unsigned excludeFromLeafCliques,
					   vector< unsigned>& leaf_cliques) 
{
}

// Separator creation, meaning create the seperator objects
// both within and between sections. Given two neighboring
// sections L and R, the separator between the interface
// cliques in L and R is contained in R.
void 
SectionScheduler::createSeparators(JT_Partition& section, vector< pair<unsigned,unsigned> >&order) {
}

void 
SectionScheduler::createSeparators() {
}

// create the virtual evidence separators
void 
SectionScheduler::createVESeparators(JT_Partition& section) {
}


// Separator iteration order and accumulated set intersection
// creation for separator driven clique potential creation, and
// also updates the seperators partial accumulator structure and
// sets up cliques other variables.
void 
SectionScheduler::computeSeparatorIterationOrder(MaxClique& clique, JT_Partition& section) {
}

void 
SectionScheduler::computeSeparatorIterationOrders(JT_Partition& section) {
}

void 
SectionScheduler::computeSeparatorIterationOrders() {
}

// Computes the preceding iterated unassigned nodes and therein the
// set of assigned nodes in each clique that should/shouldn't be
// iterated.
void 
SectionScheduler::getCumulativeUnassignedIteratedNodes() {
}

// compute the assignment order for nodes in this
// section's cliques relative to each clique's incomming separators, and while
// doing so, also set the dispositions for each of the resulting
// nodes in each clique.
void 
SectionScheduler::sortCliqueAssignedNodesAndComputeDispositions(const char *varCliqueAssignmentPrior) {
}

void 
SectionScheduler::sortCliqueAssignedNodesAndComputeDispositions(JT_Partition& section, const char *varCliqueAssignmentPrior) {
}

 
/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::setUpJTDataStructures()
 *
 *   Sets up all the data structures in a JT (other than preparing for
 *   unrolling), and calls all the routines in the necessary
 *   order. This is like the main routine of this class, as it calls
 *   what needs to be done to set up a JT.
 *
 * Preconditions:
 *   Junction tree must have been created. Should not have called 
 *   setUpDataStructures() before. The C sections in the template
 *   have variables, but the other partions (P, and E) might be empty.
 *
 * Postconditions:
 *   All data structures are set up. The JT is ready to have
 *   prepareForUnrolling() called.
 *
 * Side Effects:
 *   Changes many internal data structures in this object. 
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
SectionScheduler::setUpJTDataStructures(const char* varSectionAssignmentPrior,
					const char *varCliqueAssignmentPrior)
{
  // main() routine for this class.
  createSectionJunctionTrees(junctionTreeMSTpriorityStr);
  computeSectionInterfaces();
  createFactorCliques();
  createDirectedGraphOfCliques();
  assignRVsToCliques(varSectionAssignmentPrior,varCliqueAssignmentPrior);
  assignFactorsToCliques();
  // TODO: assignScoringFactorsToCliques();
  setUpMessagePassingOrders();
  // create seps and VE seps.
  createSeparators();
  computeSeparatorIterationOrders();

  // -- -- used only to compute weight.
  getCumulativeUnassignedIteratedNodes(); 
  // -- --

  sortCliqueAssignedNodesAndComputeDispositions(varCliqueAssignmentPrior);
}
