/*
 * GMTK_LinearSectionScheduler.cc
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "GMTK_FileSource.h"
#include "GMTK_StreamSource.h"

#include "GMTK_BoundaryTriangulate.h"

#include "GMTK_LinearSectionScheduler.h"

  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
void 
LinearSectionScheduler::setUpDataStructures(iDataStreamFile &tri_file,
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
  myjt = new JunctionTree(gm_template);
  myjt->setUpDataStructures(varSectionAssignmentPrior,varCliqueAssignmentPrior);
  myjt->prepareForUnrolling();
  infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
  ////////////////////////////////////////////////////////////////////
  algorithm->setJT(myjt);
}

  // Formerly JunctionTree::printAllJTInfo()
void 
LinearSectionScheduler::printInferencePlanSummary(char const *fileName) {
  myjt->printAllJTInfo(fileName);
}

  // Formerly GMTemplate::reportScoreStats()
void 
LinearSectionScheduler::reportScoreStats() {
}


  // NOTE:  This assumes all inference algorithms will have clique-like things they
  //        need to print posteriors of

  // Set the range of selected clique #'s in P', C', E' for printing.
  // TODO: preconditions
void 
LinearSectionScheduler::setCliquePrintRanges(char *p_range, char *c_range, char *e_range) {
  myjt->setCliquePrintRanges(p_range,c_range,e_range);
}

void 
LinearSectionScheduler::setSectionDebugRange(Range &rng) {
  myjt->setPartitionDebugRange(rng);
}

  // Print to f the order of the variables in each clique selected by setCliquePrintRanges().
void 
LinearSectionScheduler::printCliqueOrders(FILE *f) {
  myjt->printCliqueOrders(f);
}

  // Returns the size (in # of floats) of the cliques selected by setCliquePrintRanges().
void 
LinearSectionScheduler::getCliquePosteriorSize(unsigned &p_size, unsigned &c_size, unsigned &e_size) {
  myjt->cliquePosteriorSize(p_size, c_size, e_size);
}

#if 0
// now available in base class
unsigned 
LinearSectionScheduler::unroll(unsigned numFrames,
			       const UnrollTableOptions tableOption,
			       unsigned *totalNumberSections) 
{
  unsigned numSections;
  unsigned numUsableFrames = myjt->unroll(numFrames, (JunctionTree::UnrollTableOptions)tableOption, &numSections);
  myjt->sparseJoinSegementInit(numSections);
  if (totalNumberSections) *totalNumberSections = numSections;
  return numUsableFrames;
}
#endif

logpr
LinearSectionScheduler::probEvidence(unsigned *numUsableFrames,
				     unsigned *numSectionsDone,
				     const bool limitTime,
				     const bool noE, 
				     const bool cliquePosteriorNormalize,
				     const bool cliquePosteriorUnlog,
				     ObservationFile *posteriorFile)
{
  ObservationFile *observation_file = dynamic_cast<ObservationFile *>(obs_source);
  assert(observation_file);

  unsigned T; // # of sections

  // MOVE UNROLL TO SectionScheduler
  unsigned nUsableFrames = unroll(observation_file->numFrames(), ZeroTable, &T);
  if (numUsableFrames) *numUsableFrames = nUsableFrames;
  myjt->sparseJoinSegementInit(T);
  
  PartitionTables* cur_sect_tab = myjt->getSectionTables(0);
  
 // do P'
  SectionSeparator *msg = algorithm->computeForwardInterfaceSeparator(0, cur_sect_tab);

  // do C'
  unsigned t;
  for (t=1; t < T-1; t+=1) {
    // delete cur_sect_tab;
    cur_sect_tab = myjt->getSectionTables(t);
    algorithm->receiveForwardInterfaceSeparator(t, msg, cur_sect_tab);
    delete msg;
    msg = algorithm->computeForwardInterfaceSeparator(t, cur_sect_tab);
    //if (limitTime && probEvidenceTimeExpired) goto finished;
  }

  // do E'
  cur_sect_tab = myjt->getSectionTables(t);
  algorithm->receiveForwardInterfaceSeparator(t, msg, cur_sect_tab);
  delete msg; // not if E' is the first section!
  algorithm->computeForwardInterfaceSeparator(t, cur_sect_tab);

  //finished:
  
  if (numSectionsDone) *numSectionsDone = t;
  logpr probE = algorithm->probEvidence(t, cur_sect_tab);
  delete cur_sect_tab;
  return probE;
}
  
logpr 
LinearSectionScheduler::forwardBackward(unsigned *numUsableFrames,
					const bool cliquePosteriorNormalize,
					const bool cliquePosteriorUnlog,
					ObservationFile *posteriorFile)
{
  logpr result;
  return result;
}

logpr 
LinearSectionScheduler::viterbi(unsigned *numUsableFrames,
				const bool cliquePosteriorNormalize,
				const bool cliquePosteriorUnlog,
				ObservationFile *posteriorFile)
{
  logpr result;
  return result;
}


logpr 
LinearSectionScheduler::smoothing(unsigned *numUsableFrames,
				  unsigned *numSectionsDone,
				  const bool noE,
				  FILE *f,
				  const bool printObserved,
				  regex_t *preg,
				  regex_t *creg,
				  regex_t *ereg,
				  ObservationFile *posteriorFile,
				  const bool cliquePosteriorNormalize,
				  const bool cliquePosteriorUnlog)
{
  logpr result;
  return result;
}
