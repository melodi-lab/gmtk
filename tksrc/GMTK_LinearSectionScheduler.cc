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

#if HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_HG_H
#  include "hgstamp.h"
#endif

#include "GMTK_SectionTablesBase.h"
#include "GMTK_LinearSectionScheduler.h"

#include "GMTK_FileSource.h"
#include "GMTK_StreamSource.h"


VCID(HGID)

LinearSectionScheduler::~LinearSectionScheduler() {}


logpr
LinearSectionScheduler::probEvidence(SectionInferenceAlgorithm *algorithm,
				     unsigned *numUsableFrames,
				     unsigned *numSectionsDone,
				     const bool limitTime,
				     const bool noE,                        // FIXME - needed?
				     const bool cliquePosteriorNormalize,   
				     const bool cliquePosteriorUnlog,
				     ObservationFile *posteriorFile)
// FIXME - move printing here?
{
  assert(algorithm);
  FileSource *observation_source = dynamic_cast<FileSource *>(obs_source);
  assert(observation_source);

  unsigned T; // # of sections

  unsigned nUsableFrames = unroll(observation_source->numFrames(), ZeroTable, &T);
  if (numUsableFrames) *numUsableFrames = nUsableFrames;

  SectionIterator inference_it(*this,T);

  algorithm->set_inference_it(&inference_it);
  init_CC_CE_rvs(inference_it);
  
  SectionTablesBase *prev_sect_tab = NULL;
  SectionTablesBase *cur_sect_tab = algorithm->getSectionTables(inference_it.cur_jt_section());
  
  // do P'
  algorithm->prepareForwardInterfaceSeparator(cur_sect_tab);

  // do C' C' ... E'
  unsigned t;
  for (t=1; t < T; t+=1) {

    setCurrentInferenceShiftTo(inference_it, t);

    algorithm->releaseSectionTables(prev_sect_tab);
    prev_sect_tab = cur_sect_tab;
    cur_sect_tab = algorithm->getSectionTables(inference_it.cur_jt_section());

    algorithm->receiveForwardInterfaceSeparator(prev_sect_tab, cur_sect_tab);
    algorithm->prepareForwardInterfaceSeparator(cur_sect_tab);

    //if (limitTime && probEvidenceTimeExpired) goto finished;
  }

  //finished:
  
  if (numSectionsDone) *numSectionsDone = t;
  logpr probE = cur_sect_tab->probEvidence(inference_it, *this);
  algorithm->releaseSectionTables(prev_sect_tab);
  algorithm->releaseSectionTables(cur_sect_tab);
  return probE;
}
  
logpr 
LinearSectionScheduler::forwardBackward(SectionInferenceAlgorithm *algorithm,
					unsigned *numUsableFrames,
					const bool cliquePosteriorNormalize,
					const bool cliquePosteriorUnlog,
					ObservationFile *posteriorFile,
					const bool doDistributeEvidence)
{
  assert(algorithm);
  FileSource *observation_source = dynamic_cast<FileSource *>(obs_source);
  assert(observation_source);

  unsigned T; // # of sections

  unsigned nUsableFrames = unroll(observation_source->numFrames(), LongTable, &T);
  if (numUsableFrames) *numUsableFrames = nUsableFrames;

  SectionIterator inference_it(*this,T);

  algorithm->set_inference_it(&inference_it);
  init_CC_CE_rvs(inference_it);

  // Collect Evidence (forward pass)
  
  infoMsg(IM::Low,"Collecting Evidence\n");

  SectionTablesBase *prev_section = NULL, 
                    *cur_section  = algorithm->getSectionTables(0);
  // do P'

  // FIXME - move skipLISeparator check here? 
  algorithm->prepareForwardInterfaceSeparator(cur_section);
  // FIXME - Support sectionDoDist? it would go here
  // FIXME - move useLISeparator check here? 

  // do C' C' ... E'
  unsigned t;
  for (t=1; t < T; t+=1) {

    setCurrentInferenceShiftTo(inference_it, t);

    prev_section = algorithm->getSectionTables(t-1);
    cur_section  = algorithm->getSectionTables(t);

    algorithm->receiveForwardInterfaceSeparator(prev_section, cur_section);
    // FIXME - move skipLISeparator check here? 
    algorithm->prepareForwardInterfaceSeparator(cur_section);
    // FIXME - Support sectionDoDist? it would go here
    // FIXME - move useLISeparator check here? 
  }
  assert ( inference_it.at_e() );
  infoMsg(IM::Low,"Done Collecting Evidence\n");

  logpr probE = cur_section->probEvidence(inference_it, *this);
  printf("Segment %u, after CE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	 observation_source->segmentNumber(),
	 probE.val(),
	 probE.val() / observation_source->numFrames(),
	 probE.val() / nUsableFrames);

  // Distribute Evidence (backwards pass)

  if (doDistributeEvidence) {
    // FIXME - uncomment this
    // if (SectionScheduler::viterbiScore) myjt.setRootToMaxCliqueValue(); // fix #529
    infoMsg(IM::Low,"Distributing Evidence\n");
    
    assert(t=T);
    assert(inference_it.cur_st() == T-1);
    assert(inference_it.get_st_len() == T);
    
    for (t=T-1; t > 0; --t) {
      
      setCurrentInferenceShiftTo(inference_it, t);
      prev_section = algorithm->getSectionTables(t-1);
      cur_section  = algorithm->getSectionTables(t);
      
      // FIXME - move skipLISeparator check here? 
      algorithm->prepareBackwardInterfaceSeparator(cur_section);
      // FIXME - move useLISeparator check here? 
      
#if 0
      // FIXME - enable to support Viterbi/n-best
      if (viterbiScore) {
	recordPartitionViterbiValue(inference_it);
      }
#endif
      
      algorithm->sendBackwardInterfaceSeparator(prev_section, cur_section);
    }
    setCurrentInferenceShiftTo(inference_it, 0);
    cur_section = prev_section; // work on P'
    algorithm->prepareBackwardInterfaceSeparator(cur_section);
#if 0
    // FIXME - enable to support Viterbi/n-best
    if (viterbiScore) {
      recordPartitionViterbiValue(inference_it);
    }
#endif
    
    infoMsg(IM::Low,"Done Distributing Evidence\n");
  }
    
  if (viterbiScore)
    infoMsg(IM::SoftWarning,"NOTE: Clique sums will be different since viterbi option is active\n");
  if (IM::messageGlb(IM::Low)) {
    // FIXME - uncomment this
    // myjt.printProbEvidenceAccordingToAllCliques();
    logpr probe = cur_section->probEvidence(inference_it, *this);
    printf("Segment %d, after DE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
      observation_source->segmentNumber(),
      probe.val(),
      probe.val()/observation_source->numFrames(),
      probe.val()/nUsableFrames);
  }
    
  if (p_clique_print_range || c_clique_print_range || e_clique_print_range) {
      
#if 0
    // FIXME - commented out temporarily for refactor code movement
    if (cliqueOutputName && !pCliqueFile) {
      unsigned pSize, cSize, eSize;
      myjt.cliquePosteriorSize(pSize, cSize, eSize);
      unsigned cliqueSize = (pSize > cSize) ? pSize : cSize;
      cliqueSize = (cliqueSize > eSize) ? cliqueSize : eSize;
      
      if (pPartCliquePrintRange && pSize != cliqueSize) {
         error("ERROR: incompatible cliques selected for file output\n");
      }
      if (cPartCliquePrintRange && cSize != cliqueSize) {
        error("ERROR: incompatible cliques selected for file output\n");
      }
      if (ePartCliquePrintRange && eSize != cliqueSize) {
        error("ERROR: incompatible cliques selected for file output\n");
      }
      myjt.printCliqueOrders(stdout);
      pCliqueFile = instantiateWriteFile(cliqueListName, cliqueOutputName, cliquePrintSeparator,
                                         cliquePrintFormat, cliqueSize, 0, cliquePrintSwap);
    }
#endif

    // FIXME - for heterogeneous section inference algorithms, will need to loop over
    //         section_table_array telling sections to print themselves...
    printAllCliques(stdout,cliquePosteriorNormalize,cliquePosteriorUnlog,false /* FIXME - cliquePrintOnlyEntropy */, posteriorFile);
      
#if 0
    if (pCliqueFile) {
      pCliqueFile->endOfSegment();
    }
#endif
      
  }

  return probE;
}

logpr 
LinearSectionScheduler::viterbi(SectionInferenceAlgorithm *algorithm,
				unsigned nBest,
				unsigned *numUsableFrames,
				const bool cliquePosteriorNormalize,
				const bool cliquePosteriorUnlog,
				ObservationFile *posteriorFile)
{
  assert(algorithm);
  logpr result;
  return result;
}


logpr 
LinearSectionScheduler::smoothing(SectionInferenceAlgorithm *algorithm,
				  unsigned nBest,
				  unsigned *numUsableFrames,
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
  assert(algorithm);
  logpr result;
  return result;
}
