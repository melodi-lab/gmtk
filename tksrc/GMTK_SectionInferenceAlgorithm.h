/*
 * GMTK_SectionInferenceAlgorithm.h
 *   Base class for within-section inference algorithms.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_SECTIONINFERENCEALGORITHM_H
#define GMTK_SECTIONINFERENCEALGORITHM_H


// Specific section inference algorithm implementations:

//   Pedagogical 
//   Sparse join
//   LBP

#include <assert.h>

#include "fileParser.h"

#include "GMTK_SectionSeparator.h"
#include "GMTK_PartitionTables.h"

#include "GMTK_SectionScheduler.h"
#include "GMTK_SectionIterator.h"

class SectionInferenceAlgorithm {
 public:

  SectionInferenceAlgorithm(SectionScheduler *jt) : myjt(jt) {
    assert(myjt);
  }

  virtual ~SectionInferenceAlgorithm() {}

  // All message actions are named from the perspective of C_t.

  // compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  virtual SectionSeparator *computeForwardInterfaceSeparator(SectionIterator &t, PartitionTables *section_posterior) = 0; 
//virtual void computeForwardInterfaceSeparator(unsigned t, PartitionTables *currentSection, SecionSeparator *msg) = 0; 

  // recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  virtual void receiveForwardInterfaceSeparator(SectionIterator &t, SectionSeparator *msg, PartitionTables *section_posterior) = 0;


  // compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  virtual SectionSeparator *computeBackwardsInterfaceSeparator(SectionIterator &t) = 0;

  // recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
  virtual void receiveBackwardInterfaceSeparator(SectionIterator &t, SectionSeparator const &msg) = 0;



  // return P(Q_t | X_{?}), where ? depends on the messages C_t has seen so far:
  //        P(Q_t | X_{0:t}) in the forward pass
  //        P(Q_t | X_{0:T-1}) in a (full) backward pass
  //        P(Q_t | X_{0:t+\tau}) in a smoothing backward pass

  // precondition: setCurrentInferenceShiftTo(t), section_posterior has recieved necessary messages

  virtual logpr probEvidence(SectionIterator &t, PartitionTables *section_posterior) {
  // FIXME  this should probably move to SectionScheduler ?
    if (t.at_p() && myjt->P1.cliques.size() > 0) {
      return section_posterior->maxCliques[myjt->P_ri_to_C].sumProbabilities();
    } else if (t.at_c() && myjt->Co.cliques.size() > 0) {
      return section_posterior->maxCliques[myjt->C_ri_to_C].sumProbabilities();
    } else if (t.at_e()) {
      return section_posterior->maxCliques[myjt->E_root_clique].sumProbabilities();
    } else {
      return logpr();
    }
  }


  // Section subclasses can manage their own message ordering w/in a section.
  // Read/write the section's inference plan (JT, msg orders, etc)

  void createInferencePlan(); // called by gmtkTriangulate ?
  void readInferencePlan(iDataStreamFile *f);
  void writeInferencePlan(ioDataStreamFile *f);

 protected:

  SectionScheduler *myjt;
};

#endif
