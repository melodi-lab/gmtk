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

#include "GMTK_SectionTablesBase.h"

#include "GMTK_SectionScheduler.h"
#include "GMTK_SectionIterator.h"

class SectionInferenceAlgorithm {
 public:

  SectionInferenceAlgorithm(SectionScheduler *jt) : myjt(jt) {
    assert(myjt);
  }

  virtual ~SectionInferenceAlgorithm() {}

  void set_inference_it(SectionIterator *it) { 
    assert(it);
    inference_it = it; 
  }


  // allocate a new SectionTables for the indicated section - free with the below
  virtual SectionTablesBase *getSectionTables(JT_Partition& origin) = 0;
  
  // free a SectionTables allocated with the above
  virtual void releaseSectionTables(SectionTablesBase *tables) = 0;
  
  // return a pointer to the SectionTables for section t - The memory is owned by
  // this object; do not try to free it yourself
  virtual SectionTablesBase *getSectionTables(unsigned t) = 0;
  
  
  // All message actions are named from the perspective of current section C'_t.

  // Precondition: jt->setInferenceShiftTo(inference_it, t)


  // Prepare to compute forward message(s) for C'_t -> C'_{t+1} (aka gather into root)
  //   Precondition: C'_t's incoming (left interface) separators, if any, have been populated by
  //                 a previous call of receiveForwardInterfaceSeparator() for C'_{t-1} -> C'_t
  //   After this method, cur_section:
  //     - has processed any message(s) from C'_{t-1}
  //     - has processed all evidence in C'_t
  //     - has updated its CliqueTables scores to reflect P(Q_t | X_{0:t})
  //     - is ready to populate C'_{t+1}'s incoming (left interface) separators in the 
  //       subsequent receiveForwardInterfaceSeparator() call
  virtual void prepareForwardInterfaceSeparator(SectionTablesBase *cur_section) = 0; 

  // Receive forward message(s) for C'_{t-1} -> C'_t (aka sendForwardsCrossPartitions)
  //    prev_section populates cur_section's incoming (left interface) separators
  virtual void receiveForwardInterfaceSeparator(SectionTablesBase *prev_section, SectionTablesBase *cur_section) = 0;


  // Prepare to compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  virtual void prepareBackwardInterfaceSeparator(SectionTablesBase *cur_section) = 0;

  // send backward message for C'_{t-1} <- C'_t (sendBackwardCrossPartitions)
  virtual void sendBackwardInterfaceSeparator(SectionTablesBase *prev_section, SectionTablesBase *cur_section) = 0;


#if 0
  // return P(Q_t | X_{?}), where ? depends on the messages C_t has seen so far:
  //        P(Q_t | X_{0:t}) in the forward pass
  //        P(Q_t | X_{0:T-1}) in a (full) backward pass
  //        P(Q_t | X_{0:t+\tau}) in a smoothing backward pass

  // precondition: setCurrentInferenceShiftTo(t), section_posterior has recieved necessary messages

  virtual logpr probEvidence(SectionTablesBase *section_posterior) {
  // FIXME  This should probably move to SectionScheduler ?
    //      No, MaxCliqueTable (base class) ?
    //      No, subclasses of SectionTablesBase
    if (inference_it->at_p() && myjt->P1.cliques.size() > 0) {
      return section_posterior->maxCliques[myjt->P_ri_to_C].sumProbabilities();
    } else if (inference_it->at_c() && myjt->Co.cliques.size() > 0) {
      return section_posterior->maxCliques[myjt->C_ri_to_C].sumProbabilities();
    } else if (inference_it->at_e()) {
      return section_posterior->maxCliques[myjt->E_root_clique].sumProbabilities();
    } else {
      return logpr();
    }
  }
#endif


  virtual void printAllCliques(PartitionStructures& ss,
			       SectionTablesBase& st,
			       const unsigned section_num,
			       char const *const nm,
			       BP_Range *rng,
			       FILE *f,
			       const bool normalize,
			       const bool unlog,
			       const bool justPrintEntropy = false,
			       ObservationFile *obsFile = NULL)
  {}

  virtual void printAllCliques(const unsigned section,
			       SectionTablesBase *st,
			       FILE *f,
			       const bool normalize, const bool unlog,
			       const bool just_print_entropy = false,
			       ObservationFile *posterior_file = NULL)
  {}
  
  // Section subclasses can manage their own message ordering w/in a section.
  // Read/write the section's inference plan (JT, msg orders, etc)

  void createInferencePlan(); // called by gmtkTriangulate ?
  void readInferencePlan(iDataStreamFile *f);
  void writeInferencePlan(ioDataStreamFile *f);

 protected:

  SectionScheduler *myjt;
  SectionIterator  *inference_it;
};

#endif
