/*
 * GMTK_PedagogicalSectionTables.h
 *   
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2009 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_PEDAGOGICALSECTIONTABLES_H
#define GMTK_PEDAGOGICALSECTIONTABLES_H

#include "GMTK_SectionTablesBase.h"
#include "GMTK_PedagogicalCliqueTable.h"

class PedagogicalSectionTables : public SectionTablesBase {
 public:

  PedagogicalCliqueTable* maxCliques;

  MaxCliqueTable *getMaxCliques() { return NULL; }

  // Create a dummy/invalid table for arrays that can be written over.
  PedagogicalSectionTables() : maxCliques(NULL) {}

  PedagogicalSectionTables(JT_Partition& origin);

  ~PedagogicalSectionTables() { clear(); }

  logpr probEvidence(SectionIterator &inference_it, SectionScheduler &myjt) {
    logpr result;
    if (inference_it.at_p() && myjt.P1.cliques.size() > 0) {
      for (unsigned i=0; i < myjt.P_ri_to_C.size(); ++i) {
	result += maxCliques[myjt.P_ri_to_C[i]].sumProbabilities();
      }
    } else if (inference_it.at_c() && myjt.Co.cliques.size() > 0) {
      for (unsigned i=0; i < myjt.C_ri_to_C.size(); ++i) {
	result += maxCliques[myjt.C_ri_to_C[i]].sumProbabilities();
      }
    } else if (inference_it.at_e()) {
      for (unsigned i=0; i < myjt.E_root_clique.size(); ++i)
         result += maxCliques[myjt.E_root_clique[i]].sumProbabilities();
    }
    return result;
  }

  void projectToOutgoingSeparators(SectionIterator &stss_it,
				   PartitionStructures &sourceSectionStructures, 
				   ConditionalSeparatorTable *separatorTableArray,
				   ConditionalSeparatorTable::SharedLocalStructure *sepSharedStructure);

  void receiveBackwardsSeparators(SectionIterator &stss_it,
				  PartitionStructures &sourceSectionStructures, 
				  ConditionalSeparatorTable *separatorTableArray,
				  ConditionalSeparatorTable::SharedLocalStructure &sepSharedStructure);
  
  void printAllCliques(FILE *f, BP_Range *clique_print_range,
		       SectionIterator &stss_it, PartitionStructures &ss,
		       const bool normalize, const bool unlog,
		       const bool justPrintEntropy,
		       ObservationFile *obs_file = NULL);

  // EM updating.
  void emIncrement(PartitionStructures& ps,
		   const logpr probE, 
		   const bool localCliqueNormalization = false,
		   const double emTrainingBeam = -LZERO);

  void clear() {
    delete [] maxCliques;
    maxCliques = NULL;
    SectionTablesBase::clear();
  }

  void init(PartitionStructures& ps);

  // memory use reporting
  void reportMemoryUsageTo(PartitionStructures& ps,FILE *f);

};


#endif

