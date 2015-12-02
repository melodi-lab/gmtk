/*
 * GMTK_SparseJoinSectionTables.h
 *   
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2009 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_SPARSEJOINSECTIONTABLES_H
#define GMTK_SPARSEJOINSECTIONTABLES_H

#include "GMTK_SectionTablesBase.h"

class SparseJoinSectionTables : public SectionTablesBase {
 public:

  MaxCliqueTable* maxCliques;

  // Create a dummy/invalid table for arrays that can be written over.
  SparseJoinSectionTables() : maxCliques(NULL) {}

  SparseJoinSectionTables(JT_Partition& origin);

  ~SparseJoinSectionTables() { clear(); }

  logpr probEvidence(SectionIterator &inference_it, SectionScheduler &myjt) {
    if (inference_it.at_p() && myjt.P1.cliques.size() > 0) {
      return maxCliques[myjt.P_ri_to_C].sumProbabilities();
    } else if (inference_it.at_c() && myjt.Co.cliques.size() > 0) {
      return maxCliques[myjt.C_ri_to_C].sumProbabilities();
    } else if (inference_it.at_e()) {
      return maxCliques[myjt.E_root_clique].sumProbabilities();
    } else {
      return logpr();
    }
  }

  void projectToOutgoingSeparators(SectionIterator &stss_it,
				   ConditionalSeparatorTable *separatorTableArray,
				   ConditionalSeparatorTable::SharedLocalStructure &sepSharedStructureArray);
  
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
    SparseJoinSectionTables::clear();
  }

  void init(PartitionStructures& ps);

  // memory use reporting
  void reportMemoryUsageTo(PartitionStructures& ps,FILE *f);

};


#endif

