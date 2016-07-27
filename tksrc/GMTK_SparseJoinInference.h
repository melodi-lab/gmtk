/*
 * GMTK_SparseJoinInference.h
 *   Efficient Hugin-style message passing inference algorithm within sections.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_SPARSEJOININFERENCE_H
#define GMTK_SPARSEJOININFERENCE_H

#include "GMTK_SectionScheduler.h"
#include "GMTK_SectionIterator.h"
#include "GMTK_SectionInferenceAlgorithm.h"

#include "GMTK_SparseJoinSectionTables.h"

class SparseJoinInference : public SectionInferenceAlgorithm {
 public:

 SparseJoinInference(SectionScheduler *jt) : SectionInferenceAlgorithm(jt) {}


  // allocate a new SectionTables for the indicated section - free with the below
  SectionTablesBase *getSectionTables(JT_Partition& origin) { return new SparseJoinSectionTables(origin); }
  
  // free a SectionTables allocated with the above
  void releaseSectionTables(SectionTablesBase *tables) {
    if (tables) {
      SparseJoinSectionTables *t = dynamic_cast<SparseJoinSectionTables *>(tables);
      assert(t);
      delete t;
    }
  }
  
  // return a pointer to the SectionTables for section t - The memory is owned by
  // this object; do not try to free it yourself
  SectionTablesBase *getSectionTables(unsigned t);


  // All message actions are named from the perspective of C_t.

  // Prepare to compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  //   cur_section must be SparseJoinSectionTables
  void prepareForwardInterfaceSeparator(SectionTablesBase *cur_section);

  // Receive forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  //   cur_section must be SparseJoinSectionTables
  //   prev_section may be any subclass of SectionTablesBase
  void receiveForwardInterfaceSeparator(SectionTablesBase *prev_section, SectionTablesBase *cur_section);


  // Prepare to compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  //   cur_section must be SparseJoinSectionTables
  void prepareBackwardInterfaceSeparator(SectionTablesBase *cur_section);

  // Send backward message for C'_{t-1} <- C'_t (sendBackwardCrossPartitions)
  //   cur_section must be SparseJoinSectionTables
  //   prev_section may be any subclass of SectionTablesBase
  void sendBackwardInterfaceSeparator(SectionTablesBase *prev_section, SectionTablesBase *cur_section);


 private:

  void ceGatherIntoRoot(PartitionStructures &ss,
			SparseJoinSectionTables &st,
			// index of root clique in the section
			vector<unsigned> const &roots,
			// message order of the JT in this section
			vector< unsigned > &clique_order,
			// the name of the section (for debugging/status msgs)
			const char *const section_type_name,
			// number of the section in unrolled graph 
			// (for printing/debugging/status msgs only)
			const unsigned section_num,
			const bool clear_when_done = false,
			const bool also_clear_origins = false);

  void ceGatherIntoRoot(PartitionStructures &ss,
			SparseJoinSectionTables &st,
			// index of root clique in the section
			const unsigned root,
			// message order of the JT in this section
			vector< pair<unsigned,unsigned> > &message_order,
			// the name of the section (for debugging/status msgs)
			const char *const section_type_name,
			// number of the section in unrolled graph 
			// (for printing/debugging/status msgs only)
			const unsigned section_num,
			const bool clear_when_done = false,
			const bool also_clear_origins = false);

  void deScatterOutofRoot(PartitionStructures &ss,
			  SparseJoinSectionTables &st,
			  vector<unsigned> const &roots,
			  vector< pair<unsigned,unsigned> > &message_order,
			  const char *const sect_type_name,
			  const unsigned sect_num);

};

#endif

