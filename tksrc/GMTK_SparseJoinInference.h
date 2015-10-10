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

#include "GMTK_SectionSeparator.h"
#include "GMTK_SectionScheduler.h"
#include "GMTK_SectionIterator.h"
#include "GMTK_SectionInferenceAlgorithm.h"

#include "GMTK_PartitionTables.h"

class SparseJoinInference : public SectionInferenceAlgorithm {
 public:

 SparseJoinInference(SectionScheduler *jt) : SectionInferenceAlgorithm(jt) {}

  // All message actions are named from the perspective of C_t.

  // compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  SectionSeparator *computeForwardInterfaceSeparator(PartitionTables *section_posterior);

  // recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  void receiveForwardInterfaceSeparator(SectionSeparator *msg, PartitionTables *section_posterior);


  // compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  SectionSeparator *computeBackwardsInterfaceSeparator(SectionIterator &t);

  // recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
  void receiveBackwardInterfaceSeparator(SectionSeparator const &msg);

 private:

  void ceGatherIntoRoot(PartitionStructures &ss,
			PartitionTables &st,
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

  void ceSendForwardsCrossSections(// previous section
				   PartitionStructures &previous_ss,
				   PartitionTables &previous_st,
				   // root clique of the previous section (i.e., the
				   // right interface clique) 
				   const unsigned previous_sect_root,
				   // name of previous section (for debugging/status msgs)
				   const char *const previous_sect_type_name,
				   // sequence number (in unrolling) of previous section
				   // (for debugging/status msgs)
				   const unsigned previous_sect_num,

				   // next section
				   PartitionStructures &next_ss,
				   PartitionTables &next_st,
				   // leaf clique of next section (i.e., index number
				   // of the left interface clique of next section)
				   const unsigned next_sect_leaf,
				   // name (debugging/status msgs)
				   const char *const next_sect_type_name,
				   // partitiiton number (debugging/status msgs)
				   const unsigned next_sect_num);

  void deScatterOutofRoot(PartitionStructures &ss,
			  PartitionTables &st,
			  const unsigned root,
			  vector< pair<unsigned,unsigned> > &message_order,
			  const char *const sect_type_name,
			  const unsigned sect_num);

  void deSendBackwardsCrossSections(PartitionStructures &previous_ss,
				      PartitionTables &previous_st,
				      const unsigned previous_section_root,
				      const char *const previous_section_type_name,
				      const unsigned previous_section_num,
				      // 
				      PartitionStructures &next_ss,
				      PartitionTables &next_st,
				      const unsigned next_section_leaf,
				      const char *const next_section_type_name,
				      const unsigned next_section_num);
};

#endif

