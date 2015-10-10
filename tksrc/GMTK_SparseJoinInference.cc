/*
 * GMTK_SparseJoinInference.cc
 *   Efficient Hugin-style message passing inference algorithm within sections.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "GMTK_SparseJoinInference.h"
#include "GMTK_ZeroCliqueException.h"

// compute forward message for C'_t -> C'_{t+1} (aka gather into root)
SectionSeparator *
SparseJoinInference::computeForwardInterfaceSeparator(PartitionTables *section_posterior) {
  // we skip the first Co's LI separator if there is no P1
  // section, since otherwise we'll get zero probability.
  if (inference_it->at_first_c() && myjt->P1.cliques.size() == 0){
    myjt->Co.skipLISeparator();
  }
  // gather into the root of the current  section
  ceGatherIntoRoot(myjt->section_structure_array[inference_it->cur_ss()],
		   *section_posterior,
		   inference_it->cur_ri(),
		   inference_it->cur_message_order(),
		   inference_it->cur_nm(),
		   inference_it->cur_st());
#if 0
  // TODO: support this?
  if (sectionDoDist) {
    deScatterOutofRoot(section_structure_array[inference_it->cur_ss()],
		       *section_posterior,
		       inference_it->cur_ri(),
		       inference_it->cur_message_order(),
		       inference_it->cur_nm(),
		       inference_it->cur_st());
  }
#endif
  // possibly print the P or C section information
  if (inference_it->cur_section_clique_print_range() != NULL)
    printAllCliques(myjt->section_structure_array[inference_it->cur_ss()],
		    *section_posterior,
		    inference_it->cur_st(),
		    inference_it->cur_nm(),
		    inference_it->cur_section_clique_print_range(),
		    stdout,
		    true, true, //cliquePosteriorNormalize,cliquePosteriorUnlog,
		    false, NULL /*posteriorFile*/);
  // TODO: support normalize, unlog, posterior file options above

  // if the LI separator was turned off, we need to turn it back on.
  if (inference_it->at_first_c() && myjt->P1.cliques.size() == 0)
    myjt->Co.useLISeparator();
  return section_posterior;
} 

// recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
void 
SparseJoinInference::receiveForwardInterfaceSeparator(SectionSeparator *msg, PartitionTables *section_posterior) {
}


// compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
SectionSeparator *
SparseJoinInference::computeBackwardsInterfaceSeparator(SectionIterator &t) {
  return NULL;
} 


// recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
void 
SparseJoinInference::receiveBackwardInterfaceSeparator(SectionSeparator const &msg) {
}





/*-
 *-----------------------------------------------------------------------
 * SparseJoinInference::ceGatherIntoRoot
 *   
 *   Collect Evidence Gather Into Root: This routine does a collect
 *   evidence pass for this section, and gathers all messages into
 *   the root within the current section. It does so using the
 *   message order given in the argument 'message_order', and gathers
 *   into the provided root clique.  
 *
 * See Also:
 *   Dual routine: SparseJoinInference::deScatterOutofRoot()
 *
 *
 * Preconditions:
 *   It is assumed that either:
 *     1) this is the left-most section
 *  or 2) that the left interface clique within this section has had a message
 *        sent to it from the left neighbor section.
 *
 * Postconditions:
 *     All cliques in the section have all messages but one sent to it.
 *
 * Side Effects:
 *     all sections will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void 
SparseJoinInference::ceGatherIntoRoot(PartitionStructures &ss,
				      PartitionTables &st,
				      const unsigned root,
				      vector< pair<unsigned,unsigned> > &message_order,
				      const char *const section_type_name,
				      const unsigned section_num,
				      const bool clear_when_done,
				      const bool also_clear_origins)
{
  // first check that this is not an empty section.
  if (ss.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inference_debug_level = IM::glbMsgLevel(IM::Inference);
  unsigned inference_memory_debug_level = IM::glbMsgLevel(IM::InferenceMemory);

  if (! myjt->section_debug_range.contains((int)section_num)) {
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }

  bool zeroClique = false;
  try {
    // Now, do section messages.
    for (unsigned msgNo=0;msgNo < message_order.size(); msgNo ++) {
      const unsigned from = message_order[msgNo].first;
      const unsigned to = message_order[msgNo].second;
      infoMsg(IM::Inference, IM::Med+5,
	      "CE: gathering into %s,section[%d]: clique %d\n",
	      section_type_name,section_num,from);

      // this may now throw an exception on zero clique errors - RR
      st.maxCliques[from].
	ceGatherFromIncommingSeparators(ss.maxCliquesSharedStructure[from],
					st.separatorCliques,
					ss.separatorCliquesSharedStructure.ptr);
  
      infoMsg(IM::Inference, IM::Mod,
	      "CE: message %s,section[%d]: clique %d --> clique %d\n",
	      section_type_name,section_num,from,to);
      st.maxCliques[from].
	ceSendToOutgoingSeparator(ss.maxCliquesSharedStructure[from],
				  st.separatorCliques,
				  ss.separatorCliquesSharedStructure.ptr);

      // TODO: if we are just computing probE here, we should delete
      // memory in st.maxCliques[from]. Also, if we're only doing probE,
      // we should not keep the cliques around at all, only the outgoing
      // separator.
      if (clear_when_done) {
	st.maxCliques[from].
	  clearCliqueAndIncommingSeparatorMemory(ss.maxCliquesSharedStructure[from],
						 st.separatorCliques,
						 ss.separatorCliquesSharedStructure.ptr);

	if (also_clear_origins) {
	  // then clear out the origin memory used for inference.
	  ss.origin.clearCliqueAndIncommingSeparatorMemoryForClique(from); 
	}
      }
    }
  } catch (ZeroCliqueException &e) {
    zeroClique = true; // abort this section & segment
  }
  if (!zeroClique) {
    // collect to section's root clique
    infoMsg(IM::Inference, IM::Med+5,
	    "CE: gathering into section root %s,section[%d]: clique %d\n",
	    section_type_name,section_num,root);
    try {
      st.maxCliques[root].
	ceGatherFromIncommingSeparators(ss.maxCliquesSharedStructure[root],
					st.separatorCliques,
					ss.separatorCliquesSharedStructure.ptr);
    } catch (ZeroCliqueException &e) {
      zeroClique = true; // abort this section & segment
    }
    if (!zeroClique) {
      if (IM::messageGlb(IM::InferenceMemory, IM::Med+9)) {
	st.reportMemoryUsageTo(ss,stdout);
      }
    }
  }
  if (! myjt->section_debug_range.contains((int)section_num)) {
    IM::setGlbMsgLevel(IM::InferenceMemory, inference_memory_debug_level);
    IM::setGlbMsgLevel(IM::Inference, inference_debug_level);
  }
  if (zeroClique) {
    throw ZeroCliqueException(); // continue to abort segment
  }
}


/*-
 *-----------------------------------------------------------------------
 * SparseJoinInference::ceSendForwardsCrossSections
 *   
 *   Collect Evidence Send To Next Section: This routine sends a
 *   message from the right interface clique of a left (or previous)
 *   section to the left interface clique of a right (or next)
 *   section in the section series. It is assumed that the right
 *   interface clique has had all its incomming messages sent to it.
 *
 * See Also:
 *   Dual routine: SparseJoinInference::deSendBackwardsCrossSections()
 *
 *
 * Preconditions:
 *   It is assumed that:
 *     1) the right interface of the previous section must have had
 *        all messages sent to it.
 * 
 * Postconditions:
 *     the left interface of the next section is now set up.
 *
 * Side Effects:
 *     all sections will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void 
SparseJoinInference::ceSendForwardsCrossSections(PartitionStructures &previous_ss,
						 PartitionTables &previous_st,
						 const unsigned previous_sect_root,
						 const char *const previous_sect_type_name,
						 const unsigned previous_sect_num,
						 PartitionStructures &next_ss,
						 PartitionTables &next_st,
						 const unsigned next_sect_leaf,
						 const char *const next_sect_type_name,
						 const unsigned next_sect_num)
{
  // check for empty sections.
  if (previous_ss.maxCliquesSharedStructure.size() == 0 || next_ss.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inference_debug_level = IM::glbMsgLevel(IM::Inference);
  unsigned inference_memory_debug_level = IM::glbMsgLevel(IM::InferenceMemory);

  if (! myjt->section_debug_range.contains((int)next_sect_num)) {
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }


  infoMsg(IM::Inference, IM::Mod,"CE: message %s,section[%d],clique(%d) --> %s,section[%d],clique(%d)\n",
	  previous_sect_type_name,
	  previous_sect_num,
	  previous_sect_root,
	  next_sect_type_name,
	  next_sect_num,
	  next_sect_leaf);
  previous_st.maxCliques[previous_sect_root].
    ceSendToOutgoingSeparator(previous_ss.maxCliquesSharedStructure[previous_sect_root],
			      next_st.separatorCliques[next_ss.separatorCliquesSharedStructure.size()-1],
			      next_ss.separatorCliquesSharedStructure[next_ss.separatorCliquesSharedStructure.size()-1]);

  if (IM::messageGlb(IM::InferenceMemory, IM::Med+9)) {
    previous_st.reportMemoryUsageTo(previous_ss,stdout);
  }

  if (! myjt->section_debug_range.contains((int)next_sect_num)) {
    IM::setGlbMsgLevel(IM::InferenceMemory, inference_memory_debug_level);
    IM::setGlbMsgLevel(IM::Inference, inference_debug_level);
  }
}

void 
SparseJoinInference::deScatterOutofRoot(PartitionStructures &ss,
					PartitionTables &st,
					const unsigned root,
					vector< pair<unsigned,unsigned> > &message_order,
					const char *const sect_type_name,
					const unsigned sect_num)
{
}

void 
SparseJoinInference::deSendBackwardsCrossSections(PartitionStructures &previous_ss,
						  PartitionTables &previous_st,
						  const unsigned previous_section_root,
						  const char *const previous_section_type_name,
						  const unsigned previous_section_num,
						  // 
						  PartitionStructures &next_ss,
						  PartitionTables &next_st,
						  const unsigned next_section_leaf,
						  const char *const next_section_type_name,
						  const unsigned next_section_num)
{
}
