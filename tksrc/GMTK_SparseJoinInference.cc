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

#if HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_HG_H
#  include "hgstamp.h"
#endif

#include <vector>

#include "general.h"

VCID(HGID)

#include "GMTK_SparseJoinInference.h"
#include "GMTK_ZeroCliqueException.h"


SectionTablesBase *
SparseJoinInference::getSectionTables(unsigned t) { 
  if (myjt->section_table_array[t] == NULL) {
    if (t == 0) {
      myjt->section_table_array[t] = new SparseJoinSectionTables(myjt->P1);
    } else if (t < myjt->section_table_array.size() - 1) {
      myjt->section_table_array[t] = new SparseJoinSectionTables(myjt->Co);
    } else {
      myjt->section_table_array[t] = new SparseJoinSectionTables(myjt->E1);
    }
  }
  // FIXME - assert typeof(section_table_array[t] == SparseJoinInferenceSectionTable)
  return myjt->section_table_array[t];
}




// compute forward message for C'_t -> C'_{t+1} (aka gather into root)
void
SparseJoinInference::prepareForwardInterfaceSeparator(SectionTablesBase *cur_section) {
  SparseJoinSectionTables *section = dynamic_cast<SparseJoinSectionTables *>(cur_section);
  assert(section);

  // FIXME - possibly pull this up to where this method is called for performance

  // we skip the first Co's LI separator if there is no P1
  // section, since otherwise we'll get zero probability.
  if (inference_it->at_first_c() && myjt->P1.cliques.size() == 0){
    myjt->Co.skipLISeparator();
  }
  // it might be that E is the first section as well, say if this is
  // a static graph, and in this case we need in this case to skip the
  // incomming separator, which doesn't exist.
  if (!inference_it->has_c_section() && myjt->P1.cliques.size() == 0) {
    myjt->E1.skipLISeparator();
  }

  // gather into the root of the current  section
  ceGatherIntoRoot(myjt->section_structure_array[inference_it->cur_ss()],
		   *section,
		   inference_it->cur_roots(),
		   inference_it->cur_message_order(),
		   inference_it->cur_nm(),
		   inference_it->cur_st());

  // make the right interface cliques consistent if we need to
  vector< pair<unsigned,unsigned> > &subtree_reverse_messages = inference_it->cur_reverse_msg_order();

  if (subtree_reverse_messages.size() > 0) {
    deScatterOutofRoot(myjt->section_structure_array[inference_it->cur_ss()],
		       *section,
		       inference_it->cur_roots(),
		       subtree_reverse_messages,
		       inference_it->cur_nm(),
		       inference_it->cur_st());
  }

  // if the LI separator was turned off, we need to turn it back on.
  if (!inference_it->has_c_section() && myjt->P1.cliques.size() == 0) {
    myjt->E1.useLISeparator();
  }
  if (inference_it->at_first_c() && myjt->P1.cliques.size() == 0) {
    myjt->Co.useLISeparator();
  }
} 

// recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
void 
SparseJoinInference::receiveForwardInterfaceSeparator(SectionTablesBase *prev_section, SectionTablesBase *cur_section) {
  PartitionStructures &previous_ps = myjt->section_structure_array[inference_it->prev_ss()];
  vector<unsigned>     previous_part_root = inference_it->prev_ri();
  const char*const     previous_part_type_name = inference_it->prev_nm();
  unsigned             previous_part_num = inference_it->prev_st();

  PartitionStructures &next_ps = myjt->section_structure_array[inference_it->cur_ss()];
  SparseJoinSectionTables *next_st = dynamic_cast<SparseJoinSectionTables *>(cur_section);
  assert(next_st);
  vector<unsigned>     next_part_leaf = inference_it->cur_li();
  const char*const     next_part_type_name = inference_it->cur_nm();
  unsigned             next_part_num = inference_it->cur_st();

  // check for empty partitions.

  //FIXME - need API since type(prev_ps) unknown ? 
  if (previous_ps.maxCliquesSharedStructure.size() == 0 || next_ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! myjt->section_debug_range.contains( (int) next_part_num )) {
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }

#if 0
  infoMsg(IM::Inference, IM::Mod,"CE: message %s,part[%d],clique(%d) --> %s,part[%d],clique(%d)\n",
	  previous_part_type_name,
	  previous_part_num,
	  previous_part_root,
	  next_part_type_name,
	  next_part_num,
	  next_part_leaf);
#endif
  // We don't know the inference algorithm, and hence the SectionTablesBase subclass,
  // used for prev_section. But we do know that next_st is a SparseJoinSectionTables, since
  // we're running SparseJoinInference on that section. So we know how to find next_st's
  // incoming separators. Thus prev_section can figure out for itself how to project 
  // itself down onto next_st's incoming separators (all SectionInferenceAlgorithms 
  // speak the same section separator data structure).
  unsigned li_size = inference_it->cur_li().size();
  prev_section->projectToOutgoingSeparators(*inference_it, 
					    previous_ps,
					    &(next_st->separatorCliques[next_ps.separatorCliquesSharedStructure.size()-li_size]),
					    &(next_ps.separatorCliquesSharedStructure[next_ps.separatorCliquesSharedStructure.size()-li_size]));

  if (IM::messageGlb(IM::InferenceMemory, IM::Med+9)) {
    // FIXME - previous_st->reportMemoryUsageTo(previous_ps,stdout);
  }

  if (! myjt->section_debug_range.contains((int)next_part_num)) {
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }
}


// compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
void
SparseJoinInference::prepareBackwardInterfaceSeparator(SectionTablesBase *cur_section) {
  SparseJoinSectionTables *section = dynamic_cast<SparseJoinSectionTables *>(cur_section);
  assert(section);

  // FIXME - possibly pull this up to where this method is called for performance
  if (inference_it->at_first_c() && myjt->P1.cliques.size() == 0) {
    myjt->Co.skipLISeparator();    
  } else if (!inference_it->has_c_section() && myjt->P1.cliques.size() == 0) {
    myjt->E1.skipLISeparator();
  }
  
#if 0
  deScatterOutofRoot(myjt->section_structure_array[inference_it->cur_ss()],
		     *section,
		     inference_it->cur_ri(),
		     inference_it->cur_message_order(),
		     inference_it->cur_nm(),
		     inference_it->cur_st());
#endif
  if (inference_it->at_first_c() && myjt->P1.cliques.size() == 0) {
    myjt->Co.useLISeparator();
  } else if (!inference_it->has_c_section() && myjt->P1.cliques.size() == 0) {
    myjt->E1.useLISeparator();
  }

#if 0
  // FIXME - add API to get/set Viterbi values as sets of RVs
  if (myjt->viterbiScore)
    recordPartitionViterbiValue(inference_it);
#endif
} 


// send backward message for C'_{t-1} <- C'_t (sendBackwardCrossPartitions)
void 
SparseJoinInference::sendBackwardInterfaceSeparator(SectionTablesBase *prev_section, SectionTablesBase *cur_section) {
  PartitionStructures &previous_ps = myjt->section_structure_array[inference_it->prev_ss()];
  vector<unsigned>     previous_part_root = inference_it->prev_ri();
  const char*const     previous_part_type_name = inference_it->prev_nm();
  unsigned             previous_part_num = inference_it->prev_st();

  PartitionStructures &next_ps = myjt->section_structure_array[inference_it->cur_ss()];
  SparseJoinSectionTables *next_st = dynamic_cast<SparseJoinSectionTables *>(cur_section);
  assert(next_st);

  vector<unsigned>     next_part_leaf = inference_it->cur_li();
  const char*const     next_part_type_name = inference_it->cur_nm();
  unsigned             next_part_num = inference_it->cur_st();

  // check for empty partitions.

  // FIXME - see receive...
  if (previous_ps.maxCliquesSharedStructure.size() == 0 || next_ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! myjt->section_debug_range.contains((int)previous_part_num)) {
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }
#if 0
  infoMsg(IM::Inference, IM::Mod,"DE: message %s,part[%d],clique(%d) <-- %s,part[%d],clique(%d)\n",
	  previous_part_type_name,previous_part_num,previous_part_root,
	  next_part_type_name,next_part_num,next_part_leaf);
#endif
#if 0
  // FIXME - call prev_section->projectToIncomingSeporator()
  prevous_st->maxCliques[previous_part_root].
    deReceiveFromIncommingSeparator(previous_ps.maxCliquesSharedStructure[previous_part_root],
				    next_st->separatorCliques[next_ps.separatorCliquesSharedStructure.size()-1],
				    next_ps.separatorCliquesSharedStructure[next_ps.separatorCliquesSharedStructure.size()-1]);
#else
  // We don't know the inference algorithm, and hence the SectionTablesBase subclass,
  // used for prev_section. But we do know that next_st is a SparseJoinSectionTables, since
  // we're running SparseJoinInference on that section. So we know how to find next_st's
  // outgoing (for the backwards pass) separators. Thus prev_section can figure out for itself 
  // how to project receive next_st's interface separators (all SectionInferenceAlgorithms 
  // speak the same section separator data structure).
  prev_section->receiveBackwardsSeparators(*inference_it,
					   previous_ps,
					   &(next_st->separatorCliques[next_ps.separatorCliquesSharedStructure.size()-1]),
					   next_ps.separatorCliquesSharedStructure[next_ps.separatorCliquesSharedStructure.size()-1]);
#endif
  if (IM::messageGlb(IM::InferenceMemory, IM::Med+9)) {
    // FIXME - previous_st->reportMemoryUsageTo(previous_ps,stdout);
  }

  if (! myjt->section_debug_range.contains((int)next_part_num)) {
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }
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
				      SparseJoinSectionTables &st,
				      vector<unsigned> const &roots,
				      vector< pair<unsigned,unsigned> > &message_order,
				      const char *const section_type_name,
				      const unsigned section_num,
				      const bool clear_when_done,
				      const bool also_clear_origins)
{
#if 1
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
	ceSendToOutgoingSeparators(ss.maxCliquesSharedStructure[from],
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
    // collect to section's root cliques
    set<unsigned> gathered_roots;
    for (unsigned i=0; i < roots.size(); ++i) {
      unsigned root = roots[i];
      if (gathered_roots.find(root) != gathered_roots.end()) continue;
      gathered_roots.insert(root);
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
  }
  if (! myjt->section_debug_range.contains((int)section_num)) {
    IM::setGlbMsgLevel(IM::InferenceMemory, inference_memory_debug_level);
    IM::setGlbMsgLevel(IM::Inference, inference_debug_level);
  }
  if (zeroClique) {
    throw ZeroCliqueException(); // continue to abort segment
  }
#else
  // FIXME - just run through the whole message order as below once,
  //         then do final gather from incoming seps for each root

  for (unsigned i=0; i < roots.size(); ++i) {
    ceGatherIntoRoot(ss,st,roots[i],message_order,section_type_name,section_num,clear_when_done,also_clear_origins);
  }
#endif
}
void 
SparseJoinInference::ceGatherIntoRoot(PartitionStructures &ss,
				      SparseJoinSectionTables &st,
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
	ceSendToOutgoingSeparators(ss.maxCliquesSharedStructure[from],
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
 * JunctionTree::deScatterOutofRoot
 *   
 *   Distribute Evidence Scatter Outof Root: This routine does a
 *   distribute evidence pass for this section, and scatters all
 *   messages outof the root within the current section. It does so
 *   using the message order given in the argument 'message_order',
 *   and scatters out of the provided root clique (which is the right
 *   interface clique of this section). By "Scatter", I mean it
 *   sends messages from the root clique distributing everything
 *   ultimately to all leaf cliques in this section.
 *
 * See Also:
 *   Dual routine: JunctionTree::ceGatherIntoRoot()
 *
 *
 * Preconditions:
 *   It is assumed that either:
 *     1) this is the ritht-most section
 *  or 2) that the right interface clique within this section has had a message
 *        sent to it from the right neighbor section.
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
SparseJoinInference::deScatterOutofRoot(// the section
					PartitionStructures &ss, 
					SparseJoinSectionTables &st,
					// root (right interface clique) of this section
					vector<unsigned> const &roots,
					vector< pair<unsigned,unsigned> > &message_order,
					// name (debugging/status msgs)
					const char *const sect_type_name,
					// section number (debugging/status msgs)
					const unsigned sect_num)
{
  // first check that this is not an empty section.
  if (ss.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! myjt->section_debug_range.contains((int)sect_num)) {
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }


  for (unsigned i=0; i < roots.size(); ++i) {
    unsigned root = roots[i];
    infoMsg(IM::Inference, IM::Med+5,"DE: distributing out of section root %s,section[%d]: clique %d\n",
	    sect_type_name, sect_num, root);

    // FIXME - the reverse message order must be such that root[i] will deScatterToOutgoingSeparators() 
    //         before any messages are passed out of root[i]. I.e., if roots = <a,b,c>, the reverse
    //         message order <(b,x), (c,y), (a,z)> would be invalid because the message out of root b
    //         would be sent before the st.maxCliques[b].deScatterToOutgoingSeparators(). Likewise for
    //         root c. 

    st.maxCliques[root].
      deScatterToOutgoingSeparators(ss.maxCliquesSharedStructure[root],
				    st.separatorCliques,
				    ss.separatorCliquesSharedStructure.ptr);
    for (unsigned msgNoP1=message_order.size();msgNoP1 > 0; msgNoP1 --) {
      const unsigned to = message_order[msgNoP1-1].first;
      const unsigned from = message_order[msgNoP1-1].second;
      infoMsg(IM::Inference, IM::Mod,"DE: message %s,section[%d]: clique %d <-- clique %d\n",
	      sect_type_name, sect_num, to, from);
      st.maxCliques[to].
	deReceiveFromIncommingSeparator(ss.maxCliquesSharedStructure[to],
					st.separatorCliques,
					ss.separatorCliquesSharedStructure.ptr);
      
      infoMsg(IM::Inference, IM::Med+5,"DE: distributing out of %s,section[%d]: clique %d\n",
	      sect_type_name, sect_num, to);
      st.maxCliques[to].
	deScatterToOutgoingSeparators(ss.maxCliquesSharedStructure[to],
				      st.separatorCliques,
				      ss.separatorCliquesSharedStructure.ptr);
    }
  }
  if (! myjt->section_debug_range.contains((int)sect_num)) {
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }


}

