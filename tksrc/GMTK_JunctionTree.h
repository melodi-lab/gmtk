/*
 * GMTK_JunctionTree.h
 *   GMTK Junction Tree. Exact and approximate inference support for GMTK.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * $Header$
 *
 */


#ifndef GMTK_JUNCTIONTREE_H
#define GMTK_JUNCTIONTREE_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <regex.h>

#include "bp_range.h"

#include "GMTK_RV.h"
#include "GMTK_FileParser.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_MaxClique.h"
#include "GMTK_JT_Partition.h"
#include "GMTK_PartitionStructures.h"
#include "GMTK_PartitionTables.h"

#include "debug.h"

// class mention for forward references.
class GraphicalModel;
class BoundaryTriangulate;
class Partition;
class GMTemplate;
class JunctionTree;

class JunctionTree {

  friend class GMTemplate;
  friend class BoundaryTriangulate;

  // The set of base partitions from which real unrolled things are cloned from.
  // When unrolling zero time, we get:
  //   u0: P' E'
  // When unrolling 1 or more times, the method depends on
  // if the template was created using either the left or right interface
  // method.
  // If template created using left interface method, we do:
  //  u0: P' E'
  //  u1: P' C' E' 
  //  u2: P' C' C' E' 
  //  u3: P' C' C' C' E'
  //  u4: etc.
  // Note that for left interface, an E1 contains M copies of 
  // the original chunk C, so u0 is the basic template.
  //
  // If template created using right interface method, we do:
  //  u0: P' E'
  //  u1: P' C' E' 
  //  u2: P' C' C' E' 
  //  u3: P' C' C' C' E'
  //  u4: etc.
  // in the right interface case, P' contains an original P and M
  // copies of C.  The next three variables hold partitions P', C',
  // and E', where, for the standard left interface and simple
  // boundary case, we have that P' = P, C' = C , and E' = [C
  // E]. These partitions use the *same* set of C++ random variables
  // (so the intersection of the random variables in P1 and Co will be
  // the interface).
  JT_Partition P1; 
  JT_Partition Co;   // C "other", depending on if right or left interface method is used.
  JT_Partition E1; 

  // The set of random variables corresponding to the union of the rvs
  // P1, Co, E1, corresponding to the template unrolled M+S-1
  // times. In other words, the number of C repetitions is M + S, so
  // we unroll one less than that, see
  // BoundaryTriangulate::findPartitions() for more information. These
  // are determined in create_base_partitions().
  vector <RV*> partition_unrolled_rvs; 
  // mapping from name(frame) to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > partition_ppf;


  // The names of the above three partitions to use for printing,
  // debugging, etc.
  static const char* P1_n;
  static const char* Co_n;
  static const char* E1_n;
  
  // Note, while we need extra separator cliques that are between the
  // corresponding partitions interface cliques, these separators will
  // live in the partition on the right of the separator. They will be
  // pointed to by the left interface clique in the partition on the
  // right.

  // The partition structures that hold structures for a set of
  // RVs unrolled enough to cover any actual length DGM.
  sArray <PartitionStructures> partitionStructureArray;
  // The partition tables that hold the actual clique/separator tables
  // for a DGM. This might be much longer than the
  // partitionStructureArray but is certainly no shorter.
  sArray <PartitionTables> partitionTableArray;
  // the evidence probability used during island algorithm.
  logpr cur_prob_evidence;
  // the EM training beam used for island training (TODO:, move this elsewhere, perhaps in clique)
  double curEMTrainingBeam;

  ////////////////////////////////////////////////////////////////////////
  // collect/distribute/probE support variables used in various routines.
  // current set of unrolled random variables
  vector <RV*> cur_unrolled_rvs;
  // current mapping from 'name+frame' to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > cur_ppf;
  ////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // various other parameters regarding the current segment. These
  // are all computed by JunctionTree::unroll(...)
  // 
  // 1) The number of frames given, but in units of the amount by which
  // the basic template would need to be unrolled to match the given
  // number of frames.
  unsigned basicTemplateMaxUnrollAmount;
  // 2) For the number of frames given, the min amount by which we can
  // unroll and still re-use the rvs accross time.
  unsigned basicTemplateMinUnrollAmount;
  // 3) how much the modified [P' C' E'] template would need to be
  // unrolled to match the number of frames given. This
  // value could be negative.
  int modifiedTemplateMaxUnrollAmount;
  // 4) how much the modified [P' C' E'] template would need to be
  // unrolled to match the number of frames given for the structure
  // case (i.e., this would either be [P' E'], [P' C' E'], or [P' C' C' E']
  // so this value is never more than 1 and could be negative.
  int modifiedTemplateMinUnrollAmount;
  // 5) How many usable frames are there in the current segment (for
  //    some boundaries/triangulations and their resulting partitions,
  //    it might not be possible to have a segment that is a multiple
  //    of one (1) and so this number is guaranteed to be this.
  unsigned numUsableFrames;
  // 6) The offset by which we should start the observations (should be given
  // to observation matrix). Why this exists, see 5).
  unsigned frameStart;
  // 7) The current total number of real frames.
  unsigned curNumFrames;
  ////////////////////////////////////////////////////////////////////////

  // An iterator to iterate partitions in increasing order, starting
  // from the left and moving to the right, both the
  // partitionStructureArray and the partitionTableArray
  // simultaneously. This iterator does *NOT* work like STL iterators.
  // Rather, it iterates between strictly the first and last entry but
  // does not move off the end (i.e., like it.end()), so one has to
  // check specifically for end. Also, to understand this class,
  // it would be useful to look at the computeUnrollParameters()
  // function in class GMTemplate.
  class ptps_iterator {
  private:
    // the current partition table (pt) index (there can be any number of tables).
    unsigned _pt_i;
    // the current partition structure (ps) index (there are at most 4 structures)
    unsigned _ps_i; // could remove this and compute it as needed from _pt_i and _pt_len
    
    // associated jt information.
    JunctionTree& jt;
    const unsigned _pt_len;
    
  public:

    // constructor and assignment routines ...

    // constructor for linear inference.
    ptps_iterator(JunctionTree& _jt) 
      : jt(_jt),_pt_len(jt.partitionTableArray.size()) 
    {
      set_to_first_entry();
    }

    // Constructor with explicit partition tables length, used e.g.,
    // with island so that we can still use this as an iterator but
    // without the partitionTableArray needing to be so long.
    ptps_iterator(JunctionTree& _jt,unsigned arg_pt_len) 
      : jt(_jt),_pt_len(arg_pt_len)
    {
      set_to_first_entry();
    }
    // copy constructor
    ptps_iterator(ptps_iterator& tmp) 
      : jt(tmp.jt),_pt_len(tmp._pt_len)
    {
      set_to_first_entry();
    }


    // lengths
    unsigned pt_len() {
      return _pt_len;
    }
    unsigned ps_len() {
      return std::min(_pt_len,4u);
    }

    // initialization routines.
    void set_to_first_entry() { 
      // set to be at E'
      _pt_i = _ps_i = 0; 
    }
    void set_to_last_entry() {
      // set to be at E'?
      _pt_i = pt_len()-1;
      _ps_i = ps_len()-1;
    }

    void go_to_part_no (unsigned i) {
      assert ( i < pt_len() );
      if (pt_len() < 5 || i < 3) {
	_pt_i = i; _ps_i = i;
      } else {
	if (i + 1 == pt_len())
	  set_to_last_entry();
	else {
	  _pt_i = i;
	  _ps_i = 2;
	}
      }
    }

    // assignment routines
    ptps_iterator& operator=(const ptps_iterator& rhs) {
      assert (&jt == &rhs.jt); // can't do cross JT objects.
      _pt_i = rhs._pt_i;
      _ps_i = rhs._ps_i;
      return *this;
    }
    bool operator==(const ptps_iterator& o) {
      return ((_pt_i == o._pt_i) && (_ps_i == o._ps_i));
    }
    bool operator!=(const ptps_iterator& other) {
      return !(*this == other);
    }


    // location identification routines.
    // The current indices
    unsigned pt_i() { return _pt_i; }
    unsigned ps_i() { return _ps_i; }
    // The previous indices, invalid to call these
    // if we are at begin state.
    unsigned pt_prev_i() { return _pt_i-1; }
    unsigned ps_prev_i() { return _ps_i-1; }
    // The next indices, invalid to call these
    // if we are at end state.
    unsigned pt_next_i() { return _pt_i+1; }
    unsigned ps_next_i() { return _ps_i+1; }

    bool at_first_entry() {
      return (_ps_i == 0) ;
    }
    bool at_second_entry() {
      return (_ps_i == 1) ;
    }
    bool at_penultimate_entry() {
      return (_pt_i + 2  == pt_len());
    }
    bool at_last_entry() {
      // are we at E'?
      return (_pt_i + 1  == pt_len());
    }
    unsigned num_c_partitions() {
      // return the number of C partitions in the table (rather than
      // structure) array.
      return (pt_len() - 2);
    }
    bool has_c_partition() { return (num_c_partitions() > 0); }
    bool at_p() { return at_first_entry(); }
    bool at_c() { return (!at_p() && !at_e()); }
    bool at_e() { return at_last_entry(); }
    bool at_last_c() { 
      // are we at last C'?
      return (has_c_partition() && at_penultimate_entry());
    }
    bool at_first_c() { 
      // are we at the first C'?
      return (has_c_partition() && at_second_entry());
    }

    bool next_at_p() { return false; }
    bool next_at_c() { 
      const unsigned dist_from_last_entry = jt.partitionTableArray.size() - _pt_i - 1;
      return (dist_from_last_entry > 1);
    }
    bool next_at_e() { 
      return at_penultimate_entry();
    }
    bool prev_at_p() {  return (_pt_i == 1); }
    bool prev_at_c() {  return (_pt_i > 1); }
    bool prev_at_e() { return false; }

    bool at_entry(unsigned pos) {
      return (pt_i() == pos);
    }

    // movement (increment/decrement) routines 

    // prefix
    ptps_iterator& operator ++() {
      if (pt_len() <= 4) {
	assert ( ps_len() == pt_len() );
	// either [P' E'], or [P' C' E'], or [P' C1' C2' E']
	if (!at_last_entry()) {
	  // advance both table and structure to next partition
	  // since they stay in sync for this pt_len.
	  _pt_i++;
	  _ps_i++;
	} else {
	  // do nothing, already a end, we don't move off the end like
	  // with a STL iterator.
	}
      } else {
	// pt_len() > 4
	// when we increment through the partitions, the index numbers take the form:
	// (eg., when pt_len() = 9 and ps_len() = 4)
	// 
	//    pt:   P C C C C C C C E
	//    pt_i: 0 1 2 3 4 5 6 7 8
        //    ps:   P C C E
	//    ps_i: 0 1 2 2 2 2 2 2 3
	// 
	// in other words, the ps index stops at table 1 only once and
	// then sticks at partition 2 until the final E'. Therefore,
	// in the case where the structures array has 4 entries (P C C
	// E), we see that ps_i indexes the right entry of if a pair
	// of partitions (where the pairs are 0P, PC, CC, CE) rather
	// than the left (which would correspond to PC, CC, CE, E0),
	// where '0' corresponds to a null non-existent partition. We
	// note that the decrement operator in this class mimics this
	// behavior in the reverse.
	assert ( ps_len() == 4 );
	if (_pt_i < (pt_len()-2)) {
	  // advance table to next partition
	  _pt_i ++;
	  // keep structure at last C' partition
	  _ps_i = std::min (_ps_i+1,2u);
	} else if (!at_last_entry()) {
	  // advance both, allowing to continue to reach last entry
	  _pt_i++;
	  _ps_i++;
	} else {
	  // do nothing, already at last entry.
	}
      }
      return *this;
    }
    // postfix
    ptps_iterator operator ++(int) {
      ptps_iterator tmp(*this);
      ++(*this);
      return tmp;
    }


    // prefix
    ptps_iterator& operator --() {
      if (pt_len() <= 4) {
	assert ( ps_len() == pt_len() );
	// either [P' E'], or [P' C' E'], or [P' C1' C2' E']
	if (!at_first_entry()) {
	  // decrease both table and structure to next partition
	  // since they stay in sync for this pt_len().
	  _pt_i--;
	  _ps_i--;
	} else {
	  // do nothing, already a start, we don't move off the beginning.
	}
      } else {
	// pt_len() > 4
	// when we descrement through the partitions, we want to be
	// compatible with the increment, i.e., we should be able
	// to increment and descrement through the array. Therefore,
	// the index numbers are compatible.
	// (eg., when pt_len() = 9 and ps_len = 4)
	// 
	//    pt:   E C C C C C C C P
	//    pt_i: 8 7 6 5 4 3 2 1 0
        //    ps:   E C C P
	//    ps_i: 3 2 2 2 2 2 2 1 0
	// 
	// in other words, the ps index stays at partition 2 until the
	// initial C' at which it takes entry 1 only once and then
	// finally decrements to starting P'. Therefore, in the case
	// where the structures array has 4 entries (P C C E), we see
	// that ps_i indexes the right entry of if a pair of
	// partitions (where the pairs are 0P, PC, CC, CE) rather than
	// the left (which would correspond to PC, CC, CE, E0), where
	// '0' corresponds to a null non-existent partition. We note
	// that the increment operator in this class mimics this
	// behavior in the reverse.
	assert ( ps_len() == 4 );
	if (_pt_i > 2) {
	  // advance table to next partition
	  _pt_i --;
	  // keep structure at last C' partition
	  _ps_i = 2;
	} else if (!at_first_entry()) {
	  // decrement both, allowing to continue to reach first entry
	  _pt_i--;
	  _ps_i--;
	} else {
	  // do nothing, already at first entry.
	}
      }
      return *this;
    }
    // postfix
    ptps_iterator operator --(int) {
      ptps_iterator tmp(*this);
      --(*this);
      return tmp;
    }

    // the current partition name
    const char* cur_nm() {
      if (at_first_entry())
	return jt.P1_n;
      else if (at_last_entry())
	return jt.E1_n;
      else 
	return jt.Co_n;
    }

    // the previous partition name
    const char* prev_nm() {
      if (at_first_c() || (at_e() && !has_c_partition()))
	return jt.P1_n;
      else if (!at_first_entry())
	return jt.Co_n;
      else
	return "LEFT_INVALID";
    }
    // the next partition name
    const char* next_nm() {
      if (at_last_c() || (at_p() && !has_c_partition()))
	return jt.E1_n;
      else if (!at_last_entry())
	return jt.Co_n;
      else
	return "RIGHT_INVALID";
    }

    // the current partition's root clique number 
    // (equivalently, its right interface clique)
    unsigned cur_ri() {
      // also see JunctionTree::computePartitionInterfaces()
      if (at_first_entry())
	return jt.P_ri_to_C;
      else if (at_last_entry())
	return jt.E_root_clique;
      else 
	return jt.C_ri_to_C;
    }
    // the previous partition's right interface clique number
    unsigned prev_ri() {
      if (at_p())
	return ~0x0; // invalid
      else if (prev_at_p())
	return jt.P_ri_to_C;
      else 
	return jt.C_ri_to_C;
    }
    unsigned next_ri() {
      if (at_e())
	return ~0x0; // invalid
      else if (next_at_e())
	return jt.E_root_clique;
      else 
	return jt.C_ri_to_C;
    }

    // the current partition's left interface clique number
    unsigned cur_li() {
      // also see JunctionTree::computePartitionInterfaces()
      if (at_first_entry())
	return 0; // invalid entry
      else if (at_last_entry())
	return jt.E_li_to_C;
      else 
	return jt.C_li_to_C;
    }

    unsigned prev_li() {
      if (at_p())
	return ~0x0; // invalid
      else if (prev_at_p())
	return ~0x0; // invalid
      else 
	return jt.C_li_to_C;
    }

    unsigned next_li() {
      if (at_e())
	return ~0x0; // invalid
      else if (next_at_e())
	return jt.E_li_to_C;
      else 
	return jt.C_li_to_C;
    }

    vector< pair<unsigned,unsigned> > & cur_message_order() {
      if (at_first_entry())
	return jt.P1_message_order;
      else if (at_last_entry())
	return jt.E1_message_order;
      else 
	return jt.Co_message_order;
    }

    JT_Partition& cur_jt_partition() {
      if (at_first_entry())
	return jt.P1;
      else if (at_last_entry())
	return jt.E1;
      else 
	return jt.Co;
    }

    BP_Range* cur_part_clique_print_range() {
      if (at_p())
	return jt.pPartCliquePrintRange;
      else if (at_c()) 
	return jt.cPartCliquePrintRange;
      else 
	return jt.ePartCliquePrintRange;
    }

    void printState(FILE* f) {

      fprintf(f,"partition iterator: pt_len=%d,ps_len=%d,pt_i=%d,pt_prev_i=%d,pt_next_i=%d,ps_i=%d,ps_prev_i=%d,ps_next_i=%d\n",
	      pt_len(),ps_len(),
	      pt_i(),pt_prev_i(),pt_next_i(),
	      ps_i(),ps_prev_i(),ps_next_i());

      fprintf(f," at_first_entry()=%d,at_second_entry()=%d,at_penultimate_entry()=%d,at_last_entry()=%d,num_c_partitions()=%d\n",
	     at_first_entry(),at_second_entry(),at_penultimate_entry(),at_last_entry(),num_c_partitions());

      fprintf(f," at_p()=%d,at_c()=%d,at_e()=%d,at_last_c()=%d,at_first_c()=%d,next_at_p()=%d,next_at_c()=%d,next_at_e()=%d,prev_at_p()=%d,prev_at_c()=%d,prev_at_e()=%d\n",
	     at_p(),at_c(),at_e(),at_last_c(),at_first_c(),next_at_p(),next_at_c(),next_at_e(),prev_at_p(),prev_at_c(),prev_at_e());
      fprintf(f," cur_li()=%d,prev_li()=%d,next_li()=%d, cur_ri()=%d,prev_ri()=%d,next_ri()=%d\n",
	      cur_li(),prev_li(),next_li(),
	      cur_ri(),prev_ri(),next_ri());
      fprintf(f," cur_nm()=%s,prev_nm()=%s,next_nm()=%s\n",
	     cur_nm(),prev_nm(),next_nm());
    }

  };


  // this is an iterator that is used by island (it obviously
  // is not always used in the typical sequential left/right fashion.
  ptps_iterator inference_it;
  // sets of random variables for frame shifting.
  set <RV*> cur_CC_rvs;
  unsigned cur_cc_shift;
  set <RV*> cur_CE_rvs;
  unsigned cur_ce_shift;
  void shiftCCtoPosition(int);
  void shiftCCrelative(int delta) { shiftCCtoPosition(cur_cc_shift+delta); }
  void shiftCEtoPosition(int);
  void shiftCErelative(int delta) { shiftCEtoPosition(cur_cc_shift+delta); }
  void init_CC_CE_rvs(ptps_iterator& ptps_it);
  void setCurrentInferenceShiftTo(int pos);



  ////////////////////////////////////////////////////////////////////////
  // Support variables specific to the island algorithm.
  // A map that stores *only* the island partition tables in the
  // island algorithm, as a function of the absolute part number -
  map < unsigned, PartitionTables*> islandsMap;
  sArray <PartitionTables*> islandPartitionTableArray;

  ////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  // Support variables specific to Viterbi and N-best decoding
  // 
  sArray < unsigned > P_partition_values;
  sArray < unsigned > C_partition_values;
  sArray < unsigned > E_partition_values;
  void recordPartitionViterbiValue(ptps_iterator& it);


  ////////////////////////////////////////////////////////////////////////


  // Identities of cliques in junction trees: 
  // for P, 
  //    P's right  interface to C (a root in a JT partition)
  unsigned P_ri_to_C; 
  //    The next one does not exist since we currently always do CE first.
  // unsigned P_li_clique; 
  // 
  // for C
  //    C's left interface to P
  unsigned C_li_to_P;
  //    C's left interface to C (same as C_li_to_P)
  unsigned C_li_to_C;
  //    C's right interface to C (a root in a JT partition)
  unsigned C_ri_to_C;
  //    C's right interface to E (a root in a JT partition) (same as C_ri_to_C)
  unsigned C_ri_to_E;
  // 
  // for E, 
  // E's left interface to C
  unsigned E_li_to_C;
  // root inside of E.
  unsigned E_root_clique;

  // Booleans telling if the interface cliques of the two partitions
  // are the same, meaning we don't need both and can drop one (to
  // save a bit of computation). These are currently computed but are
  // not yet used for anything.
  bool P_to_C_icliques_same;
  bool C_to_C_icliques_same;
  bool C_to_E_icliques_same;

  // Message passing orders for each partition.  Increasing index
  // order is 'collect evidence' phase from left to right in direction
  // of time, and decreasing order is 'distribute evidence' phase from
  // right to left in direction of time. Note that this assumes that
  // the overal root node in the JT is on the far right within E
  // (which might not be the best order).  NOTE: These are kept here
  // rather than in the partitions, since they are re-used for all
  // cloned partitions.
  vector< pair<unsigned,unsigned> > P1_message_order;
  vector< unsigned > P1_leaf_cliques;
  vector< pair<unsigned,unsigned> > Co_message_order;
  vector< unsigned > Co_leaf_cliques;
  vector< pair<unsigned,unsigned> > E1_message_order;  
  vector< unsigned > E1_leaf_cliques;

  // do a bit of setup for the upcomming inference round.
  void prepareForNextInferenceRound();

  // memory use reporting
  void reportMemoryUsageTo(FILE *f,unsigned whichPartitions = 0x7);

  // A version of unroll that starts with the gm_template and fills up
  // base partitions.
  void create_base_partitions();
  void insertFactorClique(FactorClique& factorClique,FactorInfo& factor);


  // Helper routines that are private (only called by other member
  // functions of this class).
  static void setUpMessagePassingOrderRecurse(JT_Partition& part,
					      const unsigned root,
					      vector< pair<unsigned,unsigned> >&order,
					      const unsigned excludeFromLeafCliques,
					      vector< unsigned>& leaf_cliques);
  static void assignRVToClique(const char *const partName,
			       JT_Partition&part,
			       const unsigned root,
			       const unsigned depth,
			       RV* rv,
			       unsigned& numberOfTimesAssigned,
			       set<RV*>& parSet,
			       const bool allParentsObserved,
			       multimap< vector<double>, unsigned >& scoreSet);

  static void createDirectedGraphOfCliquesRecurse(JT_Partition& part,
					   const unsigned root,
					   vector< bool >& visited);
  static void getCumulativeAssignedNodes(JT_Partition& part,
					 const unsigned root);
  static void getCumulativeUnassignedIteratedNodes(JT_Partition& part,const unsigned root);





  void ceGatherIntoRoot(PartitionStructures& ps,
			PartitionTables& pt,
			const unsigned root,
			vector< pair<unsigned,unsigned> >& message_order,
			const char*const part_type_name,
			const unsigned part_num,
			const bool clearWhenDone = false,
			const bool alsoClearOrigins = false);

  void ceSendForwardsCrossPartitions(PartitionStructures& previous_ps,
			     PartitionTables& previous_pt,
			     const unsigned previous_part_root,
			     const char*const previous_part_type_name,
			     const unsigned previous_part_num,
			     PartitionStructures& next_ps,
			     PartitionTables& next_pt,
			     const unsigned next_part_leaf,
			     const char*const next_part_type_name,
			     const unsigned next_part_num);

  void deScatterOutofRoot(PartitionStructures& ps,
			  PartitionTables& pt,
			  const unsigned root,
			  vector< pair<unsigned,unsigned> >& message_order,
			  const char*const part_type_name,
			  const unsigned part_num);

  void deSendBackwardsCrossPartitions(PartitionStructures& previous_ps,
				      PartitionTables& previous_pt,
				      const unsigned previous_part_root,
				      const char*const previous_part_type_name,
				      const unsigned previous_part_num,
				      // 
				      PartitionStructures& next_ps,
				      PartitionTables& next_pt,
				      const unsigned next_part_leaf,
				      const char*const next_part_type_name,
				      const unsigned next_part_num
				      );

  // Support routines for island algorithm inference.
  void deleteIsland(const unsigned part);
  void storeIsland(const unsigned part,PartitionTables *pt);
  PartitionTables* createPartition(const unsigned part);

  PartitionTables* retreiveIsland(const unsigned part);

  void ceGatherIntoRoot(const unsigned part,
			PartitionTables *pt);
  void deScatterOutofRoot(const unsigned part,
			  PartitionTables* pt);
  void ceSendForwardsCrossPartitions(const unsigned lpart,
				     PartitionTables *lpt,
				     PartitionTables *rpt);

  void deSendBackwardsCrossPartitions(const unsigned left_part,
				      PartitionTables *lpt,
				      PartitionTables *rpt);
  logpr probEvidenceRoot(const unsigned part,
			 PartitionTables* pt);
  logpr setRootToMaxCliqueValue(const unsigned part,
				PartitionTables* pt);				
  void emIncrementIsland(const unsigned part,
			 PartitionTables* pt,
			 const logpr probE, 
			 const bool localCliqueNormalization);
  void printAllCliques(const unsigned part,
		       PartitionTables* pt,
		       FILE* f,
		       const bool normalize,
		       const bool justPrintEntropy = false);
  void printAllCliqueProbabilties(const unsigned part,
				  PartitionTables* pt);


  void collectDistributeIslandRecurse(const unsigned start,
				      const unsigned end,
				      const unsigned base,
				      const unsigned linear_section_threshold,
				      const bool runEMalgorithm,
				      const bool runViterbiAlgorithm,
				      const bool localCliqueNormalization);
  void collectDistributeIslandBase(const unsigned start,
				   const unsigned end,
				   const bool runEMalgorithm,
				   const bool runViterbiAlgorithm,
				   const bool localCliqueNormalization);

public:


  // Set to true if the JT should create extra separators for any
  // virtual evidence (VE) that might be usefully exploitable
  // computationally in the clique. I.e., we treat the parents
  // of an immediate observed child, where the child is a deterministic 
  // function of the parents, a a separator over the parents
  // to be intersected as normal with all the other separators.
  static unsigned useVESeparators;
  enum VESeparatorType { VESEP_PC = 0x1, VESEP_PCG = 0x2 }; 
  // booleans to indicate where ve-seps should be used.
  static unsigned veSeparatorWhere;
  enum VESeparatorWhere { VESEP_WHERE_P = 0x1, VESEP_WHERE_C = 0x2, VESEP_WHERE_E = 0x4 }; 

  
  // When doing inference if any kind, this variable determines
  // if we should clear the clique and separator value cache
  // between segments/utterances. It might be beneficial, for
  // example, to retain the value cache around between segments/utterances
  // if for example, there are many such values that are common. If
  // not, on the other hand, setting this to true will cause
  // an increase in memory use on each segment.
  static bool perSegmentClearCliqueValueCache;

  // Set to true if the JT weight that we compute should be an upper
  // bound.  It is not guaranteed to be a tight upper bound, but is
  // guaranteed to at least be *an* upper bound.
  static bool jtWeightUpperBound;

  // Set to true if the JT weight scoring mechansim should be more
  // conservative, meaning it should not underestimate the charge of a
  // node in a clique as much. See code for details.
  static bool jtWeightMoreConservative;

  // The priority string for selecting the next edge when constructing
  // the junction tree. Default is in .cc file, and see .cc file for
  // what options are supported.
  static const char* junctionTreeMSTpriorityStr;

  // The priority string for selecting which clique of a partition
  // (from the set of valid ones) should be used as the partition
  // interface clique. See .cc file in routine findInterfaceClique()
  // for documentation.
  static const char* interfaceCliquePriorityStr;

  // Set to > 0.0 if the JT weight that we compute should heavily
  // penalize any unassigned iterated nodes. Penalty = factor
  // that gets multiplied by number of unassigned iterated.
  static float jtWeightPenalizeUnassignedIterated;
  
  // scaling factors (must be > 0 and <= 1) corresponding
  // to how much the separator nodes' charge gets scaled. The more pruning
  // we do during inference, the lower this should go.
  static float jtWeightSparseNodeSepScale;
  static float jtWeightDenseNodeSepScale;


  // When doing scoring (prob(evidence)), do we make compute the 'viterbi'
  // score (meaning the score of the most probable set of variables
  // or P(evidence,best_hidden)), or the full inference score,
  // namely \sum_hidden P(evidence,hidden)
  static bool viterbiScore;

  // range of cliques within each partition to print out when doing
  // CE/DE inference. If these are NULL, then we print nothing.
  BP_Range* pPartCliquePrintRange;
  BP_Range* cPartCliquePrintRange;
  BP_Range* ePartCliquePrintRange;

  // constructor
  JunctionTree(GMTemplate& arg_gm_template)
    : curEMTrainingBeam(-LZERO),
      inference_it(*this),
      fp(arg_gm_template.fp),
      gm_template(arg_gm_template)
  {
    pPartCliquePrintRange = cPartCliquePrintRange = ePartCliquePrintRange = NULL;
  }
  ~JunctionTree() {
    delete pPartCliquePrintRange;
    delete cPartCliquePrintRange;
    delete ePartCliquePrintRange;
    clearAfterUnroll();
  }

  // the fixed file parser for this model, for RV unrolling, etc.
  FileParser& fp;

  // The fixed gm_template for this model, contains the
  // pre-triangulated graph.
  GMTemplate& gm_template;

  struct Edge {
    unsigned clique1;
    unsigned clique2;
    vector <double> weights;
  };

  // for sorting edges in decreasing weight order.
  struct EdgeCompare {  
    // sort descending
    bool operator() (const Edge& a, 
		     const Edge& b) {
      return (a.weights) > (b.weights);
    }
  };


  // Call many of the routines below in the right order.
  void setUpDataStructures(const char* varPartitionAssignmentPrior,
			   const char *varCliqueAssignmentPrior);

  // create the three junction trees for the basic partitions.
  void createPartitionJunctionTrees(const string pStr = junctionTreeMSTpriorityStr) {
    infoMsg(IM::Giga,"Creating of P partition JT\n");
    createPartitionJunctionTree(gm_template.P,pStr);
    infoMsg(IM::Giga,"Creating of C partition JT\n");
    createPartitionJunctionTree(gm_template.C,pStr);
    infoMsg(IM::Giga,"Creating of E partition JT\n");
    createPartitionJunctionTree(gm_template.E,pStr);
    infoMsg(IM::Giga,"Done creating P,C,E partition JTs\n");
  }
  // create a junction tree within a partition.
  static void createPartitionJunctionTree(Partition& part, 
					  const string pStr = junctionTreeMSTpriorityStr);

  // routine to find the interface cliques of the partitions
  void computePartitionInterfaces();
  // routine to create the factors in the appropriate partitions
  void createFactorCliques();
  // routine to find the interface cliques of a partition
  void computePartitionInterface(JT_Partition& part1,
				 unsigned int& part1_ric,
				 JT_Partition& part2,
				 unsigned int& part2_lic,
				 bool& icliques_same);


  // root the JT
  void createDirectedGraphOfCliques();
  static void createDirectedGraphOfCliques(JT_Partition& part,
				    const unsigned root);



  // Assign probability giving random variables to cliques (i.e.,
  // these are assigned only to cliques such that the random variables
  // and *all* their parents live in the clique, plus some other
  // criterion in order to make message passing as efficient as
  // possible).
  void assignRVsToCliques(const char* varPartitionAssignmentPrior,
			  const char *varCliqueAssignmentPrior);
  static void assignRVsToCliques(const char *const partName,
				 JT_Partition&part,
				 const unsigned rootClique,
				 const char* varPartitionAssignmentPrior,
				 const char *varCliqueAssignmentPrior);



  void assignFactorsToCliques();
  void assignFactorsToCliques(JT_Partition& part);


  // For the three partitions, set up the different message passing
  // orders that are to be used. This basically just does a tree
  // traversal using the previously selected root.
  void setUpMessagePassingOrders();
  static void setUpMessagePassingOrder(JT_Partition& part,
				       const unsigned root,
				       vector< pair<unsigned,unsigned> >&order,
				       const unsigned excludeFromLeafCliques,
				       vector< unsigned>& leaf_cliques);

  // Separator creation, meaning create the seperator objects
  // both within and between partitions. Given two neighboring
  // partitions L and R, the separator between the interface
  // cliques in L and R is contained in R.
  static void createSeparators(JT_Partition& part,
			       vector< pair<unsigned,unsigned> >&order);
  void createSeparators();

  // create the virtual evidence separators
  static void createVESeparators(JT_Partition& part);


  // Separator iteration order and accumulated set intersection
  // creation for separator driven clique potential creation, and
  // also updates the seperators partial accumulator structure and
  // sets up cliques other variables.
  static void computeSeparatorIterationOrder(MaxClique& clique,
					     JT_Partition& part);
  static void computeSeparatorIterationOrders(JT_Partition& part);
  void computeSeparatorIterationOrders();

  // Computes the preceding iterated unassigned nodes and therein the
  // set of assigned nodes in each clique that should/shouldn't be
  // iterated.
  void getCumulativeUnassignedIteratedNodes();

  // compute the assignment order for nodes in this
  // partition's cliques relative to each clique's incomming separators, and while
  // doing so, also set the dispositions for each of the resulting
  // nodes in each clique.
  void sortCliqueAssignedNodesAndComputeDispositions(const char *varCliqueAssignmentPrior);
  void sortCliqueAssignedNodesAndComputeDispositions(JT_Partition& part,
						     const char *varCliqueAssignmentPrior);


  // return an upper bound on the weight of the junction tree in the
  // given partition, where the JT weight is defined as the cost of
  // doing collect evidence on this JT.
  static double junctionTreeWeight(JT_Partition& part,
				   const unsigned rootClique,
				   set<RV*>* lp_nodes,
				   set<RV*>* rp_nodes);

  // Given a set of maxcliques for a partition, and an interface for
  // this (can be left right, or any set including empty, the only
  // condition is that it must be covered by at least one of the
  // cliques), compute the junction tree for this set and return the
  // estimated JT cost. This is a static routine so can be called from
  // anywhere.
  static double junctionTreeWeight(vector<MaxClique>& cliques,
				   const set<RV*>& interfaceNodes,
				   set<RV*>* lp_nodes,
				   set<RV*>* rp_nodes);
				   
  // 
  // Print all information about the JT. Must
  // have had computeSeparatorIterationOrders() called
  // already.
  void printAllJTInfo(const char* fileName);
  void printAllJTInfo(FILE* f,JT_Partition& part,const unsigned root,
		      set <RV*>* lp_nodes,set <RV*>* rp_nodes);
  void printAllJTInfoCliques(FILE* f,JT_Partition& part,const unsigned root,const unsigned treeLevel,
			     set <RV*>* lp_nodes,set <RV*>* rp_nodes);
  void printMessageOrder(FILE *f,vector< pair<unsigned,unsigned> >& message_order);
  void printCurrentRVValues(FILE* f);

  // various forms of clique printing
  void setCliquePrintRanges(char *p,char*c,char*e);
  // this one is general.
  void printAllCliques(FILE* f,const bool normalize,
		       const bool justPrintEntropy);

  void printAllCliques(PartitionStructures& ps,
		       PartitionTables& pt,
		       const unsigned partNo,
		       const char *const nm,
		       BP_Range* rng,
		       FILE* f,
		       const bool normalize,
		       const bool justPrintEntropy = false);

  // 
  // Do some last-minute data structure setup to prepare for
  // unrolling to work (such as preliminary and pre work for
  // leaving STL, etc.)
  static void prepareForUnrolling(JT_Partition& part);
  void prepareForUnrolling();

  // Set up internal structures for unrolled network k>=0 times, where
  // k is the number of times C' is duplicated (so unroll by 0 means
  // the basic template-prime, unroll by one means two C's, etc.).
  // Unrolling only affects the non-STL data structures.
  // void unroll(unsigned k);

  // unroll for frames = numFrames 
  // Unrolling only affects the non-STL data structures.
  // Returns number of frames actually used, or 0 if invalid num frames.
  enum UnrollTableOptions { LongTable, ShortTable, ZeroTable, NoTouchTable };
  unsigned unroll(unsigned numFrames,
		  const UnrollTableOptions tableOption = LongTable,
		  unsigned *totalNumberPartitions = NULL);

  // do any cleanup after an unroll. Useful to be called by
  // destructors, etc.
  void clearAfterUnroll();

  // Perhaps make different unrolls for decoding, unroll for EM
  // training unroll for viterbi training, etc.
  // ...

  // basic collect evidence phase on basic structures.
  void collectEvidence();
  void distributeEvidence();
  // compute P(E), probability of the evidence, after collect evidence has been run.
  logpr probEvidence();

  // the same, but don't unroll to length of segment.
  logpr probEvidenceFixedUnroll(const unsigned numFrames,
				unsigned* numUsableFrames = NULL,
				bool limitTime=false,
				unsigned *numPartitionsDone = NULL,
				const bool noE=false);
  // simple call
  logpr probEvidence(const unsigned numFrames, unsigned& numUsableFrames) {
    return probEvidenceFixedUnroll(numFrames,&numUsableFrames,false,NULL,false);
  }
  // version that does unrolling, and const. memory, & stops after timer interupt occurs.
  // Just like probEvidence() but keeps track of how many messages between partitions
  // have been done, and returns (even if not complete) when timer expires.  
  static bool probEvidenceTimeExpired;
  logpr probEvidenceTime(const unsigned numFrames, unsigned& numUsableFrames, 
			 unsigned &numPartitionsDone, const bool noE = false) {
    return probEvidenceFixedUnroll(numFrames,&numUsableFrames,true,&numPartitionsDone,noE);
  }

  // return the island's idea of the current prob of evidence
  logpr curProbEvidenceIsland() { return cur_prob_evidence; }
  // set island's EM beam.
  void setCurEMTrainingBeam(const double b) { curEMTrainingBeam = b; }

  // EM training increment, for use with collectEvidence and distributeEvidence.
  void emIncrement(const logpr probEvidence,
		   const bool localCliqueNormalization = false,
		   const double emTrainingBeam = -LZERO);

  // print P(E) to stdout using all cliques. After a ce,de stage,
  // all values should be the same.
  void printProbEvidenceAccordingToAllCliques();

  // Set the root to the max clique value and return that value.
  logpr setRootToMaxCliqueValue();

  // Routine that calls collect/distribute evidence using the island
  // algorithm (i.e., log-space inference).
  logpr
  collectDistributeIsland(const unsigned numFrames,
			  unsigned& numUsableFrames,
			  const unsigned base,
			  const unsigned linear_section_threshold,
			  const bool runEMalgorithm = false,
			  const bool runViterbiAlgorithm = false,
			  const bool localCliqueNormalization = false);


  // void saveViterbiValuesIsland(oDataStreamFile& vfile);
  // void saveViterbiValuesLinear(oDataStreamFile& vfile);
  // void saveViterbiValuesIsland(FILE*);
  void printSavedPartitionViterbiValues(FILE*,
					bool printObserved,
					regex_t *preg,
					char* partRangeFilter);

  void printSavedViterbiValues(FILE*,
			       bool printObserved = false,
			       regex_t *preg = NULL,
			       bool reverseOrder = false,
			       unsigned maxTriggerVars = 0,
			       const char **triggerVars = NULL,
			       const char **triggerValSets = NULL);




  // actuall message routines.
  // void collectMessage(MaxClique& from,MaxClique& to);
  // void distributeMessage(MaxClique& from,MaxClique& to);

  // Returns a good approximation of the weight of (optionally) P, C,
  // and E in the variables pWeight, cWeight, and eWeight for
  // a given number of frame unrollings.
  double junctionTreeWeight(const bool includeP,
			    const bool includeC,
			    const bool includeE,
			    const unsigned numFrames,
			    double& pWeight,
			    double& cWeight,
			    double& eWeight);

  // Compute the 'junction tree weight' (roughly, the log10(cost of
  // doing inference)) for the set of cliques given in cliques. Note,
  // cliques *must* be a valid set of maxcliques of a junction tree --
  // if they are not, unexpected results are returned.
  double junctionTreeWeight(vector<MaxClique>& cliques);


  // used to clear out hash table memory between segments
  void clearCliqueSepValueCache(bool force=false) {
    P1.clearCliqueSepValueCache(force);
    Co.clearCliqueSepValueCache(force);
    E1.clearCliqueSepValueCache(force);
  }

  // access to the current set of nodes.
  inline vector <RV*>& curNodes() { return cur_unrolled_rvs; }


};

#endif

