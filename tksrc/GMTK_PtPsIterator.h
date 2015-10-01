/*
 * GMTK_PtPsIterator.h
 *   Iterator for SectionTables and SectionStructures (formerly 
 *   Partiton, hence the Ps)
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 *
 * $Header$
 *
 */


#ifndef GMTK_PTPSITERATOR_H
#define GMTK_PTPSITERATOR_H

#include <stdio.h>
#include "GMTK_JunctionTree.h"
#include "debug.h"


  // An iterator to iterate partitions in increasing order, starting
  // from the left and moving to the right, both the
  // partitionStructureArray and the partitionTableArray
  // simultaneously. This iterator does *NOT* work like STL iterators.
  // Rather, it iterates between strictly the first and last entry but
  // does not move off the end (i.e., like it.end()), so one has to
  // check specifically for end. Also, to understand this class,
  // it would be useful to look at the computeUnrollParameters()
  // function in class GMTemplate.
  class PtPsIterator {
  protected:
    // the current partition table (pt) index (there can be any number of tables).
    unsigned _pt_i;
    // the current partition structure (ps) index (there are at most 4 structures)
    unsigned _ps_i; // could remove this and compute it as needed from _pt_i and _pt_len
    
    // associated jt information.
    JunctionTree& jt;
    unsigned _pt_len;
    
  public:

    // constructor and assignment routines ...

    // constructor for linear inference.
    PtPsIterator(JunctionTree& _jt) 
      : jt(_jt),_pt_len(jt.partitionTableArray.size()) 
    {
      set_to_first_entry();
    }

    // Constructor with explicit partition tables length, used e.g.,
    // with island so that we can still use this as an iterator but
    // without the partitionTableArray needing to be so long.
    PtPsIterator(JunctionTree& _jt,unsigned arg_pt_len) 
      : jt(_jt),_pt_len(arg_pt_len)
    {
      set_to_first_entry();
    }
    // copy constructor
    PtPsIterator(PtPsIterator& tmp) 
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


    void set_pt_len(unsigned pt_len) {
      _pt_len = pt_len;
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
      assert (pt_len()==0 || i < pt_len() );
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
    PtPsIterator& operator=(const PtPsIterator& rhs) {
      assert (&jt == &rhs.jt); // can't do cross JT objects.
      _pt_i = rhs._pt_i;
      _ps_i = rhs._ps_i;
      return *this;
    }
    bool operator==(const PtPsIterator& o) {
      return ((_pt_i == o._pt_i) && (_ps_i == o._ps_i));
    }
    bool operator!=(const PtPsIterator& other) {
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
    PtPsIterator& operator ++() {
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
    PtPsIterator operator ++(int) {
      PtPsIterator tmp(*this);
      ++(*this);
      return tmp;
    }


    // prefix
    PtPsIterator& operator --() {
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
    PtPsIterator operator --(int) {
      PtPsIterator tmp(*this);
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

#endif

