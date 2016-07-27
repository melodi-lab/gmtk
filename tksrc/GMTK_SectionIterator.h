/*
 * GMTK_SectionIterator.h Root class for the inference algorithms at
 *   the time series level.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef SECTIONITERATOR_H
#define SECTIONITERATOR_H

#include "GMTK_SectionScheduler.h"

  // An iterator to iterate sections in increasing order, starting
  // from the left and moving to the right, both the
  // section_structure_array and the section_table_array
  // simultaneously. This iterator does *NOT* work like STL iterators.
  // Rather, it iterates between strictly the first and last entry but
  // does not move off the end (i.e., like it.end()), so one has to
  // check specifically for end. Also, to understand this class,
  // it would be useful to look at the computeUnrollParameters()
  // function in class GMTemplate.
  class SectionIterator {
   protected:
    // the current section table (st) index (there can be any number of tables).
    unsigned st_i;
    // the current section structure (ss) index (there are at most 4 structures)
    unsigned ss_i; // could remove this and compute it as needed from st_i and st_len
    
    // associated junction tree information.
    SectionScheduler& jt;
    unsigned st_len;
    
  public:

    // constructor and assignment routines ...

    // constructor for linear inference.
    SectionIterator(SectionScheduler &ss) 
      : jt(ss),st_len(jt.section_table_array.size()) 
    {
      set_to_first_entry();
    }

    // Constructor with explicit section tables length, used e.g.,
    // with island so that we can still use this as an iterator but
    // without the section_table_array needing to be so long.
    SectionIterator(SectionScheduler &ss, unsigned st_len) 
      : jt(ss),st_len(st_len)
    {
      set_to_first_entry();
    }
    // copy constructor
    SectionIterator(SectionIterator& tmp) 
      : jt(tmp.jt),st_len(tmp.st_len)
    {
      set_to_first_entry();
    }


    // lengths
    unsigned get_st_len() {
      return st_len;
    }
    unsigned get_ss_len() {
      return std::min(st_len,4u);
    }


    void set_st_len(unsigned st_len) {
      this->st_len = st_len;
    }

    // initialization routines.
    void set_to_first_entry() { 
      // set to be at E'
      st_i = ss_i = 0; 
    }
    void set_to_last_entry() {
      // set to be at E'?
      st_i = get_st_len()-1;
      ss_i = get_ss_len()-1;
    }

    void go_to_section_no (unsigned i) {
      assert (get_st_len()==0 || i < get_st_len() );
      if (get_st_len() < 5 || i < 3) {
	st_i = i; ss_i = i;
      } else {
	if (i + 1 == get_st_len())
	  set_to_last_entry();
	else {
	  st_i = i;
	  ss_i = 2;
	}
      }
    }

    // assignment routines
    SectionIterator& operator=(const SectionIterator& rhs) {
      assert (&jt == &rhs.jt); // can't do cross JT objects.
      st_i = rhs.st_i;
      ss_i = rhs.ss_i;
      return *this;
    }
    bool operator==(const SectionIterator& o) {
      return ((st_i == o.st_i) && (ss_i == o.ss_i));
    }
    bool operator!=(const SectionIterator& other) {
      return !(*this == other);
    }


    // location identification routines.
    // The current indices
    unsigned cur_st() { return st_i; }
    unsigned cur_ss() { return ss_i; }
    // The previous indices, invalid to call these
    // if we are at begin state.
    unsigned prev_st() { return st_i-1; }
    unsigned prev_ss() { return ss_i-1; }
    // The next indices, invalid to call these
    // if we are at end state.
    unsigned next_st() { return st_i+1; }
    unsigned next_ss() { return ss_i+1; }

    bool at_first_entry() {
      return (ss_i == 0) ;
    }
    bool at_second_entry() {
      return (ss_i == 1) ;
    }
    bool at_penultimate_entry() {
      return (st_i + 2  == get_st_len());
    }
    bool at_last_entry() {
      // are we at E'?
      return (st_i + 1  == get_st_len());
    }
    unsigned num_c_sections() {
      // return the number of C sections in the table (rather than
      // structure) array.
      return (get_st_len() - 2);
    }
    bool has_c_section() { return (num_c_sections() > 0); }
    bool at_p() { return at_first_entry(); }
    bool at_c() { return (!at_p() && !at_e()); }
    bool at_e() { return at_last_entry(); }
    bool at_last_c() { 
      // are we at last C'?
      return (has_c_section() && at_penultimate_entry());
    }
    bool at_first_c() { 
      // are we at the first C'?
      return (has_c_section() && at_second_entry());
    }

    bool next_at_p() { return false; }
    bool next_at_c() { 
      const unsigned dist_from_last_entry = jt.section_table_array.size() - st_i - 1;
      return (dist_from_last_entry > 1);
    }
    bool next_at_e() { 
      return at_penultimate_entry();
    }
    bool prev_at_p() {  return (st_i == 1); }
    bool prev_at_c() {  return (st_i > 1); }
    bool prev_at_e() { return false; }

    bool at_entry(unsigned pos) {
      return (cur_st() == pos);
    }

    // movement (increment/decrement) routines 

    // prefix
    SectionIterator& operator ++() {
      if (get_st_len() <= 4) {
	assert ( get_ss_len() == get_st_len() );
	// either [P' E'], or [P' C' E'], or [P' C1' C2' E']
	if (!at_last_entry()) {
	  // advance both table and structure to next section
	  // since they stay in sync for this st_len.
	  st_i++;
	  ss_i++;
	} else {
	  // do nothing, already a end, we don't move off the end like
	  // with a STL iterator.
	}
      } else {
	// get_st_len() > 4
	// when we increment through the sections, the index numbers take the form:
	// (eg., when get_st_len() = 9 and get_ss_len() = 4)
	// 
	//    pt:   P C C C C C C C E
	//    pt_i: 0 1 2 3 4 5 6 7 8
        //    ps:   P C C E
	//    ps_i: 0 1 2 2 2 2 2 2 3
	// 
	// in other words, the ps index stops at table 1 only once and
	// then sticks at section 2 until the final E'. Therefore,
	// in the case where the structures array has 4 entries (P C C
	// E), we see that ps_i indexes the right entry of if a pair
	// of sections (where the pairs are 0P, PC, CC, CE) rather
	// than the left (which would correspond to PC, CC, CE, E0),
	// where '0' corresponds to a null non-existent section. We
	// note that the decrement operator in this class mimics this
	// behavior in the reverse.
	assert ( get_ss_len() == 4 );
	if (st_i < (get_st_len()-2)) {
	  // advance table to next section
	  st_i ++;
	  // keep structure at last C' section
	  ss_i = std::min (ss_i+1,2u);
	} else if (!at_last_entry()) {
	  // advance both, allowing to continue to reach last entry
	  st_i++;
	  ss_i++;
	} else {
	  // do nothing, already at last entry.
	}
      }
      return *this;
    }
    // postfix
    SectionIterator operator ++(int) {
      SectionIterator tmp(*this);
      ++(*this);
      return tmp;
    }


    // prefix
    SectionIterator& operator --() {
      if (get_st_len() <= 4) {
	assert ( get_ss_len() == get_st_len() );
	// either [P' E'], or [P' C' E'], or [P' C1' C2' E']
	if (!at_first_entry()) {
	  // decrease both table and structure to next section
	  // since they stay in sync for this get_st_len().
	  st_i--;
	  ss_i--;
	} else {
	  // do nothing, already a start, we don't move off the beginning.
	}
      } else {
	// get_st_len() > 4
	// when we descrement through the sections, we want to be
	// compatible with the increment, i.e., we should be able
	// to increment and descrement through the array. Therefore,
	// the index numbers are compatible.
	// (eg., when get_st_len() = 9 and ps_len = 4)
	// 
	//    pt:   E C C C C C C C P
	//    pt_i: 8 7 6 5 4 3 2 1 0
        //    ps:   E C C P
	//    ps_i: 3 2 2 2 2 2 2 1 0
	// 
	// in other words, the ps index stays at section 2 until the
	// initial C' at which it takes entry 1 only once and then
	// finally decrements to starting P'. Therefore, in the case
	// where the structures array has 4 entries (P C C E), we see
	// that ps_i indexes the right entry of if a pair of
	// sections (where the pairs are 0P, PC, CC, CE) rather than
	// the left (which would correspond to PC, CC, CE, E0), where
	// '0' corresponds to a null non-existent section. We note
	// that the increment operator in this class mimics this
	// behavior in the reverse.
	assert ( get_ss_len() == 4 );
	if (st_i > 2) {
	  // advance table to next section
	  st_i --;
	  // keep structure at last C' section
	  ss_i = 2;
	} else if (!at_first_entry()) {
	  // decrement both, allowing to continue to reach first entry
	  st_i--;
	  ss_i--;
	} else {
	  // do nothing, already at first entry.
	}
      }
      return *this;
    }
    // postfix
    SectionIterator operator --(int) {
      SectionIterator tmp(*this);
      --(*this);
      return tmp;
    }

    // the current section name
    const char* cur_nm() {
      if (at_first_entry())
	return jt.P1_n;
      else if (at_last_entry())
	return jt.E1_n;
      else 
	return jt.Co_n;
    }

    // the previous section name
    const char* prev_nm() {
      if (at_first_c() || (at_e() && !has_c_section()))
	return jt.P1_n;
      else if (!at_first_entry())
	return jt.Co_n;
      else
	return "LEFT_INVALID";
    }
    // the next section name
    const char* next_nm() {
      if (at_last_c() || (at_p() && !has_c_section()))
	return jt.E1_n;
      else if (!at_last_entry())
	return jt.Co_n;
      else
	return "RIGHT_INVALID";
    }


    // FIXME: These (interface index & size accessors) need to move to 
    //        someplace "inference acrchitecture"-specific


    // the current section's root clique number 
    // (equivalently, its right interface clique)
    vector<unsigned> cur_ri() {
      // also see SectionScheduler::computeSectionInterfaces()
      if (at_first_entry())
	return jt.P_ri_to_C;
      else if (at_last_entry())
	return jt.E_root_clique;
      else 
	return jt.C_ri_to_C;
    }

    vector<unsigned> cur_roots() {
      if (at_first_entry())
	return jt.P1.subtree_roots;
      else if (at_last_entry())
	return jt.E1.subtree_roots;
      else 
	return jt.Co.subtree_roots;
    }

    // the previous section's right interface clique number
    vector<unsigned> prev_ri() {
      if (at_p())
	return vector<unsigned> (1,~0x0); // invalid
      else if (prev_at_p())
	return jt.P_ri_to_C;
      else 
	return jt.C_ri_to_C;
    }

    unsigned prev_ri_size() {
      if (at_p())
	return ~0x0; // invalid
      else if (prev_at_p())
	return jt.P_ri_to_C_size;
      else 
	return jt.C_ri_to_C_size;
    }

    vector<unsigned> next_ri() {
      if (at_e())
	return vector<unsigned> (1,~0x0); // invalid
      else if (next_at_e())
	return jt.E_root_clique;
      else 
	return jt.C_ri_to_C;
    }

    // the current section's left interface clique number
    vector<unsigned> cur_li() {
      // also see SectionScheduler::computeSectionInterfaces()
      if (at_first_entry())
	return vector<unsigned> (1,0); // invalid entry
      else if (at_last_entry())
	return jt.E_li_to_C;
      else 
	return jt.C_li_to_C;
    }

    vector<unsigned> prev_li() {
      if (at_p())
	return vector<unsigned> (1,~0x0); // invalid
      else if (prev_at_p())
	return vector<unsigned> (1,~0x0); // invalid
      else 
	return jt.C_li_to_C;
    }

    vector<unsigned> next_li() {
      if (at_e())
	return vector<unsigned> (1,~0x0); // invalid
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

    vector<unsigned> & cur_collect_clique_order() {
      if (at_first_entry())
	return jt.P1_collect_clique_order;
      else if (at_last_entry())
	return jt.E1_collect_clique_order;
      else 
	return jt.Co_collect_clique_order;
    }

    vector< pair<unsigned,unsigned> > & cur_reverse_msg_order() {
assert(jt.E1_reverse_messages.size()==0);
      if (at_first_entry())
	return jt.P1_reverse_messages;
      else if (at_last_entry())
	return jt.E1_reverse_messages;
      else 
	return jt.Co_reverse_messages;
    }

    JT_Partition& cur_jt_section() {
      if (at_first_entry())
	return jt.P1;
      else if (at_last_entry())
	return jt.E1;
      else 
	return jt.Co;
    }

    JT_Partition &prev_jt_section() {
      if (at_p())
	return jt.P1; // invalid
      else if (prev_at_p())
	return jt.P1;
      else 
	return jt.Co;
    }

    BP_Range* cur_section_clique_print_range() {
      if (at_p())
	return jt.p_clique_print_range;
      else if (at_c()) 
	return jt.c_clique_print_range;
      else 
	return jt.e_clique_print_range;
    }

    void printState(FILE* f) {

      fprintf(f,"section iterator: st_len=%d,ps_len=%d,pt_i=%d,prev_st=%d,next_st=%d,ps_i=%d,prev_ss()=%d,next_ss()=%d\n",
	      get_st_len(),get_ss_len(),
	      cur_st(),prev_st(),next_st(),
	      cur_ss(),prev_ss(),next_ss());

      fprintf(f," at_first_entry()=%d,at_second_entry()=%d,at_penultimate_entry()=%d,at_last_entry()=%d,num_c_sections()=%d\n",
	     at_first_entry(),at_second_entry(),at_penultimate_entry(),at_last_entry(),num_c_sections());

      fprintf(f," at_p()=%d,at_c()=%d,at_e()=%d,at_last_c()=%d,at_first_c()=%d,next_at_p()=%d,next_at_c()=%d,next_at_e()=%d,prev_at_p()=%d,prev_at_c()=%d,prev_at_e()=%d\n",
	     at_p(),at_c(),at_e(),at_last_c(),at_first_c(),next_at_p(),next_at_c(),next_at_e(),prev_at_p(),prev_at_c(),prev_at_e());
      // FIXME - print all interface nodes, not just the first
      fprintf(f," cur_li()=%d,prev_li()=%d,next_li()=%d, cur_ri()=%d,prev_ri()=%d,next_ri()=%d\n",
	      cur_li()[0],prev_li()[0],next_li()[0],
	      cur_ri()[0],prev_ri()[0],next_ri()[0]);
      fprintf(f," cur_nm()=%s,prev_nm()=%s,next_nm()=%s\n",
	     cur_nm(),prev_nm(),next_nm());
    }

  };

#endif
