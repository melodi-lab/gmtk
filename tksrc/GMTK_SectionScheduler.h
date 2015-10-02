/*
 * GMTK_SectionScheduler.h Root class for the inference algorithms at
 *   the time series level.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 * SectionScheduler subclasses handle inference over a series of sections. 
 *
 * The specific inference task subclasses of SectionScheduler are 
 * (X is observed, Q is hidden, T is segment length):
 *
 *  - ProbEvidenceTask       compute P(X_{0:T-1}), forward pass only, O(1) memory
 *  - ViterbiTask            compute argmax_{Q_{0:T-1}} P(Q_{0:T-1} | X_{0:T-1})
 *  - ForwardBackwardTask    compute P(Q_{0:T-1} | X_{0:T-1})
 *  - SmoothingTask          compute argmax_{Q_{t-\tau}} P(Q_{t-\tau} | X_{0:t})  
 *
 * The time-series algorithm subclasses of SectionScheduler are:
 *  - LinearSectionScheduler        Just iterate sequentially over sections
 *  - IslandSectionScheduler        Skip through the time series, leaving "islands" of partial results to save memory
 *  - ArchipelagosSectionScheduler  Parallel version of Island
 * 
 * The SectionScheduler subclasses multiply inherit the *Task classes representing the
 * inference tasks the time series level inference algorithms can perform. Not all 
 * SectionSchedulers support every inference task, e.g., IslandSectionScheduler
 * can't do online filtering/smoothing (OnlineSmoothingTask).
 *
 * The Sections are responsible for doing inference within themselves via whatever
 * SectionInferenceAlgorithm the Section subclass implements. The SectionScheduler
 * classes just pass messages between Sections via instances of InterfaceSeparator.
 *
 */


#ifndef GMTK_SECTIONSCHEDULER_H
#define GMTK_SECTIONSCHEDULER_H

#include <vector>
#include <map>

#include "logp.h"

// TODO: what's the difference between range and bp_range?
#include "range.h"
#include "bp_range.h"
#include "mArray.h"

#include "fileParser.h"

#include "GMTK_ObservationSource.h"
#include "GMTK_FileParser.h"

#include "GMTK_PartitionStructures.h"
#include "GMTK_PartitionTables.h"

//#include "GMTK_SectionInferenceAlgorithm.h"
class SectionInferenceAlgorithm;

class SectionScheduler {

 public:


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

    // the current section's root clique number 
    // (equivalently, its right interface clique)
    unsigned cur_ri() {
      // also see SectionScheduler::computeSectionInterfaces()
      if (at_first_entry())
	return jt.P_ri_to_C;
      else if (at_last_entry())
	return jt.E_root_clique;
      else 
	return jt.C_ri_to_C;
    }
    // the previous section's right interface clique number
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

    // the current section's left interface clique number
    unsigned cur_li() {
      // also see SectionScheduler::computeSectionInterfaces()
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

    JT_Partition& cur_jt_section() {
      if (at_first_entry())
	return jt.P1;
      else if (at_last_entry())
	return jt.E1;
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
      fprintf(f," cur_li()=%d,prev_li()=%d,next_li()=%d, cur_ri()=%d,prev_ri()=%d,next_ri()=%d\n",
	      cur_li(),prev_li(),next_li(),
	      cur_ri(),prev_ri(),next_ri());
      fprintf(f," cur_nm()=%s,prev_nm()=%s,next_nm()=%s\n",
	     cur_nm(),prev_nm(),next_nm());
    }

  };



  SectionScheduler(GMTemplate &gm_template, FileParser &fp, SectionInferenceAlgorithm *algorithm, ObservationSource *obs_source) :
    p_clique_print_range(NULL), c_clique_print_range(NULL), e_clique_print_range(NULL),
    section_debug_range("all", 0, 0x7FFFFFFF), obs_source(obs_source), fp(fp), gm_template(gm_template), algorithm(algorithm), inference_it(*this)
  {}
    
  virtual ~SectionScheduler() {}


  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
  virtual void setUpDataStructures(iDataStreamFile &tri_file,
				   char const *varSectionAssignmentPrior,
				   char const *varCliqueAssignmentPrior,
				   bool checkTriFileCards) = 0;


  // Prepare to do inference on segment of length T. Note that this
  // is the analog of the current JunctionTree::unroll() which sets up the
  // section_structure_array for at most  P' C' C' E' (4 sections, so
  // O(1) memory), and the section_table_array for between 0 and T
  // sections depending on ZeroTable, ShortTable, or LongTable. 
  // JT::unroll() also allocates O(T) memory (optionally memory mapped)
  // to store the the Viterbi values if they aren't being written to
  // a file.

  // Since we now have the ability to write Viterbi values to files, do
  // we want to drop the ability to store them in memory? This would
  // break things, as a Viterbi file would then be a required argument.
  // It would also be slower for applications that do fit in memory
  // because of the higher overhead of file I/O operations.

  // Maybe rename this to prepareForSegment() or something to avoid confusion
  // with O(T) space graph unrolling?
  // returns number of usable frames
  enum UnrollTableOptions { LongTable, ShortTable, ZeroTable, NoTouchTable };
  virtual unsigned unroll(unsigned numFrames,
			  const UnrollTableOptions tableOption = LongTable,
			  unsigned *totalNumberSections = NULL); 


  // Formerly JunctionTree::printAllJTInfo()
  virtual void printInferencePlanSummary(char const *fileName) = 0;

  // Formerly GMTemplate::reportScoreStats()
  virtual void reportScoreStats() = 0;


  // NOTE:  This assumes all inference algorithms will have clique-like things they
  //        need to print posteriors of

  // Set the range of selected clique #'s in P', C', E' for printing.
  // TODO: preconditions
  virtual void setCliquePrintRanges(char *p_range, char *c_range, char *e_range) = 0;

  // Print to f the order of the variables in each clique selected by setCliquePrintRanges().
  virtual void printCliqueOrders(FILE *f) = 0;

  // Returns the size (in # of floats) of the cliques selected by setCliquePrintRanges().
  virtual void getCliquePosteriorSize(unsigned &p_size, unsigned &c_size, unsigned &e_size) = 0;

  // Set the range of section #s that get elevated verbosity
  // P' = 0, C' \in 1, ..., T=2, E' = T-1
  virtual void setSectionDebugRange(Range const &rng) {
    section_debug_range.SetLimits(rng.first(), rng.last()); 
    section_debug_range.SetDefStr(rng.GetDefStr()); 
  }

 protected:

  // The set of base sections from which real unrolled things are cloned from.
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
  // copies of C.  The next three variables hold sections P', C',
  // and E', where, for the standard left interface and simple
  // boundary case, we have that P' = P, C' = C , and E' = [C
  // E]. These sections use the *same* set of C++ random variables
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
  vector <RV*> section_unrolled_rvs; 
  // mapping from name(frame) to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > section_ppf;


  // The names of the above three sections to use for printing,
  // debugging, etc.
  static const char* P1_n;
  static const char* Co_n;
  static const char* E1_n;
  


  // Identities of cliques in junction trees: 
  // for P, 
  //    P's right  interface to C (a root in a JT section)
  unsigned P_ri_to_C; 
  //    The next one does not exist since we currently always do CE first.
  // unsigned P_li_clique; 
  // 
  // for C
  //    C's left interface to P
  unsigned C_li_to_P;
  //    C's left interface to C (same as C_li_to_P)
  unsigned C_li_to_C;
  //    C's right interface to C (a root in a JT section)
  unsigned C_ri_to_C;
  //    C's right interface to E (a root in a JT section) (same as C_ri_to_C)
  unsigned C_ri_to_E;
  // 
  // for E, 
  // E's left interface to C
  unsigned E_li_to_C;
  // root inside of E.
  unsigned E_root_clique;



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



  // Note, while we need extra separator cliques that are between the
  // corresponding sections interface cliques, these separators will
  // live in the section on the right of the separator. They will be
  // pointed to by the left interface clique in the section on the
  // right.

  // The section structures that hold structures for a set of
  // RVs unrolled enough to cover any actual length DGM.
  sArray <PartitionStructures> section_structure_array;
  // The section tables that hold the actual clique/separator tables
  // for a DGM. This might be much longer than the
  // section_structure_array but is certainly no shorter.
  sArray <PartitionTables> section_table_array;


  // range of cliques within each section to print out when doing
  // CE/DE inference. If these are NULL, then we print nothing.
  BP_Range* p_clique_print_range;
  BP_Range* c_clique_print_range;
  BP_Range* e_clique_print_range;

  Range section_debug_range;

  ObservationSource *obs_source; // & other common members?
  // inference task methods can dynamic cast to FileSource or StreamSource as needed?

  // the fixed file parser for this model, for RV unrolling, etc.
  FileParser   &fp;

  // The fixed gm_template for this model, contains the pre-triangulated graph.
  GMTemplate   &gm_template;


  // This guy does the actual work of inference within a section
  SectionInferenceAlgorithm *algorithm;


  // FIXME - does inference_it really need to be a member? 
  //   Maybe each *Task method has local SectionIterator(s) for its work...
  //   {Island,Archipelagos}SectionScheduler might want it as a member to
  //     avoid having to pass it as a parameter everywhere, but it could
  //     be in the subclass(es).
  //   What's the relative performance cost of passing it vs. calling this->inference_it.method() ?
  SectionIterator inference_it;
 
};

#endif
