/*-
 * GMTK_CFunctionDecisionTrees.cc
 *
 *     Code that defines and declares internal pre-defined GMTK deterministic
 *     mapping functions, that map from a set of random variables (nominally a
 *     set of parent variables and one child variable) that are non-negative integer
 *     valued down to a single integer value.
 *     This file serves two purposes:
 *         1) to define and register internal oft-used determinstic
 *            mapping functions that are useful in a variety of
 *            contexts (such as copy parent, etc.).  
 *         2) to allow the user to define their own mapping function. This is useful
 *            when a decision tree is large, complicated, and most importantly slow.
 *            By being able to define it here, at compile time, we can take advantage
 *            of the optimizing C++ comiler producing an efficient implementation of
 *            the desired formula (not to mention that now, not only formulas but
 *            one can also use loops, subroutine calls, local variables, floating
 *            point, and so on).  
 *            In this case, the user must make sure that their function doesn't have
 *            infinite loops, and so on, as otherwise inference will just stall.
 *    How to use this file?
 *         1) Search for the tag REGISTRATION_OF_USER_DEFINED_DETERMINISTIC_MAPPERS:
 *            Here you define a bunch of C-like functions that become deterministic mappers.
 *            You can define any C subroutines in this region of the code as well.
 *         2) Next, search for the tag REGISTRATION_OF_USER_DEFINED_DETERMINISTIC_MAPPERS:
 *            below and copy/uncomment necessary code to register your newly defined 
 *            deterministic mappers.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2010, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "error.h"
#include "general.h"
#include "rand.h"
#include "sArray.h"

#define DEFINE_DETERMINISTIC_MAPPER_MACROS 1
#include "GMTK_CFunctionDeterministicMappings.h"
#include "GMTK_RngDecisionTree.h"
#include "GMTK_DiscRV.h"
#include "GMTK_GMParms.h"

//DRIP constants necessary for dt mappings
#include "drip_constants.inc"


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// GMTK Internal C function deterministic mapping functions.
// DO NOT MODIFY ANYTHING IN THE NEXT BIT OF CODE STARTING HERE.
//
#define COPYPARENT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(copyParent,COPYPARENT_NUM_FEATURES)
{
  DiscRVType rv = p0; 
  return rv;
}
//
#define ALWAYSZERO_NUM_FEATURES 0
DEFINE_DETERMINISTIC_MAPPER_C_CODE(alwaysZero,ALWAYSZERO_NUM_FEATURES)
{
  return (DiscRVType) 0;
}
//
#define ALWAYSONE_NUM_FEATURES 0
DEFINE_DETERMINISTIC_MAPPER_C_CODE(alwaysOne,ALWAYSONE_NUM_FEATURES)
{
  return (DiscRVType) 1;
}
//
#define INCREMENT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(increment,INCREMENT_NUM_FEATURES)
{
  return (p0+1);
}
//
#define DECREMENT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(decrement,DECREMENT_NUM_FEATURES)
{
  return (p0-1);
}
//
#define CONDITIONAL_INCREMENT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(conditionalIncrement,CONDITIONAL_INCREMENT_NUM_FEATURES)
{
  // increments p0 if p1 is non-zero, otherwise returns p0.
  if (p1)
    return (p0+1);
  else
    return (p0);
}
//
#define CONDITIONAL_DECREMENT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(conditionalDecrement,CONDITIONAL_DECREMENT_NUM_FEATURES)
{
  // decrements p0 if p1 is non-zero, otherwise returns p0.
  if (p1)
    return (p0-1);
  else
    return (p0);
}
//
#define CONDITIONAL_LIMITED_INCREMENT_NUM_FEATURES 3
DEFINE_DETERMINISTIC_MAPPER_C_CODE(conditionalLimitedIncrement,CONDITIONAL_LIMITED_INCREMENT_NUM_FEATURES)
{
  // increments p0 if p1 is non-zero. Increment up to and including value given by p2
  // but not beyond. Otherwise returns p0.
  if (p1 && p0 < p2)
    return (p0+1);
  else
    return (p0);
}
//
#define CONDITIONAL_LIMITED_DECREMENT_NUM_FEATURES 3
DEFINE_DETERMINISTIC_MAPPER_C_CODE(conditionalLimitedDecrement,CONDITIONAL_LIMITED_DECREMENT_NUM_FEATURES)
{
  // decrements p0 if p1 is non-zero. Decrement down to and including value given by p2
  // but not below. Otherwise returns p0.
  if (p1 && p0 > p2)
    return (p0-1);
  else
    return (p0);
}


//
// The below includes routines with a variable number of parents.
// TODO: get the below working, with variable numbers of parents, we
// need to change DT code to allow variable num features.
//
DEFINE_DETERMINISTIC_MAPPER_C_CODE(allParentsEqual,CDT_VARIABLE_NUMBER_FEATURES)
{
  // returns 1 if all parents are equal, and otherwise returns zero.
  // Note that this works for any number of parents.
  if ( numParents == 0 ) {
    error("CDT allParentsEqual called with zero features. Need to have at least 1.");
  }
  DiscRVType rv = p0; 
  for (unsigned i = 1 ; i < numParents ; i++ ) {
    if (par(i) != rv)
      return (DiscRVType) 0;
  }
  return (DiscRVType) 1;
}
//
DEFINE_DETERMINISTIC_MAPPER_C_CODE(allParentsUnEqual,CDT_VARIABLE_NUMBER_FEATURES)
{
  // returns 1 if all parents are *un*equal, and otherwise returns zero.
  // Note that this works for any number of parents.
  if ( numParents == 0 ) {
    error("CDT allParentsUnEqual called with zero features. Need to have at least 1.");
  }
  // return 0 if we find any two parents that are equal. 
  // TODO: Is there a faster way than O(N^2) to do this?
  for (unsigned i = 0 ; i < numParents ; i++ ) {
    for (unsigned j = i+1 ; j < numParents ; j++ ) {
      if (par(i) == par(j))
	return (DiscRVType) 0;
    }
  }
  return (DiscRVType) 1;
}
// DO NOT MODIFY ANYTHING IN THE ABOVE BIT OF CODE ENDING HERE.
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// DEFINITION_OF_USER_DEFINED_DETERMINISTIC_MAPPERS: 
// Additional user defined DTs. A few examples are given, you can uncomment and
// modify at will, and then recompile GMTK and these deterministic functions
// will be available to you to use just like any decision tree-based deterministic
// mapping.

// note: when a, b are integers, a/b = floor(a/b).
// To get ceil and round, using only integer ops, we have that:
//    ceil((float)a/(float)b) = (a + a - 1)/b
//    round((float)a/(float)b) = (a + a/2)/b = (a + (a>>1))/b.
// One is free here to convert to/from floating point, but often one need not do that.

#define INCREMENT_POSITION_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(increment_position,INCREMENT_POSITION_NUM_FEATURES)
{ 
  // increment p0 based on value of p1, 
  // equivalent to the following dt
  // 2
  //   -1 { min(p0 + p1, PEAK_LENGTH-1) }

  const DiscRVType max_peaks=PEAK_LENGTH-1;
  if(p0+p1 < max_peaks) //   return min(p0+p1, max_peaks)
    return (p0+p1);
  else
    return (max_peaks);
}



#define LAST_DIGIT_ONE_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(last_digit_one, LAST_DIGIT_ONE_NUM_FEATURES) {
    const DiscRVType mod_ten = p0 % 10;
    if (mod_ten == 0) return (DiscRVType) 0;
    return (DiscRVType) 1;
}

#define DIVIDE_TEN_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(divide_ten, DIVIDE_TEN_NUM_FEATURES) {
    return p0 / 10;
}



#define DECREMENT_INSERTION_COUNT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(decrement_insertion_count,DECREMENT_INSERTION_COUNT_NUM_FEATURES)
{ 
  // equivalent to the following dt
  // decrement_insertion_count
  // 2
  // 1 2 0 default
  //   -1 { p0 }
  //   -1 { max(p0-p1, 0) }

  const DiscRVType min_insertions=0;
  
  if(p1) {
    if(p0-p1 > min_insertions) //    return max(p0-p1, min_insertions)
      return (p0-p1);
    else
      return (min_insertions);
  }
  else
    return p0;
}

#define CHECK_INSERTION_COUNT_MAPPING_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(check_insertion_count,CHECK_INSERTION_COUNT_MAPPING_NUM_FEATURES)
{ 
  // equivalent to the following dt
  // checkInsertionCount
  // 1
  //   -1 { ( p0 == 0 ? 1 : 0 ) } 

  if(p0)
    return (DiscRVType) 0;
  else
    return (DiscRVType) 1;
}

#define NUM_FRAMES_INSERTION_CPT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(num_frames_insertion_cpt,NUM_FRAMES_INSERTION_CPT_NUM_FEATURES)
{ 
  //equivalent to the following dt:
  // num_frames_switch_insertion_cpt %used in prologue
  //   2
  //   -1 { (p0 == 0 ? 0 : p1) }
  if(p0)
    return p1;
  else
    return (DiscRVType) 0;
}

#define DECREMENT_DELETE_COUNT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(decrement_delete_count,DECREMENT_DELETE_COUNT_NUM_FEATURES)
{ 
  //equivalent to the following dt:
  // decrement_delete_count
  // 2
  // 1 3 0 1 default
  //   -1 { p0 }
  //   -1 { p0 }
  //   -1 { max(p0-p1+1, 0) }

  // if(p1 < 2)
  //   return p0;
  // else { //return max(p0-p1+1, 0)
  //   // if(p0-p1+1 > min_deletions)
  //   if(p0-p1+1 > 0)
  //     return (p0-p1+1);
  //   else
  //     return (DiscRVType) 0;
  // }
  DiscRVType rv = p0-(p1-1);
  if(p1 < 2)
    return p0;
  if( rv >= 0)
    return rv;
  if ( rv < 0)
    return (DiscRVType) 0;
}

#define SWITCH_DELTA_DISTRIBUTION_NUM_FEATURES 4
DEFINE_DETERMINISTIC_MAPPER_C_CODE(switch_delta_distribution,SWITCH_DELTA_DISTRIBUTION_NUM_FEATURES)
{ 
// switch_delta_distribution
// 4
//   % parents: p0=num_peaks_minus_one(0), p1=peak_position(-1), p2=insertion(-1), p3=delete_count(-1)
// 2 2 0 default
//     -1 { (( p1+1 > p0 ) ? 0 : min(p3+1,DELTA_STATES-1) )} %n-1 deletes left, switch to\
//  appropriate distribution
//    -1 { 0 } %do not increment
  DiscRVType rv = DELTA_STATES-1;
  if(p2) //insertion, do not allow delta to fire
    return (DiscRVType) 0;
  else {
    if(p1+1 > p0)
      return (DiscRVType) 0;
    else {
      if(rv < p3+1)
	return (rv);
      else
	return (p3+1);
    }
  }
}

#define SWITCH_DELTA_DISTRIBUTION_CONSDELS50_NUM_FEATURES 4
DEFINE_DETERMINISTIC_MAPPER_C_CODE(switch_delta_distribution_consdels50,SWITCH_DELTA_DISTRIBUTION_CONSDELS50_NUM_FEATURES)
{ 
// switch_delta_distribution
// 4
//   % parents: p0=num_peaks_minus_one(0), p1=peak_position(-1), p2=insertion(-1), p3=delete_count(-1)
// 2 2 0 default
//     -1 { (( p1+1 > p0 ) ? 0 : min(p3+1,DELTA_STATES-1) )} %n-1 deletes left, switch to\
//  appropriate distribution
//    -1 { 0 } %do not increment
  DiscRVType rv = 51;
  if(p2) //insertion, do not allow delta to fire
    return (DiscRVType) 0;
  else {
    if(p1+1 > p0)
      return (DiscRVType) 0;
    else {
      if(rv < p3+1)
	return (rv);
      else
	return (p3+1);
    }
  }
}

#define SWITCH_DELTA_DISTRIBUTION_CONSDELS100_NUM_FEATURES 4
DEFINE_DETERMINISTIC_MAPPER_C_CODE(switch_delta_distribution_consdels100,SWITCH_DELTA_DISTRIBUTION_CONSDELS100_NUM_FEATURES)
{ 
// switch_delta_distribution
// 4
//   % parents: p0=num_peaks_minus_one(0), p1=peak_position(-1), p2=insertion(-1), p3=delete_count(-1)
// 2 2 0 default
//     -1 { (( p1+1 > p0 ) ? 0 : min(p3+1,DELTA_STATES-1) )} %n-1 deletes left, switch to\
//  appropriate distribution
//    -1 { 0 } %do not increment
  DiscRVType rv = 101;
  if(p2) //insertion, do not allow delta to fire
    return (DiscRVType) 0;
  else {
    if(p1+1 > p0)
      return (DiscRVType) 0;
    else {
      if(rv < p3+1)
	return (rv);
      else
	return (p3+1);
    }
  }
}

#define DELTA_DEL_VE_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(delta_del_ve, DELTA_DEL_VE_NUM_FEATURES)
{ 
  //p0 = DELTA(0), p1 = DELETE_COUNT(0)
    if(p0 + 1 > p1) // delta must respect the max deletion constraint
    return (DiscRVType) 0;
  else // otherwise, valid hypothesis
    return (DiscRVType) 1;
}

#define DELTA_THEO_BOUND_NUM_FEATURES 3
DEFINE_DETERMINISTIC_MAPPER_C_CODE(delta_theo_bound,DELTA_THEO_BOUND_NUM_FEATURES)
{ 
  // if an insertion occurred in the previous frame, do not let delta fire
  // if we are at the last peak in the theoretical spectrum, do not allow delta to fire
  // else, let delta be random
  // NUM_PEAKS_MINUS_ONE(0), PEAK_POSITION(-1), INSERTION(-1)
  if(p2) //insertion, do not allow delta to fire
    return (DiscRVType) 0;
  else {
    if(p1+1 > p0)
      return (DiscRVType) 0;
    else {
      return (DiscRVType) 1;
    }
  }
}

#define DELTA_THEO_BOUND_VE_NUM_FEATURES 3
DEFINE_DETERMINISTIC_MAPPER_C_CODE(delta_theo_bound_ve, DELTA_THEO_BOUND_VE_NUM_FEATURES)
{ 
  // if delta fires and we go beyond max theoretical spectrum index, return zero
  // DELTA(0), PEAK_POSITION(-1), NUM_PEAKS_MINUS_ONE(0)
  if(p0+p1 > p2)
    return (DiscRVType) 0;
  else
    return (DiscRVType) 1;
}

#define SWITCH_DELTA_DISTRIBUTION_CONSDELS24_NUM_FEATURES 4
DEFINE_DETERMINISTIC_MAPPER_C_CODE(switch_delta_distribution_consdels24,SWITCH_DELTA_DISTRIBUTION_CONSDELS24_NUM_FEATURES)
{ 
// switch_delta_distribution
// 4
//   % parents: p0=num_peaks_minus_one(0), p1=peak_position(-1), p2=insertion(-1), p3=delete_count(-1)
// 2 2 0 default
//     -1 { (( p1+1 > p0 ) ? 0 : min(p3+1,DELTA_STATES-1) )} %n-1 deletes left, switch to\
//  appropriate distribution
//    -1 { 0 } %do not increment
  DiscRVType rv = 25;
  if(p2) //insertion, do not allow delta to fire
    return (DiscRVType) 0;
  else {
    if(p1+1 > p0)
      return (DiscRVType) 0;
    else {
      if(rv < p3+1)
	return (rv);
      else
	return (p3+1);
    }
  }
}

#define CHECK_IN_RANGE_NUM_FEATURES 3
DEFINE_DETERMINISTIC_MAPPER_C_CODE(check_in_range,CHECK_IN_RANGE_NUM_FEATURES)
{ 
  //variable described: IN_RANGE (epilogue variable)
  //parents in order: PEAK_POSITION(0), NUM_PEAKS_MINUS_ONE(0), MAX_DELETES(0)
  //equivalent to the following dt:
  // checkInRange
  // 2
  //   -1 { (p1-p0 < p2+1 ? 1 : 0 )  }
  
  // const DiscRVType delete_card=9;
  // if(p1-p0 < delete_card)
  //   return (DiscRVType) 1;
  // else
  //   return (DiscRVType) 0;

  if(p1-p0 < p2+1)
    return (DiscRVType) 1;
  else
    return (DiscRVType) 0;
}

#define CONS_LIM_CHECK_INSERTION_COUNT_NUM_FEATURES 3
DEFINE_DETERMINISTIC_MAPPER_C_CODE(cons_lim_check_insertion_count,CONS_LIM_CHECK_INSERTION_COUNT_NUM_FEATURES)
{ 
  //variable described: INSERTION
  //parents in order: INSERTION_COUNT(0), CONS_INS(0), CONS_INS_LIM(0)
  //equivalent to the following dt:
// checkInsertionCount
// 3
// 0 2 0 default
//   -1 { 1 }
//   1 2 0 default
//     -1 { 0 } % insertions left, nonconsceutive insertion, fire
// 		 -1 { ( p1 > p2 ? 1 : 0 ) } %if consecutive insertion limit hit, do not allow insertions.  else,\
//  let insertion be bernoulli
  if(p0){
    if(p1){
      if(p1 > p2)
	return (DiscRVType) 1;
      else
	return (DiscRVType) 0;	
    }
    else
      return (DiscRVType) 0;
  }
  else
    return (DiscRVType) 1;
}

#define PROLOGUE_SWITCH_DELTA_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(prologue_switch_delta,PROLOGUE_SWITCH_DELTA_NUM_FEATURES)
{ 
  //variable described: delta (in the prologue)
  //parents in order: MAX_DELETES(0), CONS_DEL_LIM(0)
  //equivalent to the following dt:
// prologue_switch_delta
// 2
// 0 1 default
//   -1 { ( p0 > p1 ? p1 : p0 )} %if max deletes > max consecutive deletes, let delta be mutlinomial o\
// ver max consecutive deletes states
  if(p0 > p1)
    return p1;
  else
    return p0;
}

#define CONS_LIM_SWITCH_DELTA_DISTRIBUTION_NUM_FEATURES 6
DEFINE_DETERMINISTIC_MAPPER_C_CODE(cons_lim_switch_delta_distribution,CONS_LIM_SWITCH_DELTA_DISTRIBUTION_NUM_FEATURES)
{ 
  //variable described: delta (chunk and epilogue)
  //parents in order: NUM_PEAKS_MINUS_ONE(0), PEAK_POSITION(-1), INSERTION(-1), DELETE_COUNT(-1), CONS_DELS(0), CONS_DEL_LIM(0)
  //equivalent to the following dt:
// switch_delta_distribution
// 6
// 2 2 0 default
//   -1 { ( (p1+1 > p0) ? 0 : ( (p4>p5) ? 0 : min(p3+1, p5-p4+1)) ) } %1st condition: at last peak, 2c\
// 								      nd condition: consecutive deletes used up
// 								      -1 { 0 } %insertion in previous frame, do not fire
  if(p2)
    return (DiscRVType) 0;
  else {
    if(p1+1 > p0)
      return (DiscRVType) 0;
    else {
      if(p4 > p5)
	return (DiscRVType) 0;
      else {
	if(p3+1 < p5-p4+1)
	  return (p3+1);
	else
	  return (p5-p4+1);
      }
    }
  }
}

#define CONS_DELS_DT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(cons_dels_dt,CONS_DELS_DT_NUM_FEATURES)
{ 
  //variable described: cons_dels cpt
  //parents in order: DELTA(-1), CONS_DELS(-1)
  //equivalent to the following dt:
// cons_dels_dt
// 2
// 0 3 0 1 default
//   -1 { 0 }
//  -1 { 0 }
//  -1 { ( p0+p1-1 ) }

 if(p0 <= 1)
   return (DiscRVType) 0;
 else
   return (p0+p1-1);
}

#define CONS_INS_DT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(cons_ins_dt,CONS_INS_DT_NUM_FEATURES)
{ 
  //variable described: CONS_INS
  //parents in order: INSERTION(-1), CONS_INS(-1)
  //equivalent to the following dt:
// cons_ins_dt
// 2
// 0 2 0 default
//   -1 { 0 }
//  -1 { ( p1+1 ) }
  
  if(p0)
    return (p1+1);
  else
    return (DiscRVType) 0;
}

#define MAX_DELTA_DELS_DT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(max_delta_dels_dt,MAX_DELTA_DELS_DT_NUM_FEATURES)
{ 
  //variable described: delete_count
  //parents in order: p0 = max_dels
  //equivalent to the following dt:
  // max_delta_dels %used in prologue
  // 1
  //   %since a jump of 1 is now a deletion, knocks off two switching states
  //   -1 { min(p0, DELTA_STATES - 2) }
  
  DiscRVType rv = DELTA_STATES-2;
  if(p0 < rv)
    return (p0);
  else
    return (rv);
}

#define MAX_DELTA_DELS_CONSDELS50_DT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(max_delta_dels_consdels50_dt,MAX_DELTA_DELS_CONSDELS50_DT_NUM_FEATURES)
{ 
  //variable described: delete_count
  //parents in order: p0 = max_dels
  //equivalent to the following dt:
  // max_delta_dels %used in prologue
  // 1
  //   %since a jump of 1 is now a deletion, knocks off two switching states
  //   -1 { min(p0, DELTA_STATES - 2) }
  
  DiscRVType rv = 50;
  if(p0 < rv)
    return (p0);
  else
    return (rv);
}

#define MAX_DELTA_DELS_CONSDELS100_DT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(max_delta_dels_consdels100_dt,MAX_DELTA_DELS_CONSDELS100_DT_NUM_FEATURES)
{ 
  //variable described: delete_count
  //parents in order: p0 = max_dels
  //equivalent to the following dt:
  // max_delta_dels %used in prologue
  // 1
  //   %since a jump of 1 is now a deletion, knocks off two switching states
  //   -1 { min(p0, DELTA_STATES - 2) }
  
  DiscRVType rv = 100;
  if(p0 < rv)
    return (p0);
  else
    return (rv);
}

#define MAX_DELTA_DELS_CONSDELS24_DT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(max_delta_dels_consdels24_dt,MAX_DELTA_DELS_CONSDELS24_DT_NUM_FEATURES)
{ 
  //variable described: delete_count
  //parents in order: p0 = max_dels
  //equivalent to the following dt:
  // max_delta_dels %used in prologue
  // 1
  //   %since a jump of 1 is now a deletion, knocks off two switching states
  //   -1 { min(p0, DELTA_STATES - 2) }
  
  DiscRVType rv = 24;
  if(p0 < rv)
    return (p0);
  else
    return (rv);
}

#define EPI_INSERTION_COUNT_DT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(epi_insertion_count_dt,EPI_INSERTION_COUNT_DT_NUM_FEATURES)
{ 
// epi_insertion_count_dt
//   % parents: p0=insertion_count(0), p1=insertion(0)
// 2
//     -1 { min(p0-p1, MAX_EPI_INSERTIONS-1) }
  
  DiscRVType rv = MAX_EPI_INSERTIONS-1;
  if(p0-p1 < rv)
    return (p0-p1);
  else
    return (rv);
}

#define DEL_PRUNE_DT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(del_prune_dt,DEL_PRUNE_DT_NUM_FEATURES)
{ 
//   % parents: p0=delete_count(0), p1=delete_lower_bound(0)
// 2
//     -1 { (p0 <= p1) ? 0 : 1 }

  if(p0 <= p1)
    return (DiscRVType) 0;
  else
    return (DiscRVType) 1;
}

#define INS_PRUNE_TOL0_DT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(ins_prune_tol0_dt,INS_PRUNE_TOL0_DT_NUM_FEATURES)
{ 
  // Fit a window around the estimated insertions
  //   % parents: p0=insertion_count(0), p1=ins_bound(0), p2=ins_tol(0)
  // 2
  //     -1 { (abs(p0-p1) <= INS_TOL) ? 1 : 0 }
  if(p0==p1)
    return (DiscRVType) 1;
  else
    return (DiscRVType) 0;
}

#define INS_PRUNE_TOL1_DT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(ins_prune_tol1_dt,INS_PRUNE_TOL1_DT_NUM_FEATURES)
{ 
  // Fit a window around the estimated insertions
  //   % parents: p0=insertion_count(0), p1=ins_bound(0), p2=ins_tol(0)
  // 2
  //     -1 { (abs(p0-p1) <= INS_TOL) ? 1 : 0 }
  DiscRVType rv0 = 0;
  //  DiscRVType rv1 = INS_TOL;

  // rv0 = |p0-p1|
  if(p0 > p1)
    rv0 = p0-p1;
  else
    rv0 = p1-p0;

  //  if(rv0 <= rv1)
  if(rv0 <= 1)
    return (DiscRVType) 1;
  else
    return (DiscRVType) 0;
}

#define INS_PRUNE_TOL2_DT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(ins_prune_tol2_dt,INS_PRUNE_TOL2_DT_NUM_FEATURES)
{ 
  // Fit a window around the estimated insertions
  //   % parents: p0=insertion_count(0), p1=ins_bound(0), p2=ins_tol(0)
  // 2
  //     -1 { (abs(p0-p1) <= INS_TOL) ? 0 : 1 }
  DiscRVType rv0 = 0;
  //  DiscRVType rv1 = INS_TOL;

  // rv0 = |p0-p1|
  if(p0 > p1)
    rv0 = p0-p1;
  else
    rv0 = p1-p0;

  //  if(rv0 <= rv1)
  if(rv0 <= 2)
    return (DiscRVType) 1;
  else
    return (DiscRVType) 0;
}

#define BY_INTENSITY_DT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(by_intensity_dt,BY_INTENSITY_DT_NUM_FEATURES)
{ 
// parents: p0 == insertion(0), p1 == is_bion(0)
// dt:
// 2
// 0 2 1 default % split based on whether peak is an insertion
//   -1 2 % insertion, return constant penalty, which is state 0
//   -1 { p1 } % score the peak, p1=0 denotes a y-ion, p1=1 denotes a b-ion

  if(p0) // insertion
    return (DiscRVType) 2;
  else // not an insertion
    return (p1);
}

#define CHARGE_MZ_DT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(charge_mz_dt,CHARGE_MZ_DT_NUM_FEATURES)
{ 
// charge_mz_dt  % for now, assume charge \in {2,3}
// 2 % p0 == insertion(0), p1 = fragment_charge(0)
// 0 2 1 default
//   -1 { 2 }
//   -1 { p1 - 1 }

  if(p0) //insertion
    return (DiscRVType) 2;
  else
    return (p1-1);
}

#define BY_CHARGE_INTENSITY_DT_NUM_FEATURES 3
DEFINE_DETERMINISTIC_MAPPER_C_CODE(by_charge_intensity_dt,BY_CHARGE_INTENSITY_DT_NUM_FEATURES)
{ 
// by_charge_intensity_dt  % for now, assume charge \in {2,3}
// 3
// 0 2 1 default % insertion
//   -1 { 4 }
//   2 2 1 default % fragment_charge
//     -1 { p1 }
//     -1 { p1 + 2 }
  if(p0) //insertion
    return (DiscRVType) 4;
  else
    if(p2==1)
      return p1;
    else
      return (p1+2);
}
//
// add more deterministic mapping functions here as desired ...
// 
#ifdef USER_INTERNAL_CFUNC_DTS
#include "user_internal_cfunc_dts.cc"
#endif

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Registration code. This registers all deterministic mapping
// functions that were defined above. Once  they are registered, with a given
// name, they may be used like any decision tree.
//
// NOTE: please give all decision trees registered as such a name
// starting with "internal:" and any user defined functions
// a name starting with "user_internal:"



void
registerAllCFunctionDeterministicMappings(GMParms& gmp)
{

  ///////////////////////////////////////////////////////////////////////
  // DO NOT CHANGE ANYTHING IN THE FOLLOWING FEW LINES STARTING HERE.
  gmp.registerDeterministicCMapper("internal:copyParent",
				   COPYPARENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(copyParent));
  gmp.registerDeterministicCMapper("internal:alwaysZero",
				   ALWAYSZERO_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(alwaysZero));
  gmp.registerDeterministicCMapper("internal:alwaysOne",
				   ALWAYSONE_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(alwaysOne));
  gmp.registerDeterministicCMapper("internal:increment",
				   INCREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(increment));
  gmp.registerDeterministicCMapper("internal:decrement",
				   DECREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(decrement));
  gmp.registerDeterministicCMapper("internal:conditionalIncrement",
				   CONDITIONAL_INCREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(conditionalIncrement));
  gmp.registerDeterministicCMapper("internal:conditionalDecrement",
				   CONDITIONAL_DECREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(conditionalDecrement));
  gmp.registerDeterministicCMapper("internal:conditionalLimitedIncrement",
				   CONDITIONAL_LIMITED_INCREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(conditionalLimitedIncrement));
  gmp.registerDeterministicCMapper("internal:conditionalLimitedDecrement",
				   CONDITIONAL_LIMITED_DECREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(conditionalLimitedDecrement));
  gmp.registerDeterministicCMapper("internal:allParentsUnEqual",
				   CDT_VARIABLE_NUMBER_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(allParentsUnEqual));
  // TODO: add the variable parent deterministic mappers once the DT code can accept variable numbers of parens.
  // ...
  // DO NOT CHANGE ANYTHING IN THE ABOVE FEW LINES ENDING HERE.
  ///////////////////////////////////////////////////////////////////////




  ///////////////////////////////////////////////////////////////////////
  // REGISTRATION_OF_USER_DEFINED_DETERMINISTIC_MAPPERS: 
  // ADD USER DEFINED C FUNCTION DETERMINISTIC MAPPING REGISTRATIONS HERE.
  // Arguments are:
  //   registerDeterministicCMapper(
  //          name_of_deterministic_mapping which is of type char*,
  //          number of features (just like when defining a decision tree),
  //          C function above, use macro
  //        );

  /*
   * Here is an example mapping registration that one might want to define.
   * Uncomment to activate.
   gmp.registerDeterministicCMapper("user_internal:ajitMapping",
   AJIT_MAPPING_NUM_FEATURES,
   DETERMINISTIC_MAPPER_C_CODE_NAME(ajitMapping));
  */
  gmp.registerDeterministicCMapper("user_internal:increment_position",
				   INCREMENT_POSITION_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(increment_position));
  gmp.registerDeterministicCMapper("user_internal:decrement_insertion_count",
				   DECREMENT_INSERTION_COUNT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(decrement_insertion_count));
  gmp.registerDeterministicCMapper("user_internal:check_insertion_count",
				   CHECK_INSERTION_COUNT_MAPPING_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(check_insertion_count));
  gmp.registerDeterministicCMapper("user_internal:num_frames_insertion_cpt",
				   NUM_FRAMES_INSERTION_CPT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(num_frames_insertion_cpt));
  gmp.registerDeterministicCMapper("user_internal:decrement_delete_count",
				   DECREMENT_DELETE_COUNT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(decrement_delete_count));
  gmp.registerDeterministicCMapper("user_internal:switch_delta_distribution",
				   SWITCH_DELTA_DISTRIBUTION_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(switch_delta_distribution));
  gmp.registerDeterministicCMapper("user_internal:switch_delta_distribution_consdels24",
				   SWITCH_DELTA_DISTRIBUTION_CONSDELS24_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(switch_delta_distribution_consdels24));
  gmp.registerDeterministicCMapper("user_internal:switch_delta_distribution_consdels50",
				   SWITCH_DELTA_DISTRIBUTION_CONSDELS50_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(switch_delta_distribution_consdels50));
  gmp.registerDeterministicCMapper("user_internal:switch_delta_distribution_consdels100",
				   SWITCH_DELTA_DISTRIBUTION_CONSDELS100_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(switch_delta_distribution_consdels100));
  gmp.registerDeterministicCMapper("user_internal:check_in_range",
				   CHECK_IN_RANGE_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(check_in_range));
  gmp.registerDeterministicCMapper("user_internal:cons_lim_check_insertion_count",
				   CONS_LIM_CHECK_INSERTION_COUNT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(cons_lim_check_insertion_count));
  gmp.registerDeterministicCMapper("user_internal:prologue_switch_delta",
				   PROLOGUE_SWITCH_DELTA_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(prologue_switch_delta));
  gmp.registerDeterministicCMapper("user_internal:cons_lim_switch_delta_distribution",
				   CONS_LIM_SWITCH_DELTA_DISTRIBUTION_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(cons_lim_switch_delta_distribution));
  gmp.registerDeterministicCMapper("user_internal:cons_dels_dt",
				   CONS_DELS_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(cons_dels_dt));
  gmp.registerDeterministicCMapper("user_internal:cons_ins_dt",
				   CONS_INS_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(cons_ins_dt));
  gmp.registerDeterministicCMapper("user_internal:max_delta_dels_dt",
				   MAX_DELTA_DELS_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(max_delta_dels_dt));
  gmp.registerDeterministicCMapper("user_internal:max_delta_dels_consdels24_dt",
				   MAX_DELTA_DELS_CONSDELS24_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(max_delta_dels_consdels24_dt));
  gmp.registerDeterministicCMapper("user_internal:max_delta_dels_consdels50_dt",
				   MAX_DELTA_DELS_CONSDELS50_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(max_delta_dels_consdels50_dt));
  gmp.registerDeterministicCMapper("user_internal:max_delta_dels_consdels100_dt",
				   MAX_DELTA_DELS_CONSDELS100_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(max_delta_dels_consdels100_dt));
  gmp.registerDeterministicCMapper("user_internal:epi_insertion_count_dt",
				   EPI_INSERTION_COUNT_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(epi_insertion_count_dt));
  gmp.registerDeterministicCMapper("user_internal:del_prune_dt",
				   DEL_PRUNE_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(del_prune_dt));
  gmp.registerDeterministicCMapper("user_internal:ins_prune_tol0_dt",
				   INS_PRUNE_TOL0_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(ins_prune_tol0_dt));
  gmp.registerDeterministicCMapper("user_internal:ins_prune_tol1_dt",
				   INS_PRUNE_TOL1_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(ins_prune_tol1_dt));
  gmp.registerDeterministicCMapper("user_internal:ins_prune_tol2_dt",
				   INS_PRUNE_TOL2_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(ins_prune_tol2_dt));
  gmp.registerDeterministicCMapper("user_internal:by_intensity_dt",
				   BY_INTENSITY_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(by_intensity_dt));
  gmp.registerDeterministicCMapper("user_internal:charge_mz_dt",
				   CHARGE_MZ_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(charge_mz_dt));
  gmp.registerDeterministicCMapper("user_internal:by_charge_intensity_dt",
				   BY_CHARGE_INTENSITY_DT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(by_charge_intensity_dt));
  gmp.registerDeterministicCMapper("user_internal:delta_del_ve",
				   DELTA_DEL_VE_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(delta_del_ve));
  gmp.registerDeterministicCMapper("user_internal:delta_theo_bound",
				   DELTA_THEO_BOUND_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(delta_theo_bound));
  gmp.registerDeterministicCMapper("user_internal:delta_theo_bound_ve",
				   DELTA_THEO_BOUND_VE_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(delta_theo_bound_ve));

  gmp.registerDeterministicCMapper("user_internal:last_digit_one", 
                   LAST_DIGIT_ONE_NUM_FEATURES,
                   DETERMINISTIC_MAPPER_C_CODE_NAME(last_digit_one));

  gmp.registerDeterministicCMapper("user_internal:divide_ten", 
                   DIVIDE_TEN_NUM_FEATURES,
                   DETERMINISTIC_MAPPER_C_CODE_NAME(divide_ten));

  // Uncomment to register user defined DTs. You can defiine them in the below file if you like.
  // 

#ifdef USER_INTERNAL_CFUNC_DTS
#include "register_user_internal_cfunc_dts.cc"
#endif


}
