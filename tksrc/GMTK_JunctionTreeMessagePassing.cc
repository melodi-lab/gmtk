/*-
 * GMTK_JunctionTree.cc
 *     Junction Tree, message passing routines.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003-2009.
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */


/* The ISO C99 standard specifies that in C++ implementations these
   macros should only be defined if explicitly requested.  */
#define __STDC_LIMIT_MACROS 1
#define __STDC_CONSTANT_MACROS 1
   // The ISO C99 standard specifies that the macros in inttypes.h must
   //  only be defined if explicitly requested. 
#define __STDC_FORMAT_MACROS 1

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>
#endif
#if HAVE_STDINT_H
#  include <stdint.h>
#endif


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <ctype.h>

#include <iterator>
#include <map>
#include <set>
#include <algorithm>
#include <new>
#include <typeinfo>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_GMParms.h"
#include "GMTK_Dlinks.h"
#include "GMTK_DeepVECPT.h"
#include "GMTK_BinaryViterbiFileUtils.h"

#include "GMTK_RngDecisionTree.h"

#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_FileSource.h"
#  include "GMTK_StreamSource.h"
#endif

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)

unsigned JunctionTree::numSmoothingPartitions = 0; // used for smoothing in gmtkOnline

/*
 *  init_CC_CE_rvs(it)
 *    given an iterator representing an unrolled segment,
 *    we create the pair of random variables CC and CE that
 *    will be shifted around in time/position in order
 *    to do inference.
 *
 * Preconditions:
 *    1) partitionStructureArray must be set up and initialized before
 *    this routine is called.
 *    2) All random variables should currently be shifted to their 
 *       zero position (otherwise thigns will be totally out of sync).
 *
 */
void
JunctionTree::init_CC_CE_rvs(ptps_iterator& ptps_it)
{
  if (ptps_it.num_c_partitions() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // partitions (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = partitionStructureArray[1].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,cur_CC_rvs);
    tmp1 = partitionStructureArray[3].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,cur_CE_rvs);
  } else {
    cur_CC_rvs.clear();
    cur_CE_rvs.clear();
  } 
  cur_cc_shift = cur_ce_shift = 0;
  ptps_it.go_to_part_no(0);
}

/*
 * shiftCCtoPosition(int pos)
 *  Given an absolute partition number 'pos', we shift the pair of
 *  random variables CC so that the right position (the 2nd
 *  C of CC) is now at position 'pos'. If the CE is shifted,
 *  we first shift CE back to position 0 before shifting CC.
 */
void
JunctionTree::shiftCCtoPosition(int pos)
{
  if (cur_ce_shift != 0) {
    assert ( cur_cc_shift == 0 );
    // we need to get CE back to zero.
    adjustFramesBy (cur_CE_rvs, -cur_ce_shift
		    *gm_template.S*fp.numFramesInC(),
		    !onlineViterbi);  // gmtkOnline can't observe old values - ticket #468
    cur_ce_shift = 0;
  }
  int delta = pos - cur_cc_shift;
  adjustFramesBy (cur_CC_rvs, delta
		  *gm_template.S*fp.numFramesInC());
  cur_cc_shift = pos;
}

/*
 * shiftOriginalVarstoPosition(vector<RV*> rvs, int pos, int &prevPos)
 *  Shift the random variables in 'rvs' in time (frames) by
 *  the difference between 'prevPos' and 'pos'. This is used  in 
 *  printSavedViterbiValues to adjust the rvs' frame numbers in
 *  the unpacking buffers.
 */

void
JunctionTree::shiftOriginalVarstoPosition(vector<RV*> rvs, int pos, int &prevPos)
{
  set<RV*> uprvs(rvs.begin(),rvs.end());
  int delta = (pos - prevPos);
  adjustFramesBy(uprvs, delta);
  prevPos = pos;
}

/*
 * shiftCCtoPosition(int pos)
 *  Given an absolute partition number 'pos', we shift the pair of
 *  random variables CE so that the right position (the E
 *  of CE) is now at position 'pos'. If the CC pair is shifted,
 *  we first shift CC back to position 0 before shifting CE.
 */
void
JunctionTree::shiftCEtoPosition(int pos)
{
  if (cur_cc_shift != 0) {
    assert ( cur_ce_shift == 0 );
    // we need to get CC back to zero.
    adjustFramesBy (cur_CC_rvs, -cur_cc_shift
		    *gm_template.S*fp.numFramesInC(),
                    !onlineViterbi);  // gmtkOnline can't observe old values - ticket #468
    cur_cc_shift = 0;
  }
  int delta = pos - cur_ce_shift;
  adjustFramesBy (cur_CE_rvs, delta
		  *gm_template.S*fp.numFramesInC());
  cur_ce_shift = pos;
}



/*
 * setCurrentInferenceShiftTo(int pos)
 *
 *  This activates partitions 'pos-1' and 'pos' so the
 *  right of two successive partitions is active at position 'pos'.
 *  This routine could be called:
 *
 *      activateTwoAdjacentPartitionsWithRightPartitionAtPos(pos)
 *      alignRightSideOfPartitionPairToPos(pos)
 * 
 *  Given an absolute partition number 'pos', we shift the pair of
 *  random variables (either CC or CE) so the *right* position (the 2nd
 *  C of CC or the E of CE) is now at position 'pos'. This means that
 *  messages entirely within position 'pos' and entirely within
 *  position 'pos-1' and between positions 'pos-1' and 'pos' will be
 *  correct, but no other messages are guaranteed to be correct.
 *
 * preconditions:
 *    the class member iterator 'inference_it' must be initilzed for the
 *    current and appropriate segment length.
 *
 * side effects
 *   - modifies the random variables corresponding to a partition pair
 *   - modifies the inference_it iterator.
 *
 */
void
JunctionTree::setCurrentInferenceShiftTo(int pos)
{

  if (inference_it.at_entry(pos)) {
    // printf("== Already at shift %d\n",pos);
    return;
  }

  //   printf("========================================\n");
  //   printf("== Setting current inference shift to %d\n",pos);
  //   printf("========================================\n");

  inference_it.go_to_part_no(pos);
  if (inference_it.num_c_partitions() <= 2) {
    // then do nothing, since nothing is never shifted away from zero.
  } else {
    // need to do some work. We need to get the right
    // of the appropriate pair (either the second C of CC in PCCE,
    // or the E of CE in PCCE) to be at position pos.
    if (inference_it.at_p()) {
      // P has an intersection with the first C in PCCE, so we need
      // to get that first C into the right position. We do that
      // by getting both CC's into the right position.
      shiftCCtoPosition(0);
    } else if (inference_it.at_e()) {
      // get the CE into the right position
      shiftCEtoPosition(pos - 3);
    } else {
      // Get CC into the right position, where 'pos' corresponds
      // to the second C in PCCE.
      if (pos <= 2)
	shiftCCtoPosition(0);
      else 
	shiftCCtoPosition(pos-2);
    }
  }
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//   Standard Collect/Distribute Evidence Support
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////



/*
 *
 * Preconditions: 
 *    The linear algorithm deScatterOutofRoot( .. ) 
 *    must have *just* been called and the random variables, and the iterator it has
 *    been set appropriately to the current partition.
 * Relies on:
 *   1) That 'it' is set to the current partition.
 *   2) That all random variables in the current structure (as designated by 'it') are set to 
 *      the current viterbi value.
 *
 */
void
JunctionTree::recordPartitionViterbiValue(ptps_iterator& it)
{
  PartitionStructures& ps = partitionStructureArray[it.ps_i()];
  unsigned partitionLength = ps.packer.packedLen();
  unsigned N_best = 1;
  unsigned num_to_write = N_best * partitionLength;
  if (partitionLength > 0)  {
    // if it is not greater than zero, then the partition has
    // no hidden discrete variables.

    if (binaryViterbiFile) {
      // Do a binary dump of the Viterbi values to a file
      // in chronological order for later printing by a separate
      // program. This is O(1) memory and O(T) disk space.

      if (it.at_p()) {
	if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset, SEEK_SET) == (off_t) -1) {
	  char *err = strerror(errno);
	  error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
	}	
	ps.packer.pack(ps.hrvValuePtrs.ptr,P_partition_values.ptr);
	if (fwrite(P_partition_values.ptr, sizeof(unsigned), num_to_write, binaryViterbiFile) 
	    != num_to_write) 
	{
	  char *err = strerror(errno);
	  error("ERROR: write failed on '%s': %s\n", binaryViterbiFilename, err);
	}
      } else if (it.at_e()) {
	off_t offset = (off_t)  // P size + (T-2) * C size
	  (   (   N_best * partitionStructureArray[0].packer.packedLen()       
	        + N_best * partitionStructureArray[1].packer.packedLen() * it.num_c_partitions()  
	      ) * sizeof(unsigned)   );                         
	if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset + offset, SEEK_SET) == (off_t) -1) {
	  char *err = strerror(errno);
	  error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
	}	
	ps.packer.pack(ps.hrvValuePtrs.ptr,E_partition_values.ptr);
	if (fwrite(E_partition_values.ptr, sizeof(unsigned), num_to_write, binaryViterbiFile) 
	    != num_to_write) 
	{
	  char *err = strerror(errno);
	  error("ERROR: write failed on '%s': %s\n", binaryViterbiFilename, err);
	}
	nextViterbiOffset = gmtk_ftell(binaryViterbiFile); // remember where to start next segment
	if (nextViterbiOffset == (off_t)-1) {
	  char *err = strerror(errno);
	  error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
	}
      } else { // at a C partition
	off_t offset = (off_t)  // P size + (t-1) * C size
          (   (   N_best * partitionStructureArray[0].packer.packedLen()       
		+ N_best * partitionLength * ( it.pt_i() - 1 )
              ) * sizeof(unsigned)   );               
	if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset + offset, SEEK_SET) == (off_t) -1) {
	  char *err = strerror(errno);
	  error("seek failed on '%s': %s\n", binaryViterbiFilename, err);
	}	
	ps.packer.pack(ps.hrvValuePtrs.ptr, C_partition_values.ptr);  // stomp on old values to be O(1) memory
	if (fwrite(C_partition_values.ptr, sizeof(unsigned), num_to_write, binaryViterbiFile) 
	    != num_to_write) 
	{
	  char *err = strerror(errno);
	  error("write failed on '%s': %s\n", binaryViterbiFilename, err);
	}
      }      
    } else {
      // Store the Viterbi values in ?_partition_values in O(T)
      // memory for printing within this program.

      if (it.at_p()) {
	ps.packer.pack(ps.hrvValuePtrs.ptr,P_partition_values.ptr);
      } else if (it.at_e()) {
	ps.packer.pack(ps.hrvValuePtrs.ptr,E_partition_values.ptr);
      } else {
	ps.packer.pack(ps.hrvValuePtrs.ptr,
		       C_partition_values.ptr
		       + 
		       (it.pt_i()-1)*partitionLength);
      }
    }
  }
}



/////////////////////////////////////////////
// Reuild Viterbi trigger string without whitespace characters
/////////////////////////////////////////////
void
preProcessFormula(string& original) {
  unsigned i;
  string revised;
  for (i=0; i<original.length(); ++i) {
    if (!isspace(original[i])) {
      revised.push_back(original[i]);
    }
  }
  original = revised;
}


typedef pair<string,int>      RVKey;     // the name & offset of an RV
typedef map<RVKey,unsigned>   RVMap;     // map RV name(offset) to pX in leaf node expression
typedef vector<RVKey>         RVVec;     // vector of name(offset) pairs in "parent" order

// This replaces the variableName(offset) references in the triggerExpression with
// pX to match the actual leaf node expression syntax. The rvVec maps the X back to
// the corresponding variableName(offset). variableNames holds the set of variable
// names used in the model (used for error checking).
void
JunctionTree::parseViterbiTrigger(set<string> &variableNames, char *triggerExpression, RVVec &rvVec, string &expr) {

  RVMap rvMap;
  rvVec.clear();
  expr="";
  string trigExprStr(triggerExpression);
  preProcessFormula(trigExprStr);
  triggerExpression = (char *) trigExprStr.c_str();  // no whitespace to deal with
  unsigned varCount = 0;
  for (char *p=triggerExpression; *p; ) {
    // possible start of 'name(offset)'
    if ( ('A' <= *p && *p <= 'Z') || ('a' <= *p && *p <= 'z') ) { 
      string name;
      char  *q;
      // got {alpha}({alpha}|{digit}|_|-)*
      for (q = p; *q && ( ('A' <= *q && *q <= 'Z') || ('a' <= *q && *q <= 'z') || ('0' <= *q && *q <= '9') || *q=='_' || *q=='-'); q+=1)
	name.push_back(*q);
      if (*q == '(') {       // read 'name(' so far
	int    offset;
	string offsetStr;
	q += 1;
	if ('0' <= *q && *q <= '9') {
	  for ( ; *q && ( ('0' <= *q && *q <= '9') ); q+=1)
	    offsetStr.push_back(*q);
	  if (*q == ')') {     // read 'name(offset)'

	    p = q+1;

	    if (variableNames.find(name) != variableNames.end()) {
	      // name matches a variable in the model

	      if (RngDecisionTree::EquationClass::functionNameCollision(name)) {
		warning("WARNING: identifier '%s' in trigger formula '%s' could be either a variable or a single argument function. "
			"It will be interpretted as a variable.\n", name.c_str(), triggerExpression);
	      }
	      sscanf(offsetStr.c_str(), "%d", &offset);
	      
	      RVKey k(name,offset);
	      char parent[16]; // at most 10 digits for a 32-bit uint, plus 'p' and NULL, & round up
	      if (rvMap.find(k) == rvMap.end()) {  // new variable
		rvMap[k] = varCount;
		rvVec.push_back(k);
		sprintf(parent, "p%u", varCount);
		varCount += 1;
	      } else {                             // already encountered this variable
		sprintf(parent, "p%u", rvMap.find(k)->second);
	      }
	      expr.append(parent);
	    } else if (RngDecisionTree::EquationClass::functionNameCollision(name)) {
	      // name does not match any variable in the model, but it is a function name,
	      // so treat it as an unambiguous single-argument function

	      expr.append(name);
	      expr.push_back('(');
	      expr.append(offsetStr);
	      expr.push_back(')');
	    } else {
	      // neither a known variable name, nor a known function name. What is it?
	      error("ERROR: identifier '%s' in trigger formula '%s' is neither a variable name nor a function name.\n",
		    name.c_str(), triggerExpression);
	    }
	  } else {   // read {alpha}({alpha}|{digit}|_|-)* '(' {digit}+ but no ')'
	    expr.append(name);
	    expr.push_back('(');
	    expr.append(offsetStr);
	    p = q;
	  }
	} else {
	  // read {alpha}({alpha}|{digit}|_|-)* '(' but no offset
	    expr.append(name);
	    expr.push_back('(');
	    p = q;
	}
      } else { 
	// read {alpha}({alpha}|{digit}|_|-)* but no '('
	expr.append(name);
	p = q;
      }
    } else {
      expr.push_back(*p++);
    }
  }
}


// output rvVec will be vector of RV* corresponding to the name(offset)s in keyVec for the rvs in the current partition

// keyVec holds the variableName(offset) pairs used in a trigger expression in pX order, i.e. keyVec[X] == pX.
// rvs holds the RV instances of the current section.
// the result is that rvVec[X] is the RV * to pX in the current section

void
RVKey2RVVec(vector<RV*>rvs, RVVec keyVec, vector<RV*> &rvVec) {
  if (rvs.size() == 0) return;

  // find min frame in rvs to use as 0 offset for keyVec
  unsigned minFrame = rvs[0]->frame();
  for (unsigned i = 1; i < rvs.size(); i+=1) {
    if (rvs[i]->frame() < minFrame)
      minFrame = rvs[i]->frame();
  }

  rvVec.clear();
  for (unsigned i=0; i < keyVec.size(); i+=1) {
    bool found = false;
    for (unsigned j=0; j < rvs.size(); j+=1) {
      if (rvs[j]->name() == keyVec[i].first && (int)rvs[j]->frame() == (int)minFrame + keyVec[i].second) {
	rvVec.push_back(rvs[j]);
	found = true;
	break;
      }
    }
    if (!found)
      error("ERROR: Viterbi printing trigger variable %s(%d) not found in section\n", keyVec[i].first.c_str(), keyVec[i].second);
  }
  assert(keyVec.size() == rvVec.size());
#if 0
  for (unsigned i=0; i < keyVec.size(); i+=1) {
    printf("p%u -> %s(%d) => %s(%u)\n", i, keyVec[i].first.c_str(), keyVec[i].second, rvVec[i]->name().c_str(), rvVec[i]->frame());
  }
#endif
}


// returns true if arrays a and b (assumed to both be of length size)
// do not contain exactly the same values where regex_mask[i] is true. 
// Used for Viterbi run-length compression

bool
differentValues(unsigned size, unsigned *a, unsigned *b, vector<bool> &regex_mask) {
  for (unsigned i=0; i < size; i+=1) {
    if (regex_mask[i] && a[i] != b[i]) return true;
  }
  return false;
}


// setup data structures needed to evaluate a Viterbi printing trigger
void
JunctionTree::initializeViterbiTrigger(char *vitTrigger, set<string> &variableNames, 
				       RVVec &vitTriggerVec, string &vitTriggerExpr,
				       RngDecisionTree::EquationClass &triggerEqn,
				       char arg)
{
  if (vitTrigger) {
    parseViterbiTrigger(variableNames, vitTrigger, vitTriggerVec, vitTriggerExpr);
    try {
      triggerEqn.parseFormula(vitTriggerExpr);
    }
    catch( string const &error_message ){
      error("ERROR: In -%cVitTrigger '%s' :  %s", arg, vitTrigger, error_message.c_str());
    }
    catch( const char * const error_message ) {
      error("ERROR: In -%cVitTrigger '%s' :  %s", arg, vitTrigger, error_message);
    }
  }
}


// evaluate a Viterbi printing trigger
bool
JunctionTree::evaluateTrigger(vector<RV *> &allRVs, RVVec &vitTriggerVec, string &vitTriggerExpr, 
			      RngDecisionTree::EquationClass &triggerEqn) 
{
  vector<RV*> triggerParents;
  RVKey2RVVec(allRVs, vitTriggerVec, triggerParents);
  RngDecisionTree dt(vitTriggerExpr);
  return triggerEqn.evaluateFormula(&dt, triggerParents) > 0;
}


// Manages some state for Viterbi run-length encoding. 
// Returns true if the current chunk should be printed.

// first_C is true iff this is the first chunk in the segment => do some initialization.
// C_size is the number of (printed) RVs in the chunk.
// If printObserved is true, print observed RVs instead of just hidden.
//    Note that observed RVs only count for run-length triggering if printObserved is true.
// allrvs is the RV *s for the current chunk (including observed RVs)
// hidRVVector is just the hidden RVs in the current chunk
// previous_C_values are the values of the chunk RVs in the most recently printed chunk
// regex_mask[i] == true iff the ith RV in hidRVVector (allrvs if printObserved is true)
//    matches preg
bool
newViterbiValues(bool &first_C, unsigned &C_size, bool printObserved, 
		 vector<RV *> &allrvs, vector<RV *> hidRVVector, 
		 sArray<unsigned> &previous_C_values, 
		 vector<bool> &regex_mask, regex_t *preg)
{
  if (first_C) {
    first_C = false;
    if (printObserved) {
      C_size = allrvs.size();
    } else {
      C_size = hidRVVector.size();
    }
    previous_C_values.resize(C_size);
    regex_mask.resize(C_size);
    for (unsigned i=0; i < C_size; i+=1) {
      previous_C_values[i] = UINT32_MAX;
      if (printObserved) {
	if (preg) {
	  regex_mask[i] = !regexec(preg, allrvs[i]->name().c_str(),0,0,0);
	} else {
	  regex_mask[i] = true;
	}
      } else {
	if (preg) {
	  regex_mask[i] = !regexec(preg, hidRVVector[i]->name().c_str(),0,0,0);
	} else {
	  regex_mask[i] = true;
	}
      }
    }
  }
  sArray<unsigned> current_C_values(C_size);
  for (unsigned i=0; i < C_size; i+=1) {
    if (printObserved) {
      current_C_values[i] = ((DiscRV *)allrvs[i])->val;
    } else {
      current_C_values[i] = ((DiscRV *)hidRVVector[i])->val;
    }
  }
  bool result =  differentValues(C_size, current_C_values.ptr, previous_C_values.ptr, regex_mask);
  memcpy(previous_C_values.ptr, current_C_values.ptr, C_size * sizeof(unsigned));
  return result;
}




/*
 *
 * This routine saves the viterbi values computed by the most recent
 * linear inference run (assuming its data structures are still valid)
 * to an observation file.
 *
 * Preconditions: 
 *
 *    Assumes that distributeEvidence has just been run and all data
 *    structures (such as the compressed viterbi value array) are set
 *    up appropriately. 
 *
 *    Assumes that inference_it is currently set for the current
 *    segment.
 *  
 *    Assumes that the CC and CE partition pair random variables
 *    have been properly set up.
 * 
 *
 */

void
JunctionTree::viterbiValuesToObsFile(unsigned numFrames,
				     FILE   * binVitFile,
				     unsigned segment,
				     regex_t* preg,
				     regex_t* creg,
				     regex_t* ereg,
				     char   * frameRangeFilter)
{
  unsigned numUsableFrames;
  if (binaryViterbiFile) {
    unsigned totalNumberPartitions;
    numUsableFrames = unroll(numFrames,ZeroTable,&totalNumberPartitions);
    
    new (&inference_it) ptps_iterator(*this,totalNumberPartitions);
    init_CC_CE_rvs(inference_it);
  } else {
    numUsableFrames = this->numUsableFrames;
  }

  // Modified section -> original section initialization

  vector<RV*> unrolled_rvs;
  map<RVInfo::rvParent, unsigned> unrolled_map;

  vector<RV*> P_rvs;      // original P for printing
  vector<RV*> Pprime_rvs; // modified P' for unpacking
  vector<RV*> hidP_rvs;      // hidden subset of original P for printing
  vector<RV*> hidPprime_rvs; // hidden subset of modified P' for unpacking

  vector<vector<RV*> > C_rvs; // original Cs for printing
  vector<vector<RV*> > Cprime_rvs; // modified C's for unpacking
  vector<vector<RV*> > hidC_rvs; // hidden subset of original Cs for printing
  vector<vector<RV*> > hidCprime_rvs; // hidden subset of modified C's for unpacking

  vector<RV*> E_rvs; // ... printing
  vector<vector<RV*> > Eprime_rvs; // ... unpacking
  vector<RV*> hidE_rvs; // ... printing
  vector<vector<RV*> > hidEprime_rvs; // ... unpacking

  sArray<DiscRVType *>PprimeValuePtrs;
  vector<sArray<DiscRVType *> > CprimeValuePtrs;
  vector<sArray<DiscRVType *> > EprimeValuePtrs;

  createUnpackingMap(unrolled_rvs, unrolled_map, 
		     P_rvs, hidP_rvs, Pprime_rvs, hidPprime_rvs, 
		     C_rvs, hidC_rvs, Cprime_rvs, hidCprime_rvs,
		     E_rvs, hidE_rvs, Eprime_rvs, hidEprime_rvs,
		     PprimeValuePtrs, CprimeValuePtrs, EprimeValuePtrs);


  unsigned NP = fp.numFramesInP();
  unsigned NC = fp.numFramesInC();

  unsigned M = gm_template.M;
  unsigned S = gm_template.S;
  unsigned totalOriginalPartitions = 2 + inference_it.num_c_partitions() * S + M;

  infoMsg(IM::Printing,IM::High,"NP = %u   NC = %u   M = %u   S = %u   # orig parts = %u\n",
       NP, NC, M, S, totalOriginalPartitions);

  // Trigger initialization

  set<string> variableNames; // names of variables in the model
  for (unsigned i=0; i < partition_unrolled_rvs.size(); i+=1) {
    variableNames.insert(partition_unrolled_rvs[i]->name());
  }

  RVVec  pVitTriggerVec;
  string pVitTriggerExpr;
  RngDecisionTree::EquationClass pTriggerEqn;
  initializeViterbiTrigger(pVitTrigger, variableNames, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, 'p');

  RVVec  cVitTriggerVec;
  string cVitTriggerExpr;
  RngDecisionTree::EquationClass cTriggerEqn;
  initializeViterbiTrigger(cVitTrigger, variableNames, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, 'c');

  RVVec  eVitTriggerVec;
  string eVitTriggerExpr;
  RngDecisionTree::EquationClass eTriggerEqn;
  initializeViterbiTrigger(eVitTrigger, variableNames, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, 'e');


  Range* frameRange = NULL;
  frameRange = new Range(frameRangeFilter,0,numUsableFrames);
  if (frameRange->length() == 0) { 
    warning("WARNING: Frame range filter must specify a valid non-zero "
	    "length range within [0:%d]. Range given is %s\n",
	    numUsableFrames, frameRangeFilter);
    delete frameRange;
    frameRange = NULL;
  }

  if (frameRange == NULL)
    frameRange = new Range("all",0,numUsableFrames);

  Range::iterator* frameRange_it = new Range::iterator(frameRange->begin());

  vector<int> Cpos(C_rvs.size());
  for (unsigned int i=0; i < Cpos.size(); i+=1) 
    Cpos[i] = fp.numFramesInP() + i * fp.numFramesInC();
  int Epos = fp.numFramesInP() + C_rvs.size() * fp.numFramesInC();
  


  int primeIndex = 0;     // which of the Cprime_rvs or Eprime_rvs to unpack to
  int originalIndex = 0;  // which of the C_rvs to print from

  int minAvailableFrame = -1; // nothing unpacked yet
  int maxAvailableFrame = -1;

  while (!frameRange_it->at_end()) {
    unsigned ppp = (*frameRange_it);
    infoMsg(IM::Printing,IM::High,"frame %u ", ppp);

    // map frame to original partition
    unsigned part;
    if ( (*frameRange_it) < (int)NP ) {
      part = 0; // P
    } else if ( (*frameRange_it) >= (int)NP + ((int)totalOriginalPartitions-2)*(int)NC ) {
      part = totalOriginalPartitions - 1; // E
    } else {
      part = 1 + ((*frameRange_it) - NP) / NC; // C
    }
    infoMsg(IM::Printing,IM::High,"in original partition %u ", part);

    // already fully unpacked? if so, print it
    if (minAvailableFrame <= (*frameRange_it) && (*frameRange_it) <= maxAvailableFrame) {
      infoMsg(IM::Printing,IM::High,"is available to print:\n");
      if (part == 0) { // print P partition
	storeToObsFile((*frameRange_it), segment, P_rvs, pVitTrigger, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, preg, 'p');
      } else if (part == totalOriginalPartitions-1) { // print E partition
	int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
	shiftOriginalVarstoPosition(E_rvs, targetFrame, Epos);
	storeToObsFile((*frameRange_it), segment, E_rvs, eVitTrigger, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, ereg, 'e');
      } else {      // print C partition
	int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
	originalIndex = ((int)part - 1) % (int) C_rvs.size();
	shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
	storeToObsFile((*frameRange_it), segment, C_rvs[originalIndex], cVitTrigger, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, creg, 'c');
      }
      (*frameRange_it)++;  // move on to next frame
      continue;
    }

    /* Before we can print this original partition $C_j$ (j = part), we must 
       unpack the modified partitions 

       $$\left\{ C'_i \left| \, \max\left(-1,\left\lceil\frac{j-s-m+1}{s}\right\rceil\right) \leq i 
         \leq \left\lfloor \frac{j}{s} \right\rfloor \right. \right\}$$

       where $s$ and $m$ are the boundary algorithm parameters, $C'_{-1}=P'$, and $C'_{N_{C'}}=E'$. 
     */
    
    int numerator = ( (int)part - 1 - (int)S - (int)M + 1 );
    int firstPrimePart;
    if (numerator <= -(int)S) {
      firstPrimePart = -1; // max
    } else if (numerator <= 0){
      firstPrimePart = 0;  // ceil
    } else {
      firstPrimePart = numerator / (int)S;
      if (numerator % (int)S)
	firstPrimePart += 1;   // ceil
    }
    firstPrimePart += 1; // account for C_{-1} = P'

    unsigned lastPrimePart = (part > 0) ?  1 + ((int)part-1) / (int)S : 0;
    if (lastPrimePart >= inference_it.pt_len())  // E original partition # may > # of modified partitions
      lastPrimePart = inference_it.pt_len() - 1;

    infoMsg(IM::Printing,IM::High,"requires unpacking modified partitions %u to %u:\n  unpack:", 
	    firstPrimePart, lastPrimePart);

    for (unsigned i = (unsigned) firstPrimePart; i <= lastPrimePart; i += 1) { // unpack C'_{i} set
      infoMsg(IM::Printing,IM::High,"  %u'", i); 
      setCurrentInferenceShiftTo(i);
      PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

      if (binaryViterbiFile) { // load packed values from disk if not already in memory
	readBinaryVitPartition(ps, i);
      }
      
      // unpack
      if (inference_it.at_p()) { // P'
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(P_partition_values.ptr,PprimeValuePtrs.ptr);
      } else if (inference_it.at_e()) { // E'
	// -1 to get the preceding C', -1 to account for C'_{-1} = P'
	primeIndex = ((int)i - 2) % (int)Eprime_rvs.size(); 
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(E_partition_values.ptr,EprimeValuePtrs[primeIndex].ptr);
      } else { // C'
	assert ( inference_it.at_c() );
	primeIndex = ((int)i - 1) % (int)Cprime_rvs.size();
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(C_partition_values.ptr  + 
			   ( binaryViterbiFile ? 0 : (inference_it.pt_i()-1)*ps.packer.packedLen() ),
			   CprimeValuePtrs[primeIndex].ptr);
      }
    }

    // unpacking firstPrimePart ... lastPrimePart makes these frames available:
    if (firstPrimePart == 0) {
      minAvailableFrame = 0;
    } else {
      minAvailableFrame = NP + (  ( (firstPrimePart-1) * S + M ) * NC  );
    }

    if (lastPrimePart == inference_it.pt_len()-1) {
      maxAvailableFrame = numUsableFrames - 1;
    } else {
      maxAvailableFrame = NP + lastPrimePart * S * NC - 1;
    }
    infoMsg(IM::Printing,IM::High,"  available frames %u to %u\n", minAvailableFrame, maxAvailableFrame);
  }

  if (vitObsFile) vitObsFile->endOfSegment();

  delete frameRange;
}


void
computeVarOrder(vector<RV *> &sectionRVs, regex_t *preg, char sectionLabel, int frame, unsigned segment, vector<string> &names) {
  assert(0 <= frame);
  set<string> nameSet;
  // collect the names selected by preg
  for (unsigned i=0; i < sectionRVs.size(); i+=1) {
    string name = sectionRVs[i]->name();
    if ( sectionRVs[i]->frame() == (unsigned) frame && (!preg || !regexec(preg, name.c_str(), 0,0,0)) ) {
      if (!sectionRVs[i]->discrete()) {
	string byStr("by -");
	byStr += sectionLabel;
	byStr += "VitRegexFilter";
	string sectStr("in ");
	sectStr += sectionLabel;
	sectStr +=" section";
	error("ERROR: variable '%s(%u)' in segment %u selected for Viterbi output to observation file %s is not discrete. "
	      "Only discrete variables may be stored to Viterbi observation files. "
	      "Use -(p|c|e)VitRegexFilter to select appropriate variables for Viterbi output.\n", 
	      name.c_str(), sectionRVs[i]->frame(),  segment, preg ? byStr.c_str() : sectStr.c_str() ); 
      } else {
	nameSet.insert(name);
      }
    }
  }
  names.resize(0); // empty it
  // now add the selected names to the vector in output order
  for (unsigned i=0; i < sectionRVs.size(); i+=1) {
    if (names.size() == nameSet.size()) {
      return; // got all the allowed names
    }
    string name = sectionRVs[i]->name();
    if (nameSet.find(name) != nameSet.end()) {
      names.push_back(name);
    }
  }
  assert(names.size() == nameSet.size());
}


void 
JunctionTree::storeToObsFile(int frame, unsigned segment,
			     vector<RV *> &rvs, 
			     bool useVitTrigger,
			     RVVec  &vitTriggerVec, 
			     string &vitTriggerExpr, 
			     RngDecisionTree::EquationClass &vitTriggerEqn,
			     regex_t *reg, char sectionLabel) 
{
  assert(0 <= frame);
  bool trigger = true;
  if (useVitTrigger) 
    trigger = evaluateTrigger(rvs, vitTriggerVec, vitTriggerExpr, vitTriggerEqn);
  if (!trigger || rvs.size()== 0) return;

  // check to see if we need to instantiate the output observation file
  if (vitObsFile == NULL) {
    // we need to get the names of the variables to output in order
    computeVarOrder(rvs, reg, sectionLabel, frame, segment, vitObsVariableNames);
    // we need to inform the user of the variable order in the output file
    printf("Viterbi values will be stored in the observation file in the order:");
    for (unsigned i=0; i < vitObsVariableNames.size(); i+=1) {
      printf(" %s", vitObsVariableNames[i].c_str());
    }
    printf("\n");
    // Now we can instantiate the file. This is Viterbi output, so there are no continuous features.
    // The number of discrete features is the # of variables to output
    vitObsFile = instantiateWriteFile(vitObsListName, vitObsFileName, const_cast<char *>(vitObsNameSeparator), 
				      const_cast<char *>(vitObsFileFmt), 0, vitObsVariableNames.size(), vitObsFileSwap);
  }

  set<string> nameSet;
  for (unsigned i=0; i < rvs.size(); i+=1) {
    unsigned f = rvs[i]->frame();
    assert(f <= 2147483647);
    if (f == (unsigned)frame) {
      if (!reg || !regexec(reg,rvs[i]->name().c_str(),0,0,0)) {
	if (std::find(vitObsVariableNames.begin(), vitObsVariableNames.end(), rvs[i]->name()) == vitObsVariableNames.end()) {
	  error("ERROR: extra RV %s selected for output in section %c\n", rvs[i]->name().c_str(), sectionLabel);
	}
	nameSet.insert(rvs[i]->name());
      }
    }
  }
  if (nameSet.size() != vitObsVariableNames.size()) {
    string missingVars("");
    for (vector<string>::iterator it=vitObsVariableNames.begin(); it != vitObsVariableNames.end(); ++it) {
      if (nameSet.find(*it) == nameSet.end()) {
	missingVars.append(*it);
	missingVars.append(" ");
      }
    }
    error("ERROR: missing RVs %sin %c\n", missingVars.c_str(), sectionLabel);
  }
  // actual output. number written must be a multiple of vitObsVariableNames.size()
  for (unsigned i=0; i < vitObsVariableNames.size(); i+=1) {
    unsigned j;
    for (j=0; vitObsVariableNames[i].compare( rvs[j]->name() ); j+=1)
      ;
    vitObsFile->writeFeature(  (Data32) ( dynamic_cast<DiscRV *>(rvs[j])->val )  );
  }
}



void
JunctionTree::printModifiedSection(PartitionStructures &ps,
				   unsigned *packed_values,
				   bool useVitTrigger,
				   RVVec &vitTriggerVec,
				   string &vitTriggerExpr,
				   RngDecisionTree::EquationClass &vitTriggerEqn,
				   bool printObserved,
				   unsigned part,
				   char sectionLabel,
				   FILE *f,
				   regex_t *preg,
				   vector<bool> &regex_mask,
				   bool &first_C,
				   unsigned &C_size,
				   sArray<unsigned> &previous_values,
				   bool runLengthCompress,
				   unsigned pt_i)
{
  if (ps.packer.packedLen() > 0)
    ps.packer.unpack(packed_values
		     + 
		     (pt_i-1)*ps.packer.packedLen(),
		     ps.hrvValuePtrs.ptr);
  printUnpackedSection(ps, useVitTrigger, vitTriggerVec, vitTriggerExpr, vitTriggerEqn, printObserved, part,
		       sectionLabel, f, preg, regex_mask, first_C, C_size, previous_values, runLengthCompress, pt_i);
}


void
JunctionTree::printUnpackedSection(PartitionStructures &ps,
				   bool useVitTrigger,
				   RVVec &vitTriggerVec,
				   string &vitTriggerExpr,
				   RngDecisionTree::EquationClass &vitTriggerEqn,
				   bool printObserved,
				   unsigned part,
				   char sectionLabel,
				   FILE *f,
				   regex_t *preg,
				   vector<bool> &regex_mask,
				   bool &first_C,
				   unsigned &C_size,
				   sArray<unsigned> &previous_values,
				   bool runLengthCompress,
				   unsigned pt_i)
{
  bool trigger = true;
  if (useVitTrigger) 
    trigger = evaluateTrigger(ps.allrvs_vec, vitTriggerVec, vitTriggerExpr, vitTriggerEqn);
  if (runLengthCompress)
    trigger = trigger && newViterbiValues(first_C, C_size, printObserved, ps.allrvs_vec, 
					  ps.hidRVVector, previous_values, regex_mask, preg);
  if (trigger && printObserved && ps.allrvs.size() > 0) {
    fprintf(f,"Ptn-%d %c': ",part, sectionLabel);
    printRVSetAndValues(f,ps.allrvs,true,preg);
  } else if (trigger && !printObserved && ps.packer.packedLen() > 0) {
    fprintf(f,"Ptn-%d %c': ",part, sectionLabel);
    printRVSetAndValues(f,ps.hidRVVector,true,preg);
  }
}

void
JunctionTree::printOriginalSection(vector<RV *> sectionRVs,
				   vector<RV *> hiddenRVs,
				   bool useVitTrigger,
				   RVVec &vitTriggerVec,
				   string &vitTriggerExpr,
				   RngDecisionTree::EquationClass &vitTriggerEqn,
				   bool printObserved,
				   unsigned part,
				   char sectionLabel,
				   FILE *f,
				   regex_t *preg,
				   vector<bool> &regex_mask,
				   bool &first_C,
				   unsigned &C_size,
				   sArray<unsigned> &previous_values,
				   bool runLengthCompress,
				   int frame)
{
  bool trigger = true;
  if (useVitTrigger) 
    trigger = evaluateTrigger(sectionRVs, vitTriggerVec, vitTriggerExpr, vitTriggerEqn);
  if (runLengthCompress)
    trigger = trigger && newViterbiValues(first_C, C_size, printObserved, sectionRVs, 
					  hiddenRVs, previous_values, regex_mask, preg);
  if (trigger && printObserved && sectionRVs.size() > 0) {
    setObservedRVs(sectionRVs);
    fprintf(f,"Ptn-%d %c: ",part, sectionLabel);
    if (frame > -1) {
      printRVSetAndValues(f,sectionRVs,true,preg, frame);
    } else {
      printRVSetAndValues(f,sectionRVs,true,preg);
    }
  } else if (trigger && !printObserved && hiddenRVs.size() > 0) {
    fprintf(f,"Ptn-%d %c: ",part, sectionLabel);
    if (frame > -1) {
      printRVSetAndValues(f,hiddenRVs,true,preg, frame);
    } else {
      printRVSetAndValues(f,hiddenRVs,true,preg);
    }
  }
}



/*
 *
 * This routine prints the Viterbi values computed by the most recent
 * linear inference run (assuming its data structures are still valid)
 * in ASCII to f (typically stdout). The Viterbi values are printed
 * by modified section (P', C', E').
 *
 * Preconditions: 
 *
 *    Assumes that distributeEvidence has just been run and all data
 *    structures (such as the compressed viterbi value array) are set
 *    up appropriately. 
 *
 *    Assumes that inference_it is currently set for the current
 *    segment.
 *  
 *    Assumes that the CC and CE partition pair random variables
 *    have been properly set up.
 * 
 *
 */

void
JunctionTree::printSavedPartitionViterbiValues(FILE* f,
					       bool printObserved,
					       regex_t* preg,
					       regex_t* creg,
					       regex_t* ereg,
					       char* partRangeFilter)
{
  printSavedPartitionViterbiValues(0, NULL, f, printObserved, preg, creg, ereg, partRangeFilter);
}



/*
 * This version of the above reads the saved binary Viterbi values
 * from vitFile to populate the inference data structures rather
 * than requiring the execution of distributeEvidence().
 */

void
JunctionTree::printSavedPartitionViterbiValues(unsigned numFrames,
					       FILE* vitFile,
					       FILE* f,
					       bool printObserved,
					       regex_t* preg,
					       regex_t* creg,
					       regex_t* ereg,
					       char* partRangeFilter)
{
  fprintf(f,"Printing random variables from (P',C',E')=(%d,%d,%d) modified sections\n",
	  P_partition_values.size(),
	  C_partition_values.size(),
	  E_partition_values.size());

  if (vitFile) {
    unsigned totalNumberPartitions;
    (void) unroll(numFrames,ZeroTable,&totalNumberPartitions);
    new (&inference_it) ptps_iterator(*this,totalNumberPartitions);
    init_CC_CE_rvs(inference_it);
  }

  Range* partRange = NULL;
  if (partRangeFilter != NULL) {
    partRange = new Range(partRangeFilter,0,inference_it.pt_len());
    if (partRange->length() == 0) { 
      warning("WARNING: Part range filter must specify a valid non-zero "
	      "length range within [0:%d]. Range given is %s\n",
	      inference_it.pt_len(), partRangeFilter);
      delete partRange;
      partRange = NULL;
    }
  }
  if (partRange == NULL)
    partRange = new Range("all",0,inference_it.pt_len());

  Range::iterator* partRange_it = new Range::iterator(partRange->begin());

  set<string> variableNames; // names of variables in the model
  for (unsigned i=0; i < partition_unrolled_rvs.size(); i+=1) {
    variableNames.insert(partition_unrolled_rvs[i]->name());
  }

  RVVec  pVitTriggerVec;
  string pVitTriggerExpr;
  RngDecisionTree::EquationClass pTriggerEqn;
  initializeViterbiTrigger(pVitTrigger, variableNames, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, 'p');

  RVVec  cVitTriggerVec;
  string cVitTriggerExpr;
  RngDecisionTree::EquationClass cTriggerEqn;
  initializeViterbiTrigger(cVitTrigger, variableNames, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, 'c');

  RVVec  eVitTriggerVec;
  string eVitTriggerExpr;
  RngDecisionTree::EquationClass eTriggerEqn;
  initializeViterbiTrigger(eVitTrigger, variableNames, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, 'e');

  bool first_P = true;
  bool first_C = true;
  bool first_E = true;
  unsigned P_size = 0;
  unsigned C_size = 0;
  unsigned E_size = 0;
  sArray<unsigned> previous_P_values;
  sArray<unsigned> previous_C_values;
  sArray<unsigned> previous_E_values;
  vector<bool> pregex_mask;
  vector<bool> cregex_mask;
  vector<bool> eregex_mask;

  while (!partRange_it->at_end()) {
    unsigned part = (*partRange_it);
    setCurrentInferenceShiftTo(part);
    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

    if (vitFile) {
      readBinaryVitPartition(ps, part);
    }
    if (inference_it.at_p()) {
      printModifiedSection(ps, P_partition_values.ptr, pVitTrigger!=NULL,
			   pVitTriggerVec, pVitTriggerExpr, pTriggerEqn,
			   printObserved, part, 'P', f, preg, pregex_mask,
			   first_P, P_size, previous_P_values);
    } else if (inference_it.at_e()) {
      printModifiedSection(ps, E_partition_values.ptr, eVitTrigger!=NULL,
			   eVitTriggerVec, eVitTriggerExpr, eTriggerEqn,
			   printObserved, part, 'E', f, ereg, eregex_mask,
			   first_E, E_size, previous_E_values);
    } else {
      assert ( inference_it.at_c() );      
      printModifiedSection(ps, C_partition_values.ptr, cVitTrigger!=NULL,
			   cVitTriggerVec, cVitTriggerExpr, cTriggerEqn,
			   printObserved, part, 'C', f, creg, cregex_mask,
			   first_C, C_size, previous_C_values,
			   vitRunLength, 
			   vitFile ? 1 : inference_it.pt_i());
    }
    (*partRange_it)++;
  }
  delete partRange;
  //clearAfterUnroll();  ???
}



/* FTOC maps frame # f to chunk # given a prologue of P frames and chunks of C frames */
#define FTOC(P,C,f) \
  ( ((f) - (P)) / (C) )

/*
 * Create the printing (P,C,E _rvs) and unpacking (Pprime,Cprime,Eprime _rvs) buffers for printSavedViterbiValues
 */

void JunctionTree::createUnpackingMap(
  vector<RV*> &unrolled_rvs, map<RVInfo::rvParent, unsigned> &unrolled_map,
  vector<RV*> &P_rvs, vector<RV*> &hidP_rvs, vector<RV*> &Pprime_rvs, vector<RV *> &hidPprime_rvs,
  vector<vector<RV*> > &C_rvs, vector<vector<RV*> > &hidC_rvs, vector<vector<RV*> > &Cprime_rvs, vector<vector<RV*> > &hidCprime_rvs,
  vector<RV*> &E_rvs, vector<RV *> &hidE_rvs, vector<vector<RV*> > &Eprime_rvs, vector<vector<RV*> > &hidEprime_rvs,
  sArray<DiscRVType*> &PprimeValuePtrs, 
  vector<sArray<DiscRVType *> > &CprimeValuePtrs, 
  vector<sArray<DiscRVType *> > &EprimeValuePtrs) 
{
  unsigned M = gm_template.M;
  unsigned S = gm_template.S;
  infoMsg(IM::Printing, IM::Info,"creating C' -> C map, (M,S) = (%u,%u)    frames in (P,C,E) = (%d,%d,%d)\n", 
	  M, S, fp.numFramesInP(), fp.numFramesInC(), fp.numFramesInE());

  // the # of modified C' partitions needed  
  unsigned nCprimes = 1 + (M / S) + ( (M % S) ? 1 : 0 );

  // the # of original C partitions needed
  unsigned nCs = nCprimes * S;
  infoMsg(IM::Printing, IM::Moderate, "# C' = %u   # C = %u\n", nCprimes, nCs);

  // unroll to create RV instances shared by the printing & unpacking sets
  fp.unroll(nCs-1, unrolled_rvs, unrolled_map);

  unsigned NP = fp.numFramesInP();
  unsigned NC = fp.numFramesInC();

  P_rvs.clear(); hidP_rvs.clear();
  Pprime_rvs.clear(); hidPprime_rvs.clear();

  C_rvs.resize(nCs); hidC_rvs.resize(nCs);
  vector<set<RV*> > Calready(nCs);
  for (unsigned i=0; i < nCs; i+=1) {
    C_rvs[i].clear();
    hidC_rvs[i].clear();
    Calready[i].clear();
  }
  Cprime_rvs.resize(nCprimes); hidCprime_rvs.resize(nCprimes); 
  for (unsigned i=0; i < nCprimes; i+=1) {
    Cprime_rvs[i].clear();
    hidCprime_rvs[i].clear();
  }
  E_rvs.clear(); hidE_rvs.clear();
  Eprime_rvs.resize(nCprimes); hidEprime_rvs.resize(nCprimes);
  for (unsigned i=0; i < nCprimes; i+=1) {
    Eprime_rvs[i].clear();
    hidEprime_rvs[i].clear();
  }
  CprimeValuePtrs.resize(nCprimes);
  EprimeValuePtrs.resize(nCprimes);

  // the unpacking could end at any of the C' sets, which
  // in turn means that we don't know which of the C sets
  // unpacking E' will complete. So, we setup an E' set for
  // each possible case, and use the sentence length to 
  // determine which to use at unpacking time

  // TODO: Since the mapping is now based on partitionStructureArray,
  //       the number of C's in the unrolled model is known at map
  //       creation time. Thus we should be able to construct only
  //       the E' mapping for only the relevant C'E' transition 
  //       (Might not apply for gmtkOnline though!)

  infoMsg(IM::Printing, IM::Moderate, "\nP':\n");
  vector<RV*> P = partitionStructureArray[0].allrvs_vec;
  for (vector<RV*>::iterator it = P.begin(); it != P.end(); ++it) {
    RV *v = getRV(unrolled_rvs, unrolled_map, *it);
    if ( (*it)->frame() < NP ) { // v is in P and P'
      P_rvs.push_back(v);
      Pprime_rvs.push_back(v);
      if (v->hidden()) {
	hidP_rvs.push_back(v);
	hidPprime_rvs.push_back(v);
      }
      infoMsg(IM::Printing, IM::Moderate, "%s(%u) -> P\n", (*it)->name().c_str(), (*it)->frame());
    } else {                     // v is in P' but not P
      Pprime_rvs.push_back(v);
      unsigned t = FTOC(NP,NC,(*it)->frame()); // the unmodified C index this variable belongs in
      C_rvs[t].push_back(v);
      Calready[t].insert(v);
      if (v->hidden()) {
	hidPprime_rvs.push_back(v);
	hidC_rvs[t].push_back(v);
      }
      infoMsg(IM::Printing, IM::Moderate, "%s(%u) -> C(%u)\n", (*it)->name().c_str(), (*it)->frame(), t);
    }
  }
  PprimeValuePtrs.resize(hidPprime_rvs.size());
  for (unsigned i=0; i < hidPprime_rvs.size(); i+=1) {
    DiscRV *drv = (DiscRV *) hidPprime_rvs[i];
    PprimeValuePtrs[i] = &(drv->val);
  }

  // If the unrolling consists of just P'E' then partitionStructureArray does
  // not contain a C' to base the C' -> C mapping on. But since there are no
  // C's to unpack, we don't need to construct the C' -> C mapping. Ticket #127
  if (partitionStructureArray.size() > 2) {
    // Map the first C' variables to their corresponding C (for printing)
    // and C' (for unpacking) instances
    infoMsg(IM::Printing, IM::Moderate, "\nC':\n");
    vector<RV*> C = partitionStructureArray[1].allrvs_vec;
    for (vector<RV*>::iterator it = C.begin(); it != C.end(); ++it) {
      RV *v = getRV(unrolled_rvs, unrolled_map, *it);
      unsigned t = FTOC(NP,NC,(*it)->frame()); // the unmodified C index this variable belongs in
      C_rvs[t].push_back(v);
      Calready[t].insert(v);
      Cprime_rvs[0].push_back(v);
      if (v->hidden()) {
	hidC_rvs[t].push_back(v);
	hidCprime_rvs[0].push_back(v);
      }
      infoMsg(IM::Printing, IM::Moderate, "C(%u) C'[0] : %s(%u)\n", FTOC(NP,NC,(*it)->frame()), (*it)->name().c_str(), (*it)->frame());
    }

    // Slide each C'[0] variable forward S Cs to compute the
    // remaining C's
    for (unsigned i=1; i < nCprimes; i+=1) {
      infoMsg(IM::Printing, IM::Moderate, "\n");
      for (vector<RV*>::iterator it = Cprime_rvs[i-1].begin(); 
	   it != Cprime_rvs[i-1].end();
	   ++it)
	{
	  RV *v = *it;
	  unsigned f = NP + ( v->frame() - NP + S * NC ) % (NC * nCs); // shift v over S original Cs (mod nCs)
	  RVInfo::rvParent target(v->name(), f);
	  RV *rv = getRV(unrolled_rvs, unrolled_map, target);
	  unsigned t = FTOC(NP,NC,f);

	  // rv may already be in (hid?)C_rvs[t], resulting in excess length when printing C.
	  // We want the ordering of a vector, but the single membership of a set...
	  // Does need to be added to (hid?)Cprime_rvs[i], as that is newly constructed...
	  if (Calready[t].find(rv) == Calready[t].end()) {
	    C_rvs[t].push_back(rv);
	    if (rv->hidden()) {
	      hidC_rvs[t].push_back(rv);
	    }
	    Calready[t].insert(rv);
	  }
	  Cprime_rvs[i].push_back(rv);
	  if (rv->hidden()) {
	    hidCprime_rvs[i].push_back(rv);
	  }
	  infoMsg(IM::Printing, IM::Moderate, "C(%u) C'[%u] : %s(%u)\n", t, i, target.first.c_str(), target.second);
	}
    }
    // sort C_rvs[i] to ensure the variables end up in the same order as P_rvs and E_rvs
    for (unsigned i=0; i < nCprimes; i+=1) {
      sort(C_rvs[i].begin(), C_rvs[i].end(), rvcompare);
    }
    for (unsigned i=0; i < nCprimes; i+=1) {
      CprimeValuePtrs[i].resize(hidCprime_rvs[i].size());
      for (unsigned j=0; j < hidCprime_rvs[i].size(); j+=1) {
	DiscRV *drv = (DiscRV *) hidCprime_rvs[i][j];
	CprimeValuePtrs[i][j] = &(drv->val);
      }
    }
  }

  infoMsg(IM::Printing, IM::Moderate, "\nE':\n");
  int Eidx  = partitionStructureArray.size() - 1;
  int delta = (nCs - ( (partitionStructureArray.size() - 2) * S + M ) ) * NC;
  unsigned firstEframe = NP + nCs * NC;
  infoMsg(IM::Printing, IM::Moderate, "delta = %d    firstE = %u\n", delta, firstEframe);
  vector<RV*> E = partitionStructureArray[Eidx].allrvs_vec;
  for (vector<RV*>::iterator it = E.begin(); it != E.end(); ++it) {
    RV *v = *it;
    if (v->frame() + delta >= firstEframe) { // v is in E and E'
      RV *rv = getRV(unrolled_rvs, unrolled_map, v, delta);
      E_rvs.push_back(rv);
      if (rv->hidden()) hidE_rvs.push_back(rv);
      for (unsigned i=0; i < nCprimes; i+=1) {
	Eprime_rvs[i].push_back(rv);
	if (rv->hidden()) hidEprime_rvs[i].push_back(rv);
      }
      infoMsg(IM::Printing, IM::Moderate, "%s(%u) -> %s(%u) -> E\n", v->name().c_str(), v->frame()+delta, rv->name().c_str(), rv->frame());
    } else { // v is in E' but not E
      // The C'E' transition could occur after unpacking any of the nCprimes C' sets,
      // so Eprime_rvs[i] handles the C'[i]E' transition
      for (unsigned i=0; i < nCprimes; i+=1) {
	unsigned t = NP + (v->frame() + (S * NC * (1 - (Eidx-1))) - NP + S * NC * i) % (NC * nCs);
	unsigned originalC = FTOC(NP,NC,t);
	RVInfo::rvParent target(v->name(), t);
	RV *rv_shifted = getRV(unrolled_rvs, unrolled_map, target);
	Eprime_rvs[i].push_back(rv_shifted);
	if (rv_shifted->hidden()) 
	  hidEprime_rvs[i].push_back(rv_shifted);
	if (Calready[originalC].find(rv_shifted) == Calready[originalC].end()) {
	  C_rvs[originalC].push_back(rv_shifted);
	  if (rv_shifted->hidden()) 
	    hidC_rvs[originalC].push_back(rv_shifted);
	  Calready[originalC].insert(rv_shifted);
	  infoMsg(IM::Printing, IM::High, "Eprime_rvs[%u] += %s(%u) -> %s(%u) into C_rvs[%u]\n", 
		  i, v->name().c_str(), v->frame(), rv_shifted->name().c_str(), rv_shifted->frame(), originalC);
	} else 
	  infoMsg(IM::Printing, IM::High, "Eprime_rvs[%u] += %s(%u) -> %s(%u) already in C_rvs[%u]\n", 
		  i, v->name().c_str(), v->frame(), rv_shifted->name().c_str(), rv_shifted->frame(), originalC);
      }
      infoMsg(IM::Printing, IM::Moderate, "%s(%u) -> C(%u)\n", v->name().c_str(), v->frame()+delta, FTOC(NP,NC,v->frame()+delta));
    }
  }
  for (unsigned i=0; i < nCprimes; i+=1) {
    EprimeValuePtrs[i].resize(hidEprime_rvs[i].size());
    for (unsigned j=0; j < hidEprime_rvs[i].size(); j+=1) {
      DiscRV *drv = (DiscRV *) hidEprime_rvs[i][j];
      EprimeValuePtrs[i][j] = &(drv->val);
    }
  }

  // print mapping as debug info
  if (IM::messageGlb(IM::Printing, IM::Info)) {
    infoMsg(IM::Printing,IM::Info, "G' -> G mapping:\n\n");

    set<RV*> pset(P_rvs.begin(), P_rvs.end());
    for (vector<RV*>::iterator it = Pprime_rvs.begin(); 
	 it != Pprime_rvs.end();
	 ++it)
      {
	RV *v = *it;
	if (pset.find(v) == pset.end()) {
	  infoMsg(IM::Printing,IM::Info, "P'    C[%u]: %s(%u)\n", FTOC(NP,NC,v->frame()), v->name().c_str(), v->frame());
	} else {
	  infoMsg(IM::Printing,IM::Info, "P'    P   : %s(%u)\n", v->name().c_str(), v->frame());
	}
      }
    for (unsigned i=0; i < nCprimes; i+=1) {
      infoMsg(IM::Printing, IM::Info, "-----------------------------\n");
      for (vector<RV*>::iterator it = Cprime_rvs[i].begin();
	   it != Cprime_rvs[i].end();
	   ++it)
	{
	  RV *v = *it;
	  infoMsg(IM::Printing,IM::Info, "C'[%u] C[%u]: %s(%u)\n", i, FTOC(NP,NC,v->frame()),v->name().c_str(), v->frame());
	}
    }

    for (unsigned i=0; i < nCprimes; i+=1) {
      infoMsg(IM::Printing, IM::Info, "-----------------------------\n");
      for (vector<RV*>::iterator it = Eprime_rvs[i].begin(); 
	   it != Eprime_rvs[i].end();
	   ++it)
	{
	  RV *v = *it;
	  unsigned c = FTOC(NP,NC,v->frame());
	  if (c < nCs) {
	    infoMsg(IM::Printing,IM::Info, "E'[%u] C[%u]: %s(%u)\n", i, c, v->name().c_str(), v->frame());
	  } else {
	    infoMsg(IM::Printing,IM::Info, "E'[%u] E   : %s(%u)\n", i, v->name().c_str(), v->frame());
	  }
	}
    }
  }
}



/*
 *
 * This routine prints the Viterbi values computed by the most recent
 * linear inference run (assuming its data structures are still valid)
 * in ASCII to f (typically stdout). Unlike printSavedPartitionViterbiValues(), 
 * this method prints the values ordered by the original P, C, and E,
 * sections rather than the modified P', C', and E' sections.
 *
 * Preconditions: 
 *
 *    Assumes that distributeEvidence has just been run and all data
 *    structures (such as the compressed viterbi value array) are set
 *    up appropriately. 
 *
 *    Assumes that inference_it is currently set for the current
 *    segment.
 *  
 *    Assumes that the CC and CE partition pair random variables
 *    have been properly set up.
 * 
 *
 */
void
JunctionTree::printSavedViterbiValues(FILE* f,
				      bool printObserved,
				      regex_t* preg, regex_t* creg, regex_t* ereg)
{
  printSavedViterbiValues(0, f, NULL, printObserved, preg, creg, ereg);
}


/*
 * This version of the above reads the saved binary Viterbi values
 * from binVitFile to populate the Viterbi value data structures rather
 * than requiring the execution of distributeEvidence().
 */
void
JunctionTree::printSavedViterbiValues(unsigned numFrames,
				      FILE* f,
				      FILE* binVitFile,
				      bool printObserved,
				      regex_t* preg, regex_t* creg, regex_t* ereg)
{
  if (binVitFile) {
    unsigned totalNumberPartitions;
    (void) unroll(numFrames,ZeroTable,&totalNumberPartitions);
    new (&inference_it) ptps_iterator(*this,totalNumberPartitions);
    init_CC_CE_rvs(inference_it);
  }
  vector<RV*> unrolled_rvs;
  map<RVInfo::rvParent, unsigned> unrolled_map;

  vector<RV*> P_rvs;      // original P for printing
  vector<RV*> Pprime_rvs; // modified P' for unpacking
  vector<RV*> hidP_rvs;      // hidden subset of original P for printing
  vector<RV*> hidPprime_rvs; // hidden subset of modified P' for unpacking

  vector<vector<RV*> > C_rvs; // original Cs for printing
  vector<vector<RV*> > Cprime_rvs; // modified C's for unpacking
  vector<vector<RV*> > hidC_rvs; // hidden subset of original Cs for printing
  vector<vector<RV*> > hidCprime_rvs; // hidden subset of modified C's for unpacking

  vector<RV*> E_rvs; // ... printing
  vector<vector<RV*> > Eprime_rvs; // ... unpacking
  vector<RV*> hidE_rvs; // ... printing
  vector<vector<RV*> > hidEprime_rvs; // ... unpacking

  sArray<DiscRVType *>PprimeValuePtrs;
  vector<sArray<DiscRVType *> > CprimeValuePtrs;
  vector<sArray<DiscRVType *> > EprimeValuePtrs;

  createUnpackingMap(unrolled_rvs, unrolled_map, 
		     P_rvs, hidP_rvs, Pprime_rvs, hidPprime_rvs, 
		     C_rvs, hidC_rvs, Cprime_rvs, hidCprime_rvs,
		     E_rvs, hidE_rvs, Eprime_rvs, hidEprime_rvs,
		     PprimeValuePtrs, CprimeValuePtrs, EprimeValuePtrs);

  unsigned M = gm_template.M;
  unsigned S = gm_template.S;

  fprintf(f,"Printing random variables from (P,C,E)=(%d,%d,%d) sections\n",
	  P_partition_values.size(),
	  C_partition_values.size() * S + M,
	  E_partition_values.size());

  Range* partRange = new Range("all",0,inference_it.pt_len());

  Range::iterator* partRange_it = new Range::iterator(partRange->begin());

  vector<int> Cpos(C_rvs.size());
  for (unsigned int i=0; i < Cpos.size(); i+=1) 
    Cpos[i] = fp.numFramesInP() + i * fp.numFramesInC();
  int Epos = fp.numFramesInP() + C_rvs.size() * fp.numFramesInC();


  set<string> variableNames; // names of variables in the model
  for (unsigned i=0; i < partition_unrolled_rvs.size(); i+=1) {
    variableNames.insert(partition_unrolled_rvs[i]->name());
  }

  RVVec  pVitTriggerVec;
  string pVitTriggerExpr;
  RngDecisionTree::EquationClass pTriggerEqn;
  initializeViterbiTrigger(pVitTrigger, variableNames, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, 'p');

  RVVec  cVitTriggerVec;
  string cVitTriggerExpr;
  RngDecisionTree::EquationClass cTriggerEqn;
  initializeViterbiTrigger(cVitTrigger, variableNames, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, 'c');

  RVVec  eVitTriggerVec;
  string eVitTriggerExpr;
  RngDecisionTree::EquationClass eTriggerEqn;
  initializeViterbiTrigger(eVitTrigger, variableNames, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, 'e');

  bool first_P = true;
  bool first_C = true;
  bool first_E = true;
  unsigned P_size = 0;
  unsigned C_size = 0;
  unsigned E_size = 0;
  sArray<unsigned> previous_P_values;
  sArray<unsigned> previous_C_values;
  sArray<unsigned> previous_E_values;
  vector<bool> pregex_mask;
  vector<bool> cregex_mask;
  vector<bool> eregex_mask;

  unsigned primeIndex = 0;
  unsigned originalIndex = 0;
  unsigned Ccount = 1;

  while (!partRange_it->at_end()) {
    
    unsigned part = (*partRange_it);
    setCurrentInferenceShiftTo(part);
    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

    if (binVitFile) {
      unsigned N_best = 1;
      unsigned num_to_read = N_best * ps.packer.packedLen();
      if (inference_it.at_p()) {
	readVitIntVector(num_to_read, P_partition_values.ptr);
      } else if (inference_it.at_e()) {
	readVitIntVector(num_to_read, E_partition_values.ptr);
      } else {
	readVitIntVector(num_to_read, C_partition_values.ptr);
      }
    }
    if (inference_it.at_p()) {
      if (ps.packer.packedLen() > 0) 
        ps.packer.unpack(P_partition_values.ptr,PprimeValuePtrs.ptr);
      printOriginalSection(P_rvs, hidP_rvs, pVitTrigger != NULL,  pVitTriggerVec, 
			   pVitTriggerExpr, pTriggerEqn, printObserved, part, 'P', 
			   f, preg, pregex_mask, first_P, P_size, previous_P_values, false);
    } else if (inference_it.at_e()) {
      primeIndex = (primeIndex + Eprime_rvs.size() - 1) % Eprime_rvs.size(); // primeIndex -= 1 mod nCprimes
      if (ps.packer.packedLen() > 0) 
	ps.packer.unpack(E_partition_values.ptr,EprimeValuePtrs[primeIndex].ptr);
      // print completed C partitions
      int targetFrame = fp.numFramesInP() + (int)(part-1) * gm_template.S * fp.numFramesInC();
      for (unsigned i=0; i < gm_template.M; i+=1) { // unpacking E' completes the last M Cs
	shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
	printOriginalSection(C_rvs[originalIndex], hidC_rvs[originalIndex], cVitTrigger != NULL, cVitTriggerVec,
			     cVitTriggerExpr, cTriggerEqn, printObserved, part,'C', f, 
			     creg, cregex_mask, first_C, C_size, previous_C_values, vitRunLength);
	Ccount += 1;
	originalIndex = (originalIndex + 1) % C_rvs.size();
	targetFrame += fp.numFramesInC();
      } 
      // print E partition
      if ( (hidE_rvs.size() > 0)  || (printObserved && E_rvs.size() > 0) ) {
	shiftOriginalVarstoPosition(E_rvs, targetFrame, Epos);
	printOriginalSection(E_rvs, hidE_rvs, eVitTrigger != NULL,  eVitTriggerVec, 
			     eVitTriggerExpr, eTriggerEqn, printObserved, part, 'E', f, 
			     ereg, eregex_mask, first_E, E_size, previous_E_values, false);
      }
    } else {
      assert ( inference_it.at_c() );
      // print C partition
      {
	unsigned long offset = binVitFile ? 0 : (inference_it.pt_i()-1) * ps.packer.packedLen();
        if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(C_partition_values.ptr + offset, CprimeValuePtrs[primeIndex].ptr);
	int targetFrame = fp.numFramesInP() + (int)(part-1) * gm_template.S * fp.numFramesInC();
	for (unsigned i=0; i < gm_template.S; i+=1) { // unpacking a C' completes S Cs
	  shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
	  printOriginalSection(C_rvs[originalIndex], hidC_rvs[originalIndex], cVitTrigger != NULL, cVitTriggerVec,
			       cVitTriggerExpr, cTriggerEqn, printObserved, part, 'C', f, 
			       creg, cregex_mask, first_C, C_size, previous_C_values, vitRunLength);
	  Ccount += 1;
	  originalIndex = (originalIndex + 1) % C_rvs.size();
	  targetFrame += fp.numFramesInC();
	}
	primeIndex = (primeIndex + 1) % Cprime_rvs.size();
      }
      // previous_C = inference_it.pt_i();
    }
    (*partRange_it)++;
  }

  delete partRange;
  //clearAfterUnroll();
}


void
JunctionTree::readBinaryVitPartition(PartitionStructures& ps, unsigned part) {
  assert(part == inference_it.pt_i());
  unsigned N_best = 1;
  unsigned num_to_read = N_best * ps.packer.packedLen();
  off_t offset;
  if (inference_it.at_p()) {
    if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset, SEEK_SET) == (off_t) -1) {
      char *err = strerror(errno);
      error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
    }
    readVitIntVector(num_to_read, P_partition_values.ptr);
  } else if (inference_it.at_e()) {
    offset = (off_t)  // P' size + (T-2) * C' size
      (   (   N_best * partitionStructureArray[0].packer.packedLen()       
	    + N_best * partitionStructureArray[1].packer.packedLen() * inference_it.num_c_partitions()  
	  ) * sizeof(unsigned)   );                         
    if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset + offset, SEEK_SET) == (off_t) -1) {
      char *err = strerror(errno);
      error("ERROR: seek failed on '%s': %s\n", binaryViterbiFilename, err);
    }       
    readVitIntVector(num_to_read, E_partition_values.ptr);
  } else {
    offset = (off_t)  // P' size + (t-1) * C' size
      (   (   N_best * partitionStructureArray[0].packer.packedLen()
            + N_best * partitionStructureArray[1].packer.packedLen() * ( part - 1 )
	  ) * sizeof(unsigned)   );               
    if (gmtk_fseek(binaryViterbiFile, binaryViterbiOffset + offset, SEEK_SET) == (off_t) -1) {
      char *err = strerror(errno);
      error("seek failed on '%s': %s\n", binaryViterbiFilename, err);
    }       
    readVitIntVector(num_to_read, C_partition_values.ptr);
  } 
}


/*
 * This version prints the Viterbi values for the original (P, C, E)
 * sections specified by the partRangeFilter.
 *
 * Preconditions: 
 *
 *    If binVitFile is NULL, assumes that distributeEvidence() has just 
 *    been run and all data structures (such as the compressed Viterbi 
 *    value array) are set up appropriately. If binVitFile is non-NULL,
 *    it reads the saved binary Viterbi values from the file to populate 
 *    the Viterbi value data structures rather than requiring the execution
 *    of distributeEvidence().
 *
 *    Assumes that inference_it is currently set for the current
 *    segment.
 *  
 *    Assumes that the CC and CE partition pair random variables
 *    have been properly set up.
 * 
 */
void
JunctionTree::printSavedViterbiValues(unsigned numFrames, FILE* f,
				      FILE *binVitFile,
				      bool printObserved,
				      regex_t* preg, regex_t* creg, regex_t* ereg,
				      char *partRangeFilter)
{

  if (partRangeFilter == NULL) {
    if (binaryViterbiFile) 
      printSavedViterbiValues(numFrames, f, binaryViterbiFile, printObserved, preg, creg, ereg);
    else
      printSavedViterbiValues(f, printObserved, preg, creg, ereg);
    return;
  }

  if (binaryViterbiFile) {
    unsigned totalNumberPartitions;
    (void) unroll(numFrames,ZeroTable,&totalNumberPartitions);
    
    new (&inference_it) ptps_iterator(*this,totalNumberPartitions);
    init_CC_CE_rvs(inference_it);
  }

  vector<RV*> unrolled_rvs;
  map<RVInfo::rvParent, unsigned> unrolled_map;

  vector<RV*> P_rvs;      // original P for printing
  vector<RV*> Pprime_rvs; // modified P' for unpacking
  vector<RV*> hidP_rvs;      // hidden subset of original P for printing
  vector<RV*> hidPprime_rvs; // hidden subset of modified P' for unpacking

  vector<vector<RV*> > C_rvs; // original Cs for printing
  vector<vector<RV*> > Cprime_rvs; // modified C's for unpacking
  vector<vector<RV*> > hidC_rvs; // hidden subset of original Cs for printing
  vector<vector<RV*> > hidCprime_rvs; // hidden subset of modified C's for unpacking

  vector<RV*> E_rvs; // ... printing
  vector<vector<RV*> > Eprime_rvs; // ... unpacking
  vector<RV*> hidE_rvs; // ... printing
  vector<vector<RV*> > hidEprime_rvs; // ... unpacking

  sArray<DiscRVType *>PprimeValuePtrs;
  vector<sArray<DiscRVType *> > CprimeValuePtrs;
  vector<sArray<DiscRVType *> > EprimeValuePtrs;

  createUnpackingMap(unrolled_rvs, unrolled_map, 
		     P_rvs, hidP_rvs, Pprime_rvs, hidPprime_rvs, 
		     C_rvs, hidC_rvs, Cprime_rvs, hidCprime_rvs,
		     E_rvs, hidE_rvs, Eprime_rvs, hidEprime_rvs,
		     PprimeValuePtrs, CprimeValuePtrs, EprimeValuePtrs);
  unsigned M = gm_template.M;
  unsigned S = gm_template.S;

  fprintf(f,"Printing random variables from (P,C,E)=(%d,%d,%d) sections\n",
	  P_partition_values.size(),
	  C_partition_values.size() * S + M,
	  E_partition_values.size());

  unsigned totalOriginalPartitions = 2 + inference_it.num_c_partitions() * S + M;

  infoMsg(IM::Printing,IM::High,"M = %u   S = %u   # orig parts = %u\n",
       M, S, totalOriginalPartitions);

  Range* partRange = NULL;
  partRange = new Range(partRangeFilter,0,totalOriginalPartitions);
  if (partRange->length() == 0) { 
    warning("WARNING: Part range filter must specify a valid non-zero "
	    "length range within [0:%d]. Range given is %s\n",
	    totalOriginalPartitions, partRangeFilter);
    delete partRange;
    partRange = NULL;
  }

  if (partRange == NULL)
    partRange = new Range("all",0,totalOriginalPartitions);

  Range::iterator* partRange_it = new Range::iterator(partRange->begin());

  vector<int> Cpos(C_rvs.size());
  for (unsigned int i=0; i < Cpos.size(); i+=1) 
    Cpos[i] = fp.numFramesInP() + i * fp.numFramesInC();
  int Epos = fp.numFramesInP() + C_rvs.size() * fp.numFramesInC();
  

  set<string> variableNames; // names of variables in the model
  for (unsigned i=0; i < partition_unrolled_rvs.size(); i+=1) {
    variableNames.insert(partition_unrolled_rvs[i]->name());
  }

  RVVec  pVitTriggerVec;
  string pVitTriggerExpr;
  RngDecisionTree::EquationClass pTriggerEqn;
  initializeViterbiTrigger(pVitTrigger, variableNames, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, 'p');

  RVVec  cVitTriggerVec;
  string cVitTriggerExpr;
  RngDecisionTree::EquationClass cTriggerEqn;
  initializeViterbiTrigger(cVitTrigger, variableNames, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, 'c');

  RVVec  eVitTriggerVec;
  string eVitTriggerExpr;
  RngDecisionTree::EquationClass eTriggerEqn;
  initializeViterbiTrigger(eVitTrigger, variableNames, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, 'e');

  bool first_P = true;
  bool first_C = true;
  bool first_E = true;
  unsigned P_size = 0;
  unsigned C_size = 0;
  unsigned E_size = 0;
  sArray<unsigned> previous_P_values;
  sArray<unsigned> previous_C_values;
  sArray<unsigned> previous_E_values;
  vector<bool> pregex_mask;
  vector<bool> cregex_mask;
  vector<bool> eregex_mask;

  int primeIndex = 0;     // which of the Cprime_rvs or Eprime_rvs to unpack to
  int originalIndex = 0;  // which of the C_rvs to print from

  while (!partRange_it->at_end()) {

    unsigned part = (*partRange_it);

    /* Before we can print this original partition $C_j$ (j = part), we must 
       unpack the modified partitions 

       $$\left\{ C'_i \left| \, \max\left(-1,\left\lceil\frac{j-s-m+1}{s}\right\rceil\right) \leq i 
         \leq \left\lfloor \frac{j}{s} \right\rfloor \right. \right\}$$

       where $s$ and $m$ are the boundary algorithm parameters, $C'_{-1}=P'$, and $C'_{N_{C'}}=E'$. 
     */
    
    int numerator = ( (int)part - 1 - (int)S - (int)M + 1 );
    int firstPrimePart;
    if (numerator <= -(int)S) {
      firstPrimePart = -1; // max
    } else if (numerator <= 0){
      firstPrimePart = 0;  // ceil
    } else {
      firstPrimePart = numerator / (int)S;
      if (numerator % (int)S)
	firstPrimePart += 1;   // ceil
    }
    firstPrimePart += 1; // account for C_{-1} = P'

    unsigned lastPrimePart =  (part > 0) ?  1 + ((int)part-1) / (int)S : 0;
    if (lastPrimePart >= inference_it.pt_len())  // E original partition # may > # of modified partitions
      lastPrimePart = inference_it.pt_len() - 1;

    infoMsg(IM::Printing,IM::High,"original partition %u requires unpacking modified partitions %u to %u:\n", 
	    part, firstPrimePart, lastPrimePart);

    for (unsigned i = (unsigned) firstPrimePart; i <= lastPrimePart; i += 1) { // unpack C'_{i} set
      infoMsg(IM::Printing,IM::High,"unpack %u'\n", i); 
      setCurrentInferenceShiftTo(i);
      PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

      if (binaryViterbiFile) {
	readBinaryVitPartition(ps, i);
      }
      
      if (inference_it.at_p()) { // P'
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(P_partition_values.ptr,PprimeValuePtrs.ptr);
      } else if (inference_it.at_e()) { // E'
	// -1 to get the preceding C', -1 to account for C'_{-1} = P'
	primeIndex = ((int)i - 2) % (int)Eprime_rvs.size(); 
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(E_partition_values.ptr,EprimeValuePtrs[primeIndex].ptr);
      } else { // C'
	assert ( inference_it.at_c() );
	primeIndex = ((int)i - 1) % (int)Cprime_rvs.size();
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(C_partition_values.ptr  + 
			   ( binaryViterbiFile ? 0 : (inference_it.pt_i()-1)*ps.packer.packedLen() ),
			   CprimeValuePtrs[primeIndex].ptr);
      }
    }

    if (part == 0) { // print P partition
      printOriginalSection(P_rvs, hidP_rvs, pVitTrigger != NULL,  pVitTriggerVec, 
			   pVitTriggerExpr, pTriggerEqn, printObserved, part, 'P', f, 
			   preg, pregex_mask, first_P, P_size, previous_P_values, false);
    } else if (part == totalOriginalPartitions-1) { // print E partition
      if ( (hidE_rvs.size() > 0)  || (printObserved && E_rvs.size() > 0) ) {
	int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
	shiftOriginalVarstoPosition(E_rvs, targetFrame, Epos);
	printOriginalSection(E_rvs, hidE_rvs, eVitTrigger != NULL,  eVitTriggerVec, 
			     eVitTriggerExpr, eTriggerEqn, printObserved, part, 'E', f, 
			     ereg, eregex_mask, first_E, E_size, previous_E_values, false);
      }
    } else {      // print C partition
      int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
      originalIndex = ((int)part - 1) % (int) C_rvs.size();
      shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
      printOriginalSection(C_rvs[originalIndex], hidC_rvs[originalIndex], cVitTrigger != NULL, cVitTriggerVec,
			   cVitTriggerExpr, cTriggerEqn, printObserved, part, 'C', f, 
			   creg, cregex_mask, first_C, C_size, previous_C_values, vitRunLength);
    }
    (*partRange_it)++;
  }
  
  delete partRange;
}


/*
 * This version of the above prints the Viterbi values for the frames
 * specified by the frameRangeFilter.
 */

void
JunctionTree::printSavedViterbiFrames(unsigned numFrames, FILE* f,
				      FILE *binVitFile,
				      bool printObserved,
				      regex_t* preg, regex_t* creg, regex_t* ereg, 
				      char *frameRangeFilter)
{
  unsigned numUsableFrames;
  if (binaryViterbiFile) {
    unsigned totalNumberPartitions;
    numUsableFrames = unroll(numFrames,ZeroTable,&totalNumberPartitions);
    
    new (&inference_it) ptps_iterator(*this,totalNumberPartitions);
    init_CC_CE_rvs(inference_it);
  } else {
    numUsableFrames = this->numUsableFrames;
  }

  vector<RV*> unrolled_rvs;
  map<RVInfo::rvParent, unsigned> unrolled_map;

  vector<RV*> P_rvs;      // original P for printing
  vector<RV*> Pprime_rvs; // modified P' for unpacking
  vector<RV*> hidP_rvs;      // hidden subset of original P for printing
  vector<RV*> hidPprime_rvs; // hidden subset of modified P' for unpacking

  vector<vector<RV*> > C_rvs; // original Cs for printing
  vector<vector<RV*> > Cprime_rvs; // modified C's for unpacking
  vector<vector<RV*> > hidC_rvs; // hidden subset of original Cs for printing
  vector<vector<RV*> > hidCprime_rvs; // hidden subset of modified C's for unpacking

  vector<RV*> E_rvs; // ... printing
  vector<vector<RV*> > Eprime_rvs; // ... unpacking
  vector<RV*> hidE_rvs; // ... printing
  vector<vector<RV*> > hidEprime_rvs; // ... unpacking

  sArray<DiscRVType *>PprimeValuePtrs;
  vector<sArray<DiscRVType *> > CprimeValuePtrs;
  vector<sArray<DiscRVType *> > EprimeValuePtrs;

  createUnpackingMap(unrolled_rvs, unrolled_map, 
		     P_rvs, hidP_rvs, Pprime_rvs, hidPprime_rvs, 
		     C_rvs, hidC_rvs, Cprime_rvs, hidCprime_rvs,
		     E_rvs, hidE_rvs, Eprime_rvs, hidEprime_rvs,
		     PprimeValuePtrs, CprimeValuePtrs, EprimeValuePtrs);

  unsigned M = gm_template.M;
  unsigned S = gm_template.S;

  fprintf(f,"Printing random variables from (P,C,E)=(%d,%d,%d) sections\n",
	  P_partition_values.size(),
	  C_partition_values.size() * S + M,
	  E_partition_values.size());


  unsigned NP = fp.numFramesInP();
  unsigned NC = fp.numFramesInC();

  unsigned totalOriginalPartitions = 2 + inference_it.num_c_partitions() * S + M;

  infoMsg(IM::Printing,IM::High,"NP = %u   NC = %u   M = %u   S = %u   # orig parts = %u\n",
       NP, NC, M, S, totalOriginalPartitions);

  Range* frameRange = NULL;
  frameRange = new Range(frameRangeFilter,0,numUsableFrames);
  if (frameRange->length() == 0) { 
    warning("WARNING: Frame range filter must specify a valid non-zero "
	    "length range within [0:%d]. Range given is %s\n",
	    numUsableFrames, frameRangeFilter);
    delete frameRange;
    frameRange = NULL;
  }

  if (frameRange == NULL)
    frameRange = new Range("all",0,numUsableFrames);

  Range::iterator* frameRange_it = new Range::iterator(frameRange->begin());

  vector<int> Cpos(C_rvs.size());
  for (unsigned int i=0; i < Cpos.size(); i+=1) 
    Cpos[i] = fp.numFramesInP() + i * fp.numFramesInC();
  int Epos = fp.numFramesInP() + C_rvs.size() * fp.numFramesInC();
  

  set<string> variableNames; // names of variables in the model
  for (unsigned i=0; i < partition_unrolled_rvs.size(); i+=1) {
    variableNames.insert(partition_unrolled_rvs[i]->name());
  }

  RVVec  pVitTriggerVec;
  string pVitTriggerExpr;
  RngDecisionTree::EquationClass pTriggerEqn;
  initializeViterbiTrigger(pVitTrigger, variableNames, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, 'p');

  RVVec  cVitTriggerVec;
  string cVitTriggerExpr;
  RngDecisionTree::EquationClass cTriggerEqn;
  initializeViterbiTrigger(cVitTrigger, variableNames, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, 'c');

  RVVec  eVitTriggerVec;
  string eVitTriggerExpr;
  RngDecisionTree::EquationClass eTriggerEqn;
  initializeViterbiTrigger(eVitTrigger, variableNames, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, 'e');

  bool first_P = true;
  bool first_C = true;
  bool first_E = true;
  unsigned P_size = 0;
  unsigned C_size = 0;
  unsigned E_size = 0;
  sArray<unsigned> previous_P_values;
  sArray<unsigned> previous_C_values;
  sArray<unsigned> previous_E_values;
  vector<bool> pregex_mask;
  vector<bool> cregex_mask;
  vector<bool> eregex_mask;

  int primeIndex = 0;     // which of the Cprime_rvs or Eprime_rvs to unpack to
  int originalIndex = 0;  // which of the C_rvs to print from

  int minAvailableFrame = -1; // nothing unpacked yet
  int maxAvailableFrame = -1;

  while (!frameRange_it->at_end()) {
    unsigned ppp = (*frameRange_it);
    infoMsg(IM::Printing,IM::High,"frame %u ", ppp);

    // map frame to original partition
    unsigned part;
    if ( (*frameRange_it) < (int)NP ) {
      part = 0; // P
    } else if ( (*frameRange_it) >= (int)NP + ((int)totalOriginalPartitions-2)*(int)NC ) {
      part = totalOriginalPartitions - 1; // E
    } else {
      part = 1 + ((*frameRange_it) - NP) / NC; // C
    }
    infoMsg(IM::Printing,IM::High,"in original partition %u ", part);

    // already fully unpacked? if so, print it
    if (minAvailableFrame <= (*frameRange_it) && (*frameRange_it) <= maxAvailableFrame) {
      infoMsg(IM::Printing,IM::High,"is available to print:\n");
      if (part == 0) { // print P partition
	printOriginalSection(P_rvs, hidP_rvs, pVitTrigger != NULL,  pVitTriggerVec, 
			     pVitTriggerExpr, pTriggerEqn, printObserved, part, 'P', f, 
			     preg, pregex_mask, first_P, P_size, previous_P_values, vitRunLength, (*frameRange_it));
      } else if (part == totalOriginalPartitions-1) { // print E partition
	if ( (hidE_rvs.size() > 0)  || (printObserved && E_rvs.size() > 0) ) {
	  int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
	  shiftOriginalVarstoPosition(E_rvs, targetFrame, Epos);
	  printOriginalSection(E_rvs, hidE_rvs, eVitTrigger != NULL,  eVitTriggerVec, 
			       eVitTriggerExpr, eTriggerEqn, printObserved, part, 'E', f, 
			       ereg, eregex_mask, first_E, E_size, previous_E_values, vitRunLength, (*frameRange_it));
	}
      } else {      // print C partition
	int targetFrame = fp.numFramesInP() + (int)(part-1) * fp.numFramesInC();
	originalIndex = ((int)part - 1) % (int) C_rvs.size();
	shiftOriginalVarstoPosition(C_rvs[originalIndex], targetFrame, Cpos[originalIndex]);
	printOriginalSection(C_rvs[originalIndex], hidC_rvs[originalIndex], cVitTrigger != NULL, cVitTriggerVec,
			     cVitTriggerExpr, cTriggerEqn, printObserved, part, 'C', f, 
			     creg, cregex_mask, first_C, C_size, previous_C_values, vitRunLength, (*frameRange_it));
      }
      (*frameRange_it)++;  // move on to next frame
      continue;
    }

    /* Before we can print this original partition $C_j$ (j = part), we must 
       unpack the modified partitions 

       $$\left\{ C'_i \left| \, \max\left(-1,\left\lceil\frac{j-s-m+1}{s}\right\rceil\right) \leq i 
         \leq \left\lfloor \frac{j}{s} \right\rfloor \right. \right\}$$

       where $s$ and $m$ are the boundary algorithm parameters, $C'_{-1}=P'$, and $C'_{N_{C'}}=E'$. 
     */
    
    int numerator = ( (int)part - 1 - (int)S - (int)M + 1 );
    int firstPrimePart;
    if (numerator <= -(int)S) {
      firstPrimePart = -1; // max
    } else if (numerator <= 0){
      firstPrimePart = 0;  // ceil
    } else {
      firstPrimePart = numerator / (int)S;
      if (numerator % (int)S)
	firstPrimePart += 1;   // ceil
    }
    firstPrimePart += 1; // account for C_{-1} = P'

    unsigned lastPrimePart = (part > 0) ?  1 + ((int)part-1) / (int)S : 0;
    if (lastPrimePart >= inference_it.pt_len())  // E original partition # may > # of modified partitions
      lastPrimePart = inference_it.pt_len() - 1;

    infoMsg(IM::Printing,IM::High,"requires unpacking modified partitions %u to %u:\n  unpack:", 
	    firstPrimePart, lastPrimePart);

    for (unsigned i = (unsigned) firstPrimePart; i <= lastPrimePart; i += 1) { // unpack C'_{i} set
      infoMsg(IM::Printing,IM::High,"  %u'", i); 
      setCurrentInferenceShiftTo(i);
      PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

      if (binaryViterbiFile) { // load packed values from disk if not already in memory
	readBinaryVitPartition(ps, i);
      }
      
      // unpack
      if (inference_it.at_p()) { // P'
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(P_partition_values.ptr,PprimeValuePtrs.ptr);
      } else if (inference_it.at_e()) { // E'
	// -1 to get the preceding C', -1 to account for C'_{-1} = P'
	primeIndex = ((int)i - 2) % (int)Eprime_rvs.size(); 
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(E_partition_values.ptr,EprimeValuePtrs[primeIndex].ptr);
      } else { // C'
	assert ( inference_it.at_c() );
	primeIndex = ((int)i - 1) % (int)Cprime_rvs.size();
	if (ps.packer.packedLen() > 0) 
	  ps.packer.unpack(C_partition_values.ptr  + 
			   ( binaryViterbiFile ? 0 : (inference_it.pt_i()-1)*ps.packer.packedLen() ),
			   CprimeValuePtrs[primeIndex].ptr);
      }
    }

    // unpacking firstPrimePart ... lastPrimePart makes these frames available:
    if (firstPrimePart == 0) {
      minAvailableFrame = 0;
    } else {
      minAvailableFrame = NP + (  ( (firstPrimePart-1) * S + M ) * NC  );
    }

    if (lastPrimePart == inference_it.pt_len()-1) {
      maxAvailableFrame = numUsableFrames - 1;
    } else {
      maxAvailableFrame = NP + lastPrimePart * S * NC - 1;
    }
    infoMsg(IM::Printing,IM::High,"  available frames %u to %u\n", minAvailableFrame, maxAvailableFrame);
  }
  
  delete frameRange;
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setRootToMaxCliqueValue()
 *   
 *   This version of probEvidence() merely sets E's root clique to its
 *   current max value, and return the max value.
 *
 * See Also:
 *   probEvidence() above, the other (overloaded) version which is
 *   a const. mem version of collectEvidence().
 *
 * Preconditions:
 *   - collectEvidence() (or the other probEvidence()) must have been called
 *     for this to work.
 *   - The E partition must exist in the partitionStructureArray array
 *     and must have been called for this to work right.
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   The probability of the evidence
 *
 *-----------------------------------------------------------------------
 */
logpr
JunctionTree::setRootToMaxCliqueValue()
{
  ptps_iterator ptps_it(*this);
  ptps_it.set_to_last_entry();
  if (partitionStructureArray[ptps_it.ps_i()].maxCliquesSharedStructure.size() == 0)
    --ptps_it;

  return partitionTableArray[ptps_it.pt_i()].maxCliques[ptps_it.cur_ri()].
    maxProbability(partitionStructureArray[ptps_it.ps_i()].maxCliquesSharedStructure[ptps_it.cur_ri()],
		   true);

}



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//   Forward Messages
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::ceGatherIntoRoot
 *   
 *   Collect Evidence Gather Into Root: This routine does a collect
 *   evidence pass for this partition, and gathers all messages into
 *   the root within the current partition. It does so using the
 *   message order given in the argument 'message_order', and gathers
 *   into the provided root clique.  
 *
 * See Also:
 *   Dual routine: JunctionTree::deScatterOutofRoot()
 *
 *
 * Preconditions:
 *   It is assumed that either:
 *     1) this is the left-most partition
 *  or 2) that the left interface clique within this partition has had a message
 *        sent to it from the left neighbor partition.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::ceGatherIntoRoot(PartitionStructures& ps,
			       PartitionTables& pt,
			       // index of root clique in the partition
			       const unsigned root,
			       // message order of the JT in this partition
			       vector< pair<unsigned,unsigned> >& message_order,
			       // the name of the partition (for debugging/status msgs)
			       const char*const part_type_name,
			       // number of the partition in unrolled graph 
			       // (for printing/debugging/status msgs only)
			       const unsigned part_num,
			       const bool clearWhenDone,
			       const bool alsoClearOrigins)
{
  // first check that this is not an empty partition.
  if (ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! partitionDebugRange.contains((int)part_num)) {
#if 0
    printf("ceGather [part %u]: lowering inference level to %d\n", 
	   part_num, IM::glbMsgLevel(IM::DefaultModule));
#endif
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }

  bool zeroClique = false;
  try {
    // Now, do partition messages.
    for (unsigned msgNo=0;msgNo < message_order.size(); msgNo ++) {
      const unsigned from = message_order[msgNo].first;
      const unsigned to = message_order[msgNo].second;
      infoMsg(IM::Inference, IM::Med+5,
	      "CE: gathering into %s,part[%d]: clique %d\n",
	      part_type_name,part_num,from);

      // this may now throw an exception on zero clique errors - RR
      pt.maxCliques[from].
	ceGatherFromIncommingSeparators(ps.maxCliquesSharedStructure[from],
					pt.separatorCliques,
					ps.separatorCliquesSharedStructure.ptr);
  
      infoMsg(IM::Inference, IM::Mod,
	      "CE: message %s,part[%d]: clique %d --> clique %d\n",
	      part_type_name,part_num,from,to);
      pt.maxCliques[from].
	ceSendToOutgoingSeparator(ps.maxCliquesSharedStructure[from],
				  pt.separatorCliques,
				  ps.separatorCliquesSharedStructure.ptr);

      // TODO: if we are just computing probE here, we should delete
      // memory in pt.maxCliques[from]. Also, if we're only doing probE,
      // we should not keep the cliques around at all, only the outgoing
      // separator.
      if (clearWhenDone) {
	pt.maxCliques[from].
	  clearCliqueAndIncommingSeparatorMemory(ps.maxCliquesSharedStructure[from],
						 pt.separatorCliques,
						 ps.separatorCliquesSharedStructure.ptr);

	if (alsoClearOrigins) {
	  // then clear out the origin memory used for inference.
	  ps.origin.clearCliqueAndIncommingSeparatorMemoryForClique(from); 
	}
      }
    }
  } catch (ZeroCliqueException const &e) {
    zeroClique = true; // abort this partition & segment
  }
  if (!zeroClique) {
    // collect to partition's root clique
    infoMsg(IM::Inference, IM::Med+5,
	    "CE: gathering into partition root %s,part[%d]: clique %d\n",
	    part_type_name,part_num,root);
    try {
      pt.maxCliques[root].
	ceGatherFromIncommingSeparators(ps.maxCliquesSharedStructure[root],
					pt.separatorCliques,
					ps.separatorCliquesSharedStructure.ptr);
    } catch (ZeroCliqueException const &e) {
      zeroClique = true; // abort this partition & segment
    }
    if (!zeroClique) {
      if (IM::messageGlb(IM::InferenceMemory, IM::Med+9)) {
	pt.reportMemoryUsageTo(ps,stdout);
      }
    }
  }
  if (! partitionDebugRange.contains((int)part_num)) {
#if 0
    printf("ceGather [part %u]: raising inference level to %d\n", 
	   part_num, inferenceDebugLevel);
#endif
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }
  if (zeroClique) {
    throw ZeroCliqueException(); // continue to abort segment
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::ceSendForwardsCrossPartitions
 *   
 *   Collect Evidence Send To Next Partition: This routine sends a
 *   message from the right interface clique of a left (or previous)
 *   partition to the left interface clique of a right (or next)
 *   partition in the partition series. It is assumed that the right
 *   interface clique has had all its incomming messages sent to it.
 *
 * See Also:
 *   Dual routine: JunctionTree::deSendBackwardsCrossPartitions()
 *
 *
 * Preconditions:
 *   It is assumed that:
 *     1) the right interface of the previous partition must have had
 *        all messages sent to it.
 * 
 * Postconditions:
 *     the left interface of the next partition is now set up.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::ceSendForwardsCrossPartitions(// previous partition
				    PartitionStructures& previous_ps,
				    PartitionTables& previous_pt,
				    // root clique of the previous partition (i.e., the
				    // right interface clique) 
				    const unsigned previous_part_root,
				    // name of previous partition (for debugging/status msgs)
				    const char*const previous_part_type_name,
				    // sequence number (in unrolling) of previous partition
				    // (for debugging/status msgs)
				    const unsigned previous_part_num,
				    // next partition
				    PartitionStructures& next_ps,
				    PartitionTables& next_pt,
				    // leaf clique of next partition (i.e., index number
				    // of the left interface clique of next partition)
				    const unsigned next_part_leaf,
				    // name (debugging/status msgs)
				    const char*const next_part_type_name,
				    // partitiiton number (debugging/status msgs)
				    const unsigned next_part_num)
{
  // check for empty partitions.
  if (previous_ps.maxCliquesSharedStructure.size() == 0 || next_ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! partitionDebugRange.contains((int)next_part_num)) {
#if 0
    printf("ceGather [part %u]: lowering inference level to %d\n", 
	   next_part_num, IM::glbMsgLevel(IM::DefaultModule));
#endif
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }


  infoMsg(IM::Inference, IM::Mod,"CE: message %s,part[%d],clique(%d) --> %s,part[%d],clique(%d)\n",
	  previous_part_type_name,
	  previous_part_num,
	  previous_part_root,
	  next_part_type_name,
	  next_part_num,
	  next_part_leaf);
  previous_pt.maxCliques[previous_part_root].
    ceSendToOutgoingSeparator(previous_ps.maxCliquesSharedStructure[previous_part_root],
			      next_pt.separatorCliques[next_ps.separatorCliquesSharedStructure.size()-1],
			      next_ps.separatorCliquesSharedStructure[next_ps.separatorCliquesSharedStructure.size()-1]);

  if (IM::messageGlb(IM::InferenceMemory, IM::Med+9)) {
    previous_pt.reportMemoryUsageTo(previous_ps,stdout);
  }

  if (! partitionDebugRange.contains((int)next_part_num)) {
#if 0
    printf("ceGather [part %u]: raising inference level to %d\n", 
	   next_part_num, inferenceDebugLevel);
#endif
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectEvidence()
 *   
 *   Collect Evidence: This routine performes a complete collect
 *   evidence pass for this series of partitions that have been
 *   unrolled a certain amount, given by the partitionStructureArray array. It
 *   sets the appropriate names, etc. of the partitions depending on
 *   if this is a left-interface or right-interface form of inference.
 *   the root within the current partition. 
 *  
 *   This routine demonstrates a simple version of linear space
 *   collect evidence inference. It keeps everything in memory at the
 *   same time when doing going forward, so it is suitable for use in
 *   a collect evidence/distribute evidence (forward/backward)
 *   framework.  A companion routine (aptly named
 *   'distributeEvidence') will, assuming collect evidence has been
 *   called, do the distribute evidence stage, thereby leaving all
 *   cliques in all partitions locally and therefore globally
 *   consistent, ready to be used say in EM training.
 *
 * See Also:
 *    1) contant memory (not dept. on time T) version of collect evidence
 *       JunctionTree::probEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *    3) O(Tn) memory version of collect/distribute evidence
 *       JunctionTree::collectEvidenceOnlyKeepSeps / JunctionTree::distributeEvidenceOnlyKeepSeps
 *
 * Preconditions:
 *   The parititons must have been created and placed in the array partitionStructureArray.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::collectEvidence()
{
  // This routine handles all of:
  // 
  //    unrolled 0 times: (so there is a single P1, and E1)  
  //    unrolled 1 time: so there is a P1, C1, C3, E1
  //    unrolled 2 or more times: so there is a P1 C1 [C2 ...] C3, E1
  //    etc.

  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this);

  init_CC_CE_rvs(inference_it);

  // we skip the first Co's LI separator if there is no P1
  // partition, since otherwise we'll get zero probability.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.skipLISeparator();
  // gather into the root of the current  partition
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   partitionTableArray[inference_it.pt_i()],
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i());

  if (sectionDoDist) {
  // Send messages from the root clique to the rest of the cliques
  // in this partition so that they are consistant with the observations
  // in this partition. We originally wanted to send messages only to
  // the cliques actually being printed, but deScatterToOutgoingSeparators()
  // sends messages to all of a clique's outgoing separators (rather
  // than just those on the path to a printing clique) and we decided
  // not to implement a "subset" scatter. We think that in the common
  // cases there won't be much extra work from the full scatter.
    deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		       partitionTableArray[inference_it.pt_i()],
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());
  }

  // if the LI separator was turned off, we need to turn it back on.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.useLISeparator();

  for (unsigned part=1;part < inference_it.pt_len(); part ++ ) {

    setCurrentInferenceShiftTo(part);

    // inference_it.printState(stdout);

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[inference_it.ps_prev_i()],
			  partitionTableArray[inference_it.pt_prev_i()],
			  inference_it.prev_ri(),
			  inference_it.prev_nm(),
			  inference_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[inference_it.ps_i()],
			  partitionTableArray[inference_it.pt_i()],
			  inference_it.cur_li(),
			  inference_it.cur_nm(),
			  inference_it.pt_i());


    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // it might be that E is the first partition as well, say if this is
    // a static graph, and in this case we need in this case to skip the
    // incomming separator, which doesn't exist.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();
    // next, gather into the root of the final E partition
    ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		     partitionTableArray[inference_it.pt_i()],
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());

    if (sectionDoDist) {
      // Send messages from the root clique to the rest of the cliques
      // in this partition so that they are consistant with the observations
      // in this partition. We originally wanted to send messages only to
      // the cliques actually being printed, but deScatterToOutgoingSeparators()
      // sends messages to all of a clique's outgoing separators (rather
      // than just those on the path to a printing clique) and we decided
      // not to implement a "subset" scatter. We think that in the common
      // cases there won't be much extra work from the full scatter.
      deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
			 partitionTableArray[inference_it.pt_i()],
			 inference_it.cur_ri(),
			 inference_it.cur_message_order(),
			 inference_it.cur_nm(),
			 inference_it.pt_i());
    }
    
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();

    // if the LI separator was turned off, we need to turn it back on.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();

  }
  assert ( inference_it.at_e() );


#if 0
  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this);

  init_CC_CE_rvs(inference_it);
  setCurrentInferenceShiftTo(0);

  // we skip the first Co's LI separator if there is no P1
  // partition, since otherwise we'll get zero probability.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.skipLISeparator();
  // gather into the root of the current  partition
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   partitionTableArray[inference_it.pt_i()],
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i());
  // if the LI separator was turned off, we need to turn it back on.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.useLISeparator();

  for (unsigned part=1;part < inference_it.pt_len(); part ++ ) {
    setCurrentInferenceShiftTo(part);

    // inference_it.printState(stdout);

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[inference_it.ps_prev_i()],
			  partitionTableArray[inference_it.pt_prev_i()],
			  inference_it.prev_ri(),
			  inference_it.prev_nm(),
			  inference_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[inference_it.ps_i()],
			  partitionTableArray[inference_it.pt_i()],
			  inference_it.cur_li(),
			  inference_it.cur_nm(),
			  inference_it.pt_i());


    // it might be that E is the first partition as well, say if this is
    // a static graph, and in this case we need in this case to skip the
    // incomming separator, which doesn't exist.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();
    // next, gather into the root of the final E partition
    ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		     partitionTableArray[inference_it.pt_i()],
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();
  }
  assert ( inference_it.at_e() );

#endif

}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectEvidenceOnlyKeepSeps()
 *   
 *   Collect Evidence: This routine performes a collect evidence pass 
 *   similar to the constant memory JunctionTree::probEvidenceFixedUnroll(),
 *   but it keeps the separator cliques around for use in the corresponding
 *   JunctionTree::distributeEvidenceOnlyKeepSeps() backwards pass.
 *
 * See Also:
 *    1) contant memory (not dept. on time T) version of collect evidence
 *       JunctionTree::probEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *   The parititons must have been created and placed in the array partitionStructureArray.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *     The separator cliques are remembered in interfaceTemp for later use
 *     in JunctionTree::distributeEvidenceOnlyKeepSeps()
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */

logpr
JunctionTree::collectEvidenceOnlyKeepSeps(const unsigned int numFrames,
				    unsigned* numUsableFrames)
{
  FileSource *gomFS;
  // This should be safe since gmtkOnline is the only program
  // that does inference and doesn't use FileSource and gmtkOnline
  // only uses onlineFixedUnroll
  gomFS= static_cast<FileSource *>(globalObservationMatrix);
  assert(typeid(*globalObservationMatrix) == typeid(*gomFS));

  // Unroll, but do not use the long table array (2nd parameter is
  // false) for allocation, but get back the long table length
  // in a local variable for use in our iterator.
  unsigned totalNumberPartitions;
  {
    unsigned tmp = unroll(numFrames,ZeroTable,&totalNumberPartitions);
    gomFS->justifySegment(tmp);
    if (numUsableFrames) 
      *numUsableFrames = tmp;
    // limit scope of tmp.
  }

  // This routine handles all of:
  // 
  //    unrolled 0 times: (so there is a single P1, and E1)  
  //    unrolled 1 time: so there is a P1, C1, C3, E1
  //    unrolled 2 or more times: so there is a P1 C1 [C2 ...] C3, E1
  //    etc.

  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this, totalNumberPartitions);

  init_CC_CE_rvs(inference_it);

  PartitionTables* prev_part_tab = NULL;
  PartitionTables* cur_part_tab = new PartitionTables(inference_it.cur_jt_partition());

  interfaceTemp.resize(inference_it.pt_len());

  // we skip the first Co's LI separator if there is no P1
  // partition, since otherwise we'll get zero probability.
  if (inference_it.at_first_c() && P1.cliques.size() == 0) {
    Co.skipLISeparator();
  }

  // gather into the root of the current  partition
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   *cur_part_tab,
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i());

  // if the LI separator was turned off, we need to turn it back on.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.useLISeparator();

  unsigned part;
  for (part=1; part < inference_it.pt_len(); part ++ ) {

    delete prev_part_tab;
    prev_part_tab = cur_part_tab;

    setCurrentInferenceShiftTo(part);
    cur_part_tab = new PartitionTables(inference_it.cur_jt_partition());
      
    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[inference_it.ps_prev_i()],
			  *prev_part_tab,
			  inference_it.prev_ri(),
			  inference_it.prev_nm(),
			  inference_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[inference_it.ps_i()],
			  *cur_part_tab,
			  inference_it.cur_li(),
			  inference_it.cur_nm(),
			  inference_it.pt_i());

    // remember the incoming separator so the forward pass can be re-done in the backward pass
    unsigned interfaceCliqueNum = 
      partitionStructureArray[inference_it.ps_i()].separatorCliquesSharedStructure.size()-1;
    cur_part_tab->separatorCliques[ interfaceCliqueNum ].preserve = true;
    interfaceTemp[ inference_it.pt_i() ] = cur_part_tab->separatorCliques[ interfaceCliqueNum ];

    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // it might be that E is the first partition as well, say if this is
    // a static graph, and in this case we need in this case to skip the
    // incomming separator, which doesn't exist.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();
    // next, gather into the root of the partition
    ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		     *cur_part_tab,
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());

    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();
    // if the LI separator was turned off, we need to turn it back on.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();
  }
  logpr probE;
  if (viterbiScore) {
    probE = cur_part_tab->maxCliques[E_root_clique].maxProb();
  } else {
    probE = cur_part_tab->maxCliques[E_root_clique].sumProbabilities();
  }
  delete cur_part_tab;
  delete prev_part_tab;
  return probE;
}




///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//   Backwards Messages
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::deScatterOutofRoot
 *   
 *   Distribute Evidence Scatter Outof Root: This routine does a
 *   distribute evidence pass for this partition, and scatters all
 *   messages outof the root within the current partition. It does so
 *   using the message order given in the argument 'message_order',
 *   and scatters out of the provided root clique (which is the right
 *   interface clique of this partition). By "Scatter", I mean it
 *   sends messages from the root clique distributing everything
 *   ultimately to all leaf cliques in this partition.
 *
 * See Also:
 *   Dual routine: JunctionTree::ceGatherIntoRoot()
 *
 *
 * Preconditions:
 *   It is assumed that either:
 *     1) this is the ritht-most partition
 *  or 2) that the right interface clique within this partition has had a message
 *        sent to it from the right neighbor partition.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::deScatterOutofRoot(// the partition
				 PartitionStructures& ps,
				 PartitionTables& pt,
				 // root (right interface clique) of this partition
				 const unsigned root,
				 // message order
				 vector< pair<unsigned,unsigned> >& message_order,
				 // name (debugging/status msgs)
				 const char*const part_type_name,
				 // partition number (debugging/status msgs)
				 const unsigned part_num)
{
  if (ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! partitionDebugRange.contains((int)part_num)) {
#if 0
    printf("deScatter [part %u]: lowering inference level to %d\n", part, IM::glbMsgLevel(IM::DefaultModule));
#endif
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }

  infoMsg(IM::Inference, IM::Med+5,"DE: distributing out of partition root %s,part[%d]: clique %d\n",
	  part_type_name,part_num,root);
  pt.maxCliques[root].
    deScatterToOutgoingSeparators(ps.maxCliquesSharedStructure[root],
				  pt.separatorCliques,
				  ps.separatorCliquesSharedStructure.ptr);
  for (unsigned msgNoP1=message_order.size();msgNoP1 > 0; msgNoP1 --) {
    const unsigned to = message_order[msgNoP1-1].first;
    const unsigned from = message_order[msgNoP1-1].second;
    infoMsg(IM::Inference, IM::Mod,"DE: message %s,part[%d]: clique %d <-- clique %d\n",
	    part_type_name,part_num,to,from);
    pt.maxCliques[to].
      deReceiveFromIncommingSeparator(ps.maxCliquesSharedStructure[to],
				      pt.separatorCliques,
				      ps.separatorCliquesSharedStructure.ptr);

    infoMsg(IM::Inference, IM::Med+5,"DE: distributing out of %s,part[%d]: clique %d\n",
	    part_type_name,part_num,to);
    pt.maxCliques[to].
      deScatterToOutgoingSeparators(ps.maxCliquesSharedStructure[to],
				    pt.separatorCliques,
				    ps.separatorCliquesSharedStructure.ptr);
  }

  if (! partitionDebugRange.contains((int)part_num)) {
#if 0
    printf("deScatter [part %u]: raising inference level to %d\n", part, inferenceDebugLevel);
#endif
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }


}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::deSendBackwardsCrossPartitions()
 *   
 *   Distribute Evidence Receive To Previous Partition: This routine
 *   sends a message from the left interface clique of a right (or
 *   next) partition to the right interface clique of a left (or
 *   previous) partition in the partition series. It is assumed that
 *   the left interface clique has had all its messages sent to it.
 *
 * See Also:
 *   Dual routine: JunctionTree::ceSendForwardsCrossPartitions()
 *
 * Preconditions:
 *   It is assumed that:
 *     1) the left interface of the next partition must have had
 *        all messages sent to it.
 * 
 * Postconditions:
 *     the right interface of the previous partition is now set up.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::deSendBackwardsCrossPartitions(// previous partition
					     PartitionStructures& previous_ps,
					     PartitionTables& previous_pt,
					     // root (right interface clique) of previous partition 
					     const unsigned previous_part_root,
					     // name of prev part (debugging/status msgs)
					     const char*const previous_part_type_name,
					     // number of prev part (debugging/status msgs)
					     const unsigned previous_part_num,

					     // the next partition
					     PartitionStructures& next_ps,
					     PartitionTables& next_pt,
					     // leaf (left interface cliuqe) of next partition
					     const unsigned next_part_leaf,
					     // name of next part (debugging/status msgs)
					     const char*const next_part_type_name,
					     // number of next part (debugging/status msgs)
					     const unsigned next_part_num
					     )
{
  // check for empty partitions.
  if (previous_ps.maxCliquesSharedStructure.size() == 0 || next_ps.maxCliquesSharedStructure.size() == 0)
    return;

  unsigned inferenceDebugLevel = IM::glbMsgLevel(IM::Inference);
  unsigned inferenceMemoryDebugLevel = IM::glbMsgLevel(IM::InferenceMemory);

  if (! partitionDebugRange.contains((int)previous_part_num)) {
#if 0
    printf("deScatter [part %u]: lowering inference level to %d\n", 
	   previous_part_num, IM::glbMsgLevel(IM::DefaultModule));
#endif
    IM::setGlbMsgLevel(IM::Inference, IM::glbMsgLevel(IM::DefaultModule));
    IM::setGlbMsgLevel(IM::InferenceMemory, IM::glbMsgLevel(IM::DefaultModule));
  }

  infoMsg(IM::Inference, IM::Mod,"DE: message %s,part[%d],clique(%d) <-- %s,part[%d],clique(%d)\n",
	  previous_part_type_name,previous_part_num,previous_part_root,
	  next_part_type_name,next_part_num,next_part_leaf);
  previous_pt.maxCliques[previous_part_root].
    deReceiveFromIncommingSeparator(previous_ps.maxCliquesSharedStructure[previous_part_root],
				    next_pt.separatorCliques[next_ps.separatorCliquesSharedStructure.size()-1],
				    next_ps.separatorCliquesSharedStructure[next_ps.separatorCliquesSharedStructure.size()-1]);

  if (! partitionDebugRange.contains((int)previous_part_num)) {
#if 0
    printf("deScatter [part %u]: raising inference level to %d\n", 
	   previous_part_num, inferenceDebugLevel);
#endif
    IM::setGlbMsgLevel(IM::InferenceMemory, inferenceMemoryDebugLevel);
    IM::setGlbMsgLevel(IM::Inference, inferenceDebugLevel);
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::distributeEvidence()
 *   
 *   Distribute Evidence: This routine performs a complete distribute
 *   evidence pass for this series of partitions that have been
 *   unrolled a certain amount, given by the partitionStructureArray array. It
 *   sets the appropriate names, etc. of the partitions depending on
 *   if this is a left-interface or right-interface form of inference.
 *   the root within the current partition.
 *  
 *   This routine demonstrates a simple version of linear space
 *   distribute evidence inference. It uses data kept in memory when
 *   doing going backwards, and it is assumed that it will be used in
 *   tandem with a previous collect evidence call.  A companion
 *   routine ('collectEvidence') will, do collect evidence appropriately
 *   leaving all data structures set up for distributeEvidence() to be called.
 *
 *   After distributeEvidence() is called, all cliques in all
 *   partitions will be locally (& globally if it is a JT) consistant, and so will
 *   be ready for EM training, etc.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) contant memory (not dept. on time T) version of collect evidence
 *       JunctionTree::probEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *   - collectEvidence() must have been called right before this.
 *   - The parititons must have been created and placed in the array partitionStructureArray.
 *   - inference_it must be initialized to the current segment
 *   - the pair of partition RVS set up correct.
 * 
 *
 * Postconditions:
 *     All cliques are now locally consistant, and will be globally consistant
 *     if we have a junction tree. 
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::distributeEvidence()
{
  for (unsigned part= (inference_it.pt_len()-1) ; part > 0 ; part -- ) {

    setCurrentInferenceShiftTo(part);
    
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();    
    else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();

    deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		       partitionTableArray[inference_it.pt_i()],
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();
    else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();

    if (viterbiScore)
      recordPartitionViterbiValue(inference_it);

    // send backwads message to previous partition
    deSendBackwardsCrossPartitions(partitionStructureArray[inference_it.ps_prev_i()],
				   partitionTableArray[inference_it.pt_prev_i()],
				   inference_it.prev_ri(),
				   inference_it.prev_nm(),
				   inference_it.pt_prev_i(),
				   //
				   partitionStructureArray[inference_it.ps_i()],
				   partitionTableArray[inference_it.pt_i()],
				   inference_it.cur_li(),
				   inference_it.cur_nm(),
				   inference_it.pt_i());
  }

  setCurrentInferenceShiftTo(0);
  // do the final scatter out of root of initial P partition.

  deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		     partitionTableArray[inference_it.pt_i()],
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());

  if (viterbiScore)
    recordPartitionViterbiValue(inference_it);



#if 0

  ptps_iterator ptps_it(*this);
  ptps_it.set_to_last_entry();

  // the set of rvs that need to be adjusted, if any.
  set <RV*> two_part_rvs;
  if (ptps_it.num_c_partitions() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // partitions (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = partitionStructureArray[3].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,two_part_rvs);

    // we start by adjusting C E up to appropriate values.
    adjustFramesBy (two_part_rvs, (ptps_it.num_c_partitions()-2)
		    *gm_template.S*fp.numFramesInC());    

  } 


  while (!ptps_it.at_first_entry()) {

    // ptps_it.printState(stdout);

    if (ptps_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();    
    else if (!ptps_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();

    deScatterOutofRoot(partitionStructureArray[ptps_it.ps_i()],
		       partitionTableArray[ptps_it.pt_i()],
		       ptps_it.cur_ri(),
		       ptps_it.cur_message_order(),
		       ptps_it.cur_nm(),
		       ptps_it.pt_i());
    if (ptps_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();
    else if (!ptps_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();
    if (viterbiScore)
      recordPartitionViterbiValue(ptps_it);

    // send backwads message to previous partition
    deSendBackwardsCrossPartitions(partitionStructureArray[ptps_it.ps_prev_i()],
				   partitionTableArray[ptps_it.pt_prev_i()],
				   ptps_it.prev_ri(),
				   ptps_it.prev_nm(),
				   ptps_it.pt_prev_i(),
				   //
				   partitionStructureArray[ptps_it.ps_i()],
				   partitionTableArray[ptps_it.pt_i()],
				   ptps_it.cur_li(),
				   ptps_it.cur_nm(),
				   ptps_it.pt_i());


    if (ptps_it.num_c_partitions() > 2) {
      // then we need to do some RV shifting.
      if (ptps_it.at_e()) {
	// we need to shift the two C',E' parts back to normal, and then shift
	// the (C,C) part forward appropriately

	// First, shift C',E' rvs back to normal frame
	adjustFramesBy (two_part_rvs, -(ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

	set <RV*> tmp1,tmp2;
	tmp1 = partitionStructureArray[1].returnRVsAndTheirObservedParentsAsSet();
	tmp2 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
	// We're done now with two_part_rvs, so we make a new one. The
	// next routine clears out the original set values.
	unionRVs(tmp1,tmp2,two_part_rvs);
	// now adjust C',C' variables to match.
	adjustFramesBy (two_part_rvs, (ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

      } else if (ptps_it.pt_i() > 2) {
	// shift the two C partitions over by one more
	adjustFramesBy (two_part_rvs, -1*gm_template.S*fp.numFramesInC());
      }
    }

    --ptps_it;
  }

  // do the final scatter out of root of initial P partition.
  deScatterOutofRoot(partitionStructureArray[ptps_it.ps_i()],
		     partitionTableArray[ptps_it.pt_i()],
		     ptps_it.cur_ri(),
		     ptps_it.cur_message_order(),
		     ptps_it.cur_nm(),
		     ptps_it.pt_i());



  if (viterbiScore)
    recordPartitionViterbiValue(ptps_it);

#endif
}

/*-
 *-----------------------------------------------------------------------
 * JunctionTree::distributeEvidence()
 *   
 *   Distribute Evidence: This routine performs a complete distribute
 *   evidence pass for this series of partitions that have run
 *   JunctionTree::collectEvidenceOnlyKeepSeps()
 *
 *   After distributeEvidence() is called, all cliques in all
 *   partitions will be locally (& globally if it is a JT) consistant, and so will
 *   be ready for EM training, etc.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) contant memory (not dept. on time T) version of collect evidence
 *       JunctionTree::probEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *   - collectEvidenceOnlyKeepSeps() must have been called right before this.
 *   - The parititons must have been created and placed in the array partitionStructureArray.
 *   - inference_it must be initialized to the current segment
 *   - the pair of partition RVS set up correct.
 * 
 *
 * Postconditions:
 *     All cliques are now locally consistant, and will be globally consistant
 *     if we have a junction tree. 
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */

void
JunctionTree::distributeEvidenceOnlyKeepSeps()
{
  setCurrentInferenceShiftTo(inference_it.pt_len()-1);

  PartitionTables* prev_part_tab = NULL;
  PartitionTables* cur_part_tab  = new PartitionTables(inference_it.cur_jt_partition());

  unsigned interfaceCliqueNum;

  // re-do E' gather into root using remembered separator clique
  if (!inference_it.at_first_entry()) {
    interfaceCliqueNum = partitionStructureArray[inference_it.ps_i()].separatorCliquesSharedStructure.size()-1;
    cur_part_tab->separatorCliques[ interfaceCliqueNum ] = interfaceTemp[ inference_it.pt_i() ];
    interfaceTemp[ inference_it.pt_i() ].preserve = false;
  }

  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.skipLISeparator();
  else  if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
    E1.skipLISeparator();
  
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   *cur_part_tab,
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i(), false, false);

  cur_part_tab->maxCliques[inference_it.cur_ri()].
    maxProbability(partitionStructureArray[inference_it.ps_i()].maxCliquesSharedStructure[inference_it.cur_ri()], true);

  for (unsigned part= (inference_it.pt_len()-1) ; part > 0 ; part -- ) {

    setCurrentInferenceShiftTo(part);

    // distribute evidence to the current partition
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();
    else  if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();

    deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		       *cur_part_tab,
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());

    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();
    else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();

    if (viterbiScore)
      recordPartitionViterbiValue(inference_it);

    // reconstruct previous partition

    setCurrentInferenceShiftTo(part-1);
    prev_part_tab = new PartitionTables(inference_it.cur_jt_partition());

    if (!inference_it.at_first_entry()) {
      interfaceCliqueNum = partitionStructureArray[inference_it.ps_i()].separatorCliquesSharedStructure.size()-1;
      prev_part_tab->separatorCliques[ interfaceCliqueNum ] = interfaceTemp[ inference_it.pt_i() ];
      interfaceTemp[ inference_it.pt_i() ].preserve = false;
    }
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();
    else  if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();

    ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		     *prev_part_tab,
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());

    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();
    else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();
    setCurrentInferenceShiftTo(part);

    // send backwads message to previous partition
    deSendBackwardsCrossPartitions(partitionStructureArray[inference_it.ps_prev_i()], // previous
                                   *prev_part_tab,
				   inference_it.prev_ri(),
				   inference_it.prev_nm(),
				   inference_it.pt_prev_i(),
				   // current partition
				   partitionStructureArray[inference_it.ps_i()],
				   *cur_part_tab,
				   inference_it.cur_li(),
				   inference_it.cur_nm(),
				   inference_it.pt_i());
    //    partitionTableArray[inference_it.pt_i()] = *cur_part_tab;

    delete cur_part_tab;
    cur_part_tab = prev_part_tab;
  }

  setCurrentInferenceShiftTo(0);

  // do the final scatter out of root of initial P partition.
  deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		     *cur_part_tab,
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());

  if (viterbiScore)
    recordPartitionViterbiValue(inference_it);
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::emIncrement()
 *
 *    A version of emIncrement that works with the data structures set up by 
 *    collectEvidence() and distributeEvidence() above. Note, EM training
 *    can also be performed with the island algorithm.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) distributeEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *    collectEvidence() AND distributeEvidence() must have just been
 *    called setting up the data structures. All cliques must exist and
 *    must have their clique tables filled out and ready.
 *    Also, all parametr accumulators must be set up and ready to go.
 *
 * Postconditions:
 *   The accumulators are increment accordingly.
 *
 * Side Effects:
 *   This will update the accumulators of all trainable parameter objects.
 *
 * Results:
 *   None
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::emIncrement(const logpr probE,
			  const bool localCliqueNormalization,
			  const double emTrainingBeam)
{

  // Quite simply, just iterate through all partitions and call emIncrement
  // therein.

  // forward order
  // 
  // for (unsigned part=0;part < inference_it.pt_len(); part ++ ) {
  // 
  // reverse order increment (which is the same as island algorithm order,
  // and this can make slight numeric differences). To avoid
  // questions about this ("why does island EM produce small differences
  // than normal EM???") we use the reverse order here.
  for (unsigned part=inference_it.pt_len();part > 0 ; part --) {
    setCurrentInferenceShiftTo(part-1);
    infoMsg(IM::Training, IM::High-1,
	    "EM: accumulating stats for %s,part[%d]\n",
	    inference_it.cur_nm(),inference_it.pt_i());
    partitionTableArray[inference_it.pt_i()].
      emIncrement(partitionStructureArray[inference_it.ps_i()],
		  probE,
		  localCliqueNormalization,
		  emTrainingBeam);
  }



#if 0
  // Quite simply, just iterate through all partitions and call emIncrement
  // therein.
  ptps_iterator ptps_it(*this);
  ptps_it.set_to_first_entry();

  // the set of rvs that need to be adjusted, if any.
  set <RV*> two_part_rvs;
  if (ptps_it.num_c_partitions() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // partitions (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = partitionStructureArray[1].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,two_part_rvs);
  }

  while (1) {
    infoMsg(IM::Training, IM::High-1,
	    "EM: accumulating stats for %s,part[%d]\n",
	    ptps_it.cur_nm(),ptps_it.pt_i());

    if (ptps_it.num_c_partitions() > 2) {
      // then we need to do some RV shifting.
      if (ptps_it.at_e()) {
	// we need to shift the two C1',C2' parts back to normal, and then shift
	// the (C,E) part forward.

	// First, shift C1',C2' rvs back to normal frame
	adjustFramesBy (two_part_rvs, -(ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

	set <RV*> tmp1,tmp2;
	tmp1 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
	tmp2 = partitionStructureArray[3].returnRVsAndTheirObservedParentsAsSet();
	// We're done now with two_part_rvs, so we make a new one. The
	// next routine clears out the original set values.
	unionRVs(tmp1,tmp2,two_part_rvs);
	// now adjust C2',E' variables to match.
	adjustFramesBy (two_part_rvs, (ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

      } else if (ptps_it.pt_i() > 2) {
	// shift the two C partitions over by one more
	adjustFramesBy (two_part_rvs, 1*gm_template.S*fp.numFramesInC());
      }
    }

    partitionTableArray[ptps_it.pt_i()].
      emIncrement(partitionStructureArray[ptps_it.ps_i()],
		  probE,
		  localCliqueNormalization,
		  emTrainingBeam);
    if (ptps_it.at_last_entry())
      break;
    else
      ++ptps_it;
  }


#endif
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::probEvidence()
 *   
 *   This version of probEvidence() merely sums up the entries in E's
 *   root clique and return the result. This, if collect evidence has
 *   been called all the way to E's root clique, this function will
 *   return the probabilty of the evidence prob(E).
 *
 *   Note that if we are in viterbi mode, then rather than summing
 *   we return the score of the maximum value clique entry.
 *
 * See Also:
 *   probEvidence() above, the other (overloaded) version which is
 *   a const. mem version of collectEvidence().
 *
 * Preconditions:
 *   - collectEvidence() (or the other probEvidence()) must have been called
 *   - The E partition must exist in the partitionStructureArray array
 *     and must have been called for this to work right.
 *   - both the partitionStructureArray and the partitionTableArray must
 *     be in a consistent state (i.e., unroll() must have been called
 *     and one of the above mentioned routines should have been successfully called).
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   The probability of the evidence
 *
 *-----------------------------------------------------------------------
 */
logpr
JunctionTree::probEvidence()
{

  inference_it.set_to_last_entry();
  // first check to see if there are any cliques in the final entry,
  // and if not we use the last C rather than the last E.
  if (partitionStructureArray[inference_it.ps_i()].maxCliquesSharedStructure.size() == 0)
    --inference_it;
  
  // this next routine is not necessary since the max and sum routines
  // do not require access to random variables, all the info for the
  // scores lives in the tables.
  // setCurrentInferenceShiftTo(inference_it.pt_i());

  if (viterbiScore) {
    return partitionTableArray[inference_it.pt_i()].maxCliques[inference_it.cur_ri()].maxProb();
  } else {
    return partitionTableArray[inference_it.pt_i()].maxCliques[inference_it.cur_ri()].sumProbabilities();
  }


#if 0
  ptps_iterator ptps_it(*this);
  ptps_it.set_to_last_entry();

  if (partitionStructureArray[ptps_it.ps_i()].maxCliquesSharedStructure.size() == 0)
    --ptps_it;

  if (viterbiScore) {
    return partitionTableArray[ptps_it.pt_i()].maxCliques[ptps_it.cur_ri()].maxProb();
  } else {
    return partitionTableArray[ptps_it.pt_i()].maxCliques[ptps_it.cur_ri()].sumProbabilities();
  }
#endif

}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::printProbEvidenceAccordingToAllCliques()
 *   
 *   This routine will cycle through all cliques of all partitions and will
 *   print out the sums of the probs. of each cliuqe. Therefore, if 
 *   collect/distribute evidence has been called (and it is working)
 *   this routine should print out exactly the same value for all
 *   cliques (to within numerical precision and logp.h table error roundoff).
 *
 *   This routine really is only used for debugging/status messages.
 *
 * Preconditions:
 *   - collectEvidence() and distributeEvidence() should have been called
 *     and the partitionStructureArray array left set up.
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::printProbEvidenceAccordingToAllCliques()
{
  for (unsigned part=0;part < inference_it.pt_len(); part ++ ) {
    // next call is not needed here.
    // setCurrentInferenceShiftTo(part);
    for (unsigned cliqueNo=0;cliqueNo<partitionStructureArray[inference_it.ps_i()].maxCliquesSharedStructure.size();cliqueNo++) {
      printf("Part no %d: clique no %d: log probE = %f\n",
	     inference_it.ps_i(),cliqueNo,
	     partitionTableArray[inference_it.pt_i()].maxCliques[cliqueNo].sumProbabilities().valref());
    }
  }

#if 0
  ptps_iterator ptps_it(*this);
  ptps_it.set_to_first_entry();

  while (1) {
    for (unsigned cliqueNo=0;cliqueNo<partitionStructureArray[ptps_it.ps_i()].maxCliquesSharedStructure.size();cliqueNo++) {
      printf("Part no %d: clique no %d: log probE = %f\n",
	     ptps_it.ps_i(),cliqueNo,
	     partitionTableArray[ptps_it.pt_i()].maxCliques[cliqueNo].sumProbabilities().valref());
    }
    if (ptps_it.at_last_entry())
      break;
    else
      ++ptps_it;    
  }
#endif

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fixed memory versions of computing the probability of evidence
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::probEvidenceFixedUnroll()
 *
 *    A constant memory (i.e., indep. of T), combination of unroll and
 *    collectEvidence() above. It is constnat memory in that
 *    it keeps no more than two partitions in memory simultaneously.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) distributeEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *    same as unroll() above.
 *
 * Postconditions:
 *   The probability of the evidence is returned to the caller.
 *
 * Side Effects:
 *   This will update the space managers averages and statistics.
 *
 *   TODO: hash table usage, sharingetc. so that hash tables are not re-used 
 *         for the entire length.
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
/*
 * a version of the above routine but that only unrolls by
 * a fixed amount, there is only a constant (rather than a linear cost)
 * for the graph data structures.
 */

logpr 
JunctionTree::probEvidenceFixedUnroll(const unsigned int numFrames,
				      unsigned* numUsableFrames,
				      bool limitTime,
				      unsigned* numPartitionsDone,
				      const bool noE, 
				      const bool cliquePosteriorNormalize,
				      const bool cliquePosteriorUnlog,
				      ObservationFile *posteriorFile)
{

  FileSource *gomFS;
  // This should be safe since gmtkOnline is the only program
  // that does inference and doesn't use FileSource and gmtkOnline
  // only uses onlineFixedUnroll
  gomFS= static_cast<FileSource *>(globalObservationMatrix);
  assert(typeid(*globalObservationMatrix) == typeid(*gomFS));

  // Unroll, but do not use the long table array (2nd parameter is
  // false) for allocation, but get back the long table length
  // in a local variable for use in our iterator.
  unsigned totalNumberPartitions;
  {
    unsigned tmp = unroll(numFrames,ZeroTable,&totalNumberPartitions);
    gomFS->justifySegment(tmp);
    if (numUsableFrames) 
      *numUsableFrames = tmp;
    // limit scope of tmp.
  }
  if (numPartitionsDone)
    *numPartitionsDone = 0;

  // set up our iterator
  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this,totalNumberPartitions);

  init_CC_CE_rvs(inference_it);

  PartitionTables* prev_part_tab = NULL;
  PartitionTables* cur_part_tab
    = new PartitionTables(inference_it.cur_jt_partition());

  // we skip the first Co's LI separator if there is no P1
  // partition, since otherwise we'll get zero probability.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.skipLISeparator();
  // gather into the root of the current  partition
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   *cur_part_tab,
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i());
  if (sectionDoDist) {
    deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		       *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());
  }
  
  // possibly print the P or C partition information
  if (inference_it.cur_part_clique_print_range() != NULL)
    printAllCliques(partitionStructureArray[inference_it.ps_i()],
		    *cur_part_tab,
		    inference_it.pt_i(),
		    inference_it.cur_nm(),
		    inference_it.cur_part_clique_print_range(),
		    stdout,
		    cliquePosteriorNormalize,cliquePosteriorUnlog,
		    false, posteriorFile);
  // if the LI separator was turned off, we need to turn it back on.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.useLISeparator();

  for (unsigned part=1;part < inference_it.pt_len(); part ++ ) {
    delete prev_part_tab;
    prev_part_tab = cur_part_tab;

    setCurrentInferenceShiftTo(part);
    cur_part_tab
      = new PartitionTables(inference_it.cur_jt_partition());

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[inference_it.ps_prev_i()],
			  *prev_part_tab,
			  inference_it.prev_ri(),
			  inference_it.prev_nm(),
			  inference_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[inference_it.ps_i()],
			  *cur_part_tab,
			  inference_it.cur_li(),
			  inference_it.cur_nm(),
			  inference_it.pt_i());

    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // it might be that E is the first partition as well, say if this is
    // a static graph, and in this case we need in this case to skip the
    // incomming separator, which doesn't exist.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();
    if (!(inference_it.at_e() && noE)) {  
      // we only do this if we're either not at an E, or if we are at
      // an E and noE is false.

      // next, gather into the root of the final E partition
      ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		       *cur_part_tab,
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());
      if (sectionDoDist) {
	deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
			   *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
			   inference_it.cur_ri(),
			   inference_it.cur_message_order(),
			   inference_it.cur_nm(),
			   inference_it.pt_i());
      }

      // possibly print the P or C partition information
      if (inference_it.cur_part_clique_print_range() != NULL)
	printAllCliques(partitionStructureArray[inference_it.ps_i()],
			*cur_part_tab,
			inference_it.pt_i(),
			inference_it.cur_nm(),
			inference_it.cur_part_clique_print_range(),
			stdout,
			cliquePosteriorNormalize, cliquePosteriorUnlog,
			false, posteriorFile);
    }
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();

    // if the LI separator was turned off, we need to turn it back on.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();


    if (limitTime && probEvidenceTimeExpired)
      goto finished;
  }
  assert ( inference_it.at_e() );

 finished:

  logpr rc;
  if (inference_it.at_e()) {
    // then we finished.
    rc = cur_part_tab->maxCliques[E_root_clique].sumProbabilities();
  }
  if (numPartitionsDone)
    *numPartitionsDone = inference_it.pt_i();

  delete cur_part_tab;

  return rc;


#if 0

  // Unroll, but do not use the long table array (2nd parameter is
  // false) for allocation, but get back the long table length
  // in a local variable for use in our iterator.
  unsigned totalNumberPartitions;
  {
    unsigned tmp = unroll(numFrames,ShortTable,&totalNumberPartitions);
    if (numUsableFrames) 
      *numUsableFrames = tmp;
    // limit scope of tmp.
  }
  if (numPartitionsDone)
    *numPartitionsDone = 0;

  // set up our iterator
  ptps_iterator ptps_it(*this,totalNumberPartitions);
  ptps_it.set_to_first_entry();

  // the set of rvs that need to be adjusted, if any.
  set <RV*> two_part_rvs;
  if (ptps_it.num_c_partitions() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // partitions (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = partitionStructureArray[1].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,two_part_rvs);
  } 

  // adjustment indices to swap the two C tables for
  // reuse. These are either zero (when we have
  // less than 2 C partitions, or swap between 1 and -1.
  // the only effec the table indices.
  int prev_i_adjust = 0;
  int cur_i_adjust = 0;

  while (!ptps_it.at_last_entry()) {

    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (ptps_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // gather into the root of the current  partition
    ceGatherIntoRoot(partitionStructureArray[ptps_it.ps_i()],
		     partitionTableArray[ptps_it.ps_i()+cur_i_adjust],
		     ptps_it.cur_ri(),
		     ptps_it.cur_message_order(),
		     ptps_it.cur_nm(),
		     ptps_it.pt_i());

    // possibly print the P or C partition information
    if (ptps_it.cur_part_clique_print_range() != NULL)
      printAllCliques(partitionStructureArray[ptps_it.ps_i()],
		      partitionTableArray[ptps_it.ps_i()+cur_i_adjust],
		      ptps_it.pt_i(),
		      ptps_it.cur_nm(),
		      ptps_it.cur_part_clique_print_range(),
		      stdout,
		      normalizePrintedCliques);

    if (limitTime && probEvidenceTimeExpired)
      goto finished;

    // if the LI separator was turned off, we need to turn it back on.
    if (ptps_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();

    // advance to next partition
    ++ptps_it;

    if (ptps_it.num_c_partitions() > 2) {
      // then we need to do some RV shifting.
      if (ptps_it.at_e()) {
	// we need to shift the two C1',C2' parts back to normal, and then shift
	// the (C,E) part forward.

	// First, shift C1',C2' rvs back to normal frame
	adjustFramesBy (two_part_rvs, -(ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

	set <RV*> tmp1,tmp2;
	tmp1 = partitionStructureArray[2].returnRVsAndTheirObservedParentsAsSet();
	tmp2 = partitionStructureArray[3].returnRVsAndTheirObservedParentsAsSet();
	// We're done now with two_part_rvs, so we make a new one. The
	// next routine clears out the original set values.
	unionRVs(tmp1,tmp2,two_part_rvs);
	// now adjust C2',E' variables to match.
	adjustFramesBy (two_part_rvs, (ptps_it.num_c_partitions()-2)
			*gm_template.S*fp.numFramesInC());

	// previous adjust needs to be what current just was
	prev_i_adjust = cur_i_adjust;
	// and turn off the current one since current is now an E
	cur_i_adjust = 0;
      } else if (ptps_it.pt_i() > 2) {

	assert ( ptps_it.at_c() );

	// shift the two C partitions over by one more
	adjustFramesBy (two_part_rvs, 1*gm_template.S*fp.numFramesInC());

	// swap settings
	if (cur_i_adjust == 0) {
	  cur_i_adjust = -1;
	  prev_i_adjust = 1;
	} else {
	  cur_i_adjust = 0;
	  prev_i_adjust = 0;
	}
      }
    }

    partitionTableArray[ptps_it.ps_i()+cur_i_adjust].init(partitionStructureArray[ptps_it.ps_i()]);

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[ptps_it.ps_prev_i()],
			  partitionTableArray[ptps_it.ps_prev_i()+prev_i_adjust],
			  ptps_it.prev_ri(),
			  ptps_it.prev_nm(),
			  ptps_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[ptps_it.ps_i()],
			  partitionTableArray[ptps_it.ps_i()+cur_i_adjust],
			  ptps_it.cur_li(),
			  ptps_it.cur_nm(),
			  ptps_it.pt_i());


  }

  assert ( ptps_it.at_e() );

  // it might be that E is the first partition as well, say if this is
  // a static graph, and in this case we need in this case to skip the
  // incomming separator, which doesn't exist.
  if (!ptps_it.has_c_partition() && P1.cliques.size() == 0)
    E1.skipLISeparator();
  // next, gather into the root of the final E partition
  if (!noE) 
    ceGatherIntoRoot(partitionStructureArray[ptps_it.ps_i()],
		     partitionTableArray[ptps_it.ps_i()],
		     ptps_it.cur_ri(),
		     ptps_it.cur_message_order(),
		     ptps_it.cur_nm(),
		     ptps_it.pt_i());
  if (!ptps_it.has_c_partition() && P1.cliques.size() == 0)
    E1.useLISeparator();

  // possibly print E partition clique info.
  if (ptps_it.cur_part_clique_print_range() != NULL)
      printAllCliques(partitionStructureArray[ptps_it.ps_i()],
		      partitionTableArray[ptps_it.ps_i()],
		      ptps_it.pt_i(),
		      ptps_it.cur_nm(),
		      ptps_it.cur_part_clique_print_range(),
		      stdout,
		      normalizePrintedCliques);

  // root clique of last partition did not do partition, since it
  // never sent to next separator (since there is none). We explicitly
  // call pruning on the root clique of the last partition.
  // curPart->maxCliques[E_root_clique].ceDoAllPruning();

 finished:

  logpr rc = partitionTableArray[ptps_it.ps_i()].maxCliques[E_root_clique].sumProbabilities();
  if (numPartitionsDone)
    *numPartitionsDone = ptps_it.pt_i();


  // lastly, shift variables back for distribute evidence.
  adjustFramesBy (two_part_rvs, -(ptps_it.num_c_partitions()-2)
		  *gm_template.S*fp.numFramesInC());

  return rc;

#endif

}


void
JunctionTree::enqueueMoreFrames(unsigned &truePtLen, unsigned part, StreamSource *globalObservationMatrix, 
				unsigned numNewFrames, unsigned *numUsableFrames)
{
  // ticket #468 - skip enqueue on first iteration if P is empty to keep PCCE first
  // C observation data available in the queue
  if (truePtLen == 0 && (part > 1 || fp.numFramesInP() > 0)) {  // if truePtLen is known, no need to read more frames
    // read in the same # of frames that we're about to consume to maintain
    // enough queued frames to be sure we don't overshoot the C'->E' transition
    (void) globalObservationMatrix->enqueueFrames(numNewFrames);
  }
  // update the ptps_iterator if we just found out the true length of the segment
  if (truePtLen == 0 && globalObservationMatrix->numFrames() != 0) {
    unsigned basicTempMaxUnrollAmnt;
    unsigned basicTempMinUnrollAmnt;
    int      modTempMaxUnrollAmnt;
    int      modTempMinUnrollAmnt;
    unsigned numUsableFrm;
    unsigned frmStart;
    unsigned T;
    
    if ( globalObservationMatrix->numFrames() <= globalObservationMatrix->minFutureFrames() )
      error("Segment of %d frames too short as model requires at least %d frames\n", 
	    globalObservationMatrix->numFrames() + globalObservationMatrix->minPastFrames(), 
	    globalObservationMatrix->minPastFrames() + fp.numFramesInP() + fp.numFramesInE() + globalObservationMatrix->minFutureFrames());
    
    T  = globalObservationMatrix->numFrames() - globalObservationMatrix->minFutureFrames();
    if (!gm_template.computeUnrollParameters(T,
					     basicTempMaxUnrollAmnt,
					     basicTempMinUnrollAmnt,
					     modTempMaxUnrollAmnt,
					     modTempMinUnrollAmnt,
					     numUsableFrm,
					     frmStart))
      error("Segment of %d frames too short with current GMTK template of length [P=%d,C=%d,E=%d] %d frames, and M=%d,S=%d boundary parameters. "
	    "Use longer utterances, different template, or decrease M,S if >1.\n",
	    T,
	    fp.numFramesInP(),fp.numFramesInC(),fp.numFramesInE(),
	    fp.numFrames(),
	    gm_template.M,gm_template.S);
    
    truePtLen = modTempMaxUnrollAmnt + 3;
    inference_it.set_pt_len(truePtLen);
    if (numUsableFrames) 
      *numUsableFrames = numUsableFrm;
    
  } // if just discovered segment length
}


// DBN online filtering/smoothing
logpr 
JunctionTree::onlineFixedUnroll(StreamSource *globalObservationMatrix,
				unsigned *numUsableFrames,
				unsigned *numPartitionsDone,
				const bool noE,
				FILE *f,
				const bool printObserved,
				regex_t *preg,
				regex_t *creg,
				regex_t *ereg,
				char *partRangeFilter,
				ObservationFile *posteriorFile,
				const bool cliquePosteriorNormalize,
				const bool cliquePosteriorUnlog)
{

  if (fp.numFramesInC() == 0) {
    error("ERROR: gmtkOnline does not support empty chunks\n");
  }

  if (!f) f = stdout; // just for console progress logging

  // Unroll, but do not use the long table array (2nd parameter is
  // false) for allocation, but get back the long table length
  // in a local variable for use in our iterator.
  unsigned totalNumberPartitions;

  // We don't know how many frames are in the segment yet,
  // so just assume it's very long. Should be OK, since the
  // ZeroTable option just uses the minimum number of partitions.
  // We'll reset the ptps_iterator's partition length later
  // when we actually know it.

#if 1
#  define FRAME_NUMBER_LIMIT (1073741824U)
#else
  // to force frequent overflows for testing
#  define FRAME_NUMBER_LIMIT (60U)
#endif

  (void) unroll(FRAME_NUMBER_LIMIT, ZeroTable, &totalNumberPartitions);


  // Try to read in enough frames for P' C'_1 C'_2 ... C'_\tau
  //                               or P' C'_1 C'_2 ... C'_{\tau-1} E'
  unsigned M = gm_template.M;
  unsigned S = gm_template.S;

  // Here |X| is the number of frames in partition X.

  // Assuming the segment is infinitely long, in the general case
  // we need |P| + \tau s|C| + \delta_C' frames. 

  // \delta_C' is the maximum # of frames into the C'_t, C'_{t+1} frame overlap
  // in which C'_t contains an observed variable. That's how many
  // frames beyond |P| + \tau s |C| we must read before we can compute the 
  // posterior of P'.

  // Let t_0 be the frame # of the latest observed variable in C'_1
  // \delta_C' = max(0, t_0 - |P| - s|C| + 1) <= m|C|

  // If |E|=0, we can safely treat the last section as C', so in that case
  // we don't need to read any extra frames to handle finite length segments.
  // If E is not empty, we need to preload |P| + \tau s |C| + m|C| + |E| frames
  // to detect the C' to E' transition in time for inference to work correctly.

  // We also need to preload startSkip frames so that they can be skipped :)
  // and the max Dlink lag so that the last modified section has sufficient future available.

  // GMParms has checked that -startSkip is big enough

  // The StreamSource doesn't implement -endSkip since the
  // stream's length is unknown, but ObservationSource::minFutureFrames()
  // has the correct number of frames to skip at the end

  unsigned extraInterfaceFrames; // \delta_{C'}
  {
    if (0 == fp.numFramesInE()) {
      unsigned maxObservedCPrime1Frame = 0; // frame # of latest observed variable in C'_1
      unsigned lastFrameBeforeCPrime2 = fp.numFramesInP() + S * fp.numFramesInC() - 1; // last frame # before C'_2
      vector<RV*> C = partitionStructureArray[1].allrvs_vec;
      for (unsigned i=0; i < C.size(); i+=1) {
	if (C[i]->observed() && C[i]->frame() > maxObservedCPrime1Frame) {
	  maxObservedCPrime1Frame = C[i]->frame();
	}
      }
      infoMsg(IM::ObsStream, IM::Med+1, "frame # of latest observed variable in C'_1 is %u\n", maxObservedCPrime1Frame);
      if (maxObservedCPrime1Frame > lastFrameBeforeCPrime2) {
	extraInterfaceFrames = maxObservedCPrime1Frame - lastFrameBeforeCPrime2;
      } else {
	extraInterfaceFrames = 0;
      }
      assert(extraInterfaceFrames <= M * fp.numFramesInC());
      infoMsg(IM::ObsStream, IM::Med, "Adding %u extra preload frames due to observed variables in C'C' interface\n", extraInterfaceFrames);
    } else {
      extraInterfaceFrames = M * fp.numFramesInC();
      infoMsg(IM::ObsStream, IM::Med, "Adding %u extra preload frames to detect C' to E' transition\n", extraInterfaceFrames + fp.numFramesInE());
    }
  }

  unsigned tau = numSmoothingPartitions; // # of "future" C's for smoothing

  unsigned numPreloadFrames = 
    globalObservationMatrix->startSkip() + 
    fp.numFramesInP() + 
    tau * S * fp.numFramesInC() + 
    fp.numFramesInE() + 
    extraInterfaceFrames +
    globalObservationMatrix->minFutureFrames();
  // Assume the above won't over-flow the frame # with so few frames
  assert(FRAME_NUMBER_LIMIT > numPreloadFrames);

  // setCurrentInferenceShiftTo(C'_t) may need the frames for C'_{t-1}: (S+M)|C| extra frames
  // Note that these don't count towards the lag - we're just making extra space in the
  // buffer to hold the old C' frames, we're not preloading them.
  globalObservationMatrix->setActiveFrameCount(numPreloadFrames + (S+M) * fp.numFramesInC());

  bool overflow = false; // true iff we had to reset the frame # due to more than ~ 2^30 frames

  // NOTE: we assume \tau << 2^30  i.e. if the frame # had to overflow, we've definitely
  // processed more than \tau partitions, and we can instantiate the first 2(\tau + 1) partitions
  // without the frame number overflowing

  unsigned numNewFrames = fp.numFramesInC() * S;  // read in this many frames for the next C'

  infoMsg(IM::ObsStream, IM::Info, "preloading %u frames to smooth with %u modified sections\n", numPreloadFrames, tau);
  globalObservationMatrix->preloadFrames(numPreloadFrames);

  if (globalObservationMatrix->EOS()) return logpr(0.0);

  fprintf(f,"========\nSegment %u\n", globalObservationMatrix->segmentNumber());

  unsigned truePtLen = 0; // 0 until we know the true number of modified partitions
  {
    unsigned T = globalObservationMatrix->numFrames();

    bool rememberedViterbiScore = viterbiScore;
    viterbiScore = false; // avoid allocating space for O(T) viterbi values in unroll()

    if (T > 0) {
      if (T <= globalObservationMatrix->minFutureFrames()) {
	warning("WARNING: observation stream segment of %u frames is too short. Model requires at least %u frames\n", 
		T,  globalObservationMatrix->minFutureFrames());
	return logpr(0.0);
      }
      T -= globalObservationMatrix->minFutureFrames(); // fake -endSkip
    }
    infoMsg(IM::Inference, IM::Info, "onlineFixedUnroll: total # partitions %u\n", totalNumberPartitions);
    
    viterbiScore = rememberedViterbiScore;  // do compute viterbi values in deScatterOutofRoot()? (max-product semiring)
  }
  if (numPartitionsDone)
    *numPartitionsDone = 0;
  

  // output printing state variables

  set<string> variableNames; // names of variables in the model
  for (unsigned i=0; i < partition_unrolled_rvs.size(); i+=1) {
    variableNames.insert(partition_unrolled_rvs[i]->name());
  }

  RVVec  pVitTriggerVec;
  string pVitTriggerExpr;
  RngDecisionTree::EquationClass pTriggerEqn;
  initializeViterbiTrigger(pVitTrigger, variableNames, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, 'p');

  RVVec  cVitTriggerVec;
  string cVitTriggerExpr;
  RngDecisionTree::EquationClass cTriggerEqn;
  initializeViterbiTrigger(cVitTrigger, variableNames, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, 'c');

  RVVec  eVitTriggerVec;
  string eVitTriggerExpr;
  RngDecisionTree::EquationClass eTriggerEqn;
  initializeViterbiTrigger(eVitTrigger, variableNames, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, 'e');

  bool first_P = true;
  bool first_C = true;
  bool first_E = true;

  unsigned P_size = 0;
  unsigned C_size = 0;
  unsigned E_size = 0;

  sArray<unsigned> previous_P_values;
  sArray<unsigned> previous_C_values;
  sArray<unsigned> previous_E_values;

  vector<bool> pregex_mask;
  vector<bool> cregex_mask;
  vector<bool> eregex_mask;
  
  // FIXME: it looks like we do support writing posteriors to files - this should be symmetric...
  FILE *vitFile = NULL; // I don't think we support writing output to files since the size is unbounded


  unsigned numBufferedPartitions = tau > 0 ? tau + 1 : 2;
  vector<PartitionTables *> partitionBuffer(numBufferedPartitions, NULL);


  // C' (or E') must end at frame <= FRAME_NUMBER_LIMIT to avoid overflow
  
  // We also want to only trigger overflow after processing a section such that
  //    part % (\tau + 1) = \tau - 1 
  // so that the new section lands in the last slot of the section buffer and doesn't 
  // stomp on the \tau previous sections
  
  // The i^th C' (i starts at 1 since P' is 0) 
  //                       begins at frame |P| + (i-1)S|C|                     and 
  //                         ends at frame |P| +   i  S|C| +       M|C| - 1
  // An E' following i-1 C's ends at frame |P| + (i-1)S|C| + |E| + M|C| - 1

  // So we want to find the modified section with largest last frame # x <= FRAME_NUMBER_LIMIT such that 
  //      x = |P| + (i-1)S|C| + max(S|C|, |E|) + M|C| - 1
  // and  i % (\tau + 1) = \tau - 1                            (only if we're smoothing, ie \tau > 0)

  // HOWEVER, the above analysis is irrelevant, since the preloaded frames' numbers are
  // running ahead of the last frame of the section we're currently running inference on...
  // The difference between the frame # inference is working on and the last frame # in the
  // active area of the frame queue can be up to globalObservationMatrix->getActiveFrameCount()
  assert(FRAME_NUMBER_LIMIT > globalObservationMatrix->getActiveFrameCount());
  unsigned partition_number_limit = (FRAME_NUMBER_LIMIT - globalObservationMatrix->getActiveFrameCount() ) / (S * fp.numFramesInC());
 
  if (numSmoothingPartitions > 0) {
    if ( partition_number_limit % numBufferedPartitions != (tau-1) ) {
      partition_number_limit = ( partition_number_limit / numBufferedPartitions - 1 ) * numBufferedPartitions + numSmoothingPartitions - 1;
    }
    //    assert(partition_number_limit > 2*tau+1); // resetting to C'_{2\tau+1}, so limit must be larger than that
    assert(partition_number_limit > 3*tau+2); // resetting to C'_{3\tau+2}, so limit must be larger than that
    assert(partition_number_limit % numBufferedPartitions == tau - 1);
  } else {
    // overflow C'_\theta to C'_2
    // resetting to C'_2, so limit must be > 2. Might need to subtract 1, so must be > 3. Must be even, so must be > 4.
    assert(partition_number_limit > 4); 
    if ( partition_number_limit % 2 != 0 ) {
      partition_number_limit -= 1;
    }
    assert(partition_number_limit > 2);
    assert(partition_number_limit % 2 == 0);
  }
  // Let \theta = partition_mumber_limit. Instead of incrementing from 
  // C'_{\theta-1} to C'_{theta), we're going to wrap section C'_{\theta} to C'_{2\tau + 1}
  // and renumber the frames to match. We go to 2\tau+1 because if we went
  // to C'_{\tau}, we'd try to print P' again during the smoothing pass...

  // Set up our iterator, write over the member island iterator since
  // we assume the member does not have any dynamc sub-members.
  new (&inference_it) ptps_iterator(*this,totalNumberPartitions);

  init_CC_CE_rvs(inference_it);

  PartitionTables* prev_part_tab = NULL;
  PartitionTables* cur_part_tab = new PartitionTables(inference_it.cur_jt_partition());
  partitionBuffer[0] = cur_part_tab;

  // we skip the first Co's LI separator if there is no P1
  // partition, since otherwise we'll get zero probability.
  if (inference_it.at_first_c() && P1.cliques.size() == 0)
    Co.skipLISeparator();

  // it might be that E is the first partition as well, say if this is
  // a static graph, and in this case we need in this case to skip the
  // incomming separator, which doesn't exist.
  if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
    E1.skipLISeparator();

  // gather into the root of the first partition
  ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		   *cur_part_tab,
		   inference_it.cur_ri(),
		   inference_it.cur_message_order(),
		   inference_it.cur_nm(),
		   inference_it.pt_i());

  if (numSmoothingPartitions == 0) {
    // must do scatter & print for P' here since it won't happen in the smoothing loop
    PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
    if (viterbiScore && P1.cliques.size() > 0 ) { // nothing to set if P is empty - ticket #468
      cur_part_tab->maxCliques[inference_it.cur_ri()].maxProbability(ps.maxCliquesSharedStructure[inference_it.cur_ri()], true);
    }
    deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		       *cur_part_tab,
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());
    

    // if the LI separator was turned off, we need to turn it back on.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();
    
    // print P'
    if (viterbiScore) {
      if (inference_it.at_p()) {
	printUnpackedSection(ps, pVitTrigger!=NULL, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, printObserved, inference_it.pt_i(),
			     'P', f, preg, pregex_mask, first_P, P_size, previous_P_values);
      } else if (inference_it.at_e()) {
	printUnpackedSection(ps, eVitTrigger!=NULL, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, printObserved, inference_it.pt_i(),
			     'E', f, ereg, eregex_mask, first_E, E_size, previous_E_values);
      } else {
	assert ( inference_it.at_c() );      
	printUnpackedSection(ps, cVitTrigger!=NULL, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, printObserved, inference_it.pt_i(),
			     'C', f, creg, cregex_mask, first_C, C_size, previous_C_values, vitRunLength,   vitFile ? 1 : inference_it.pt_i());
      }
    } else {
      // possibly print the P or C partition information
      if (inference_it.cur_part_clique_print_range() != NULL) {
	printAllCliques(partitionStructureArray[inference_it.ps_i()],
			*cur_part_tab,
			inference_it.pt_i(),
			inference_it.cur_nm(),
			inference_it.cur_part_clique_print_range(),
			stdout,
			cliquePosteriorNormalize, cliquePosteriorUnlog,
			false, posteriorFile);			
      }
    } 
    enqueueMoreFrames(truePtLen, 0, globalObservationMatrix, numNewFrames, numUsableFrames);
  } // handle scatter & print P' if filtering
  fflush(f); if (!posteriorFile) fflush(stdout);
  
  unsigned part, cur_part_idx, prev_part_idx;
  for (part = 1; part < inference_it.pt_len(); part += 1) {  // main inference loop over sections
    cur_part_idx = part % numBufferedPartitions;
    prev_part_idx = (cur_part_idx + numBufferedPartitions - 1) % numBufferedPartitions;
    delete partitionBuffer[cur_part_idx];
    
    setCurrentInferenceShiftTo(part);
    cur_part_tab = new PartitionTables(inference_it.cur_jt_partition());
    prev_part_tab = partitionBuffer[prev_part_idx];
    partitionBuffer[cur_part_idx] = cur_part_tab;

    // send from previous to current
    ceSendForwardsCrossPartitions(// previous partition
			  partitionStructureArray[inference_it.ps_prev_i()],
			  *prev_part_tab,
			  inference_it.prev_ri(),
			  inference_it.prev_nm(),
			  inference_it.pt_prev_i(),
			  // current partition
			  partitionStructureArray[inference_it.ps_i()],
			  *cur_part_tab,
			  inference_it.cur_li(),
			  inference_it.cur_nm(),
			  inference_it.pt_i());

    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();

    // it might be that E is the first partition as well, say if this is
    // a static graph, and in this case we need in this case to skip the
    // incomming separator, which doesn't exist.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();

    // next, gather into the root of the current partition
    ceGatherIntoRoot(partitionStructureArray[inference_it.ps_i()],
		     *cur_part_tab,
		     inference_it.cur_ri(),
		     inference_it.cur_message_order(),
		     inference_it.cur_nm(),
		     inference_it.pt_i());

    // if the LI separator was turned off, we need to turn it back on.
    if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();


    if (part >= numSmoothingPartitions || overflow) {    // now smooth if there's enough future

      if (viterbiScore) {
	// The current partition is acting as E' - need to set its back_max_cvn here as a precondition
	//   to the smoothing DE pass...
	PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
	cur_part_tab->maxCliques[inference_it.cur_ri()].maxProbability(ps.maxCliquesSharedStructure[inference_it.cur_ri()], true);
      }

      for (unsigned i=1; i <= numSmoothingPartitions; i+=1) {    // length \tau backwards pass (smoothing)

	// skip unconnected separators
	if (inference_it.at_first_c() && P1.cliques.size() == 0)
	  Co.skipLISeparator();    
	else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
	  E1.skipLISeparator();

	deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
			   *cur_part_tab,
			   inference_it.cur_ri(),
			   inference_it.cur_message_order(),
			   inference_it.cur_nm(),
			   inference_it.pt_i());
 
	// restore unconnected separators
	if (inference_it.at_first_c() && P1.cliques.size() == 0)
	  Co.useLISeparator();
	else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
	  E1.useLISeparator();
	

	// send backwads message to previous partition
	deSendBackwardsCrossPartitions(partitionStructureArray[inference_it.ps_prev_i()],
				       *prev_part_tab,
				       inference_it.prev_ri(),
				       inference_it.prev_nm(),
				       inference_it.pt_prev_i(),
				       //
				       partitionStructureArray[inference_it.ps_i()],
				       *cur_part_tab,
				       inference_it.cur_li(),
				       inference_it.cur_nm(),
				       inference_it.pt_i());
	
	setCurrentInferenceShiftTo(part-i);
	cur_part_tab  = prev_part_tab;
	prev_part_tab = partitionBuffer[(part + numBufferedPartitions - 1 - i) % numBufferedPartitions];
      }  // length \tau backward pass if smoothing

      // scatter out of root for the printing partition
      
      // skip unconnected separators
      if (inference_it.at_first_c() && P1.cliques.size() == 0)
	Co.skipLISeparator();    
      else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
	E1.skipLISeparator();
      
      deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
			 *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
			 inference_it.cur_ri(),
			 inference_it.cur_message_order(),
			 inference_it.cur_nm(),
			 inference_it.pt_i());
      
      // restore unconnected separators
      if (inference_it.at_first_c() && P1.cliques.size() == 0)
	Co.useLISeparator();
      else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
	E1.useLISeparator();
      
      // print part-\tau

      PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
 
      if (viterbiScore) {
	if (inference_it.at_p()) {
	  printUnpackedSection(ps, pVitTrigger!=NULL, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, printObserved, inference_it.pt_i(),
			       'P', f, preg, pregex_mask, first_P, P_size, previous_P_values);
	} else if (inference_it.at_e()) {
	  printUnpackedSection(ps, eVitTrigger!=NULL, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, printObserved, inference_it.pt_i(),
			       'E', f, ereg, eregex_mask, first_E, E_size, previous_E_values);
	} else {
	  assert ( inference_it.at_c() );      
	  printUnpackedSection(ps, cVitTrigger!=NULL, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, printObserved, inference_it.pt_i(),
			       'C', f, creg, cregex_mask, first_C, C_size, previous_C_values, vitRunLength,   vitFile ? 1 : inference_it.pt_i());
	}
      } else {
	// possibly print the P or C partition information
	if (inference_it.cur_part_clique_print_range() != NULL) {
	  printAllCliques(partitionStructureArray[inference_it.ps_i()],
			  *cur_part_tab,
			  inference_it.pt_i(),
			  inference_it.cur_nm(),
			  inference_it.cur_part_clique_print_range(),
			  stdout,
			  cliquePosteriorNormalize, cliquePosteriorUnlog,
			  false, posteriorFile);			
	}
      }
      fflush(f); if (!posteriorFile) fflush(stdout);

   
      // equeue more frames

      if (part == partition_number_limit && truePtLen == 0) {  // frame number is about to overflow
	// There's no point trying to reset the frame #'s when the segment length is known,
	// since all of the frames are already in the queue, presumably with non-overflowed numbers.
	overflow = true;
	
	if (numSmoothingPartitions > 0) {
	  assert( (part % numBufferedPartitions) == (tau - 1) );
	  //	  assert( part > 2 * tau + 1 );
	  assert(part > 3 * tau + 2);
	  unsigned firstFrameOfCurrentPart = fp.numFramesInP() + (part-1) * S * fp.numFramesInC();
	  assert( globalObservationMatrix->firstFrameInQueue() <= firstFrameOfCurrentPart && firstFrameOfCurrentPart <= globalObservationMatrix->lastFrameInQueue() );
	  unsigned firstFrameOfResetPart = fp.numFramesInP() + (3 * tau) * S * fp.numFramesInC(); // current section's first frame # after reset to section C'_{3 \tau+ 1}
	  assert( firstFrameOfCurrentPart > firstFrameOfResetPart );
	  unsigned delta = firstFrameOfCurrentPart - firstFrameOfResetPart;
	  assert( globalObservationMatrix->firstFrameInQueue() >= delta);
	  unsigned firstFrame = globalObservationMatrix->firstFrameInQueue() - delta;
	  infoMsg(IM::Inference, IM::Info, "resetting first queued frame number to %u\n", firstFrame);
	  globalObservationMatrix->resetFrameNumbers(firstFrame);
	  infoMsg(IM::Inference, IM::Info, "resetting ptps to partition %u to %u\n", part, 2 * numSmoothingPartitions);
	  //	  part = 2 * numSmoothingPartitions; // continue @ partition 2 \tau + 1, about to increment part at top of loop
	  part = 3 * tau + 1; // continue @ section 3 \tau + 2, about to increment part at top of loop
	  // don't reset to partition \tau or it'll print P' on every overflow
	  // don't reset to partition 2\tau+1 or reset frame #s may underflow for small \tau
	} else {
	  assert( part % 2 == 0 );
	  assert( part > 3); // P' C'_1 can't be reused, resetting to C'_2, C'_3, so current section must be later than that
	  unsigned firstFrameOfCurrentPart = fp.numFramesInP() + (part-1) * S * fp.numFramesInC();
	  assert( globalObservationMatrix->firstFrameInQueue() <= firstFrameOfCurrentPart && firstFrameOfCurrentPart <= globalObservationMatrix->lastFrameInQueue() );
	  unsigned firstFrameOfResetPart = fp.numFramesInP() + S * fp.numFramesInC(); // current section's first frame # after reset to section C'_2
	  assert( firstFrameOfCurrentPart > firstFrameOfResetPart );
	  unsigned delta = firstFrameOfCurrentPart - firstFrameOfResetPart;
	  assert( globalObservationMatrix->firstFrameInQueue() >= delta);
	  unsigned firstFrame = globalObservationMatrix->firstFrameInQueue() - delta;
	  infoMsg(IM::Inference, IM::Info, "resetting first queued frame number to %u\n", firstFrame);
	  globalObservationMatrix->resetFrameNumbers(firstFrame);
	  infoMsg(IM::Inference, IM::Info, "resetting ptps section %u to 3\n", part);
	  part = 2; // continue @ section 3, about to increment part at top of loop
	}
      }

      enqueueMoreFrames(truePtLen, part, globalObservationMatrix, numNewFrames, numUsableFrames);

    } // if handled the first \tau partitions
        
  } // loop over partitions

  if (part <= numSmoothingPartitions && !overflow) { 
    // do normal DE pass, since there weren't enough frames for smoothing

    unsigned shortPart = part-1; // partition # we're handling in the short segment

    if (viterbiScore) {
      // The E' didn't get maxProb'd yet because it's a short segment
      setCurrentInferenceShiftTo(shortPart);

      PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];
      cur_part_tab->maxCliques[inference_it.cur_ri()].maxProbability(ps.maxCliquesSharedStructure[inference_it.cur_ri()], true);
    }

    for ( ; shortPart > 0; --shortPart) {

      setCurrentInferenceShiftTo(shortPart);
      cur_part_tab = partitionBuffer[shortPart];
      prev_part_tab = partitionBuffer[shortPart-1];
      
      // skip unconnected separators
      if (inference_it.at_first_c() && P1.cliques.size() == 0)
	Co.skipLISeparator();    
      else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
	E1.skipLISeparator();

      deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
			 *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
			 inference_it.cur_ri(),
			 inference_it.cur_message_order(),
			 inference_it.cur_nm(),
			 inference_it.pt_i());
      
      // restore unconnected separators
      if (inference_it.at_first_c() && P1.cliques.size() == 0)
	Co.useLISeparator();
      else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
	E1.useLISeparator();
      
      // send backwads message to previous partition
      deSendBackwardsCrossPartitions(partitionStructureArray[inference_it.ps_prev_i()],
				     *prev_part_tab,
				     inference_it.prev_ri(),
				     inference_it.prev_nm(),
				     inference_it.pt_prev_i(),
				     //
				     partitionStructureArray[inference_it.ps_i()],
				     *cur_part_tab,
				     inference_it.cur_li(),
				     inference_it.cur_nm(),
				     inference_it.pt_i());
    }  // DE pass
    
    // DE pass P'
    setCurrentInferenceShiftTo(shortPart);
    cur_part_tab = partitionBuffer[shortPart];
      
    // skip unconnected separators
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.skipLISeparator();    
    else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.skipLISeparator();

    deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
		       *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
		       inference_it.cur_ri(),
		       inference_it.cur_message_order(),
		       inference_it.cur_nm(),
		       inference_it.pt_i());
    
    // restore unconnected separators
    if (inference_it.at_first_c() && P1.cliques.size() == 0)
      Co.useLISeparator();
    else if (!inference_it.has_c_partition() && P1.cliques.size() == 0)
      E1.useLISeparator();

  } // short segment - normal DE pass


  // print left-overs: the last \tau partitions, or all of them if T < \tau + 1
  
  if (numSmoothingPartitions > 0) { // none left-over if filtering (\tau = 0)
    unsigned offset = ( part > numSmoothingPartitions+1 || overflow )  ?  numSmoothingPartitions : part;
    unsigned unprintedPartitionCount  = ( part > numSmoothingPartitions+1 || overflow )  ?  numSmoothingPartitions : part; 
    for (unsigned i=0; i < unprintedPartitionCount; i+=1) {
      setCurrentInferenceShiftTo(part - offset + i);
      cur_part_idx = (part - offset + i) % numBufferedPartitions;
      cur_part_tab = partitionBuffer[cur_part_idx];
      
      PartitionStructures& ps = partitionStructureArray[inference_it.ps_i()];

      // Why is the deScatterOutofRoot() here? The C's printed during the main loop are printed immediately
      // after their deScatterOutofRoot(), so the C' RVs are left in the max probability values. However,
      // this code is printing the last \tau C's, which have experienced a "normal" backwards pass - 
      // C'[ T-\tau] will be the last C' to have been deScatterOutofRoot()ed, so the C' shared structure RVs
      // will be "stuck" at that partition's max probability values for the last \tau C's. Calling 
      // deScatterOutofRoot() here right before printing the C' ensures that it is set to the correct values.
      // This does mean that the last \tau partitions get scattered twice. This could possibly be avoided 
      // by having a circular buffer of length \tau+1 for the packed RV values, but that would add work for
      // every partition instead of just the last \tau.
      deScatterOutofRoot(partitionStructureArray[inference_it.ps_i()],
			 *cur_part_tab, //partitionTableArray[inference_it.pt_i()],
			 inference_it.cur_ri(),
			 inference_it.cur_message_order(),
			 inference_it.cur_nm(),
			 inference_it.pt_i());
      if (viterbiScore) {
	if (inference_it.at_p()) {
	  printUnpackedSection(ps, pVitTrigger!=NULL, pVitTriggerVec, pVitTriggerExpr, pTriggerEqn, printObserved, inference_it.pt_i(),
			       'P', f, preg, pregex_mask, first_P, P_size, previous_P_values);
	} else if (inference_it.at_e()) {
	  printUnpackedSection(ps, eVitTrigger!=NULL, eVitTriggerVec, eVitTriggerExpr, eTriggerEqn, printObserved, inference_it.pt_i(),
			       'E', f, ereg, eregex_mask, first_E, E_size, previous_E_values);
	} else {
	  assert ( inference_it.at_c() );      
	  printUnpackedSection(ps, cVitTrigger!=NULL, cVitTriggerVec, cVitTriggerExpr, cTriggerEqn, printObserved, inference_it.pt_i(),
			       'C', f, creg, cregex_mask, first_C, C_size, previous_C_values, vitRunLength,   vitFile ? 1 : inference_it.pt_i());
	}
      } else {
	// possibly print the partition information
	if (inference_it.cur_part_clique_print_range() != NULL) {
	  printAllCliques(partitionStructureArray[inference_it.ps_i()],
			  *cur_part_tab,
			  inference_it.pt_i(),
			  inference_it.cur_nm(),
			  inference_it.cur_part_clique_print_range(),
			  stdout,
			  cliquePosteriorNormalize, cliquePosteriorUnlog,
			  false, posteriorFile);			
	}
      }
      fflush(f); if (!posteriorFile) fflush(stdout);
    }   // printing left-over partitions
  }

  logpr rc;
  rc = cur_part_tab->maxCliques[E_root_clique].sumProbabilities();

  // delete partitionBuffer[i]
  for (unsigned i=0; i < numBufferedPartitions; i+=1) {
    delete partitionBuffer[i];
  }

  return rc;
}



/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
