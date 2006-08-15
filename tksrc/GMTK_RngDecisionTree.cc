/*-
 * GMTK_RngDecisionTree.cc
 *     General Class to map from vectors of integers to some basic type.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu> & Chris Bartels <bartels@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
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

#include "GMTK_DiscRV.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_PackCliqueValue.h"
#include "GMTK_RngDecisionTree.h"


VCID("$Header$")


/////////////////////////////////////////////////////////////////////
// File extension for compiled DT files 
/////////////////////////////////////////////////////////////////////
const string DTFileExtension = ".index";

/////////////////////////////////////////////////////////////////////
// Static computation stack for formula evaluation 
/////////////////////////////////////////////////////////////////////
RngDecisionTree::sArrayStack<RngDecisionTree::EquationClass::stack_element_t> 
  RngDecisionTree::EquationClass::stack(0);

/////////////////////////////////////////////////////////////////////
// Arrays to classify formula tokens and map their strings to their 
// enumerations 
/////////////////////////////////////////////////////////////////////

map<string, RngDecisionTree::EquationClass::tokenEnum> 
  RngDecisionTree::EquationClass::delimiter;
map<string, RngDecisionTree::EquationClass::tokenEnum> 
  RngDecisionTree::EquationClass::function;
map<string, RngDecisionTree::EquationClass::tokenEnum> 
  RngDecisionTree::EquationClass::variable;

map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::infixToken;
map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::unaryToken;
map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::functionToken;
map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::oneValFunctionToken;
map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::twoValFunctionToken;
map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::manyValFunctionToken;
map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::variableToken;

map<RngDecisionTree::EquationClass::tokenEnum, unsigned> 
  RngDecisionTree::EquationClass::tokenPriority;

#ifdef MAIN
RAND rnd(true);
#endif 

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                    CONSTRUCTOR/DESTRUCTOR
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/*-
 *-----------------------------------------------------------------------
 * ~RngDecisionTree
 *      destructor
 * 
 * Preconditions:
 *      object must be "filled" in.
 *
 * Postconditions:
 *      object is not filled in.
 *
 * Side Effects:
 *      all memory associated with this object is reclaimed.
 *
 * Results:
 *      nothing
 *
 *-----------------------------------------------------------------------
 */
RngDecisionTree::~RngDecisionTree()
{
  if (root != NULL) {
    destructorRecurse(root);
  }
  delete root;
  delete dtFile;
}


/*-
 *-----------------------------------------------------------------------
 * destructorRecurse
 *      destructor helper
 * 
 * Preconditions:
 *      object must be "filled" in.
 *
 * Postconditions:
 *      object is not filled in.
 *
 * Side Effects:
 *      all memory associated with this object is reclaimed.
 *
 * Results:
 *      nothing
 *
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::destructorRecurse(RngDecisionTree::Node* node)
{
  if (node == NULL)
    return;


  switch (node->nodeType) {
  case NonLeafNodeArray: 
    // printf("deleting NonLeafNodeArray\n");
    for (unsigned i=0;i<node->nln_a().children.size(); i++) {
      destructorRecurse(&(node->nln_a().children[i]));
    }
    node->nln_a().children.clear();
    // node->nln_a().children.~sArray < Node >();
    node->nln_a().~NonLeafNodeArrayStruct();
    break;

  case NonLeafNodeHash:
    // printf("deleting NonLeafNodeHash\n");
    for (unsigned i=0;i<node->nln_h().children.size(); i++) {
      destructorRecurse(&(node->nln_h().children[i]));
    }
    node->nln_h().children.clear();      
    // try to call destructors
    // node->nln_h().children.~sArray < Node >();
    // node->nln_h().nodeMapper.~shash_map < unsigned, Node* >();
    node->nln_h().~NonLeafNodeHashStruct();
    break;
	
  case NonLeafNodeRngs:
    // printf("deleting NonLeafNodeRngs\n");
    for (unsigned i=0;i<node->nln_r().children.size(); i++) {
      destructorRecurse(&(node->nln_r().children[i].nd));
    }
    destructorRecurse(node->nln_r().def);
    node->nln_r().children.clear();
    delete node->nln_r().def;
    // node->nln_r().children.~sArray < RangeNode >();
    node->nln_r().~NonLeafNodeRngsStruct();
    break;

  case LeafNodeVal:
    // printf("deleting LeafNodeVal\n");
    // nothing to delete here
    node->ln_v().~LeafNodeValStruct();
    break;

  case LeafNodeEquation:
    // printf("deleting LeafNodeEquation\n");
    // make sure anything in an equation gets deleted
    node->ln_e().~LeafNodeEquationStruct();
    break;

  default:
    coredump("INTERNAL ERROR: unknown DT node type.\n");
    break;
  }

}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                    Reading/Writing
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * read:
 *      reads in from a file.
 * 
 * Preconditions:
 *      object should not be filled on (data will be lost if so).
 *
 * Postconditions:
 *      object will be filled in if this is a plain DT. This
 *      might, however, be a pointer to another file that contains
 *      a list of DTs. In that case, the object is not filled
 *      in and nextIterableDT() must be called for
 *      each DT in the file.
 *
 * Side Effects:
 *      Totally changes the object, might kill program if error occurs.
 *
 * Results:
 *      none 
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::read(iDataStreamFile& is)
{
  NamedObject::read(is);

  //////////////////////////////////////////////////////////////////////
  // Read the beginning of the next decision tree 
  //////////////////////////////////////////////////////////////////////
  is.read(dtFileName,"Can't read Decision Tree's either file or numFeatures");

  //////////////////////////////////////////////////////////////////////
  // Entry is the name of a file containing multiple instances of the 
  // decision tree.
  //////////////////////////////////////////////////////////////////////
  if (!strIsInt(dtFileName.c_str(), (int*)&_numFeatures)) {

    //////////////////////////////////////////////////////////////////////
    // Make sure this instance already hasn't been initialized as an
    // iterable decision tree 
    //////////////////////////////////////////////////////////////////////
    if (iterable()) {
      error("ERROR: in DT named '%s' in file '%s' line %d, can't have DTs defined recursively in files",
        name().c_str(),is.fileName(),is.lineNo());
    }

    initializeIterableDT(dtFileName); 
  }
  //////////////////////////////////////////////////////////////////////
  // Entry is an integer indicating an inline decision tree 
  //////////////////////////////////////////////////////////////////////
  else {

    if (_numFeatures < 0) {
      error("ERROR: in DT named '%s', file '%s' line %d, decision tree must have >= 0 features",
      name().c_str(), is.fileName(),is.lineNo());
    }

    if (root != NULL) {
      destructorRecurse(root);
      delete root;
      root = NULL;
    }
    root = new Node;
    readRecurse(is,*root);
  }

}


/*-
 *-----------------------------------------------------------------------
 * initializeIterableDT 
 * 
 * Preconditions:
 *      Decision tree should not have been previously initialized. 
 *
 * Postconditions:
 *      Object is initialized to read in per-utterance decision trees 
 *
 * Side Effects:
 *      None
 *
 * Results:
 *      None
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::initializeIterableDT(
  string& fileName
  ) 
{
  //////////////////////////////////////////////////////////////////////
  // Make sure this instance already hasn't been initialized as an
  // iterable decision tree 
  //////////////////////////////////////////////////////////////////////
  if (iterable()) {
    error("ERROR: in DT named '%s', can't have DTs defined recursively in files",
      fileName.c_str()); 
  }

  //////////////////////////////////////////////////////////////////////
  // Initialize variables 
  //////////////////////////////////////////////////////////////////////
  dtFileName = fileName;
  dtFile     = new iDataStreamFile(dtFileName.c_str(), false, false);
  numDTs     = 0;
  dtNum      = -1;

  //////////////////////////////////////////////////////////////////////
  // Read in the first decision tree 
  //////////////////////////////////////////////////////////////////////
  beginIterableDT();
}


/*-
 *-----------------------------------------------------------------------
 * seek 
 *      Seek to a particular decision tree 
 *
 * Preconditions:
 *      Object should be filled in using 'read'.  
 *
 * Postconditions:
 *      The decision tree file pointer will point to the appropriate tree  
 *
 * Side Effects:
 *      If the index file has not been opened, it will be. 
 *
 * Results:
 *      none
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::seek(
  unsigned dt_nmbr 
  )
{
  string   indexFileName;
  unsigned numIndexDTs;
  unsigned position; 
  unsigned result; 

  if (!iterable()) {
    error("ERROR: trying to seek in non-iterable DT, '%s'\n",  
      curName.c_str() ); 
  }

  //////////////////////////////////////////////////////////////////////////
  // Open index file if necessary
  //////////////////////////////////////////////////////////////////////////
  if (indexFile == NULL) {
    indexFileName = dtFileName + DTFileExtension;
    indexFile = new iDataStreamFile(indexFileName.c_str(), true, false);
  }

  //////////////////////////////////////////////////////////////////////////
  // Check that the number of decision trees matches in the index file 
  //////////////////////////////////////////////////////////////////////////
  indexFile->rewind();
  indexFile->read(numIndexDTs, "num index DTs");

  if (numIndexDTs != numDTs) {
      error("ERROR: Index file '%s' lists '%d' decision trees does match '%s' which lists '%d'\n",
        indexFileName.c_str(), numIndexDTs, dtFileName.c_str(), numDTs );
  }
 
  //////////////////////////////////////////////////////////////////////////
  // Read in the position of the tree in the DT file 
  //////////////////////////////////////////////////////////////////////////
  indexFile->fseek( sizeof(unsigned)*(dt_nmbr), SEEK_CUR ); 
  indexFile->read(position, "DT offset");

  result = dtFile->fseek(position, SEEK_SET); 
  dtNum = dt_nmbr - 1; 
}

/*-
 *-----------------------------------------------------------------------
 * readRecurse:
 *      continues (from read) reading in from a file.
 * 
 * Preconditions:
 *      object should not be filled on (data will be lost if so).
 *
 * Postconditions:
 *      object will be filled in.
 *
 * Side Effects:
 *      Totally changes the object, might kill program if error occurs.
 *
 * Results:
 *      the current node it read. This routine also builds up the
 *      left (or prev) pointers in a doubly linked list of leaf
 *      nodes. The currently right-most leaf node is returned
 *      in the argument prevLeaf.
 *
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::readRecurse(iDataStreamFile& is,Node& node)
{
  int curFeat;
  is.read(curFeat,"Can't read DecisionTree's feature value");
  if (curFeat < -1) 
    error("ERROR: DT '%s', file '%s' line %d, feature number (=%d) must be non-negative",
	  name().c_str(),is.fileName(),is.lineNo(),curFeat);
  if (curFeat >= (int)_numFeatures) 
    error("ERROR: DT '%s', file '%s' line %d, feature number (=%d) must be < numFeatures (=%d)",
	  name().c_str(),is.fileName(),is.lineNo(),curFeat,
	  _numFeatures);

  if (curFeat == -1) {
    //////////////////////////////////////////////////////////////////////
    // Definitely is a leaf node, either value or equation.
    //////////////////////////////////////////////////////////////////////

    string leafNodeVal;
    int    val;

    //////////////////////////////////////////////////////////////////////
    // Skip leading spaces, then read until a '{' or a space is found 
    //////////////////////////////////////////////////////////////////////
    is.readStringUntil(leafNodeVal, '{', true,
        "Can't read RngDecisionTree's '{' character");

    //////////////////////////////////////////////////////////////////////
    // Check if leaf node is a formula
    //////////////////////////////////////////////////////////////////////
    if (leafNodeVal.size() == 0) {

      //////////////////////////////////////////////////////////////////////
      // Read until '}' is found (will not stop at spaces) 
      //////////////////////////////////////////////////////////////////////
      is.readStringUntil(leafNodeVal, '}', false, 
        "Can't read RngDecisionTree's '}' character");

      node.nodeType = LeafNodeEquation;
      new (&node.ln_e()) LeafNodeEquationStruct();
      
      try {
	node.ln_e().equation.parseFormula(leafNodeVal);
      }
      catch( string error_message ){
        error("ERROR: In file '%s' line %d, DT '%s', equation '%s':  %s", 
          is.fileName(), is.lineNo(),name().c_str(), leafNodeVal.c_str(), 
	  error_message.c_str());
      }
      catch( const char * const error_message ) {
        error("ERROR: In file '%s' line %d, DT '%s', equation '%s':  %s", 
          is.fileName(), is.lineNo(),name().c_str(), leafNodeVal.c_str(), 
	  error_message );
      }
    }
    //////////////////////////////////////////////////////////////////////
    // Check for integer
    //////////////////////////////////////////////////////////////////////
    else if (strIsInt(leafNodeVal.c_str(),&val)) {
      node.nodeType = LeafNodeVal;
      new (&node.ln_v()) LeafNodeValStruct();
      node.ln_v().value = val;
    }
    else {
      error("ERROR: In file '%s' line %d, DT '%s', invalid leaf node value '%s'.",
        is.fileName(), is.lineNo(),name().c_str(), leafNodeVal.c_str() ); 
    }
  }
  else {
    //////////////////////////////////////////////////////////////////////
    // Definitely a non leaf node, either array, hash, or rngs
    //////////////////////////////////////////////////////////////////////

    unsigned numSplits;
    is.read(numSplits,"Can't read DecisionTree's number of splits");
    if (numSplits < 1)
      error("ERROR: DT '%s', file '%s' line %d, can't have < 1 node splits",
	    name().c_str(),is.fileName(),is.lineNo());

    // we need to start reading the splits before we can determine what type
    // of node this is. Here are the cases.
    //
    // Array: if the splits are a list of consecutive integers of the form:
    //        3 4 5 6 7 8
    // or if it is a range of the form 
    //        3 ... 100
    // either of the appropriate length, then this can be handled by an array
    // of nodes, so we use the array form.

    // Hash: if the splits are a list of non-consequtive integers, then
    // a hash table is used.

    // Lastly, if it is a list of BP_Rngs, then the Rng array (the slowest
    // but most general of the three) is used.

    // in all cases, we end with the string 'default'

    // mode/state variable. 
    // 0 for array sequence of integers
    // 1 for array range [ start ... end ]
    // 2 for hash, and
    // 3 for ranges
    unsigned mode = 0;
    int firstUnsigned=~0x0;
    int previousUnsigned=~0x0;

    // keep a temporary copy of the split int values in case this is
    // mode 0, 1, or 2 since if we move to 3, we can't allocate 3's
    // data structures until we've deleted 2's.
    sArray < unsigned > splitIntVals(numSplits-1);

    // read (numSplits-1) ints or ranges.
    string token;
    for (unsigned i=0; (i<(numSplits-1)) || (mode == 1) ;i++) {

      is.read(token,"Can't read DecisionTree's integer value/range");

      if (mode == 0) {
	// we are still in array mode, search for an int.

	int tokenLen;
	int tokenVal;
	if (strIsInt(token.c_str(),&tokenVal,&tokenLen) && (tokenLen == (int)token.size())) {
	  if (tokenVal < 0) 
	    error("ERROR: DT '%s', file '%s' line %d: got negative int (%d) but must be non-negative",
		  name().c_str(),is.fileName(),is.lineNo(),tokenVal);
	  // then we got a valid unsigned.
	  if (i == 0) {
	    firstUnsigned = tokenVal;
	  } else {
	    // we've seen an int before.
	    if (tokenVal != previousUnsigned + 1) {
	      // then we've got a non-consecutive int, we need to
	      // move into hash mode (and hope we stay there).

	      // nothing to free because we haven't yet
	      // set up the array case.

	      // set up hash table and change mode.
	      mode = 2;
	      node.nodeType = NonLeafNodeHash;
	      new (& node.nln_h() ) NonLeafNodeHashStruct(numSplits-1);
	      node.nln_h().ftr = (leafNodeValType)curFeat;
	      node.nln_h().children.resize(numSplits);
	      // new (& (node.nln_h().children) ) sArray < Node >(numSplits);
	      // new (& (node.nln_h().nodeMapper) ) shash_map < unsigned, Node * > ();

	      // now insert all the integers we've got so far, all of which
	      // are consecutive except for this last one.
	      for (unsigned j=0;j<i;j++) {
		// since these are consecutive, we need not check the hash found condition.
		node.nln_h().nodeMapper.insert(firstUnsigned+j,&(node.nln_h().children[j]));
	      }
	      // insert the newly found non-consecutive int.
	      bool foundp;
	      node.nln_h().nodeMapper.insert(tokenVal,&(node.nln_h().children[i]),foundp);
	      if (foundp) 
		error("ERROR: DT '%s', file '%s' line %d: dupilcate integer entry %d",
		      name().c_str(),is.fileName(),is.lineNo(),tokenVal);
	    }
	  }
	  previousUnsigned = tokenVal;
	  splitIntVals[i] = tokenVal;
	} else if (token == "...") {
	  if (i != 1)
	    error("ERROR: DT '%s', file '%s' line %d: array range indicateor '...' can only occur at 2nd position between two ints",
		  name().c_str(),is.fileName(),is.lineNo());
	  // then this is an int array in "int ... int" form. 
	  mode = 1;
	} else {
	  // then this is perhaps a range, and we've moving from mode
	  // 0 to mode 3. Need to convert the (consecutive) ints we've
	  // seen so far into ranges, and move into most general but
	  // slowest range mode.

	  // Set up range array and change mode.
	  mode = 3;
	  node.nodeType = NonLeafNodeRngs;
	  new (& node.nln_r() ) NonLeafNodeRngsStruct();
	  node.nln_r().ftr = (leafNodeValType)curFeat;
	  // assume it is ordered, but we'll need to make sure
	  // later on.
	  node.nln_r().ordered = true;
	  // new (& (node.nln_r().children) ) sArray < RangeNode >(numSplits-1);
	  node.nln_r().children.resize(numSplits-1);
	  // allocate one more for the default.
	  node.nln_r().def = new Node;

	  // now insert all the integers we've got so far, all of which
	  // are consecutive integers except for this last one which should be a range.
	  for (unsigned j=0;j<i;j++) {
	    char buff[1024];
	    sprintf(buff,"%d",firstUnsigned+j);
	    // this will die if buff is not a rnage
	    new (&node.nln_r().children[j].rng) BP_Range(buff,
							 0,
							 MAX_BP_RANGE_VALUE);
	  }
	  // insert the newly found presumably range spec.
	  new (&node.nln_r().children[i].rng) BP_Range(token.c_str(),
						       0,
						       MAX_BP_RANGE_VALUE);

	}
      } else if (mode == 1) {
	// we are in array range mode. We must have already seen an an
	// "int ..." and we expect an "int default" for this to be
	// valid. We finish up here.

	int tokenLen;
	int tokenVal;
	if (!strIsInt(token.c_str(),&tokenVal,&tokenLen) || (tokenLen != (int)token.size()) || (tokenVal < 0)) {
	    error("ERROR: DT '%s', file '%s' line %d: expecting non negative int at split position %d, but got '%s'",
		  name().c_str(),is.fileName(),is.lineNo(),i,token.c_str());
	}

	// tokenVal - firstUnsigned + 1 + 1 == numSplits
	if ((unsigned)tokenVal + 2 != numSplits + firstUnsigned) {
	    error("ERROR: DT '%s', file '%s' line %d: array range '%d ... %d' must be %d long",
		  name().c_str(),is.fileName(),is.lineNo(),firstUnsigned,tokenVal,numSplits-1);
	}

	// everything checks out ok. Allocate.
	node.nodeType = NonLeafNodeArray;
	new (& node.nln_a() ) NonLeafNodeArrayStruct();
	node.nln_a().ftr = (leafNodeValType)curFeat;
	node.nln_a().base = (leafNodeValType)firstUnsigned;
	node.nln_a().children.resize(numSplits);

	// leave the loop. 
	goto doneWithSplits;

      } else if (mode == 2) {
	// we are in hash mode, having not survived array mode in
	// at least first iteration.

	assert ( i > 0 );
	int tokenLen;
	int tokenVal;
	if (strIsInt(token.c_str(),&tokenVal,&tokenLen) && (tokenLen == (int)token.size())) {
	  // then we got a valid int.
	  if (tokenVal < 0) 
	    error("ERROR: DT '%s', file '%s' line %d: got negative int (%d) but must be non-negative",
		  name().c_str(),is.fileName(),is.lineNo(),tokenVal);

	  // insert the newly found int token
	  bool foundp;
	  node.nln_h().nodeMapper.insert(tokenVal,&(node.nln_h().children[i]),foundp);
	  if (foundp) 
	    error("ERROR: DT '%s', file '%s' line %d: dupilcate integer entry %d",
		  name().c_str(),is.fileName(),is.lineNo(),tokenVal);
	  splitIntVals[i] = tokenVal;
	} else {
	  // need to move to range mode (if we've got a valid range)

	  // then this is perhaps a range. Need to convert
	  // the ints we've seen so far into ranges, and move into
	  // most general but slowest range mode. Need to delete
	  // all existing hash-mode structures as well.  

	  // delete hash mode structures since we came from  hash mode.
	  // node.nln_h().children.~sArray < Node >();
	  // node.nln_h().nodeMapper.~shash_map < unsigned, Node * > ();
	  node.nln_h().~NonLeafNodeHashStruct();

	  // set up range array and change mode.
	  mode = 3;
	  node.nodeType = NonLeafNodeRngs;
	  new (& node.nln_r() ) NonLeafNodeRngsStruct();
	  node.nln_r().ftr = (leafNodeValType)curFeat;
	  // assume it is ordered, but we'll need to make sure
	  // later on.
	  node.nln_r().ordered = true;
	  node.nln_r().children.resize(numSplits-1);

	  // allocate one more for the default.
	  node.nln_r().def = new Node;

	  // now insert all the integers we've got so far, all of which
	  // are consecutive integers except for this last one which should be a range.
	  for (unsigned j=0;j<i;j++) {
	    char buff[1024];
	    sprintf(buff,"%d",splitIntVals[j]);
	    // this will die if buff is not a range
	    new (&node.nln_r().children[j].rng) BP_Range(buff,
							 0,
							 MAX_BP_RANGE_VALUE);
	  }
	  // insert the newly found presumably range spec.
	  new (&node.nln_r().children[i].rng) BP_Range(token.c_str(),
						       0,
						       MAX_BP_RANGE_VALUE);
	}
      } else {
	// we are in rngs mode, having not survived either array or hash mode.
	assert ( mode == 3 );

	// insert the newly found presumably range spec.
	new (&node.nln_r().children[i].rng) BP_Range(token.c_str(),
						     0,
						     MAX_BP_RANGE_VALUE);

      }
    }
  doneWithSplits:
    if (mode == 2) {
      // TODO: we're in hash mode, but if the user
      // has consecutive values, we should move to mode == 0
      // to speed up querying. We sort the stored values to
      // check if they are consecutive. If they are,
      // we use array mode. If not, we use hash mode.
    }

    if (mode == 0) {
      // still haven't done this.
      node.nodeType = NonLeafNodeArray;
      new (& node.nln_a()) NonLeafNodeArrayStruct();
      node.nln_a().ftr = (leafNodeValType)curFeat;
      node.nln_a().base = (leafNodeValType)firstUnsigned;
      node.nln_a().children.resize(numSplits);
    }

    
    // ranges are ok, now make sure we get default string.
    is.read(token,"Can't read DecisionTree's 'default'");
    if (strcmp(RNG_DECISION_TREE_DEF_STR,token.c_str()))
      error("ERROR: DT '%s', file '%s' line %d: expecting default str (%s) got (%s)",
	    name().c_str(),is.fileName(),is.lineNo(),RNG_DECISION_TREE_DEF_STR,token.c_str());

    // free up the memory since we don't need it any longer.
    splitIntVals.clear();

    if (node.nodeType == NonLeafNodeRngs) {

      //////////////////////////////////////////////
      // check for overlap and order errors in the strings
      // TODO: this is way inefficient, see if it is possible not to do
      // an N^2 algorithm. How? Do own version quicksort that backs out when
      // it finds a non-comparable entry.

      for (unsigned i=0;i<numSplits-1;i++) {
	for (unsigned j=i+1;j<numSplits-1;j++) {
	  if (node.nln_r().children[i].rng.overlapP(node.nln_r().children[j].rng))
	    error("ERROR: DT '%s', file '%s' line %d: range %d (%s) and %d (%s) have a non-empty intersection.",
		  name().c_str(),is.fileName(),is.lineNo(),
		  i,node.nln_r().children[i].rng.rangeStr(),
		  j,node.nln_r().children[j].rng.rangeStr()
		  );

#if 0	
	  ////////////////////////////////////////////////////////
	  // this next check is required to ensure there is
	  // an ordering of the ranges, so sorting makes sense
	  // when we do binary search. An unordered list of ranges would, for
	  // example, have things like 1,3,5 and 2,4 where you can't say
	  // that either (1,3,5) < (2,4) nor (2,4) < (1,3,5). In these
	  // cases, we do not sort the lists, and have to resort to linear search
	  // when doing a query.
	  if (node.nln_r().ordered
	      &&
	      !(node.nln_r().children[i].rng < node.nln_r().children[j].rng)
	      && 
	      !(node.nln_r().children[j].rng < node.nln_r().children[i].rng)
	      ) {
	    // then it isn't ordered, and we'll have to do linear search
	    // on a query. 
	    // 
	    // TODO: check if it is ordered alread, and if
	    // so, don't do the sort below.
	    node.nln_r().ordered = false;
	  }
#endif
	}
      }
    }

    if (node.nodeType == NonLeafNodeArray) {
      for (unsigned i=0; i<numSplits; i++) {
	readRecurse(is,node.nln_a().children[i]);
      }
    } else if (node.nodeType == NonLeafNodeHash) {
      for (unsigned i=0; i<numSplits; i++) {
	readRecurse(is,node.nln_h().children[i]);
      }
    } else if (node.nodeType == NonLeafNodeRngs) {
      for (unsigned i=0; i<numSplits - 1; i++) {
	readRecurse(is,node.nln_r().children[i].nd);
      }
      // default case is special here.
      readRecurse(is,*(node.nln_r().def));
    } else {
      // shouldn't happen.
      assert ( 0 );
    }


    if (node.nodeType == NonLeafNodeRngs) { 
      // when ordered, do the sort right now.
      if (numSplits >= DT_SPLIT_SORT_THRESHOLD && node.nln_r().ordered) 
	{
	  //////////////////////////////////////////////////////////
	  // sort the entries in ascending so we can do bin-like search later. Sort only
	  // if there are a sufficient number to warrant a sort (condition above).
	  // note that even if some of the entries aren't comparable (i.e., are "equal",
	  // the binary search should check this condition.
	  // printf("---------- Sorting\n");
	  node.nln_r().children.sort();
	  // printf("-------- Done Sorting\n");
	} else {
	  // didn't sort so can't be ordered.
	  node.nln_r().ordered = false;
	}
    }

  }
}


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::EquationClass
 *   Constructor 
 *
 * Preconditions:
 *   none 
 *
 * Postconditions:
 *   Initializes static members if not already done 
 *
 * Side Effects:
 *   none   
 *
 * Results:
 *   none   
 *-----------------------------------------------------------------------
 */
RngDecisionTree::EquationClass::EquationClass()
{
  //////////////////////////////////////////////////////////////////////////
  // Initialize the static maps of tokens 
  //////////////////////////////////////////////////////////////////////////
  if (delimiter.size() == 0) {

    ////////////////////////////////////////////////////////////////////////
    // 'delimiter' and 'function' map text to the tokenEnum 
    ////////////////////////////////////////////////////////////////////////
    delimiter["&"]  = TOKEN_BITWISE_AND;
    delimiter["|"]  = TOKEN_BITWISE_OR;
    delimiter["~"]  = TOKEN_BITWISE_NOT;
    delimiter[":"]  = TOKEN_COLON;
    delimiter[","]  = TOKEN_COMMA;
    delimiter["./"] = TOKEN_DIVIDE_CEIL;
    delimiter["/"]  = TOKEN_DIVIDE_FLOOR;
    delimiter["/."] = TOKEN_DIVIDE_FLOOR;
    delimiter["~/"] = TOKEN_DIVIDE_ROUND;
    delimiter["=="] = TOKEN_EQUALS;
    delimiter["^"]  = TOKEN_EXPONENT;
    delimiter[">"]  = TOKEN_GREATER_THAN;
    delimiter[">="] = TOKEN_GREATER_THAN_EQ;
    delimiter["<"]  = TOKEN_LESS_THAN;
    delimiter["<="] = TOKEN_LESS_THAN_EQ;
    delimiter["&&"] = TOKEN_LOGICAL_AND;
    delimiter["||"] = TOKEN_LOGICAL_OR;
    delimiter["-"]  = TOKEN_MINUS;
    delimiter["!"]  = TOKEN_NOT;
    delimiter["!="] = TOKEN_NOT_EQUAL;
    delimiter["+"]  = TOKEN_PLUS;
    delimiter["<<"] = TOKEN_SHIFT_LEFT;
    delimiter[">>"] = TOKEN_SHIFT_RIGHT;
    delimiter[" "]  = TOKEN_SPACE;
    delimiter["\t"] = TOKEN_SPACE;
    delimiter["("]  = TOKEN_LEFT_PAREN;
    delimiter["*"]  = TOKEN_TIMES;
    delimiter[")"]  = TOKEN_RIGHT_PAREN;
    delimiter["?"]  = TOKEN_QUESTION_MARK;

    function["abs"] = TOKEN_ABSOLUTE_VALUE;
    function["ceil_divide"]  = TOKEN_DIVIDE_CEIL_FUNCTION;
    function["floor_divide"] = TOKEN_DIVIDE_FLOOR_FUNCTION;
    function["round_divide"] = TOKEN_DIVIDE_ROUND_FUNCTION;
    function["max"]     = TOKEN_MAX;
    function["median"]  = TOKEN_MEDIAN;
    function["min"]     = TOKEN_MIN;
    function["mod"]     = TOKEN_MOD;
    function["all_diff"] = TOKEN_ALL_DIFFERENT;
    function["rotate"]  = TOKEN_ROTATE;
    function["xor"]     = TOKEN_BITWISE_XOR;

    variable["cardinality_child"]   = TOKEN_CARDINALITY_CHILD;
    variable["cc"]                  = TOKEN_CARDINALITY_CHILD;
    variable["cardinality_parent_"] = TOKEN_CARDINALITY_PARENT;
    variable["cp"]                  = TOKEN_CARDINALITY_PARENT;
    variable["parent_"]             = TOKEN_PARENT_VALUE;
    variable["p"]                   = TOKEN_PARENT_VALUE;
    variable["parent_minus_one_"]   = TOKEN_PARENT_VALUE_MINUS_ONE;
    variable["pmo"]                 = TOKEN_PARENT_VALUE_MINUS_ONE;
    variable["parent_plus_one_"]    = TOKEN_PARENT_VALUE_PLUS_ONE;
    variable["ppo"]                 = TOKEN_PARENT_VALUE_PLUS_ONE;
    variable["max_value_child"]     = TOKEN_MAX_VALUE_CHILD;
    variable["mc"]                  = TOKEN_MAX_VALUE_CHILD;
    variable["max_value_parent_"]   = TOKEN_MAX_VALUE_PARENT;
    variable["mp"]                  = TOKEN_MAX_VALUE_PARENT;

    ////////////////////////////////////////////////////////////////////////
    //  Maps tokens to commands
    ////////////////////////////////////////////////////////////////////////
    infixToken[TOKEN_BITWISE_AND]     = COMMAND_BITWISE_AND;
    infixToken[TOKEN_BITWISE_OR]      = COMMAND_BITWISE_OR;
    infixToken[TOKEN_BITWISE_XOR]     = COMMAND_BITWISE_XOR;
    infixToken[TOKEN_DIVIDE_CEIL]     = COMMAND_DIVIDE_CEIL;
    infixToken[TOKEN_DIVIDE_FLOOR]    = COMMAND_DIVIDE_FLOOR;
    infixToken[TOKEN_DIVIDE_ROUND]    = COMMAND_DIVIDE_ROUND;
    infixToken[TOKEN_EXPONENT]        = COMMAND_EXPONENT;
    infixToken[TOKEN_EQUALS]          = COMMAND_EQUALS;
    infixToken[TOKEN_NOT_EQUAL]       = COMMAND_NOT_EQUAL;
    infixToken[TOKEN_GREATER_THAN]    = COMMAND_GREATER_THAN;
    infixToken[TOKEN_GREATER_THAN_EQ] = COMMAND_GREATER_THAN_EQ;
    infixToken[TOKEN_LESS_THAN]       = COMMAND_LESS_THAN;
    infixToken[TOKEN_LESS_THAN_EQ]    = COMMAND_LESS_THAN_EQ;
    infixToken[TOKEN_LOGICAL_AND]     = COMMAND_LOGICAL_AND;
    infixToken[TOKEN_LOGICAL_OR]      = COMMAND_LOGICAL_OR; 
    infixToken[TOKEN_MINUS]           = COMMAND_MINUS;
    infixToken[TOKEN_PLUS]            = COMMAND_PLUS;
    infixToken[TOKEN_QUESTION_MARK]   = COMMAND_BRANCH_IF_FALSE;
    infixToken[TOKEN_SHIFT_LEFT]      = COMMAND_SHIFT_LEFT;
    infixToken[TOKEN_SHIFT_RIGHT]     = COMMAND_SHIFT_RIGHT;
    infixToken[TOKEN_TIMES]           = COMMAND_TIMES;

    unaryToken[TOKEN_BITWISE_NOT]     = COMMAND_BITWISE_NOT;
    unaryToken[TOKEN_MINUS]           = COMMAND_NEGATE;
    unaryToken[TOKEN_NOT]             = COMMAND_NOT;

    functionToken[TOKEN_MAX] = COMMAND_MAX;
    functionToken[TOKEN_MIN] = COMMAND_MIN;

    oneValFunctionToken[TOKEN_ABSOLUTE_VALUE] = COMMAND_ABSOLUTE_VALUE;

    twoValFunctionToken[TOKEN_BITWISE_XOR]  = COMMAND_BITWISE_XOR;
    twoValFunctionToken[TOKEN_DIVIDE_CEIL_FUNCTION]  = COMMAND_DIVIDE_CEIL;
    twoValFunctionToken[TOKEN_DIVIDE_FLOOR_FUNCTION] = COMMAND_DIVIDE_FLOOR;
    twoValFunctionToken[TOKEN_DIVIDE_ROUND_FUNCTION] = COMMAND_DIVIDE_ROUND;
    twoValFunctionToken[TOKEN_MOD]          = COMMAND_MOD;

    manyValFunctionToken[TOKEN_ALL_DIFFERENT] = COMMAND_ALL_DIFFERENT;

    // Special cases
    // TOKEN_MEDIAN => COMMAND_MEDIAN
    // TOKEN_ROTATE => COMMAND_ROTATE

    variableToken[TOKEN_CARDINALITY_CHILD]  = COMMAND_PUSH_CARDINALITY_CHILD;
    variableToken[TOKEN_CARDINALITY_PARENT] = COMMAND_PUSH_CARDINALITY_PARENT;
    variableToken[TOKEN_PARENT_VALUE]       = COMMAND_PUSH_PARENT_VALUE;
    variableToken[TOKEN_PARENT_VALUE_MINUS_ONE] = 
      COMMAND_PUSH_PARENT_VALUE_MINUS_ONE;
    variableToken[TOKEN_PARENT_VALUE_PLUS_ONE]  = 
      COMMAND_PUSH_PARENT_VALUE_PLUS_ONE;
    variableToken[TOKEN_MAX_VALUE_CHILD]    = COMMAND_PUSH_MAX_VALUE_CHILD;
    variableToken[TOKEN_MAX_VALUE_PARENT]   = COMMAND_PUSH_MAX_VALUE_PARENT;

    ////////////////////////////////////////////////////////////////////////
    //  Maps tokens to priority levels 
    ////////////////////////////////////////////////////////////////////////
    tokenPriority[TOKEN_BITWISE_AND]     = BITWISE_AND_PRCDNC;
    tokenPriority[TOKEN_BITWISE_NOT]     = UNARY_PRCDNC;
    tokenPriority[TOKEN_BITWISE_OR]      = BITWISE_OR_PRCDNC;
    tokenPriority[TOKEN_DIVIDE_CEIL]     = MULT_PRCDNC;
    tokenPriority[TOKEN_DIVIDE_FLOOR]    = MULT_PRCDNC;
    tokenPriority[TOKEN_DIVIDE_ROUND]    = MULT_PRCDNC;
    tokenPriority[TOKEN_EXPONENT]        = EXPONENT_PRCDNC;
    tokenPriority[TOKEN_EQUALS]          = EQUALITY_PRCDNC;
    tokenPriority[TOKEN_GREATER_THAN]    = RELATIONAL_PRCDNC;
    tokenPriority[TOKEN_GREATER_THAN_EQ] = RELATIONAL_PRCDNC; 
    tokenPriority[TOKEN_LEFT_PAREN]      = PAREN_PRCDNC;
    tokenPriority[TOKEN_LESS_THAN]       = RELATIONAL_PRCDNC;
    tokenPriority[TOKEN_LESS_THAN_EQ]    = RELATIONAL_PRCDNC;
    tokenPriority[TOKEN_LOGICAL_AND]     = LOGICAL_AND_PRCDNC;
    tokenPriority[TOKEN_LOGICAL_OR]      = LOGICAL_OR_PRCDNC;
    tokenPriority[TOKEN_MINUS]           = ADDITIVE_PRCDNC;
    tokenPriority[TOKEN_NOT]             = UNARY_PRCDNC;
    tokenPriority[TOKEN_NOT_EQUAL]       = EQUALITY_PRCDNC;
    tokenPriority[TOKEN_PLUS]            = ADDITIVE_PRCDNC;
    tokenPriority[TOKEN_QUESTION_MARK]   = CONDITIONAL_PRCDNC;
    tokenPriority[TOKEN_RIGHT_PAREN]     = PAREN_PRCDNC;
    tokenPriority[TOKEN_SHIFT_LEFT]      = SHIFT_PRCDNC;
    tokenPriority[TOKEN_SHIFT_RIGHT]     = SHIFT_PRCDNC;
    tokenPriority[TOKEN_TIMES]           = MULT_PRCDNC;
  }

}


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::evaluateFormula
 *   Calculate an equation's value given a set parent variables 
 *
 * Preconditions:
 *   none 
 *
 * Postconditions:
 *   none 
 *
 * Side Effects:
 *   none   
 *
 * Results:
 *   Returns a single value 
 *-----------------------------------------------------------------------
 */
leafNodeValType 
RngDecisionTree::EquationClass::evaluateFormula(
	const vector< RV* >& variables,
	const RV* const rv
)
{
  unsigned crrnt_cmnd, end_cmnd;
  unsigned last; 
  unsigned command, operand; 
  int      i;
  stack_element_t value;
  unsigned val, number, position, bitwidth;
  unsigned mask_1, mask_2;

  stack.clear();

  for (crrnt_cmnd = 0, 
       end_cmnd   = commands.size();
       crrnt_cmnd != end_cmnd;
       ++crrnt_cmnd) {

    command = GET_COMMAND(commands[crrnt_cmnd]); 

    char* const missingParentErrorString = "ERROR: Reference to non-existant parent variable in formula. Asking for parent %d but only %d parents are available.\n";

    switch (command) {
	
      case COMMAND_PUSH_CARDINALITY_CHILD:
        stack.push_back( (rv->discrete() ? RV2DRV(rv)->cardinality : 0) );
        break;

      case COMMAND_PUSH_CARDINALITY_PARENT:	
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        if (operand >= variables.size()) {	
          error(missingParentErrorString,operand,variables.size());
        }
        stack.push_back( (variables[operand]->discrete() ? 
          RV2DRV(variables[operand])->cardinality : 0) );
        break;	

      case COMMAND_PUSH_PARENT_VALUE:	
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        if (operand >= variables.size()) {	
          error(missingParentErrorString,operand,variables.size());
        }
        stack.push_back( RV2DRV(variables[operand])->discrete() ? 
          RV2DRV(variables[operand])->val : 0 );	
        break;	

      case COMMAND_PUSH_PARENT_VALUE_MINUS_ONE:	
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        if (operand >= variables.size()) {	
          error(missingParentErrorString,operand,variables.size());
        }
        stack.push_back( RV2DRV(variables[operand])->discrete() ? 
          (RV2DRV(variables[operand])->val-1) : 0 );	
        break;	

      case COMMAND_PUSH_PARENT_VALUE_PLUS_ONE:	
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        if (operand >= variables.size()) {	
          error(missingParentErrorString,operand,variables.size());
        }
        stack.push_back( RV2DRV(variables[operand])->discrete() ? 
          (RV2DRV(variables[operand])->val+1) : 0 );	
        break;	

      case COMMAND_PUSH_MAX_VALUE_CHILD:	
        stack.push_back( (rv->discrete() ? (RV2DRV(rv)->cardinality - 1) : 0) );
        break;	

      case COMMAND_PUSH_MAX_VALUE_PARENT:	
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        if (operand >= variables.size()) {	
          error(missingParentErrorString,operand,variables.size());
        }
        stack.push_back( (variables[operand]->discrete() ? 
          (RV2DRV(variables[operand])->cardinality - 1) : 0) );
        break;	

      case COMMAND_PUSH_CONSTANT:
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        stack.push_back(operand); 
        break;

      case COMMAND_ABSOLUTE_VALUE: 
        last = stack.stackSize() - 1;
        if (stack[last] < 0) { 
          stack[last] = -stack[last];
        }
        break;

      case COMMAND_ALL_DIFFERENT:
        int  index_1;
        int  index_2;
        bool all_different; 
    
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        all_different = true;

        for( index_1=(stack.stackSize()-operand); 
             (index_1<stack.stackSize()) && (all_different); 
             index_1++)
        {
          for(index_2=index_1+1; index_2<stack.stackSize(); index_2++)
          {
            if (stack[index_1] == stack[index_2]) 
            {
              all_different = false;
              break;
            }
          }
        }

        stack[stack.stackSize()-operand] = all_different;
        stack.pop_back(operand-1);
        break;
 
      case COMMAND_BITWISE_AND: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] & stack[last];
        stack.pop_back();
        break;
 
      case COMMAND_BITWISE_OR: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] | stack[last];
        stack.pop_back();
        break;

      case COMMAND_BITWISE_NOT:
        last = stack.stackSize() - 1;
        stack[last] = ~(stack[last]); 
        break;
 
      case COMMAND_BITWISE_XOR: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] ^ stack[last];
        stack.pop_back();
        break;

      case COMMAND_BRANCH: 
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        crrnt_cmnd += operand;
        break;

      case COMMAND_BRANCH_IF_FALSE: 
        if (!stack[stack.stackSize()-1]) {
          operand = GET_OPERAND(commands[crrnt_cmnd]); 
          crrnt_cmnd += operand;
        }
        stack.pop_back();
        break;
 
      case COMMAND_DIVIDE_CEIL:
        last = stack.stackSize() - 1;
        if (stack[last] == 0) {
          error("ERROR:  Divide by zero error\n"); 
        }
        stack[last-1] = (stack[last-1] + stack[last] - 1) / stack[last];
        stack.pop_back();
        break;
 
      case COMMAND_DIVIDE_FLOOR:
        last = stack.stackSize() - 1;
        if (stack[last] == 0) {
          error("ERROR:  Divide by zero error\n"); 
        }
        stack[last-1] = stack[last-1] / stack[last];
        stack.pop_back();
        break;
 
      case COMMAND_DIVIDE_ROUND:
        last = stack.stackSize() - 1;
        if (stack[last] == 0) {
          error("ERROR:  Divide by zero error\n"); 
        }
        stack[last-1] = (stack[last-1] + (stack[last]>>1)) / stack[last];
        stack.pop_back();
        break;

      case COMMAND_EQUALS: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] == stack[last];
        stack.pop_back();
        break;
 
      case COMMAND_EXPONENT: 
        last = stack.stackSize() - 1;
        value = 1; 
        for (i=0; i<stack[last]; ++i) {
          value = value*stack[last-1]; 
        }
        stack[last-1] = value;
        stack.pop_back();
        break;

      case COMMAND_GREATER_THAN: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] > stack[last];
        stack.pop_back();
        break;

      case COMMAND_GREATER_THAN_EQ: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] >= stack[last];
        stack.pop_back();
        break;

      case COMMAND_LESS_THAN: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] < stack[last];
        stack.pop_back();
        break;

      case COMMAND_LESS_THAN_EQ: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] <= stack[last];
        stack.pop_back();
        break;
 
      case COMMAND_LOGICAL_AND: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] && stack[last];
        stack.pop_back();
        break;
 
      case COMMAND_LOGICAL_OR: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] || stack[last];
        stack.pop_back();
        break;

      case COMMAND_MAX:
        last = stack.stackSize() - 1;
        if (stack[last] > stack[last-1]) {
          stack[last-1] = stack[last];
        }
        stack.pop_back();
        break;

      //////////////////////////////////////////////////////////////////////
      // Calculate the median
      //////////////////////////////////////////////////////////////////////
      case COMMAND_MEDIAN:

        stack_element_t pivot_value;
        stack_element_t high_value;
        unsigned pivot_index;
        unsigned index;
        unsigned first;
        unsigned lower_bound_index, upper_bound_index;
        unsigned pivot_lb, pivot_ub;  
        unsigned not_eq_pivot_index;  
        unsigned median_index;  
        bool     odd_nmbr_vls;
        bool     high_value_found;

        operand = GET_OPERAND(commands[crrnt_cmnd]); 

        first = stack.stackSize() - operand;
        last  = stack.stackSize() - 1;

        median_index = (first+last)/2;

        if (operand & 1) {
          odd_nmbr_vls = true;
        }
        else {
          odd_nmbr_vls = false;
        }

        // Entries with index <  lower_bound_index are <= pivot_value
        // Entries with index >= upper_bound_index are >= pivot_value
        lower_bound_index = first;
        upper_bound_index = last; 

        high_value = 0;
        high_value_found = false;
        do {

          //////////////////////////////////////////////////////////////////// 
          // Partial sort:  One of the elements is randomly chosen as a pivot,
          // the array is then shuffled so that all values <= pivot_value are
          // at the bottom of the array, and all values > pivot_value are
          // at the top of the array.
          ////////////////////////////////////////////////////////////////////  

          // Randomly choose a pivot 
          pivot_index = rnd.uniform(lower_bound_index, upper_bound_index);
          pivot_value = stack[pivot_index];

          // Put pivot at the lower_bound_index  
          swap(stack[lower_bound_index], stack[pivot_index]);

          // Entries with an index < pivot_lb have values <= pivot_value
          // Entries with an index > pivot_ub have values >  pivot_value
          // Entries with index s.t. (pivot_lb>=index>= pivot_ub) are unknown 
          // 
          // Move entries that have values>pivot_value to spots above the 
          // upper bound.
          pivot_lb = lower_bound_index+1;
          pivot_ub = upper_bound_index;

          while (pivot_lb<pivot_ub) { 
            if (stack[pivot_lb] <= pivot_value) {
              pivot_lb++;
            }
            else {          
              swap(stack[pivot_lb], stack[pivot_ub]);
              pivot_ub--; 
            }
          }

          // Now pivot_lb=pivot_ub, so look at the value at pivot_lb to
          // determine which side of the pivot it should be on.
          // pivot_ub might be an invalid index
          if (stack[pivot_lb] <= pivot_value) {
            pivot_ub++; 
          }
          else {
            pivot_lb--;
          }

          ////////////////////////////////////////////////////////////////////  
          // Deal with ties:  When there are values equal to the pivot this
          // can cause problems because the lower_bound_index and 
          // upper_bound_index will never converge.  This moves the 
          // pivot_value and all entries equal to it so they sit together 
          // at the pivot boundary 
          ////////////////////////////////////////////////////////////////////  

          // Put pivot itself (still stored at lower_bound_index) at pivot 
          // boundary 
          swap(stack[lower_bound_index],stack[pivot_lb]);

          // Now find and move other entries that equal pivot_value to the 
          // center.  If the pivot_lb==lower_bound_index we know that there
          // are not any ties (or any values <= to it).
          // The variable not_eq_pivot_index will store the largest index of 
          // a value strictly less than the pivot 
          if (pivot_lb <= lower_bound_index) {
            not_eq_pivot_index = pivot_lb;
          }
          else {
            not_eq_pivot_index = pivot_lb-1;   
            index = lower_bound_index;
            while ((not_eq_pivot_index > 0) && (index <= not_eq_pivot_index)) {
              if (stack[index] == pivot_value) {
                swap(stack[index], stack[not_eq_pivot_index]);
                if (not_eq_pivot_index > 0 ) {
                  not_eq_pivot_index--;
                }
              }
              else {
                index++;
              }
            }
          }

          ////////////////////////////////////////////////////////////////////  
          // Move upper and lower bounds and check stoppping conditions
          ////////////////////////////////////////////////////////////////////  

          //////////////////////////////////////////////////////////////////  
          // Note that:
          //   (pivot_lb+1)=pivot_ub
          //   Entries with an index <= pivot_lb have values <= pivot_value
          //   Entries s.t. (not_eq_pivot_index<=index<pivot_ub) = pivot_value
          //   Entries with an index >= pivot_ub have values >  pivot_value
          //////////////////////////////////////////////////////////////////  

          if (median_index <= not_eq_pivot_index) {
            upper_bound_index = not_eq_pivot_index; 
          }
          else if (median_index < pivot_ub) {
            upper_bound_index = lower_bound_index = median_index;
          }
          else { 
            lower_bound_index = pivot_ub;

            // Should never be the case that median_index > not_eq_pivot_index 
            // and pivot_ub is invalid, but protect this with an assert 
            assert(lower_bound_index<=last);
          }

          //////////////////////////////////////////////////////////////////  
          // When there are an even number of entries, do a second search
          // to find the entry with the median_index+1'th value 
          //////////////////////////////////////////////////////////////////  
          if ((lower_bound_index >= upper_bound_index) && 
              (high_value_found == false) &&
              (!odd_nmbr_vls)) {
            high_value_found = true;
            high_value = stack[median_index];
            median_index++;
            lower_bound_index = median_index;
            upper_bound_index = last;
          }

        } while (lower_bound_index < upper_bound_index);

        /////////////////////////////////////////////////////////////////////
        // With an odd number of values median is middle value
        // With an even number of values, median is average of the two values
        //   straddling the middle
        /////////////////////////////////////////////////////////////////////
        if (odd_nmbr_vls) {
          stack[first] = stack[median_index];
        }
        else {
          stack[first] = (high_value+stack[median_index]+1)>>1;
        }

        stack.pop_back(operand-1);

        break;
        ////////////////////////////////////////////////////////////////////  
        // End of MEDIAN 
        ////////////////////////////////////////////////////////////////////  

      case COMMAND_MIN:
        last = stack.stackSize() - 1;
        if (stack[last] < stack[last-1]) {
          stack[last-1] = stack[last];
        }
        stack.pop_back();
        break;
 
      case COMMAND_MINUS:
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] - stack[last];
        stack.pop_back();
        break;
 
      case COMMAND_MOD:
        last = stack.stackSize() - 1;
        if (stack[last] == 0) {
          error("ERROR:  Mod by zero error\n"); 
        }
        stack[last-1] = stack[last-1] % stack[last];
        stack.pop_back();
        break;

      case COMMAND_NEGATE:
        last = stack.stackSize() - 1;
        stack[last] = -(stack[last]); 
        break;

      case COMMAND_NOT:
        last = stack.stackSize() - 1;
        stack[last] = !(stack[last]); 
        break;

      case COMMAND_NOT_EQUAL: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] != stack[last];
        stack.pop_back();
        break;
 
      case COMMAND_PLUS:
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] + stack[last];
        stack.pop_back();
        break;

      case COMMAND_ROTATE:
        //////////////////////////////////////////////////////////////////////
        // rotate(val, number, position, bitwidth)
        //////////////////////////////////////////////////////////////////////
        last     = stack.stackSize() - 1;
        bitwidth = stack[last]; 
        position = stack[last-1];
        value    = (int)stack[last-2];

        if (value < 0) {
          number   = (-value) % bitwidth; 

          mask_1 = ((1<<(bitwidth-number))-1) << (position); 
          mask_2 = ((1<<number)-1) << (position+bitwidth-number); 

          val    = stack[last-3] & ~mask_1 & ~mask_2; 
          val    |= ((stack[last-3] & mask_1) << number); 
          val    |= ((stack[last-3] & mask_2) >> (bitwidth-number)); 
        }
        else {
          number   = value % bitwidth; 

          mask_1 = ((1<<number)-1) << (position); 
          mask_2 = ((1<<(bitwidth-number))-1) << (position+number); 

          val    = stack[last-3] & ~mask_1 & ~mask_2; 
          val    |= ((stack[last-3] & mask_1) << (bitwidth-number)); 
          val    |= ((stack[last-3] & mask_2) >> number); 
        }

        stack[last-3] = val;          
        stack.pop_back(3);
        break;

      case COMMAND_SHIFT_LEFT:
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] << stack[last];
        stack.pop_back();
        break;

      case COMMAND_SHIFT_RIGHT:
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] >> stack[last];
        stack.pop_back();
        break;

      case COMMAND_TIMES: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] * stack[last];
        stack.pop_back();
        break;

      default:
        assert(0);
        break;
    }
  }

  assert(stack.stackSize() == 1);
  return((unsigned)stack[0]);
}    


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::parseFormula
 *   Parse a string containing an equation.  The equation is converted 
 *   into a series of integer commands that are executed by evaluateFormula.
 * 
 * Preconditions:
 *   none      
 *
 * Postconditions:
 *   The decision tree's equation data structure is filled in.  The 
 *   computation stack will grow large enough to evaluate this equation. 
 *      
 * Side Effects:
 *   none      
 *
 * Results:
 *   none 
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::EquationClass::parseFormula(
  string formula
  )
{
  parsingCommandContainer new_commands;
  tokenStruct             token;
  unsigned depth = 0; 
  unsigned i; 

  assert( (LAST_COMMAND_INDEX & COMMAND_MASK) == LAST_COMMAND_INDEX );

  preProcessFormula(formula);

  getToken(formula, token);

  parseExpression(token, formula, new_commands, LOWEST_PRECEDENCE, depth);

  if (token.token != TOKEN_END) {
    string error_message = "Invalid token at " + formula;
    throw(error_message); 
  }

  assert(formula.length() == 0);

  if (depth != 1) {
    string error_message = "Not enough factors for given operators";
    throw(error_message); 
  }

  commands.resize( new_commands.size() ); 
  for (i=0; i<new_commands.size(); ++i) {
    commands[i] = new_commands[i]; 
  }
}


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::preProcessFormula
 * 
 * Preconditions:
 *   none      
 *
 * Postconditions:
 *   Takes a formula as given by the user and outputs one which the 
 *   parser can handle.  Specifically, it removes whitespace and changes
 *   all characters to lower case.
 *      
 * Side Effects:
 *   none      
 *
 * Results:
 *   none 
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::EquationClass::preProcessFormula(
  string& original 
  )
{
  unsigned i;
  string revised;

  //////////////////////////////////////////////////////////////////////////
  // Build string without whitespace characters, and covert to lower case 
  //////////////////////////////////////////////////////////////////////////
  for (i=0; i<original.length(); ++i) {
    if (!isspace(original[i])) {
      revised.push_back(tolower(original[i]));
    }
  }

  original = revised;
}


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::parseExpression
 *   Parses a formula, recursively calling itself to process pieces of 
 *   the equation with higher precedence that the current instance. 
 * 
 * Preconditions:
 *   none      
 *
 * Postconditions:
 *   Text that is parsed is erased from formula, and the commands created
 *   are placed in cmmnds.  'depth' is updated along with the computation 
 *   stack size
 *
 * Side Effects:
 *   none      
 *
 * Results:
 *   none 
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::EquationClass::parseExpression(
  tokenStruct&             token,
  string&                  formula, 
  parsingCommandContainer& cmmnds, 
  unsigned                 prvs_precedence, 
  unsigned&                depth  
  )
{
  formulaCommand new_command;
  tokenStruct    next_token;
  unsigned       next_precedence;

  //////////////////////////////////////////////////////////////////////////
  // Recursively call until a factor is found 
  ///////////////////////////////////////////////////////////////////////////
  if (prvs_precedence == HIGHEST_PRECEDENCE) {
    parseFactor(token, formula, cmmnds, depth);
  }
  else if ((prvs_precedence == UNARY_PRCDNC) && 
           (unaryToken[token.token] != COMMAND_INVALID)) { 
    next_token = token;
    getToken(formula, token);
    parseFactor(token, formula, cmmnds, depth);

    new_command = MAKE_COMMAND( unaryToken[next_token.token], 0 );
    cmmnds.push_back(new_command);
  }
  else {
    parseExpression( token, formula, cmmnds, prvs_precedence-1, depth );
  }
  
  //////////////////////////////////////////////////////////////////////////
  // Process if this is in infix token with the current precedence level 
  //////////////////////////////////////////////////////////////////////////
  while((infixToken[token.token] != COMMAND_INVALID) &&
        (tokenPriority[token.token] == prvs_precedence)) {

    if (token.token == TOKEN_QUESTION_MARK) {
      parseQuestionMark( token, formula, cmmnds, depth );
    }
    else {
      next_token = token;
      getToken(formula, token);
      next_precedence = tokenPriority[token.token];

      parseExpression( token, formula, cmmnds, prvs_precedence-1, depth );

      switch (next_token.token) {
        case TOKEN_COLON:
          break;

        default:
          new_command = MAKE_COMMAND( infixToken[next_token.token], 0 );
          cmmnds.push_back(new_command);
          changeDepth( -1, depth );
          break;
      }
    }
  }

}


/*-
 *-----------------------------------------------------------------------
 * parseQuestionMark
 *   Support function for parseExpression.  Parses the right side of the
 *   conditional (?) operator.  
 * 
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   Text that is parsed is erased from formula, and the commands created
 *   are placed in cmmnds.  'depth' is updated along with the computation 
 *   stack size.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::EquationClass::parseQuestionMark(
  tokenStruct&             token,
  string&                  formula,  
  parsingCommandContainer& commands,
  unsigned&                depth  
  )
{
  parsingCommandContainer true_commands;
  parsingCommandContainer false_commands;
  tokenStruct             next_token;
  formulaCommand          new_command;
  unsigned                start_depth;

  //////////////////////////////////////////////////////////////////////////
  // Parse the 'true' section of the conditional 
  //////////////////////////////////////////////////////////////////////////
  start_depth = depth; 
  next_token = token;
  getToken(formula, token);
  parseExpression( token, formula, true_commands, LOWEST_PRECEDENCE, depth );

  if (token.token != TOKEN_COLON) {
    string error_message = "Expecting ':' after '?' operator";
    throw(error_message); 
  }

  if (start_depth+1 != depth) {
    string error_message = "Not enough factors in true section of conditional";
    throw(error_message); 
  }
  changeDepth( -1, depth );

  //////////////////////////////////////////////////////////////////////////
  // Parse the 'false' section of the conditional 
  //////////////////////////////////////////////////////////////////////////
  next_token = token;
  getToken(formula, token);
  parseExpression( token, formula, false_commands, LOWEST_PRECEDENCE, depth );

  if (start_depth+1 != depth) {
    string error_message = "Not enough factors in false section of conditional";
    throw(error_message); 
  }

  //////////////////////////////////////////////////////////////////////////
  // Add the branch commands and the commands for the true and false clauses 
  //////////////////////////////////////////////////////////////////////////
  new_command = MAKE_COMMAND( COMMAND_BRANCH_IF_FALSE, true_commands.size()+1 );
  commands.push_back(new_command);
  changeDepth( -1, depth );
  copy( true_commands.begin(), true_commands.end(), inserter(commands, 
        commands.end()) );
  new_command = MAKE_COMMAND( COMMAND_BRANCH, false_commands.size() );
  commands.push_back(new_command);
  copy( false_commands.begin(), false_commands.end(), inserter(commands, 
    commands.end()) );

}


/*-
 *-----------------------------------------------------------------------
 * parseFactor
 *   Support function for parseFormula.  A factor evaluates to a value
 *   (possibly recursively) that operators can evaluate.  Examples are
 *   integers, functions, 'p0', and 'c0'.
 * 
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   Text that is parsed is erased from formula, and the commands created
 *   are placed in cmmnds.  'depth' is updated along with the computation 
 *   stack size.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
bool
RngDecisionTree::EquationClass::parseFactor(
  tokenStruct&             token,
  string&                  formula,  
  parsingCommandContainer& commands,
  unsigned&                depth  
  )
{
  formulaCommand new_command;
  tokenEnum      next_token;
  bool           isFactor = true;

  switch (token.token) {

    case TOKEN_LEFT_PAREN:  
      getToken(formula, token); 
      parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
      if (token.token == TOKEN_RIGHT_PAREN) {
        getToken(formula, token);
      }
      else {
        string error_message = "Expecting right parenthesis at '" + formula + 
          "'";
        throw(error_message); 
      }
      break;
 
    case TOKEN_INTEGER:
      new_command = MAKE_COMMAND( COMMAND_PUSH_CONSTANT, token.number ); 
      commands.push_back(new_command);
      changeDepth( 1, depth );
      getToken(formula, token); 
      break;

    default:
      //////////////////////////////////////////////////////////////////////
      // Is it a variable token 
      //////////////////////////////////////////////////////////////////////
      if (variableToken[token.token] != COMMAND_INVALID) {
        new_command = MAKE_COMMAND( variableToken[token.token], token.number );
        commands.push_back(new_command);
        changeDepth( 1, depth );
        getToken(formula, token); 
      }

      //////////////////////////////////////////////////////////////////////
      // Functions which take two or more operands 
      //////////////////////////////////////////////////////////////////////
      else if (functionToken[token.token] != COMMAND_INVALID) {
  
        next_token = token.token;

        getToken(formula, token); 
        if (token.token != TOKEN_LEFT_PAREN) {
          string error_message = "Expecting left parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        getToken(formula , token); 
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        if (token.token != TOKEN_COMMA) {
          string error_message = "Function requires at least two operands"; 
          throw(error_message);
        }

        getToken(formula, token);
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);

        while (token.token == TOKEN_COMMA) {
          new_command = MAKE_COMMAND( functionToken[next_token], 0 );
          commands.push_back(new_command);
          changeDepth( -1, depth );

          getToken(formula, token);
          parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        }

        if (token.token != TOKEN_RIGHT_PAREN) {
          string error_message = "Expecting right parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        new_command = MAKE_COMMAND( functionToken[next_token], 0 );
        commands.push_back(new_command);
        changeDepth( -1, depth );
        getToken(formula , token); 
      }
      //////////////////////////////////////////////////////////////////////
      // Functions which take exactly one operand 
      //////////////////////////////////////////////////////////////////////
      else if (oneValFunctionToken[token.token] != COMMAND_INVALID) {
  
        next_token = token.token;

        getToken(formula, token); 
        if (token.token != TOKEN_LEFT_PAREN) {
          string error_message = "Expecting left parenthesis at '" + 
          formula + "'";
          throw(error_message);
        }

        getToken(formula, token);
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);

        if (token.token != TOKEN_RIGHT_PAREN) {
          string error_message = "Expecting right parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        new_command = MAKE_COMMAND( oneValFunctionToken[next_token], 0 );
        commands.push_back(new_command);
        getToken(formula, token); 
      }
      //////////////////////////////////////////////////////////////////////
      // Functions which take exactly two operands 
      //////////////////////////////////////////////////////////////////////
      else if (twoValFunctionToken[token.token] != COMMAND_INVALID) {
  
        next_token = token.token;

        getToken(formula, token); 
        if (token.token != TOKEN_LEFT_PAREN) {
          string error_message = "Expecting left parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        getToken(formula, token); 
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        if (token.token != TOKEN_COMMA) {
          string error_message = "Function requires exactly two operands"; 
          throw(error_message);
        }

        getToken(formula, token);
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);

        if (token.token != TOKEN_RIGHT_PAREN) {
          string error_message = "Expecting right parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        new_command = MAKE_COMMAND( twoValFunctionToken[next_token], 0 );
        commands.push_back(new_command);
        changeDepth( -1, depth );
        getToken(formula, token); 
      }
      //////////////////////////////////////////////////////////////////////
      // Functions which takes >=2 operands and the function 
      // evaluator evaluates all operands at once
      //////////////////////////////////////////////////////////////////////
      else if (manyValFunctionToken[token.token] != COMMAND_INVALID) {

        unsigned nmbr_operands = 0;

        next_token = token.token;

        getToken(formula, token); 
        if (token.token != TOKEN_LEFT_PAREN) {
          string error_message = "Expecting left parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        getToken(formula , token); 
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        nmbr_operands++; 
        if (token.token != TOKEN_COMMA) {
          string error_message = "Function requires at least two operands"; 
          throw(error_message);
        }

        getToken(formula, token);
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        nmbr_operands++; 

        while (token.token == TOKEN_COMMA) {
          getToken(formula, token);
          parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
          nmbr_operands++; 
        }

        if (token.token != TOKEN_RIGHT_PAREN) {
          string error_message = "Expecting right parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        new_command = MAKE_COMMAND( manyValFunctionToken[next_token], 
          nmbr_operands );
        commands.push_back(new_command);
        changeDepth( -(nmbr_operands-1), depth );
        getToken(formula , token); 
      }
      //////////////////////////////////////////////////////////////////////
      // Handle to median function separately
      //    median(a,b,...)
      //////////////////////////////////////////////////////////////////////
      else if (token.token == TOKEN_MEDIAN) {

        unsigned nmbr_operands = 0;

        next_token = token.token;

        getToken(formula, token); 
        if (token.token != TOKEN_LEFT_PAREN) {
          string error_message = "Expecting left parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        getToken(formula , token); 
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        nmbr_operands++; 
        if (token.token != TOKEN_COMMA) {
          string error_message = "Median requires at least three operands"; 
          throw(error_message);
        }

        getToken(formula , token); 
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        nmbr_operands++; 
        if (token.token != TOKEN_COMMA) {
          string error_message = "Median requires at least three operands"; 
          throw(error_message);
        }

        getToken(formula, token);
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        nmbr_operands++; 

        while (token.token == TOKEN_COMMA) {

          getToken(formula, token);
          parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
          nmbr_operands++; 
        }

        if (token.token != TOKEN_RIGHT_PAREN) {
          string error_message = "Expecting right parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        new_command = MAKE_COMMAND( COMMAND_MEDIAN, nmbr_operands );
        commands.push_back(new_command);
        changeDepth( -(nmbr_operands-1), depth );
        getToken(formula , token); 
      }
      //////////////////////////////////////////////////////////////////////
      // Handle to rotate function separately
      //    rotate(val,num,pos,length)
      //////////////////////////////////////////////////////////////////////
      else if (token.token == TOKEN_ROTATE) {

        next_token = token.token;

        getToken(formula, token); 
        if (token.token != TOKEN_LEFT_PAREN) {
          string error_message = "Expecting left parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        getToken(formula, token); 
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        if (token.token != TOKEN_COMMA) {
          string error_message = "Function requires four operands:  rotate(value, number, position, bitwidth)";
          throw(error_message);
        }

        getToken(formula, token);
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        if (token.token != TOKEN_COMMA) {
          string error_message = "Function requires four operands:  rotate(value, number, position, bitwidth)";
          throw(error_message);
        }

        getToken(formula, token);
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);
        if (token.token != TOKEN_COMMA) {
          string error_message = "Function requires four operands:  rotate(value, number, position, bitwidth)";
          throw(error_message);
        }

        getToken(formula, token);
        parseExpression(token, formula, commands, LOWEST_PRECEDENCE, depth);

        if (token.token != TOKEN_RIGHT_PAREN) {
          string error_message = "Expecting right parenthesis at '" + 
            formula + "'";
          throw(error_message);
        }

        new_command = MAKE_COMMAND( COMMAND_ROTATE, 0 );
        commands.push_back(new_command);

        changeDepth( -3, depth );
        getToken(formula , token); 
      }     
      //////////////////////////////////////////////////////////////////////
      // Not a valid factor 
      //////////////////////////////////////////////////////////////////////
      else {
        isFactor = false;
      }

      break;
  }

  return(isFactor);
}


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::getToken
 *   Support function for parseFormula.  Parses a string and gives the   
 *   next token in the string.
 * 
 * Preconditions:
 *   none 
 *
 * Postconditions:
 *   The token found is written to the token structure, and the 
 *   characters making up the token are removed from the expression.
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   none 
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::EquationClass::getToken(
  string&      expression, 
  tokenStruct& token 
  )
{
  //////////////////////////////////////////////////////////////////////////
  // Local variables 
  //////////////////////////////////////////////////////////////////////////
  string   token_string;

  map<string, tokenEnum>::iterator crrnt_dlmtr; 
  map<string, tokenEnum>::iterator end_dlmtr; 
  size_t dlmtr_lctn; 
  size_t minimum_dlmtr_lctn;
 
  vector<map<string, tokenEnum>::iterator> found_dlmtr_cntnr; 
  vector<map<string, tokenEnum>::iterator>::iterator crrnt_found_dlmtr; 
  vector<map<string, tokenEnum>::iterator>::iterator end_found_dlmtr; 
  map<string, tokenEnum>::iterator found_dlmtr; 

  //////////////////////////////////////////////////////////////////////////
  // Check for end of string 
  //////////////////////////////////////////////////////////////////////////
  if (expression.length() == 0) {
    token.token = TOKEN_END;
    minimum_dlmtr_lctn = 0; 
  }
  else {
    //////////////////////////////////////////////////////////////////////////
    // Find the earliest delimiter token 
    //////////////////////////////////////////////////////////////////////////
    token.token = LAST_TOKEN_INDEX;
    minimum_dlmtr_lctn = expression.size();

    for( crrnt_dlmtr = delimiter.begin(),
         end_dlmtr   = delimiter.end();
         crrnt_dlmtr != end_dlmtr;  
         ++crrnt_dlmtr ) { 

      dlmtr_lctn = expression.find( (*crrnt_dlmtr).first, 0 );    
      if (dlmtr_lctn != string::npos) { 
        if (dlmtr_lctn < minimum_dlmtr_lctn) {
          minimum_dlmtr_lctn = dlmtr_lctn;
          found_dlmtr_cntnr.clear();
        }
        if (dlmtr_lctn <= minimum_dlmtr_lctn) {
          found_dlmtr_cntnr.push_back(crrnt_dlmtr);
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // If next token is a delimiter, determine the longest delimiter which 
    // matches. 
    //////////////////////////////////////////////////////////////////////////
    if (minimum_dlmtr_lctn == 0) {
      found_dlmtr = *found_dlmtr_cntnr.begin();
      for( crrnt_found_dlmtr = found_dlmtr_cntnr.begin(),
           end_found_dlmtr  = found_dlmtr_cntnr.end();
           crrnt_found_dlmtr != end_found_dlmtr;  
           ++crrnt_found_dlmtr) { 

        if ((*(*crrnt_found_dlmtr)).first.size()>(*found_dlmtr).first.size()) {
          found_dlmtr = *crrnt_found_dlmtr;
        }
      }

      token.token = (*found_dlmtr).second;
      minimum_dlmtr_lctn += (*found_dlmtr).first.size();
    }
    //////////////////////////////////////////////////////////////////////////
    // The next token is not a delimiter, the token is contained in the string
    // between the beginning of the expression and delimiter.  Now parse 
    // this string for a token.
    //////////////////////////////////////////////////////////////////////////
    else {

      token_string = expression.substr(0, minimum_dlmtr_lctn);

      //////////////////////////////////////////////////////////////////////////
      // Check for an integer 
      //////////////////////////////////////////////////////////////////////////
      map<string, tokenEnum>::iterator found_string;
      int    number;    
      bool   is_integer;

      is_integer = getInteger( token_string, number ); 
      if (is_integer) {
        token.token  = TOKEN_INTEGER;
        token.number = number; 
      }
      else {

        ////////////////////////////////////////////////////////////////////////
        // Check if token string is a function 
        ////////////////////////////////////////////////////////////////////////
        found_string = function.find( token_string );
        if (found_string != function.end()) {
          token.token = (*found_string).second; 
        }

        ////////////////////////////////////////////////////////////////////////
        // Check if token string is a variable name (which doesn't require a
        // following integer)
        ////////////////////////////////////////////////////////////////////////
        found_string = variable.find( token_string );
        if (found_string != variable.end()) {  
          token.token = (*found_string).second;
          token.number = number; 
        }
   
        if (token.token == LAST_TOKEN_INDEX) {
          //////////////////////////////////////////////////////////////////////
          // Next, check for variable name (such as "p0" or "c2")
          //////////////////////////////////////////////////////////////////////

          string variable_name; 
          int integer_location = -1;

          for (unsigned i=0; i<token_string.size(); ++i) {
            if (isdigit(token_string[i])) { 
              integer_location = i;
              break;
            }
          }

          if (integer_location >= 0) {
         
            is_integer = getInteger( token_string.substr(integer_location,
              token_string.length()), number );

            if (is_integer) {
              found_string = variable.find( token_string.substr(0, 
                integer_location) );
              if (found_string != variable.end()) {
                token.token = (*found_string).second;
                token.number = number; 
              }
            }
          } 
        } 
      } 
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Exit with error message if no token was found 
  //////////////////////////////////////////////////////////////////////////
  if (token.token == LAST_TOKEN_INDEX) {
    string error_message = "Invalid symbol at '" + expression + "'";
    throw(error_message); 
  }
 
  //////////////////////////////////////////////////////////////////////////
  // Erase portion of string that was parsed 
  //////////////////////////////////////////////////////////////////////////
  expression.erase(0, minimum_dlmtr_lctn ); 

}


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::getInteger
 *   Determine if string is an integer and gives the integer if
 *   one is found. 
 * 
 * Preconditions:
 *   none     
 *
 * Postconditions:
 *   If found, the integer is stored in number otherwise number is 
 *   undetermined. 
 *
 * Side Effects:
 *   none     
 *
 * Results:
 *   Returns true if given expression is an integer, false if not 
 *-----------------------------------------------------------------------
 */
bool 
RngDecisionTree::EquationClass::getInteger(
  const string& expression, 
  int& number 
  )
{
  string   nmbr;
  unsigned index = 0;
  bool     is_number = true;
  bool     number_found = false;

  while ((index < expression.size()) && (is_number)) {

    if ((expression[index]>='0') && (expression[index]<='9')) {
      nmbr += expression[index];
    }
    else {
      is_number = false;
    }

    ++index;
  }
 
  if (is_number) {
    number_found = true;
    number = atoi(nmbr.c_str());
  }

  return(number_found);
}


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::changeDepth
 *   The depth is the size of the formula computation stack.  This 
 *   proceedure updates the depth and increases the size of the stack, if
 *   needed.  
 * 
 * Preconditions:
 *   none     
 *
 * Postconditions:
 *   detph+=change, stack size may be increased 
 *
 * Side Effects:
 *   none     
 *
 * Results:
 *   none     
 *-----------------------------------------------------------------------
 */
void 
RngDecisionTree::EquationClass::changeDepth(
  const int change, 
  unsigned& depth 
  )
{
  if (((int)depth+change)<0)
  {
    throw("Not enough factors for given operators");
  }

  depth += change;
  stack.growIfNeeded(depth);
}


/*-
 *-----------------------------------------------------------------------
 * beginIterableDT
 *      rewinds the DT file to the beginning.
 * 
 * Preconditions:
 *      This must be a DT object that gets its DT instantiations
 *      from a file.
 *
 * Postconditions:
 *      next DT is read in.
 *
 * Side Effects:
 *      changes internal structures.
 *
 * Results:
 *      returns nil
 *
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::beginIterableDT()
{
  //////////////////////////////////////////////////////////////////////
  // first make sure this is a DT from file object
  //////////////////////////////////////////////////////////////////////
  if (!iterable()) {
    error("ERROR: can't call beginIterableDT() for non-file DT");
  }

  //////////////////////////////////////////////////////////////////////
  // Rewind file pointer to the beginning of the files and read the 
  // number decision trees 
  //////////////////////////////////////////////////////////////////////
  dtFile->rewind();
  dtFile->read(numDTs, "num DTs");
  dtNum = -1;

  //////////////////////////////////////////////////////////////////////
  // If necessary, use index file to seek to first utterance 
  //////////////////////////////////////////////////////////////////////
  if (firstDT != 0) {
    seek( firstDT ); 
  }
 
  //////////////////////////////////////////////////////////////////////
  // Read in the first tree 
  //////////////////////////////////////////////////////////////////////
  nextIterableDT();
}


/*-
 *-----------------------------------------------------------------------
 * nextIterableDT
 *      reads in the next DT from the file.
 * 
 * Preconditions:
 *      This must be a DT object that gets its DT instantiations
 *      from a file.
 *
 * Postconditions:
 *      next DT is read in.
 *
 * Side Effects:
 *      changes internal structures.
 *
 * Results:
 *      returns nil
 *
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::nextIterableDT()
{
  int    readDTNum;

  //////////////////////////////////////////////////////////////////////
  // first make sure this is a DT from file object
  //////////////////////////////////////////////////////////////////////
  if (!iterable()) {
    error("ERROR: can't call nextIterableDT() for non-file DT");
  }

  //////////////////////////////////////////////////////////////////////
  // Increment tree index and make sure the file matches 
  //////////////////////////////////////////////////////////////////////
  dtNum++;

  dtFile->read(readDTNum, "DT num");
  if (readDTNum != dtNum) {
    error("ERROR: reading from file '%s', expecting DT number %d but got number %d\n",
      dtFileName.c_str(), dtNum, readDTNum);
  }

  //////////////////////////////////////////////////////////////////////
  // next, delete old DT if it exists.
  //////////////////////////////////////////////////////////////////////
  if (root != NULL) {
    destructorRecurse(root);
    delete root;
    root = NULL;
  }

  // read in the rest of the DT.
  dtFile->read(curName,"cur name");
  dtFile->read(_numFeatures,"num feats");
  if (_numFeatures <= 0)
    error("ERROR: reading dynamic decision tree '%s' with current name '%s' from file '%s', but decision tree must have > 0 features",
	  name().c_str(),
	  curName.c_str(),
	  dtFile->fileName());
  root = new Node;
  readRecurse(*dtFile,*root);
}


/*-
 *-----------------------------------------------------------------------
 * writeIndexFile 
 *   Writes an index file for a iterable DT 
 * 
 * Preconditions:
 *   The file containing the per-utterance DTs should have already been 
 *   opened using read
 *
 * Postconditions:
 *   File will be written
 *
 * Side Effects:
 *   File pointer for per-utterance DTs is moved to the end of the file.
 *   If it is to be used again beginIterableDT must be called. 
 *
 * Results:
 *   none
 *
 * TODO: Move DT parsing into a subroutine
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::writeIndexFile()
{
  int      numDTs;
  unsigned position;
  int      i;

  assert(iterable()); 

  //////////////////////////////////////////////////////////////////////////
  // Open index file for writing 
  //////////////////////////////////////////////////////////////////////////
  string outputFileName = dtFile->fileName() + DTFileExtension;
  oDataStreamFile outFile(outputFileName.c_str(), true);

  printf("Writing file %s\n", outputFileName.c_str() ); 
 
  //////////////////////////////////////////////////////////////////////////
  // Write number of trees 
  //////////////////////////////////////////////////////////////////////////
  dtFile->rewind();
  dtFile->read(numDTs, "num DTs");
  outFile.writeInt(numDTs, "num DTs");

  //////////////////////////////////////////////////////////////////////////
  // Write indices 
  //////////////////////////////////////////////////////////////////////////
  for(i=0; i<numDTs; i++) {

    position = dtFile->ftell();

    if (position < 0) {
      error("ERROR: problem reading from file '%s'\n", dtFileName.c_str() );
    }

    outFile.writeInt(position, "Decision tree index");

    ////////////////////////////////////////////////////////////////////////
    // Parse file to the beginning of the next tree 
    ////////////////////////////////////////////////////////////////////////
    dtFile->read(dtNum, "DT num");
    if (i != dtNum) {
      error("ERROR: reading from file '%s', expecting DT number %d but got number %d\n",
        dtFileName.c_str(), dtNum, i);
    }

    // next, delete old DT if it exists.
    if (root != NULL) {
      destructorRecurse(root);
      delete root;
      root = NULL;
    }

    // read in the rest of the DT.
    dtFile->read(curName,"cur name");
    dtFile->read(_numFeatures,"num feats");
    if (_numFeatures <= 0) {
      error("ERROR: reading dynamic decision tree '%s' with current name '%s' from file '%s', but decision tree must have > 0 features", 
        name().c_str(), curName.c_str(), dtFile->fileName());
    } 


    root = new Node;
    readRecurse(*dtFile,*root);
  }

}


/*-
 *-----------------------------------------------------------------------
 * write:
 *      writes to from a file.
 * 
 * Preconditions:
 *      object should be filled in
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      nothing.
 *
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.nl();
  if (!iterable()) {
    os.write(_numFeatures,"RngDecisionTree:: write numFeatures");
    os.nl();
    writeRecurse(os,root,0);
  } else {
    // just write out the file name
    os.write(dtFileName,"RngDecisionTree:: write numFeatures");
    os.nl();
  }
}


/*-
 *-----------------------------------------------------------------------
 * writeRecurse:
 *      continues (from write) writing out to a file 
 * 
 * Preconditions:
 *      object should be filled in
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      the current node it read
 *
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::writeRecurse(oDataStreamFile& os,
			      RngDecisionTree::Node *n,
			      const int depth)
{
  if (n->nodeType == LeafNodeVal) {
    os.space(depth*2);
    os.write(-1);
    os.write(n->ln_v().value);
    os.nl();
  } else if (n->nodeType == LeafNodeEquation) {
    os.space(depth*2);
    os.write(-1);
    // TODO: equation guy should be able to write itself.
    os.write("equation");
    os.nl();
  } else if (n->nodeType == NonLeafNodeArray) {
    os.space(depth*2);
    os.write(n->nln_a().ftr,"writeRecurse, ftr");
    os.write(n->nln_a().children.size(),"writeRecurse, numsplits");
    if (n->nln_a().children.size() > 1) {
      os.write(n->nln_a().base);
      os.write("...");
      os.write(n->nln_a().base+n->nln_a().children.size()-2);
    }
    os.write(RNG_DECISION_TREE_DEF_STR);
    os.nl();
    for (unsigned i=0;i<n->nln_a().children.size();i++) 
      writeRecurse(os,&(n->nln_a().children[i]),depth+1);

  } else if (n->nodeType == NonLeafNodeHash) {
    os.space(depth*2);
    os.write(n->nln_h().ftr,"writeRecurse, ftr");
    os.write(n->nln_h().children.size(),"writeRecurse, numsplits");
    os.write("hashed_int_vals");
    os.write(RNG_DECISION_TREE_DEF_STR);    
    os.nl();
    for (unsigned i=0;i<n->nln_h().children.size();i++) 
      writeRecurse(os,&(n->nln_h().children[i]),depth+1);
  } else {
    os.space(depth*2);
    os.write(n->nln_r().ftr,"writeRecurse, ftr");
    os.write(n->nln_r().children.size()+1,"writeRecurse, numsplits");
    for (unsigned i=0;i<n->nln_r().children.size();i++) {
      os.write(n->nln_r().children[i].rng.rangeStr());
    }
    os.write(RNG_DECISION_TREE_DEF_STR);
    os.nl();
    for (unsigned i=0;i<n->nln_r().children.size();i++) 
      writeRecurse(os,&(n->nln_r().children[i].nd),depth+1);
    writeRecurse(os,n->nln_r().def,depth+1);
  }
}



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                    Querying the decision tree
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 *  queryRecurse
 *      Exact same routine as above but for arrays of random variables.
 *      This code is duplicated here for simplicity and efficiency.
 *      NOTE: Any update here should also be done in the
 *      version that takes arrays of integers above.
 *
 *   Note that this routine assumes that all RV's may be safely cast
 *   into DiscRVs. If this is not the case, arbitrary results will
 *   follow. If 'drv' is non-NULL, then any DT formulas will have
 *   access to features of the current child random variable (such
 *   as child's cardinality, frame number, and so on). Otherwise, if
 *   'drv' is NULL, then both child cardinality and child frame will
 *   evaluate to 0.
 *
 * Preconditions:
 *      same as above
 *
 * Postconditions:
 *      same as above
 *
 * Side Effects:
 *      same as above
 *
 * Results:
 *      as above
 *
 *-----------------------------------------------------------------------
 */
leafNodeValType RngDecisionTree::queryRecurse(const vector < RV* >& arr,
					      RngDecisionTree::Node *n,
					      const RV* rv)
{

  // include label here to implement tail recursion with gotos.
 topOfRoutine:

  if (n->nodeType == LeafNodeVal) {
    return n->ln_v().value;
  } else if  (n->nodeType == LeafNodeEquation) {
    leafNodeValType answer;
    answer = n->ln_e().equation.evaluateFormula( arr, rv );
    return(answer);
  } else if (n->nodeType == NonLeafNodeArray) {
    assert ( n->nln_a().ftr < int(arr.size()) );
    const int val = RV2DRV(arr[n->nln_a().ftr])->val - n->nln_a().base;
    
    if (val >= 0 &&  (val+1) < (int)n->nln_a().children.size()) {
      n = & n->nln_a().children[val];
    } else {
      // default case
      n = & n->nln_a().children[n->nln_a().children.size()-1];
    }
    goto topOfRoutine;
  } else if (n->nodeType == NonLeafNodeHash) {
    assert ( n->nln_h().ftr < int(arr.size()) );
    const int val = RV2DRV(arr[n->nln_h().ftr])->val;
    Node** tmp = n->nln_h().nodeMapper.find(val);
    if (tmp != NULL)
      n = *tmp;
    else
      n = & n->nln_h().children[n->nln_h().children.size()-1];
    goto topOfRoutine;
  } else if (n->nodeType == NonLeafNodeRngs) {
    assert ( n->nln_r().ftr < int(arr.size()) );
    const int val = RV2DRV(arr[n->nln_r().ftr])->val;

    // use a switch to knock off the short cases for
    // which we use simple linear search with little bookkeeping.
    if (n->nln_r().children.size() < DT_SPLIT_SORT_THRESHOLD
	|| 
	!n->nln_r().ordered) {
      for (unsigned i=0;i<n->nln_r().children.size();i++) {
	// Do a linear search.
	if (n->nln_r().children[i].rng.contains(val)) {
	  // return queryRecurse(arr,n->nonLeafNode.children[i],rv);
	  n = &(n->nln_r().children[i].nd);
	  goto topOfRoutine;
	}
      }
    } else {
      // eliminate simple boundary conditions.
      const int maxRngNum = n->nln_r().children.size()-1;
      if ((n->nln_r().children[0].rng) > val)
	goto failedQuery;
      if ((n->nln_r().children[maxRngNum].rng) < val)
	goto failedQuery;

      // do binary search
      int l=0,u=maxRngNum;
      int m=0;
      while (l<=u) {
	// all these are conditional on the val being contained in the rng.
	// rngs[l] <= val && val <= rngs[u]  
	m = (l+u)/2; 
	if ((n->nln_r().children[m].rng) > val) 
	  // rngs[l] <= val && val < rngs[m]
	  u = m-1;
	// rngs[l] <= val && val <= rngs[u]
	else if ((n->nln_r().children[m].rng) < val)
	  // rngs[m] < val && val < rngs[u]
	  l=m+1;
	// rngs[l] < val && val < rngs[u]
	else {
	  // found potential range that might contain value
	  // since neither val < rng nor val > rng. 
	  if (n->nln_r().children[m].rng.contains(val)) {
	    // return queryRecurse(arr,n->nonLeafNode.children[m],rv);
	    n = &(n->nln_r().children[m].nd);
	    goto topOfRoutine;
	  }
	  // linearly search up and down for all ranges that have overlapping
	  // boundaries, using boundariesOverlap() and query them as well since
	  // they are treated as equal.
	  if (m < maxRngNum) {
	    int im=m;
	    do {
	      if (n->nln_r().children[im].rng.boundariesOverlap(n->nln_r().children[im+1].rng)) {
		if (n->nln_r().children[im+1].rng.contains(val)) {
		  n = &(n->nln_r().children[im+1].nd);
		  goto topOfRoutine;
		}
	      } else {
		// boundaries don't ovelrap so no chance of finding it here.
		break;
	      }
	    } while (++im < maxRngNum);
	  }
	  if (m > 0) {
	    int im=m;
	    do {
	      if (n->nln_r().children[im].rng.boundariesOverlap(n->nln_r().children[im-1].rng)) {
		if (n->nln_r().children[im-1].rng.contains(val)) {
		  n = &(n->nln_r().children[im-1].nd);
		  goto topOfRoutine;
		}
	      } else {
		// boundaries don't ovelrap so no chance of finding it here.
		break;
	      }
	    } while (--im > 0);
	  }
	  break;
	}
      }
    }
    
    // failed lookup, so return must be the default one.
  failedQuery:

    // failed lookup, so return must be the default one. 
    // return queryRecurse(arr,n->nonLeafNode.children[n->nonLeafNode.rngs.size()],rv);
    n = n->nln_r().def;
    goto topOfRoutine;
  }
  assert ( 0 ) ;
  // this should never happen.
  return ~0x0;
}



/*-
 *-----------------------------------------------------------------------
 *  computeParentsSatisfyingChild()
 *      Count the number of parent assigments will satisfy the current observed child.
 *
 * Preconditions:
 *      Child is assumed to be assigned to (observed at)  an appropriate value.
 *      If any parents are observed, we assume they are already set to their observed
 *      values.
 *
 * Postconditions:
 *      child must be observed
 *
 * Side Effects:
 *      changes values of rvs
 *
 * Results:
 *      count in reference 'num'
 *
 *-----------------------------------------------------------------------
 */
void RngDecisionTree::computeParentsSatisfyingChild(
	    // input arguments
	    unsigned par, // parent number
	    vector <RV*> & parents, 
	    vector <RV*> & hiddenParents,
	    PackCliqueValue& hiddenParentPacker,
	    sArray < DiscRVType*>& hiddenNodeValPtrs,
	    RV* child,
	    // output arguments
	    sArray < unsigned >& packedParentVals,
	    unsigned& num)
{
  if (par == parents.size()) {
    unsigned val = query(child->allParents,child);
    if (val == RV2DRV(child)->val) {
      // we've got a hit.
      packedParentVals.growByNIfNeededAndCopy(2,(num+1)*hiddenParentPacker.packedLen());
      // printf("before pack"); printRVSetAndValues(stdout,parents);
      hiddenParentPacker.pack(hiddenNodeValPtrs.ptr,&packedParentVals.ptr[num*hiddenParentPacker.packedLen()]);
      // hiddenParentPacker.unpack(&packedParentVals.ptr[num*hiddenParentPacker.packedLen()],hiddenNodeValPtrs.ptr);
      // printf("after unpack"); printRVSetAndValues(stdout,parents);
      // number of entries increases 
      num++;
    }

    if (message(Max+5)) {
      // optionally print out diagnostics
      printf("Pr[%s(%d)=%d|",child->name().c_str(),child->frame(),
	     RV2DRV(child)->val);
      printRVSetAndValues(stdout,parents,false);
      printf("]=%d\n",(val == RV2DRV(child)->val));
    }
  } else {
    DiscRV*drv = RV2DRV(parents[par]);
    if (!drv->discreteObservedImmediate()) {
      HidDiscRV*hdrv = (HidDiscRV*)drv;
      for (hdrv->val = 0; hdrv->val < hdrv->cardinality; hdrv->val ++) {
	computeParentsSatisfyingChild(par+1,parents,hiddenParents,hiddenParentPacker,hiddenNodeValPtrs,
				      child,packedParentVals,num);
      }
    } else {
      // We assume any observed parents are already set to observed value.
      // ObsDiscRV*odrv = (ObsDiscRV*)drv;
      // odrv->setToObservedValue();
      computeParentsSatisfyingChild(par+1,parents,hiddenParents,hiddenParentPacker,hiddenNodeValPtrs,
				    child,packedParentVals,num);
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 *  computeParentsChildSatisfyingGrandChild
 *      Count the number of parent assigments giving a child that satisfies the observed grand child.
 *
 * Preconditions:
 *      grandChild is assumed to be assigned to (observed at)  an appropriate value.
 *      If any parents are observed, we assume they are already set to their observed
 *      values.
 *
 * Postconditions:
 *      child must be observed
 *
 * Side Effects:
 *      changes values of rvs
 *
 * Results:
 *      count in reference 'num'
 *
 *
 *-----------------------------------------------------------------------
 */
void RngDecisionTree::computeParentsChildSatisfyingGrandChild(
	    // input arguments
	    unsigned par, // parent number
	    vector <RV*> & parents, 
	    vector <RV*> & hiddenParents,
	    PackCliqueValue& hiddenParentPacker,
	    sArray < DiscRVType*>& hiddenNodeValPtrs,
	    RV* child,
	    RV* grandChild,
	    // output arguments
	    sArray < unsigned >& packedParentVals,
	    unsigned& num)
{
  if (par == parents.size()) {
    unsigned val = query(child->allParents,child);
    RV2DRV(child)->val = val;
    logpr cur_p;
    // parent of granGhild is now Guaranteed to be set. 
    grandChild->probGivenParents(cur_p);
    if (!cur_p.zero()) {
      // we've got a hit.
      packedParentVals.growByNIfNeededAndCopy(2,(num+1)*hiddenParentPacker.packedLen());
      hiddenParentPacker.pack(hiddenNodeValPtrs.ptr,&packedParentVals.ptr[num*hiddenParentPacker.packedLen()]);
      num++;
    }

    if (message(Max+5)) {
      // optionally print out diagnostics
      printf("Pr[%s(%d)=%d|%s(%d)=%d]=%f;",
	     grandChild->name().c_str(),grandChild->frame(),RV2DRV(grandChild)->val,
	     child->name().c_str(),child->frame(),RV2DRV(child)->val,
	     cur_p.valref());
      printf("Pr[%s(%d)=%d|",child->name().c_str(),child->frame(),
	     RV2DRV(child)->val);
      printRVSetAndValues(stdout,parents,false);
      printf("]=%d\n",(val == RV2DRV(child)->val));
    }
  } else {
    DiscRV*drv = RV2DRV(parents[par]);
    if (!drv->discreteObservedImmediate()) {
      HidDiscRV*hdrv = (HidDiscRV*)drv;
      for (hdrv->val = 0; hdrv->val < hdrv->cardinality; hdrv->val ++) {
	computeParentsChildSatisfyingGrandChild(par+1,parents,hiddenParents,hiddenParentPacker,
						hiddenNodeValPtrs,
						child,grandChild,
						packedParentVals,num);
      }
    } else {
      // We assume any observed parents are already set to observed value.
      // ObsDiscRV*odrv = (ObsDiscRV*)drv;
      // odrv->setToObservedValue();
      computeParentsChildSatisfyingGrandChild(par+1,parents,hiddenParents,hiddenParentPacker,
					      hiddenNodeValPtrs,
					      child,grandChild,
					      packedParentVals,num);
    }
  }
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                    Misc
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
//        test code Support
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#include <algorithm>

////////////////////////////////////////////////////////////////////
// Test formulas 
////////////////////////////////////////////////////////////////////

class TestRandomVariable : public DiscRV
{
  public:

    TestRandomVariable( RVInfo new_info, string new_name, unsigned new_card ) 
       : DiscRV(new_info) { cardinality = new_card; return; }
    void begin(logpr& p) { return; }
    bool next(logpr& p) { return(false); }
    RV* create() { return(NULL); } //*
};

bool RngDecisionTree::verifyParseError(
  string               formula
  )
{
  Node     node;
  bool     correct;

  new (&node.ln_e()) LeafNodeEquationStruct();

  try {
    node.ln_e().equation.parseFormula(formula);
    printf("Error not correctly caught in formula:'%s'\n", formula.c_str());
    correct = false;
  }
  catch (...) {
    printf("   Parse error correctly caught\n");
    correct = true;
  }
 
  return(correct);
}


bool RngDecisionTree::testFormula(
  string               formula,
  const vector< RV* >& variables, 
  RV*                  child, 
  unsigned             desired_answer 
  )
{
  Node     node;
  unsigned answer;
  bool     correct;

  new (&node.ln_e()) LeafNodeEquationStruct();

  try {
    node.ln_e().equation.parseFormula(formula);
  }
  catch( string error_message ){
    error("   PARSE ERROR: %s", error_message.c_str());
  }
  catch( const char * const error_message ) {
    error("   PARSE ERROR: %s", error_message );
  }

  answer = node.ln_e().equation.evaluateFormula( variables, child );
  printf("   Answer: %d   0x%x\n", answer, answer);

  if (answer == desired_answer) {
    correct = true;
  }
  else { 
    printf("====>**** NOT CORRECT!!!\n");
    correct = false;
  }

  return(correct);
}

void test_formula()
{
  RngDecisionTree dt;
  vector<RV*>     vars;
  string          formula;
  bool            correct;

  RVInfo::FeatureRange tmp_fr;
  RVInfo::ListIndex    tmp_li;
  vector< RVInfo::rvParent > tmp_switchingParents;
  vector<vector< RVInfo::rvParent > > tmp_conditionalParents;
  vector< CPT::DiscreteImplementaton > tmp_discImplementations;
  vector< MixtureCommon::ContinuousImplementation > tmp_contImplementations;
  vector< RVInfo::ListIndex > tmp_listIndices;
  vector< RVInfo::WeightInfo > tmp_rvWeightInfo;

  RVInfo dummy(
    0, 1, 0, 10, "test", "aa_info", RVInfo::t_discrete, RVInfo::d_hidden, 10,  
    tmp_fr, false, false, NULL, tmp_li, tmp_switchingParents, 
    tmp_conditionalParents, tmp_discImplementations, tmp_contImplementations,
    tmp_listIndices, tmp_rvWeightInfo 
    );

  TestRandomVariable p0(dummy, "p0", 8 );
  TestRandomVariable p1(dummy, "p1", 40 );
  TestRandomVariable p2(dummy, "p2", 2 );
  TestRandomVariable child(dummy, "child", 1000000 );

  p0.val  = 7;
  p1.val  = 3;
  p2.val  = 0;
  vars.push_back(&p0);
  vars.push_back(&p1);
  vars.push_back(&p2);

  correct = true;

  formula = "+";
  correct &= dt.verifyParseError(formula);

  formula = "()*()";
  correct &= dt.verifyParseError(formula);

  formula = "()/3";
  correct &= dt.verifyParseError(formula);

  formula = "(3!=4)";
  correct &= dt.testFormula( formula, vars, &child, 1);

  formula = "p0 != 7";
  correct &= dt.testFormula( formula, vars, &child, 0);

  formula = "!0!=!1";
  correct &= dt.testFormula( formula, vars, &child, 1);

  formula = "cc";
  correct &= dt.testFormula( formula, vars, &child, 1000000 );

  formula = "cardinality_child";
  correct &= dt.testFormula( formula, vars, &child, 1000000 );

  formula = "cardinality_parent_0+cp0";
  correct &= dt.testFormula( formula, vars, &child, 16 );

  formula = "parent_1+p2";
  correct &= dt.testFormula( formula, vars, &child, 3 );

  formula = "parent_minus_one_0+pmo1";
  correct &= dt.testFormula( formula, vars, &child, 8 );

  formula = "parent_plus_one_0+ppo1";
  correct &= dt.testFormula( formula, vars, &child, 12 );

  formula = "max_value_child+mc";
  correct &= dt.testFormula( formula, vars, &child, (2*1000000)-2 );

  formula = "max_value_parent_2+mp0";
  correct &= dt.testFormula( formula, vars, &child, 8 );

  formula = "round_divide(11,7)";
  correct &= dt.testFormula( formula, vars, &child, 2 );

  formula = "10~/7";
  correct &= dt.testFormula( formula, vars, &child, 1 );

  formula = "floor_divide(11,7)";
  correct &= dt.testFormula( formula, vars, &child, 1 );

  formula = "11/.7";
  correct &= dt.testFormula( formula, vars, &child, 1 );

  formula = "ceil_divide (13,7)";
  correct &= dt.testFormula( formula, vars, &child, 2 );

  formula = "ceil_divide (14,7)";
  correct &= dt.testFormula( formula, vars, &child, 2 );

  formula = "15./7";
  correct &= dt.testFormula( formula, vars, &child, 3 );

  formula = "1+1";
  correct &= dt.testFormula( formula, vars, &child, 2 );

  formula = "(1+1)";
  correct &= dt.testFormula( formula, vars, &child, 2 );

  formula = "2*3+4";
  correct &= dt.testFormula( formula, vars, &child, 10 );

  formula = "2+3*4";
  correct &= dt.testFormula( formula, vars, &child, 14 );

  formula = "1+2+3+4+5";
  correct &= dt.testFormula( formula, vars, &child, 15 );

  formula = "3+(4==4)*5^2|6";
  correct &= dt.testFormula( formula, vars, &child, 30 );

  formula = "((p0&p1)|cp0)";
  correct &= dt.testFormula( formula, vars, &child, 11 );

  formula = "(cp1/cp0)*p1";
  correct &= dt.testFormula( formula, vars, &child, 15 );

  formula = "cp2^5";
  correct &= dt.testFormula( formula, vars, &child, 32 );

  formula = "7>3";
  correct &= dt.testFormula( formula, vars, &child, 1 );

  formula = "(7<3)";
  correct &= dt.testFormula( formula, vars, &child, 0 );

  formula = "(3>=3)&&(3==3)&&(3<=3)";
  correct &= dt.testFormula( formula, vars, &child, 1 );

  formula = "p0 || (27<2) ||  (229<=1)";
  correct &= dt.testFormula( formula, vars, &child, 1 );

  formula = "(MIN(CP0,CP1,CP2))";
  correct &= dt.testFormula( formula, vars, &child, 2 );

  formula = "2+max(p0,  P1,p2)*3 + 4";
  correct &= dt.testFormula( formula, vars, &child, 27 );

  formula = "  mod  (  8 , 3 ) ";
  correct &= dt.testFormula( formula, vars, &child, 2 );

  formula = "((12>3)?(p0*p1):(cp0+p2))";
  correct &= dt.testFormula( formula, vars, &child, 21 );

  formula = "12<3 ? p1/0 : cp0+p2";
  correct &= dt.testFormula( formula, vars, &child, 8 );

  formula = " 1 ? 27==p0 ? cp0+cp2 :  p0+2 :7/0 ";
  correct &= dt.testFormula( formula, vars, &child, 9 );

  formula = "  (  max  (  ( 12/  4) , 14,p2  ) + ( cp1   -p1))*   p0+27";
  correct &= dt.testFormula( formula, vars, &child, 384 );

  formula = " Mod(  cp1  - p1 , 17 + P0  ) - 4";
  correct &= dt.testFormula( formula, vars, &child, 9 );

  formula = "3<<4";
  correct &= dt.testFormula( formula, vars, &child, 48 );

  formula = "238>>3";
  correct &= dt.testFormula( formula, vars, &child, 29 );

  formula = "abs((p2-p1-p0))";
  correct &= dt.testFormula( formula, vars, &child, 10 );

  formula = "!0";
  correct &= dt.testFormula( formula, vars, &child, 1 );

  formula = "~10&15";
  correct &= dt.testFormula( formula, vars, &child, 5 );

  formula = "-(3-12)";
  correct &= dt.testFormula( formula, vars, &child, 9 );

  formula = "-3--12-4*-2";
  correct &= dt.testFormula( formula, vars, &child, 17 );

  ////////////////////////////////////////////////////////////////////// 
  // rotate(val,num,pos,bitwidth)
  ////////////////////////////////////////////////////////////////////// 
  formula = "rotate(1234, 0, 3, 7)";
  correct &= dt.testFormula( formula, vars, &child, 1234 );

  formula = "rotate(255, 3, 0, 16)";
  correct &= dt.testFormula( formula, vars, &child, 57375 );

  formula = "rotate(255, 3, 3, 12)";
  correct &= dt.testFormula( formula, vars, &child, 28703 );

  formula = "rotate(204, 0-3, 2, 6)";
  correct &= dt.testFormula( formula, vars, &child, 120 );

  formula = "rotate(4044, 0-3, 2, 6)";
  correct &= dt.testFormula( formula, vars, &child, 3960 );

  formula = "rotate(255, 0-3, 0, 16)";
  correct &= dt.testFormula( formula, vars, &child, 2040 );

  formula = "rotate(255, 0-3, 2, 16)";
  correct &= dt.testFormula( formula, vars, &child, 2019 );

  ////////////////////////////////////////////////////////////////////// 
  // median 
  ////////////////////////////////////////////////////////////////////// 
  formula = "median(8,12,37)";
  correct &= dt.testFormula( formula, vars, &child, 12 );

  formula = "median(7,4,4,7,4,7)";
  correct &= dt.testFormula( formula, vars, &child, 6 );

  formula = "median(4,4,4,4,4,4)";
  correct &= dt.testFormula( formula, vars, &child, 4 );

  formula = "median(4,4,4,2,4)";
  correct &= dt.testFormula( formula, vars, &child, 4 );

  formula = "median(p0,max(p1*8-12,0)+3, 2*p2)";
  correct &= dt.testFormula( formula, vars, &child, 7 );

  for (int nmbr_tests=0; nmbr_tests<10; nmbr_tests++) {
    vector<int> some_numbers;
    int nmbr_operands;
    int random;
    int answer;
    char number_string[128]; 

    some_numbers.clear();
    nmbr_operands = rnd.uniform(3, 20);
    formula = "median(";
    for (int i=0; i<nmbr_operands; i++) {
      if (i != 0) {
        formula = formula + ",";
      }
      random = rnd.uniform(0, 15);
      some_numbers.push_back(random);
      sprintf(number_string, "%d", random);
      formula = formula + number_string;
    }
    formula = formula + ")"; 

    sort(some_numbers.begin(), some_numbers.end());
    
    if (nmbr_operands & 1) {
      answer = some_numbers[nmbr_operands/2];
    }
    else {
      answer = (some_numbers[nmbr_operands/2] + 
                some_numbers[nmbr_operands/2-1] + 1) / 2;
    }

    correct &= dt.testFormula( formula, vars, &child, answer );
  }

  formula = "all_diff( 0, 0 )";
  correct &= dt.testFormula( formula, vars, &child, 0 );

  formula = "all_diff( 0, 1 )";
  correct &= dt.testFormula( formula, vars, &child, 1 );

  formula = "all_diff( 7, 3+4, p0-1, 7, 7, 7, 7 ) + 4";
  correct &= dt.testFormula( formula, vars, &child, 4 );

  formula = "4+all_diff( 7, 3+5, p0-2, 9, 10, 11, 12, 13, 14, 14 )";
  correct &= dt.testFormula( formula, vars, &child, 4 );

  formula = "all_diff( 7, 3+5, p0-2, 9, 10, 11, 12, 13, 14, 15, 16 )";
  correct &= dt.testFormula( formula, vars, &child, 1 );

  ////////////////////////////////////////////////////////////////////// 
  // Finish 
  ////////////////////////////////////////////////////////////////////// 
  if (correct) {
    printf("All formulas gave the correct answer\n"); 
  }
  else {
    error("A formula is not giving the correct answer\n");
  }
}

////////////////////////////////////////////////////////////////////
// Test decision trees 
////////////////////////////////////////////////////////////////////
/*
 * A decision tree looks something like:
 * <numFtrs>
 * <ftr> <num_splits> r1 r2 ... rs 
 * <ftr> <num_splits> r1 r2 ... rs # split 1
 * <ftr> <num_splits> r1 r2 ... rs # split 2
 * -1 value
 * <ftr> <num_splits> r1 r2 ... rs # split s
 */


// Some example DTs to write to a file, read back in,
// and then test out.

char *dtStr1 =
"% this is a decision tree file\n"
"%\n"
"dt_name 3  % number of features\n"
//    0   1  2  3  4   5  6   7  8   9
"0 10 10 11 12 13 0:5 6:9 14,43 50:10:150 55:10:155 default\n"
"  1 2 0:10 default\n"
"    2 2 0:10 default\n"
"      -1 42\n"
"      -1 {  (  p0  \n  +   1    )    }\n"
"    2 2 0:5 default\n"
"      -1 {(cp0+1)}\n"
"      -1 {(00+1)}\n"
"  1 2 0:10 default\n"
"    2 2 0:10 default\n"
"      -1 {(p0+p1+5)}\n"
"      -1 6\n"
"    2 2 0:5 default\n"
"      -1 7\n"
"      -1 8\n"
"  -1 10 % when feature[0] = 12, map to 10 regardless of all else\n"
"  -1 11 % when feature[0] = 13, map to 11 regardless of all else\n"
"  -1 12 % when feature[0] = 0:5, map to 12 regardless of all else\n"
"  -1 13 % when feature[0] = 6:9, map to 13 regardless of all else\n"
"  -1 14 % when feature[0] = 4, map to 14 regardless of all else\n"
"  -1 15 % when feature[0] >=50, map to 15 regardless of all else\n"
"  -1 16 % when feature[0] = 15, map to 16 regardless of all else\n"
"  1 3 0:10 11 default  % none of the above, formula again\n"   
"    2 7 0 ... 5 default\n"
"      -1 99\n"
"      -1 101\n"
"      -1 202\n"
"      -1 303\n"
"      -1 404\n"
"      -1 505\n"
"      -1 606\n"
"    2 10 100 5 3 4 9 303 404 40 50 default\n"
"      -1 11\n"
"      -1 12\n"
"      -1 13\n"
"      -1 14\n"
"      -1 15\n"
"      -1 16\n"
"      -1 17\n"
"      -1 18\n"
"      -1 19\n"
"      -1 20\n"
"    2 10 100 5 3 4 9 303 404 40 51 default\n"
"      -1 1100\n"
"      -1 1200\n"
"      -1 1300\n"
"      -1 1400\n"
"      -1 1500\n"
"      -1 1600\n"
"      -1 1700\n"
"      -1 1800\n"
"      -1 1900\n"
"      -1 2000\n"
"\n\n";



char *dtStr2 =
"% this is a decision tree file\n"
"%\n"
"dt_name 1  % number of features\n"
"0 14 10 17,19,21 11 12 13 0:5 6:9 14 50: 15 18,20 23,25,27 22,24,26,28 default\n"
// "0 14 0 ... 12 default\n"
"  -1 2\n"
"  -1 3\n"
"  -1 4\n"
"  -1 5\n"
"  -1 0\n"
"  -1 1\n"
"  -1 6\n"
"  -1 50\n"
"  -1 7\n"
"  -1 8\n"
"  -1 9\n"
"  -1 10\n"
"  -1 11\n"
"  0 16 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 default\n"
"     -1 0\n"
"     -1 1\n"
"     -1 2\n"
"     -1 3\n"
"     -1 4\n"
"     -1 5\n"
"     -1 6\n"
"     -1 7\n"
"     -1 8\n"
"     -1 9\n"
"     -1 10\n"
"     -1 11\n"
"     -1 12\n"
"     -1 13\n"
"     -1 14\n"
"     -1 15\n"
;


char *dtStr3 = 
"wordPosition % name \n"
"3 % parents: wordTransition(-1), phoneTransition(-1), wordPosition(-1)\n"
"1 2 0 default\n"
"  -1 {(p2)}\n"
"   0 2 1 default\n"
"      -1 0\n"
"       2 2 7 default\n"
"           -1 7\n"
"	    -1 {(p2+1)}\n"
"\n\n";


char *dtStr4 = 
"eou13A_0 % DT name\n"
"2 % number of parents\n"
"0 2 1 default\n"
"         1 3 4 3 default\n"
"                 -1 1\n"
"                 -1 0\n"
"                 -1 0\n"
"\n\n";

char *dtStr5 = 
"foobar % DT name\n"
"3 % number of parents\n"
"0 2 0:4,8,10 default\n"
"  -1 1\n"
"  -1 0\n"
"\n\n";



void test_dts()
{

  {
    oDataStreamFile dtfile ("/tmp/foo.dt",false);
    dtfile.write(dtStr5);
    dtfile.write('\n',"",false);
  }
  iDataStreamFile is (file,false);
  RngDecisionTree dt;
  dt.read(is);




#ifdef REMOVE_THIS_CODE_FOR_NOW 

  // first write out the file
  
  printf("Writing to /tmp/foo.dt\n");



  {
    char *file = "/tmp/foo.dt";
    iDataStreamFile is (file,false);
    RngDecisionTree dt;
    dt.read(is);
    oDataStreamFile dtfile ("/tmp/foo_output.dt",false);
    dt.write(dtfile);
    dtfile.write('\n',"",false);
  }


  // a test for memory leaks 
  while (1) {
    char *file = "/tmp/foo.dt";
    iDataStreamFile is (file,false);
    RngDecisionTree dt;
    dt.read(is);
  }

  char *file = "/tmp/foo.dt";
  iDataStreamFile is (file,false);
  RngDecisionTree dt;
  dt.read(is);


  printf("Found decision tree\n");

  RVInfo::FeatureRange tmp_fr;
  RVInfo::ListIndex    tmp_li;
  vector< RVInfo::rvParent > tmp_switchingParents;
  vector<vector< RVInfo::rvParent > > tmp_conditionalParents;
  vector< CPT::DiscreteImplementaton > tmp_discImplementations;
  vector< MixtureCommon::ContinuousImplementation > tmp_contImplementations;
  vector< RVInfo::ListIndex > tmp_listIndices;
  vector< RVInfo::WeightInfo > tmp_rvWeightInfo;

  RVInfo dummy(
    0.0, 1, 0, 10, "test", "aa_info", RVInfo::t_discrete, RVInfo::d_hidden, 100, 
    tmp_fr, NULL, tmp_li, tmp_switchingParents, 
    tmp_conditionalParents, tmp_discImplementations, tmp_contImplementations,
    tmp_listIndices, tmp_rvWeightInfo 
    );

  TestRandomVariable p0(dummy, "p0",  8 );
  TestRandomVariable p1(dummy, "p1", 40 );
  TestRandomVariable p2(dummy, "p2",  2 );
  TestRandomVariable child(dummy, "child",  1000000 );

  p0.val  = 7;
  p1.val  = 3;
  p2.val  = 0;

  vector<RV*> vars;
  vars.push_back(&p0);
  vars.push_back(&p1);
  vars.push_back(&p2);

  assert ( dt.numFeatures() == 3 );

  iDataStreamFile stin ("-",false,false);

  vector<int> vec;
  vec.resize(dt.numFeatures());

  while (1) {
    printf("Enter a length %d intvec:",dt.numFeatures());
    fflush(stdout);
    stin.read(vec,dt.numFeatures());

    p0.val = vec[0];
    p1.val = vec[1];
    p2.val = vec[2];

    printf("Querying with vector and cards: ");
    fflush(stdout);
    printf(" %d:%d",p0.val,p0.cardinality);
    printf(" %d:%d",p1.val,p1.cardinality);
    printf(" %d:%d",p2.val,p2.cardinality);
    printf(" chd_crd = %d",child.cardinality);
    printf("\n");
    printf("### RESULT ==> %d\n",dt.query(vars,(RV*)&child));
  }


#endif 
}

int
main(int argc,char *argv[])
{


  printf("sizeof RngDecisionTree::Node = %d\n",sizeof(RngDecisionTree::Node));
  /*
  printf("sizeof EquationClass = %d\n",sizeof(RngDecisionTree::EquationClass));
  printf("sizeof BP_Range = %d\n",sizeof(BP_Range));

  printf("sizeof NonLeafNodeArrayStruct = %d\n",sizeof(RngDecisionTree::NonLeafNodeArrayStruct));
  printf("sizeof NonLeafNodeHashStruct = %d\n",sizeof(RngDecisionTree::NonLeafNodeHashStruct));
  printf("sizeof NonLeafNodeRngsStruct = %d\n",sizeof(RngDecisionTree::NonLeafNodeRngsStruct));
  printf("sizeof LeafNodeValStruct = %d\n",sizeof(RngDecisionTree::LeafNodeValStruct));
  printf("sizeof LeafNodeEquationStruct = %d\n",sizeof(RngDecisionTree::LeafNodeEquationStruct));
*/

  test_dts();

  // Test the formula parser 
  // test_formula();



}

#endif

