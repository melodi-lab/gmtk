/*-
 * GMTK_RngDecisionTree.cc
 *     General Class to map from vectors of integers to some basic type.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
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
#include "sArray.h"

#include "GMTK_RandomVariable.h"
#include "GMTK_RngDecisionTree.h"

VCID("$Header$");

/////////////////////////////////////////////////////////////////////
// File extension for compiled DT files 
/////////////////////////////////////////////////////////////////////
const string DTFileExtension = ".index";

/////////////////////////////////////////////////////////////////////
// Static computation stack for formula evaluation 
/////////////////////////////////////////////////////////////////////
RngDecisionTree::sArrayStack<int> RngDecisionTree::EquationClass::stack(0);

/////////////////////////////////////////////////////////////////////
// Arrays to classify formula tokens and map their strings to their 
// enumerations 
/////////////////////////////////////////////////////////////////////

map<string, RngDecisionTree::EquationClass::tokenEnum> 
  RngDecisionTree::EquationClass::delimiter;
map<string, RngDecisionTree::EquationClass::tokenEnum> 
  RngDecisionTree::EquationClass::function;

map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::termToken;
map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::expressionToken;
map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::functionToken;
map<RngDecisionTree::EquationClass::tokenEnum, 
  RngDecisionTree::EquationClass::formulaCommand> 
  RngDecisionTree::EquationClass::twoValFunctionToken;
			  
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
  destructorRecurse(root);
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
  if (node->nodeType != NonLeafNode) {
    // do nothing since this is a leaf
  } else {
    for (unsigned i=0;i<node->nonLeafNode.children.size();i++) {
      destructorRecurse(node->nonLeafNode.children[i]);
      delete node->nonLeafNode.children[i];
      if (i<node->nonLeafNode.children.size()-1)
	delete node->nonLeafNode.rngs[i];
    }
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
 *      in and clampNextDecisionTree() must be called for
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
  is.read(dtFileName,"RngDecisionTree:: read file/numFeatures");

  //////////////////////////////////////////////////////////////////////
  // Entry is the name of a file containing multiple instances of the 
  // decision tree.
  //////////////////////////////////////////////////////////////////////
  if (!strIsInt(dtFileName.c_str(), (int*)&_numFeatures)) {

    //////////////////////////////////////////////////////////////////////
    // Make sure this instance already hasn't been initialized as an
    // iterable decision tree 
    //////////////////////////////////////////////////////////////////////
    if (clampable()) {
      error("ERROR: in DT named '%s' in file '%s', can't have DTs defined recursively in files",
        name().c_str(),is.fileName());
    }

    initializeIterableDT(dtFileName); 
  }
  //////////////////////////////////////////////////////////////////////
  // Entry is an integer indicating an inline decision tree 
  //////////////////////////////////////////////////////////////////////
  else {

    if (_numFeatures < 0) {
      error("ERROR: in DT named '%s', file '%s', decision tree must have >= 0 features",
      name().c_str(), is.fileName());
    }

    rightMostLeaf = NULL;
    root = readRecurse(is,rightMostLeaf);

    // set up the 'right' pointers in the doubly 
    // linked lists.
    Node *tmp = rightMostLeaf;
    Node *rightLeaf = NULL;
    while (tmp != NULL) {
      tmp->leafNode.nextLeaf = rightLeaf;
      rightLeaf = tmp;
      tmp = tmp->leafNode.prevLeaf;
    }
    leftMostLeaf = rightLeaf;
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
  if (clampable()) {
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
  clampFirstDecisionTree();
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

  if (!clampable()) {
    error("ERROR: trying to seek in non-clampable DT, '%s'\n",  
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
RngDecisionTree::Node* 
RngDecisionTree::readRecurse(iDataStreamFile& is,
			     Node* &prevLeaf)
{

  int curFeat;
  is.read(curFeat,"RngDecisionTree:: readRecurse curFeat");
  if (curFeat < -1) 
    error("ERROR: DT '%s', file '%s', feature number (=%d) must be non-negative",name().c_str(),is.fileName(),curFeat);
  if (curFeat >= (int)_numFeatures) 
    error("ERROR: DT '%s', file '%s', feature number (=%d) must be < numFeatures (=%d)",name().c_str(),is.fileName(),curFeat,
	  _numFeatures);
  Node *node = new Node;
  if (curFeat == -1) {
    // then this is a leaf node.
    node->leafNode.value = 0;
    // leaf node, now read the next bit of leaf node as a string to determine what it is.
    string leafNodeVal;
    is.read(leafNodeVal,"RngDecisionTree:: readRecurse leaf node value");
    int val;
    if (strIsInt(leafNodeVal.c_str(),&val)) {
      node->nodeType = LeafNodeVal;
      node->leafNode.value = val;
      node->leafNode.leafNodeString = leafNodeVal;
    } else if (leafNodeVal == "expand") {
      node->nodeType = LeafNodeFullExpand;
      node->leafNode.leafNodeString = leafNodeVal;
    } else if (leafNodeVal[0] == '(') {
      node->nodeType = LeafNodeFormula;
      // TODO: remove exceptions
      try {
        node->leafNode.equation.parseFormula(leafNodeVal);
      }
      catch( string error_message ){
        error("ERROR: In file '%s', DT '%s', equation '%s':  %s", 
          is.fileName(), name().c_str(), leafNodeVal.c_str(), 
	  error_message.c_str());
      }
      catch( const char * const error_message ) {
        error("ERROR: In file '%s', DT '%s', equation '%s':  %s", 
          is.fileName(), name().c_str(), leafNodeVal.c_str(), 
	  error_message );
      }
    } else {
      error("ERROR: In file '%s', DT '%s', invalid leaf node value '%s':  %s", 
	    is.fileName(), name().c_str(), leafNodeVal.c_str() ); 
    }

    node->leafNode.prevLeaf = prevLeaf;
    prevLeaf = node;
  } else {
    node->nodeType = NonLeafNode;
    node->nonLeafNode.ftr = (leafNodeValType)curFeat;
    node->nonLeafNode.ordered = true;
    unsigned numSplits;
    is.read(numSplits,"RngDecisionTree:: readRecurse numSplits");
    if (numSplits < 1)
      error("ERROR: DT '%s', file '%s', can't have < 1 node splits",name().c_str(),is.fileName());
    if (numSplits > RNG_DECISION_TREE_MAX_ARY)
      error("ERROR: DT '%s', file '%s': can't have > %d splits",
	    name().c_str(),is.fileName(),RNG_DECISION_TREE_MAX_ARY);
    node->nonLeafNode.children.resize(numSplits);
    // rngs is smaller since last string is always the default catch-all.
    node->nonLeafNode.rngs.resize(numSplits-1);
    for (unsigned i=0;i<numSplits;i++) {
      char *str;
      is.read(str,"RngDecisionTree:: readRecurse, reading range");
      if (i==(numSplits-1)) {
	if (strcmp(RNG_DECISION_TREE_DEF_STR,str))
	  error("ERROR: DT '%s', file '%s': expecting default str (%s) got (%s)",
		name().c_str(),is.fileName(),RNG_DECISION_TREE_DEF_STR,str);
      } else {
	// note: ideally, we would limit the maximum range
	// value to be equal to the cardinality of the random
	// variable here. We can't do that, however, because
	// the DT is generic, and could be used with multiple
	// different RVs with different cardinalities.
	node->nonLeafNode.rngs[i]
	  = new BP_Range(str,
			 0,
			 MAX_BP_RANGE_VALUE);
      }
      ///////////////////////////////////////////////////////////////
      // WARNING: We assume here that BP_Range will make
      // its own copy of the string, so we delete it here.
      // If the above class changes, we might need to change this code.
      delete [] str;
      ///////////////////////////////////////////////////////////////
    }

    //////////////////////////////////////////////
    // check for overlap and order errors in the strings
    for (unsigned i=0;i<numSplits-1;i++) {
      for (unsigned j=i+1;j<numSplits-1;j++) {
	if (node->nonLeafNode.rngs[i]
	    ->overlapP(node->nonLeafNode.rngs[j]))
	  error("ERROR: DT '%s', file '%s': range %d (%s) and %d (%s) have a non-empty intersection.",name().c_str(),is.fileName(),
		i,
		node->nonLeafNode.rngs[i]->rangeStr(),j,
		node->nonLeafNode.rngs[j]->rangeStr());
	
	////////////////////////////////////////////////////////
	// this next check is required to ensure there is
	// an ordering of the ranges, so sorting makes sense
	// when we do binary search. An unordered list of ranges would, for
	// example, have things like 1,3,5 and 2,4 where you can't say
	// that either (1,3,5) < (2,4) nor (2,4) < (1,3,5). In these
	// cases, we do not sort the lists, and have to resort to linear search
	// when doing a query.
	if ( 
	    (!((*node->nonLeafNode.rngs[i]) < (*node->nonLeafNode.rngs[j])))
	    &&
	    (!((*node->nonLeafNode.rngs[j]) < (*node->nonLeafNode.rngs[i])))
	    ) {
	  // then it isn't ordered, and we'll have to do linear search
	  // on a query. 
	  // 
	  // TODO: check if it is ordered alread, and if
	  // so, don't do the sort below.
	  node->nonLeafNode.ordered = false;
	}
      }
    }

    for (unsigned i=0; i<numSplits; i++)
      node->nonLeafNode.children[i] =
	readRecurse(is,prevLeaf);

    if (numSplits >= DT_SPLIT_SORT_THRESHOLD && node->nonLeafNode.ordered) {
      //////////////////////////////////////////////////////////
      // sort the entries so we can do binsearch later. Sort only
      // if 1) there are a sufficient number to warrant a sort and 2)
      // the list of ranges is sortable.

      // copy in
      vector< pair<BP_Range*,Node*> > arr;
      arr.resize(numSplits-1);
      for (unsigned i=0;i<(numSplits-1);i++) {
	arr[i].first = node->nonLeafNode.rngs[i];
	arr[i].second = node->nonLeafNode.children[i];
      }
      //'sort
      sort(arr.begin(),
	   arr.end(),
	   RngCompare());
      // copy out
      for (unsigned i=0;i<(numSplits-1);i++) {
	node->nonLeafNode.rngs[i] = arr[i].first;
	node->nonLeafNode.children[i] = arr[i].second;
      }
    }

  }
  return node;
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

    delimiter["&"]  = TOKEN_BITWISE_AND;
    delimiter["|"]  = TOKEN_BITWISE_OR;
    delimiter[":"]  = TOKEN_COLON;
    delimiter[","]  = TOKEN_COMMA;
    delimiter["/"]  = TOKEN_DIVIDE;
    delimiter["=="] = TOKEN_EQUALS;
    delimiter["^"]  = TOKEN_EXPONENT;
    delimiter[">"]  = TOKEN_GREATER_THAN;
    delimiter[">="] = TOKEN_GREATER_THAN_EQ;
    delimiter["<"]  = TOKEN_LESS_THAN;
    delimiter["<="] = TOKEN_LESS_THAN_EQ;
    delimiter["&&"] = TOKEN_LOGICAL_AND;
    delimiter["||"] = TOKEN_LOGICAL_OR;
    delimiter["-"]  = TOKEN_MINUS;
    delimiter["+"]  = TOKEN_PLUS;
    delimiter[" "]  = TOKEN_SPACE;
    delimiter["\t"] = TOKEN_SPACE;
    delimiter["("]  = TOKEN_LEFT_PAREN;
    delimiter["*"]  = TOKEN_TIMES;
    delimiter[")"]  = TOKEN_RIGHT_PAREN;
    delimiter["?"]  = TOKEN_QUESTION_MARK;

    function["max"] = TOKEN_MAX;
    function["min"] = TOKEN_MIN;
    function["mod"] = TOKEN_MOD;
    function["xor"] = TOKEN_BITWISE_XOR;

    expressionToken[TOKEN_MINUS] = COMMAND_MINUS;
    expressionToken[TOKEN_PLUS]  = COMMAND_PLUS;

    termToken[TOKEN_BITWISE_AND]     = COMMAND_BITWISE_AND;
    termToken[TOKEN_BITWISE_OR]      = COMMAND_BITWISE_OR;
    termToken[TOKEN_BITWISE_XOR]     = COMMAND_BITWISE_XOR;
    termToken[TOKEN_DIVIDE]          = COMMAND_DIVIDE;
    termToken[TOKEN_EXPONENT]        = COMMAND_EXPONENT;
    termToken[TOKEN_EQUALS]          = COMMAND_EQUALS;
    termToken[TOKEN_GREATER_THAN]    = COMMAND_GREATER_THAN;
    termToken[TOKEN_GREATER_THAN_EQ] = COMMAND_GREATER_THAN_EQ;
    termToken[TOKEN_LESS_THAN]       = COMMAND_LESS_THAN;
    termToken[TOKEN_LESS_THAN_EQ]    = COMMAND_LESS_THAN_EQ;
    termToken[TOKEN_LOGICAL_AND]     = COMMAND_LOGICAL_AND;
    termToken[TOKEN_LOGICAL_OR]      = COMMAND_LOGICAL_OR; 
    termToken[TOKEN_TIMES]           = COMMAND_TIMES;
    termToken[TOKEN_QUESTION_MARK]   = COMMAND_QUESTION_MARK;

    functionToken[TOKEN_MAX] = COMMAND_MAX;
    functionToken[TOKEN_MIN] = COMMAND_MIN;

    twoValFunctionToken[TOKEN_BITWISE_XOR] = COMMAND_BITWISE_XOR;
    twoValFunctionToken[TOKEN_MOD]         = COMMAND_MOD; 
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
  const vector< RandomVariable* >& variables 
  )
{
  unsigned crrnt_cmnd, end_cmnd;
  unsigned last; 
  unsigned command, operand; 
  int      value, i;

  stack.clear();

  for (crrnt_cmnd = 0, 
       end_cmnd   = commands.size();
       crrnt_cmnd != end_cmnd;
       ++crrnt_cmnd) {

    command = GET_COMMAND(commands[crrnt_cmnd]); 

    switch (command) {
	
      case COMMAND_PUSH_PARENT:	
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        if (operand >= variables.size()) {	
          error("ERROR:  Reference to non-existant parent\n"); 	
        }
        stack.push_back( variables[operand]->val );	
        break;	

      case COMMAND_PUSH_CARDINALITY:
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        if (operand >= variables.size()) {
          error("ERROR:  Reference to non-existant parent cardinality\n"); 
        }
        stack.push_back( variables[operand]->cardinality ); 
        break;

      case COMMAND_PUSH_CONSTANT:
        operand = GET_OPERAND(commands[crrnt_cmnd]); 
        stack.push_back(operand); 
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
 
      case COMMAND_BITWISE_XOR: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] ^ stack[last];
        stack.pop_back();
        break;
 
      case COMMAND_DIVIDE:
        last = stack.stackSize() - 1;
        if (stack[last] == 0) {
          error("ERROR:  Divide by zero error\n"); 
        }
        stack[last-1] = stack[last-1] / stack[last];
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

      case COMMAND_PLUS:
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] + stack[last];
        stack.pop_back();
        break;

      case COMMAND_TIMES: 
        last = stack.stackSize() - 1;
        stack[last-1] = stack[last-1] * stack[last];
        stack.pop_back();
        break;

      case COMMAND_QUESTION_MARK: 
        last = stack.stackSize() - 1;
        stack[last-2] = stack[last-2] ? stack[last-1] : stack[last]; 
        stack.pop_back();
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
 * RngDecisionTree::EquationClass::evaluateFormula
 *   Stub function to support some leftover tests
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
  const vector<unsigned>& value,
  const vector<unsigned>& cardinality
  )
{
  assert(0);
  return(0);
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
 *   The decision tree's equation data structure is filled in. 
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
  string leafNodeVal
  )
{
  tokenStruct token;
  unsigned    depth;

  assert( (LAST_COMMAND_INDEX & COMMAND_MASK) == LAST_COMMAND_INDEX );
 
  depth = 0; 
  do {
    getToken(leafNodeVal, token); 
    parseExpression(token, leafNodeVal, depth ); 
  } while (token.token != TOKEN_END);

}


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::parseExpression
 *   Support function for parseFormula.  An expression includes terms
 *   separated by + and -. 
 * 
 * Preconditions:
 *   none 
 *
 * Postconditions:
 *   The decision tree's equation data structure is partially filled in. 
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
  tokenStruct& token,
  string&      leafNodeVal,
  unsigned&    depth  
  )
{
  formulaCommand new_command;
  tokenEnum      next_token;

  parseTerm(token, leafNodeVal, depth); 

  while( expressionToken[token.token] != COMMAND_INVALID ) { 

    next_token  = token.token;
    getToken(leafNodeVal, token); 
    parseTerm( token, leafNodeVal, depth );

    new_command = MAKE_COMMAND( expressionToken[next_token], 0 );
    commands.resizeAndCopy( commands.size()+1 );
    commands[commands.size()-1] = new_command;
    changeDepth( -1, depth );
  } 
}


/*-
 *-----------------------------------------------------------------------
 * RngDecisionTree::EquationClass::parseTerm
 *   Support function for parseFormula.  An term contains factors 
 *   separated by operators that are evaluated in the order listed. 
 * 
 * Preconditions:
 *   none 
 *
 * Postconditions:
 *   The decision tree's equation data structure is partially filled in. 
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   none 
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::EquationClass::parseTerm(
  tokenStruct& token,
  string&      leafNodeVal,  
  unsigned&    depth  
  )
{
  formulaCommand new_command;
  tokenEnum      next_token;
 
  parseFactor(token, leafNodeVal, depth);

  while( termToken[token.token] != COMMAND_INVALID ) { 
 
    switch (token.token) {

      case TOKEN_QUESTION_MARK:  
        getToken(leafNodeVal, token); 
        parseExpression(token, leafNodeVal, depth);

        if (token.token != TOKEN_COLON) {
          string error_message = "Expecting colon after '?' operator"; 
          throw(error_message); 
        }

        getToken(leafNodeVal, token); 
        parseExpression(token, leafNodeVal, depth);

        new_command = MAKE_COMMAND( COMMAND_QUESTION_MARK, 0 );
        commands.resizeAndCopy( commands.size()+1 );
        commands[commands.size()-1] = new_command;
        changeDepth( -1, depth );
        break;

      default:
        next_token  = token.token;
        getToken(leafNodeVal, token); 
        parseFactor(token, leafNodeVal, depth);

        new_command = MAKE_COMMAND( termToken[next_token], 0 );
        commands.resizeAndCopy( commands.size()+1 );
        commands[commands.size()-1] = new_command;
        changeDepth( -1, depth );
        break;
    }
  }
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
 *   The decision tree's equation data structure is partially filled in. 
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
RngDecisionTree::EquationClass::parseFactor(
  tokenStruct& token,
  string&      leafNodeVal,  
  unsigned&    depth  
  )
{
  formulaCommand new_command;
  tokenEnum      next_token;

  switch (token.token) {

    case TOKEN_LEFT_PAREN:  
      getToken(leafNodeVal, token); 
      parseExpression(token, leafNodeVal, depth);

      if (token.token == TOKEN_RIGHT_PAREN) {

        getToken(leafNodeVal, token); 
      }
      else {
        string error_message = "Expecting right parenthesis at '" + leafNodeVal 
          + "'";
        throw(error_message); 
      }
      break;
 
    case TOKEN_INTEGER:
      new_command = MAKE_COMMAND( COMMAND_PUSH_CONSTANT, token.number ); 
      commands.resizeAndCopy( commands.size()+1 );
      commands[commands.size()-1] = new_command;
      changeDepth( 1, depth );
      getToken(leafNodeVal, token); 
      break;
 
    case TOKEN_PARENT:
      new_command = MAKE_COMMAND( COMMAND_PUSH_PARENT, token.number );
      commands.resizeAndCopy( commands.size()+1 );
      commands[commands.size()-1] = new_command;
      changeDepth( 1, depth );
      getToken(leafNodeVal, token); 
      break;
 
    case TOKEN_CARDINALITY:
      new_command = MAKE_COMMAND( COMMAND_PUSH_CARDINALITY, token.number );
      commands.resizeAndCopy( commands.size()+1 );
      commands[commands.size()-1] = new_command;
      changeDepth( 1, depth );
      getToken(leafNodeVal, token); 
      break;

    default:
      //////////////////////////////////////////////////////////////////////
      // Functions which take two or more operands 
      //////////////////////////////////////////////////////////////////////
      if (functionToken[token.token] != COMMAND_INVALID) {
  
        next_token = token.token;

        getToken(leafNodeVal, token); 
        if (token.token != TOKEN_LEFT_PAREN) {
          string error_message = "Expecting left parenthesis at '" + 
	    leafNodeVal + "'";
          throw(error_message);
        }

        getToken(leafNodeVal, token); 
        parseExpression(token, leafNodeVal, depth);
        if (token.token != TOKEN_COMMA) {
          string error_message = "Function requires at least two operands"; 
          throw(error_message);
        }

        getToken(leafNodeVal, token);
        parseExpression(token, leafNodeVal, depth);

        while (token.token == TOKEN_COMMA) {
          new_command = MAKE_COMMAND( functionToken[next_token], 0 );
          commands.resizeAndCopy( commands.size()+1 );
          commands[commands.size()-1] = new_command;
          changeDepth( -1, depth );

          getToken(leafNodeVal, token);
          parseExpression(token, leafNodeVal, depth);
        }

        if (token.token != TOKEN_RIGHT_PAREN) {
          string error_message = "Expecting right parenthesis at '" + 
            leafNodeVal + "'";
          throw(error_message);
        }

        new_command = MAKE_COMMAND( functionToken[next_token], 0 );
        commands.resizeAndCopy( commands.size()+1 );
        commands[commands.size()-1] = new_command;
        changeDepth( -1, depth );

        getToken(leafNodeVal, token); 
      }
      //////////////////////////////////////////////////////////////////////
      // Functions which take exactly two operands 
      //////////////////////////////////////////////////////////////////////
      else if (twoValFunctionToken[token.token] != COMMAND_INVALID) {
  
        next_token = token.token;

        getToken(leafNodeVal, token); 
        if (token.token != TOKEN_LEFT_PAREN) {
          string error_message = "Expecting left parenthesis at '" + 
	    leafNodeVal + "'";
          throw(error_message);
        }

        getToken(leafNodeVal, token); 
        parseExpression(token, leafNodeVal, depth);
        if (token.token != TOKEN_COMMA) {
          string error_message = "Function requires exactly two operands"; 
          throw(error_message);
        }

        getToken(leafNodeVal, token);
        parseExpression(token, leafNodeVal, depth);

        if (token.token != TOKEN_RIGHT_PAREN) {
          string error_message = "Expecting right parenthesis at '" + 
            leafNodeVal + "'";
          throw(error_message);
        }

        new_command = MAKE_COMMAND( twoValFunctionToken[next_token], 0 );
        commands.resizeAndCopy( commands.size()+1 );
        commands[commands.size()-1] = new_command;
        changeDepth( -1, depth );
        getToken(leafNodeVal, token); 
      }

      break;
  }
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
  unsigned index;

  map<string, tokenEnum>::iterator crrnt_dlmtr; 
  map<string, tokenEnum>::iterator end_dlmtr; 
  size_t dlmtr_lctn; 
  size_t minimum_dlmtr_lctn;
 
  vector<map<string, tokenEnum>::iterator> found_dlmtr_cntnr; 
  vector<map<string, tokenEnum>::iterator>::iterator crrnt_found_dlmtr; 
  vector<map<string, tokenEnum>::iterator>::iterator end_found_dlmtr; 
  map<string, tokenEnum>::iterator found_dlmtr; 

  //////////////////////////////////////////////////////////////////////////
  // Skip leading spaces 
  //////////////////////////////////////////////////////////////////////////
  index = 0;
  while ((index < expression.length()) &&   
         ((expression[index] == ' ') || (expression[index] == '\t'))) {
    ++index;
  }

  //////////////////////////////////////////////////////////////////////////
  // Find the earliest delimiter token 
  //////////////////////////////////////////////////////////////////////////
  token.token = LAST_TOKEN_INDEX;
  minimum_dlmtr_lctn = expression.size();

  for( crrnt_dlmtr = delimiter.begin(),
       end_dlmtr   = delimiter.end();
       crrnt_dlmtr != end_dlmtr;  
       ++crrnt_dlmtr ) { 

    dlmtr_lctn = expression.find( (*crrnt_dlmtr).first, index );    
    if (dlmtr_lctn >= 0) {
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
  // If no tokens were found, return the end token 
  //////////////////////////////////////////////////////////////////////////
  if (found_dlmtr_cntnr.size() == 0) {
    token.token = TOKEN_END;
  }
  //////////////////////////////////////////////////////////////////////////
  // If next token is a delimiter, determine the longest delimiter which 
  // matches. 
  //////////////////////////////////////////////////////////////////////////
  else if (minimum_dlmtr_lctn == 0) {

    found_dlmtr = *found_dlmtr_cntnr.begin();
    for( crrnt_found_dlmtr = found_dlmtr_cntnr.begin(),
         end_found_dlmtr  = found_dlmtr_cntnr.end();
         crrnt_found_dlmtr != end_found_dlmtr;  
         ++crrnt_found_dlmtr) { 

      if ( (*(*crrnt_found_dlmtr)).first.size() > (*found_dlmtr).first.size() )       {
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
    string rest;
    int    number;    
    bool   is_integer;

    is_integer = getInteger( token_string, number ); 
    if (is_integer) {
      token.token  = TOKEN_INTEGER;
      token.number = number; 
    } 
    //////////////////////////////////////////////////////////////////////////
    // Next, check for parent or cardinality (such as "p0" or "c2")
    //////////////////////////////////////////////////////////////////////////
    else if ( token_string[0] == 'p' ) {
  
      rest = token_string.substr(1, token_string.length());
      is_integer = getInteger( rest,  number ); 
      if (is_integer) {
        token.token  = TOKEN_PARENT;
        token.number = number; 
      } 
    }
    else if ( token_string[0] == 'c' ) {
  
      rest = token_string.substr(1, token_string.length());
      is_integer = getInteger( rest,  number ); 
      if (is_integer) {
        token.token  = TOKEN_CARDINALITY;
        token.number = number; 
      } 
    }

    //////////////////////////////////////////////////////////////////////////
    // Check if token string is a multiple character key word 
    //////////////////////////////////////////////////////////////////////////
    else {

      map<string, tokenEnum>::iterator found_string;
      found_string = function.find( token_string );
      if (found_string != function.end()) {
        token.token = (*found_string).second; 
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
 *   Determine if string begins with an integer and gives the integer if
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
 *   Returns true if given expression begins with a integer, false if not 
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

  while ((index < expression.size()) && (is_number)) {

    if ((expression[index]>='0') && (expression[index]<='9')) {
      nmbr += expression[index]; 
    }
    else {
      is_number = false;
    }

    ++index;
  }

  if (index == 0) {
    is_number = false; 
  }

  number = atoi(nmbr.c_str());
  return(is_number);
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
  assert( ((change>0) || (((int)depth+change)>=0)) );  
  depth += change;
  stack.growIfNeeded(depth);
}


/*-
 *-----------------------------------------------------------------------
 * clampFirstDecisionTree
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
RngDecisionTree::clampFirstDecisionTree()
{
  //////////////////////////////////////////////////////////////////////
  // first make sure this is a DT from file object
  //////////////////////////////////////////////////////////////////////
  if (!clampable()) {
    error("ERROR: can't call clampFirstDecisionTree() for non-file DT");
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
  clampNextDecisionTree();
}


/*-
 *-----------------------------------------------------------------------
 * clampNextDecisionTree
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
RngDecisionTree::clampNextDecisionTree()
{
  int    readDTNum;

  //////////////////////////////////////////////////////////////////////
  // first make sure this is a DT from file object
  //////////////////////////////////////////////////////////////////////
  if (!clampable()) {
    error("ERROR: can't call clampNextDecisionTree() for non-file DT");
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
  }

  // read in the rest of the DT.
  dtFile->read(curName,"cur name");
  dtFile->read(_numFeatures,"num feats");
  if (_numFeatures <= 0)
    error("ERROR: reading dynamic decision tree '%s' with current name '%s' from file '%s', but decision tree must have > 0 features",
	  name().c_str(),
	  curName.c_str(),
	  dtFile->fileName());
  rightMostLeaf = NULL;
  root = readRecurse(*dtFile,rightMostLeaf);
  Node *tmp = rightMostLeaf;
  Node *rightLeaf = NULL;
  while (tmp != NULL) {
    tmp->leafNode.nextLeaf = rightLeaf;
    rightLeaf = tmp;
    tmp = tmp->leafNode.prevLeaf;
  }
  leftMostLeaf = rightLeaf;
}


/*-
 *-----------------------------------------------------------------------
 * writeIndexFile 
 *   Writes an index file for a clampable DT 
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
 *   If it is to be used again clampFirstDecisionTree must be called. 
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

  assert(clampable()); 

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
    }

    // read in the rest of the DT.
    dtFile->read(curName,"cur name");
    dtFile->read(_numFeatures,"num feats");
    if (_numFeatures <= 0) {
      error("ERROR: reading dynamic decision tree '%s' with current name '%s' from file '%s', but decision tree must have > 0 features", 
        name().c_str(), curName.c_str(), dtFile->fileName());
    } 

    rightMostLeaf = NULL;
    root = readRecurse(*dtFile,rightMostLeaf);
    Node *tmp = rightMostLeaf;
    Node *rightLeaf = NULL;
    while (tmp != NULL) {
      tmp->leafNode.nextLeaf = rightLeaf;
      rightLeaf = tmp;
      tmp = tmp->leafNode.prevLeaf;
    }
    leftMostLeaf = rightLeaf;
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
  if (!clampable()) {
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
  if (n->nodeType != NonLeafNode) {
    os.space(depth*2);
    os.write(-1);
    os.write(n->leafNode.leafNodeString);
    os.nl();
  } else {
    os.space(depth*2);
    os.write(n->nonLeafNode.ftr,"writeRecurse, ftr");
    os.write(n->nonLeafNode.children.size(),"writeRecurse, numsplits");
    for (unsigned i=0;i<n->nonLeafNode.rngs.size();i++) {
      os.write(n->nonLeafNode.rngs[i]->rangeStr());
    }
    os.write(RNG_DECISION_TREE_DEF_STR);
    os.nl();
    for (unsigned i=0;i<n->nonLeafNode.children.size();i++) 
      writeRecurse(os,n->nonLeafNode.children[i],depth+1);
  }
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                    Querying the decision tree
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * query
 *      queries the decision tree
 *      NOTE: Any update here should also be done in the
 *      version that takes arrays of random variables below.
 * 
 * Preconditions:
 *      Must be a filled in DC structure.
 *
 * Postconditions:
 *      same as before.
 *
 * Side Effects:
 *      none.
 *
 * Results:
 *      The resulting value of the query.
 *
 *-----------------------------------------------------------------------
 */
leafNodeValType RngDecisionTree::query(const vector <int >& arr,
				       const vector <int > &cards)
{
  assert(0);
  assert ( unsigned(arr.size()) == _numFeatures );
  return queryRecurse(arr,cards,root);
}


/*-
 *-----------------------------------------------------------------------
 * queryRecurse
 *      support for querying the decision tree
 *      NOTE: Any update here should also be done in the
 *      version that takes arrays of random variables below.
 * 
 * Preconditions:
 *      Must be a filled in DC structure.
 *
 * Postconditions:
 *      same as query
 *
 * Side Effects:
 *      none.
 *
 * Results:
 *      The resulting value of the query starting at node.
 *
 *-----------------------------------------------------------------------
 */
leafNodeValType RngDecisionTree::queryRecurse(const vector < int >& arr,
					      const vector < int >& cards,
					      RngDecisionTree::Node *n)
{
  assert(0);
  leafNodeValType dummy;
  return(dummy); 
}


/*-
 *-----------------------------------------------------------------------
 * query & queryRecurse
 *      Exact same routine as above but for arrays of random variables.
 *      This code is duplicated here for simplicity and efficiency.
 *      NOTE: Any update here should also be done in the
 *      version that takes arrays of integers above.
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
leafNodeValType RngDecisionTree::query(const vector < RandomVariable* >& arr)
{
  assert ( unsigned(arr.size()) == _numFeatures );
  return queryRecurse(arr,root);
}

leafNodeValType RngDecisionTree::queryRecurse(const vector < RandomVariable* >& arr,
					      RngDecisionTree::Node *n)
{

  if (n->nodeType == LeafNodeVal) {
    return n->leafNode.value;
  } else if (n->nodeType == NonLeafNode) {

    assert ( n->nonLeafNode.ftr < int(arr.size()) );

    /*
     * assert ( arr[n->nonLeafNode.ftr]->val < 
     * RNG_DECISION_TREE_MAX_CARDINALITY );
     */

    const int val = arr[n->nonLeafNode.ftr]->val;

    // use a switch to knock off the short cases for
    // which we use simple linear search with little bookkeeping.
    if (n->nonLeafNode.rngs.size() < DT_SPLIT_SORT_THRESHOLD
	|| 
	!n->nonLeafNode.ordered) {
      for (unsigned i=0;i<n->nonLeafNode.rngs.size();i++) {
	// Do a linear search.
	if (n->nonLeafNode.rngs[i]->contains(val))
	  return queryRecurse(arr,n->nonLeafNode.children[i]);
      }
    } else {
      // eliminate simple boundary conditions.
      const int maxRngNum = n->nonLeafNode.rngs.size()-1;
      if (*(n->nonLeafNode.rngs[0]) > val)
	goto failedQuery;
      if (*(n->nonLeafNode.rngs[maxRngNum]) < val)
	goto failedQuery;

      // do binary search
      int l=0,u=maxRngNum;
      int m=0;
      while (l<=u) {
	// all these are conditional on the val being contained in the rng.
	// rngs[l] <= val && val <= rngs[u]  
	m = (l+u)/2; 
	if (*(n->nonLeafNode.rngs[m]) > val) 
	  // rngs[l] <= val && val < rngs[m]
	  u = m-1;
	// rngs[l] <= val && val <= rngs[u]
	else if (*(n->nonLeafNode.rngs[m]) < val)
	  // rngs[m] < val && val < rngs[u]
	  l=m+1;
	// rngs[l] < val && val < rngs[u]
	else {
	  // found potential range that might contain value
	  // since neither val < rng nor val > rng. 
	  if (n->nonLeafNode.rngs[m]->contains(val))
	    return queryRecurse(arr,n->nonLeafNode.children[m]);
	  break;
	}
      }
    }
    
    // failed lookup, so return must be the default one.
  failedQuery:

    // failed lookup, so return must be the default one.
    return queryRecurse(arr,n->nonLeafNode.children[n->nonLeafNode.rngs.size()]);
  } else if (n->nodeType == LeafNodeFullExpand) {

    leafNodeValType res = 0;
    for (unsigned i=0;i<arr.size()-1;i++) {
      unsigned val = arr[i]->val;
      unsigned card = arr[i]->cardinality;
      res = card*(res + val);
    }
    res += arr[arr.size()-1]->val;
    return res;

  } else {
    leafNodeValType answer;

    answer = n->leafNode.equation.evaluateFormula( arr );
    return(answer);
  }

}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                    Misc
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////





#ifdef MAIN

////////////////////////////////////////////////////////////////////
//        test code Support
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


char *dtStr1 =
"% this is a decision tree file\n"
"%\n"
"dt_name 3  % number of features\n"
"0 10 10 11 12 13 0:5 6:9 14 50: 15 default\n"
"  1 2 0:10 default\n"
"    2 2 0:10 default\n"
"      -1 expand\n"
"      -1 (p0+1)\n"
"    2 2 0:5 default\n"
"      -1 ( c0 + 1 )\n"
"      -1 (m0 +1)\n"
"  1 2 0:10 default\n"
"    2 2 0:10 default\n"
"      -1 (p0+p1+ 5)\n"
"      -1 6\n"
"    2 2 0:5 default\n"
"      -1 7\n"
"      -1 8\n"
"  -1 10 % when feature[0] = 10, map to 10 regardless of all else\n"
"  -1 11 % when feature[0] = 11, map to 11 regardless of all else\n"
"  -1 12 % when feature[0] = 12, map to 12 regardless of all else\n"
"  -1 13 % when feature[0] = 13, map to 13 regardless of all else\n"
"  -1 14 % when feature[0] = 14, map to 14 regardless of all else\n"
"  -1 15 % when feature[0] = 15, map to 15 regardless of all else\n"
"  -1 16 % when feature[0] >= 50, map to 16 regardless of all else\n"
"  1 2 0:10 default\n"
"    2 2 0:10 default\n"
"      -1 9\n"
"      -1 10\n"
"    2 2 0:5 default\n"
"      -1 11\n"
"      -1 12\n";


char *dtStr2 =
"% this is a decision tree file\n"
"%\n"
"dt_name 1  % number of features\n"
"0 14 10 17,19,21 11 12 13 0:5 6:9 14 50: 15 18,20 23,25,27 22,24,26,28 default\n"
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


int
main(int argc,char *argv[])
{

  // first write out the file
  if (argc == 1)
    {
      oDataStreamFile dtfile ("/tmp/foo.dt",false);
      dtfile.write(dtStr1);
    }

  char *file = "/tmp/foo.dt";
  if (argc > 1)
    // read it in again
    file = argv[1];
  iDataStreamFile is (file,false);

  RngDecisionTree dt;
  dt.read(is);

  printf("Found decision tree\n");
  oDataStreamFile os("-",false);
  dt.write(os);

  vector<int> vec;
  vector<int> card;
  vec.resize(dt.numFeatures());
  card.resize(dt.numFeatures());
  iDataStreamFile stin ("-",false,false);

  // first test iterating through all leaf values.
#if 0
  for (RngDecisionTree::iterator it = dt.begin();
       it != dt.end(); it++) {
    printf("leaf value = %d\n",it.value());
  }
#endif

  printf("Enter a length %d set of cardinalities:",dt.numFeatures());
  fflush(stdout);
  stin.read(card,dt.numFeatures());

  while (1) {
    printf("Enter a length %d intvec:",dt.numFeatures());
    fflush(stdout);
    stin.read(vec,dt.numFeatures());

    printf("Querying with vector and cards: ");
    fflush(stdout);
    for (unsigned i=0;i<dt.numFeatures();i++) {
      printf(" %d:%d",vec[i],card[i]);
    }
    printf("\n");
    printf("### RESULT ==> %d\n",dt.query(vec,card));
  }

}



#endif
