/*-
 * GMTK_RngDecisionTree
 *      General class to map from vectors of integers to some 
 *      basic type (such as int, float, etc.). Uses bp_ranges
 *      to form the queries.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */



#ifndef GMTK_RNG_DECISION_TREE_H
#define GMTK_RNG_DECISION_TREE_H


#include "logp.h"
#include "sArray.h"
#include "fileParser.h"
#include "bp_range.h"

#include "GMTK_RandomVariable.h"


/////////////////////////////////////////////////
// The maximum branching factor on any decision tree node.
// This may be safely increased, and is here just for
// debugging during reading to help enforce sanity.
#define RNG_DECISION_TREE_MAX_ARY 20

/////////////////////////////////////////////////
// The maximum "cardinality" of any of the integers 
// used for making decisions on. This may
// safely be increased.
#define RNG_DECISION_TREE_MAX_CARDINALITY 512

#define RNG_DECISION_TREE_DEF_STR "default"

template <class T = int>
class RngDecisionTree {

private:

protected:

  struct Node {
    bool leaf;
    struct NonLeafNode {
      int ftr;
      sArray< Node* > children;
      sArray< BP_Range* > rngs;
    } nonLeafNode;
    struct LeafNode {
      T value;
    } leafNode;
  };

  ///////////////////////////////////////////////////////////    
  // the size of the integer vector that we are making
  // decisions based on. 
  int _numFeatures;

  ///////////////////////////////////////////////////////////    
  // The root of the decision tree.
  Node *root;

  ///////////////////////////////////////////////////////////    
  // support for reading
  Node* readRecurse(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  // support for querying
  T queryRecurse(const sArray < int >& arr,
		 Node *n);


  ///////////////////////////////////////////////////////////    
  // support for destructor
  void destructorRecurse(Node *n);

public:

  RngDecisionTree() {}
  ~RngDecisionTree();

  ///////////////////////////////////////////////////////////    
  // read in the basic parameters, assuming file pointer 
  // is located at the correct position.
  void read(iDataStreamFile& is);
  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  void write(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  // Make a query and return the value corresponding to
  // the array of integers.
  T query(const sArray < int >& arr);


  ///////////////////////////////////////////////////////////    
  // Make a query and return the value corresponding to
  // the array of integers, from an array of random variables.
  T query(const sArray < RandomVariable* >& arr);


  ///////////////////////////////////////////////////////////    
  // make available the number of featurse.
  int numFeatures() { return _numFeatures; }

};



/*-
 *-----------------------------------------------------------------------
 * ~RngDecisionTree
 *      destructor
 * 
 * Preconditions:
 *      object must be "filled" in.
 *
 * Postconditions:
 *      sky is the limit
 *
 * Side Effects:
 *      all memory is gone.
 *
 * Results:
 *      nothing
 *
 *-----------------------------------------------------------------------
 */
template <class T>
RngDecisionTree<T>::~RngDecisionTree()
{
  destructorRecurse(root);
  delete root;
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
 *      sky is the limit
 *
 * Side Effects:
 *      all memory is gone.
 *
 * Results:
 *      nothing
 *
 *-----------------------------------------------------------------------
 */
template <class T>
void
RngDecisionTree<T>::destructorRecurse(RngDecisionTree<T>::Node* node)
{
  if (node->leaf) {
    // do nothing
  } else {
    for (int i=0;i<node->nonLeafNode.children.len();i++) {
      destructorRecurse(node->nonLeafNode.children[i]);
      delete node->nonLeafNode.children[i];
      if (i<node->nonLeafNode.children.len()-1)
	delete node->nonLeafNode.rngs[i];
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * read:
 *      reads in from a file.
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
 *      nothing.
 *
 *-----------------------------------------------------------------------
 */
template <class T> 
void
RngDecisionTree<T>::read(iDataStreamFile& is)
{
  is.read(_numFeatures,"RngDecisionTree:: read numFeatures");
  if (_numFeatures <= 0)
    error("RngDecisionTree::read decision tree must have >= 0 features");
  root = readRecurse(is);
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
 *      the current node it read
 *
 *-----------------------------------------------------------------------
 */
template <class T> 
RngDecisionTree<T>::Node* 
RngDecisionTree<T>::readRecurse(iDataStreamFile& is)
{

  int curFeat;
  is.read(curFeat,"RngDecisionTree:: readRecurse curFeat");
  if (curFeat < -1) 
    error("RngDecisionTree::readRecurse, feature number (=%d) must be non-negative",curFeat);
  if (curFeat >= _numFeatures) 
    error("RngDecisionTree::readRecurse, feature number (=%d) must be < numFeatures (=%d)",
	  _numFeatures);
  Node *node = new Node;
  if (curFeat == -1) {
    // leaf node
    node->leaf = true;
    is.read(node->leafNode.value,"RngDecisionTree:: readRecurse value");
  } else {
    node->leaf = false;
    node->nonLeafNode.ftr = curFeat;
    int numSplits;
    is.read(numSplits,"RngDecisionTree:: readRecurse numSplits");
    if (numSplits <= 1)
      error("RngDecisionTree:: readRecurse, can't have < 1 split");
    if (numSplits > RNG_DECISION_TREE_MAX_ARY)
      error("RngDecisionTree:: readRecurse, can't have > %d splits",
	    RNG_DECISION_TREE_MAX_ARY);
    node->nonLeafNode.children.resize(numSplits);
    // rngs is smaller since last string is always the default catch-all.
    node->nonLeafNode.rngs.resize(numSplits-1);
    for (int i=0;i<numSplits;i++) {
      char *str;
      is.read(str,"RngDecisionTree:: readRecurse, reading range");
      if (i==(numSplits-1)) {
	if (strcmp(RNG_DECISION_TREE_DEF_STR,str))
	  error("RngDecisionTree::readRecurse, expecting default str (%s) got (%s)",
		RNG_DECISION_TREE_DEF_STR,str);
      } else {
	node->nonLeafNode.rngs[i]
	  = new BP_Range(str,
			 0,
			 RNG_DECISION_TREE_MAX_CARDINALITY);
      }
      ///////////////////////////////////////////////////////////////
      // WARNING: We assume here that BP_Range will make
      // its own copy of the string, so we delete it here.
      // If the above class changes, we might need to change this code.
      delete [] str;
      ///////////////////////////////////////////////////////////////
    }
    // check for overlap errors in the strings
    for (int i=0;i<numSplits-1;i++) {
      for (int j=i+1;j<numSplits-1;j++) {
	if (node->nonLeafNode.rngs[i]
	    ->overlapP(
		       node->nonLeafNode.rngs[j]))
	  error("RngDecisionTree:: readRecurse, range %d (%s) and %d (%s) have a non-empty intersection.",i,
		node->nonLeafNode.rngs[i]->rangeStr(),j,
		node->nonLeafNode.rngs[j]->rangeStr());
      }
    }

    for (int i=0;i<numSplits;i++)
      node->nonLeafNode.children[i] =
	readRecurse(is);
  }
  return node;
}


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
template <class T> 
T RngDecisionTree<T>::query(const sArray < int >& arr)
{
  assert ( arr.len() == _numFeatures );
  return queryRecurse(arr,root);
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
 * Results:
 *      The resulting value of the query starting at node.
 *
 *-----------------------------------------------------------------------
 */
template <class T> 
T RngDecisionTree<T>::queryRecurse(const sArray < int >& arr,
				RngDecisionTree<T>::Node *n)
{
  if (n->leaf)
    return n->leafNode.value;

  assert ( n->nonLeafNode.ftr < arr.len() );
  assert ( arr[n->nonLeafNode.ftr] >= 0 &&
	   arr[n->nonLeafNode.ftr] <= 
	   RNG_DECISION_TREE_MAX_CARDINALITY );

  const int val = arr[n->nonLeafNode.ftr];
  for (int i=0;i<n->nonLeafNode.rngs.len();i++ ) {
    if (n->nonLeafNode.rngs[i]->contains(val))
      return queryRecurse(arr,n->nonLeafNode.children[i]);
  }
  // failed lookup, so return must be the default one.
  return queryRecurse(arr,n->nonLeafNode.children[n->nonLeafNode.rngs.len()]);
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
template <class T> 
T RngDecisionTree<T>::query(const sArray < RandomVariable* >& arr)
{
  assert ( arr.len() == _numFeatures );
  return queryRecurse(arr,root);
}
template <class T> 
T RngDecisionTree<T>::queryRecurse(const sArray < RandomVariable* >& arr,
				   RngDecisionTree<T>::Node *n)
{
  if (n->leaf)
    return n->leafNode.value;

  assert ( n->nonLeafNode.ftr < arr.len() );
  assert ( arr[n->nonLeafNode.ftr]->val >= 0 &&
	   arr[n->nonLeafNode.ftr]->val <= 
	   RNG_DECISION_TREE_MAX_CARDINALITY );

  const int val = arr[n->nonLeafNode.ftr]->val;
  for (int i=0;i<n->nonLeafNode.rngs.len();i++ ) {
    if (n->nonLeafNode.rngs[i]->contains(val))
      return queryRecurse(arr,n->nonLeafNode.children[i]);
  }
  // failed lookup, so return must be the default one.
  return queryRecurse(arr,n->nonLeafNode.children[n->nonLeafNode.rngs.len()]);
}



#endif // defined GMTK_CPT


