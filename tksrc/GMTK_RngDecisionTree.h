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

/*
 *
 * TODO: modify so that the resulting leaf node values can
 * just copy one of the parent values (making it easy
 * to have one of the parent values just specify the index).
 * TODO: make it such that a leaf node is some simple integer formula
 * on the parent values (or the result is an int that can be
 * type cast to a T).
 * TODO: make this not a template class as it is only used with
 * unsigned.
 * 
 */



#ifndef GMTK_RNG_DECISION_TREE_H
#define GMTK_RNG_DECISION_TREE_H

#include "fileParser.h"
#include "logp.h"
#include "vector.h"
#include "bp_range.h"
#include "GMTK_NamedObject.h"

/////////////////////////////////////////////////
// The maximum branching factor on any decision tree node.
// This may be safely increased, and is here just for
// debugging during reading to help enforce sanity.
#define RNG_DECISION_TREE_MAX_ARY 100

/////////////////////////////////////////////////
// The maximum "cardinality" of any of the integers 
// used for making decisions on. This may
// safely be increased.
#define RNG_DECISION_TREE_MAX_CARDINALITY 512


/////////////////////////////////////////////////
// The string that is used to specify the 'fail' case
// in a decision tree. This condition will
// succeed only if none of the other rng conditions succeed.
#define RNG_DECISION_TREE_DEF_STR "default"

/////////////////////////////////////////////////
// Forward declare a RandomVariable so we can
// use pointers to it in this class.
class RandomVariable;

template <class T = int>
class RngDecisionTree : public NamedObject {

private:

protected:

  ////////////////////////////////////
  // this is used if we are to obtain multiple
  // decision trees (using the same name) from
  // a file
  iDataStreamFile* dtFile; // the file pointer
  string dtFileName;       // the file name
  unsigned numDTs;         // number of DTs in this file
  int dtNum;               // the current DT number
  string curName;          // the current DT name

  enum NodeType { NonLeafNode, LeafNodeVal, LeafNodeFormula, LeafNodeFullExpand };

  enum FormulaType { Integer, ParentValue, ParentCardinality, ParentValCard };

  struct FormEntry {
    FormulaType fType;
    union {
      unsigned parentIndex;
      T value;
    };
  };

  struct Node {
    NodeType nodeType;
    // Ideally, this should be a union of
    // the following two structs, but
    // either the compiler or C++ doesn't let
    // us keep a union of structs with templates as below.
    struct NonLeafNode {
      // This is when nodeType == NonLeafNode.
      int ftr;
      vector< Node* > children;
      vector< BP_Range* > rngs;
    } nonLeafNode;
    struct LeafNode {
      // The string form of the leaf node.
      string leafNodeString;
      // The value, if it is just a value,
      // i.e., when nodeType == LeafNodeVal
      T value;
      // The formula if it is a formula,
      // i.e., when nodeType == LeafNodeFormula
      vector <FormEntry> formula;
      // Note: when nodeType == LeafNodeFullExpand,
      // we just compute the one-to-one mapping
      // from joint state space of the parents
      // to the integers, so the above two things
      // are not used.

      // In order to have easy access to all of the
      // leaf nodes, we also keep a 
      // doubly linked list of these nodes.
      Node* prevLeaf;
      Node* nextLeaf;
    } leafNode;
  };

  ///////////////////////////////////////////////////////////    
  // the size of the integer vector that we are making
  // decisions based on. 
  unsigned _numFeatures;

  ///////////////////////////////////////////////////////////    
  // The root of the decision tree.
  Node *root;

  ///////////////////////////////////////////////////////////    
  // The leaf nodes.
  Node *rightMostLeaf;
  Node *leftMostLeaf;

  ///////////////////////////////////////////////////////////    
  // support for reading
  Node* readRecurse(iDataStreamFile& is,
		    Node* & prevLeaf);

  ///////////////////////////////////////////////////////////    
  // support for writing
  void writeRecurse(oDataStreamFile& os,
		     Node *n,
		     const int depth);

  ///////////////////////////////////////////////////////////    
  // support for querying
  T queryRecurse(const vector < int >& arr,
		 const vector < int >& cards,
		 Node *n);
  T queryRecurse(const vector < RandomVariable* >& arr,
		 Node *n);


  ///////////////////////////////////////////////////////////    
  // support for destructor
  void destructorRecurse(Node *n);

public:

  RngDecisionTree() { dtFile = NULL; root = NULL;  }
  ~RngDecisionTree();
  bool clampable() { return (dtFile != NULL); }


  ///////////////////////////////////////////////////////////    
  // read in the basic parameters, assuming file pointer 
  // is located at the correct position. Returns true
  // if this is a variable DT (that comes from a file).
  void read(iDataStreamFile& is);
  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  void write(oDataStreamFile& os);
  ///////////////////////////////////////////////////////////    
  // read in the next DT assuming this is a file.a
  void clampFirstDecisionTree();
  void clampNextDecisionTree();

  ///////////////////////////////////////////////////////////    
  // iterators for iterating through leaf values.

  class iterator {
    friend class RngDecisionTree;
    Node *leaf;
    const RngDecisionTree& mydt;
  public:
    iterator(RngDecisionTree& dt) : mydt(dt) { 
      leaf = dt.leftMostLeaf;
      // We could also start with rightMostLeaf and
      // move leftwards.
    }
    // prefix
    iterator& operator ++() { leaf = leaf->leafNode.nextLeaf; return *this; }
    // postfix
    iterator operator ++(int) { iterator tmp=*this; ++*this; return tmp; }
    T value() { return leaf->leafNode.value; }
    bool valueNode() { return (leaf->nodeType == LeafNodeVal); }
    bool operator == (const iterator &it) { return it.leaf == leaf; } 
    bool operator != (const iterator &it) { return it.leaf != leaf; } 
  };
  friend class iterator;
  
  iterator begin() { return iterator(*this); }
  iterator end() { iterator it(*this); it.leaf = NULL; return it; }


  ///////////////////////////////////////////////////////////    
  // Make a query and return the value corresponding to
  // the array of integers.
  T query(const vector < int >& arr,
	  const vector <int > & cards);


  ///////////////////////////////////////////////////////////    
  // Make a query and return the value corresponding to
  // the array of integers, from an array of random variables.
  T query(const vector < RandomVariable* >& arr);


  ///////////////////////////////////////////////////////////    
  // make available the number of featurse.
  unsigned numFeatures() { return _numFeatures; }

  ///////////////////////////////////////////////////////////    
  // return the number of leaves.
  int numLeaves();

};


#include "GMTK_RandomVariable.h"

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
template <class T>
RngDecisionTree<T>::~RngDecisionTree()
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
template <class T>
void
RngDecisionTree<T>::destructorRecurse(RngDecisionTree<T>::Node* node)
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
 *      Returns false if everything is fine. If it returns
 *      true, then that means that this DT is not ready to be used,
 *      and clampNextDecisionTree() must be called (once for each
 *      DT that is given in the DT file).
 *
 *-----------------------------------------------------------------------
 */
template <class T> 
void
RngDecisionTree<T>::read(iDataStreamFile& is)
{
  NamedObject::read(is);

  // read the next entry which is either an integer, or
  // a string which gives the name of the DT file to obtain
  // multiple instances of the internals of this decision tree.
  is.read(dtFileName,"RngDecisionTree:: read file/numFeatures");
  if (!strIsInt(dtFileName.c_str(),&(int)_numFeatures)) {
    if (clampable())
      error("ERROR: in DT named '%s', can't have DTs defined recursively from files",name().c_str());
    // then this must be a file name, turn off cpp pipe since we need to
    // rewind later.
    dtFile = new iDataStreamFile(dtFileName.c_str(),is.binary(),false);
    dtNum = -1;
    numDTs = 0;
    // pre-read here, so we always have a valid DT, 
    // even before clamping.
    dtFile->read(numDTs,"num DTs");
    if (numDTs == 0)
      error("ERROR: in DT named '%s', File '%s' specifies an invalid number of DTs\n",name().c_str(),dtFileName.c_str());
    clampNextDecisionTree();
    return;
  }
  if (_numFeatures < 0)
    error("ERROR: in DT named '%s', file '%s', decision tree must have >= 0 features",name().c_str(),is.fileName());
  rightMostLeaf = NULL;
  root = readRecurse(is,rightMostLeaf);
  Node *tmp = rightMostLeaf;
  Node *rightLeaf = NULL;
  while (tmp != NULL) {
    tmp->leafNode.nextLeaf = rightLeaf;
    rightLeaf = tmp;
    tmp = tmp->leafNode.prevLeaf;
  }
  leftMostLeaf = rightLeaf;
  return;
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
RngDecisionTree<T>::readRecurse(iDataStreamFile& is,
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
    node->leafNode.value = 0;
    // leaf node, now read the next bit as a string to determine what it is.
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
    } else if (leafNodeVal[0] = '(') {
      node->nodeType = LeafNodeFormula;
      //
      // try and parse it. It should be
      // a string of the form (A1+A2+A3+A4+...An) 
      // where Ai is either:
      //    1) a single integer
      //    2) the string "pj" where j is the j'th parent's value
      //    3) the string "cj" where j is the j'th parent's cardinality
      //    4) the string "mj" where j is the j'th parent's cardinality*value (i.e., the product of value * card)
      // 
      //  So, an example might be m0+m1+p2+1
      //    which will go from the joint product state space of parents 0, 1, and 2 but will add one 
      //    to the current state.
      // 


      // start parsing the string
      unsigned pos = 1;
      bool expectingPlus = false;
      node->leafNode.leafNodeString = leafNodeVal;
      while (1) {
	while (pos < leafNodeVal.size() && (leafNodeVal[pos] == ' ' || leafNodeVal[pos] == '\t'))
	  pos++;
	if (pos == leafNodeVal.size()) {
	  // get next portion
	  is.read(leafNodeVal,"RngDecisionTree:: readRecurse leaf node value");
	  node->leafNode.leafNodeString += leafNodeVal;
	  pos = 0;
	  continue;
	}
	int val;
	int len;
	if (leafNodeVal[pos] == 'p' || 
	    leafNodeVal[pos] == 'c' || 
	    leafNodeVal[pos] == 'm') {
	  if (expectingPlus) {
	    error("ERROR: DT '%s', file '%s': Expecting '+' at position %d in string '%s'\n",
		  name().c_str(),is.fileName(),pos,leafNodeVal.c_str());
	  }
	  FormEntry fe;
	  if (leafNodeVal[pos] == 'p')
	    fe.fType = ParentValue;
	  else if (leafNodeVal[pos] == 'c')
	    fe.fType = ParentCardinality;
	  else
	    fe.fType = ParentValCard;

	  pos++;
	  if (!strIsInt(leafNodeVal.c_str()+pos,&val,&len)) {
	    error("ERROR: DT '%s', file '%s': Expecting integer at position %d in string '%s'\n",
		  name().c_str(),is.fileName(),pos,leafNodeVal.c_str());
	  }
	  // we have an int value
	  pos += len;
	  fe.parentIndex = val;
	  if (fe.parentIndex >= _numFeatures)
	    error("ERROR: DT '%s' in file '%s' specifies a parent index (%c%d) larger than number of parents (%d)\n",
		  name().c_str(),is.fileName(),
		  fe.parentIndex,
		  (fe.fType == ParentValue?'p':(fe.fType == ParentCardinality?'c':'m')),
		  _numFeatures);
	  node->leafNode.formula.push_back(fe);
	  expectingPlus = true;
	} else if (leafNodeVal[pos] == '+') {
	  if (!expectingPlus)
	    error("ERROR: DT '%s', file '%s': Not expecting '+' at position %d in string '%s'\n",
		  name().c_str(),is.fileName(),pos,leafNodeVal.c_str());
	  // consume token
	  pos++;
	  expectingPlus = false;
	} else if (strIsInt(leafNodeVal.c_str()+pos,&val,&len)) {
	  if (expectingPlus) {
	    error("ERROR: DT '%s', file '%s': Expecting '+' at position %d in string '%s'\n",
		  name().c_str(),is.fileName(),pos,leafNodeVal.c_str());
	  }
	  FormEntry fe;
	  fe.fType = Integer;
	  fe.value = val;
	  node->leafNode.formula.push_back(fe);
	  pos += len;
	  expectingPlus = true;
	} else if (leafNodeVal[pos] == ')') {
	  // the end of the formula
	  break;
	} else {
	  error("ERROR: DT '%s', file '%s': Parse error at position %d in string '%s'\n",
		name().c_str(),is.fileName(),pos,leafNodeVal.c_str());
	}
      }
      // make sure we didnt' have string such as "p1+p2+"
      if (!expectingPlus) 
	error("ERROR: DT '%s', file '%s': Must not end with a '+' at position %d in string '%s'\n",
	      name().c_str(),is.fileName(),pos,leafNodeVal.c_str());
      if (node->leafNode.formula.size() == 0) {
	error("ERROR: DT '%s', file '%s': Must have formula between '()' in '%s'\n",
	      name().c_str(),is.fileName(),leafNodeVal.c_str());
      }
    }

    node->leafNode.prevLeaf = prevLeaf;
    prevLeaf = node;
  } else {
    node->nodeType = NonLeafNode;
    node->nonLeafNode.ftr = (T)curFeat;
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
    for (unsigned i=0;i<numSplits-1;i++) {
      for (unsigned j=i+1;j<numSplits-1;j++) {
	if (node->nonLeafNode.rngs[i]
	    ->overlapP(
		       node->nonLeafNode.rngs[j]))
	  error("ERROR: DT '%s', file '%s': range %d (%s) and %d (%s) have a non-empty intersection.",name().c_str(),is.fileName(),
		i,
		node->nonLeafNode.rngs[i]->rangeStr(),j,
		node->nonLeafNode.rngs[j]->rangeStr());
      }
    }

    for (unsigned i=0;i<numSplits;i++)
      node->nonLeafNode.children[i] =
	readRecurse(is,prevLeaf);
  }
  return node;
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
template <class T> 
void
RngDecisionTree<T>::clampFirstDecisionTree()
{
  // first make sure this is a DT from file object
  if (!clampable())
    error("ERROR: can't call clampFirstDecisionTree() for non-file DT");
  dtNum = -1;
  dtFile->rewind();
  dtFile->read(numDTs,"num DTs");
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
template <class T> 
void
RngDecisionTree<T>::clampNextDecisionTree()
{
  // first make sure this is a DT from file object
  if (!clampable())
    error("ERROR: can't call clampNextDecisionTree() for non-file DT");
  dtNum++;
  int i;
  dtFile->read(i,"DT num");
  if (i != dtNum)
    error("ERROR: reading from file '%s', expecting DT number %d but got number %d\n",dtFileName.c_str(),dtNum,i);
  // next, delete old DT if it exists.
  if (root != NULL) {
    destructorRecurse(root);
    delete root;
  }
  // read in the rest of the DT.
  dtFile->read(curName,"cur name");
  dtFile->read(_numFeatures,"num feats");
  if (_numFeatures <= 0)
    error("ERROR: reading decision tree from file '%s', but decision tree must have >= 0 features",
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
template <class T> 
void
RngDecisionTree<T>::write(oDataStreamFile& os)
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
template <class T> 
void
RngDecisionTree<T>::writeRecurse(oDataStreamFile& os,
				 RngDecisionTree<T>::Node *n,
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
template <class T> 
T RngDecisionTree<T>::query(const vector < int >& arr,
			    const vector <int > &cards)
{
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
template <class T> 
T RngDecisionTree<T>::queryRecurse(const vector < int >& arr,
				   const vector < int >& cards,
				   RngDecisionTree<T>::Node *n)
{
  assert ( arr.size() == cards.size() );

  if (n->nodeType == LeafNodeVal) {
    return n->leafNode.value;
  } else if (n->nodeType == NonLeafNode) {

    assert ( n->nonLeafNode.ftr < int(arr.size()) );
    assert ( arr[n->nonLeafNode.ftr] >= 0 &&
	     arr[n->nonLeafNode.ftr] <= 
	     RNG_DECISION_TREE_MAX_CARDINALITY );

    const int val = arr[n->nonLeafNode.ftr];
    for (unsigned i=0;i<n->nonLeafNode.rngs.size();i++ ) {
      // TODO: turn this into a log(n) operation
      // rather than linear. To do this, we'll need
      // to make sure that the ranges are in order.
      // To do this, define an operator < for bp_ranges.
      if (n->nonLeafNode.rngs[i]->contains(val))
	return queryRecurse(arr,cards,n->nonLeafNode.children[i]);
    }
    // failed lookup, so return must be the default one.
    return queryRecurse(arr,
			cards,
			n->nonLeafNode.children[n->nonLeafNode.rngs.size()]);
  } else if (n->nodeType == LeafNodeFullExpand) {

    T res = 0;
    for (unsigned i=0;i<arr.size()-1;i++) {
      res = cards[i]*(res + arr[i]);
    }
    res += arr[arr.size()-1];
    return res;

  } else {
    // formula

    T res = 0;
    for (unsigned i=0;i<n->leafNode.formula.size();i++) {
      if (n->leafNode.formula[i].fType == Integer)
	res += n->leafNode.formula[i].value;
      else if (n->leafNode.formula[i].fType == ParentValue) {
	if (n->leafNode.formula[i].parentIndex >= arr.size())
	  error("ERROR: Dynamic error, DT '%s' specifies a parent index (p%d) larger than number of parents (%d)\n",
		name().c_str(),
		n->leafNode.formula[i].parentIndex,
		arr.size());
	res += arr[n->leafNode.formula[i].parentIndex];
      } else if (n->leafNode.formula[i].fType == ParentCardinality) {
	if (n->leafNode.formula[i].parentIndex >= arr.size())
	  error("ERROR: Dynamic error, DT '%s' specifies a card index (c%d) larger than number of parents (%d)\n",
		name().c_str(),
		n->leafNode.formula[i].parentIndex,
		arr.size());
	res += cards[n->leafNode.formula[i].parentIndex];
      } else {
	// we must have that (n->leafNode.formula[i].fType == ParentValCard)
	if (n->leafNode.formula[i].parentIndex >= arr.size())
	  error("ERROR: Dynamic error, DT '%s' specifies an index (m%d) larger than number of parents (%d)\n",
		name().c_str(),
		n->leafNode.formula[i].parentIndex,
		arr.size());
	res += cards[n->leafNode.formula[i].parentIndex]*
	       arr[n->leafNode.formula[i].parentIndex];
      }
    }
    return res;
  }
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
T RngDecisionTree<T>::query(const vector < RandomVariable* >& arr)
{
  assert ( unsigned(arr.size()) == _numFeatures );
  return queryRecurse(arr,root);
}
template <class T> 
T RngDecisionTree<T>::queryRecurse(const vector < RandomVariable* >& arr,
				   RngDecisionTree<T>::Node *n)
{

  if (n->nodeType == LeafNodeVal) {
    return n->leafNode.value;
  } else if (n->nodeType == NonLeafNode) {

    assert ( n->nonLeafNode.ftr < int(arr.size()) );
    assert ( arr[n->nonLeafNode.ftr]->val >= 0 &&
	     arr[n->nonLeafNode.ftr]->val <= 
	     RNG_DECISION_TREE_MAX_CARDINALITY );

    const int val = arr[n->nonLeafNode.ftr]->val;
    for (unsigned i=0;i<n->nonLeafNode.rngs.size();i++ ) {
      if (n->nonLeafNode.rngs[i]->contains(val))
	return queryRecurse(arr,n->nonLeafNode.children[i]);
    }
    // failed lookup, so return must be the default one.
    return queryRecurse(arr,n->nonLeafNode.children[n->nonLeafNode.rngs.size()]);
  } else if (n->nodeType == LeafNodeFullExpand) {

    T res = 0;
    for (unsigned i=0;i<arr.size()-1;i++) {
      unsigned val = arr[i]->val;
      unsigned card = arr[i]->cardinality;
      res = card*(res + val);
    }
    res += arr[arr.size()-1]->val;
    return res;

  } else {
    // formula

    T res = 0;
    for (unsigned i=0;i<n->leafNode.formula.size();i++) {
      if (n->leafNode.formula[i].fType == Integer)
	res += n->leafNode.formula[i].value;
      else if (n->leafNode.formula[i].fType == ParentValue) {
	// This check is redundant.
	// if (n->leafNode.formula[i].parentIndex >= arr.size())
	//  error("ERROR: Dynamic error, DT '%s' specifies a parent index (p%d) larger than number of parents (%d)\n",
	//	name().c_str(),
	//	n->leafNode.formula[i].parentIndex,
	//	arr.size());
	res += arr[n->leafNode.formula[i].parentIndex]->val;
      } else if (n->leafNode.formula[i].fType == ParentCardinality) {
	// This check is redundant.
	// if (n->leafNode.formula[i].parentIndex >= arr.size())
	//  error("ERROR: Dynamic error, DT '%s' specifies a card index (c%d) larger than number of parents (%d)\n",
	//	name().c_str(),
	//	n->leafNode.formula[i].parentIndex,
	//	arr.size());
	res += arr[n->leafNode.formula[i].parentIndex]->cardinality;
      } else {
	// we must have that (n->leafNode.formula[i].fType == ParentValCard)
	// This check is redundant.
	// if (n->leafNode.formula[i].parentIndex >= arr.size())
	//  error("ERROR: Dynamic error, DT '%s' specifies an index (m%d) larger than number of parents (%d)\n",
	//	name().c_str(),
	//	n->leafNode.formula[i].parentIndex,
	//	arr.size());
	res += arr[n->leafNode.formula[i].parentIndex]->cardinality*
	       arr[n->leafNode.formula[i].parentIndex]->val;
      }
    }
    return res;
  }

}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                    Misc
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



#endif 


