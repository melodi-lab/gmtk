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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"

#include "GMTK_RngDecisionTree.h"


VCID("$Header$");


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
 *      Returns false if everything is fine. If it returns
 *      true, then that means that this DT is not ready to be used,
 *      and clampNextDecisionTree() must be called (once for each
 *      DT that is given in the DT file).
 *
 *-----------------------------------------------------------------------
 */
void
RngDecisionTree::read(iDataStreamFile& is)
{
  NamedObject::read(is);

  // read the next entry which is either an integer, or
  // a string which gives the name of the DT file to obtain
  // multiple instances of the internals of this decision tree.
  is.read(dtFileName,"RngDecisionTree:: read file/numFeatures");
  if (!strIsInt(dtFileName.c_str(),(int*)&_numFeatures)) {
    if (clampable())
      error("ERROR: in DT named '%s' in file '%s', can't have DTs defined recursively in files",
	    name().c_str(),is.fileName());
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
void
RngDecisionTree::clampNextDecisionTree()
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
  assert ( arr.size() == cards.size() );

  if (n->nodeType == LeafNodeVal) {
    return n->leafNode.value;
  } else if (n->nodeType == NonLeafNode) {

    assert ( n->nonLeafNode.ftr < int(arr.size()) );

    assert ( arr[n->nonLeafNode.ftr] >= 0 );
    /*
      assert ( arr[n->nonLeafNode.ftr] >= 0 &&
      arr[n->nonLeafNode.ftr] <= 
      RNG_DECISION_TREE_MAX_CARDINALITY );
    */

    const int val = arr[n->nonLeafNode.ftr];

    // use a switch to knock off the short cases for
    // which we use simple linear search with little bookkeeping.
    if (n->nonLeafNode.rngs.size() < DT_SPLIT_SORT_THRESHOLD
	||
	!n->nonLeafNode.ordered) {
      for (unsigned i=0;i<n->nonLeafNode.rngs.size();i++) {
	// Do a linear search.
	if (n->nonLeafNode.rngs[i]->contains(val))
	  return queryRecurse(arr,cards,n->nonLeafNode.children[i]);
      }
    } else {
      // do a binary search.
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
	    return queryRecurse(arr,cards,n->nonLeafNode.children[m]);
	  break;
	}
      }
    }

    // failed lookup, so return must be the default one.
  failedQuery:
    return queryRecurse(arr,
			cards,
			n->nonLeafNode.children[n->nonLeafNode.rngs.size()]);
  } else if (n->nodeType == LeafNodeFullExpand) {

    leafNodeValType res = 0;
    for (unsigned i=0;i<arr.size()-1;i++) {
      res = cards[i]*(res + arr[i]);
    }
    res += arr[arr.size()-1];
    return res;

  } else {
    // formula

    leafNodeValType res = 0;
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
    // formula

    leafNodeValType res = 0;
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
