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
#include "bp_range.h"
#include "GMTK_NamedObject.h"

#include <vector>
#include <algorithm>

/////////////////////////////////////////////////
// The maximum branching factor on any decision tree node.
// This may be safely increased, and is here just for
// debugging during reading to help enforce sanity.
#define RNG_DECISION_TREE_MAX_ARY 400000

/////////////////////////////////////////////////
// The maximum "cardinality" of any of the integers 
// used for making decisions on. This may
// safely be increased if needed. 
#define RNG_DECISION_TREE_MAX_CARDINALITY ((1<<(sizeof(unsigned)*8-1))-1)

/////////////////////////////////////////////////
// The maximum range value
#define MAX_BP_RANGE_VALUE 100000000


/////////////////////////////////////////////////
// The string that is used to specify the 'fail' case
// in a decision tree. This condition will
// succeed only if none of the other rng conditions succeed.
#define RNG_DECISION_TREE_DEF_STR "default"


/////////////////////////////////////////////////////////////
// The threshold beyond which the children in each node
// of a DT are sorted when they are read in. I.e., 
// we sort if the number of splits is greater than or equal to this.
#define DT_SPLIT_SORT_THRESHOLD 3

/////////////////////////////////////////////////
// Forward declare a RandomVariable so we can
// use pointers to it in this class.
class RandomVariable;

/////////////////////////////////////////////////
typedef unsigned leafNodeValType;

class RngDecisionTree : public NamedObject {

private:

protected:

  ////////////////////////////////////
  // this is used if we are to obtain multiple
  // decision trees (using the same name) from
  // a file
  iDataStreamFile* indexFile;  // file pointer to the index file
  iDataStreamFile* dtFile;     // the file pointer
  string           dtFileName; // the file name
  unsigned         numDTs;     // number of DTs in this file
  int              dtNum;      // the current DT number
  string           curName;    // the current DT name
  unsigned         firstDT;    // index of the first decision tree
 
  ///////////////////////////////////////////////////////////////////////////    
  // Equation Parsing 
  ///////////////////////////////////////////////////////////////////////////    

  class EquationClass
  {

    public:

      void parseFormula(
        string& leafNodeVal 
        );

      leafNodeValType evaluateFormula(
        const vector< RandomVariable* >& variables 
      );

      leafNodeValType evaluateFormula(
        const vector<unsigned>& value, 
        const vector<unsigned>& cardinality
      );

    protected:

      struct commandStruct { 
        unsigned command;
        unsigned immediate;
      };

      // Vector of commands with immediates 
      vector<commandStruct> commands;

      ///////////////////////////////////////////////////////////////////////
      // Tokens used by the equation parsing function  
      ///////////////////////////////////////////////////////////////////////    
      typedef enum {

        TOKEN_SPACE, 
        TOKEN_PARENT, 
        TOKEN_CARDINALITY, 
        TOKEN_INTEGER, 
        TOKEN_PLUS, 
        TOKEN_MINUS, 
        TOKEN_TIMES, 
        TOKEN_DIVIDE, 
        TOKEN_LEFT_PAREN, 
        TOKEN_RIGHT_PAREN, 
        TOKEN_COMMA, 
        TOKEN_EXPONENT, 
        TOKEN_BITWISE_AND, 
        TOKEN_BITWISE_OR, 
        TOKEN_XOR,
        TOKEN_LOGICAL_AND, 
        TOKEN_LOGICAL_OR, 
        TOKEN_MOD,
        TOKEN_MIN,
        TOKEN_MAX,

        TOKEN_END, 
        LAST_TOKEN_INDEX

      } tokenEnum;

      /////////////////////////////////////////////////////////////////////// 
      // Structure used for storing tokens with associated values 
      /////////////////////////////////////////////////////////////////////// 
      struct tokenStruct {
        tokenEnum token;

        string text; 
        int    number; 
      }; 

      ///////////////////////////////////////////////////////////////////////
      // Commands used in equation evaluation 
      ///////////////////////////////////////////////////////////////////////    
      enum {

        COMMAND_PUSH_PARENT,      
        COMMAND_PUSH_CARDINALITY,
        COMMAND_PUSH_CONSTANT, 
        COMMAND_PLUS, 
        COMMAND_MINUS, 
        COMMAND_TIMES, 
        COMMAND_DIVIDE, 
        COMMAND_EXPONENT, 
        COMMAND_BITWISE_AND, 
        COMMAND_BITWISE_OR, 
        COMMAND_XOR, 
        COMMAND_LOGICAL_AND, 
        COMMAND_LOGICAL_OR, 
        COMMAND_MOD, 
        COMMAND_MIN, 
        COMMAND_MAX, 

        LAST_COMMAND_INDEX
      };

      void parseExpression(
        tokenStruct& token,
        string&      leafNodeVal 
      );
 
      void parseTerm(
        tokenStruct& token,
        string&      leafNodeVal 
      );

      void parseFactor(
        tokenStruct& token,
        string&      leafNodeVal 
      );

      void getToken(
        string&      expression, 
        tokenStruct& item
      );

      bool getInteger(
        const string& expression, 
        int&   number 
      );
  };

  ///////////////////////////////////////////////////////////////    
  // Node Structures 
  ///////////////////////////////////////////////////////////////    
  enum NodeType{NonLeafNode, LeafNodeVal, LeafNodeFormula, LeafNodeFullExpand};

  struct Node;

  ///////////////////////////////////////////////////////////    
  // Tree structure for when node is not a leaf node 
  ///////////////////////////////////////////////////////////    
  struct NonLeafNodeStruct {
      // This is when nodeType == NonLeafNode.
      int ftr;
      bool ordered;
      vector< Node* > children;
      vector< BP_Range* > rngs;
   } nonLeafNode;

  ///////////////////////////////////////////////////////////    
  // Structure for various types of leaf nodes 
  ///////////////////////////////////////////////////////////    
  struct LeafNodeStruct {

    // The string form of the leaf node.
    string leafNodeString;

    // The value, if it is just a value,
    // i.e., when nodeType == LeafNodeVal
    leafNodeValType value;

    // Note: when nodeType == LeafNodeFullExpand,
    // we just compute the one-to-one mapping
    // from joint state space of the parents
    // to the integers, so the above two things
    // are not used.

    // Used when the leaf node is an equation 
    EquationClass equation;

    // In order to have easy access to all of the
    // leaf nodes, we also keep a 
    // doubly linked list of these nodes.
    Node* prevLeaf;
    Node* nextLeaf;
  } ;

  ///////////////////////////////////////////////////////////    
  // Node structure (nonLeafNode and leafNodes can not be 
  //  unioned because they have entries with constructors)  
  ///////////////////////////////////////////////////////////    
  struct Node {
    NodeType nodeType;
    NonLeafNodeStruct nonLeafNode;
    LeafNodeStruct    leafNode;
  };

  ////////////////////////////////////////////////////////////////
  // structure for comparing range pointers. 
  struct RngCompare {  
    RngCompare() {}
    bool operator() (const pair<BP_Range*,Node*>& a, 
		     const pair<BP_Range*,Node*>& b) {
      return (*a.first) < (*b.first);
    }
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
  leafNodeValType queryRecurse(const vector < int >& arr,
		 const vector < int >& cards,
		 Node *n);
  leafNodeValType queryRecurse(const vector < RandomVariable* >& arr,
		 Node *n);


  ///////////////////////////////////////////////////////////    
  // support for destructor
  void destructorRecurse(Node *n);

  ///////////////////////////////////////////////////////////    
  // seek through decision tree file to a particular tree 
  void seek(
    unsigned dt_nmbr 
  );

public:

  RngDecisionTree() : indexFile(NULL), dtFile(NULL), firstDT(0), root(NULL) {}; 
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
  // Set the index of the first decision tree
  void setFirstDecisionTree(
    unsigned first_index
  ) { firstDT = first_index; } 
 
  ///////////////////////////////////////////////////////////    
  // Set the file pointer and read ino the first DT 
  void clampFirstDecisionTree();

  ///////////////////////////////////////////////////////////    
  // read in the next DT 
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
    leafNodeValType value() { return leaf->leafNode.value; }
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
  leafNodeValType query(const vector < int >& arr,
	  const vector <int > & cards);


  ///////////////////////////////////////////////////////////    
  // Make a query and return the value corresponding to
  // the array of integers, from an array of random variables.
  leafNodeValType query(const vector < RandomVariable* >& arr);


  ///////////////////////////////////////////////////////////    
  // make available the number of featurse.
  unsigned numFeatures() { return _numFeatures; }

  ///////////////////////////////////////////////////////////    
  // return the number of leaves.
  int numLeaves();

  ///////////////////////////////////////////////////////////    
  // ********************* 
  void initializeIterableDT(
    string& fileName
  ); 

  ///////////////////////////////////////////////////////////    
  // write index file for a clampable DT 
  void writeIndexFile();

};


#endif 


