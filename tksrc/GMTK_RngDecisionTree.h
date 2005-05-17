/*-
 * GMTK_RngDecisionTree
 *      General class to map from vectors of integers to some 
 *      basic type (such as int, float, etc.). Uses bp_ranges
 *      to form the queries.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu> & Chris Bartels <bartels@ee.washington.edu>
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

#include "bp_range.h"
#include "fileParser.h"
#include "GMTK_NamedObject.h"
#include "logp.h"
#include "sArray.h"
#include "shash_map.h"

#include "GMTK_DiscRV.h"

#include <algorithm>
#include <map>
#include <vector>

class PackCliqueValue;

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

//////////////////////////////////////////////////////////////////
// Computes the max over N objects, for various N.
#define MAX_OF_2(a,b)          ((a)>(b)?(a):(b))
#define MAX_OF_3(a,b,c)        MAX_OF_2((a),MAX_OF_2((b),(c)))
#define MAX_OF_4(a,b,c,d)      MAX_OF_3((a),(b),MAX_OF_2((c),(d)))
#define MAX_OF_5(a,b,c,d,e)    MAX_OF_4((a),(b),(c),MAX_OF_2((d),(e)))
#define MAX_OF_6(a,b,c,d,e,f)  MAX_OF_5((a),(b),(c),(d),MAX_OF_2((e),(f)))


/////////////////////////////////////////////////
// Forward declare a RandomVariable so we can
// use pointers to it in this class.
class RV;
class DiscRV;

// unsigned cnt_NonLeafNodeRngsStruct = 0;
// unsigned cnt_NonLeafNodeArrayStruct = 0;
// unsigned cnt_Node = 0;
// unsigned cnt_RangeNode = 0;

/////////////////////////////////////////////////
typedef unsigned leafNodeValType;

class RngDecisionTree : public NamedObject, IM {

#ifdef MAIN
  friend int main(int,char**);
#endif

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

  template <class T>
  class sArrayStack : public sArray<T>
  {
    public:
      sArrayStack(int size) {
        stack_top = -1;
        this->resize(size);
      }
      
      inline void push_back(T item) {
        ++stack_top;
        this->ptr[stack_top] = item;
      }

      inline void pop_back() {
        --stack_top;
        assert(stack_top >= -1);
      }

      inline void pop_back(unsigned nmbr) {
        stack_top -= nmbr;
        assert(stack_top >= -1);
      }

      inline void clear() {
        stack_top = -1;
      }

      inline T stackSize() const {
        return(stack_top+1);
      }

    private:
      int stack_top;
  };

  class EquationClass
  {
  public:

    typedef unsigned formulaCommand;
    typedef sArray<formulaCommand> formulaCommandContainer;

    EquationClass();

    void parseFormula(
		      string leafNodeVal 
		      );

    leafNodeValType evaluateFormula(const vector< RV* >& variables,
				    const RV* rv = NULL);

  protected:

    typedef int stack_element_t;
    static sArrayStack<stack_element_t> stack;

    // Vector of commands 
    formulaCommandContainer commands;

    ///////////////////////////////////////////////////////////////////////
    // Tokens used by the equation parsing function  
    ///////////////////////////////////////////////////////////////////////    
    typedef enum {

      TOKEN_ABSOLUTE_VALUE, 
      TOKEN_ALL_DIFFERENT,
      TOKEN_BITWISE_AND, 
      TOKEN_BITWISE_OR, 
      TOKEN_BITWISE_NOT, 
      TOKEN_BITWISE_XOR,
      TOKEN_COLON, 
      TOKEN_COMMA, 
      TOKEN_DIVIDE_CEIL, 
      TOKEN_DIVIDE_FLOOR, 
      TOKEN_DIVIDE_ROUND, 
      TOKEN_DIVIDE_CEIL_FUNCTION, 
      TOKEN_DIVIDE_FLOOR_FUNCTION, 
      TOKEN_DIVIDE_ROUND_FUNCTION, 
      TOKEN_EQUALS, 
      TOKEN_EXPONENT, 
      TOKEN_GREATER_THAN, 
      TOKEN_GREATER_THAN_EQ, 
      TOKEN_INTEGER, 
      TOKEN_LEFT_PAREN, 
      TOKEN_LESS_THAN, 
      TOKEN_LESS_THAN_EQ, 
      TOKEN_LOGICAL_AND, 
      TOKEN_LOGICAL_OR, 
      TOKEN_MAX,
      TOKEN_MEDIAN,
      TOKEN_MIN,
      TOKEN_MINUS, 
      TOKEN_MOD,
      TOKEN_NOT,
      TOKEN_PLUS, 
      TOKEN_QUESTION_MARK, 
      TOKEN_RIGHT_PAREN, 
      TOKEN_ROTATE, 
      TOKEN_SHIFT_LEFT, 
      TOKEN_SHIFT_RIGHT, 
      TOKEN_SPACE, 
      TOKEN_TIMES, 

      TOKEN_CARDINALITY_CHILD, 
      TOKEN_CARDINALITY_PARENT, 
      TOKEN_MAX_VALUE_CHILD, 
      TOKEN_MAX_VALUE_PARENT, 
      TOKEN_PARENT_VALUE, 
      TOKEN_PARENT_VALUE_MINUS_ONE, 
      TOKEN_PARENT_VALUE_PLUS_ONE, 

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
      OPERAND_SHIFT = 6,
      COMMAND_MASK  = 0x3f, 
      OPERAND_MASK  = ~COMMAND_MASK,

      COMMAND_INVALID = 0,
      COMMAND_PUSH_CARDINALITY_CHILD, 
      COMMAND_PUSH_CARDINALITY_PARENT, 
      COMMAND_PUSH_PARENT_VALUE, 
      COMMAND_PUSH_PARENT_VALUE_MINUS_ONE, 
      COMMAND_PUSH_PARENT_VALUE_PLUS_ONE, 
      COMMAND_PUSH_MAX_VALUE_CHILD, 
      COMMAND_PUSH_MAX_VALUE_PARENT, 
      COMMAND_PUSH_CONSTANT,    
      COMMAND_ABSOLUTE_VALUE, 
      COMMAND_ALL_DIFFERENT, 
      COMMAND_BITWISE_AND, 
      COMMAND_BITWISE_NOT, 
      COMMAND_BITWISE_OR, 
      COMMAND_BITWISE_XOR, 
      COMMAND_BRANCH,
      COMMAND_BRANCH_IF_FALSE,
      COMMAND_DIVIDE_CEIL, 
      COMMAND_DIVIDE_FLOOR, 
      COMMAND_DIVIDE_ROUND, 
      COMMAND_EQUALS, 
      COMMAND_EXPONENT, 
      COMMAND_GREATER_THAN, 
      COMMAND_GREATER_THAN_EQ, 
      COMMAND_LESS_THAN, 
      COMMAND_LESS_THAN_EQ, 
      COMMAND_LOGICAL_AND, 
      COMMAND_LOGICAL_OR, 
      COMMAND_MAX,
      COMMAND_MEDIAN, 
      COMMAND_MIN, 
      COMMAND_MINUS, 
      COMMAND_MOD, 
      COMMAND_NEGATE, 
      COMMAND_NOT, 
      COMMAND_PLUS, 
      COMMAND_ROTATE,
      COMMAND_SHIFT_LEFT, 
      COMMAND_SHIFT_RIGHT, 
      COMMAND_TIMES, 

      LAST_COMMAND_INDEX
    };

#define MAKE_COMMAND(command, operand)  \
        (command | (operand << OPERAND_SHIFT))

#define GET_COMMAND(command) (command & COMMAND_MASK) 
#define GET_OPERAND(command) ((command & OPERAND_MASK) >> OPERAND_SHIFT)

    typedef enum {
      HIGHEST_PRECEDENCE = 0,
      PAREN_PRCDNC,  
      UNARY_PRCDNC,  
      EXPONENT_PRCDNC,  
      MULT_PRCDNC,  
      ADDITIVE_PRCDNC,  
      SHIFT_PRCDNC,  
      RELATIONAL_PRCDNC,  
      EQUALITY_PRCDNC,  
      BITWISE_AND_PRCDNC,  
      BITWISE_OR_PRCDNC,  
      LOGICAL_AND_PRCDNC,  
      LOGICAL_OR_PRCDNC,  
      CONDITIONAL_PRCDNC,
      LAST,
      LOWEST_PRECEDENCE = (LAST-1) 
    } precedence_t;

    static map<string, tokenEnum> delimiter;
    static map<string, tokenEnum> function;
    static map<string, tokenEnum> variable;

    static map<tokenEnum, formulaCommand> infixToken;
    static map<tokenEnum, formulaCommand> unaryToken;
    static map<tokenEnum, formulaCommand> functionToken;
    static map<tokenEnum, formulaCommand> oneValFunctionToken;
    static map<tokenEnum, formulaCommand> twoValFunctionToken;
    static map<tokenEnum, formulaCommand> manyValFunctionToken;
    static map<tokenEnum, formulaCommand> variableToken;
    static map<tokenEnum, unsigned>       tokenPriority;

    typedef vector<formulaCommand> parsingCommandContainer;

    void preProcessFormula(
			   string& original 
			   );

    void parseExpression(
			 tokenStruct&             token,
			 string&                  formula, 
			 parsingCommandContainer& commands,
			 unsigned                 precedence_level, 
			 unsigned&                depth  
			 );

    void parseQuestionMark(
			   tokenStruct&             token,
			   string&                  formula,  
			   parsingCommandContainer& commands,
			   unsigned&                depth  
			   );

    bool parseFactor(
		     tokenStruct&             token,
		     string&                  formula, 
		     parsingCommandContainer& commands,
		     unsigned&                depth  
		     );

    void getToken(
		  string&      expression, 
		  tokenStruct& item
		  );

    bool getInteger(
		    const string& expression, 
		    int&   number 
		    );

    void changeDepth(
		     const int change, 
		     unsigned& depth 
		     );

  };

  ///////////////////////////////////////////////////////////////    
  // Node Structure Types. 
  ///////////////////////////////////////////////////////////////    
  enum NodeType{NonLeafNodeArray, NonLeafNodeHash, NonLeafNodeRngs, 
		LeafNodeVal, LeafNodeEquation };

  // forward declaration.
  struct Node;

  // Note: none of the following structures can contain either
  // pointers to themselves or members which contain pointers to
  // themselves.


  ///////////////////////////////////////////////////////////    
  // Caess when node is not a leaf node 
  ///////////////////////////////////////////////////////////    
  struct NonLeafNodeArrayStruct {
     // Which feature (parent) are we querying
     int ftr;
     // value to subtract from parent to get index into array
     unsigned base;
     // last entry of children is the "default" node.
     sArray < Node > children;

    // NonLeafNodeArrayStruct() { printf("::NonLeafNodeArrayStruct, cnt = %d\n",cnt_NonLeafNodeArrayStruct++); }
    // ~NonLeafNodeArrayStruct() { printf("::~NonLeafNodeArrayStruct, cnt = %d\n",--cnt_NonLeafNodeArrayStruct); }

  };

  struct NonLeafNodeHashStruct {
    // Which feature (parent) are we querying
    int ftr;
    // last entry of children is the "default" case
    sArray < Node > children;
    // hash table from parent int value to pointer to node within children.
    // TODO: the hash member is big, try to use something smaller.
    shash_map < unsigned, Node* > nodeMapper;

    NonLeafNodeHashStruct(unsigned starting_size)
      : nodeMapper(starting_size) {}

  };


  // forward declaration
  struct RangeNode;

  struct NonLeafNodeRngsStruct {
    // Which feature (parent) are we querying
    int ftr;
    // true if the ranges are ordered, so we can do binsearch.
    bool ordered;
    // here, 'children' structure does not include "default" case.
    // so this is of size numSplits - 1.
    sArray < RangeNode > children;
    // default special case if all ranges don't match.
    Node* def;

    // NonLeafNodeRngsStruct() { printf("::NonLeafNodeRngsStruct, cnt = %d\n",cnt_NonLeafNodeRngsStruct++); }
    // ~NonLeafNodeRngsStruct() { printf("::~NonLeafNodeRngsStruct, cnt = %d\n",--cnt_NonLeafNodeRngsStruct); }
  };


  ///////////////////////////////////////////////////////////    
  // Structure for various types of leaf nodes 
  ///////////////////////////////////////////////////////////    
  struct LeafNodeValStruct {

    // The value, if it is just a value,
    // i.e., when nodeType == LeafNodeVal
    leafNodeValType value;
  } ;

  struct LeafNodeEquationStruct {
    // Used when the leaf node is an equation 
    EquationClass equation;
  } ;

  ///////////////////////////////////////////////////////////    
  // Node structure (nonLeafNode and leafNodes can not be 
  // unioned because they have entries with constructors)  
  // Therefore, we hack it using 'unsigned char foo[max(sizeof,sizeof)]' 
  // constructs and C++'s placement new/delete constructs. 
  // This significantly reduces memory requirements and
  // also allows for contiguous allocation of Nodes in
  // the hierarchy (further reducing space requirements).
  ///////////////////////////////////////////////////////////    
  struct Node {
    NodeType nodeType;

    unsigned char _buffer[MAX_OF_5(sizeof(NonLeafNodeArrayStruct),sizeof(NonLeafNodeHashStruct),sizeof(NonLeafNodeRngsStruct),sizeof(LeafNodeValStruct),sizeof(LeafNodeEquationStruct))];

    inline NonLeafNodeArrayStruct& nln_a()  { return (*((NonLeafNodeArrayStruct*)&_buffer)); }
    inline NonLeafNodeHashStruct&  nln_h()  { return (*((NonLeafNodeHashStruct*)&_buffer)); }
    inline NonLeafNodeRngsStruct&  nln_r()  { return (*((NonLeafNodeRngsStruct*)&_buffer)); }
    inline LeafNodeValStruct&      ln_v()   { return (*((LeafNodeValStruct*)&_buffer)); }
    inline LeafNodeEquationStruct& ln_e()   { return (*((LeafNodeEquationStruct*)&_buffer)); }


#if 0
    Node(const Node&nd) { 
      *this = nd;
    }
    Node& operator=(const Node&nd) {
      nodeType = nd.nodeType;
      ::memcpy(_buffer,nd._buffer,sizeof(_buffer));
      printf("--- Node(nd), op = %d\n",cnt_Node); 
      return *this;
    }

    Node() { printf("::Node(), cnt = %d\n",cnt_Node++); }
    ~Node() { printf("::~Node(), cnt = %d\n",--cnt_Node); }
#endif

  };

  // complete declaration of RangeNode mentioned above.
  struct RangeNode {
    BP_Range rng;
    Node nd;


    bool operator<(const RangeNode& rn) const { return (rng < rn.rng); }

#if 0
    RangeNode(const RangeNode& rn)
      : rng(rn.rng),nd(rn.nd) 
    {
      printf("::RangeNode(rn), cnt = %d\n",cnt_RangeNode++);
    }

    RangeNode& operator=(const RangeNode& rn) {
      rng = rn.rng; nd = rn.nd; return *this; 
      printf("-- RangeNode op  = %d\n",cnt_RangeNode);
    }
    RangeNode()
    {
      printf("::RangeNode(), cnt = %d\n",cnt_RangeNode++);
    }
    ~RangeNode()
    {
      printf("::~RangeNode(), cnt = %d\n",--cnt_RangeNode);
    }
#endif

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
  // support for reading
  void readRecurse(iDataStreamFile& is,Node& node);


  ///////////////////////////////////////////////////////////    
  // support for writing
  void writeRecurse(oDataStreamFile& os,
		    Node *n,
		    const int depth); 



  ///////////////////////////////////////////////////////////    
  // support for querying:
  // This routine assumes that all parent RV's may be safely cast into
  // DiscRVs. If this is not the case, arbitrary results will
  // follow. The child random variable "rv" may be any type of random
  // variable, but if it is actually continous, then any DT formula
  // asking for the child's cardinality will evaluate to 0.
  leafNodeValType queryRecurse(const vector < RV* >& arr,
			       Node *n,
			       const RV* rv);


  ///////////////////////////////////////////////////////////    
  // support for destructor
  void destructorRecurse(Node *n);

public:

  RngDecisionTree() : indexFile(NULL), dtFile(NULL), firstDT(0), root(NULL) {}; 
  ~RngDecisionTree();

  // return true if this DT changes from one segment to the next. We
  // know this by if the dtFile is available (if it is, presumably
  // this DT is iterable).
  bool iterable() { return (dtFile != NULL); }

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
  // seek through decision tree file to a particular tree 
  void seek(
    unsigned dt_nmbr 
  );
 
  ///////////////////////////////////////////////////////////    
  // Set the file pointer and read ino the first DT 
  void beginIterableDT();

  ///////////////////////////////////////////////////////////    
  // read in the next DT 
  void nextIterableDT();

  ///////////////////////////////////////////////////////////    
  // Function only used for testing 
  bool testFormula(
    string               formula,
    const vector< RV* >& variables, 
    RV*                  child, 
    unsigned             desired_answer 
    );

  ///////////////////////////////////////////////////////////    
  ///////////////////////////////////////////////////////////    
  // Make a query and return the value corresponding to the array of
  // integers, from an array of random variables.  Note that this
  // routine assumes that all parent RVs may be safely cast into
  // DiscRVs. If this is not the case, arbitrary results will
  // follow. If the child 'rv' is discrete, then any DT formulas will
  // have access to features of the current child random variable
  // (such as child's cardinality, frame number, and so
  // on). Otherwise, if 'drv' is continous, then the child cardinality
  // will evaluate to 0.
  leafNodeValType query(const vector < RV* >& arr,
			const RV* rv) {
    assert ( unsigned(arr.size()) == _numFeatures );
    return queryRecurse(arr,root,rv);
  }


  ///////////////////////////////////////////////////////////    
  // make available the number of features.
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
  // write index file for a iterable DT 
  void writeIndexFile();

  ///////////////////////////////////////////////////////////    
  // count the number of parent assignments that satisfy a 
  // observed child. All variables are discrete.
  void computeParentsSatisfyingChild(
	    // input arguments
	    unsigned par, // parent number
	    vector <RV*> & parents, 
	    vector <RV*> & hiddenParents,
	    PackCliqueValue& hiddenParentPacker,
	    sArray < DiscRVType*>& hiddenNodeValPtrs,
	    RV* child,
	    // output arguments
	    sArray < unsigned >& packedParentVals,
	    unsigned& num);
  void computeParentsChildSatisfyingGrandChild(
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
	    unsigned& num);


};


#endif 


