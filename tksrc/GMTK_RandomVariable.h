/*
 * "Copyright 2001, University of Washington and International Business Machines Corporation. All Rights Reserved
 *
 *    Written by Jeff Bilmes and Geoffrey Zweig
 *
 * NO WARRANTY
 * THE PROGRAM IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED INCLUDING, WITHOUT
 * LIMITATION, ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT,
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Each Recipient is
 * solely responsible for determining the appropriateness of using the Program
 * and assumes all risks associated with such use, including but not limited
 * to the risks and costs of program errors, compliance with applicable laws,
 * damage to or loss of data, programs or equipment, and unavailability or
 * interruption of operations.

 * DISCLAIMER OF LIABILITY
 * THE UNIVERSITY OF WASHINGTON, INTERNATIONAL BUSINESS MACHINES CORPORATION,
 * JEFF BILMES AND GEOFFREY ZWEIG SHALL NOT HAVE ANY LIABILITY FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING WITHOUT LIMITATION LOST PROFITS), HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE  OF
 * THE PROGRAM, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES."
*/

/* 
 * GMTK_RandomVariable.h
 * Defines the functions that all random variables must satisfy in order
 * for the inference, sampling, and other generic routines to work.
 */

#ifndef GMTK_RANDOMVARIABLE_H
#define GMTK_RANDOMVARIABLE_H

#include <vector>
#include <string>

#include "logp.h"
#include "GMTK_RngDecisionTree.h"

// Variables come in these basic varieties.
enum vartype {Continuous,Discrete};

/////////////////////////////////////////////////////////////////////////////
// For doing EM with simulated data, it is necessary to be able to store 
// the values of
// all the variables in an array of type VariableValue. The easiest way
// of accomodating both discrete and continuous variables is to let 
// VariableValue objects store either kind of data.
struct VariableValue
{
    int ival;  // in the case of integer values, this gets used

    // for continuous variables, the values can be stored in the fval array
    vector<float> fval;   
};

/*
 *
 * The random variable class defines the basic functions that random
 * variables must implement. 
 *
 * Switching parents are guaranteed to be discrete.
 * A nodes conditioning parents must be lower indexed than it.
 *
 ******************************************************************
 * TODO: RandomVariable and it's descendants needs to be re-thought
 * as right now there is much that really belongs in children
 * (DiscreteRandomVariable, etc.).
 ******************************************************************
 *
 */
class RandomVariable
{
public:


  //////////////////////////////////////////////////////////////////////////////
  // This is the integer type of the values that a discrete random variable
  // may take on. Possibilities include unsigned char, char, short, int, 
  // unsigned long, and so on. Note, if this changes,
  // then might need to change Decision tree code as well.
  typedef short DiscreteVariableType;

  ////////////////////////////////////////////////////////////////////////
  // What is my name?
  string label;
  string name() { return label; }

  /////////////////////////////////////////////////////////////////////////
  // Is the node hidden, or is it an observation.
  bool hidden;

  /////////////////////////////////////////////////////////////////////////
  // Is the variable discrete?
  // Inference conditions on this; cliques keep track of the values of 
  // their discrete members.
  bool discrete;

  /////////////////////////////////////////////////////////////////////////
  // What time frame does the node belong to?
  // Counting starts from 0.
  int timeIndex;

  /////////////////////////////////////////////////////////////////////////
  // Initialize with the variable type.
  // The default timeIndex value of -1 indicates a static network.
  // The default value of "hidden" is true.
  // Discrete nodes must be specified with their cardinalities.
  RandomVariable(string _label, vartype vt, int card=0) {
    hidden=true; 
    discrete=(vt==Discrete); 
    timeIndex=-1;
    cardinality=card; 
    if (discrete) 
      assert(card!=0);
    label=_label;
    dtMapper=NULL;
  }

  virtual ~RandomVariable() {;}

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // PARENT MEMBER SUPPORT                      ////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // The set of switching parents, if any.
  vector<RandomVariable *> switchingParents;

  ////////////////////////////////////////////////////////////////////////
  // allPossibleParents is the union of the switchingParents and the
  // all possible conditionalParents.
  // Used to determine topological orderings
  vector<RandomVariable *> allPossibleParents;

  ////////////////////////////////////////////////////////////////////////
  // The set of variables that use this variable either as a switching
  // parent or as a (possible) conditional parent.
  vector<RandomVariable *> allPossibleChildren;

  ////////////////////////////////////////////////////////////////////////
  // For each possible different list of
  // conditional parents that might exist for all possible
  // values of the switching parents, this array gives that list of
  // appropriate conditional parents. For example, 
  // suppose that S is the set of conditional parents, and that
  // 0 <= S <= 5 corresponds to one set of conditional parents,
  // and 6 <= S < = 10 corresponds to another set of conditional
  // parents, and those are the only two set of conditional parents
  // that exist for all values of the switching parents, this
  // list is of size two. 
  vector< vector < RandomVariable* > > conditionalParentsList;

  ////////////////////////////////////////////////////////////////////////
  // This is set to the current set of conditional parents,
  // which is dependent on the current value of the switching parents.
  // Note that this points to one of the entries in conditionalParentsList
  vector<RandomVariable *> *curConditionalParents;

  ////////////////////////////////////////////////////////////////////////
  // This decision tree is used by the intFromSwitchingState() routine
  // to map from the current set of switching parent values
  // to an integer, which indicates which set of conditional
  // parents should be active for those switching parent values.
  // This variable is NULL if there are no Switching 
  // parents of this node.
  RngDecisionTree *dtMapper;

  ////////////////////////////////////////////////////////////////////////
  // A method to make it easy to set up the graph topology
  // The probability of a variable taking a particular value depends on
  // the values of its parents. There come in two varieties: switching 
  // parents and conditional parents. The values of the switching parents
  // determines an additional set of conditional parents. The variable
  // is conditioned on these conditional parents. Note that a switching 
  // parent may also be a conditional parent.
  // A node's switching and conditional parents must all be lower numbered 
  // than the node itself, thus allowing a topological ordering  
  // This sets the following arrays:
  // switchingParents, conditionalParentsList, allPossibleParents, and
  // (for the parents of the variable) it updates allPossibleChildren
  void setParents(vector<RandomVariable *> &sparents,
		  vector<vector<RandomVariable *> > &cpl);


  // Cached value of findConditionalParents(). Can reuse
  // this value w/o needing to do the integer map lookup
  // again.
  unsigned cachedIntFromSwitchingState;

  ////////////////////////////////////////////////////////////////////////
  // Looks at the currently (presumably clamped) values of the
  // switching parents and returns an 1D integer index. 
  // This indexes the set of conditional parents to use.
  unsigned intFromSwitchingState() 
  {
    // if dtMapper is NULL, we assume that
    // there are no switching parrents, which means
    // that there is only one set of conditional parents.
    // Return the index in that case.
    if (dtMapper == NULL)
      return 0;
    else
      return dtMapper->query(switchingParents);
  }

  ////////////////////////////////////////////////////////////////////////
  // Looks at the values of the switching parents, and sets the 
  // conditionalParents array pointer appropriately.
  virtual void findConditionalParents() = 0;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // SUPPORT ONLY FOR DISCRETE RANDOM VARIABLES ////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  /*
   * TODO: move this stuff to the child class.  
   */

  /////////////////////////////////////////////////////////////////
  // in the discrete case, the actual value of the variable.
  // Cliques keep track of the values of their discrete members.
  DiscreteVariableType val;
  
  /////////////////////////////////////////////////////////////////////////
  // Again in the discrete case.
  // The maximum number of possible values this RV may take in, irrespective
  // of the values of the parents. So, we must
  // have that 0 <= val < cardinality
  int cardinality;

  /////////////////////////////////////////////////
  // return true only if this node is discrete and
  // if *all* CPTs are MTCPT.
  virtual bool deterministic() { return false; }  


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // END OF SUPPORT ONLY FOR DISCRETE RANDOM VARIABLES /////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // Printing: show all the pertinent information about the node
  virtual void reveal(bool show_vals);
  // print the switching parent names/timeframes along with their currently clamped values.
  void printSwitchingParentInfo(FILE* stream);
  // print the current conditional parent name/times along with their currently clamped values.
  void printCurConditionalParentInfo(FILE* stream);

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Probability evaluation                            /////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // The inference algorithm guarantees that when this is called, the
  // variable and its parents will be clamped to appropriate values.
  // This also assumes that findConditionalParents has been called.
  virtual logpr probGivenParents() = 0;


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Parameter and Structure adjustment/checking       /////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  // Sets the parameters determining probGivenParents() to random values.
  virtual void makeRandom() = 0;

  ////////////////////////////////////////////////////////////////////////
  // Sets the parameters determining probGivenParents() to uniform values.
  virtual void makeUniform() = 0;

  ////////////////////////////////////////////////////////////////////////
  // Ties the parameters of 'this' with whatever those of 'other' are. 
  // 'other' and 'this' must be identical structuraly.
  virtual void tieParametersWith(RandomVariable*const other) = 0;

  ////////////////////////////////////////////////////////////////////////
  // Returns true of 'this' and 'other' are structurally identical.
  bool identicalStructureWith(const RandomVariable& other);

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Value iteration                                   /////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // Sets a variable to the first value that is possible.
  // Values with 0 probability can be ignored, as long as some value is set.
  // This function must always do a clamping.
  // Observation variables (already clamped do nothing).
  // It is assumed that the parent values are clamped at this point.
  virtual void clampFirstValue() = 0;

  ////////////////////////////////////////////////////////////////////////
  // Sets a variable to the next possible value.
  // Returns false when there are no more values or the variable is an
  // observation.
  // It is assumed that the parent values are clamped at this point.
  virtual bool clampNextValue() = 0;

  ////////////////////////////////////////////////////////////////////////
  // Sets the variable according to the probability distribution determined
  // by its parent's values. (i.e., sample from this distribution)
  // Used for simulation. 
  virtual void instantiate() = 0;


  ////////////////////////////////////////////////////////////////////////
  // It can be useful to tell a variable to cache the value it's currently
  // set to. For example, wheen keeping track of the best values seen
  // so far in enumerativeViterbi(). Calling this function tells a variable
  // to store its current value (wherever it wants).
  // Not expected to be called many times, so just leave it virtual.
  // This like a push operation, onto a one element size stack.
  virtual void cacheValue() = 0;

  ////////////////////////////////////////////////////////////////////////
  // Sets the variable's value to its cached value.
  // This like a pop operation, off of the stack.
  virtual void restoreCachedValue() = 0;

  /////////////////////////////////////////////////////////////////////////
  // stores a variable's value elsewhere
  virtual void storeValue(VariableValue &vv) = 0;

  /////////////////////////////////////////////////////////////////////////
  // sets a variables value as specified
  virtual void setValue(VariableValue &vv) = 0;

  /////////////////////////////////////////
  /////////////////////////////////////////
  // SUPPORT FOR EM  VARIABLES ////////////
  /////////////////////////////////////////
  /////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // On the backward pass, this will be called with a posterior probability
  // for the variable and its parents with all possible sets of values.
  // In inference, the GM is used as a global memory, with variable values
  // dynamically changed as inference proceeds. The inference algorithm
  // guarantees that the variable and its parents will have been set to
  // the the appropriate values when increment() is called, so that the
  // variable can just "look" at its parents to see what their values are.
  // Then the count of seeing the variable and its parents with those values
  // is incremented by the posterior amount.
  // For continuous variables, statistics are accumulated.
  virtual void emIncrement(logpr posterior) = 0;


  /////////////////////////////////////////
  /////////////////////////////////////////
  // END OF SUPPORT FOR EM  VARIABLES /////
  /////////////////////////////////////////
  /////////////////////////////////////////

 
  /////////////////////////////////////////
  /////////////////////////////////////////
  // SUPPORT FOR REPRODUCTION /////////////
  /////////////////////////////////////////
  /////////////////////////////////////////

  /////////////////////////////////////////  
  // to unroll the network, the variables in the repreating segments are 
  // cloned. Cloning must:
  // - copy the parent arrays (child arrays will be reset)
  // - copy the member data - cardinality, dtMapper, hidden, discrete
  virtual RandomVariable *clone();

  ////////////////////////////////////////////////////////
  // create an empty object of the same type as this (or
  // for whatever derived class this routine lives in).
  virtual RandomVariable *create() = 0;

};

#endif
