
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

#ifndef GMTK_DISCRETERANDOMVARIABLE_H
#define GMTK_DISCRETERANDOMVARIABLE_H

#include <vector>

#include "GMTK_RandomVariable.h"
#include "GMTK_CPT.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_ObservationMatrix.h"

//////////////////////////////////////////////////////////////////////////
// A constant consisting of a special feature element value to indicate
// that this discrete RV should always use the value that it is currently
// assigned to (i.e., it is a fixed clamped value).
#define DRV_USE_FIXED_VALUE_FEATURE_ELEMENT (unsigned)(1 << (sizeof(unsigned)*8-1))

class DiscreteRandomVariable : public RandomVariable
{
private:
  friend class CPT;
  friend class MDCPT;
  friend class MSCPT;
  friend class MTCPT;
  friend class FileParser;

  //////////////////////////////////////////////////////////////////////
  // CPT array, one for each set of possible parents we might
  // have (so size of this array is the number of different
  // possible conditional parents).
  vector < CPT* > conditionalCPTs;

  ////////////////////////////////////////////////////////////////
  // The current CPT after findConditionalParents() is called.
  // It is "current" in the sence that it valid for the set
  // of parent values that are clamped. If the parent values
  // change, this CPT will no longer be valid until another
  // findConditionalParents() is called.
  CPT* curCPT;

  // iterator used between clamp functions.
  CPT::iterator it;

  // the feature file element corresponding to this RV.
  unsigned featureElement;

public:

  DiscreteRandomVariable(RVInfo&,
			 string _label, 
			 int card);


  ////////////////////////////////////////////////
  // Assuming the parents have been allocated, this forces
  // the internal CPT structures to 1) be allocated and 
  // 2) match the cardinalities of the parents.
  void allocateProbabiltyTables();

  void setCpts(vector<CPT*> &cpts);

  // returns true if all cpts are determinstic.
  bool deterministic() {
    // first check if we have any cpts. If not
    // we be conservative and return false (this
    // case might arrise when doing elimination
    // experiments without needing any parameters).
    // if (conditionalCPTs.size() == 0)
    // return false;

    // TODO: make this a member of the rv_info (so that we don't
    // re-compute this over and over).
    for (unsigned i=0;i<rv_info.discImplementations.size();i++) {
      if (rv_info.discImplementations[i] != CPT::di_MTCPT)
	return false; 
    }
    //    for (unsigned i=0;i<conditionalCPTs.size();i++) {
    //      if (conditionalCPTs[i]->cptType != CPT::di_MTCPT)
    //	return false; 
    //    }
    return true;
  }

  // returns true if all cpts are sparse or determinstic.
  bool sparse() {
    // TODO: make this a member of the rv_info (so that we don't
    // re-compute this over and over).
    for (unsigned i=0;i<rv_info.discImplementations.size();i++) {
      if ((rv_info.discImplementations[i] != CPT::di_MTCPT)
	  &&
	  (rv_info.discImplementations[i] != CPT::di_MSCPT))
	return false; 
    }
    return true;
  }



  // This routine returns the average cardinality (average number of
  // possible child random variable values with non-zero probability)
  // for the current random variable.  This routines assume that the
  // RV's parents are in the same clique (or an earlier clique
  // relative to the JT root) so that we only consider the child
  // values with non-zero probability.
  unsigned averageCardinality() {
    unsigned res;
    for (unsigned i=0;i<rv_info.discImplementations.size();i++) {
      if (rv_info.discImplementations[i] == CPT::di_MDCPT)
	res += cardinality;
      else if (rv_info.discImplementations[i] == CPT::di_MTCPT)
	res += 1;
      else { // sparce CPT
	// we only get the sparse CPT cardinality if it has
	// been read in, otherwise, we just take card.
	unsigned tmp;
	if (i < conditionalCPTs.size())
	  tmp = conditionalCPTs[i]->averageCardinality();
	else 
	  tmp = cardinality;
	res += tmp;
      }
    }
    if (res == 0) {
      // if there are no CPTs (i.e., if we didn't load in the
      // parameters during say a triangulation procedure, then we just
      // return the cardinality.
      return cardinality;
    } else {
      // otherwise, we return what we've got.
      return res/rv_info.discImplementations.size();
    }
  }

  // This routine returns the maximum possible cardinality (max number
  // of possible child random variable values with non-zero
  // probability) for the current random variable.  This routines
  // assume that the RV's parents are in the same clique (or an
  // earlier clique relative to the JT root) so that we only consider
  // the child values with non-zero probability.
  unsigned maxCardinality() {
    unsigned res = 0;
    for (unsigned i=0;i<rv_info.discImplementations.size();i++) {
      if (rv_info.discImplementations[i] == CPT::di_MDCPT) {
	res = cardinality;
	// it'll never get bigger than this so we 
	// might as well leave the loop
	break;
      } else if (rv_info.discImplementations[i] == CPT::di_MTCPT) {
	if (res == 0) res = 1;
      } else { // sparce CPT
	// we only get the sparse CPT cardinality if it has
	// been read in, otherwise, we just take card.
	unsigned tmp;
	if (i < conditionalCPTs.size())
	  tmp = conditionalCPTs[i]->maxCardinality();
	else 
	  tmp = cardinality;
	if (tmp > res) res = tmp;
      }
    }
    if (res == 0) {
      // if there are no CPTs (i.e., if we didn't load in the
      // parameters during say a triangulation procedure, then we just
      // return the cardinality.
      return cardinality;
    } else {
      // otherwise, we return the maximum.
      return res;
    }
  }

  // various other routines want a cardinality to 'use'.  they use the
  // routine useCardinality() which can be defined either as the avg
  // or max cardinality above.
  unsigned useCardinality() { return maxCardinality(); }  


  ////////////////////////////////////////////////////////////////
  // Set up conditional parents pointers and other tables.
  void findConditionalParents();

  ////////////////////////////////////////////////////////////////
  // These next several routines use the switching state and
  // parent set determined when clampFirstValue() is called.
  // clampFirstValue() invokes findConditionalParents(), and  
  // calls to clampNextValue() and probGivenParents() rely on being
  // called in the appropriate context. i.e. the values of the switching
  // parents must not have changed since the last call to clampFirstValue().
  // This imposes some constraints -- which do happen to be satisfied by
  // the inference loops.
  // 
  // compute the probability
  logpr probGivenParents() {
    logpr _cachedProb = curCPT->probGivenParents(*curConditionalParents,val);
    // 
    // TODO: make this weight stuff part of a subclass with a new version
    // of this function, same speed since probGivenParents is virtual anyway.
    if (wtStatus != wt_NoWeight) {
      if (wtStatus == wt_Constant)
	_cachedProb.valref() *= wtWeight;
      else // get weight from observation matrix at current time frame
	_cachedProb.valref() *= 
	  (*globalObservationMatrix.floatVecAtFrame((unsigned)timeIndex, wtFeatureElement));
    }
    return _cachedProb;
  }

  logpr probGivenParentsWSetup() {
    findConditionalParents();
    // this isn't necessary.
    // curCPT->becomeAwareOfParentValues(*curConditionalParents);
    return probGivenParents();
  }


  // set the RV (which must be 
  void setToObservedValue() {
    // assert (!hidden);
    // observed, so set value from observation matrix
    if (globalObservationMatrix.active() && featureElement != DRV_USE_FIXED_VALUE_FEATURE_ELEMENT) {
      // printf("getting value of random variable '%s', time index %d, el %d\n",
      // label.c_str(),timeIndex,featureElement);
      unsigned tmp = globalObservationMatrix.unsignedAtFrame(timeIndex,featureElement);

      // TODO: do this check outside of this loop.
      if (tmp >= (unsigned)cardinality) 
	error("ERROR: RV '%s' at time index %d has cardinality %d, but feature element position %d in observation file (time %d of segment %d) has value %u.\n",
	      label.c_str(),
	      timeIndex,
	      cardinality,
	      featureElement,
	      timeIndex,
	      globalObservationMatrix.segmentNumber(),
	      tmp);
      val = tmp;
    }
    // otherwise, we keep the value set to what it was before.
    return;
  }




  // clamp this RV to its "first" value
  void clampFirstValue() {
    findConditionalParents();
    if (!hidden) {
      setToObservedValue();
      return;
    }
    // a hidden variable, so we set up the iterator.
    curCPT->becomeAwareOfParentValues(*curConditionalParents);
    it = curCPT->begin(); 
    val = it.val();
    // assert ( val >= 0 );
  }

  // continue on
  bool clampNextValue() { 
    if (!hidden) 
      return false;
    it++; 
    if (it!=curCPT->end()) 
      val = it.val();
    // assert ( val >= 0 );
    return (it != curCPT->end()); 
  }


  // clamp this RV to its "first" value
  void begin() {
    findConditionalParents();
    if (!hidden) {
      // must be done externally
      // setToObservedValue();
      return;
    }
    // a hidden variable, so we set up the iterator.
    curCPT->becomeAwareOfParentValues(*curConditionalParents);
    it = curCPT->begin(); 
    val = it.val();
    // assert ( val >= 0 );
  }

  // continue on
  bool next() { 
    if (!hidden) 
      return false;
    it++; 
    if (it!=curCPT->end()) 
      val = it.val();
    // assert ( val >= 0 );
    return (it != curCPT->end()); 
  }





  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // Value caching support.
  RandomVariable::DiscreteVariableType cached_val;
  void cacheValue() {cached_val=val;}
  void restoreCachedValue() {val=cached_val;}

  /////////////////////////////////////////////////////////////////////////
  // stores a variable's value elsewhere
  void storeValue(VariableValue &vv) {vv.ival = val;}

  /////////////////////////////////////////////////////////////////////////
  // sets a variables value as specified
  void setValue(VariableValue &vv) {val = vv.ival;}

  void makeRandom() { 
    for (unsigned i=0;i<conditionalCPTs.size();i++) 
      conditionalCPTs[i]->makeRandom();
  }

  void makeUniform()  { 
    for (unsigned i=0;i<conditionalCPTs.size();i++) 
      conditionalCPTs[i]->makeUniform();
  }

  void tieParametersWith(RandomVariable*const other,
			 bool checkStructure=true);

  ////////////////////////////////////////////////////////////////
  // Sample, set value.
  void instantiate() { 
    findConditionalParents(); 
    curCPT->becomeAwareOfParentValues(*curConditionalParents);
    val = curCPT->randomSample(); 
  }

  ///////////////////////////////////////////////////
  // EM Support
  void emIncrement(logpr posterior) { 
    findConditionalParents();
    curCPT->emIncrement(posterior,this);
  }

  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  // reproduction routines.

  RandomVariable *create() { 
    return new DiscreteRandomVariable(rv_info,label,cardinality);
  }
  RandomVariable *clone();
  RandomVariable *cloneWithoutParents();

  ///////////////////////////////////////////////////

};

#endif
