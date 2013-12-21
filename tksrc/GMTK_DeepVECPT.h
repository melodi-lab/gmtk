/*-
 * GMTK_DeepVECPT.h
 *
 *   A "deep virtual evidence" CPT class. This CPT gives GMTK the ability
 *   to apply deep models as virtual evidence, where the probabilties come
 *   from a deep model applied to a window of observations.
 *
 *   Specifically, this CPT can be used for a binary variable C that
 *   has exactly one parent A.  The CPT uses the output of the deep model
 *   as the probabilities P(C=1|A=a). 
 *   Note that C need not be observed, although it typically is to provide
 *   virtual online evidence to variable A.  In this way, the actual
 *   deep model outputs need not even be probabilties, they can be
 *   arbitrary scores. If C is hidden, however, then the values of
 *   P(C=0|A=a) are taken to be 1-P(C=1|A=a), so in such a case (when
 *   C is hidden) it is more sensible for the scores for each a to be
 *   actual values between 0 and 1.
 *
 *  Written by Richard Rogers <rprogers@ee.washington.edu>
 * 
 * 
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_DEEPVECPT_H
#define GMTK_DEEPVECPT_H

#include "fileParser.h"
#include "logp.h"

#include "sArray.h"

#include "GMTK_CPT.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRVType.h"
#include "GMTK_NamedObject.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_FileSource.h"
#include "GMTK_DeepNN.h"

// we need to interface to the external global observation
// matrix object to get some parametes (such as start skip, end skip,
// and the length of each utterance to ensure lengths are the same).
extern ObservationSource *globalObservationMatrix;

class DeepVECPT : public CPT {

  ObservationSource *obs;

  // TODO: redo all this so that it works well
  // with obs.openFile(). Either change to char* and make
  // proper destructor, or something else.

  // use frames t-window_radius to t+window_radius as input
  // to the deep model for frame t
  unsigned window_radius; 

  // number of floats per frame taken from the file
  unsigned nfs;

  // these variables are used in case that we
  // get the VECPT information from the global
  // observation matrix directly, rather than
  // from a local observation file here. Otherwise,
  // this variable is set to zero. We have
  // two offsets, one for floats and one for ints.
  unsigned obs_file_foffset;

  ////////////////
  // current parent value
  unsigned curParentValue;

  ////////////////
  // the value
  int _val;

  ////////////////
  // the Deep Multi-Layer Perceptron
  DeepNN *dmlp;

  ////////////////
  // remember the computed CPT so we don't have to recompute it
  unsigned  cached_segment; 
  unsigned  cached_frame;
  double   *cached_CPT;

  float    *input_vector;

  // Apply the deep neural network to get the probability
  logpr applyNN(DiscRVType parentValue, DiscRV * drv);

public:

  ///////////////////////////////////////////////////////////  
  // General constructor, 
  // VECPTs always have one parent, and a binary child.
  DeepVECPT() : CPT(di_DeepVECPT), dmlp(NULL),
    cached_segment(0xFFFFFFFF), cached_frame(0xFFFFFFFF), cached_CPT(NULL), input_vector(NULL)
  { 
    _numParents = 1; _card = 2; cardinalities.resize(_numParents); 
  }

  ~DeepVECPT() { 
    if (cached_CPT) delete[] cached_CPT; 
    if (input_vector) delete[] input_vector;
  }


  DeepNN *getDeepNN() { return dmlp; }

  // a VECPT is considered iterable since its implementation can
  // change not only from segment to segment, but even within a
  // segment.
  bool iterable() { return true; }


  // Number of input features taken from a single obs file frame
  unsigned numFeaturesPerFrame() { return nfs; }


  // Number of past/future frames included in the input
  unsigned windowRadius() { return window_radius; }


  unsigned obsOffset() { return obs_file_foffset; }
  

  ///////////////////////////////////////////////////////////    
  // Semi-constructors: useful for debugging.
  // See parent class for further documention.
  void setNumParents(const int _nParents) {
    assert ( _nParents == 1 );
    CPT::setNumParents(_nParents);
    bitmask &= ~bm_basicAllocated;
  }
  void setNumCardinality(const int var, const int card) {
    if (var == 1) {
      // setting child card.
      assert ( card == 2);
      CPT::setNumCardinality(var,card);
    } else if (var == 0) {
      // setting parent cardinality, which can be anything, so no check.
      CPT::setNumCardinality(var,card);
    } else {
      // we're never allowed to have more than 1 parent in this case.
      assert (0);
    }
  }
  void allocateBasicInternalStructures() {
    error("DeepVECPT::allocateBasicInternalStructures() not implemented");
  }

  ///////////////////////////////////////////////////////////  
  // Probability evaluation, compute Pr( child | parents ), and
  // iterator support. See GMTK_CPT.h for documentation.
  void becomeAwareOfParentValues( vector <RV *>& parents,
				  const RV* rv );
  void begin(iterator& it,DiscRV* drv, logpr& p);
  void becomeAwareOfParentValuesAndIterBegin(vector < RV *>& parents, 
					     iterator &it, 
					     DiscRV* drv,
					     logpr& p);
  logpr probGivenParents(vector <RV *>& parents,
			 DiscRV * drv);
  bool next(iterator &it,logpr& p);

  ///////////////////////////////////
  // random sample given current parent value.
  int randomSample(DiscRV* drv);

  ////////////////////////////////////////////////////////////////////
  // How many parameters does this consume? We return 0 since the
  // virtual evidence does not constitute parameters in the normal
  // sense.
  unsigned totalNumberParameters() { return 0; }

  // The above is true for EM training, but gmtkDMLPtrain actually
  // learns the deep VE CPT parameters, so count them
  unsigned totalNumberDMLPParameters() {
    return dmlp->totalNumberParameters();
  }

  ///////////////////////////////////////////////////////////  
  // These routines are no-ops in this case since all
  // distributions come from files.
  void normalize() {}
  // set all values to random values.
  void makeRandom() {}
  // set all values to uniform values.
  void makeUniform() {}

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration() {}
  void emIncrement(logpr p,vector <RV*>& parents,RV*rv) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}

  // parallel training (DeepVECPTs are never EM trained)
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false) {}
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  void emZeroOutObjectsAccumulators() {}
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  const string typeName() { return "DeepVirtualEvidenceCPT"; }
  //////////////////////////////////

  // support to change the segment number
  void setSegment(const unsigned segNo) {
    obs->openSegment(segNo);
  }
  unsigned numFrames() {
    // return the number of frames in the current segment.
    return obs->numFrames();
  }

};



#endif // defined DeepVECPT
