/*-
 * GMTK_EMable
 *        Objects that inherit from this are "EMable" in the sense
 *        that we can use them in EM. This class provides some support
 *        for such objects.
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


#ifndef GMTK_EMABLE_H
#define GMTK_EMABLE_H

#include "fileParser.h"
#include "logp.h"
#include "machine-dependent.h"
#include "sArray.h"
#include "debug.h"

#include "GMTK_NamedObject.h"

class RV;

class EMable : public NamedObject, public IM {

protected:

  //////////////////////////////////////////////////////////////////
  // A bitmask giving the "state" of the object. Used for
  //   1) error checking
  //   2) swapping parameters when EMables are shared (so
  //      that parameters aren't swapped twice)
  //   3) Loading/Storing accumulators, to make sure they aren't
  //      loaded/stored more than one time.
  //   and a number of other things (see below).
  enum {
    // Basic data structures allocated (after a read)
    bm_basicAllocated  = (1 << 0),
    // Are the em structures allocated
    bm_emAllocated     = (1 << 1),
    // True when an EM epoch is ongoing, so that we don't call
    // startEmEpoch multiple times when things are tied together.
    bm_emEpochOnGoing  = (1 << 2),
    // May the current and next parameters been swapped.
    bm_swappable         = (1 << 3),
    // Has the accumulators been loaded/stored.
    bm_accLoadStore    = (1 << 4),
    // Are we training this object.
    bm_amTraining    = (1 << 5),
    // True if  the parameters of this object shared by
    // multiple users. NOTE: this might not be used
    // by all child classes of EMable, used only by
    // those objects who keep track of if they are being
    // shared.
    bm_isShared    = (1 << 6),
    // bm_isUsed == is this object "used" by anyone, in that
    // is a mean used by a component, or is a component used
    // by a mixture. Again, this field might not be needed
    // by all EM objects.
    bm_isUsed    = (1 << 7),
    // some objects will initialize their accumlators at
    // the beginning of an EM epoch, while others will
    // initialize them at the end (since the real accumulators
    // live in other objects). This bit is meant for those
    // objects who initializet their accumulators at the end.
    bm_accInitialized    = (1 << 8),
    // some emable objects, such as a mixture of gaussians,
    // might want to cache the computation of probabilities
    // after it evaluates them for inference, so that
    // they can be used again for learning. This next bit
    // is for those such objects. If the bit is set,
    // the probabilities will be cached. 
    bm_cacheProbs    = (1 << 9),

    // Initial State, where no data structures have been allocated.
    bm_initState      = 0x0
  };

  unsigned int bitmask;

  ////////////////////////////////////////////////
  // the sum of the EM accumulated probability for this
  // object.
  logpr accumulatedProbability;

  ////////////////////////////////////////////////
  // If the posterior probability falls below this
  // value, we do not increment.
  static logpr minIncrementProbabilty;
  static unsigned long missedIncrementCount;


  ///////////////////////////////////////////////////////
  // The minimum continuous accumulated probability of mean and covariance -like
  // objects. If the accumulated probability falls below this
  // value, then the mean or variance like object will not
  // update its values.
  static logpr _minContAccumulatedProbability;
  // same thing as above but for discrete objects.
  static logpr _minDiscAccumulatedProbability;

 
public:

  // global variable to determine if this process includes priors in
  // the counts. We keep this here since when using parallelism, we
  // only want to add the priors to one process, not all of them.
  static bool useDirichletPriors;


  // global variable to determine if, when we do the accumulators, we
  // actually want to accumulate according to the Fisher kernel
  // scores.  I.e., rather than actually accumulating the EM expected
  // counts, what we do is accumulate so that the result is the
  // derivative of the log likelihood.
  static bool accumulateFisherKernelScores;

  
  EMable() { bitmask = bm_amTraining; }
  virtual ~EMable() {}

  static logpr minContAccumulatedProbability()
  { return _minContAccumulatedProbability; }
  static logpr minDiscAccumulatedProbability()
  { return _minContAccumulatedProbability; }
  static logpr setMinContAccumulatedProbability(const logpr floor);
  static logpr setMinDiscAccumulatedProbability(const logpr floor);

  inline logpr get_accumulatedProbability() { return accumulatedProbability; }

  // returns the type of the sub-object in string
  // form that is suitable for printing and identifying
  // the type of the object.
  virtual const string typeName() = 0;

  /////////////////////////////////////////////////////////////////
  // clear the swap bit, needed for sharing.
  void clearBasicAllocatedBit() { bitmask &= ~bm_basicAllocated; }
  void setBasicAllocatedBit() { bitmask |= bm_basicAllocated; }
  bool basicAllocatedBitIsSet() { return (bitmask & bm_basicAllocated); }


  //////////////////////////////////
  // EM training                  //
  //////////////////////////////////


  /*
   *
   ***************************************************************
   * The following dummy routines describe the
   * generic EM Loop functionality for an EMable.
   * An "iteration" is an entire pass through the data.
   * There might, however, be variations to this theme
   * (for example, for mean- and variance- like objects)
   * to support either parameter tying, or to improve
   * computational efficiency.
   *
   *emStartIeration() {
   *  // First check if this type object (means, variances,
   *  // dlink matrices, etc) is currently set to be trained
   *  // or is it fixed for this EM.
   *  if (!object of this type being trained with EM)
   *    return;
   *
   *  // Next check if EM has been started for this object.
   *  // If it has, then we don't do anything.
   *  if(emOnGoingBitIsSet())
   *    return; // already done
   *  // Next check to make sure any internal EM structures
   *  // have been allocated, using the allocated bit. If
   *  // not, allocate them and set the bit.
   *  if (!emEmAllocatedBitIsSet()) {
   *    // Allocate interal structurse
   *    // ... 
   *    // and set the bit
   *    emSetEmAllocatedBit();
   *  }
   *  // Next, make it such that we set the ongoing bit.
   *  emSetOnGoingBit();
   *  // Next, set the swappable bit so that this object's parameters
   *  // will get swapped when an emEndIteration is called later.
   *  emSetSwappableBit();
   *
   *  // Finally, do any object specific initialization of
   *  // the next parameters, such as setting to zero (or
   *  // to prior values), and so on, and calling
   *  // any member objects emStartIteration();
   *  // ...
   *}
   *
   *
   *emIncrement(logpr prob, ...) {
   *  // First check if this type object (means, variances,
   *  // dlink matrices, etc) is currently set to be trained
   *  // or is it fixed for this EM.
   *  if (!object of this type being trained with EM)
   *    return;
   *
   *  // Next, call emStartIteration. The idea here
   *  // is that we do not know, during each EM iteration,
   *  // what objects are going to be actively trained using EM.
   *  // We just start incrementing and in a lazy fashion
   *  // start the EM iterations going. If an iteration is
   *  // already going, this ends up being a nop, otherwise
   *  // things are initialized appropriately.
   *
   *  emStartIteration();
   *
   *  // Next, perhaps check that this prob is large enough and if
   *  // it is not we don't bother to increment anything.
   *  if (prob < minIncrementProbabilty) {
   *     missedIncrementCount++;
   *     return;
   *  }
   *
   *  // Next, accumulate the probabilities.
   *  accumulatedProbability+= prob;
   *
   *  // Finally, do object specific incrementing which depends
   *  // also on the remaining arguments after prob.
   *  // ..
   *
   *}
   *
   *emEndIteration() {
   *  // First check if this type object (means, variances,
   *  // dlink matrices, etc) is currently set to be trained
   *  // or is it fixed for this EM.
   *  if (!object of this type being trained with EM)
   *    return;
   *
   *  // Next check if EM is on going for this object, and if it
   *  // is not, return making this a nop. This way, endIteration
   *  // may be called from a collection of global EMable objects
   *  // w/o respect to who is using them via parameter tying.
   *
   *  if ( !emOnGoingBitIsSet() )
   *    return; 
   *
   *  // Next, check if the accumulator is zero, if it is
   *  // then this means that emIncrement() either was not
   *  // called or it was called with very low or zero values
   *  // for prob. Issue a warning for now.
   *
   *  if (accumulatedProbability.zero()) {
   *    warning("...");
   *    // TODO: also check if this when inverted will explode.
   *  }
   *
   *  // Next, do any obect specific ending, which might include
   *  // ending for sub-objects.
   *
   *  // Finaly, stop the EM iteration. That way emEndIteration
   *  // may be called multiple times (via shared objects for 
   *  // parameter tieing) and anything other than the first
   *  // time will be a nop.
   *  emClearOnGoingBit();
   *
   *}
   *
   *
   *emSwapCurAndNew() {
   *  // First check if this type object (means, variances,
   *  // dlink matrices, etc) is currently set to be trained
   *  // or is it fixed for this EM.
   *  if (!object of this type being trained with EM)
   *    return;
   *
   *
   *  // Next, check if this object is swappable and if
   *  // not this is a nop.
   *  if (!emSwappableBitIsSet())
   *    return;
   *
   *  // Next, do the stuff to swap the current and next parameters,
   *  // possibly calling this routine on member variables.
   *
   *  // Finally, reset the swappable bit so that this is only
   *  // called once.
   *  emClearSwappableBit();
   *}
   *
   *
   *The global EM loop should look something like:
   *
   *  for each training sample
   *      for each object and relevant probability
   *           ob.emIncrement(prob)
   *  GlobalCollection.emEndIteration();
   *  if converged
   *      GlobalCollection.emSwapCurAndNew();     
   *
   */
  
  

  ////////////////////////////////////////////////////////////////////
  // begins a new epoch. Also ensures that data for EM is allocated
  // and possibly changes the alocated bit, and ongoing bit
  virtual void emStartIteration()  { assert(0); }
  // The next few are for Gaussian component objects.
  virtual void emStartIteration(sArray<float>& componentsNextMeans)  { assert(0); }
  // The next few are for Gaussian component objects.
  virtual void emStartIteration(sArray<float>& xzAccumulators,
				sArray<float>& zzAccumulators,
				sArray<float>& zAccumulators)
  { assert(0); }

  ////////////////////////////////////////////////////////////////////
  // Accumulate new data into the internal structures for EM.
  // assumes that the ongoing bit it set. The various 
  // versions are for differnet objects that need to utilize
  // different information when they are being called.
  // E.g., a discrete random variable that uses a CPT, the CPT needs to know 
  // the current values of the parents. A continuous random variable
  // with a mean object, the mean object needs to know the current values
  // of the observations, say. These should never be called directly,
  // and should be re-defined in child classes. They are not abstract
  // though (i.e., "=0") because we want the option to need not define them
  // for a particular child class.
  virtual void emIncrement(logpr prob) { assert(0); }
  virtual void emIncrement(logpr prob,RV*r) { assert(0); }
  virtual void emIncrement(logpr prob,vector < RV* >& parents, RV*r) { assert(0); }
  // for real-valued things.
  virtual void emIncrement(const logpr prob,
			   const float*f,
			   const Data32* const base,
			   const int stride)
                        { assert(0); }
  virtual void emIncrement(const logpr prob,
			   const float fprob, // precomputed prob.unlog() for speed
			   const float*f,
			   const Data32* const base,
			   const int stride)
                        { assert(0); }
  virtual void emIncrement(logpr prob,
			   const unsigned frameIndex, 
			   const unsigned firstFeatureElement) 
                        { assert(0); }
  virtual void emIncrement(logpr prob,sArray<logpr>& a) { assert(0); }  
  virtual void emIncrement(logpr prob,const int val) { assert(0); }  


  ////////////////////////////////////////////////////////////////////
  // Accumulate new data into the internal structures for this 
  // em iteration, clears the ongoing bit.
  virtual void emEndIteration() { assert(0); }
  // This second two is for mean/covariance objects who need the
  // corresponding accumulated stats from its shared components.
  // virtual void emEndIteration(const float*const m) { assert(0); }
  // virtual void emEndIteration(const logpr prob,const float*const m,const float*const v) { assert(0); }

  ////////////////////////////////////////////////////////////////////////////
  // if swap bit not set, swaps the current and new parameters, set swap bit.
  // otherwise does nothing.
  virtual void emSwapCurAndNew() = 0;

  /////////////////////////////////////////////////////////////////
  // clear the swap bit, needed for sharing.
  void emClearSwappableBit() { bitmask &= ~bm_swappable; }
  void emSetSwappableBit() { bitmask |= bm_swappable; }
  bool emSwappableBitIsSet() { return (bitmask & bm_swappable); }

  /////////////////////////////////////////////////////////////////
  // basic allocated bit
  void emClearEmAllocatedBit() { bitmask &= ~bm_emAllocated; }
  void emSetEmAllocatedBit() { bitmask |= bm_emAllocated; }
  bool emEmAllocatedBitIsSet() { return (bitmask & bm_emAllocated); }

  /////////////////////////////////////////////////////////////////
  // clear the swap bit, needed for sharing.
  void emClearOnGoingBit() { bitmask &= ~bm_emEpochOnGoing; }
  void emSetOnGoingBit() { bitmask |= bm_emEpochOnGoing; }
  bool emOnGoingBitIsSet() { return (bitmask & bm_emEpochOnGoing); }

  /////////////////////////////////////////////////////////////////
  void emClearAmTrainingBit() { bitmask &= ~bm_amTraining; }
  void emSetAmTrainingBit() { bitmask |= bm_amTraining; }
  bool emAmTrainingBitIsSet() { return (bitmask & bm_amTraining); }

  /////////////////////////////////////////////////////////////////
  void emClearSharedBit() { bitmask &= ~bm_isShared; }
  void emSetSharedBit() { bitmask |= bm_isShared; }
  bool emSharedBitIsSet() { return (bitmask & bm_isShared); }

  /////////////////////////////////////////////////////////////////
  void emClearUsedBit() { bitmask &= ~bm_isUsed; }
  void emSetUsedBit() { bitmask |= bm_isUsed; }
  bool emUsedBitIsSet() { return (bitmask & bm_isUsed); }

  /////////////////////////////////////////////////////////////////
  void emClearAccInitializedBit() { bitmask &= ~bm_accInitialized; }
  void emSetAccInitializedBit() { bitmask |= bm_accInitialized; }
  bool emAccInitializedBitIsSet() { return (bitmask & bm_accInitialized); }

  /////////////////////////////////////////////////////////////////
  void emClearCacheProbsBit() { bitmask &= ~bm_cacheProbs; }
  void emSetCacheProbsBit() { bitmask |= bm_cacheProbs; }
  bool emCacheProbsBitIsSet() { return (bitmask & bm_cacheProbs); }


  // return the number of parameters for object.
  virtual unsigned totalNumberParameters() = 0;

  /////////////////////////////////////////////////////////////
  // For parallel EM training, a number of routines
  // to support the loading and storing of objects accumulators.
  // These routines are a bit tricky only because they
  // need to work with 
  //    1) parameter sharing
  //    2) objects which are not trained for an iteration
  //    3) parallelism.

  ///////////////////////////////////////////////////////////////
  // store the current set of accumulators to a file.
  virtual void emStoreAccumulators(oDataStreamFile& ofile);

  ///////////////////////////////////////////////////////////////
  // load the current set of accumulators from a file.
  virtual void emLoadAccumulators(iDataStreamFile& ifile);

  //////////////////////////////////////////////////////////////////////
  // accumulate (add to) the current set of accumulators from a file.
  virtual void emAccumulateAccumulators(iDataStreamFile& ifile);

  ///////////////////////////////////////////////////////////////
  // initialize (say to zero, or to appropriate starting values) the current set of accumulators.
  virtual void emInitAccumulators();

  ///////////////////////////////////////////////////////////////
  // write accumulators in a simple, easy, unencoded, and fixed
  // length way to read by humans, or to be used as a input to feature
  // transform for a kernel machine.
  virtual void emWriteUnencodedAccumulators(oDataStreamFile& ofile,
					    bool writeLogVals = true);


  ///////////////////////////////////////////////////////////////
  // virtual functions for objects to do the actual work.
  virtual void emStoreObjectsAccumulators(oDataStreamFile& ofile,  // file to which to write
					  bool writeLogVals = true, // write log(prob) vals
					  bool writeZeros = false // write just a vector of zeros of appropriate length
					  ) = 0;


  virtual void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) = 0;
  virtual void emZeroOutObjectsAccumulators() = 0;
  virtual void emLoadObjectsAccumulators(iDataStreamFile& ifile) = 0;
  virtual void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) = 0;

  //////////////////////////////////////////////



};

#endif

