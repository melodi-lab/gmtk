/*-
 * GMTK_MixtureCommon.h
 *        Common elements that are shared by all Gaussiaan type clases.
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


#ifndef GMTK_MIXTURECOMMON_H
#define GMTK_MIXTURECOMMON_H

#include <map>
#include <set>

#include "fileParser.h"
#include "logp.h"
#include "machine-dependent.h"

#include "GMTK_NamedObject.h"
#include "GMTK_EMable.h"


class Dense1DPMF;
class MeanVector;
class DiagCovarVector;
class DlinkMatrix;
class Component;


class MixtureCommon : public EMable {

protected:

  ///////////////////////////////////////////////////////
  // dim = dimensionality of this Mixture
  // Typically, this will be the number of features over which
  // this mixture applies. All components of this
  // mixture must be of the same dimension of course.
  const unsigned _dim;

  ////////////////////////////////////////////////////////////////////
  // The number of mixture components for this mixture.
  unsigned numComponents;
  
public:

  /////////////////////////////////////////////////////////
  // support for component vanishing & splitting with multiple
  // potentially shared objects.
  static set<pair<Dense1DPMF*,unsigned> > vanishingComponentSet;
  static set<pair<Dense1DPMF*,unsigned> > splittingComponentSet;

  /////////////////////////////////////////////////////////
  // Support for cloning of portitions of components with shared objects.
  // These maps are used to remember when, in a given iteration,
  // a portion of a component has been cloned. If it already
  // has, rather than cloning again, we use the clone
  // that has already occured.
  static map<MeanVector*,MeanVector*> meanCloneMap;
  static map<DlinkMatrix*,DlinkMatrix*> dLinkMatCloneMap;
  static map<DiagCovarVector*,DiagCovarVector*> diagCovarCloneMap;
  static map<Component*,Component*> mcCloneMap;

  ///////////////////////////////////////////////////////////
  // This object must be aware of types of its derived
  // classes for reading from file purposes. 
  enum ContinuousImplementation { ci_unknown=-1,
				  ci_mixture=0,
				  ci_gausSwitchMixture=1,
				  ci_logitSwitchMixture=2,
				  ci_mlpSwitchMixture=3,
				  ci_zeroScoreMixture=4,
				  ci_unityScoreMixture=5
				  };
  const ContinuousImplementation mixType;

  MixtureCommon(const int dim, 
		     const ContinuousImplementation _mixType) 
    : _dim(dim),mixType(_mixType) {}
  virtual ~MixtureCommon() {}

  //////////////////////////////////////////////////////
  // Support for mixture vanishing
  static double mixCoeffVanishRatio;
  // The value that, if any mixture coef gets smaller than,
  // causes that component to vanish. I.e., if a mixture 
  // coeff is such that 
  // 
  //     p_component < 1/(curNumComponents*mixCoeffVanishRatio).
  // 
  // then it should vanish. Set to very large values (1e20) to
  // effectively turn off.

  // Support for mixture splitting, if a mixture coefficient prob.
  // gets greater or equal to mixCoeffSplitRatio/curNumComponents,
  // then the mixture will split. Typical values
  // might be something in the range 1.0 - 10.0. To
  // split single component, set it to unity. In other words, if
  // a mixture coeff is such that;:
  // 
  //     p_component >= mixCoeffSplitRatio/curNumComponents;
  // 
  // then it should split. Set to very large values (e.g., 1e10)
  // to effectively turn off.
  static double mixCoeffSplitRatio;
  // 
  // NOTE: must have the following for sanity:
  //  mixCoeffSplitRatio/curNumComponents > 
  //                         1/(curNumComponents*mixCoeffVanishRatio)
  static void checkForValidRatioValues();

  ///////////////////////////////////////////////////////////
  // set to true if we are supposed to cache the component probabilities
  // during EM training.
  static bool cacheComponentsInEmTraining;

  //////////////////////////////////////////////////////
  // force splitting of the number of top mixture componets
  // regardless of all else. Zero to turn off.
  static unsigned numTopToForceSplit;
  //////////////////////////////////////////////////////
  // Force vanishing of the number of bottom  mixture componets
  // regardless of all else. Zero to turn off.
  static unsigned numBottomToForceVanish;
  ////////////////////////////////////////////////////////////////
  // structure for comparing, needed for splitting/vanishing  
  struct LogpUnsignedPairCompare {  
    bool operator() (const pair<logpr,unsigned>& a, 
		     const pair<logpr,unsigned>& b) {
      return (a.first) < (b.first);
    }
  };

  /////////////////////////////////////
  // return the dimensionality
  unsigned dim() { return _dim; }

  //////////////////////////////////////////////
  // read/write basic parameters
  virtual void read(iDataStreamFile& is) = 0;
  virtual void write(oDataStreamFile& os) = 0;
  //////////////////////////////////////////////

  //////////////////////////////////
  // set all current parameters to valid but random values
  virtual void makeRandom() = 0;
  // set all current parameters to valid but "uniform" values 
  // (for Gaussians this means N(0,1))
  virtual void makeUniform() = 0;
  //////////////////////////////////


  //////////////////////////////////
  // probability evaluation
  virtual logpr log_p(const float *const x,    // real-valued scoring obs at time t
		      const Data32* const base, // ptr to base obs at time t
		      const int stride) = 0;       // stride
  virtual logpr log_p(const unsigned frameIndex, 
		      const unsigned firstFeatureElement) = 0;
  //////////////////////////////////

  //////////////////////////////////
  // EM Support                   //
  //////////////////////////////////
  virtual void emIncrement(logpr prob,
			   const unsigned frameIndex, 
			   const unsigned firstFeatureElement) = 0;

  //////////////////////////////////
  // Sample Generation            //
  //////////////////////////////////
  virtual void sampleGenerate(float *sample,
			      const Data32* const base,
			      const int stride) = 0;
  //////////////////////////////////

};


#endif
