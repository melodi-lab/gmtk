/*-
 * GMTK_DlinkMatrix.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */


#ifndef GMTK_DLINKMATRIX_H
#define GMTK_DLINKMATRIX_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"
#include "GMTK_Dlinks.h"

class LinMeanCondDiagGaussian;
class MeanVector;
class DiagCovarVector;

class DlinkMatrix : public EMable  {

  friend class Dlinks;
  friend class LinMeanCondDiagGaussian;
  friend class MeanVector;
  friend class DiagCovarVector;

  ///////////////////////////////////////////////////////
  // The actual dlink structure
  Dlinks* dLinks;

  ///////////////////////////////////////////////////////////
  // The acutal matrix data values, packed
  // into one 1D array
  sArray< float > arr;
  ///////////////////////////////////////////////////////////  

  ///////////////////////////////////////////////////////////
  // Data structures support for EM
  sArray< float > nextArr;

  /////////////////////////////////////////////////////////
  // the denominator used in the case of shared dlinks
  sArray<double> sharedZZDenominator;

  /////////////////////////////////////////////////
  // counts the number of gaussian components
  // that are sharing this mean at all. This is a static
  // count, and is computed as any object that
  // might use a mean (such as a Gaussian component)
  // is read in.
  unsigned numTimesShared;

  /////////////////////////////////////////////////
  // counts the number of gaussian components
  // that are sharing this mean at EM training time. This is a dynamic
  // count, and is computed as EM training is run. This
  // value does not necessarily equal the number of
  // objects that have specified this object in
  // the object files.
  unsigned refCount;

public:

  ///////////////////////////////////////////////////////////
  // General constructor
  DlinkMatrix();
  ~DlinkMatrix() { }


  // When noisy cloning an object, this gives
  // the fraction to multiply to get the STD of the noise.
  static double cloneSTDfrac;
  static void checkForValidValues();

  //////////////////////////////////
  // set all current parameters to random/uniform values
  void makeRandom();
  void makeUniform();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  // create a copy of self, but with slightly perturbed
  // means values.
  DlinkMatrix* noisyClone();

  /////////////////////////////////////////////////
  // create a copy of self, with entirely new parameters with
  // identical values; NOTHING is shared
  DlinkMatrix* identicalIndependentClone();

  ///////////////////////////////////
  unsigned totalNumberParameters() { return arr.len(); }
  void recursivelyClearUsedBit() {  emClearUsedBit();  }
  void recursivelySetUsedBit() { emSetUsedBit();  }

  ///////////////////////////////////////////////////////////  
  // num number of features (dimensionality) for this 
  int dim() const { return dLinks->dim(); }
  unsigned totalNumberLinks() { return dLinks->totalNumberLinks(); }
  unsigned zzAccumulatorLength() { return dLinks->zzAccumulatorLength(); }

  ///////////////////////////////////////////////////////////  
  // numLinks: return the number of links for the ith
  // feature.
  int numLinks(const int i) const { 
    assert ( i >=0 && i < dim() );
    return dLinks->numLinks(i);
  }

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration(sArray<float>& xzAccumulators,
			sArray<float>& zzAccumulators,
			sArray<float>& zAccumulators);
  void emIncrement(const logpr prob,
		   const float fprob,
		   const float* const f,
		   const Data32* const base,
		   const int stride,
		   float* xzAccumulators,
		   float* zzAccumulators,
		   float* zAccumulators);
  void emIncrement(const logpr prob,
		   const float fprob,
		   const float* const f,
		   const Data32* const base,
		   const int stride);
  void emEndIterationNoSharingAlreadyNormalized(const float*const xzAccumulators);
  void emEndIterationSharedMeansCovarsDlinks(const float*const xzAccumulators,
					     const float*const zAccumulators,
					     const float*const zzAccumulators,
					     const MeanVector* mean,
					     const DiagCovarVector* covar);

  void emSwapCurAndNew();

  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false) {};
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {};
  void emZeroOutObjectsAccumulators() {};
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {};
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {};
  const string typeName() { return "dlink matrix"; }
  // need to override parent class's routine in this case.
  void emStoreAccumulators(oDataStreamFile& ofile);
  //////////////////////////////////


};



#endif
