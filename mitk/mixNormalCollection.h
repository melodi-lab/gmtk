/**
 *: mixnormalcollection.h
 *
 * @author Gang Ji
 * @author gang@ee.washington.edu
 * @version 1.0
 * $Id$
 */


#ifndef MIX_NORMAL_COLLECTION_H
#define MIX_NORMAL_COLLECTION_H

#include <cstdio>


#include "mixNormal.h"
#include "range.h"
#include "readRange.h"
#include "GMTK_ObservationMatrix.h"

/**
 * collection of MixNormal objects
 */
class MixNormalCollection {

public:
  // constructors and destructor
  MixNormalCollection(RangeSetCollection rangeSetCol, unsigned numComps, unsigned maxIterKMeans,
		      bool fullCovar=false,float covAddConst=0, float covAddEpsilon=0,double clampCov=1e-10);
  ~MixNormalCollection(){delete [] _ftrMI;}

  // initialization of parameters
  //void setInitParams(unsigned dimension,
  //		     unsigned numComps,
  //	     unsigned maxIterKMeans,
  //	     bool fullCoVar = false);

  static void setVarianceFloor(double vf);
  static void setDetFloor(double df);
  static void setMixtureCoeffVanishNumber(float mcvf);
  static void setReRandOnlyOneComp(bool b);
  static void setNoReRandOnDrop(bool b);
  static void setReRandsPerMixCompRedux(unsigned r);
  ///////////////////////////////////////////////// initialization

  bool noSamplesFound(int& rangeSpecNum);
  void printNumSamples();

  // interface to EM learning
  void startEpoch();
  void addToEpoch(ObservationMatrix* obsMat,
		  size_t featureVecDim,
		  size_t totalNumFramesInSentence,
		  size_t numFramesToProcess,
		  unsigned firstFrameToProcess,
		  RangeSetCollection &tupleCol);
  void endEpoch();
  ////////////////////////////////////////////////// EM


  // compute active and MI
  int reComputeNumActive(double &maxDist, double &aveDist,
  			 double &minDist, const double termDist,
			 unsigned numEMIter);

  void computeMI(FILE *out_fp, unsigned nSamples, 
		 const RangeSetCollection &rangeCol, FILE* rangeFileFP,bool marginalizeFirstParentOut);


  // MI computation using data //////////////////////////////////////////
  void computeMIUsingData(ObservationMatrix& obsMat,
			  RangeSetCollection rangeSetCol,
			  Range &sentenceRange, 
			  const bool quiet,
			  FILE *outFileMI, Range &lrrng,int labpos,
			  FILE* rangeFileFP);


  void startEpochMI(const RangeSetCollection &rangeCol);

  void addToEpochMI(ObservationMatrix* obsMat, size_t featureSize, 
		    size_t numberOfSample, size_t n_samps, unsigned firstFrame, 
		    const RangeSetCollection &rangeSetCol);

  void endEpochMI(FILE *outFileMI, FILE* rangeFileFP);
  //////////////////////////////////////////////////// MI computation using data


  // kMeans //////////////////////////////////////////////////////
  void kmeans(ObservationMatrix* obsMat,
	      RangeSetCollection rangeSetCol,
	      Range &lrrng,
	      Range &kMeansRange,
	      unsigned numMixtures,
	      unsigned maxIter,
	      int labpos,
	      const bool quiet);

  void startKMeansEpoch(bool randLabel, bool estCov);
  void addToKMeansEpoch(ObservationMatrix* obsMat,
				  size_t featureVecDim,
				  size_t totalNumFramesInSentence,
				  size_t numFramesToProcess,
				  unsigned firstFrameToProcess,
		       RangeSetCollection &tupleCol);
  unsigned endKMeansEpoch();
  //////////////////////////////////////////////////////// kMeans

  void calcB();
  void checkB();

  // I/O ///////////////////////////////////////////////////////
  void writeCurKMeansParams(FILE *const fp, bool isBin=true);
  void readCurKMeansParams(FILE *fp, bool isBin=true);
  int readCurParams(FILE *fp, const bool forceAllActive = false, bool isBin=true);
  void writeCurParams(FILE *const fp,bool isBin=true);
  void dumpParameters(FILE* ofp, 
	              unsigned tupleNum);
  void dumpCurIterParams(int iter);
  ///////////////////////////////////////////////////////// I/O

  void generateData(FILE* ofp,unsigned numSamples,unsigned mixtureNum);
  void generateDataUsingCov(FILE* ofp,unsigned numSamples, unsigned mixtureNum);

private:
  const unsigned _numComps;
  const unsigned _numMis;  
  MixNormal * const _ftrMI;
  MixNormal * const _ftrMI_endp;

  float covAddConst; // adds a small amount of noise to the diagonal entries of covariance matrices to prevent numerical issues.
  float covAddEpsilon; // we add a random number between covAddConst-covAddEpsilon and covAddConst+Epsilon

  double clampCov; // Clamps the covariance to avoid sigularities when it gets too small 

};

#endif
