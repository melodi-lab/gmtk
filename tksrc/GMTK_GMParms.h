/*-
 * GMTK_GMParms.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *    The GM world, all aspects associated with a GM
 *    as read in from a file.
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


#ifndef GMTK_GMPARMS_H
#define GMTK_GMPARMS_H


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

/////////////////////
// forward refs
class Dense1DPMF;
class Sparse1DPMF;
class MeanVector;
class DiagCovarVector;
class RealMatrix;
class PackedSparseRealMatrix;
class DlinkMatrix;
class WeightMatrix;
class MDCPT;
class MSCPT;

class DiagGaussian;
class LinMeanCondDiagGaussian;
class NLinMeanCondDiagGaussian;

class MixGaussians;
class GausSwitchingMixGaussians;
class LogitSwitchingMixGaussians;
class MLPSwitchingMixGaussians;

template <class T>
class RngDecisionTree;
class Dlinks;
class GMTK_GM;


class GMParms {
public:

  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // BASIC SHARED LOW-LEVEL PARAMETERS: These are the objects that 
  //  higher level objects (see below) might share together.
  //  All of these objects are "EMable" in the sense that
  //  they may be trained using EM (plus possibly some other gradient
  //  based training method).
  //////////////////////////////////////////////////////////////////

  /////////////////////////////  
  // Collection of dense prob. mass functions
  sArray< Dense1DPMF* > dPmfs;

  /////////////////////////////  
  // Collection of sparse prob. mass functions
  sArray< Sparse1DPMF* > sPmfs;

  /////////////////////////////
  // Collection of means
  sArray< MeanVector* > means;

  ////////////////////////////////
  // Collection of diag. covariances
  sArray< DiagCovarVector* > covars;

  ////////////////////////////////
  // Collection of objects
  // used for linear dependencies via
  // a dlink topology structure.
  sArray< DlinkMatrix* > dLinkMats;

  ////////////////////////////////
  // Collection of 2D Dense matrices, used
  // for weight matrices of MLPs, or 
  // for logistic regression.
  sArray< WeightMatrix* > weightMats;

  ////////////////////////////////
  // Collection of multi-dimensional dense CPTs
  sArray< MDCPT* > mdCpts;

  ///////////////////////////////////
  // Collection of multi-dimensional sparse CPTs (transition matrices, etc.)
  sArray< MSCPT* > msCpts;


  //////////////////////////////////////////////////////////////////
  // Basic Gaussian Components
  //////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////
  // Collection of diag. covariance Gaussians
  sArray< DiagGaussian* > diagGaussians;

  ////////////////////////////////
  // Collection of diagonal covariance Gaussians with linear mean 
  // dependency links (these cover the case full Covariance 
  // Gaussians, plus other forms such as banded or block diagonal, 
  // factored sparse inverse covariances, and so on.
  sArray< LinMeanCondDiagGaussian* > linMeanCondGaussians;

  ////////////////////////////////
  // Collection of diagonal covariance Gaussians with linear and/or 
  // non-linear mean dependency links (these cover the case 
  // of "non-linear" Gaussians, and so on).
  sArray< NLinMeanCondDiagGaussian* > nLinMeanCondGaussians;


  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // OBSERVATION DENSITIES
  //////////////////////////////////////////////////////////////////

  ////////////////////////////////
  // Mixtures of Gaussians (could be a heterogeneous mixutre of
  // different types above)
  sArray < MixGaussians* > mixGaussians;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with Gaussians.
  sArray < GausSwitchingMixGaussians* > gausSwitchMixGaussians;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with logistic regression (i.e., 1 layer MLP)
  sArray< LogitSwitchingMixGaussians* > logitSwitchMixGaussians;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with 2 layer (2 weight matrix) MLP
  sArray< MLPSwitchingMixGaussians* > mlpSwitchMixGaussians;


  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // A global collection of decision trees mapping vectors
  // of integers to integers. These might be used for
  // different purposes, and other strucures might index
  // into this array for a variety of purposes.
  //////////////////////////////////////////////////////////////////

  sArray< RngDecisionTree<int>* > dts;

  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // Structure between observations as a collection of
  // DLINKs
  //////////////////////////////////////////////////////////////////

  sArray< Dlinks* > dlinks;

  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // Structure of model between hidden variables
  // and between hidden and observed variables.
  //////////////////////////////////////////////////////////////////

  GMTK_GM* gm;

  /********************************************************************/

  ///////////////////////////////////////////////////////////  
  // General constructor
  GMParms(); 

  ///////////////////////////////////////////////////////////    
  // read in all the basic parameters, assuming file pointer 
  // is located at the correct position.
  void readBasic(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  void writeBasic(oDataStreamFile& os);

};

////////////////////////////////////////////////
// The global GM parameter object, must be
// defined in a mainprogram.
extern GMParms GM_Parms;

#endif
