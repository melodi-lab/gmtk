/*-
 * GMTK_World.h
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


#ifndef GMTK_WORLD_H
#define GMTK_WORLD_H


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"


class GM_World {

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
  // Collection of dense real matrices
  sArray< RealMatrix* > realMats;

  ////////////////////////////////
  // Collection of packed sparse real matrices
  sArray< PackedSparseRealMatrix* > psRealMats;

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
  sArray< SMDCPT* > sMdCpts;


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
  sArray < GMTK_MixGaussians* > mixGaussians;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with Gaussians.
  sArray < GMTK_GausSwitchingMixGaussians* > gausSwitchMixGaussians;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with logistic regression (i.e., 1 layer MLP)
  sArray< GMTK_LogitSwitchingMixGaussians* > logitSwitchMixGaussians;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with 2 layer (2 weight matrix) MLP
  sArray< GMTK_LogitSwitchingMixGaussians* > mlpSwitchMixGaussians;


  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // A global collection of decision trees mapping vectors
  // of integers to integers. These might be used for
  // different purposes, and other strucures might index
  // into this array for a variety of purposes.
  //////////////////////////////////////////////////////////////////

  sArray< RngDecisionTree* > dts;

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

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  World(); 

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
extern World world;

#endif // defined GMTK_CPT
