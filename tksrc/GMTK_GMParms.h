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


#include <map>
#include <string>
#include <vector>

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
  vector< Dense1DPMF* > dPmfs;
  map< string, Dense1DPMF* > dPmfsMap;

  /////////////////////////////  
  // Collection of sparse prob. mass functions
  vector< Sparse1DPMF* > sPmfs;
  map< string, Sparse1DPMF* > sPmfsMap;

  /////////////////////////////
  // Collection of means
  vector< MeanVector* > means;
  map< string, MeanVector* > meansMap;

  ////////////////////////////////
  // Collection of diag. covariances
  vector< DiagCovarVector* > covars;
  map< string, DiagCovarVector* > covarsMap;

  ////////////////////////////////
  // Collection of objects
  // used for linear dependencies via
  // a dlink topology structure.
  vector< DlinkMatrix* > dLinkMats;
  map< string, DlinkMatrix* > dLinkMatsMap;

  ////////////////////////////////
  // Collection of 2D Dense matrices, used
  // for weight matrices of MLPs, or 
  // for logistic regression.
  vector< WeightMatrix* > weightMats;
  map< string, WeightMatrix* > weightMatsMap;

  ////////////////////////////////
  // Collection of multi-dimensional dense CPTs
  vector< MDCPT* > mdCpts;
  map< string, MDCPT* > mdCptsMap;

  ///////////////////////////////////
  // Collection of multi-dimensional sparse CPTs (transition matrices, etc.)
  vector< MSCPT* > msCpts;
  map< string, MSCPT* > msCptsMap;

  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // Basic Gaussian Components
  //////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////
  // Collection of diag. covariance Gaussians
  vector< DiagGaussian* > diagGaussians;
  map< string, DiagGaussian* > diagGaussiansMap;

  ////////////////////////////////
  // Collection of diagonal covariance Gaussians with linear mean 
  // dependency links (these cover the case full Covariance 
  // Gaussians, plus other forms such as banded or block diagonal, 
  // factored sparse inverse covariances, and so on.
  vector< LinMeanCondDiagGaussian* > linMeanCondGaussians;
  map< string, LinMeanCondDiagGaussian* > linMeanCondGaussiansMap;

  ////////////////////////////////
  // Collection of diagonal covariance Gaussians with linear and/or 
  // non-linear mean dependency links (these cover the case 
  // of "non-linear" Gaussians, and so on).
  vector< NLinMeanCondDiagGaussian* > nLinMeanCondGaussians;
  map< string, NLinMeanCondDiagGaussian* > nLinMeanCondGaussiansMap;

  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // OBSERVATION DENSITIES
  //////////////////////////////////////////////////////////////////

  ////////////////////////////////
  // Mixtures of Gaussians (could be a heterogeneous mixutre of
  // different types above)
  vector < MixGaussians* > mixGaussians;
  map< string, MixGaussians* > mixGaussiansMap;


  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with Gaussians.
  vector < GausSwitchingMixGaussians* > gausSwitchMixGaussians;
  map< string, GausSwitchingMixGaussians* > gausSwitchMixGaussiansMap;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with logistic regression (i.e., 1 layer MLP)
  vector< LogitSwitchingMixGaussians* > logitSwitchMixGaussians;
  map< string, LogitSwitchingMixGaussians* > logitSwitchMixGaussiansMap;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with 2 layer (2 weight matrix) MLP
  vector< MLPSwitchingMixGaussians* > mlpSwitchMixGaussians;
  map< string, MLPSwitchingMixGaussians* > mlpSwitchMixGaussiansMap;


  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // A global collection of decision trees mapping vectors
  // of integers to integers. These might be used for
  // different purposes, and other strucures might index
  // into this array (or keep their own pointers to individual dts
  // that live here).
  //////////////////////////////////////////////////////////////////

  vector< RngDecisionTree<int>* > dts;
  map< string,  RngDecisionTree<int>* > dtsMap;


  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // Structure between observations as a collection of
  // DLINKs
  //////////////////////////////////////////////////////////////////

  vector< Dlinks* > dlinks;
  map< string, Dlinks* > dlinksMap;

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
  // read/write in all the basic parameters, assuming file pointer 
  // is located at the correct position.
  void readBasic(iDataStreamFile& is);
  void writeBasic(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  // read and write all the DTs.
  void readDTs(iDataStreamFile& is);
  void writeDTs(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  // read and write basic Gaussians
  void readGaussians(iDataStreamFile& is);
  void writeGaussians(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  // read and write Observation Params
  void readObservationParams(iDataStreamFile& is);
  void writeObservationParams(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  // read and write Observation Structures (DLINKS)
  void readObservationStructures(iDataStreamFile& is);
  void writeObservationStructures(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  // read and write structure
  void readStructure(iDataStreamFile& is);
  void writeStructure(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  // read/write an entire GM (params + structure, i.e.,
  // all of the above) from a single file.
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  ////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////    
  // Return the total number of parameters used by this
  // program. Use this for reporting number of parameters
  // for papers, etc.
  unsigned totalNumberParameters();


};

////////////////////////////////////////////////
// The global GM parameter object, must be
// actually defined near where main() is defined.
extern GMParms GM_Parms;

#endif
