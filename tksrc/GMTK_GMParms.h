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
class MTCPT;

class GaussianComponent;
class DiagGaussian;
class LinMeanCondDiagGaussian;
class NLinMeanCondDiagGaussian;

class MixGaussiansCommon;
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
  typedef map< string, unsigned > DPmfsMapType;
  DPmfsMapType dPmfsMap;

  /////////////////////////////  
  // Collection of sparse prob. mass functions
  vector< Sparse1DPMF* > sPmfs;
  typedef map< string, unsigned  > SPmfsMapType;
  SPmfsMapType sPmfsMap;

  /////////////////////////////
  // Collection of means
  vector< MeanVector* > means;
  typedef map< string, unsigned > MeansMapType;
  MeansMapType meansMap;

  ////////////////////////////////
  // Collection of diag. covariances
  vector< DiagCovarVector* > covars;
  typedef map< string, unsigned > CovarsMapType;
  CovarsMapType covarsMap;

  ////////////////////////////////
  // Collection of objects
  // used for linear dependencies via
  // a dlink topology structure.
  vector< DlinkMatrix* > dLinkMats;
  typedef map< string, unsigned > DLinkMatsMapType;
  DLinkMatsMapType dLinkMatsMap;

  ////////////////////////////////
  // Collection of 2D Dense matrices, used
  // for weight matrices of MLPs, or 
  // for logistic regression.
  vector< WeightMatrix* > weightMats;
  typedef map< string, unsigned > WeightMatsMapType;
  WeightMatsMapType weightMatsMap;

  ////////////////////////////////
  // Collection of multi-dimensional dense CPTs
  vector< MDCPT* > mdCpts;
  typedef map< string, unsigned > MdCptsMapType;
  MdCptsMapType mdCptsMap;

  ///////////////////////////////////
  // Collection of multi-dimensional sparse CPTs (transition matrices, etc.)
  vector< MSCPT* > msCpts;
  typedef map< string, unsigned > MsCptsMapType;
  MsCptsMapType msCptsMap;

  ///////////////////////////////////
  // Collection of deterministic "CPTs" 
  vector< MTCPT* > mtCpts;
  typedef map< string, unsigned > mtCptsMapType;
  mtCptsMapType mtCptsMap;

  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // Basic Gaussian Components
  //////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////
  // Collection of diag. covariance Gaussians

  vector< GaussianComponent* > gaussianComponents;
  typedef map< string, unsigned > GaussianComponentsMapType;
  GaussianComponentsMapType gaussianComponentsMap;

  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // OBSERVATION DENSITIES
  //////////////////////////////////////////////////////////////////

  ////////////////////////////////
  // Mixtures of Gaussians (could be a heterogeneous mixutre of
  // different types above)
  vector < MixGaussiansCommon* > mixGaussians;
  typedef map< string, unsigned > MixGaussiansMapType;
  MixGaussiansMapType mixGaussiansMap;

  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // A global collection of decision trees mapping vectors
  // of integers to integers. These might be used for
  // different purposes, and other strucures might index
  // into this array (or keep their own pointers to individual dts
  // that live here).
  //////////////////////////////////////////////////////////////////

  vector< RngDecisionTree<unsigned>* > dts;
  typedef map< string,  unsigned > DtsMapType;
  DtsMapType dtsMap;

  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // Structure between observations as a collection of
  // DLINKs
  //////////////////////////////////////////////////////////////////

  vector< Dlinks* > dlinks;
  typedef map< string, unsigned > DlinksMapType;
  DlinksMapType dlinksMap;

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


  /////////////////////////////////////////////////
  // Support for "accumulative" reading into 
  // the various data structures. If 'reset' is true
  // then we clear out the arrays and maps before we
  // read them in.
  void readDPmfs(iDataStreamFile& is,bool reset = false);
  void readSPmfs(iDataStreamFile& is,bool reset = false);
  void readMeans(iDataStreamFile& is,bool reset = false);
  void readCovars(iDataStreamFile& is,bool reset = false);
  void readDLinkMats(iDataStreamFile& is,bool reset = false);
  void readWeightMats(iDataStreamFile& is,bool reset = false);
  void readDdCpts(iDataStreamFile& is,bool reset = false);
  void readMsCpts(iDataStreamFile& is,bool reset = false);
  void readMtCpts(iDataStreamFile& is,bool reset = false);

  void readGaussianComponents(iDataStreamFile& is,bool reset = false);
  void readDiagGaussians(iDataStreamFile& is,bool reset = false);
  void readLinMeanCondGaussians(iDataStreamFile& is,bool reset = false);
  void readNLinMeanCondGaussians(iDataStreamFile& is,bool reset = false);

  void readGaussianMixtures(iDataStreamFile& is,bool reset = false);
  void readMixGaussians(iDataStreamFile& is,bool reset = false);
  void readGausSwitchMixGaussians(iDataStreamFile& is,bool reset = false);
  void readLogitSwitchMixGaussians(iDataStreamFile& is,bool reset = false);
  void readMlpSwitchMixGaussians(iDataStreamFile& is,bool reset = false);

  void readDts(iDataStreamFile& is,bool reset = false);
  void readDlinks(iDataStreamFile& is,bool reset = false);

  void writeDPmfs(oDataStreamFile& os);
  void writeSPmfs(oDataStreamFile& os);
  void writeMeans(oDataStreamFile& os);
  void writeCovars(oDataStreamFile& os);
  void writeDLinkMats(oDataStreamFile& os);
  void writeWeightMats(oDataStreamFile& os);
  void writeDdCpts(oDataStreamFile& os);
  void writeMsCpts(oDataStreamFile& os);
  void writeMtCpts(oDataStreamFile& os);
  void writeDiagGaussians(oDataStreamFile& os);
  void writeLinMeanCondGaussians(oDataStreamFile& os);
  void writeNLinMeanCondGaussians(oDataStreamFile& os);
  void writeMixGaussians(oDataStreamFile& os);
  void writeGausSwitchMixGaussians(oDataStreamFile& os);
  void writeLogitSwitchMixGaussians(oDataStreamFile& os);
  void writeMlpSwitchMixGaussians(oDataStreamFile& os);
  void writeDts(oDataStreamFile& os);
  void writeDlinks(oDataStreamFile& os);


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
  // read and write the DLINKS
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
