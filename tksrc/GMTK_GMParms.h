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
class Dlinks;
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

  //////////////////////////////////////////////////////////////////
  // A bitmask giving a flag if each type of EMable object is
  // currently being trained. Can only have 8*sizeof(unsigned) of these
  // for now.
  enum {
    emTrainDense1DPMFs                    = (1 << 0),
    emTrainSparse1DPMFs                   = (1 << 1),
    emTrainMeans                          = (1 << 2),
    emTrainCovars                         = (1 << 3),
    emTraindLinkMats                      = (1 << 4),
    emTrainWeightMats                     = (1 << 5),
    emTrainDiagGaussians                  = (1 << 6),
    emTrainMDCPTs                         = (1 << 7),
    emTrainMSCPTs                         = (1 << 8),
    emTrainMixGaussians                   = (1 << 9),
    emTrainGausSwitchingMixGaussians      = (1 << 10),
    emTrainLogitSwitchingMixGaussians     = (1 << 11),
    emTrainMLPSwitchingMixGaussians       = (1 << 12),
    emDefaultState                        = ~0x0
  };

  unsigned int emTrainBitmask;
  
  bool amTrainingDense1DPMFs() { return (emTrainBitmask & emTrainDense1DPMFs); }
  bool amTrainingSparse1DPMFs() { return (emTrainBitmask & emTrainSparse1DPMFs); }
  bool amTrainingMeans() { return (emTrainBitmask & emTrainMeans); }
  bool amTrainingCovars() { return (emTrainBitmask & emTrainCovars); }
  bool amTrainingdLinkMats() { return (emTrainBitmask & emTraindLinkMats); }
  bool amTrainingWeightMats() { return (emTrainBitmask & emTrainWeightMats); }
  bool amTrainingDiagGaussians() { return (emTrainBitmask & emTrainDiagGaussians); }
  bool amTrainingMDCPTs() { return (emTrainBitmask & emTrainMDCPTs); }
  bool amTrainingMSCPTs() { return (emTrainBitmask & emTrainMSCPTs); }
  bool amTrainingMixGaussians() { return (emTrainBitmask & emTrainMixGaussians); }
  bool amTrainingGausSwitchingMixGaussians() 
  { return (emTrainBitmask & emTrainGausSwitchingMixGaussians); }
  bool amTrainingLogitSwitchingMixGaussians() 
  { return (emTrainBitmask & emTrainLogitSwitchingMixGaussians); }
  bool amTrainingMLPSwitchingMixGaussians() 
  { return (emTrainBitmask & emTrainMLPSwitchingMixGaussians); }


  class ParamFile {
  public:
    string fileName;
    // number of items contained in this file.
    unsigned count;
  };
  typedef vector < ParamFile >  ParamFileVector;

  /********************************************************************/
  /********************************************************************/
  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // BASIC SHARED LOW-LEVEL SUPPORT PARAMETERS: These are the objects that 
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


  //////////////////////////////////////////////////////////////////
  // Basic Gaussian Components of various types. 
  //////////////////////////////////////////////////////////////////

  vector< GaussianComponent* > gaussianComponents;
  typedef map< string, unsigned > GaussianComponentsMapType;
  GaussianComponentsMapType gaussianComponentsMap;

  /********************************************************************/
  /********************************************************************/
  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // OBJECTS USED DIRECTLY BY DISCRETE RANDOM VARIABLES
  //////////////////////////////////////////////////////////////////

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
  /********************************************************************/
  /********************************************************************/


  //////////////////////////////////////////////////////////////////
  // OBJECTS USED DIRECTLY BY CONTINUOUS RANDOM VARIABLES
  //////////////////////////////////////////////////////////////////

  ////////////////////////////////
  // Mixtures of Gaussians (could be a heterogeneous mixutre of
  // different types above)
  vector < MixGaussians* > mixGaussians;
  typedef map< string, unsigned > MixGaussiansMapType;
  MixGaussiansMapType mixGaussiansMap;

  vector < GausSwitchingMixGaussians* > gausSwitchingMixGaussians;
  typedef map< string, unsigned > GausSwitchingMixGaussiansMapType;
  GausSwitchingMixGaussiansMapType gausSwitchingMixGaussiansMap;

  vector < LogitSwitchingMixGaussians* > logitSwitchingMixGaussians;
  typedef map< string, unsigned > LogitSwitchingMixGaussiansMapType;
  LogitSwitchingMixGaussiansMapType logitSwitchingMixGaussiansMap;

  vector < MLPSwitchingMixGaussians* > mlpSwitchingMixGaussians;
  typedef map< string, unsigned > MlpSwitchingMixGaussiansMapType;
  MlpSwitchingMixGaussiansMapType  mlpSwitchingMixGaussiansMap;

  /********************************************************************/
  /********************************************************************/
  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // STRUCTURAL OBJECTS (i.e., objects having to do with underlying
  // structure). These are not changed during parameter training
  // and are presumably learned external to this program.
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // First, A global collection of decision trees mapping vectors
  // of integers to integers. These might be used for
  // different purposes, and other strucures might index
  // into this array (or keep their own pointers to individual dts
  // that live here).
  //////////////////////////////////////////////////////////////////

  vector< RngDecisionTree<unsigned>* > dts;
  typedef map< string,  unsigned > DtsMapType;
  DtsMapType dtsMap;

  //////////////////////////////////////////////////////////////////
  // Next, the structure between observations as a collection of
  // DLINKs
  //////////////////////////////////////////////////////////////////

  vector< Dlinks* > dLinks;
  typedef map< string, unsigned > DLinksMapType;
  DLinksMapType dLinksMap;

  /********************************************************************/
  /********************************************************************/
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
  ~GMParms(); 

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
  void readDLinks(iDataStreamFile& is,bool reset = false);
  void readWeightMats(iDataStreamFile& is,bool reset = false);
  void readMdCpts(iDataStreamFile& is,bool reset = false);
  void readMsCpts(iDataStreamFile& is,bool reset = false);
  void readMtCpts(iDataStreamFile& is,bool reset = false);
  void readDTs(iDataStreamFile& is,bool reset = false);
  void readGaussianComponents(iDataStreamFile& is,bool reset = false);
  void readMixGaussians(iDataStreamFile& is,bool reset = false);
  void readGausSwitchMixGaussians(iDataStreamFile& is,bool reset = false);
  void readLogitSwitchMixGaussians(iDataStreamFile& is,bool reset = false);
  void readMlpSwitchMixGaussians(iDataStreamFile& is,bool reset = false);

  ////
  void writeDPmfs(oDataStreamFile& os);
  void writeSPmfs(oDataStreamFile& os);
  void writeMeans(oDataStreamFile& os);
  void writeCovars(oDataStreamFile& os);
  void writeDLinkMats(oDataStreamFile& os);
  void writeDLinks(oDataStreamFile& os);
  void writeWeightMats(oDataStreamFile& os);
  void writeMdCpts(oDataStreamFile& os);
  void writeMsCpts(oDataStreamFile& os);
  void writeMtCpts(oDataStreamFile& os);
  void writeDTs(oDataStreamFile& os);
  void writeGaussianComponents(oDataStreamFile& os);
  void writeMixGaussians(oDataStreamFile& os);
  void writeGausSwitchMixGaussians(oDataStreamFile& os);
  void writeLogitSwitchMixGaussians(oDataStreamFile& os);
  void writeMlpSwitchMixGaussians(oDataStreamFile& os);


  ///////////////////////////////////////////////////////////    
  // read/write in all the basic parameters, assuming file pointer 
  // is located at the correct position.
  void readBasic(iDataStreamFile& is);
  void writeBasic(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  // read/write in all the all parameters for simple
  // file format
  void readAll(iDataStreamFile& is);
  void writeAll(oDataStreamFile& os);


  ///////////////////////////////////////////////////////////    
  // read/write an entire GM (params + structure, i.e.,
  // all of the above) from a single file consisting
  // of sets of <keyword,fileName> pairs
  void read(iDataStreamFile& is,bool dataFilesAreBinary=false);
  void write(oDataStreamFile& os);

  ////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////    
  // Return the total number of parameters used by this
  // program. Use this for reporting number of parameters
  // for papers, etc.
  unsigned totalNumberParameters();


  ///////////////////////////////////////////////////////////    
  // Support for EM, applies to all EMable objects contained herein.
  ///////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  // calls the end EM on all objects EM epochs
  void emEndIteration();

  ////////////////////////////////////////////////////////////////////////////
  // calls the swap routine on all objects current and next parameters.
  void emSwapCurAndNew();


  ////////////////////////////////////////////////////////////////////////////
  void makeRandom();

  ////////////////////////////////////////////////////////////////////////////
  void makeUniform();


};

////////////////////////////////////////////////
// The global GM parameter object, must be
// actually defined near where main() is defined.
extern GMParms GM_Parms;

#endif
