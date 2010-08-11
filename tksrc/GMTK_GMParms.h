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

#include "general.h"
#include "fileParser.h"
#include "logp.h"
#include "sArray.h"
#include "debug.h"

/////////////////////
// forward refs
class Dense1DPMF;
class Sparse1DPMF;
class MeanVector;
class DiagCovarVector;
class DlinkMatrix;
class Dlinks;
class RealMatrix;
class DirichletTable;
class MDCPT;
class MSCPT;
class MTCPT;
class NGramCPT;
class FNGramCPT;
class FNGramImp;
class VECPT;
class Vocab;
class LatticeADT;
class LatticeNodeCPT;
class LatticeEdgeCPT;

class Component;
class GaussianComponent;
class DiagGaussian;
class LinMeanCondDiagGaussian;
class NLinMeanCondDiagGaussian;

class MixtureCommon;
class Mixture;
class GausSwitchingMixture;
class LogitSwitchingMixture;
class MLPSwitchingMixture;

class RngDecisionTree;
class Dlinks;
class NameCollection;
class GMTK_GM;


class GMParms : public IM  {
public:

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

  typedef map< string, unsigned > ObjectMapType;

  template <class T> void add(T* ob, vector<T*>& arr, ObjectMapType& map) {
    if (map.find(ob->name()) != map.end())
      error("ERROR: name '%s' already exists in map.",ob->name().c_str());
    arr.push_back(ob);
    map[ob->name()] = arr.size()-1;
  }

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
  ObjectMapType dPmfsMap;
  void add( Dense1DPMF* ob );

  /////////////////////////////  
  // Collection of sparse prob. mass functions
  vector< Sparse1DPMF* > sPmfs;
  ObjectMapType sPmfsMap;
  void add( Sparse1DPMF* ob );

  /////////////////////////////
  // Collection of means
  vector< MeanVector* > means;
  ObjectMapType meansMap;
  void add(MeanVector* ob);

  ////////////////////////////////
  // Collection of diag. covariances
  vector< DiagCovarVector* > covars;
  ObjectMapType covarsMap;
  void add(DiagCovarVector* ob);

  ////////////////////////////////
  // Collection of objects
  // used for linear dependencies via
  // a dlink topology structure.
  vector< DlinkMatrix* > dLinkMats;
  ObjectMapType dLinkMatsMap;
  void add(DlinkMatrix* ob);

  ////////////////////////////////
  // Collection of 2D Dense matrices, used
  // for weight matrices of MLPs, or 
  // for logistic regression, or anytime
  // an object needs a shared scalar, vector, or matrix of 
  // real values.
  vector< RealMatrix* > realMats;
  ObjectMapType realMatsMap;
  void add(RealMatrix*ob);


  ////////////////////////////////
  // Collection of N-D Dense count (Dirichlet hyperparameter) matrices, used
  // for meta counts for DenseCPTs and Dense1DPMFs. 
  vector< DirichletTable* > dirichletTabs;
  ObjectMapType dirichletTabsMap;
  void add(DirichletTable*ob);



  //////////////////////////////////////////////////////////////////
  // Basic Gaussian Components of various types. 
  //////////////////////////////////////////////////////////////////

  vector< Component* > components;
  ObjectMapType componentsMap;
  void add(Component*);

  /********************************************************************/
  /********************************************************************/
  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // OBJECTS USED DIRECTLY BY DISCRETE RANDOM VARIABLES
  //////////////////////////////////////////////////////////////////

  ////////////////////////////////
  // Collection of multi-dimensional dense CPTs
  vector< MDCPT* > mdCpts;
  ObjectMapType mdCptsMap;
  void add(MDCPT*);

  ///////////////////////////////////
  // Collection of multi-dimensional sparse CPTs (transition matrices, etc.)
  vector< MSCPT* > msCpts;
  ObjectMapType msCptsMap;
  void add(MSCPT*);

  ///////////////////////////////////
  // Collection of deterministic "CPTs" 
  vector< MTCPT* > mtCpts;
  ObjectMapType mtCptsMap;
  void add(MTCPT*);

  ///////////////////////////////////
  // Collection of Vocab objects
  vector<Vocab*> vocabs;
  ObjectMapType vocabsMap;
  void add(Vocab*);

  ///////////////////////////////////
  // Collection of n-gram "CPTs"
  vector<NGramCPT*> ngramCpts;
  ObjectMapType ngramCptsMap;
  void add(NGramCPT*);

  ///////////////////////////////////
  // Collection of FNgram "CPTs"
  vector<FNGramCPT*> fngramCpts;
  ObjectMapType fngramCptsMap;
  void add(FNGramCPT*);
  vector<FNGramImp*> fngramImps;
  ObjectMapType fngramImpsMap;
  void add(FNGramImp*);

  ///////////////////////////////////
  // Collection of Lattice CPTs
  vector<LatticeADT*> latticeAdts;
  vector<LatticeADT*> iterableLatticeAdts;
  ObjectMapType latticeAdtsMap;
  void add(LatticeADT*);
  vector<LatticeNodeCPT*> latticeNodeCpts;
  ObjectMapType latticeNodeCptsMap;
  void add(LatticeNodeCPT*);
  vector<LatticeEdgeCPT*> latticeEdgeCpts;
  ObjectMapType latticeEdgeCptsMap;
  void add(LatticeEdgeCPT*);

  ///////////////////////////////////
  // Collection of Virtual Evidence "CPTs"
  vector<VECPT*> veCpts;
  ObjectMapType veCptsMap;
  void add(VECPT*);

  /********************************************************************/
  /********************************************************************/
  /********************************************************************/


  //////////////////////////////////////////////////////////////////
  // OBJECTS USED DIRECTLY BY CONTINUOUS RANDOM VARIABLES
  //////////////////////////////////////////////////////////////////

  ////////////////////////////////
  // Mixture distributions (could be a heterogeneous mixutre of
  // different types above)
  vector < Mixture* > mixtures;
  ObjectMapType mixturesMap;
  void add(Mixture*);

  vector < GausSwitchingMixture* > gausSwitchingMixtures;
  ObjectMapType gausSwitchingMixturesMap;
  void add(GausSwitchingMixture*);

  vector < LogitSwitchingMixture* > logitSwitchingMixtures;
  ObjectMapType logitSwitchingMixturesMap;
  void add(LogitSwitchingMixture*);

  vector < MLPSwitchingMixture* > mlpSwitchingMixtures;
  ObjectMapType  mlpSwitchingMixturesMap;
  void add(MLPSwitchingMixture*);

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

  vector< RngDecisionTree* > dts;
  vector< RngDecisionTree* > iterableDts;
  ObjectMapType dtsMap;
  void add(RngDecisionTree*);
  

  //////////////////////////////////////////////////////////////////
  // Next, the structure between observations as a collection of
  // DLINKs
  //////////////////////////////////////////////////////////////////

  vector< Dlinks* > dLinks;
  ObjectMapType dLinksMap;
  void add(Dlinks*);
  
  /********************************************************************/
  /********************************************************************/
  /********************************************************************/

  //////////////////////////////////////////////////////////////////
  // OTHER OBJECTS (i.e., objects not in any of the above categories)
  //////////////////////////////////////////////////////////////////

  vector< NameCollection* > ncls;
  ObjectMapType nclsMap;
  void add(NameCollection*);

  //////////////////////////////////////////////////////////////////
  // sorting the name collections will break any mappings that use
  // them, so they must always be unsorted afterwards
  void sort_name_collections();
  void unsort_name_collections();

  void commit_nc_changes();


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
  void readRealMats(iDataStreamFile& is,bool reset = false);
  void readDirichletTabs(iDataStreamFile& is,bool reset = false);
  void readMdCpts(iDataStreamFile& is,bool reset = false);
  void readMsCpts(iDataStreamFile& is,bool reset = false);
  void readMtCpts(iDataStreamFile& is,bool reset = false);
  void readVocabs(iDataStreamFile& is, bool reset = false);
  void readNgramCpts(iDataStreamFile& is, bool reset = false);
  void readFNgramImps(iDataStreamFile& is, bool reset = false);
  void readLatticeAdts(iDataStreamFile& is, bool reset = false);
  void readVECpts(iDataStreamFile& is, bool reset = false);
  void readDTs(iDataStreamFile& is,bool reset = false);
  void readComponents(iDataStreamFile& is,bool reset = false);
  void readMixtures(iDataStreamFile& is,bool reset = false);
  void readNameCollections(iDataStreamFile& is,bool reset = false);
  void readGausSwitchMixtures(iDataStreamFile& is,bool reset = false);
  void readLogitSwitchMixtures(iDataStreamFile& is,bool reset = false);
  void readMlpSwitchMixtures(iDataStreamFile& is,bool reset = false);

  //
  void writeDPmfs(oDataStreamFile& os);
  void writeSPmfs(oDataStreamFile& os);
  void writeMeans(oDataStreamFile& os);
  void writeCovars(oDataStreamFile& os);
  void writeDLinkMats(oDataStreamFile& os);
  void writeDLinks(oDataStreamFile& os);
  void writeRealMats(oDataStreamFile& os);
  void writeDirichletTabs(oDataStreamFile& os);
  void writeMdCpts(oDataStreamFile& os);
  void writeMsCpts(oDataStreamFile& os);
  void writeMtCpts(oDataStreamFile& os);
  //void writeVocabs(oDataStreamFile& os);
  //void writeNgramCpts(oDataStreamFile& os);
  //void writeFNgramCpts(oDataStreamFile& os);
  //void writeVECpts(oDataStreamFile& os);
  void writeDTs(oDataStreamFile& os);
  void writeComponents(oDataStreamFile& os);
  void writeMixtures(oDataStreamFile& os);
  void writeNameCollections(oDataStreamFile& os);
  void writeGausSwitchMixtures(oDataStreamFile& os);
  void writeLogitSwitchMixtures(oDataStreamFile& os);
  void writeMlpSwitchMixtures(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  // read/write in all the all parameters for simple file format; if
  // remove_unnamed is true, then Mixtures not named in a
  // NameCollection will be removed
  void readAll(iDataStreamFile& is);
  void writeAll(oDataStreamFile& os, bool remove_unnamed=false);


  ///////////////////////////////////////////////////////////    
  // read/write an entire GM (params + structure, i.e.,
  // all of the above) from a single file consisting
  // of sets of <keyword,fileName> pairs
  void read(iDataStreamFile& is);
  void write(const char *const outputFileFormat, 
	     const char * const cppCommandOptions,
	     const int intTag=CSWT_EMPTY_TAG,
	     bool remove_unnamed=false);

  ///////////////////////////////////////////////////////////    
  // read/write the trainable parameters, i.e., the ones that this program
  // might modify
  void readTrainable(iDataStreamFile& is);
  void writeTrainable(oDataStreamFile& os, bool remove_unnamed=false);

  ///////////////////////////////////////////////////////////    
  // read/write the non-trainable parameters, i.e., the ones that this program
  // will not modify (e.g., DTs, DLINK structures, etc)
  void readNonTrainable(iDataStreamFile& is);
  void writeNonTrainable(oDataStreamFile& os, bool remove_unnamed=false);

  ////////////////////////////////////////////////////////////////
  // load internal global objects, after all other parameters
  // have been read in.
  void finalizeParameters();

  ////////////////////////////////////////////////////////////////
  void begin();
  void next();
  unsigned setSegment(const unsigned segmentNo);

  // some parameter objects require the stride in the observation
  // matrix. Set this here.
  void setStride(const unsigned stride);

  // check that dlinks etc. are good w.r.t. observation stream
  void checkConsistentWithGlobalObservationStream();

  ///////////////////////////////////////////////////////////    
  // Return the total number of parameters used by this
  // program. Use this for reporting number of parameters
  // for papers, etc.
  unsigned totalNumberParameters();

  ///////////////////////////////////////////////////////////    
  // Remove mixtures that are not listed in a named collection
  void markOnlyNamedMixtures();

  ///////////////////////////////////////////////////////////    
  // Remove parameters that are not used by anyone. Return
  // true if any were removed.
  void markUsedMixtureComponents(bool remove_unnamed=false);


  ///////////////////////////////////////////////////////////    
  // using the file given by name file, set
  // all the objects given in that file so that
  // they are not trained.
  void markObjectsToNotTrain(const char*const fileName,
			     const char*const cppOptions=NULL);


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
  // EM parallelism and Kernels
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
  void emInitAccumulators(bool startEMIteration = false);
  void emWriteUnencodedAccumulators(oDataStreamFile& ofile,bool writeLogVals = true);

  ////////////////////////////////////////////////////////////////////////////
  void makeRandom();

  ////////////////////////////////////////////////////////////////////////////
  void makeUniform();

  ////////////////////////////////////////////////////////////////////////////
  // Read all decision trees listed in a master file and create support files
  void compileDecisionTrees(iDataStreamFile& is);

  ////////////////////////////////////////////////////////////////////////////
  // Set the index of the first utterance  
  void setFirstUtterance( unsigned first_index ); 

  ////////////////////////////////////////////////////////////////////////////
  // Write index files for iterable decision trees 
  void writeDecisionTreeIndexFiles();

private:

  unsigned firstUtterance; 
};

////////////////////////////////////////////////
// The global GM parameter object, must be
// actually defined near where main() is defined.
extern GMParms GM_Parms;

#endif

