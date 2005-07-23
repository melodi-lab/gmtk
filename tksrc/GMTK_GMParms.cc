/*-
 * GMTK_GMParms.cc
 *        All parameters associated with a GM.
 *
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"

#include "GMTK_GMParms.h"

#include "GMTK_Dense1DPMF.h"
#include "GMTK_Sparse1DPMF.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_RealMatrix.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_Dlinks.h"
#include "GMTK_WeightMatrix.h"
#include "GMTK_DirichletTable.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_USCPT.h"
#include "GMTK_NGramCPT.h"
#include "GMTK_FNGramCPT.h"
#include "GMTK_VECPT.h"
#include "GMTK_Vocab.h"
#include "GMTK_LatticeADT.h"
#include "GMTK_LatticeNodeCPT.h"
#include "GMTK_LatticeEdgeCPT.h"
#include "GMTK_NameCollection.h"

// particular Gaussian components
#include "GMTK_GaussianComponent.h"
#include "GMTK_DiagGaussian.h"
#include "GMTK_LinMeanCondDiagGaussian.h"

#include "GMTK_MixtureCommon.h"
#include "GMTK_Mixture.h"
#include "GMTK_ZeroScoreMixture.h"
#include "GMTK_UnityScoreMixture.h"

#include "GMTK_ObservationMatrix.h"

VCID("$Header$")

/////////////////////////////////
// an integer that specifies the maximum number of objects (such
// as means, covariances, DTs, etc.) that may be specified at
// one time in a file. This can be safely increased (to the
// extend that memory on the machine exists), but is here
// many for checking obviously invalid values.
const unsigned GMPARMS_MAX_NUM = 900000000;

////////////////////////////////
// Magic String definitions
#define MAGIC_DT_FILE "GMTK_DT_FILE"
#define MAGIC_PRM_FILE "GMTK_PRM_FILE"

////////////////////////////////////////////////
// Special name for global collection name. I.e.,
// with these names the arrays will index the global 
// entries.
#define NAMED_COLLECTION_GLOBAL_NAME "global"

// possibly use this at some point for ascii files.
// static const char* sectionSeparator = "############################################################################";


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


GMParms::GMParms()
{
  emTrainBitmask = emDefaultState;
}


GMParms::~GMParms()
{
  deleteObsInVector(dPmfs);
  deleteObsInVector(sPmfs);
  deleteObsInVector(means);
  deleteObsInVector(covars);
  deleteObsInVector(dLinkMats);
  deleteObsInVector(weightMats);
  deleteObsInVector(dirichletTabs);
  deleteObsInVector(components);
  deleteObsInVector(mdCpts);
  deleteObsInVector(msCpts);
  deleteObsInVector(mtCpts);
  deleteObsInVector(vocabs);
  deleteObsInVector(ngramCpts);
  deleteObsInVector(fngramCpts);
  deleteObsInVector(latticeAdts);
  deleteObsInVector(latticeNodeCpts);
  deleteObsInVector(veCpts);
  deleteObsInVector(mixtures);
  // deleteObsInVector(gausSwitchingMixtures);
  // deleteObsInVector(logitSwitchingMixtures);
  // deleteObsInVector(mlpSwitchingMixtures);
  deleteObsInVector(dts);
  deleteObsInVector(dLinks);
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Adding  Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void GMParms::add(Dense1DPMF* ob) { add(ob,dPmfs,dPmfsMap); }
void GMParms::add(Sparse1DPMF* ob) { add(ob,sPmfs,sPmfsMap); }
void GMParms::add(MeanVector* ob){ add(ob,means,meansMap); }
void GMParms::add(DiagCovarVector* ob){ add(ob,covars,covarsMap); }
void GMParms::add(DlinkMatrix* ob) { add(ob,dLinkMats,dLinkMatsMap); }
void GMParms::add(WeightMatrix*ob) { add(ob,weightMats,weightMatsMap); }
void GMParms::add(DirichletTable*ob) { add(ob,dirichletTabs,dirichletTabsMap); }
void GMParms::add(Component*ob){ add(ob,
				     components,
				     componentsMap); }
void GMParms::add(MDCPT*ob) { add(ob,mdCpts,mdCptsMap); }
void GMParms::add(MSCPT*ob) { add(ob,msCpts,msCptsMap); }
void GMParms::add(MTCPT*ob) { add(ob,mtCpts,mtCptsMap); }
void GMParms::add(Vocab* ob) { add(ob, vocabs, vocabsMap); }
void GMParms::add(NGramCPT* ob) { add(ob,ngramCpts,ngramCptsMap); }
void GMParms::add(FNGramCPT* ob) { add(ob,fngramCpts,fngramCptsMap); }
void GMParms::add(FNGramImp* ob) { add(ob,fngramImps,fngramImpsMap); }
void GMParms::add(LatticeADT* ob) {
	add(ob,latticeAdts,latticeAdtsMap);
	if ( ob->iterable() )
		iterableLatticeAdts.push_back(ob);
}
void GMParms::add(LatticeNodeCPT* ob) { add(ob,latticeNodeCpts,latticeNodeCptsMap); }
void GMParms::add(LatticeEdgeCPT* ob) { add(ob,latticeEdgeCpts,latticeEdgeCptsMap); }
void GMParms::add(VECPT*ob) { add(ob,veCpts,veCptsMap); }
void GMParms::add(Mixture*ob) { add(ob,mixtures,mixturesMap); }
void GMParms::add(GausSwitchingMixture*ob) { assert (0); }
void GMParms::add(LogitSwitchingMixture*ob) { assert (0); }
void GMParms::add(MLPSwitchingMixture*ob) { assert (0); }
void GMParms::add(RngDecisionTree*ob) {
  add(ob,dts,dtsMap);
  if (ob->iterable())
    iterableDts.push_back(ob);
}
void GMParms::add(Dlinks* ob) { add(ob,dLinks,dLinksMap); }


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        READING Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * readDPmfs
 *      Read in the dense PMF functions from a file.
 *
 * Preconditions:
 *      No conditions.
 *
 * Postconditions:
 *      All dense PMFs from file are now contained in internal arrays.
 *
 * Side Effects:
 *      Modifies internal arrays of object.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
GMParms::readDPmfs(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num Dense PMFs");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of dense PMFs (%d) in file '%s' line %d exceeds maximum",num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    dPmfs.resize(num);
  } else {
    start = dPmfs.size();
    dPmfs.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    Dense1DPMF* ob;

    is.read(cnt,"Can't read DPMF number");
    if (cnt != i) 
      error("ERROR: dense PMF count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new Dense1DPMF;
    ob->read(is);
    if (dPmfsMap.find(ob->name()) != dPmfsMap.end())
      error("ERROR: dense PMF named '%s' already defined but is specified for a second time in file '%s' line %d",ob->name().c_str(),is.fileName(),is.lineNo());
    dPmfs[i+start] = ob;
    dPmfsMap[ob->name()] = i+start;
  }
}


/*-
 *-----------------------------------------------------------------------
 * readSPmfs
 *      Read in the sparse PMF functions from a file.
 * 
 * Preconditions:
 *      No conditions.
 *
 * Postconditions:
 *      All sparse PMFs from file are now contained in internal arrays.
 *
 * Side Effects:
 *      Modifies internal arrays of object.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::readSPmfs(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num Sparse PMFs");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of sparse PMFs (%d) in file '%s' line %d exceeds maximum",num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    sPmfs.resize(num);
  } else {
    start = sPmfs.size();
    sPmfs.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    Sparse1DPMF* ob;

    is.read(cnt,"Can't read SPMF number");
    if (cnt != i) 
      error("ERROR: dense PMF count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new Sparse1DPMF;
    ob->read(is);
    if (sPmfsMap.find(ob->name()) != sPmfsMap.end())
      error("ERROR: sparse PMF named '%s' already defined but is specified for a second time in file '%s' line %d",ob->name().c_str(),is.fileName(),is.lineNo());
    sPmfs[i+start] = ob;
    sPmfsMap[ob->name()] = i+start;
  }
}


/*-
 *-----------------------------------------------------------------------
 * readMeans
 *      Read in the mean from a file.
 * 
 * Preconditions:
 *      No conditions.
 *
 * Postconditions:
 *      All means from file are now contained in internal arrays.
 *
 * Side Effects:
 *      Modifies internal arrays of object.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::readMeans(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num means");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of means (%d) in file '%s' line %d exceeds maximum",num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    means.resize(num);
  } else {
    start = means.size();
    means.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    MeanVector* ob;

    is.read(cnt,"Can't read mean number");
    if (cnt != i) 
      error("ERROR: mean count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new MeanVector();
    ob->read(is);
    // printf("New mean, ob's name = %s\n",ob->name().c_str());
    if (meansMap.find(ob->name()) != meansMap.end())
      error("ERROR: mean named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    means[i+start] = ob;
    meansMap[ob->name()] = i+start;
  }
}


void 
GMParms::readCovars(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num covars");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of covars (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    covars.resize(num);
  } else {
    start = covars.size();
    covars.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    DiagCovarVector* ob;

    is.read(cnt,"Can't read covar num");
    if (cnt != i) 
      error("ERROR: covar count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new DiagCovarVector;
    ob->read(is);
    if (covarsMap.find(ob->name()) != covarsMap.end())
      error("ERROR: covar named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    covars[i+start] = ob;
    covarsMap[ob->name()] = i+start;
  }
}


void 
GMParms::readDLinkMats(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num dlink matrices");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of dlink matrices (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    dLinkMats.resize(num);
  } else {
    start = dLinkMats.size();
    dLinkMats.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    DlinkMatrix* ob;

    is.read(cnt,"Can't read dlink matrix number");
    if (cnt != i) 
      error("ERROR: dlink matrix count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new DlinkMatrix;
    ob->read(is);
    if (dLinkMatsMap.find(ob->name()) != dLinkMatsMap.end())
      error("ERROR: dlink matrix named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    dLinkMats[i+start] = ob;
    dLinkMatsMap[ob->name()] = i+start;
  }
}


void 
GMParms::readDLinks(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't rad num dlinks");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of dlinks (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    dLinks.resize(num);
  } else {
    start = dLinks.size();
    dLinks.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    Dlinks* ob;

    is.read(cnt,"Can't read dlinks num");
    if (cnt != i) 
      error("ERROR: dlink count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new Dlinks;
    ob->read(is);
    if (dLinksMap.find(ob->name()) != dLinksMap.end())
      error("ERROR: dlink structure named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    dLinks[i+start] = ob;
    dLinksMap[ob->name()] = i+start;
  }
}


void 
GMParms::readWeightMats(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Cant' read num weight matrices");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of weight matrices (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    weightMats.resize(num);
  } else {
    start = weightMats.size();
    weightMats.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    WeightMatrix* ob;

    is.read(cnt,"Can't read weight mat num");
    if (cnt != i) 
      error("ERROR: weight matrix count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new WeightMatrix;
    ob->read(is);
    if (weightMatsMap.find(ob->name()) != weightMatsMap.end())
      error("ERROR: weight matrix named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    weightMats[i+start] = ob;
    weightMatsMap[ob->name()] = i+start;
  }
}



void 
GMParms::readDirichletTabs(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Cant' read num dirichlet tables");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of dirichlet tables (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    dirichletTabs.resize(num);
  } else {
    start = dirichletTabs.size();
    dirichletTabs.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    DirichletTable* ob;

    is.read(cnt,"Can't read dirichlet table num");
    if (cnt != i) 
      error("ERROR: dirichlet table count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new DirichletTable();
    ob->read(is);
    if (dirichletTabsMap.find(ob->name()) != dirichletTabsMap.end())
      error("ERROR: dirichlet table named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    dirichletTabs[i+start] = ob;
    dirichletTabsMap[ob->name()] = i+start;
  }
}


void 
GMParms::readMdCpts(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num DenseCPTs");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of Dense CPTs (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    mdCpts.resize(num);
  } else {
    start = mdCpts.size();
    mdCpts.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    MDCPT* ob;

    is.read(cnt,"Can't read DenseCPT num");
    if (cnt != i) 
      error("ERROR: Dense CPT count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new MDCPT;
    ob->read(is);
    if (mdCptsMap.find(ob->name()) != mdCptsMap.end()) {
      if (ob->name() == string(USMDCPT_NAME))
	error("ERROR: special internal unity score DenseCPT named '%s' must not be used in parameter files, as it is used internally",USMDCPT_NAME);
      else
	error("ERROR: Dense CPT named '%s' already defined but is specified for a second time in file '%s' line %d",
	      ob->name().c_str(),is.fileName(),is.lineNo());
    }
    mdCpts[i+start] = ob;
    mdCptsMap[ob->name()] = i+start;
  }
}


void 
GMParms::readMsCpts(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num SparseCPTs");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of Sparse CPTs (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    msCpts.resize(num);
  } else {
    start = msCpts.size();
    msCpts.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    MSCPT* ob;

    is.read(cnt,"Can't read SparseCPT num");
    if (cnt != i) 
      error("ERROR: Sparse CPT count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new MSCPT;
    ob->read(is);
    if (msCptsMap.find(ob->name()) != msCptsMap.end())
      error("ERROR: Sparse CPT named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    msCpts[i+start] = ob;
    msCptsMap[ob->name()] = i+start;
  }
}



void 
GMParms::readMtCpts(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num DeterministicCPTs");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of deterministic CPTs (%d) exceeds maximum",num);
  if (reset) {
    start = 0;
    mtCpts.resize(num);
  } else {
    start = mtCpts.size();
    mtCpts.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    MTCPT* ob;

    is.read(cnt,"Can't read DeterministicCPT num");
    if (cnt != i) 
      error("ERROR: deterministic CPT count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new MTCPT;
    ob->read(is);
    if (mtCptsMap.find(ob->name()) != mtCptsMap.end())
      error("ERROR: deterministic CPT named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    mtCpts[i+start] = ob;
    mtCptsMap[ob->name()] = i+start;
  }
}


void GMParms::readVocabs(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num Vocabs");
  if ( num > GMPARMS_MAX_NUM )
    error("ERROR: number of Vocabs (%d) exceeds maximum", num);
  if ( reset ) {
    start = 0;
    vocabs.resize(num);
  } else {
    start = vocabs.size();
    vocabs.resize(start + num);
  }
  for ( unsigned i = 0; i <num;i++ ) {
    // first read the count
    Vocab* ob;

    is.read(cnt, "Can't read Vocab num");
    if ( cnt != i ) 
      error("ERROR: Vocab count (%d), out of order in file '%s' line %d, expecting %d", 
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new Vocab();
    ob->read(is);
    if ( vocabsMap.find(ob->name()) != vocabsMap.end() )
      error("ERROR: Vocab named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    vocabs[i + start] = ob;
    vocabsMap[ob->name()] = i + start;
  }
}


void GMParms::readNgramCpts(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Cant' read num NgramCPTs");
  if ( num > GMPARMS_MAX_NUM )
    error("ERROR: number of NGram CPTs (%d) exceeds maximum", num);
  if ( reset ) {
    start = 0;
    ngramCpts.resize(num);
  } else {
    start = ngramCpts.size();
    ngramCpts.resize(start + num);
  }
  for ( unsigned i = 0; i <num;i++ ) {
    // first read the count
    NGramCPT* ob;

    is.read(cnt, "Can't read NgramCPT num");
    if ( cnt != i )
      error("ERROR: NGramCPT count (%d), out of order in file '%s' line %d, expecting %d", 
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new NGramCPT();
    ob->read(is);
    if ( ngramCptsMap.find(ob->name()) != ngramCptsMap.end() )
      error("ERROR: NGramCPT named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    ngramCpts[i + start] = ob;
    ngramCptsMap[ob->name()] = i + start;
  }
}


void GMParms::readFNgramImps(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num, "Cant' read num FNgramCPTs");
  if ( num > GMPARMS_MAX_NUM )
    error("ERROR: number of FNGram CPTs (%d) exceeds maximum", num);
  if ( reset ) {
    start = 0;
    fngramImps.resize(num);
  } else {
    start = fngramImps.size();
    fngramImps.resize(start + num);
  }
  for ( unsigned i = 0; i <num; i++ ) {
    // first read the count
    FNGramImp* ob;

    is.read(cnt, "Can't read FNGramCPT num");
    if ( cnt != i )
      error("ERROR: FNGramCPT count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt, is.fileName(),is.lineNo(), i);

    ob = new FNGramImp();
    ob->read(is);
    if ( fngramImpsMap.find(ob->name()) != fngramImpsMap.end() )
      error("ERROR: FNGramCPT named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(), is.fileName(),is.lineNo());
    fngramImps[i + start] = ob;
    fngramImpsMap[ob->name()] = i + start;
  }
}


void GMParms::readLatticeAdts(iDataStreamFile& is, bool reset) {
	unsigned num;
	unsigned cnt;
	unsigned start = 0;

	is.read(num, "Can't read num Lattice CPTs");
	if ( num > GMPARMS_MAX_NUM )
		error("ERROR: number of Lattice CPTs (%d) exceeds maximum", num);

	if ( reset ) {
		start = 0;
		latticeAdts.resize(num);
		latticeNodeCpts.resize(num);
		latticeEdgeCpts.resize(num);
	} else {
		start = latticeAdts.size();
		latticeAdts.resize(start + num);
		latticeNodeCpts.resize(start + num);
		latticeEdgeCpts.resize(start + num);
	}

	for ( unsigned i = 0; i < num; i++ ) {
		is.read(cnt, "Can't read Lattice CPT index");
		if ( cnt != i )
		  error("ERROR: Lattice CPT count (%d), out of order in file '%s' line %d, expecting %d", cnt, is.fileName(),is.lineNo(), i);
		LatticeADT* ob = new LatticeADT();
		ob->read(is);
		if ( latticeAdtsMap.find(ob->name()) != latticeAdtsMap.end() )
			error("ERROR: Lattice CPT named '%s' already defined but is specified for a second time in file '%s' line %d", ob->name().c_str(), is.fileName(),is.lineNo());
		latticeAdts[i+start] = ob;
		latticeAdtsMap[ob->name()] = i + start;

		// add if it iterable
		if ( ob->iterable() )
			iterableLatticeAdts.push_back(ob);

		// add node CPT
		LatticeNodeCPT* ndCpt = new LatticeNodeCPT();
		ndCpt->setLatticeADT(*ob);
		ndCpt->setName(ob->name());
		latticeNodeCpts[i+start] = ndCpt;
		latticeNodeCptsMap[ob->name()] = i + start;

		// add edge CPT
		LatticeEdgeCPT* edgeCpt = new LatticeEdgeCPT();
		edgeCpt->setLatticeADT(*ob);
		edgeCpt->setName(ob->name());
		latticeEdgeCpts[i+start] = edgeCpt;
		latticeEdgeCptsMap[ob->name()] = i + start;
	}
}


void GMParms::readVECpts(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num, "Can't read num VirtualEvidenceCPTs");
  if ( num > GMPARMS_MAX_NUM )
    error("ERROR: number of Ve CPTs (%d) exceeds maximum", num);
  if ( reset ) {
    start = 0;
    veCpts.resize(num);
  } else {
    start = veCpts.size();
    veCpts.resize(start + num);
  }
  for ( unsigned i = 0; i <num; i++ ) {
    // first read the count
    VECPT* ob;

    is.read(cnt, "Can't read VirtualEvidenceCPT num");
    if ( cnt != i )
      error("ERROR: VECPT count (%d), out of order in file '%s' line %d, expecting %d", 
	    cnt, is.fileName(), is.lineNo(),i);

    ob = new VECPT();
    ob->read(is);
    if ( veCptsMap.find(ob->name()) != veCptsMap.end() )
      error("ERROR: VECPT named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(), is.fileName(),is.lineNo());
    veCpts[i + start] = ob;
    veCptsMap[ob->name()] = i + start;
  }
}




void
GMParms::readDTs(
  iDataStreamFile& is, 
  bool reset 
  )
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num Decision Trees");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of DTs (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    dts.resize(num);
  } else {
    start = dts.size();
    dts.resize(start+num);
  }

  for (unsigned i=0;i<num;i++) {
    // first read the count
    RngDecisionTree* ob;

    is.read(cnt,"Can't read DecisionTree num");
    if (cnt != i) 
      error("ERROR: DT count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    ob = new RngDecisionTree;
    ob->read(is);
    if (dtsMap.find(ob->name()) != dtsMap.end())
      error("ERROR: DT named '%s' already defined but is specified for a second time in file '%s' line %d",
	    ob->name().c_str(),is.fileName(),is.lineNo());
    dts[i+start] = ob;
    dtsMap[ob->name()] = i+start;
    if (ob->iterable()) {
      iterableDts.push_back(ob);
    } 
  }   
}


void 
GMParms::readComponents(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;

  unsigned start = 0;
  is.read(num,"Can't read num components");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of components (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    start = 0;
    components.resize(num);
  } else {
    start = components.size();
    components.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    Component* gc = NULL;

    is.read(cnt,"Can't read component num");
    if (cnt != i) 
      error("ERROR: component count (%d), out of order in file '%s' line %d , expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    // next read the dimension of this component
    int dim;
    is.read(dim,"Can't read component dimension");

    // read the Gaussian type, note that when
    // this is written, the object itself will write the type.
    int t;
    is.read(t,"Can't read component type");
    if (t == Component::DiagGaussian) {
      gc = new DiagGaussian(dim);
    } else if (t == Component::LinMeanCondDiagGaussian) {
      gc = new LinMeanCondDiagGaussian(dim);
    } else if (t == Component::PolyNLinMeanCondDiagGaussian) {
      error("PolyNLinMeanCondDiag not implemented");
      // gc = new PolyNLinMeanCondDiagGaussian(dim);
    } else {
      error("Error: reading file %s line %d, unknown component type %d in file",
	    is.fileName(),is.lineNo(),t);
    }
    gc->read(is);

    // this next check is redundant since the dim > 0 check is
    // already done by the mean, variance, etc. objects. We leave
    // it here, however, since 1) it costs almost nothing, and 2) as new object types 
    // are added, we might need such a check here.
    if (dim <= 0)
      error("ERROR: component named '%s' in file '%s' line %d specifies a non-positive dimension (%d). Must be > 0.",
	    gc->name().c_str(),is.fileName(),is.lineNo(),dim);

    if (componentsMap.find(gc->name()) != componentsMap.end())
      error("ERROR: component named '%s' already defined but is specified for a second time in file '%s' line %d",
	    gc->name().c_str(),is.fileName(),is.lineNo());
    components[i+start] = gc;
    componentsMap[gc->name()] = i+start;
  }
}


void 
GMParms::readMixtures(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num mixtures");
  if (num > GMPARMS_MAX_NUM) error("ERROR: number of mixtures (%d) in file '%s' line %d exceeds maximum",
				   num,is.fileName(),is.lineNo());
  if (reset) {
    // this isn't implemented at the moment.
    assert(0);
    // start = 0;
    // mixtures.resize(num);
  } else {
    start = mixtures.size();
    mixtures.resize(start+num);
  }


  for (unsigned i=0;i<num;i++) {
    // first read the count
    Mixture* gm;

    is.read(cnt,"Can't read mixture num");
    if (cnt != i) 
      error("ERROR: reading file %s line %d, mixture count (%d), out of order, expecting %d",
	    is.fileName(),is.lineNo(),cnt,i);


    // next read the dimension of this mixture
    int dim;
    is.read(dim,"Can't read mixture dimension");

    gm = new Mixture(dim);
    gm->read(is);

    // this next check is redundant since the dim > 0 check is
    // already done by the mean, variance, and component objects. We leave
    // it here, however, since 1) it costs almost nothing, and 2) as new object types 
    // are added, we might need such a check here.
    if (dim <= 0)
      error("ERROR: mixture named '%s' in file '%s' line %d specifies a non-positive dimension (%d). Must be > 0.",
	    gm->name().c_str(),is.fileName(),is.lineNo(),dim);

    if (mixturesMap.find(gm->name()) != mixturesMap.end()) {
      if (gm->name() == string(ZEROSCOREMIXTURE_NAME))
	error("ERROR: special internal mixture named '%s' must not be used in parameter files, as it is used internally",ZEROSCOREMIXTURE_NAME);
      else if (gm->name() == string(UNITYSCOREMIXTURE_NAME))
	error("ERROR: special internal mixture named '%s' must not be used in parameter files, as it is used internally",UNITYSCOREMIXTURE_NAME);
      else
	error("ERROR: mixture named '%s' already defined but is specified for a second time in file '%s' line %d",
	      gm->name().c_str(),is.fileName(),is.lineNo());
    }


    mixtures[i+start] = gm;
    mixturesMap[gm->name()] = i+start;
  }
}


void 
GMParms::readNameCollections(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;
  unsigned start = 0;

  is.read(num,"Can't read num NamedCollections");
  if (reset) {
    start = 0;
    ncls.resize(num);
  } else {
    start = ncls.size();
    ncls.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    NameCollection* nc;

    is.read(cnt,"Can't read NameCollection num");
    if (cnt != i) 
      error("ERROR: collection order count (%d), out of order in file '%s' line %d, expecting %d",
	    cnt,is.fileName(),is.lineNo(),i);

    nc = new NameCollection();
    nc->read(is);

    if (nc->name() == NAMED_COLLECTION_GLOBAL_NAME) {
      error("ERROR: special internal collection name '%s' can not be defined in file '%s' line %d",
	    nc->name().c_str(),is.fileName(),is.lineNo());

    }

    if (nclsMap.find(nc->name()) != nclsMap.end()) {
      if (nc->name() == string(NAMED_COLLECTION_GLOBAL_NAME))
	error("ERROR: special internal collection named '%s' must not be used in parameter files, as it is used internally to refer to global table.",NAMED_COLLECTION_GLOBAL_NAME);
      else
	error("ERROR: collection named '%s' already defined but is specified for a second time in file '%s' line %d",
	      nc->name().c_str(),is.fileName(),is.lineNo());
    }

    ncls[i+start] = nc;
    nclsMap[nc->name()] = i+start;
  }
}



void 
GMParms::readGausSwitchMixtures(iDataStreamFile& is, bool reset)
{
  unsigned num;
  is.read(num,"Can't read num GSMGs");
  if (num > 0)
    error("ERROR: reading file '%s' line %d, GausSwitchMixtures not implemented just yet",
	  is.fileName(),is.lineNo());
}

void 
GMParms::readLogitSwitchMixtures(iDataStreamFile& is, bool reset)
{
  unsigned num;
  is.read(num,"Can't read num LSMGs");
  if (num > 0)
    error("ERROR: reading file '%s' line %d, LogitSwitchMixtures not implemented just yet",
	  is.fileName(),is.lineNo());
}

void 
GMParms::readMlpSwitchMixtures(iDataStreamFile& is, bool reset)
{
  unsigned num;
  is.read(num,"Can't read num MSMGs");
  if (num > 0)
    error("ERROR: reading file '%s' line %d, MLPSwitchMixtures not implemented just yet",
	  is.fileName(),is.lineNo());
}



/*-
 *-----------------------------------------------------------------------
 * readAll
 *   read in everything (all parameters) in one go.
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      files read in
 *
 * Side Effects:
 *      changes internal arrays.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
GMParms::readAll(iDataStreamFile& is)
{
  // just read everything in one go.

  // first read structural items
  readDTs(is);
  readDLinks(is);

  // then read basic numeric items.
  readDPmfs(is);
  readSPmfs(is);
  readMeans(is);
  readCovars(is);
  readDLinkMats(is);
  readWeightMats(is);
  readDirichletTabs(is);
  readMdCpts(is);
  readMsCpts(is);
  readMtCpts(is);
  readVocabs(is);
  readNgramCpts(is);
  readFNgramImps(is);
  readLatticeAdts(is);
  readVECpts(is);

  // next read definitional items
  readComponents(is);
  readMixtures(is);
  readGausSwitchMixtures(is);
  readLogitSwitchMixtures(is);
  readMlpSwitchMixtures(is);  
}


/*-
 *-----------------------------------------------------------------------
 * readTrainable
 *   read in just the trainable parameters, i.e.,
 *   the things that might be modified by the
 *   program when training is occuring.
 *
 *   Note that the union of this routine
 *   and readNonTrainable should be the same as 
 *   readAll
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      files read in
 *
 * Side Effects:
 *      changes internal arrays.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::readTrainable(iDataStreamFile& is)
{
  infoMsg(Med+8,"Reading Dense Prob. Mass Functions\n");
  readDPmfs(is);
  infoMsg(Med+8,"Reading Sparse Prob. Mass Functions\n");
  readSPmfs(is);
  infoMsg(Med+8,"Reading Means Functions\n");
  readMeans(is);
  infoMsg(Med+8,"Reading Covars Functions\n");
  readCovars(is);
  infoMsg(Med+8,"Reading Dlink Mats\n");
  readDLinkMats(is);
  infoMsg(Med+8,"Reading Weight  Mats\n");
  readWeightMats(is);  
  infoMsg(Med+8,"Reading Dense CPTs\n");
  readMdCpts(is);

  // next read definitional items
  infoMsg(Med+8,"Reading Components\n");
  readComponents(is);
  infoMsg(Med+8,"Reading Mixtures\n");
  readMixtures(is);
  readGausSwitchMixtures(is);
  readLogitSwitchMixtures(is);
  readMlpSwitchMixtures(is);  
}



/*-
 *-----------------------------------------------------------------------
 * readnonTrainable
 *   read in just the *NON*-trainable parameters, i.e.,
 *   the things that, guaranteed, will not be
 *   modified when the program is training.
 *
 *   Note that the union of this routine
 *   and readTrainable should be the same as 
 *   readAll
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      files read in
 *
 * Side Effects:
 *      changes internal arrays.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::readNonTrainable(iDataStreamFile& is)
{
  // first read structural items
  infoMsg(Med+8,"Reading Decision Trees\n");
  readDTs(is);
  infoMsg(Med+8,"Reading Dlinks\n");
  readDLinks(is);

  infoMsg(Med+8,"Reading SparseCPTs\n");
  readMsCpts(is);
  infoMsg(Med+8,"Reading Deterministic CPTs\n");
  readMtCpts(is);
  infoMsg(Med+8,"Reading Vocabs\n");
  readVocabs(is);
  infoMsg(Med+8,"Reading NgramCPTs\n");
  readNgramCpts(is);
  infoMsg(Med+8,"Reading FNgramCPTs\n");
  readFNgramImps(is);
  infoMsg(Med+8,"Reading VirtualEvidenceCPTs\n");
  readVECpts(is);
}


void 
GMParms::read(
  iDataStreamFile& is
  )
{
  // read a file consisting of a list of keyword,filename
  // pairs. the keyword says which structure to read in,
  // and the filename says where to get it.
  string keyword;
  string fileName;
  string binStatus;

  const string INLINE_FILE_KEYWORD("inline");
  
  map<string,iDataStreamFile*> fileNameMap;
  fileNameMap[INLINE_FILE_KEYWORD] = &is;

  while (is.readString(keyword)) {

    if (!is.readString(fileName)) {
      error("ERROR: while reading file '%s' line %d , got keyword '%s' without a filename",
	    is.fileName(),is.lineNo(),keyword.c_str());
    }

    bool binary_p = is.binary();
    if (fileName != INLINE_FILE_KEYWORD) {
      // read binary status of file if this is not an inline declarator
      if (!is.readString(binStatus)) {
	error("ERROR: while reading file '%s' line %d, got keyword '%s' and filename '%s' without a binary status",
	      is.fileName(),is.lineNo(),keyword.c_str(),fileName.c_str());
      }
      if (binStatus == "ascii" || binStatus == "ASCII")
	binary_p = false;
      else if (binStatus == "binary" || binStatus == "BINARY")
	binary_p = true;
      else {
	error("ERROR: while reading file '%s' line %d, got string '%s' when expecting 'ascii'/'binary' keyword",
	      is.fileName(),is.lineNo(),binStatus.c_str());
      }
      infoMsg(Med+8,"Reading keyword '%s' from %s file '%s'.\n",
	      keyword.c_str(),(binary_p?"binary":"ASCII"),
	      fileName.c_str());
    } else {
      infoMsg(Med+8,"Reading keyword '%s' from inline.\n",keyword.c_str());
    }

    map<string,iDataStreamFile*>::iterator it = fileNameMap.find(fileName);
    if (it == fileNameMap.end()) {
      fileNameMap[fileName] = new iDataStreamFile(fileName.c_str(),binary_p);
      it = fileNameMap.find(fileName);
    } else if (((*it).second)->binary() != binary_p) {
      error("ERROR: reading '%s'. File '%s' line %d had binary status = %d, now binary status = %d. Can't mix binary and ASCII files",
	    is.fileName(),is.lineNo(),fileName.c_str(),((*it).second)->binary(),binary_p);
    }

    if (keyword == "DPMF_IN_FILE") {
      readDPmfs(*((*it).second),false);

    } else if (keyword == "SPMF_IN_FILE") {
      readSPmfs(*((*it).second),false);

    } else if (keyword == "MEAN_IN_FILE") {
      readMeans(*((*it).second),false);

    } else if (keyword == "COVAR_IN_FILE") {
      readCovars(*((*it).second),false);

    } else if (keyword == "DLINK_MAT_IN_FILE") {
      readDLinkMats(*((*it).second),false);

    } else if (keyword == "DLINK_IN_FILE") {
      readDLinks(*((*it).second),false);

    } else if (keyword == "WEIGHT_MAT_IN_FILE") {
      readWeightMats(*((*it).second),false);

    } else if (keyword == "DIRICHLET_TAB_IN_FILE") {
      readDirichletTabs(*((*it).second),false);

    } else if (keyword == "DENSE_CPT_IN_FILE") {
      readMdCpts(*((*it).second),false);

    } else if (keyword == "SPARSE_CPT_IN_FILE") {
      readMsCpts(*((*it).second),false);

    } else if (keyword == "DETERMINISTIC_CPT_IN_FILE") {
      readMtCpts(*((*it).second),false);

    } else if (keyword == "VOCAB_IN_FILE") {
      readVocabs(*((*it).second),false);

    } else if (keyword == "NGRAM_CPT_IN_FILE") {
      readNgramCpts(*((*it).second),false);

    } else if (keyword == "FNGRAM_CPT_IN_FILE") {
      readFNgramImps(*((*it).second),false);

    } else if (keyword == "LATTICE_CPT_IN_FILE") {
	    readLatticeAdts(*((*it).second),false);

    } else if (keyword == "VE_CPT_IN_FILE") {
      readVECpts(*((*it).second),false);

    } else if (keyword == "DT_IN_FILE") {
      readDTs(*((*it).second),false);

    } else if (keyword == "MC_IN_FILE") {
      readComponents(*((*it).second),false);

    } else if (keyword == "MX_IN_FILE") {
      readMixtures(*((*it).second),false);

    } else if (keyword == "NAME_COLLECTION_IN_FILE") {
      readNameCollections(*((*it).second),false);

    } else if (keyword == "GSMG_IN_FILE") {
      error("GSMG_IN_FILE in file '%s' line %d, not implemented",
	    is.fileName(),is.lineNo());

    } else if (keyword == "LSMG_IN_FILE") {
      error("LSMG_IN_FILE in file '%s' line %d, not implemented",
	    is.fileName(),is.lineNo());

    } else if (keyword == "MSMG_IN_FILE") {
      error("MSMG_IN_FILE in file '%s' line %d, not implemented",
	    is.fileName(),is.lineNo());

    } else {
      error("ERROR: encountered unknown file type '%s' in file '%s' line %d",
	    keyword.c_str(),is.fileName(),is.lineNo());
    }
  }

  // now go through and delete all the input files
  for (map<string,iDataStreamFile*>::iterator it = fileNameMap.begin();
       it != fileNameMap.end(); it++) {
    if ((*it).first != INLINE_FILE_KEYWORD) {
      // don't delete is
      delete ((*it).second);
    }
  }

}


/*-
 *-----------------------------------------------------------------------
 * GMParms::writeDecisionTreeIndexFiles()
 *
 * Preconditions:
 *    The iterableDts vector has been filled in using readDTs 
 *
 * Postconditions:
 *    An index file is written for each iterable DT
 *
 * Side Effects:
 *    The first decison tree is set to 0 for all trees. 
 *
 * Results:
 *    none
 *-----------------------------------------------------------------------
 */
void 
GMParms::
writeDecisionTreeIndexFiles()
{
  vector<RngDecisionTree*>::iterator crrnt_tree; 
  vector<RngDecisionTree*>::iterator end_tree; 

  for (crrnt_tree = iterableDts.begin(), 
       end_tree   = iterableDts.end();
       crrnt_tree != end_tree;
       crrnt_tree++) {
 
    (*crrnt_tree)->setFirstDecisionTree(0);
    (*crrnt_tree)->writeIndexFile();
  } 
}


/*-
 *-----------------------------------------------------------------------
 * GMParms::finalizeParameters()
 *   load internal global objects, after all other parameters
 *   have been read in.
 *
 * Preconditions:
 *      Should be called after all internal objects have been read in. This routine
 *      should not be called before other GMTK objects have
 *      been loaded. Note also this should be done *before* any association
 *      with the RV structure, so that it doesn't have any undefined values.
 *
 * Postconditions:
 *      New global objects have been done.
 *
 * Side Effects:
 *      changes internal GMTK object arrays. Note that
 *      this routine will add to the internal GMKT object arrays
 *      by appending at the end. This routine should be called last,
 *      after all other objects have been allocated. Also, when
 *      writing out any of these objects, we should be sure
 *      not to write out any of the objects which are being
 *      stored here.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::finalizeParameters()
{
  ///////////////////////////////////////////////////////
  // Now that presumably everything has been read in,
  // we insert the global internal objects:
  //     1) a named collection which references the global arrays.
  //     2) two special scoring mixtures (unity and zero)
  //     3) one special scoring MDCPT (unity)

  // Load the global named collection.  first, make sure that the name
  // hasn't already been defined, meaning either that this routine has
  // not been called or nobody else defined a collection with this name.
  assert (nclsMap.find(string(NAMED_COLLECTION_GLOBAL_NAME)) == nclsMap.end());

  NameCollection* nc = new NameCollection();
  nc->_name = NAMED_COLLECTION_GLOBAL_NAME;
  // copy the tables:
  // TODO: figure out a better way than copying the entire collection
  // to global. Perhaps, however, STL does the right thing and will
  // share the internal array, but should check this.
  nc->mxTable = mixtures;
  nc->spmfTable = sPmfs;
  ncls.push_back(nc);
  nclsMap[nc->name()] = ncls.size()-1;

  /////////////////////////////////////////////////////////////////////
  // now we load 2 extra mixtures. 
  // IMPORTANT: when making changes here, also see routine:  GMParms::writeMixtures(oDataStreamFile& os)

  // Load the zero scoring Mixture
  assert (mixturesMap.find(string(ZEROSCOREMIXTURE_NAME)) == mixturesMap.end());
  ZeroScoreMixture* zs = new ZeroScoreMixture();
  mixtures.push_back(zs);
  mixturesMap[zs->name()] = mixtures.size()-1;

  // Load the zero scoring Mixture
  assert (mixturesMap.find(string(UNITYSCOREMIXTURE_NAME)) == mixturesMap.end());
  UnityScoreMixture* us = new UnityScoreMixture();
  mixtures.push_back(us);
  mixturesMap[us->name()] = mixtures.size()-1;


  /////////////////////////////////////////////////////////////////////
  // and we load 1 extra MDCPT
  // IMPORTANT: when making changes here, also see routine:  GMParms::writeMdCpts(oDataStreamFile& os)

  // load the unity scoring MDCPT.  Note, this CPT corresponds to
  // a random variable:
  //   1) of any cardinality
  //   2) with no switching or conditional parents
  //   3) that is observed
  //   4) that has probability = 1 regardless of the
  //      observed value, so that this is a conditional rather than
  //      a scoring observation, similar to the way conditional
  //      floating point observations work.
  assert (mdCptsMap.find(string(USMDCPT_NAME)) == mdCptsMap.end());
  USCPT *uscpt = new USCPT();
  mdCpts.push_back(uscpt);
  mdCptsMap[uscpt->name()] = mdCpts.size()-1;

}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        WRITING Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * write<OBJECTS>
 *      Write the <OBJECT> to a file specified by os
 * 
 * Preconditions:
 *      None, as if the arrays are zero sized, nothing will be written
 *      other than the number 0
 *
 * Postconditions:
 *      All objects have been written.
 *      
 * Side Effects:
 *      None.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */


void 
GMParms::writeDPmfs(oDataStreamFile& os)
{
  os.nl(); os.writeComment("dense PMFs");os.nl();
  os.write(dPmfs.size(),"num dPMFs"); os.nl();
  for (unsigned i=0;i<dPmfs.size();i++) {
    // first write the count
    os.write(i,"dDPMF cnt");
    os.nl();
    dPmfs[i]->write(os);
  }
  os.nl();
}



void 
GMParms::writeSPmfs(oDataStreamFile& os)
{
  os.nl(); os.writeComment("sparse PMFs");os.nl();
  os.write(sPmfs.size(),"num sPMFs"); os.nl();
  for (unsigned i=0;i<sPmfs.size();i++) {
    // first write the count
    os.write(i,"sPMFs cnt");
    os.nl();
    sPmfs[i]->write(os);
  }
  os.nl();
}


/*-
 *-----------------------------------------------------------------------
 * writeMeans
 *      writes out the means 
 * 
 * Preconditions:
 *      markUsedMixtureComponents() should have been called immediately before.
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeMeans(oDataStreamFile& os)
{
  os.nl(); os.writeComment("means");os.nl();
  // first scan through and find the number that are used.
  unsigned used = 0;
  for (unsigned i=0;i<means.size();i++)
    if (means[i]->emUsedBitIsSet())
      used++;
  if (used != means.size())
    warning("NOTE: saving only %d used means out of a total of %d",
	    used,means.size());
  os.write(used,"num Means"); os.nl();
  unsigned index=0;
  for (unsigned i=0;i<means.size();i++) {
    if (means[i]->emUsedBitIsSet()) {
      // first write the count
      os.write(index++,"means cnt");
      os.nl();
      means[i]->write(os);
    }
  }


  assert ( used == index );
  os.nl();
}



/*-
 *-----------------------------------------------------------------------
 * writeMeans
 *      writes out the means 
 * 
 * Preconditions:
 *      markUsedMixtureComponents() should have been called immediately before.
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeCovars(oDataStreamFile& os)
{
  os.nl(); os.writeComment("diagonal covariance matrices");os.nl();
  unsigned used = 0;
  for (unsigned i=0;i<covars.size();i++)
    if (covars[i]->emUsedBitIsSet())
      used++;
  if (used != covars.size())
    warning("NOTE: saving only %d used covariances out of a total of %d",
	    used,covars.size());
  os.write(used,"num covars"); os.nl();
  unsigned index=0;
  for (unsigned i=0;i<covars.size();i++) {
    if (covars[i]->emUsedBitIsSet()) {
      // first write the count
      os.write(index++,"covar cnt");
      os.nl();
      covars[i]->write(os);
    }
  }
  assert ( used == index );
  os.nl();
}



/*-
 *-----------------------------------------------------------------------
 * writeDlinkmats
 *      writes out the dlink mats
 * 
 * Preconditions:
 *      markUsedMixtureComponents() should have been called immediately before.
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeDLinkMats(oDataStreamFile& os)
{
  os.nl(); os.writeComment("dlink matrices");os.nl();
  unsigned used = 0;
  for (unsigned i=0;i<dLinkMats.size();i++) {
    if (dLinkMats[i]->emUsedBitIsSet())
      used++;
  }
  if (used != dLinkMats.size())
    warning("NOTE: saving only %d used dlink matrices out of a total of %d",
	    used,dLinkMats.size());
  os.write(used,"num dlink mats"); os.nl();
  unsigned index=0;
  for (unsigned i=0;i<dLinkMats.size();i++) {
    if (dLinkMats[i]->emUsedBitIsSet()) {
      // first write the count
      os.write(index++,"dlink mat cnt");
      os.nl();
      dLinkMats[i]->write(os);
    }
  }
  assert ( used == index );
  os.nl();
}



/*-
 *-----------------------------------------------------------------------
 * writeDlinks
 *      writes out the dlinks
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeDLinks(oDataStreamFile& os)
{
  os.nl(); os.writeComment("dlink structures");os.nl();
  os.write(dLinks.size(),"num dlinks"); os.nl();
  for (unsigned i=0;i<dLinks.size();i++) {
    // first write the count
    os.write(i,"dlink cnt");
    os.nl();
    dLinks[i]->write(os);
  }
  os.nl();
}



/*-
 *-----------------------------------------------------------------------
 * writeWeightMats
 *      writes out the weight mats
 * 
 * Preconditions:
 *      markUsedMixtureComponents() should have been called immediately before.
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeWeightMats(oDataStreamFile& os)
{
  os.nl(); os.writeComment("weight matrices");os.nl();
  unsigned used = 0;
  for (unsigned i=0;i<weightMats.size();i++) {
    if (weightMats[i]->emUsedBitIsSet())
      used++;
  }
  if (used != weightMats.size())
    warning("NOTE: saving only %d used weight  matrices out of a total of %d",
	    used,weightMats.size());
  os.write(used,"num weight mats"); os.nl();
  unsigned index = 0;
  for (unsigned i=0;i<weightMats.size();i++) {
    if (weightMats[i]->emUsedBitIsSet()) {
      // first write the count
      os.write(index++,"weight mat cnt");
      os.nl();
      weightMats[i]->write(os);
    }
  }
  assert ( used == index );
  os.nl();
}




/*-
 *-----------------------------------------------------------------------
 * writeDiricletTabs
 *      writes out the Dirichlet Tables
 * 
 * Preconditions:
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeDirichletTabs(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Dirichlet Tables");os.nl();
  os.write(weightMats.size(),"num dirichlet tabs"); os.nl();
  for (unsigned i=0;i<dirichletTabs.size();i++) {
    // first write the count
    os.write(i,"diriclet tab cnt");
    os.nl();
    dirichletTabs[i]->write(os);
  }
  os.nl();
}



/*-
 *-----------------------------------------------------------------------
 * writeMdCpts
 *
 *   Writes out all the MDCPTs in the global object. If there are
 *   any internal pre-defined CPTs in this set, then they are
 *   assumed to be kept at the beginning of the array, and are 
 *   not written out (since they are automatically created anew each
 *   time the program loads). See routine loadGlobal() for more details.
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeMdCpts(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Dense CPTs");os.nl();
  // leave out the 1st one (ie., the -1) as it is an internal
  // object. See routine loadGlobal()
  os.write(mdCpts.size()-1,"num Dense CPTs"); os.nl();

  // Next, get a pointer to the unity score CPT that we should not
  // write out.  Note that it potentially might not be at the end of
  // the array since we might have done automatic allocation of MDCPTs
  // when reading in the .str file.
  const string usname = string(USMDCPT_NAME);
  assert (mdCptsMap.find(usname) != mdCptsMap.end());
  const unsigned idx = mdCptsMap[string(USMDCPT_NAME)];
  assert ( idx < mdCpts.size() );
  USCPT *uscpt = (USCPT*) mdCpts[idx];
  unsigned cnt = 0;
  for (unsigned i=0;i<mdCpts.size();i++) {
    if (mdCpts[i] == uscpt)
      continue;
    // first write the count
    os.write(cnt++,"Dense CPT cnt");
    os.nl();
    mdCpts[i]->write(os);
  }
  os.nl();
}


/*-
 *-----------------------------------------------------------------------
 * writemscpts
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeMsCpts(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Sparse CPTs");os.nl();
  os.write(msCpts.size(),"num Sparse CPTs"); os.nl();
  for (unsigned i=0;i<msCpts.size();i++) {
    // first write the count
    os.write(i,"Sparse CPT cnt");
    os.nl();
    msCpts[i]->write(os);
  }
  os.nl();
}



/*-
 *-----------------------------------------------------------------------
 * writemtcpts
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeMtCpts(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Deterministic CPTs");os.nl();
  os.write(mtCpts.size(),"num deterministic CPTs"); os.nl();
  for (unsigned i=0;i<mtCpts.size();i++) {
    // first write the count
    os.write(i,"deterministic CPT cnt");
    os.nl();
    mtCpts[i]->write(os);
  }
  os.nl();
}


/*-
 *-----------------------------------------------------------------------
 * writeDTs
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeDTs(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Decision Trees");os.nl();
  os.write(dts.size(),"num DTS"); os.nl();
  for (unsigned i=0;i<dts.size();i++) {
    // first write the count
    os.write(i,"DTS cnt");
    os.nl();
    dts[i]->write(os);
  }
  os.nl();
}



/*-
 *-----------------------------------------------------------------------
 * writeComponents 
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeComponents(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Components");os.nl();
  unsigned used = 0;
  for (unsigned i=0;i<components.size();i++) {
    if (components[i]->emUsedBitIsSet())
      used++;
  }
  if (used != components.size())
    warning("NOTE: saving only %d used components out of a total of %d",
	    used,components.size());
  os.write(used,"num GCs"); os.nl();
  unsigned index = 0;  
  for (unsigned i=0;i<components.size();i++) {
    if (components[i]->emUsedBitIsSet()) {
      // first write the count
      os.write(index++,"GC cnt");
      os.nl();

      // next write the dimension of this component
      os.write(components[i]->dim(),"GC dim");
      os.nl();

      ////////////////////////////////////////////////////////////
      // Assume that the GC's write routine will 
      // itself write the component type
      components[i]->write(os);
    }
  }
  os.nl();
}


/*-
 *-----------------------------------------------------------------------
 * writeMixtures
 *
 *   Writes out all the mixtures in the global object. If there are
 *   any internal pre-defined mixtures in this set, then they are assumed
 *   to be kept at the beginning of the array, and are not written out
 *   (since they are automatically created anew each time the program
 *   loads). See routine loadGlobal() for more details.
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out, so moves file pointers.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeMixtures(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Mixtures of components");os.nl();
  // Leave out the first two (ie., the -2) as they are internal
  // objects. See routine finalizeParameters().
  os.write(mixtures.size()-2,"num MIXCOMPONENTS"); os.nl();

  for (unsigned i=0;i<mixtures.size()-2;i++) {
    // first write the count
    os.write(i,"MIXCOMPONENTS cnt");
    os.nl();

    // next write the dimension of this mixture
    os.write(mixtures[i]->dim(),"MG dim");
    os.nl();

    mixtures[i]->write(os);
  }

  os.nl();
}





/*-
 *-----------------------------------------------------------------------
 * writeNameCollections
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeNameCollections(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Collection of Names");os.nl();
  unsigned numToWrite = ncls.size();
  if (nclsMap.find(string(NAMED_COLLECTION_GLOBAL_NAME)) != nclsMap.end()) {
    numToWrite --;
  }
  os.write(numToWrite,"num collections"); os.nl();
  for (unsigned i=0;i<ncls.size();i++) {
    if (ncls[i]->name() == NAMED_COLLECTION_GLOBAL_NAME)
      continue;
    // write the count
    os.write(i,"NCLS cnt");
    os.nl();
    ncls[i]->write(os);
  }
  os.nl();
}





/*-
 *-----------------------------------------------------------------------
 * writeGausSwitchMixtures
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeGausSwitchMixtures(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Gaussian Switching Mixtures");os.nl();
  os.write(0,"num GausSwitchMIXTURE"); os.nl();
}




/*-
 *-----------------------------------------------------------------------
 * writeLogitSwitchMixtures
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeLogitSwitchMixtures(oDataStreamFile& os)
{
  os.nl(); os.writeComment("Logistic-Regression-based Switching Mixtures");os.nl();
  os.write(0,"num GausSwitchMIXTURE"); os.nl();
}




/*-
 *-----------------------------------------------------------------------
 * writeMlpSwitchMixtures
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      one
 *
 * Side Effects:
 *      all "used" parameters are written out.
 *
 * Results:
 *      nil.
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeMlpSwitchMixtures(oDataStreamFile& os)
{
  os.nl(); os.writeComment("MLP-based Switching Mixtures");os.nl();
  os.write(0,"num GausSwitchMIXTURE"); os.nl();
}


/*
 *  
 * ==================================================================
 *			 GROUP WRITE ROUTINES
 * ==================================================================
 *
 */



/*-
 *-----------------------------------------------------------------------
 * writeAll
 *   write out everything (all parameters) in one go.
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      files write out are written.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeAll(oDataStreamFile& os)
{
  markUsedMixtureComponents();

  // just write everything in one go.

  // write structural items
  writeDTs(os);
  writeDLinks(os);

  // then write basic numeric items.
  writeDPmfs(os);
  writeSPmfs(os);
  writeMeans(os);
  writeCovars(os);
  writeDLinkMats(os);
  writeWeightMats(os);  
  writeMdCpts(os);
  writeMsCpts(os);
  writeMtCpts(os);

  // next write definitional items
  writeComponents(os);
  writeMixtures(os);
  writeGausSwitchMixtures(os);
  writeLogitSwitchMixtures(os);
  writeMlpSwitchMixtures(os);

}


/*-
 *-----------------------------------------------------------------------
 * writeTrainable
 *   write out trainable parameters, see readTrainable for docs. 
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      files write out are written.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeTrainable(oDataStreamFile& os)
{

  markUsedMixtureComponents();

  // just write everything in one go.

  // then write basic numeric items.
  writeDPmfs(os);
  writeSPmfs(os);
  writeMeans(os);
  writeCovars(os);
  writeDLinkMats(os);
  writeWeightMats(os);  
  writeMdCpts(os);

  // next write definitional items
  writeComponents(os);
  writeMixtures(os);
  writeGausSwitchMixtures(os);
  writeLogitSwitchMixtures(os);
  writeMlpSwitchMixtures(os);
}


/*-
 *-----------------------------------------------------------------------
 * writeNonTrainable
 *   write out non-trainable parameters, see readNonTrainable for docs. 
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      files write out are written.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void 
GMParms::writeNonTrainable(oDataStreamFile& os)
{
  // just write everything in one go.
  markUsedMixtureComponents();

  // write structural items
  writeDTs(os);
  writeDLinks(os);

  writeMsCpts(os);
  writeMtCpts(os);
}



void 
GMParms::write(const char *const outputFileFormat, const char * const cppCommandOptions, const int intTag)
{
  //
  // read a file consisting of a list of 
  // 
  //      <keyword>  <filename> {ascii|binary}
  //
  // triplets. 
  // The keyword says which structure to write out,
  // The filename says where to write it.
  // And the {ascii|binary} tag says if we should write it in ascii or
  //   binary.
  //

  if (outputFileFormat == NULL)
    return;

  string keyword;
  string fileName;
  string binStatus;

  map<string,oDataStreamFile*> fileNameMap;
  char buff[2048];
  iDataStreamFile is(outputFileFormat,false,true,cppCommandOptions);

  markUsedMixtureComponents();

  while (is.readString(keyword)) {

    if (!is.readString(fileName)) {
      error("ERROR: while reading file '%s' line %d, got keyword '%s' without a filename",
	    is.fileName(),is.lineNo(),keyword.c_str());
    }
    

    bool binary_p = is.binary();
    // read binary status of file
    if (!is.readString(binStatus)) {
      error("ERROR: while reading file '%s' line %d, got keyword '%s' and filename '%s' without a binary status",
	    is.fileName(),is.lineNo(),keyword.c_str(),fileName.c_str());
    }
    if (binStatus == "ascii" || binStatus == "ASCII")
      binary_p = false;
    else if (binStatus == "binary" || binStatus == "BINARY")
      binary_p = true;
    else {
      error("ERROR: while reading file '%s' line %d, got string '%s' when expecting 'ascii'/'binary' keyword",
	    is.fileName(),is.lineNo(),binStatus.c_str());
    }

    copyStringWithTag(buff,fileName.c_str(),intTag,sizeof(buff));
    fileName = buff;

    map<string,oDataStreamFile*>::iterator it = fileNameMap.find(fileName);
    if (it == fileNameMap.end()) {
      // then this is a previously unencountered filename.
      fileNameMap[fileName] = new oDataStreamFile(fileName.c_str(),binary_p);
      it = fileNameMap.find(fileName);
    } else if (((*it).second)->binary() != binary_p) {
      // need to close the file and re-open it with the appropriate binary status, and
      // in append mode
      // delete ((*it).second);
      // fileNameMap[fileName] = new oDataStreamFile(fileName.c_str(),binary_p,true);
      // it = fileNameMap.find(fileName);
      error("ERROR: reading '%s' line %d. File '%s' had binary status = %d, now binary status = %d. Can't mix binary and ASCII files",
	    is.fileName(),is.lineNo(),fileName.c_str(),((*it).second)->binary(),binary_p);
    }

    if (keyword == "DPMF_OUT_FILE") {
      writeDPmfs(*((*it).second));

    } else if (keyword == "SPMF_OUT_FILE") {
      writeSPmfs(*((*it).second));

    } else if (keyword == "MEAN_OUT_FILE") {
      writeMeans(*((*it).second));

    } else if (keyword == "COVAR_OUT_FILE") {
      writeCovars(*((*it).second));

    } else if (keyword == "DLINK_MAT_OUT_FILE") {
      writeDLinkMats(*((*it).second));

    } else if (keyword == "DLINK_OUT_FILE") {
      writeDLinks(*((*it).second));

    } else if (keyword == "WEIGHT_MAT_OUT_FILE") {
      writeWeightMats(*((*it).second));

    } else if (keyword == "DIRICHLET_TAB_OUT_FILE") {
      writeDirichletTabs(*((*it).second));

    } else if (keyword == "DENSE_CPT_OUT_FILE") {
      writeMdCpts(*((*it).second));

    } else if (keyword == "SPARSE_CPT_OUT_FILE") {
      writeMsCpts(*((*it).second));

    } else if (keyword == "DETERMINISTIC_CPT_OUT_FILE") {
      writeMtCpts(*((*it).second));

    } else if (keyword == "DT_OUT_FILE") {
      writeDTs(*((*it).second));

    } else if (keyword == "MC_OUT_FILE") {
      writeComponents(*((*it).second));

    } else if (keyword == "MX_OUT_FILE") {
      writeMixtures(*((*it).second));

    } else if (keyword == "NAME_COLLECTION_OUT_FILE") {
      writeNameCollections(*((*it).second));

    } else if (keyword == "GSMG_OUT_FILE") {
      error("GSMG_OUT_FILE in file '%s' line %d, not implemented",
	    is.fileName(),is.lineNo());

    } else if (keyword == "LSMG_OUT_FILE") {
      error("LSMG_OUT_FILE in file '%s' line %d, not implemented",
	    is.fileName(),is.lineNo());

    } else if (keyword == "MSMG_OUT_FILE") {
      error("MSMG_OUT_FILE in file '%s' line %d, not implemented",
	    is.fileName(),is.lineNo());

    } else {
      error("ERROR: encountered unknown file type '%s' in file '%s' line %d",
	    keyword.c_str(),is.fileName(),is.lineNo());
    }
  }

  // now go through and delete all the input files
  for (map<string,oDataStreamFile*>::iterator it = fileNameMap.begin();
       it != fileNameMap.end(); it++) {
      delete ((*it).second);
  }
}





////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Misc Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * begin()
 *      do any necessary bookkeeping work with regard to the
 *      parameters in order that we properly instantiate the first example.
 * 
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      nil
 *
 * Side Effects:
 *      might change internal objects.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
GMParms::begin()
{
  for(unsigned i = 0; i<iterableDts.size(); i++) {
    iterableDts[i]->setFirstDecisionTree(firstUtterance);
    iterableDts[i]->beginIterableDT();
  }
  for ( unsigned i = 0; i < iterableLatticeAdts.size(); i++ ) {
	  iterableLatticeAdts[i]->beginIterableLattice();
  }
  for (unsigned i=0;i<dLinks.size();i++) {
    dLinks[i]->clearArrayCache();
  }
  // clear any mixtures component cache
  for (unsigned i=0;i<mixtures.size();i++) {
    mixtures[i]->emptyComponentCache();
  }
}



/*-
 *-----------------------------------------------------------------------
 * next()
 *      Move on to the next example.
 * 
 * Preconditions:
 *      begin() must have been called.
 *
 * Postconditions:
 *      nil
 *
 * Side Effects:
 *      might change internal objects.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
GMParms::next()
{
  for(unsigned i = 0; i<iterableDts.size(); i++) {
    iterableDts[i]->nextIterableDT();
  }
  for ( unsigned i = 0; i < iterableLatticeAdts.size(); i++ ) {
	  iterableLatticeAdts[i]->nextIterableLattice();
  }
  for (unsigned i=0;i<dLinks.size();i++) {
    dLinks[i]->clearArrayCache();
  }
  // clear any mixtures component cache
  for (unsigned i=0;i<mixtures.size();i++) {
    mixtures[i]->emptyComponentCache();
  }

}




/*-
 *-----------------------------------------------------------------------
 * setSegment
 *      do any necessary bookkeeping work with regard to the
 *      parameters in order that we properly move to the segment
 *      given by the argument to this function.
 *
 * Preconditions:
 *      nil
 *
 * Postconditions:
 *      nil
 *
 * Side Effects:
 *      might change internal objects.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
unsigned
GMParms::setSegment(const unsigned segmentNo)
{
  globalObservationMatrix.loadSegment(segmentNo);
  const unsigned numFrames = globalObservationMatrix.numFrames();

  for(unsigned i = 0; i<iterableDts.size(); i++) {
    iterableDts[i]->seek(segmentNo);
    iterableDts[i]->nextIterableDT();
  }
  for (unsigned i=0; i< veCpts.size(); i++) {
    veCpts[i]->setSegment(segmentNo);
    if (veCpts[i]->numFrames() != numFrames) 
      error("ERROR: number of frames in segment %d for main observation matrix is %d, but VirtualEvidenceCPT '%s' observation matrix has %d frames in that segment",segmentNo,numFrames,veCpts[i]->name().c_str(),veCpts[i]->numFrames());
  }
  // set the lattice CPT frame indices
  for (unsigned i=0; i < latticeAdts.size(); i++ ) {
	  // I cannot call iterableLatticeAdts here because
	  // there is overlap of the object and reset frame indices
	  // needs to be done to all lattices
	  if ( latticeAdts[i]->iterable() ) {
		  latticeAdts[i]->seek(segmentNo);
		  latticeAdts[i]->nextIterableLattice();
	  }
	  latticeAdts[i]->resetFrameIndices(numFrames);
  }
  for (unsigned i=0;i<dLinks.size();i++) {
    dLinks[i]->clearArrayCache();
  }
  // clear any mixtures component cache
  for (unsigned i=0;i<mixtures.size();i++) {
    mixtures[i]->emptyComponentCache();
  }

  return numFrames;
}




/*-
 *-----------------------------------------------------------------------
 * setStride()
 *      Inform all objects what the current stride will be from the
 *      observation matrix that will be used.
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      internal objects are set to use stride.
 *
 * Side Effects:
 *      internal objects are set to use stride.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
GMParms::setStride(const unsigned stride)
{
  for (unsigned i=0;i<dLinks.size();i++) {
    dLinks[i]->preCompute(stride);
  }
}




/*-
 *-----------------------------------------------------------------------
 * checkConsistentWithGlobalObservationStream()
 *      Check that the parameters are consistent
 *      with the observation matrix
 * 
 * Preconditions:
 *      global observation matrix must have been read in.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      possible error if problems occur.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
GMParms::checkConsistentWithGlobalObservationStream()
{
  if ((int)globalObservationMatrix.startSkip() < -(int)Dlinks::_globalMinLag)
    error("ERROR: a start skip of %d is invalid for a minimum dlink lag of %d\n",
	  globalObservationMatrix.startSkip(),
	  Dlinks::_globalMinLag);

  if ((int)globalObservationMatrix.endSkip() < (int)Dlinks::_globalMaxLag)
    error("ERROR: an end skip of %d is invalid for a maximum dlink lag of %d\n",
	  globalObservationMatrix.endSkip(),
	  Dlinks::_globalMaxLag);

  if ((int)globalObservationMatrix.numContinuous() <= (int)Dlinks::_globalMaxOffset)
    error("ERROR: there is a dlink ofset of value %d which is too large for the observation matrix with only %d continuous features.",
	  Dlinks::_globalMaxOffset,
	  globalObservationMatrix.numContinuous());
}


/*-
 *-----------------------------------------------------------------------
 * totalNumberParameters
 *      return the total number of parameters used by this system.
 *      NOTE: The code in this routine is very dependent on the
 *      way that the internal objects are implemented, which means
 *      that not all objects are queried for their number of parameters.
 *      For example, sparsePMFs use dense PMFs, so we only query
 *      all the dense ones to get the total number of parameters.
 * 
 * Preconditions:
 *      parameters should be read in to have non zero value
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      obvious
 *
 *-----------------------------------------------------------------------
 */
unsigned GMParms::totalNumberParameters()
{
  unsigned sum=0;

  // gets dpmfs, and covers spmfs also. This covers mixture coefficients.
  for (unsigned i=0;i<dPmfs.size();i++)
    sum += dPmfs[i]->totalNumberParameters();

  // components, gets means, vars, dlinks, weight mats
  for (unsigned i=0;i<components.size();i++)
    sum += components[i]->totalNumberParameters();

  // for discrete RVs
  for (unsigned i=0;i<mdCpts.size();i++)
    sum += mdCpts[i]->totalNumberParameters();
  for (unsigned i=0;i<msCpts.size();i++)
    sum += msCpts[i]->totalNumberParameters();
  for (unsigned i=0;i<mtCpts.size();i++)
    sum += mtCpts[i]->totalNumberParameters();
  for (unsigned i=0;i<ngramCpts.size();i++)
    sum += ngramCpts[i]->totalNumberParameters();
  for (unsigned i=0;i<fngramCpts.size();i++)
    sum += fngramCpts[i]->totalNumberParameters();
  for (unsigned i=0;i<veCpts.size();i++)
    sum += veCpts[i]->totalNumberParameters();
  return sum;

}




/*-
 *-----------------------------------------------------------------------
 * markUsedMixtureComponents()
 *      Remove all mixture component objects that are not used by anyone
 *      (possibly because of vanishing).
 * 
 * Preconditions:
 *      parameters should be read in to have non zero value
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      will remove all paramters not used.
 *
 * Results:
 *      obvious
 *
 *-----------------------------------------------------------------------
 */
void GMParms::markUsedMixtureComponents()
{

  ///////////////////////////////////////////////////
  // First, go through *all* component related parameters
  // and mark as not used.
  for (unsigned i=0;i<means.size();i++)
    means[i]->recursivelyClearUsedBit();
  for (unsigned i=0;i<covars.size();i++)
    covars[i]->recursivelyClearUsedBit();
  for (unsigned i=0;i<dLinkMats.size();i++)
    dLinkMats[i]->recursivelyClearUsedBit();
  for (unsigned i=0;i<weightMats.size();i++)
    weightMats[i]->recursivelyClearUsedBit();
  for (unsigned i=0;i<components.size();i++)
    components[i]->recursivelyClearUsedBit();
  for (unsigned i=0;i<mixtures.size();i++)
    mixtures[i]->recursivelyClearUsedBit();

  ///////////////////////////////////////////////
  // Now, set only those bits that are used
  // by some mixture.
  for (unsigned i=0;i<mixtures.size();i++)
    mixtures[i]->recursivelySetUsedBit();

}




/*-
 *-----------------------------------------------------------------------
 *  markObjectsToNotTrain(fn)
 *      Using the file given by file name, mark all objects
 *      listed in that file so that they are not trained by EM. 
 *  The format of the file is a list of pairs, where
 *  the first element in the pair specifies an object type
 *  and the 2nd gives an object name.
 *
 *
 * 
 * Preconditions:
 *      parameters should be read in to have non zero value
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      will set all objects listed in file so that they are not trained.
 *
 * Results:
 *      il
 *
 *-----------------------------------------------------------------------
 */
void GMParms::markObjectsToNotTrain(const char*const fileName,
				    const char *const cppOptions)
{
  if (fileName == NULL)
    return;

  iDataStreamFile is(fileName, // name
			false,    // not binary 
			true,     // run cpp
			cppOptions);

  string objType;
  string objName;
  while (is.readString(objType)) {
    if (!is.readString(objName)) {
      error("ERROR: while reading file '%s' line %d, got object type '%s' without an object name",
	    is.fileName(),is.lineNo(),
	    objType.c_str());
    }

#define EMCLEARAMTRAININGBIT_CODE(mapName,name) \
      if (objName == "*") { \
         for (unsigned i=0;i<name.size();i++)  \
            name[i]->emClearAmTrainingBit(); \
      } else { \
	if ((it=mapName.find(objName)) == mapName.end()) \
	  error("ERROR: can't find object '%s' of type '%s' listed in file '%s' line %d of objects to not train.", \
	      objName.c_str(), \
	      objType.c_str(), \
	      is.fileName(),is.lineNo()); \
        name[(*it).second]->emClearAmTrainingBit(); }


    ObjectMapType::iterator it;
    if (objType == "DPMF") {
      EMCLEARAMTRAININGBIT_CODE(dPmfsMap,dPmfs);
    } else if (objType == "SPMF") {
      EMCLEARAMTRAININGBIT_CODE(sPmfsMap,sPmfs);
    } else if (objType == "MEAN") {
      EMCLEARAMTRAININGBIT_CODE(meansMap,means);
    } else if (objType == "COVAR") {
      EMCLEARAMTRAININGBIT_CODE(covarsMap,covars);
    } else if (objType == "DLINKMAT") {
      EMCLEARAMTRAININGBIT_CODE(dLinkMatsMap,dLinkMats);
    } else if (objType == "WEIGHTMAT") {
      EMCLEARAMTRAININGBIT_CODE(weightMatsMap,weightMats);
    } else if (objType == "COMPONENT") {
      EMCLEARAMTRAININGBIT_CODE(componentsMap,components);
    } else if (objType == "DENSECPT") {
      EMCLEARAMTRAININGBIT_CODE(mdCptsMap,mdCpts);
    } else if (objType == "SPARSECPT") {
      EMCLEARAMTRAININGBIT_CODE(msCptsMap,msCpts);
    } else if (objType == "DETERMINISTICCPT") {
      EMCLEARAMTRAININGBIT_CODE(mtCptsMap,mtCpts);
    } else if (objType == "MIXTURE") {
      EMCLEARAMTRAININGBIT_CODE(mixturesMap,mixtures);
    } else {
      error("ERROR: bad object type name '%s' in file '%s' line %d of objects to not train.",
	    objType.c_str(),
	    is.fileName(),is.lineNo());
    }

  }
}



 /*-
  *-----------------------------------------------------------------------
  * makeRandom
  *      make everyone random
  * 
  * Preconditions:
  *      all must be read in.
  *
  * Postconditions:
  *      What is true after the function is called.
  *
  * Side Effects:
  *      all variable parameters are changed
  *
  * Results:
  *      nil
  *
  *-----------------------------------------------------------------------
  */
 void
 GMParms::makeRandom()
 {
   // go through all EMable objects possibly
   // used by any RV and make the call

   // first do the basic objects
   // which also does everything for continuous RVs
   for (unsigned i=0;i<dPmfs.size();i++)
     dPmfs[i]->makeRandom();
   for (unsigned i=0;i<sPmfs.size();i++)
     sPmfs[i]->makeRandom();
   for (unsigned i=0;i<means.size();i++)
     means[i]->makeRandom();
   for (unsigned i=0;i<covars.size();i++)
     covars[i]->makeRandom();
   for (unsigned i=0;i<dLinkMats.size();i++)
     dLinkMats[i]->makeRandom();
   for (unsigned i=0;i<weightMats.size();i++)
     weightMats[i]->makeRandom();

   // components
   for (unsigned i=0;i<components.size();i++)
     components[i]->makeRandom();

   // for discrete RVs
   for (unsigned i=0;i<mdCpts.size();i++)
     mdCpts[i]->makeRandom();
   for (unsigned i=0;i<msCpts.size();i++)
     msCpts[i]->makeRandom();
   for (unsigned i=0;i<mtCpts.size();i++)
     mtCpts[i]->makeRandom();

 }



void
GMParms::makeUniform()
{
  // go through all EMable objects possibly
  // used by any RV and make the call

   // first do the basic objects
   // which also does everything for continuous RVs
  for (unsigned i=0;i<dPmfs.size();i++)
    dPmfs[i]->makeUniform();
  for (unsigned i=0;i<sPmfs.size();i++)
    sPmfs[i]->makeUniform();
  for (unsigned i=0;i<means.size();i++)
    means[i]->makeUniform();
  for (unsigned i=0;i<covars.size();i++)
    covars[i]->makeUniform();
  for (unsigned i=0;i<dLinkMats.size();i++)
    dLinkMats[i]->makeUniform();
  for (unsigned i=0;i<weightMats.size();i++)
    weightMats[i]->makeUniform();

   // components
  for (unsigned i=0;i<components.size();i++)
    components[i]->makeUniform();

   // for discrete RVs
  for (unsigned i=0;i<mdCpts.size();i++)
    mdCpts[i]->makeUniform();
  for (unsigned i=0;i<msCpts.size();i++)
    msCpts[i]->makeUniform();
  for (unsigned i=0;i<mtCpts.size();i++)
    mtCpts[i]->makeUniform();

}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        EM Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * emEndIteration()
 *      Call emEndIteration() for all top level objects (top level
 *      objects are those who will not be pointed to by other objects).
 *
 * Preconditions:
 *      em should be running.
 *
 * Postconditions:
 *      EM iteration is finished.
 *
 * Side Effects:
 *      Changes almost all of the parameter objects.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
GMParms::emEndIteration()
{
  // go through all EMable objects possibly
  // used by any RV and make the call

   /////////////////////////////////////////////////////////////////
   // First, do the mixtures. This will recursively call
   // this for all mean-like objects, covariance-like objects, and so
   // on, so there is no need to do those here.  
   //
   // Mean-like objects include:
   //
   //       components    
   //       means
   //       linCondMeans
   //       nonLinCondMeans
   //       ...
   // 
   // Variance-like objects include
   // 
   //       covars
   // 
  for (unsigned i=0;i<mixtures.size();i++)
    mixtures[i]->emEndIteration();

   //////////////////////////////////////////////////////////////
   // We don't do the following code  here because those objects use
   // ref counts needed for sharing. The components will
   // themselves end the EM iteration for those objects, and will call
   // it the appropriate number of times so that the ref counts are set
   // properly.
   //
   //    for (unsigned i=0;i<means.size();i++)
   //      means[i]->emEndIteration();
   //    for (unsigned i=0;i<covars.size();i++)
   //      covars[i]->emEndIteration();
   //    for (unsigned i=0;i<dLinkMats.size();i++)
   //      dLinkMats[i]->emEndIteration();
   //    for (unsigned i=0;i<weightMats.size();i++)
   //      weightMats[i]->emEndIteration();
   //
   //////////////////////////////////////////////////////////////

   // finish up the spmfs and dpmfs, some of which
   // might already have been ended by the
   // above code via the mixture component densities. 
  for (unsigned i=0;i<dPmfs.size();i++)
    dPmfs[i]->emEndIteration();
  for (unsigned i=0;i<sPmfs.size();i++)
    sPmfs[i]->emEndIteration();

   // finally, end all the CPT. for discrete RVs
  for (unsigned i=0;i<mdCpts.size();i++)
    mdCpts[i]->emEndIteration();
  for (unsigned i=0;i<msCpts.size();i++)
    msCpts[i]->emEndIteration();
  for (unsigned i=0;i<mtCpts.size();i++)
    mtCpts[i]->emEndIteration();


#if 0
  // uncomment this when the objects get written.
  for (unsigned i=0;i<gausSwitchingMixtures.size();i++)
    gausSwitchingMixtures[i]->emEndIteration();
  for (unsigned i=0;i<logitSwitchingMixtures.size();i++)
    logitSwitchingMixtures[i]->emEndIteration();
  for (unsigned i=0;i<mlpSwitchingMixtures.size();i++)
    mlpSwitchingMixtures[i]->emEndIteration();
#endif

}




/*-
 *-----------------------------------------------------------------------
 * emSwapCurAndNew()
 *      Call emSwapCurAndNew() for all top level objects (top level
 *      objects are those who will not be pointed to by other objects).
 *
 * Preconditions:
 *      EM should not be running.
 *
 * Postconditions:
 *      EM iteration is finished.
 *
 * Side Effects:
 *      Changes almost all of the parameter objects.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
GMParms::emSwapCurAndNew()
{
  // go through all EMable objects possibly
  // used by any RV and make the call

   // for continuous RVs, this will recursively
   // call swap for all sub objects.
  for (unsigned i=0;i<mixtures.size();i++)
    mixtures[i]->emSwapCurAndNew();

   // make sure that all dpmfs and spmfs have been swapped.
  for (unsigned i=0;i<dPmfs.size();i++)
    dPmfs[i]->emSwapCurAndNew();
  for (unsigned i=0;i<sPmfs.size();i++)
    sPmfs[i]->emSwapCurAndNew();

  // Don't swap these since it should have
  // been done by the mixture swapping above.
  //    for (unsigned i=0;i<means.size();i++)
  //      means[i]->emSwapCurAndNew();
  //    for (unsigned i=0;i<covars.size();i++)
  //      covars[i]->emSwapCurAndNew();
  //    for (unsigned i=0;i<dLinkMats.size();i++)
  //      dLinkMats[i]->emSwapCurAndNew();
  //    for (unsigned i=0;i<weightMats.size();i++)
  //      weightMats[i]->emSwapCurAndNew();
  //    for (unsigned i=0;i<components.size();i++)
  //      components[i]->emSwapCurAndNew();

   // for discrete RVs
  for (unsigned i=0;i<mdCpts.size();i++)
    mdCpts[i]->emSwapCurAndNew();
  for (unsigned i=0;i<msCpts.size();i++)
    msCpts[i]->emSwapCurAndNew();
  for (unsigned i=0;i<mtCpts.size();i++)
    mtCpts[i]->emSwapCurAndNew();


#if 0
  // uncomment this when the objects get written.
  for (unsigned i=0;i<gausSwitchingMixtures.size();i++)
    gausSwitchingMixtures[i]->emSwapCurAndNew();
  for (unsigned i=0;i<logitSwitchingMixtures.size();i++)
    logitSwitchingMixtures[i]->emSwapCurAndNew();
  for (unsigned i=0;i<mlpSwitchingMixtures.size();i++)
    mlpSwitchingMixtures[i]->emSwapCurAndNew();
#endif

  // clear out the cloning maps, under the assumption
  // that all clones become 'real' objects.
  MixtureCommon::vanishingComponentSet.clear();
  MixtureCommon::splittingComponentSet.clear();

  MixtureCommon::meanCloneMap.clear();
  MixtureCommon::dLinkMatCloneMap.clear();
  MixtureCommon::diagCovarCloneMap.clear();
  MixtureCommon::mcCloneMap.clear();

}

void
GMParms::emStoreAccumulators(oDataStreamFile& ofile)
{
  // first do the basic objects
  for (unsigned i=0;i<dPmfs.size();i++)
    dPmfs[i]->emStoreAccumulators(ofile);
  for (unsigned i=0;i<sPmfs.size();i++)
    sPmfs[i]->emStoreAccumulators(ofile);
  for (unsigned i=0;i<means.size();i++)
    means[i]->emStoreAccumulators(ofile);
  for (unsigned i=0;i<covars.size();i++)
    covars[i]->emStoreAccumulators(ofile);
  for (unsigned i=0;i<dLinkMats.size();i++)
    dLinkMats[i]->emStoreAccumulators(ofile);
  for (unsigned i=0;i<weightMats.size();i++)
    weightMats[i]->emStoreAccumulators(ofile);

   // components
  for (unsigned i=0;i<components.size();i++)
    components[i]->emStoreAccumulators(ofile);

   // for discrete RVs
  for (unsigned i=0;i<mdCpts.size();i++)
    mdCpts[i]->emStoreAccumulators(ofile);
  for (unsigned i=0;i<msCpts.size();i++)
    msCpts[i]->emStoreAccumulators(ofile);
  for (unsigned i=0;i<mtCpts.size();i++)
    mtCpts[i]->emStoreAccumulators(ofile);

   // for continuous RVs
  for (unsigned i=0;i<mixtures.size();i++)
    mixtures[i]->emStoreAccumulators(ofile);
#if 0
  // uncomment this when the objects get written.
  for (unsigned i=0;i<gausSwitchingMixtures.size();i++)
    gausSwitchingMixtures[i]->emStoreAccumulators(ofile);
  for (unsigned i=0;i<logitSwitchingMixtures.size();i++)
    logitSwitchingMixtures[i]->emStoreAccumulators(ofile);
  for (unsigned i=0;i<mlpSwitchingMixtures.size();i++)
    mlpSwitchingMixtures[i]->emStoreAccumulators(ofile);
#endif

}

void
GMParms::emLoadAccumulators(iDataStreamFile& ifile)
{
  /////////////////////////////////////////////////////////////
  // First, make sure the EM iterations have been
  // started. Like in the emEndIteration function above,
  // we do not call this for all objects. This call
  // ensure that all objects internal EM data structures
  // have been allocated, so that they can be accumulated to.

  // components
  for (unsigned i=0;i<components.size();i++)
    components[i]->emStartIteration();

  // do the basic discrete objects
  for (unsigned i=0;i<dPmfs.size();i++)
    dPmfs[i]->emStartIteration();
  for (unsigned i=0;i<sPmfs.size();i++)
    sPmfs[i]->emStartIteration();

  // for discrete RVs
  for (unsigned i=0;i<mdCpts.size();i++)
    mdCpts[i]->emStartIteration();
  for (unsigned i=0;i<msCpts.size();i++)
    msCpts[i]->emStartIteration();
  for (unsigned i=0;i<mtCpts.size();i++)
    mtCpts[i]->emStartIteration();

  // for continuous RVs
  for (unsigned i=0;i<mixtures.size();i++)
    mixtures[i]->emStartIteration();
  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Now, we actually load the accumulators.

  // first do the basic objects
  for (unsigned i=0;i<dPmfs.size();i++)
    dPmfs[i]->emLoadAccumulators(ifile);
  for (unsigned i=0;i<sPmfs.size();i++)
    sPmfs[i]->emLoadAccumulators(ifile);
  for (unsigned i=0;i<means.size();i++)
    means[i]->emLoadAccumulators(ifile);
  for (unsigned i=0;i<covars.size();i++)
    covars[i]->emLoadAccumulators(ifile);
  for (unsigned i=0;i<dLinkMats.size();i++)
    dLinkMats[i]->emLoadAccumulators(ifile);
  for (unsigned i=0;i<weightMats.size();i++)
    weightMats[i]->emLoadAccumulators(ifile);

  // components
  for (unsigned i=0;i<components.size();i++)
    components[i]->emLoadAccumulators(ifile);

  // for discrete RVs
  for (unsigned i=0;i<mdCpts.size();i++)
    mdCpts[i]->emLoadAccumulators(ifile);
  for (unsigned i=0;i<msCpts.size();i++)
    msCpts[i]->emLoadAccumulators(ifile);
  for (unsigned i=0;i<mtCpts.size();i++)
    mtCpts[i]->emLoadAccumulators(ifile);

  // for continuous RVs
  for (unsigned i=0;i<mixtures.size();i++)
    mixtures[i]->emLoadAccumulators(ifile);
#if 0
  // uncomment this when the objects get written.
  for (unsigned i=0;i<gausSwitchingMixtures.size();i++)
    gausSwitchingMixtures[i]->emLoadAccumulators(ifile);
  for (unsigned i=0;i<logitSwitchingMixtures.size();i++)
    logitSwitchingMixtures[i]->emLoadAccumulators(ifile);
  for (unsigned i=0;i<mlpSwitchingMixtures.size();i++)
    mlpSwitchingMixtures[i]->emLoadAccumulators(ifile);
#endif
  /////////////////////////////////////////////////////////////
}

void
GMParms::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  // first do the basic objects
  for (unsigned i=0;i<dPmfs.size();i++)
    dPmfs[i]->emAccumulateAccumulators(ifile);
  for (unsigned i=0;i<sPmfs.size();i++)
    sPmfs[i]->emAccumulateAccumulators(ifile);
  for (unsigned i=0;i<means.size();i++)
    means[i]->emAccumulateAccumulators(ifile);
  for (unsigned i=0;i<covars.size();i++)
    covars[i]->emAccumulateAccumulators(ifile);
  for (unsigned i=0;i<dLinkMats.size();i++)
    dLinkMats[i]->emAccumulateAccumulators(ifile);
  for (unsigned i=0;i<weightMats.size();i++)
    weightMats[i]->emAccumulateAccumulators(ifile);

  // components
  for (unsigned i=0;i<components.size();i++)
    components[i]->emAccumulateAccumulators(ifile);

  // for discrete RVs
  for (unsigned i=0;i<mdCpts.size();i++)
    mdCpts[i]->emAccumulateAccumulators(ifile);
  for (unsigned i=0;i<msCpts.size();i++)
    msCpts[i]->emAccumulateAccumulators(ifile);
  for (unsigned i=0;i<mtCpts.size();i++)
    mtCpts[i]->emAccumulateAccumulators(ifile);

  // for continuous RVs
  for (unsigned i=0;i<mixtures.size();i++)
    mixtures[i]->emAccumulateAccumulators(ifile);
#if 0
  // uncomment this when the objects get written.
  for (unsigned i=0;i<gausSwitchingMixtures.size();i++)
    gausSwitchingMixtures[i]->emAccumulateAccumulators(ifile);
  for (unsigned i=0;i<logitSwitchingMixtures.size();i++)
    logitSwitchingMixtures[i]->emAccumulateAccumulators(ifile);
  for (unsigned i=0;i<mlpSwitchingMixtures.size();i++)
    mlpSwitchingMixtures[i]->emAccumulateAccumulators(ifile);
#endif
}


/*-
 *-----------------------------------------------------------------------
 * setFirstUtterance
 *      Set the index of the first utterance so that the parameters to
 *      identify which per-utterance decision trees to use. 
 *
 * Preconditions:
 *      none 
 *
 * Postconditions:
 *      firstUtterance is set 
 *
 * Side Effects:
 *      none 
 *
 * Results:
 *      nil
 *-----------------------------------------------------------------------
 */
void
GMParms::setFirstUtterance(
  unsigned first_index 
  )
{
  firstUtterance = first_index;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////





#ifdef MAIN

#include "rand.h"
#include "GMTK_ObservationMatrix.h"

RAND rnd(false);
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;



int
main()
{

  // read in a file of decision trees
  {
  iDataStreamFile isdt("ex.means",false);
  GM_Parms.readMeans(isdt);

  oDataStreamFile os("-");
  GM_Parms.writeMeans(os);  
  }
  {
  iDataStreamFile isdt("ex.covars",false);
  GM_Parms.readCovars(isdt);

  oDataStreamFile os("-");
  GM_Parms.writeCovars(os);  
  }



#if 0
  // read in basic structures
  iDataStreamFile is("dataFiles/test1.gmb",false);

  // write both out again to stdout (i.e., "-")
  oDataStreamFile os("-");
  GM_Parms.writeDTs(os);
#endif

}


#endif
