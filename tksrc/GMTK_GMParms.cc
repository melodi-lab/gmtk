/*-
 * GMTK_World.cc
 *        All aspects about a GM.
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
#include <ieeefp.h>
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
#include "GMTK_PackedSparseRealMatrix.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_WeightMatrix.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"

#include "GMTK_GaussianComponent.h"
#include "GMTK_DiagGaussian.h"

#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_MixGaussians.h"


VCID("$Header$");


////////////////////////////////
// Magic String definitions
#define MAGIC_DT_FILE "GMTK_DT_FILE"
#define MAGIC_PRM_FILE "GMTK_PRM_FILE"


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////

GMParms::GMParms()
{}

void GMParms::readDPmfs(iDataStreamFile& is,bool reset)
{

}


void 
GMParms::readBasic(iDataStreamFile& is)
{

  unsigned num;
  unsigned cnt;
  string str;

  is.read(str,"GMTK_GMParms::readBasic, magic");
  if (str != MAGIC_PRM_FILE)
    error("GMTK_GMParms::readBasic. Expecting basic param file, got (%s) in file (%s)",str.c_str(),is.fileName());


  is.read(num,"GMTK_GMParms::readBasic, dpmfs");
  if (num < 0) error("GMTK_GMParms::readBasic num dpmfs = %d",num);
  dPmfs.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt dpmfs");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,dpmfs, out of order count",cnt);
    dPmfs[i] = new Dense1DPMF;
    dPmfs[i]->read(is);
    if (dPmfsMap.find(dPmfs[i]->name()) != dPmfsMap.end())
      error("GMTK_GMParms::readBasic,dpmfs, multiple use of name '%s'",
	    dPmfs[i]->name().c_str());
    dPmfsMap[dPmfs[i]->name()] = i;
  }


  is.read(num,"GMTK_GMParms::readBasic, spmfs");
  if (num < 0) error("GMTK_GMParms::readBasic num spmfs = %d",num);
  sPmfs.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt spmfs");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,spmfs, out of order count",cnt);
    sPmfs[i] = new Sparse1DPMF;
    sPmfs[i]->read(is);
    sPmfsMap[sPmfs[i]->name()] = i;
  }

  is.read(num,"GMTK_GMParms::readBasic, means");
  if (num < 0) error("GMTK_GMParms::readBasic num means = %d",num);
  means.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt means");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,means,  out of order count",cnt);
    means[i] = new MeanVector;
    means[i]->read(is);
    meansMap[means[i]->name()] = i;
  }

  is.read(num,"GMTK_GMParms::readBasic, covars");
  if (num < 0) error("GMTK_GMParms::readBasic num covars = %d",num);
  covars.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt covars");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,covars, out of order count",cnt);
    covars[i] = new DiagCovarVector;
    covars[i]->read(is);
    covarsMap[covars[i]->name()] = i;
  }

  is.read(num,"GMTK_GMParms::readBasic, DlinkMatrix");
  if (num < 0) error("GMTK_GMParms::readBasic num DlinkMatrix = %d",num);
  dLinkMats.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt dlinks");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,dlinks, out of order count",cnt);
    dLinkMats[i] = new DlinkMatrix;
    dLinkMats[i]->read(is);
    dLinkMatsMap[dLinkMats[i]->name()] = i;
  }

  is.read(num,"GMTK_GMParms::readBasic, WeightMatrix");
  if (num < 0) error("GMTK_GMParms::readBasic num WeightMatrix = %d",num);
  weightMats.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt weights");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic, weights, out of order count",cnt);
    weightMats[i] = new WeightMatrix;
    weightMats[i]->read(is);
    weightMatsMap[weightMats[i]->name()] = i;
  }


  is.read(num,"GMTK_GMParms::readBasic, MDCPT");
  if (num < 0) error("GMTK_GMParms::readBasic num MDCPT = %d",num);
  mdCpts.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt MDCPTs");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic, MDCPTs, out of order count",cnt);
    mdCpts[i] = new MDCPT;
    mdCpts[i]->read(is);
    mdCptsMap[mdCpts[i]->name()] = i;
  }

  is.read(num,"GMTK_GMParms::readBasic, MSCPT");
  if (num < 0) error("GMTK_GMParms::readBasic num MSCPT = %d",num);
  msCpts.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt MSCPTS");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic, MSCPTs, out of order count",cnt);
    msCpts[i] = new MSCPT;
    msCpts[i]->read(is);
    msCptsMap[msCpts[i]->name()] = i;
  }


  is.read(num,"GMTK_GMParms::readBasic, MTCPT");
  if (num < 0) error("GMTK_GMParms::readBasic num MTCPT = %d",num);
  mtCpts.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt MTCPTS");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic, MTCPTs, out of order count",cnt);
    mtCpts[i] = new MTCPT;
    mtCpts[i]->read(is);
    mtCptsMap[mtCpts[i]->name()] = i;
  }

}

void 
GMParms::writeBasic(oDataStreamFile& os)
{

  os.write(MAGIC_PRM_FILE,"GMTK_GMParms::writeBasic, magic");
  os.nl();

  os.nl(); os.writeComment("dense PMFs");os.nl();
  os.write(dPmfs.size(),"GMTK_GMParms::writeBasic, dpmfs");
  os.nl();
  for (unsigned i=0;i<dPmfs.size();i++) {
    os.write(i);
    dPmfs[i]->write(os);
  }
  os.nl();

  os.nl(); os.writeComment("sparse PMFs");os.nl();
  os.write(sPmfs.size(),"GMTK_GMParms::writeBasic, spmfs");
  os.nl();
  for (unsigned i=0;i<sPmfs.size();i++) {
    os.write(i);
    sPmfs[i]->write(os);
  }
  os.nl();


  os.nl(); os.write(means.size(),"GMTK_GMParms::writeBasic, means");
  os.nl();
  for (unsigned i=0;i<means.size();i++) {
    os.write(i);
    means[i]->write(os);
  }
  os.nl();

  os.nl(); os.writeComment("covars");os.nl();
  os.write(covars.size(),"GMTK_GMParms::writeBasic, covars");
  os.nl();
  for (unsigned i=0;i<covars.size();i++) {
    os.write(i);
    covars[i]->write(os);
  }
  os.nl();

  os.nl(); os.writeComment("dlink matrices");os.nl();
  os.write(dLinkMats.size(),"GMTK_GMParms::writeBasic, DlinkMatrix");
  os.nl();
  for (unsigned i=0;i<dLinkMats.size();i++) {
    os.write(i); os.nl();
    dLinkMats[i]->write(os);
  }
  os.nl();

  os.nl(); os.writeComment("weight matrices");os.nl();
  os.write(weightMats.size(),"GMTK_GMParms::writeBasic, WeightMatrix");
  os.nl();
  for (unsigned i=0;i<weightMats.size();i++) {
    os.write(i); os.nl();
    weightMats[i]->write(os);
  }
  os.nl();

  os.nl(); os.writeComment("MDCPTs");os.nl();
  os.write(mdCpts.size(),"GMTK_GMParms::writeBasic, MDCPT");
  os.nl();
  for (unsigned i=0;i<mdCpts.size();i++) {
    os.write(i); os.nl();
    mdCpts[i]->write(os);
  }
  os.nl();

  os.nl();  os.writeComment("MSCPTs");os.nl();
  os.write(msCpts.size(),"GMTK_GMParms::writeBasic, MSCPT");
  os.nl();
  for (unsigned i=0;i<msCpts.size();i++) {
    os.write(i); os.nl();
    msCpts[i]->write(os);
  }
  os.nl();


  os.nl(); os.writeComment("MTCPTs");os.nl();
  os.write(mtCpts.size(),"GMTK_GMParms::writeBasic, MTCPT");
  os.nl();
  for (unsigned i=0;i<mtCpts.size();i++) {
    os.write(i); os.nl();
    mtCpts[i]->write(os);
  }
  os.nl();

}



void 
GMParms::readDTs(iDataStreamFile& is)
{
  unsigned num;
  unsigned cnt;
  char *str;

  is.read(str,"GMTK_GMParms::readDTs, magic");
  if (strcmp(str,MAGIC_DT_FILE))
    error("GMTK_GMParms::readDTs. Expecting DT file, got (%s)",str);
  delete [] str;

  is.read(num,"GMTK_GMParms::readDTs, dpmfs");
  if (num < 0) error("GMTK_GMParms::readDTs num dpmfs = %d",num);
  dts.resize(num);
  for (unsigned i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readDTs, cnt");
    if (cnt != i) 
      error("GMTK_GMParms::readDTs, out of order count",cnt);

    dts[i] = new RngDecisionTree<unsigned>;
    dts[i]->read(is);
    dtsMap[dts[i]->name()] = i;
  }
}



void 
GMParms::writeDTs(oDataStreamFile& os)
{
  os.write(MAGIC_DT_FILE,"GMTK_GMParms::writeDTs, magic");
  os.nl();
  os.write(dts.size(),"GMTK_GMParms::writeDTs, dpmfs");
  os.nl();
  for (unsigned i=0;i<dts.size();i++) {
    os.write(i); os.nl();
    dts[i]->write(os);
  }
}


void 
GMParms::readGaussianComponents(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;

  unsigned start = 0;
  is.read(num,"GMTK_GMParms::readGaussianComponents");
  if (num < 0) error("GMTK_GMParms::readGaussianComponents num = %d",num);
  if (reset) {
    start = 0;
    gaussianComponents.resize(num);
  } else {
    start = gaussianComponents.size();
    gaussianComponents.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    GaussianComponent* gc;

    is.read(cnt,"GMTK_GMParms::readGaussianComponents, cnt");
    if (cnt != i) 
      error("GMTK_GMParms::readGaussianComponents, out of order count",cnt);

    // next read the dimension of this Gaussian
    int dim;
    is.read(dim,"GMTK_GMParms::readGaussianComponents, dim");

    // read the Gaussian type
    int t;
    is.read(t,"GMTK_GMParms::readGaussianComponents, type");
    if (t == GaussianComponent::Diag) {
      gc = new DiagGaussian(dim);
    } else if (t == GaussianComponent::LinMeanCondDiag) {
      error("LinMeanCondDiag not implemented");
      // gc = new LinMeanCondDiagGaussian(dim);
    } else if (t == GaussianComponent::NLinMeanCondDiag) {
      error("NLinMeanCondDiag not implemented");
      // gc = new NLinMeanCondDiagGaussian(dim);
    } else {
      error("Error: unknown gaussian component type in file");
    }
    gc->read(is);
    gaussianComponents[i+start] = gc;
    gaussianComponentsMap[gc->name()] = i+start;
  }
}




void 
GMParms::readGaussianMixtures(iDataStreamFile& is, bool reset)
{
  unsigned num;
  unsigned cnt;

  unsigned start = 0;
  is.read(num,"GMTK_GMParms::readGaussianMixtures");
  if (num < 0) error("GMTK_GMParms::readGaussianMixtures num = %d",num);
  if (reset) {
    start = 0;
    mixGaussians.resize(num);
  } else {
    start = mixGaussians.size();
    mixGaussians.resize(start+num);
  }
  for (unsigned i=0;i<num;i++) {
    // first read the count
    MixGaussiansCommon* gm;

    is.read(cnt,"GMTK_GMParms::readGaussianMixtures, cnt");
    if (cnt != i) 
      error("GMTK_GMParms::readGaussianMixtures, out of order count",cnt);

    // next read the dimension of this Gaussian
    int dim;
    is.read(dim,"GMTK_GMParms::readGaussianMixtures, dim");

    // read the Gaussian type
    int t;
    is.read(t,"GMTK_GMParms::readGaussianMixtures, type");
    if (t == MixGaussiansCommon::ci_mixGaussian) {
      gm = new MixGaussians(dim);
    } else if (t == MixGaussiansCommon::ci_gausSwitchMixGaussian) {
      error("GausSwitchMix not implemented");
      // gm = new GausSwitchingMixGaussians();
    } else if (t == MixGaussiansCommon::ci_logitSwitchMixGaussian) {
      error("LogitSwitchMix not implemented");
      // gm = new LogitSwitchingMixGaussians();      
    } else if (t == MixGaussiansCommon::ci_mlpSwitchMixGaussian) {
      error("MlpSwitchMix not implemented");
      // gm = new MLPSwitchingMixGaussians(dim);
    } else {
      error("Error: unknown gaussian mixture type in file");
    }
    gm->read(is);
    mixGaussians[i+start] = gm;
    mixGaussiansMap[gm->name()] = i+start;
  }
}






#ifdef MAIN

GMParms GM_Parms;

int
main()
{

  // read in a file of decision trees
  iDataStreamFile isdt("dataFiles/foo.dt",false);
  GM_Parms.readDTs(isdt);

  // read in basic structures
  iDataStreamFile is("dataFiles/test1.gmb",false);
  GM_Parms.readBasic(is);

  // write both out again to stdout (i.e., "-")
  oDataStreamFile os("-");
  GM_Parms.writeDTs(os);
  GM_Parms.writeBasic(os);

}


#endif
