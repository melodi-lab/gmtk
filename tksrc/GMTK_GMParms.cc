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



void 
GMParms::readBasic(iDataStreamFile& is)
{

  int num;
  int cnt;
  char *str;

  is.read(str,"GMTK_GMParms::readBasic, magic");
  if (strcmp(str,MAGIC_PRM_FILE))
    error("GMTK_GMParms::readBasic. Expecting basic param file, got (%s) in file (%s)",str,is.fileName());
  delete [] str;


  is.read(num,"GMTK_GMParms::readBasic, dpmfs");
  if (num < 0) error("GMTK_GMParms::readBasic num dpmfs = %d",num);
  dPmfs.resize(num);
  for (int i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt dpmfs");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,dpmfs, out of order count",cnt);
    dPmfs[i] = new Dense1DPMF;
    dPmfs[i]->read(is);
  }


  is.read(num,"GMTK_GMParms::readBasic, spmfs");
  if (num < 0) error("GMTK_GMParms::readBasic num spmfs = %d",num);
  sPmfs.resize(num);
  for (int i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt spmfs");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,spmfs, out of order count",cnt);
    sPmfs[i] = new Sparse1DPMF;
    sPmfs[i]->read(is);
  }

  is.read(num,"GMTK_GMParms::readBasic, means");
  if (num < 0) error("GMTK_GMParms::readBasic num means = %d",num);
  means.resize(num);
  for (int i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt means");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,means,  out of order count",cnt);
    means[i] = new MeanVector;
    means[i]->read(is);
  }

  is.read(num,"GMTK_GMParms::readBasic, covars");
  if (num < 0) error("GMTK_GMParms::readBasic num covars = %d",num);
  covars.resize(num);
  for (int i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt covars");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,covars, out of order count",cnt);
    covars[i] = new DiagCovarVector;
    covars[i]->read(is);
  }

  is.read(num,"GMTK_GMParms::readBasic, DlinkMatrix");
  if (num < 0) error("GMTK_GMParms::readBasic num DlinkMatrix = %d",num);
  dLinkMats.resize(num);
  for (int i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt dlinks");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic,dlinks, out of order count",cnt);
    dLinkMats[i] = new DlinkMatrix;
    dLinkMats[i]->read(is);
  }

  is.read(num,"GMTK_GMParms::readBasic, WeightMatrix");
  if (num < 0) error("GMTK_GMParms::readBasic num WeightMatrix = %d",num);
  weightMats.resize(num);
  for (int i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt weights");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic, weights, out of order count",cnt);
    weightMats[i] = new WeightMatrix;
    weightMats[i]->read(is);
  }


  is.read(num,"GMTK_GMParms::readBasic, MDCPT");
  if (num < 0) error("GMTK_GMParms::readBasic num MDCPT = %d",num);
  mdCpts.resize(num);
  for (int i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt MDCPTs");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic, MDCPTs, out of order count",cnt);
    mdCpts[i] = new MDCPT;
    mdCpts[i]->read(is);
  }

  is.read(num,"GMTK_GMParms::readBasic, MSCPT");
  if (num < 0) error("GMTK_GMParms::readBasic num MSCPT = %d",num);
  msCpts.resize(num);
  for (int i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readBasic, cnt MSCPTS");
    if (cnt != i) 
      error("GMTK_GMParms::readBasic, MSCPTs, out of order count",cnt);
    msCpts[i] = new MSCPT;
    msCpts[i]->read(is);
  }

}


void 
GMParms::writeBasic(oDataStreamFile& os)
{

  os.write(MAGIC_PRM_FILE,"GMTK_GMParms::writeBasic, magic");
  os.nl();

  os.write(dPmfs.len(),"GMTK_GMParms::writeBasic, dpmfs");
  os.nl();
  for (int i=0;i<dPmfs.len();i++) {
    os.write(i);
    dPmfs[i]->write(os);
  }
  os.nl();

  os.write(sPmfs.len(),"GMTK_GMParms::writeBasic, spmfs");
  os.nl();
  for (int i=0;i<sPmfs.len();i++) {
    os.write(i);
    sPmfs[i]->write(os);
  }
  os.nl();

  os.write(means.len(),"GMTK_GMParms::writeBasic, means");
  os.nl();
  for (int i=0;i<means.len();i++) {
    os.write(i);
    means[i]->write(os);
  }
  os.nl();

  os.write(covars.len(),"GMTK_GMParms::writeBasic, covars");
  os.nl();
  for (int i=0;i<covars.len();i++) {
    os.write(i);
    covars[i]->write(os);
  }
  os.nl();

  os.write(dLinkMats.len(),"GMTK_GMParms::writeBasic, DlinkMatrix");
  os.nl();
  for (int i=0;i<dLinkMats.len();i++) {
    os.write(i); os.nl();
    dLinkMats[i]->write(os);
  }
  os.nl();

  os.write(weightMats.len(),"GMTK_GMParms::writeBasic, WeightMatrix");
  os.nl();
  for (int i=0;i<weightMats.len();i++) {
    os.write(i); os.nl();
    weightMats[i]->write(os);
  }
  os.nl();


  os.write(mdCpts.len(),"GMTK_GMParms::writeBasic, MDCPT");
  os.nl();
  for (int i=0;i<mdCpts.len();i++) {
    os.write(i); os.nl();
    mdCpts[i]->write(os);
  }
  os.nl();

  os.write(msCpts.len(),"GMTK_GMParms::writeBasic, MSCPT");
  os.nl();
  for (int i=0;i<msCpts.len();i++) {
    os.write(i); os.nl();
    msCpts[i]->write(os);
  }
  os.nl();
}



void 
GMParms::readDTs(iDataStreamFile& is)
{
  int num;
  int cnt;
  char *str;

  is.read(str,"GMTK_GMParms::readDTs, magic");
  if (strcmp(str,MAGIC_DT_FILE))
    error("GMTK_GMParms::readDTs. Expecting DT file, got (%s)",str);
  delete [] str;

  is.read(num,"GMTK_GMParms::readDTs, dpmfs");
  if (num < 0) error("GMTK_GMParms::readDTs num dpmfs = %d",num);
  dts.resize(num);
  for (int i=0;i<num;i++) {
    is.read(cnt,"GMTK_GMParms::readDTs, cnt");
    if (cnt != i) 
      error("GMTK_GMParms::readDTs, out of order count",cnt);

    dts[i] = new RngDecisionTree<int>;
    dts[i]->read(is);
  }
}



void 
GMParms::writeDTs(oDataStreamFile& os)
{
  os.write(MAGIC_DT_FILE,"GMTK_GMParms::writeDTs, magic");
  os.nl();
  os.write(dts.len(),"GMTK_GMParms::writeDTs, dpmfs");
  os.nl();
  for (int i=0;i<dts.len();i++) {
    os.write(i); os.nl();
    dts[i]->write(os);
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
