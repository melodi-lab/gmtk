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


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////

GMParms::GMParms()
{}



void 
GMParms::readBasic(iDataStreamFile& is)
{
  int tmp;

  is.read(tmp,"GMTK_GMParms::readBasic, dpmfs");
  if (tmp < 0) error("GMTK_GMParms::readBasic num dpmfs = %d",tmp);
  dPmfs.resize(tmp);
  for (int i=0;i<tmp;i++) {
    dPmfs[i] = new Dense1DPMF;
    dPmfs[i]->read(is);
  }


  is.read(tmp,"GMTK_GMParms::readBasic, spmfs");
  if (tmp < 0) error("GMTK_GMParms::readBasic num spmfs = %d",tmp);
  sPmfs.resize(tmp);
  for (int i=0;i<tmp;i++) {
    sPmfs[i] = new Sparse1DPMF;
    sPmfs[i]->read(is);
  }

  is.read(tmp,"GMTK_GMParms::readBasic, means");
  if (tmp < 0) error("GMTK_GMParms::readBasic num means = %d",tmp);
  means.resize(tmp);
  for (int i=0;i<tmp;i++) {
    means[i] = new MeanVector;
    means[i]->read(is);
  }

  is.read(tmp,"GMTK_GMParms::readBasic, covars");
  if (tmp < 0) error("GMTK_GMParms::readBasic num covars = %d",tmp);
  covars.resize(tmp);
  for (int i=0;i<tmp;i++) {
    covars[i] = new DiagCovarVector;
    covars[i]->read(is);
  }

  is.read(tmp,"GMTK_GMParms::readBasic, DlinkMatrix");
  if (tmp < 0) error("GMTK_GMParms::readBasic num DlinkMatrix = %d",tmp);
  dLinkMats.resize(tmp);
  for (int i=0;i<tmp;i++) {
    dLinkMats[i] = new DlinkMatrix;
    dLinkMats[i]->read(is);
  }

  is.read(tmp,"GMTK_GMParms::readBasic, WeightMatrix");
  if (tmp < 0) error("GMTK_GMParms::readBasic num WeightMatrix = %d",tmp);
  weightMats.resize(tmp);
  for (int i=0;i<tmp;i++) {
    weightMats[i] = new WeightMatrix;
    weightMats[i]->read(is);
  }


  is.read(tmp,"GMTK_GMParms::readBasic, MDCPT");
  if (tmp < 0) error("GMTK_GMParms::readBasic num MDCPT = %d",tmp);
  mdCpts.resize(tmp);
  for (int i=0;i<tmp;i++) {
    mdCpts[i] = new MDCPT;
    mdCpts[i]->read(is);
  }

  is.read(tmp,"GMTK_GMParms::readBasic, MSCPT");
  if (tmp < 0) error("GMTK_GMParms::readBasic num MSCPT = %d",tmp);
  msCpts.resize(tmp);
  for (int i=0;i<tmp;i++) {
    msCpts[i] = new MSCPT;
    msCpts[i]->read(is);
  }

}


void 
GMParms::writeBasic(oDataStreamFile& os)
{

  os.write(dPmfs.len(),"GMTK_GMParms::writeBasic, dpmfs");
  os.nl();
  for (int i=0;i<dPmfs.len();i++) {
    dPmfs[i]->write(os);
  }
  os.nl();

  os.write(sPmfs.len(),"GMTK_GMParms::writeBasic, spmfs");
  os.nl();
  for (int i=0;i<sPmfs.len();i++) {
    sPmfs[i]->write(os);
  }
  os.nl();

  os.write(means.len(),"GMTK_GMParms::writeBasic, means");
  os.nl();
  for (int i=0;i<means.len();i++) {
    means[i]->write(os);
  }
  os.nl();

  os.write(covars.len(),"GMTK_GMParms::writeBasic, covars");
  os.nl();
  for (int i=0;i<covars.len();i++) {
    covars[i]->write(os);
  }
  os.nl();

  os.write(dLinkMats.len(),"GMTK_GMParms::writeBasic, DlinkMatrix");
  os.nl();
  for (int i=0;i<dLinkMats.len();i++) {
    dLinkMats[i]->write(os);
  }
  os.nl();

  os.write(weightMats.len(),"GMTK_GMParms::writeBasic, WeightMatrix");
  os.nl();
  for (int i=0;i<weightMats.len();i++) {
    weightMats[i]->write(os);
  }
  os.nl();


  os.write(mdCpts.len(),"GMTK_GMParms::writeBasic, MDCPT");
  os.nl();
  for (int i=0;i<mdCpts.len();i++) {
    mdCpts[i]->write(os);
  }
  os.nl();

  os.write(msCpts.len(),"GMTK_GMParms::writeBasic, MSCPT");
  os.nl();
  for (int i=0;i<msCpts.len();i++) {
    msCpts[i]->write(os);
  }
  os.nl();
}



void 
GMParms::readDTs(iDataStreamFile& is)
{
  int tmp;
  is.read(tmp,"GMTK_GMParms::readDTs, dpmfs");
  if (tmp < 0) error("GMTK_GMParms::readDTs num dpmfs = %d",tmp);
  dts.resize(tmp);
  for (int i=0;i<tmp;i++) {
    dts[i] = new RngDecisionTree<int>;
    dts[i]->read(is);
  }
}



void 
GMParms::writeDTs(oDataStreamFile& os)
{
  os.write(dts.len(),"GMTK_GMParms::writeDTs, dpmfs");
  os.nl();
  for (int i=0;i<dts.len();i++) {
    dts[i]->write(os);
  }
}




#ifdef MAIN

GMParms GM_Parms;

int
main()
{
  iDataStreamFile is("dataFiles/test1.gmb",false);
  GM_Parms.readBasic(is);

  // oDataStreamFile os("dataFiles/test1_out.gmb");
  oDataStreamFile os("-");
  GM_Parms.writeBasic(os);
}


#endif
