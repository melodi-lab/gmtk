/*-
 *
 * GMTK_VECPT.cc
 *     Virtual Evidence CPT.
 *     A CPT for P(C|A=a) when C is binary, and where the
 *     scores come from an observed file (pfile, htk file, ascii, binary).
 *
 *     This CPT allows GMTK to be used either as a 
 *         1) Hybrid DBN/ANN (Dynamic Bayesian Network/Artificial Neural Network), or a
 *         2) Hybrid DBN/SVM, or a
 *         3) Hybrid DBN/KM (Kernal machine) decoder.
 *     or anything else, depending on the types of scores given in the observatioon file.
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
#include "rand.h"

#include "GMTK_VECPT.h"
#include "GMTK_GMParms.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ObservationMatrix.h"

VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * VECPT::read(is)
 *      read in the table.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the 'mscpt' member function in the object.
 *
 *-----------------------------------------------------------------------
 */

void
VECPT::read(iDataStreamFile& is)
{
  string str;

  // read the name of the object.
  NamedObject::read(is);

  // read the cardinalities.
  is.read(_numParents,"Can't read VirtualEvidenceCPT's num parents");
  if (_numParents != 1) 
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' must have exactly one(1) rather than %d parents",
	  is.fileName(),is.lineNo(),name().c_str(),_numParents);
  cardinalities.resize(_numParents);
  // read the parent cardinality
  is.read(cardinalities[0],"Can't read VirtualEvidenceCPT's parent cardinality");

  // read the self cardinality, must be binary
  is.read(_card,"Can't read VirtualEvidenceCPT's self cardinality");
  if (_card != 2)
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s', child cardinality must be two(2) rather than %d.",
	  is.fileName(),is.lineNo(),name().c_str(),_card);




  // next read the file name from which this is going to get its probabilties.
  is.read(obsFileName,"Can't read VirtualEvidenceCPT's obs file name");
  // read number of floats in file
  is.read(nfs,"Can't read VirtualEvidenceCPT's obs file num floats");
  // read number of ints in file
  is.read(nis,"Can't read VirtualEvidenceCPT's obs file num ints");
  // check valid number of ints. 
  if (!(nis == 0 || (nis == nfs))) {
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' just before reading observation file '%s' with %d ints and %d floats. Must have either 0 ints, or same number of ints as floats",
	  is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs);
  }
  if (nfs == 0) {
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' just before reading observation file '%s' with %d ints and %d floats. Must have either > 0 floats",
	  is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs);
  }
  // read float range in file
  is.read(frs,"Can't read VirtualEvidenceCPT's obs file float range");
  // read int range in file
  is.read(irs,"Can't read VirtualEvidenceCPT's obs file int range");
  // per observation range in file
  is.read(pr_rs,"Can't read VirtualEvidenceCPT's obs file per-utterance range string");
  // read format string
  is.read(fmt,"Can't read VirtualEvidenceCPT's obs file file format");
  // endian swapping condition
  is.read(str,"Can't read VirtualEvidenceCPT's obs file swapping condition");
  // use first letter to determine condition.
  char c = (str.c_str())[0];
  if (c == 't' || c == 'T')
    iswp = true;
  else if (c == 'f' || c == 'F')
    iswp = false;
  else {
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' just before reading observation file '%s' with %d ints and %d floats. Endian swap condition is '%c', must be T or F",
	  is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs,c);
  }

  // Read Observation Matrix tranform string
  is.read(preTransforms,"Can't read VirtualEvidenceCPT's obs file pre-transform string");
  is.read(postTransforms,"Can't read VirtualEvidenceCPT's obs file post-transform string");

  const char*pr_str=pr_rs.c_str();

  // Now try opening the file:
  obs.openFile(obsFileName.c_str(),
	       frs.c_str(),
	       irs.c_str(),
	       nfs,
	       nis,
	       obs.formatStrToNumber(fmt.c_str()),
	       iswp,
	       globalObservationMatrix.startSkip(),
	       globalObservationMatrix.endSkip(),
	       false,  // ); // do not run CPP if ascii file.
	       NULL,
	       (const char**)&pr_str,
	       NULL,
	       NULL,
	       &preTransforms,
	       postTransforms
	       );


  // ultimately add this:  Added - Karim
  //      NULL, // CPP command options
  // pr_rs.c_str());
  
  // still here? Do more error checking.

  if (obs.numContinuous() == 0) {
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' and reading observation file '%s' with %d ints and %d floats. Range string '%s'. Must have or specify > 0 floats",
	  is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs,frs.c_str());
  }

  if (!((obs.numDiscrete() == 0) || (obs.numDiscrete() == obs.numContinuous()))) {
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' and reading observation file '%s' with %d ints and %d floats. Float range string is '%s', int range string is '%s'. Must have or specify either 0 ints, or same number of ints as floats",
	  is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs,frs.c_str(),irs.c_str());
  }

  // make sure we have same number of segments as global observation case.
  if (globalObservationMatrix.numSegments() != obs.numSegments()) {
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' and reading observation file '%s' with %d ints and %d floats. Number of segments %d must match that of the global observation file, which has %d segments",
	  is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs,obs.numSegments(),globalObservationMatrix.numSegments());
  }

  // we check that the cardinality is compatible with the VECPT.
  // we need either that
  //    1. obs.numDiscrete() == 0 and parent_random_var_card == obs.numContinuous()
  // or
  //    2. (obs.numDiscrete() == obs.numContinuous()), and (obs.numDiscrete() <= parent_random_var_card), and
  //       all discrete values are in the range 0 <= val <= (parent_random_var_card-1)
  // 

  if (obs.numDiscrete() == 0 && cardinalities[0] == obs.numContinuous()) {
    veMode = VE_Dense;
  } else if ( obs.numDiscrete() == obs.numContinuous() && 
	       obs.numDiscrete() <= cardinalities[0] )  {
    veMode = VE_Sparse;
  } else {
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' with parent cardinality %d and reading observation file '%s' with %d ints and %d floats. Either num floats in obs file must match parent cardinality, or alternatively obs file must have same number floats as ints, and both must be no more than parent cardinality.",
	  is.fileName(),is.lineNo(),name().c_str(),cardinalities[0],
	  obsFileName.c_str(),nis,nfs);
  }

  setBasicAllocatedBit();
}


/*-
 *-----------------------------------------------------------------------
 * VECPT::write(os)
 *      write out data to file 'os'. 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effects other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
VECPT::write(oDataStreamFile& os)
{
  // write without any formatting.
  NamedObject::write(os);

  os.write(1);  // num parents
  os.write(cardinalities[0]);
  os.write(2); // self card

  os.write(obsFileName);
  os.write(nfs);
  os.write(nis);
  os.write(frs);
  os.write(irs);
  os.write(pr_rs);
  os.write(fmt);
  os.write((iswp?"T":"F"));
  os.write(preTransforms);
  os.write(postTransforms);
  os.nl();
  os.nl();
}



////////////////////////////////////////////////////////////////////
//        Probability routines (see GMTK_CPT.h for docs)
////////////////////////////////////////////////////////////////////


void VECPT::becomeAwareOfParentValues( vector <RV *>& parents,
				       const RV* rv ) 
{
  assert ( parents.size() == 1 );
  curParentValue = RV2DRV(parents[0])->val;
}

void VECPT::begin(iterator& it,DiscRV* drv, logpr& p) 
{
  assert ( bitmask & bm_basicAllocated );
  it.setCPT(this);
  it.drv = drv;
  it.uInternalState = curParentValue;
  drv->val = 0;

  if (veMode == VE_Dense) {
    p.valref() = (*obs.floatVecAtFrame(drv->frame(),curParentValue));
    p = 1.0 - p;
  } else {
    // do a slow linear search since order can be anything.
    unsigned *base = obs.unsignedAtFrame(drv->frame());
    unsigned i;
    for (i=0;i<obs.numDiscrete();i++) {
      if (base[i] == curParentValue)
	break;
    }
    if (i<obs.numDiscrete()) {
      p.valref() = (*obs.floatVecAtFrame(drv->frame(),i));
      p = 1.0 - p;
    } else {
      p.set_to_zero();
    }
  }
}

void VECPT::becomeAwareOfParentValuesAndIterBegin(vector < RV *>& parents, 
						  iterator &it,
						  DiscRV* drv,
						  logpr& p) 
{
  becomeAwareOfParentValues(parents,drv);
  it.uInternalState = curParentValue;
  begin(it,drv,p);
}

logpr VECPT::probGivenParents(vector <RV *>& parents,
			      DiscRV * drv) 
{
  assert ( bitmask & bm_basicAllocated );
  curParentValue = RV2DRV(parents[0])->val;
  register DiscRVType val = drv->val;

  logpr p((void*)NULL);
  if (veMode == VE_Dense) {
    p.valref() = (*obs.floatVecAtFrame(drv->frame(),curParentValue));
    // The obseved value (1) is the one corresponding
    // to the value in the file.
    if (val == 0) {
      p = 1.0 - p;
    }
  } else {
    // do a slow linear search since order can be anything.
    // TODO: make sorted assumption and do binary search.
    unsigned *base = obs.unsignedAtFrame(drv->frame());
    unsigned i;
    for (i=0;i<obs.numDiscrete();i++) {
      if (base[i] == curParentValue)
	break;
    }
    if (i<obs.numDiscrete()) {
      p.valref() = (*obs.floatVecAtFrame(drv->frame(),i));
      // The obseved value (1) is the one corresponding
      // to the value in the file.
      if (val == 0) {
	p = 1.0 - p;
      }
    } else {
      if (val == 0) {
	p.set_to_zero();
      } else {
	p.set_to_one();
      }
    }
  }
  return p;
}

bool VECPT::next(iterator &it,logpr& p)
{
  if (it.drv->val == 1)
    return false;
  it.drv->val = 1;
  if (veMode == VE_Dense) {
    p.valref() = (*obs.floatVecAtFrame(it.drv->frame(),it.uInternalState));
  } else {
    // do a slow linear search since order can be anything.
    unsigned *base = obs.unsignedAtFrame(it.drv->frame());
    unsigned i;
    for (i=0;i<obs.numDiscrete();i++) {
      if (base[i] == it.uInternalState)
	break;
    }
    if (i<obs.numDiscrete()) {
      p.valref() = (*obs.floatVecAtFrame(it.drv->frame(),i));
    } else {
      p.set_to_one();
    }
  }
  return true;
}



/*-
 *-----------------------------------------------------------------------
 * VECPT::random sample(is)
 *      random sample according to current 'distribution' 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes val
 *
 *-----------------------------------------------------------------------
 */
int VECPT::randomSample(DiscRV* drv)
{
  // sum up the observations in current frame.

  logpr sum;
  for (unsigned i=0;i<obs.numContinuous();i++) {
    logpr tmp((void*)NULL);
    tmp.valref() = (*obs.floatVecAtFrame(drv->frame(),i));
    sum += tmp;
  }
  logpr uniform = (rnd.drand48()*sum.unlog());
  sum.set_to_zero();
  unsigned i=0;
  do {
    logpr tmp((void*)NULL);
    tmp.valref() = (*obs.floatVecAtFrame(drv->frame(),i));
    sum += tmp;
    if (uniform <= sum)
      break;
  } while (++i < obs.numContinuous());

  if (veMode == VE_Dense) {
    return i;
  } else {
    // this assumes that ints in sparse case are ordered
    return *(obs.unsignedAtFrame(drv->frame()) + i); 
  }

}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#endif
