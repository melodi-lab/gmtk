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
 * An example of how you'd specify VECPTs:
 (in a masterfile:)

VE_CPT_IN_FILE inline

2  % num VECPTs 

0
VECPT0  % name of VECPT
1 % num par
2 % par card
2 % self card
VECPT0_FILE
nfs:2 % nfloats
nis:0 % nints
frs:all % float range
irs:all % int range
pr:all % must be all
fmt:ascii
swap:F % endian swapping condition
preTransforms:X
postTransforms:X
sentRange:all
EOF

1
VECPT1
1 % num par
4 % par card
2 % self card
VECPT1_FILE
nfs:4 % nfloats
fmt:ascii
EOF

Note the "EOF" at the end of each VECPT.  Also, parsing of the VECPTs is
rudimentary and no spaces are allowed around the ":".

Also the order in which the arguments are written is not important.


The default values are:

nfs:0  % will cause an error if not updated
nis:0
frs:all
irs:all
pr:all
fmt:pfile
swap:F 
preTransforms:X
postTransforms:X
sentRange:all


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

  //Initialize defaults
  // _numParents and _card already initialized in VECPT constructor
  cardinalities[0]=0; // triggers an error if not updated
  nfs=0;  // idem
  nis=0;
  frs="all";
  irs="all";
  pr_rs=NULL;
  fmt="pfile";
  iswp=false;
  sentRange=NULL;
  preTransforms=NULL;
  postTransforms=NULL;

  string str;
  string option_name;
  string option_value;
  unsigned lineNum=0;
  // read the name of the object.
  NamedObject::read(is);

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
  
  //  bool done=false;
  is.read(str);
  while(! (is.isEOF() || str=="EOF")) {
    lineNum++;
    // parse the string we have just read
    string::size_type len=str.length();
    string::size_type pos=str.find(":",0);
    if (pos == string::npos) error("ERROR: Invalid option format at line %u while reading the VECPT.  The format should be x:y where x is the option name and y the option value",is.lineNo());
    option_value=str.substr(pos+1,len-pos-1);
    option_name=str.substr(0,pos);
    char * tmpString=new char[len];
    strcpy(tmpString,option_value.c_str());
    if(option_name == "nfs") {
      sscanf(option_value.c_str(),"%u",&nfs);
    }  
    else if(option_name == "nis") {
      sscanf(option_value.c_str(),"%u",&nis);
    }  
    else if(option_name == "frs") {
      frs=option_value;
    }  
    else if(option_name == "irs") {
      irs=option_value;
    }  
    else if(option_name == "pr") {
      pr_rs=new char[option_value.length()];
      strcpy(pr_rs,option_value.c_str());
    }  
    else if(option_name == "fmt") {
      fmt=option_value;
    }  
    else if(option_name == "swap") {
      char c = (option_value.c_str())[0];
      if (c == 't' || c == 'T')
	iswp = true;
      else if (c == 'f' || c == 'F')
	iswp = false;
      else {
	error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' just before reading observation file '%s' with %d ints and %d floats. Endian swap condition is '%c', must be T or F",
	      is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs,c);
      }
    }  
    else if(option_name == "preTransforms") {
      sscanf(option_value.c_str(),"%s",preTransforms);
    }  
    else if(option_name == "postTransforms") {
      sscanf(option_value.c_str(),"%s",postTransforms);
    }  
    else if(option_name == "sentRange") {
      sentRange=new char[option_value.length()];
      strcpy(sentRange,option_value.c_str());
    }  
    else {
      error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s'.  Unknown option line '%s'",
	    is.fileName(),is.lineNo(),name().c_str(),str.c_str());
    }

    delete [] tmpString;

    // is.read() returns (via errorReturn() in
    // miscSupport/fileParser.{h,cc}) false when a problem occurs, in
    // particular, when we reach the EOF.  Alternatively we could
    // check for the end of file directly, feof(fh), which I am doing here
    is.read(str);
  } 

  // check valid number of ints. 
  if (!(nis == 0 || (nis == nfs))) {
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' just before reading observation file '%s' with %d ints and %d floats. Must have either 0 ints, or same number of ints as floats",
	  is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs);
  }
  if (nfs == 0) {
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' just before reading observation file '%s' with %d ints and %d floats. Must have > 0 floats",
	  is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs);
  }



  //  if (pr_rs.length() != 0) pr_str=(char*)pr_rs.c_str();

  //  if(sentRange.length() !=0) sentRangeStr=(char*)sentRange.c_str();

  if(obsFileName.length()==0) {
    error("ERROR: need to supply the observation file name");
  }

  if(cardinalities[0]==0) {
    error("ERROR: need to supply parent cardinality");
  }
  
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
	       //	       (const char**)&pr_str,
	       (const char**)&pr_rs,
	       NULL,
	       NULL,
	       &preTransforms,
	       postTransforms,
	       0,  // FTROP_NONE
	       //	       (const char**)&sentRangeStr
	       (const char**)&sentRange
	       );


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

#define VECPT_TMP_OUTPUT_STR_LEN 500

  char* tmp=new char[VECPT_TMP_OUTPUT_STR_LEN];

  // write without any formatting.
  NamedObject::write(os);

  sprintf(tmp,"numParents:%d",1);
  os.write("numParents:1");  // num parents
  sprintf(tmp,"parentCard:%u",cardinalities[0]);
  os.write(tmp);
  os.write("selfCard:2"); // self card
  sprintf(tmp,"of:%s",obsFileName.c_str());
  os.write(tmp);
  sprintf(tmp,"nfs:%u",nfs);
  os.write(tmp);
  sprintf(tmp,"nis:%u",nis);
  os.write(tmp);
  sprintf(tmp,"frs:%s",frs.c_str());
  os.write(tmp);
  sprintf(tmp,"nfs:%s",irs.c_str());
  os.write(tmp);
  sprintf(tmp,"pr_rs:%s",pr_rs);
  os.write(tmp);
  sprintf(tmp,"fmt:%s",fmt.c_str());
  os.write(fmt);
  sprintf(tmp,"swap:%s",iswp?"T":"F");
  os.write((iswp?"T":"F"));
  sprintf(tmp,"preTransforms:%s",preTransforms);
  os.write(tmp);
  sprintf(tmp,"postTransforms:%s",postTransforms);
  os.write(tmp);
  sprintf(tmp,"sentRange:%s",sentRange);
  os.write(tmp);
  os.write("EOF");
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
