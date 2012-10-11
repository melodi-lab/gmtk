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
VECPT0_FILE % file to read in.
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
END

1
VECPT1
1 % num par
4 % par card
2 % self card
VECPT1_FILE % file to read in
nfs:4 % nfloats
fmt:ascii
END

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
#include <typeinfo>

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_VECPT.h"
#include "GMTK_GMParms.h"
#include "GMTK_DiscRV.h"
#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_Filter.h"
#  include "GMTK_FilterFile.h"
#  include "GMTK_ObservationFile.h"
#  include "GMTK_FileSource.h"
#endif

VCID("$Header$")

#define FILE_NAME_WHEN_USING_GLOBAL_OBSERVATION_FILE "**global**"

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
  obs_file_foffset = 0;
  obs_file_ioffset = 0;
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

  ///////////////////////////////////////////////////////////////////////
  // The first set of options match that of any other CPT. Namely,
  // we have
  //   1) a name, 
  //   2) num parents (which must be 1 in this case), 
  //   3) parent cardinality
  //   4) self cardinaltiy (which must be 2 in this case).
  //   5) file observation name for the VECPT

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

  try {

    if(obsFileName.length()==0) {
      string error_message;
      stringprintf(error_message,"need to supply observation file name");
      throw(error_message);
    }

    if(cardinalities[0]==0) {
      string error_message;
      stringprintf(error_message,"must supply parent cardinality");
      throw(error_message);
    }


    // We have two options: The VECPT values can either come from a
    // separate observation file, or can come from the global
    // observation matrix. In this later case, then the options change
    // (and we can specify an optional float and int offset into this
    // file).

    if (obsFileName == FILE_NAME_WHEN_USING_GLOBAL_OBSERVATION_FILE) {
      // we take observations from the global observation matrix.
      obs = static_cast<FileSource *>(globalObservationMatrix);
      // I think the above should be safe, as there's no chance of
      // trying to read a VECPT from anything other than a FileSource
    
      is.read(str);
      while(! (is.isEOF() || str=="END")) {
	lineNum++;

	// parse the string we have just read
	string::size_type len=str.length();
	string::size_type pos=str.find(":",0);
	if (pos == string::npos || len < 3 || pos == 0) {
	  string error_message;
	  stringprintf(error_message,"Invalid format '%s' which should be of form 'flag:value' where 'flag' is the option name and 'value' the option value",
		       str.c_str());
	  throw(error_message);
	}

	option_name = str.substr(0,pos);
	option_value = str.substr(pos+1,len-pos-1);

	if (option_name == "f_offset") {
	  if (!strIsInt(option_value.c_str(),&obs_file_foffset)) {
	    string error_message;
	    stringprintf(error_message,"Invalid value '%s' for float offset 'f_offset', must be integer.",option_value.c_str());
	    throw(error_message);
	  }
	}
	else if (option_name == "i_offset") {
	  if (!strIsInt(option_value.c_str(),&obs_file_ioffset)) {
	    string error_message;
	    stringprintf(error_message,"Invalid value '%s' for ints offset 'i_offset', must be integer.",option_value.c_str());
	    throw(error_message);
	  }
	}
	else if(option_name == "nfs") {
	  if (!strIsInt(option_value.c_str(),&nfs)) {
	    string error_message;
	    stringprintf(error_message,"Invalid value '%s' for number floats option 'nfs'. Must be integer.",option_value.c_str());
	    throw(error_message);
	  }
	}  
	else if(option_name == "nis") {
	  if (!strIsInt(option_value.c_str(),&nis)) {
	    string error_message;
	    stringprintf(error_message,"Invalid value '%s' for number int option 'nis'. Must be integer.",option_value.c_str());
	    throw(error_message);
	  }
	}  
	else {
	  string error_message;
	  stringprintf(error_message,"Unknown option:value pair '%s'",str.c_str());
	  throw(error_message);
	}

	// is.read() returns (via errorReturn() in
	// miscSupport/fileParser.{h,cc}) false when a problem occurs, in
	// particular, when we reach the EOF.  Alternatively we could
	// check for the end of file directly, feof(fh), which I am doing here
	is.read(str);
      } 

      // only dense case is supported when VECPT scores come from global observation matrix.
      // veMode = VE_Dense;

      // check valid number of ints. 
      if (!(nis == 0 || (nis == nfs))) {
	string error_message;
	stringprintf(error_message,"specifies %d ints and %d floats, but must have either 0 ints, or same number of ints as floats",
		     nis,nfs);
	throw(error_message);
      }
      if (nfs == 0) {
	string error_message;
	stringprintf(error_message,"specifies %d floats, but must have > 0 floats",
		     nfs);
	throw(error_message);
      }

      if (nis == 0) {
	// check that it works with the current global observation matrix.
	if (nfs + obs_file_foffset > globalObservationMatrix->numContinuous()) {
	  string error_message;
	  stringprintf(error_message,"specifies %d floats and offset %d, but global observation matrix only has %d",
		       nfs,obs_file_foffset,globalObservationMatrix->numContinuous());
	  throw(error_message);
	}
	if (cardinalities[0] != nfs) {
	  string error_message;
	  stringprintf(error_message,"specifies %d floats and offset %d, but random variable has cardinality %d",
		       nfs,obs_file_foffset,cardinalities[0]);
	  throw(error_message);
	}
	// must mean dense mode.
	veMode = VE_Dense;
      } else if (nis > 0) {
	if (nis != nfs) {
	  error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' with %d ints and %d floats. Must be equal in sparse case.",
		is.fileName(),is.lineNo(),name().c_str(),nis,nfs);
	}
	// check that it works with the current global observation matrix.
	if (nfs + obs_file_foffset > globalObservationMatrix->numContinuous()) {
	  string error_message;
	  stringprintf(error_message,"specifies %d floats and offset %d, but global observation matrix only has %d real entries per frame",
		       nfs,obs_file_foffset,globalObservationMatrix->numContinuous());
	  throw(error_message);
	}
	if (nis + obs_file_ioffset > globalObservationMatrix->numDiscrete()) {
	  string error_message;
	  stringprintf(error_message,"specifies %d ints and offset %d, but global observation matrix only has %d discrete entries per frame",
		       nis,obs_file_ioffset,globalObservationMatrix->numDiscrete());
	  throw(error_message);
	}

	if (cardinalities[0] < nfs) {
	  string error_message;
	  stringprintf(error_message,"specifies %d floats/ints but random variable has cardinality %d",
		       nfs,cardinalities[0]);
	  throw(error_message);
	}	

	// must mean sparse mode.
	veMode = VE_Sparse;
      } 

    } else {

      ///////////////////////////////////////////////////////////////////////
      // Next, we have a set of optional arguments that the user may give
      // to the observation code. If an optional argument is not given,
      // then the default value may be used.  These optional arguments
      // consist of lines of a "flag : value" syntax, where "flag"
      // indicates the current argument, and "value" is its value.
      // The optional arguments must end with the string "END"

      //  bool done=false;

      is.read(str);
      while(! (is.isEOF() || str=="END")) {
	lineNum++;

	// parse the string we have just read
	string::size_type len=str.length();
	string::size_type pos=str.find(":",0);
	if (pos == string::npos || len < 3 || pos == 0) {
	  string error_message;
	  stringprintf(error_message,"Invalid format '%s' which should be of form 'flag:value' where 'flag' is the option name and 'value' the option value",
		       str.c_str());
	  throw(error_message);
	}

	option_name = str.substr(0,pos);
	option_value = str.substr(pos+1,len-pos-1);


	if(option_name == "nfs") {
	  if (!strIsInt(option_value.c_str(),&nfs)) {
	    string error_message;
	    stringprintf(error_message,"Invalid value '%s' for number floats option 'nfs'. Must be integer.",option_value.c_str());
	    throw(error_message);
	  }
	}
	else if(option_name == "nis") {
	  if (!strIsInt(option_value.c_str(),&nis)) {
	    string error_message;
	    stringprintf(error_message,"Invalid value '%s' for number int option 'nis'. Must be integer.",option_value.c_str());
	    throw(error_message);
	  }
	}  
	else if(option_name == "frs") {
	  // range code needs to check this syntax.
	  frs=option_value;
	}  
	else if(option_name == "irs") {
	  // range code needs to check this syntax.
	  irs=option_value;
	}  
	else if(option_name == "pr") {
	  // range code needs to check this syntax.
	  pr_rs=new char[option_value.length()];
	  strcpy(pr_rs,option_value.c_str());
	}  
	else if(option_name == "fmt") {
	  // observation file code needs to check this syntax.
	  fmt=option_value;
	}  
	else if(option_name == "swap") {
	  char c = (option_value.c_str())[0];
	  if (c == 't' || c == 'T')
	    iswp = true;
	  else if (c == 'f' || c == 'F')
	    iswp = false;
	  else {
	    string error_message;
	    stringprintf(error_message,"Endian swap condition given as '%c', must be T or F",c);
	    throw(error_message);
	  }
	}
	else if(option_name == "preTransforms") {
	  preTransforms=new char[option_value.length()];
	  strcpy(preTransforms,option_value.c_str());
	}  
	else if(option_name == "postTransforms") {
	  postTransforms=new char[option_value.length()];
	  strcpy(postTransforms,option_value.c_str());
	}
	else if(option_name == "sentRange") {
	  sentRange=new char[option_value.length()];
	  strcpy(sentRange,option_value.c_str());
	}  
	else {
	  string error_message;
	  stringprintf(error_message,"Unknown option:value pair '%s'",str.c_str());
	  throw(error_message);
	}

	// is.read() returns (via errorReturn() in
	// miscSupport/fileParser.{h,cc}) false when a problem occurs, in
	// particular, when we reach the EOF.  Alternatively we could
	// check for the end of file directly, feof(fh), which I am doing here
	is.read(str);
      } 

      // check valid number of ints. 
      if (!(nis == 0 || (nis == nfs))) {
	string error_message;
	stringprintf(error_message,"specifies %d ints and %d floats, but must have either 0 ints, or same number of ints as floats",
		     nis,nfs);
	throw(error_message);
      }
      if (nfs == 0) {
	string error_message;
	stringprintf(error_message,"specifies %d ints and %d floats, but must have > 0 floats",
		     nis,nfs);
	throw(error_message);
      }

#if 0
      obs = new ObservationMatrix;
      // Now try opening the file:
      obs->openFile(obsFileName.c_str(),
		    frs.c_str(),
		    irs.c_str(),
		    nfs,
		    nis,
		    obs->formatStrToNumber(fmt.c_str()),
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
#else
      assert(fmt=="ascii");
      int ifmt = formatStrToNumber(fmt.c_str());
      if (ifmt < 0) {
	string error_message;
	stringprintf(error_message,"unknown file format '%s' for VE CPT", fmt.c_str());
	throw(error_message);
      }
      ObservationFile *obsFile = 
	instantiateFile((unsigned)ifmt, (char *)obsFileName.c_str(), nfs, nis,
			0, iswp, false, NULL, NULL, NULL, NULL,	sentRange);

      Filter *preFilt = instantiateFilters(preTransforms, obsFile->numLogicalContinuous(),
                                                          obsFile->numLogicalDiscrete());
      obsFile = new FilterFile(preFilt, obsFile, frs.c_str(), irs.c_str(), pr_rs);
      assert(obsFile);
      Filter *postTrans = instantiateFilters(postTransforms, obsFile->numLogicalContinuous(),
                                                             obsFile->numLogicalDiscrete());
      obsFile = new FilterFile(postTrans, obsFile);
      obs = new FileSource(obsFile, DEFAULT_BUFFER_SIZE,
			   globalObservationMatrix->startSkip(), 
			   globalObservationMatrix->endSkip());
#endif

      // still here? Do more error checking.

      if (obs->numContinuous() == 0) {
	error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' and reading observation file '%s' with %d ints and %d floats. Range string '%s'. Must have or specify > 0 floats",
	      is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs,frs.c_str());
      }

      if (!((obs->numDiscrete() == 0) || (obs->numDiscrete() == obs->numContinuous()))) {
	error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' and reading observation file '%s' with %d ints and %d floats. Float range string is '%s', int range string is '%s'. Must have or specify either 0 ints, or same number of ints as floats",
	      is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs,frs.c_str(),irs.c_str());
      }

      FileSource fs;
      if ( typeid(*globalObservationMatrix) == typeid(fs) &&
	   typeid(*obs) == typeid(fs) )
      {
	FileSource *gomFS = static_cast<FileSource *>(globalObservationMatrix);
	FileSource *obsFS = static_cast<FileSource *>(obs);
	// make sure we have same number of segments as global observation case.
	if (gomFS->numSegments() != obsFS->numSegments()) {
	  error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' and reading observation file '%s' with %d ints and %d floats. Number of segments %d must match that of the global observation file, which has %d segments",
		is.fileName(),is.lineNo(),name().c_str(),obsFileName.c_str(),nis,nfs,obsFS->numSegments(),gomFS->numSegments());
	}
      } else {
	warning("WARNING: StreamSource cannot verify that VECPT and observations have the same number of segments");
      }

      // we check that the cardinality is compatible with the VECPT.
      // we need either that
      //    1. obs->numDiscrete() == 0 and parent_random_var_card == obs->numContinuous()
      // or
      //    2. (obs->numDiscrete() == obs->numContinuous()), and (obs->numDiscrete() <= parent_random_var_card), and
      //       all discrete values are in the range 0 <= val <= (parent_random_var_card-1)
      // 

      if (obs->numDiscrete() == 0 && cardinalities[0] == obs->numContinuous()) {
	veMode = VE_Dense;
      } else if ( obs->numDiscrete() == obs->numContinuous() && 
		  obs->numDiscrete() <= cardinalities[0] )  {
	veMode = VE_Sparse;
      } else {
	error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s' with parent cardinality %d and reading observation file '%s' with %d ints and %d floats. Either num floats in obs file must match parent cardinality, or alternatively obs file must have same number floats as ints, and both must be no more than parent cardinality.",
	      is.fileName(),is.lineNo(),name().c_str(),cardinalities[0],
	      obsFileName.c_str(),nis,nfs);
      }

    }
  } catch ( string error_message ) {
    error("ERROR: reading file '%s' line %u, of VirtualEvidenceCPT spec '%s': %s",
	  is.fileName(),is.lineNo(),name().c_str(),error_message.c_str());
  } catch( const char * const error_message ) {
    error("ERROR: reading file '%s' line %u, of VirtualEvidenceCPT spec '%s': %s",
	  is.fileName(),is.lineNo(),name().c_str(),error_message);
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
  if (obsFileName == FILE_NAME_WHEN_USING_GLOBAL_OBSERVATION_FILE) {
    sprintf(tmp,"f_offset:%u",obs_file_foffset);
    os.write(tmp);
    sprintf(tmp,"i_offset:%u",obs_file_ioffset);
    os.write(tmp);
    sprintf(tmp,"nfs:%u",nfs);
    os.write(tmp);
  } else {
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
  }
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


logpr VECPT::probGivenParents(vector <RV *>& parents,
			      DiscRV * drv) 
{
  assert ( bitmask & bm_basicAllocated );
  curParentValue = RV2DRV(parents[0])->val;
  register DiscRVType val = drv->val;

  logpr p((void*)NULL);
  if (veMode == VE_Dense) {
    p.valref() = (*obs->floatVecAtFrame(drv->frame(),curParentValue+obs_file_foffset));
    // The obseved value (1) is the one corresponding
    // to the value in the file. I.e., the score 
    // in the file corresponds to Pr(child = 1 | parent = j) = f_t(j), 
    // and f_t(j) is the value stored in the file.  
    if (val == 0) {
      // if the child RV has zero value, then we invert the probability.
      // see comment 'zero valued child case' elsewhere in this file.
      p = 1.0 - p;
    }
  } else {
    // do a slow linear search since order can be anything.
    // TODO: make sorted assumption and do binary search.
    unsigned *base = obs->unsignedVecAtFrame(drv->frame());
    unsigned i;
    for (i=obs_file_ioffset;i<obs->numDiscrete();i++) {
      if (base[i] == curParentValue)
	break;
    }
    if (i<obs->numDiscrete()) {
      p.valref() = (*obs->floatVecAtFrame(drv->frame(),i));
      // The obseved value (1) is the one corresponding
      // to the value in the file.
      if (val == 0) {
	p = 1.0 - p;
      }
    } else {
      // parent value is not in the table of parent values in the
      // current frame's int observation vector.
      if (val == 0) {
	// see comment 'zero valued child case' elsewhere in this file.
	p.set_to_one();
      } else {
        // default: If the child is == 1, and there are any parent
	// values that occur during inference that are not given in the observation
	// file (which can happen in this case since the parent values are sparses),
	// then for this case we give any of those parent values a zero score.
	// I.e., this corresponds to Pr( child = 1 | parents = j) = f_t(j)
	//  unspecified entries in the file for f_t(j) are presumed to be zero.
	p.set_to_zero();
      }
    }
  }
  return p;
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

void VECPT::begin(iterator& it,DiscRV* drv, logpr& p) 
{
  assert ( bitmask & bm_basicAllocated );
  it.setCPT(this);
  it.drv = drv;
  it.uInternalState = curParentValue;
  drv->val = 0;

  if (veMode == VE_Dense) {
    p.valref() = (*obs->floatVecAtFrame(drv->frame(),curParentValue+obs_file_foffset));
    // zero valued child case: since this is the zero value of the
    // child, and the child is presumed to be binary, and since the
    // values stored in the observaion file correspond to Pr(child =
    // 1| parents), we need to invert the value in this case (assuming
    // that it is a probability). Of course, if the values in the
    // observaion file are not probabilities (which is valid if we
    // want to use the VECPT as an arbitrary score), then this next
    // operation does not make sense, but then in that case the child
    // should not be hidden in which case this code is not run.
    p = 1.0 - p;
  } else {
    // Sparse case.
    // 
    // do a slow linear search since order can be anything.
    unsigned *base = obs->unsignedVecAtFrame(drv->frame());
    unsigned i;
    for (i=obs_file_ioffset;i<obs->numDiscrete();i++) {
      if (base[i] == curParentValue)
	break;
    }
    // since this is a begin, this case corresponds to where
    // the child is hidden and equal to value 0, so we use the zero-value case.
    if (i<obs->numDiscrete()) {
      // then we found an entry in the observation file corresponding to the current parent value,
      // and we use the score in the observaion file at this position.
      p.valref() = (*obs->floatVecAtFrame(drv->frame(),i));
      // zero valued child case: see comment above.
      p = 1.0 - p;
    } else {
      // we did not find a parent value for this case, so this
      // means that 
      // Pr(child = 0 | parent) = 1 for the zero valued child case.
      p.set_to_one();
    }
  }
}

bool VECPT::next(iterator &it,logpr& p)
{
  if (it.drv->val == 1)
    return false;
  // here, the RV child now gets value 1.
  it.drv->val = 1;
  if (veMode == VE_Dense) {
    p.valref() = (*obs->floatVecAtFrame(it.drv->frame(),it.uInternalState + obs_file_foffset));
    // we don't invert the score since this is the Pr(child = 1 | parent) case.
  } else {
    // do a slow linear search since order can be anything.
    unsigned *base = obs->unsignedVecAtFrame(it.drv->frame());
    unsigned i;
    for (i=obs_file_ioffset;i<obs->numDiscrete();i++) {
      if (base[i] == it.uInternalState)
	break;
    }
    // since this is next, the child is hidden but
    // here it is the case that the child is at iteration where it
    // it is instantiated to 1 (unity). So we take the scores for
    // the =1 case.
    if (i<obs->numDiscrete()) {
      p.valref() = (*obs->floatVecAtFrame(it.drv->frame(),i));
    } else {
      // we did not find a parent value for this case, so this
      // means that 
      // Pr(child = 0 | parent) = 0 for the unity valued child case.
      p.set_to_zero();
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
  // TODO: this routine is currently not correct (but it is not used).

  logpr sum;
  for (unsigned i=0;i<obs->numContinuous();i++) {
    logpr tmp((void*)NULL);
    tmp.valref() = (*obs->floatVecAtFrame(drv->frame(),i));
    sum += tmp;
  }
  logpr uniform = (rnd.drand48()*sum.unlog());
  sum.set_to_zero();
  unsigned i=0;
  do {
    logpr tmp((void*)NULL);
    tmp.valref() = (*obs->floatVecAtFrame(drv->frame(),i));
    sum += tmp;
    if (uniform <= sum)
      break;
  } while (++i < obs->numContinuous());

  if (veMode == VE_Dense) {
    return i;
  } else {
    // this assumes that ints in sparse case are ordered
    return *(obs->unsignedVecAtFrame(drv->frame()) + i); 
  }

}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#endif
