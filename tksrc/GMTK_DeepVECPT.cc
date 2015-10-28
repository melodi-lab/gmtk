/*-
 *
 * GMTK_DeepVECPT.cc
 *     Deep Virtual Evidence CPT.
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes 
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 * A Virtual Evidence CPT that uses a deep neural network to
 * compute P(parent | observations)
 *
 * An example of how you'd specify DeepVECPTs:
 * (in a masterfile:)



DEEP_VE_CPT_IN_FILE inline 1
0
cptName
1                    % number of parents (always 1)
<parent card>        % must match final layer's # outputs
2                    % self cardinality (always 2)
deepNNName           % DeepNN to use to compute probabilities
prior_1dpmf:priorPMFName  % optional prior for NN output
f_offset:0           % starting index in observation matrix for floats
nfs:0                % # of floats to take from observation matrix
radius:<d>           % # of frames to take from observation matrix = 2d + 1
END


Note: parsing of the DeepVECPTs is rudimentary and no spaces are allowed around the ":".

The default values are:

nfs:0      % will cause an error if not updated
f_offset:0 %
radius:0   % use a single frame as input to the model

 *
 *
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif


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
#include "phipac.h"

#include "GMTK_DeepVECPT.h"
#include "GMTK_GMParms.h"
#include "GMTK_DiscRV.h"
#include "GMTK_FileSource.h"

VCID(HGID)

unsigned DeepVECPT::globalMinSkip = 0;

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * DeepVECPT::read(is)
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
DeepVECPT::read(iDataStreamFile& is)
{

  //Initialize defaults
  // _numParents and _card already initialized in DeepVECPT constructor
  cardinalities[0]=0; // triggers an error if not updated
  nfs=0;  // idem
  obs_file_foffset = 0;
  window_radius = 0;

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
  //   5) DeepNN name to use to get probabilities

  // read the name of the object.
  NamedObject::read(is);
  is.read(_numParents,"Can't read DeepVirtualEvidenceCPT's num parents");
  if (_numParents != 1)
    error("ERROR: reading file '%s' line %d, DeepVirtualEvidenceCPT '%s' must have exactly one (1) rather than %d parents",
          is.fileName(),is.lineNo(),name().c_str(),_numParents);
  cardinalities.resize(_numParents);
  // read the parent cardinality
  is.read(cardinalities[0],"Can't read DeepVirtualEvidenceCPT's parent cardinality");
  cached_CPT = new double[cardinalities[0]];

  // read the self cardinality, must be binary
  is.read(_card,"Can't read DeepVirtualEvidenceCPT's self cardinality");
  if (_card != 2)
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s', child cardinality must be two (2) rather than %d.",
          is.fileName(),is.lineNo(),name().c_str(),_card);

  is.read(str, "Can't read DeepVirtualEvidenceCPT's DeepNN name");
  if (GM_Parms.deepNNsMap.find(str) == GM_Parms.deepNNsMap.end()) {
    error("ERROR: reading file '%s' line %d,  DeepVirtualEvidenceCPT '%s' uses undefined DeepNN '%s'\n", is.fileName(), is.lineNo(), name().c_str(), str.c_str());
  }
  dmlp = GM_Parms.deepNNs[ GM_Parms.deepNNsMap[str] ];
  assert(dmlp);
  input_vector = new float[dmlp->numInputs() + 1];

  try {

    if(cardinalities[0]==0) {
      string error_message;
      stringprintf(error_message,"must supply parent cardinality > 0");
      throw(error_message);
    }
    obs = globalObservationMatrix;
    
    ///////////////////////////////////////////////////////////////////////
    // Next, we have a set of optional arguments that the user may give
    // to the observation code. If an optional argument is not given,
    // then the default value may be used.  These optional arguments
    // consist of lines of a "flag : value" syntax, where "flag"
    // indicates the current argument, and "value" is its value.
    // The optional arguments must end with the string "END"
    
    //  bool done=false;
    
    prior = NULL;
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
      else if(option_name == "nfs") {
	if (!strIsInt(option_value.c_str(),&nfs)) {
	  string error_message;
	  stringprintf(error_message,"Invalid value '%s' for number floats option 'nfs'. Must be integer.",option_value.c_str());
	  throw(error_message);
	}
      }
      else if(option_name == "radius") {
	if (!strIsInt(option_value.c_str(),&window_radius)) {
	  string error_message;
	  stringprintf(error_message,"Invalid value '%s' for window radius option 'radius'. Must be integer.",option_value.c_str());
	  throw(error_message);
	}
	if (obs->minPastFrames() < window_radius) {
	  obs->setMinPastFrames(window_radius);
	}
	if (obs->minFutureFrames() < window_radius) {
	  obs->setMinFutureFrames(window_radius);
	}
	if (window_radius > globalMinSkip) globalMinSkip = window_radius;
      }
      else if(option_name == "prior_1dpmf") {
	if (prior) {
	  error("ERROR: reading file '%s' line %d, DeepVirtualEvidenceCPT '%s' has multiple prior_1dpmf specifications\n",
		is.fileName(), is.lineNo(), name().c_str());
	}
	if (GM_Parms.dPmfsMap.find(option_value) == GM_Parms.dPmfsMap.end()) {
	  error("ERROR: reading file '%s' line %d,  DeepVirtualEvidenceCPT '%s' uses undefined Dense1DPMF '%s'\n", 
		is.fileName(), is.lineNo(), name().c_str(), option_value.c_str());
	}
	prior = GM_Parms.dPmfs[ GM_Parms.dPmfsMap[option_value] ];
	assert(prior);
	if (prior->card() != cardinalities[0]) {
	  error("ERROR: reading file '%s' line %d,  DeepVirtualEvidenceCPT '%s' needs a Dense1DPMF with cardinality %u, but '%s' has cardinality %u\n",
		is.fileName(), is.lineNo(), name().c_str(), cardinalities[0], prior->card());
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

    if (nfs == 0) {
      string error_message;
      stringprintf(error_message,"specifies %d floats, but must have > 0 floats", nfs);
      throw(error_message);
    }

    // check that it works with the current global observation matrix.
    if (nfs + obs_file_foffset > globalObservationMatrix->numContinuous()) {
      string error_message;
      stringprintf(error_message,"specifies %d floats at offset %d, but global observation matrix only has %d",
		   nfs,obs_file_foffset,globalObservationMatrix->numContinuous());
      throw(error_message);
    }


    // FXIME - check startSkip,endSkip >= radius

    // still here? Do more error checking.

    // input vector from the observation matrix is 
    // (1 + 2*window_radius) frames times nfs floats elements

    // final output layer size must match parent cardinality

    // check that DeepNN inputs matches the input vector size
    if (dmlp->numInputs() != (1 + 2 * window_radius) * nfs) {
      string error_message;
      stringprintf(error_message, "DeepVirtualEvidenceCPT '%s' requires DeepNN '%s' to have %u inputs, but it has %u", 
		   _name.c_str(), dmlp->name().c_str(), (1+2*window_radius)*nfs, dmlp->numInputs());
      throw(error_message);
    }

    // check that the final # outputs = parent cardinality
    if (dmlp->numOutputs() != cardinalities[0]) {
      string error_message;
      stringprintf(error_message,"DeepVirtualEvidenceCPT '%s' requires DeepNN '%s' to have %u outputs, but it has %u",
		   _name.c_str(), dmlp->name().c_str(), cardinalities[0], dmlp->numOutputs());
      throw(error_message);
    }

  } catch ( string const & error_message ) {
    error("ERROR: reading file '%s' line %u, of DeepVirtualEvidenceCPT spec '%s': %s",
	  is.fileName(),is.lineNo(),name().c_str(),error_message.c_str());
  } catch( const char * const error_message ) {
    error("ERROR: reading file '%s' line %u, of DeepVirtualEvidenceCPT spec '%s': %s",
	  is.fileName(),is.lineNo(),name().c_str(),error_message);
  }
  setBasicAllocatedBit();
}



/*-
 *-----------------------------------------------------------------------
 * DeepVECPT::write(os)
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
DeepVECPT::write(oDataStreamFile& os)
{

#define DeepVECPT_TMP_OUTPUT_STR_LEN 500

  char* tmp=new char[DeepVECPT_TMP_OUTPUT_STR_LEN];

  // write without any formatting.
  NamedObject::write(os);

  os.write("numParents:1");  // num parents
  sprintf(tmp,"parentCard:%u",cardinalities[0]);
  os.write(tmp);
  os.write("selfCard:2"); // self card
  sprintf(tmp,"deepNN:%s", dmlp->name().c_str());
  os.write(tmp);
  sprintf(tmp,"f_offset:%u",obs_file_foffset);
  os.write(tmp);
  sprintf(tmp,"nfs:%u",nfs);
  os.write(tmp);
  sprintf(tmp,"radius:%u",window_radius);
  os.write(tmp);
  os.write("END");
  os.nl();
  os.nl();
}


////////////////////////////////////////////////////////////////////
//        Probability routines (see GMTK_CPT.h for docs)
////////////////////////////////////////////////////////////////////


void 
DeepVECPT::becomeAwareOfParentValues( vector <RV *>& parents, const RV* rv ) {
  assert ( parents.size() == 1 );
  curParentValue = RV2DRV(parents[0])->val;
  assert(curParentValue <  cardinalities[0]);
}


logpr
DeepVECPT::applyNN(DiscRVType parentValue, DiscRV * drv) {
  REGISTER DiscRVType val = drv->val;

  logpr p((void*)NULL);

  // check if the CPT is cached
  unsigned frame = drv->frame();
  unsigned segment = obs->segmentNumber();
  if (frame == cached_frame && segment == cached_segment) {
    p.setFromP(cached_CPT[parentValue]);
    if (prior) {
      logpr priorP(prior->p(curParentValue));
      if (priorP.essentially_zero()) {
	infoMsg(IM::Inference, IM::Info, "DeepVECPT '%s' prior '%s' has zero value for %u\n",
	  name().c_str(), prior->name().c_str(), curParentValue);
	p.setFromLogP(LZERO);
      } else {
	p /= priorP;
      }
    }
    // The obseved value (1) is the one corresponding
    // to the value in the file. I.e., the score 
    // in the file corresponds to Pr(child = 1 | parent = j) = f_t(j), 
    // and f_t(j) is the value stored in the file.  
    if (val == 0) {
      // if the child RV has zero value, then we invert the probability.
      // see comment 'zero valued child case' elsewhere in this file.
      p = 1.0 - p;
    }
    return p;
  }

  // Not in the cache, so compute & cache it
  cached_frame = frame;
  cached_segment = segment;

  // assemble input vector
  assert(input_vector);
  input_vector[dmlp->numInputs()] = 1.0; // homogeneous coordinates
  float *dest = input_vector;
  // guarantees [frame - window_radius, frame + window_radius] are in cache
  unsigned stride = obs->stride();
  float *src  = obs->floatVecAtFrame(frame) - window_radius * stride + obs_file_foffset; 
  unsigned diameter = 1 + 2 * window_radius;
  for (unsigned i = 0; i < diameter; i+=1, src += stride, dest += nfs) {
    memcpy(dest, src, nfs * sizeof(float));
  }
  memcpy(cached_CPT, dmlp->applyDeepModel(input_vector), cardinalities[0] * sizeof(double)) ;

  // logpr the CPT entry
  assert(parentValue < cardinalities[0]);
  p.setFromP(cached_CPT[parentValue]);

  if (prior) {
    logpr priorP(prior->p(curParentValue));
    if (priorP.essentially_zero()) {
      infoMsg(IM::Inference, IM::Info, "DeepVECPT '%s' prior '%s' has zero value for %u\n",
        name().c_str(), prior->name().c_str(), curParentValue);
      p.setFromLogP(LZERO);
    } else {
      p /= priorP;
    }
  }
  // The obseved value (1) is the one corresponding
  // to the value in the file. I.e., the score 
  // in the file corresponds to Pr(child = 1 | parent = j) = f_t(j), 
  // and f_t(j) is the value stored in the file.  
  if (val == 0) {
    // if the child RV has zero value, then we invert the probability.
    // see comment 'zero valued child case' elsewhere in this file.
    p = 1.0 - p;
  }
  return p;
}

logpr 
DeepVECPT::probGivenParents(vector <RV *>& parents, DiscRV * drv) {
  assert ( bitmask & bm_basicAllocated );
  curParentValue = RV2DRV(parents[0])->val;
  assert(curParentValue <  cardinalities[0]);
  return applyNN(curParentValue, drv);
}


void DeepVECPT::becomeAwareOfParentValuesAndIterBegin(vector < RV *>& parents, 
						  iterator &it,
						  DiscRV* drv,
						  logpr& p) 
{
  becomeAwareOfParentValues(parents,drv);
  it.uInternalState = curParentValue;
  begin(it,drv,p);
}

void DeepVECPT::begin(iterator& it,DiscRV* drv, logpr& p) 
{
  assert ( bitmask & bm_basicAllocated );
  it.setCPT(this);
  it.drv = drv;
  it.uInternalState = curParentValue;
  drv->val = 0;
  assert(curParentValue <  cardinalities[0]);
  p = applyNN(curParentValue, drv);
}

bool DeepVECPT::next(iterator &it,logpr& p)
{
  if (it.drv->val == 1)
    return false;
  // here, the RV child now gets value 1.
  it.drv->val = 1;
  p = applyNN(curParentValue, it.drv);
  return true;
}



/*-
 *-----------------------------------------------------------------------
 * DeepVECPT::random sample(is)
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
int DeepVECPT::randomSample(DiscRV* drv)
{
  // sum up the observations in current frame.
  // TODO: this routine is currently not correct (but it is not used).
#if 0
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
#endif
  return 0;
}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#endif
