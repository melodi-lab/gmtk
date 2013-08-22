/*-
 *
 * GMTK_DeepVECPT.cc
 *     Deep Virtual Evidence CPT.
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes 
 * Licensed under the Open Software License version 3.0
 *
 * An example of how you'd specify DeepVECPTs:
 * (in a masterfile:)


DOUBLE_MAT_IN_FILE inline <m>
0
matrixName0
<rows>              % # of outputs
<cols>              % # of inputs +1 for bias vector
<value>...<value>
...

DEEP_VE_CPT_IN_FILE inline 1
0
modelName
1                    % number of parents (always 1)
<parent card>        % must match final layer's # outputs
2                    % self cardinality (always 2)
f_offset:0           % starting index in observation matrix for floats
nfs:0                % # of floats to take from observation matrix
radius:<d>           % # of frames to take from observation matrix = 2d + 1
matrices:m M_0 M_1 ... M_{m-1}  % m is # of weight layers (matrices)
                                % matrix0 must be M_0 rows by 1 + nfs(2 radius + 1) columns
                                % matrixX must be M_X rows by 1 + M_{X-1} columns
                                % M_{m-1} = <parent card> rows by 1 + M_{m-2} columns 
matrix0:matrixName0  
squash0:<squash fn>  % softmax, logistic(b), tanh, odd root, etc.
matrix1:matrixName1
squash1:<squash fn>
...
matrix<m-1>:matrixName<m-1>
squash<m-1>:<squash fn>
END

Note the "EOF" at the end of each DeepVECPT.  Also, parsing of the DeepVECPTs is
rudimentary and no spaces are allowed around the ":".

Also the order in which the arguments are written is not important.

The default values are:

matrices:0 % will cause an error if not updated
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

  // read the name of the object.
  NamedObject::read(is);
  is.read(_numParents,"Can't read DeepVirtualEvidenceCPT's num parents");
  if (_numParents != 1)
    error("ERROR: reading file '%s' line %d, DeepVirtualEvidenceCPT '%s' must have exactly one(1) rather than %d parents",
          is.fileName(),is.lineNo(),name().c_str(),_numParents);
  cardinalities.resize(_numParents);
  // read the parent cardinality
  is.read(cardinalities[0],"Can't read DeepVirtualEvidenceCPT's parent cardinality");
  cached_CPT = new double[cardinalities[0]];
  // read the self cardinality, must be binary
  is.read(_card,"Can't read DeepVirtualEvidenceCPT's self cardinality");
  if (_card != 2)
    error("ERROR: reading file '%s' line %d, VirtualEvidenceCPT '%s', child cardinality must be two(2) rather than %d.",
          is.fileName(),is.lineNo(),name().c_str(),_card);

  try {

    if(cardinalities[0]==0) {
      string error_message;
      stringprintf(error_message,"must supply parent cardinality");
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
      }
      else if(option_name == "matrices") {
	if (!strIsInt(option_value.c_str(),&num_matrices)) {
	  string error_message;
	  stringprintf(error_message,"Invalid value '%s' for number of matrices option 'matrices. Must be integer.", option_value.c_str());
	  throw(error_message);
	}
	if (num_matrices == 0) {
	  string error_message;
	  stringprintf(error_message,"DeepVirtualEvidenceCPT requires at least 1 matrix");
	  throw(error_message);
	}
	layer_output_count.resize(num_matrices);
	layer_matrix_name.resize(num_matrices);
	layer_squash_name.resize(num_matrices);
	layer_squash_func.resize(num_matrices);
	layer_logistic_beta.resize(num_matrices);
	layer_matrix.resize(num_matrices);

	if (!is.read(layer_output_count, num_matrices)) {
	  string error_message;
	  stringprintf(error_message,"Invalid matrix row counts for option 'matrices'. Must be %u integers.", num_matrices);
	  throw(error_message);
	}
	for (unsigned layer=0; layer < num_matrices; layer+=1) {
	  is.read(str);
	  // parse the string we have just read
	  string::size_type len=str.length();
	  string::size_type pos=str.find(":",0);
	  if (pos == string::npos || len < 7 || pos == 0) {
	    string error_message;
	    stringprintf(error_message,"Invalid format '%s' which should be of form 'matrix%u:name' where 'name' is the name of a double matrix.",
			 str.c_str(), layer);
	    throw(error_message);
	  }
	  option_name = str.substr(0,pos);
	  option_value = str.substr(pos+1,len-pos-1);
	  string matrixNum;
	  stringprintf(matrixNum,"matrix%u", layer);
	  if (option_name != matrixNum) {
	    string error_message;
	    stringprintf(error_message, "Expected 'matrix%u:name' but got '%s'",
			 layer, layer, str.c_str());
	    throw(error_message);
	  }
	  layer_matrix_name[layer] = option_value;
	  if (GM_Parms.doubleMatsMap.find(option_value) == GM_Parms.doubleMatsMap.end()) {
	    string error_message;
	    stringprintf(error_message,"DeepVirtualEvidenceCPT '%s' specifies double matrix name '%s' that does not exist",
			 _name.c_str(), option_value.c_str());
	    throw(error_message);
	  }
	  unsigned matrix_index = GM_Parms.doubleMatsMap[option_value];
	  layer_matrix[layer] = GM_Parms.doubleMats[matrix_index];

	  is.read(str);
	  // parse the string we have just read
	  len=str.length();
	  pos=str.find(":",0);
	  if (pos == string::npos || len < 7 || pos == 0) {
	    string error_message;
	    stringprintf(error_message,"Invalid format '%s' which should be of form 'matrix%u:name' where 'name' is the name of a double matrix.", 
			 str.c_str(), layer);
	    throw(error_message);
	  }
	  option_name = str.substr(0,pos);
	  option_value = str.substr(pos+1,len-pos-1);
	  string squashNum;
	  stringprintf(squashNum,"squash%u", layer);
	  if (option_name != squashNum) {
	    string error_message;
	    stringprintf(error_message, "Expected 'squash%u:function but got '%s'",
			 layer, str.c_str());
	    throw(error_message);
	  }
	  layer_squash_name[layer] = option_value;
	  if (option_value == "softmax") {
	    layer_squash_func[layer] = SOFTMAX;
	  } else if (option_value.compare(0, 8, "logistic") == 0) {
	    layer_squash_func[layer] = LOGISTIC;
	    if (option_value[8] == '(') {
	      string beta_str = option_value.substr(9, option_value.find(')', 9)-9);
	      int float_len;
	      if (!strIsFloat(beta_str.c_str(),&(layer_logistic_beta[layer]), &float_len) ||
		  float_len != (int)beta_str.length()) 
	      {
		string error_message;
		stringprintf(error_message,"Invalid value '%s' for logistic beta. Must be float.", beta_str.c_str());
		throw(error_message);
	      }
	    } else {
	      layer_logistic_beta[layer] = 1.0;
	    }
	  } else if (option_value == "tanh") {
	    layer_squash_func[layer] = TANH;
	  } else if (option_value == "oddroot") {
	    layer_squash_func[layer] = ODDROOT;
	  } else if (option_value == "linear") {
	    layer_squash_func[layer] = LINEAR;
	  } else if (option_value == "rectlin") {
	    layer_squash_func[layer] = RECTLIN;
	  } else {
	    string error_message;
	    stringprintf(error_message, "Invalid squash function '%s', must be one of 'softmax', 'logistic', 'tanh', "
			 " 'oddroot', 'linear', or 'rectlin'", option_value.c_str());
	    throw(error_message);
	  }
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


    // matrix rows = # of outputs for the layer
    // matrix cols = # of inputs for the layer (+ 1 for bias vector)
    
    // input vector from the observation matrix is 
    // (1 + 2*window_radius) frames times nfs floats elements

    // final output layer size must match parent cardinality

    // check that layer 0 matches the input vector size (+1 for bias element)
    if (layer_matrix[0]->cols() != (int)((1 + 2 * window_radius) * nfs + 1)) {
      string error_message;
      stringprintf(error_message, "DeepVirtualEvidenceCPT '%s' requires matrix0 to have %u columns, but double matrix '%s' has %d columns", 
		   _name.c_str(), (1+2*window_radius)*nfs + 1, layer_matrix_name[0].c_str(), layer_matrix[0]->cols());
      throw(error_message);
    }

    // check that the # rows for each layer's matrix matches the declared
    // output size in the matrix:... option   Also find the maximum # of outputs
    max_outputs = 0;
    for (unsigned i=0; i < num_matrices; i+=1) {
      if (layer_matrix[i]->rows() != (int)layer_output_count[i]) {
	string error_message;
	stringprintf(error_message,"DeepVirtualEvicenceCPT '%s' expects double matrix '%s' to have %u rows, but it has %d",
		     _name.c_str(), layer_matrix_name[i].c_str(), layer_output_count[i], layer_matrix[i]->rows());
	throw(error_message);
      }
      if (layer_output_count[i] > max_outputs)
	max_outputs = layer_output_count[i];
    }

    // now check that the input vector size for layer i matches the
    // output vector size from layer i-1
    for (unsigned i=1; i < num_matrices; i+=1) {
      if (layer_matrix[i]->cols() != (int)layer_output_count[i-1] + 1) {
	string error_message;
	stringprintf(error_message,"DeepVirtualEvidenceCPT '%s' expects double matrix '%s' to have %u cols, but it has %d",
		     _name.c_str(), layer_matrix_name[i].c_str(), layer_output_count[i-1]+1, layer_matrix[i]->cols());
	throw(error_message);
      }
    }
    // check that the final # outputs = parent cardinality
    if (layer_matrix[num_matrices-1]->rows() != (int)cardinalities[0]) {
      string error_message;
      stringprintf(error_message,"DeepVirtualEvidenceCPT '%s' final output matrix '%s' must have %u rows to match the parent variable's cardinality, but it has %d",
		   _name.c_str(), layer_matrix_name[num_matrices-1].c_str(), cardinalities[0], layer_matrix[num_matrices-1]->rows());
      throw(error_message);
    }

  } catch ( string error_message ) {
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
  sprintf(tmp,"f_offset:%u",obs_file_foffset);
  os.write(tmp);
  sprintf(tmp,"nfs:%u",nfs);
  os.write(tmp);
  sprintf(tmp,"radius:%u",window_radius);
  os.write(tmp);
  // FIXME - write matrices: stanza
  os.write("END");
  os.nl();
  os.nl();
}


////////////////////////////////////////////////////////////////////
//        Probability routines (see GMTK_CPT.h for docs)
////////////////////////////////////////////////////////////////////


void DeepVECPT::becomeAwareOfParentValues( vector <RV *>& parents,
				       const RV* rv ) 
{
  assert ( parents.size() == 1 );
  curParentValue = RV2DRV(parents[0])->val;
}

void
rectlin(double *q, unsigned len) {
  for (unsigned i=0; i < len; i+=1)
    if (q[i] < 0.0) q[i] = 0.0;
}

void
softmax(double *q, unsigned len) {
  double k  = *q;
#if 1
  unsigned n;
  if (len % 2) {
    n = 1;
  } else {
    if (*q < q[1]) {
      k = q[1];
    }
    n = 2;
  }
  double k2 = k;
  for (unsigned i=n; i < len; i+=2) {
    k  = (k  > q[i])   ? k  : q[i];
    k2 = (k2 > q[i+1]) ? k2 : q[i+1];
  }
  k = (k > k2) ? k : k2;
#else
  for (unsigned i=1; i < len; i+=1)
    if (k < q[i])
      k = q[i];
#endif
  double sum = 0.0;
  for (unsigned i=0; i < len; i+=1)
    sum += exp(q[i] - k);
  sum = log(sum);
  for (unsigned i=0; i < len; i+=1)
    q[i] = exp ( q[i] - k - sum );
}

void
logistic(double *q, unsigned len, float beta=1.0) {
  for (unsigned i=0; i < len; i+=1)
    if (q[i] < -30)
      q[i] = 0;
    else if (q[i] > 30)
      q[i] = 1;
    else
      q[i] = 1.0 / (1.0 + exp( -beta * q[i] ));
}

void
hyptan(double *q, unsigned len) {
  for (unsigned i=0; i < len; i+=1) 
    q[i] = tanh(q[i]);
}

#if 0
#define CUBE_ROOT_OF_2 1.25992104989
void
oddroot(double *q, unsigned len) {
  for (unsigned i=0; i < len; i+=1) {
    double x = pow( 3 * q[i] + sqrt( 4.0 + 9.0 * q[i] * q[i] ), 1.0/3.0 );
    q[i] = x / CUBE_ROOT_OF_2 - CUBE_ROOT_OF_2 / x;
  }
}
#else
void
oddroot(double *q, unsigned len) {
  for (unsigned i=0; i < len; i+=1) {
    double y = q[i];
    bool negate = false;
    if (y < 0) { negate = true; y = -y; }
    double x = (y <= 20) ? y : pow(3.0 * y, 1.0/3.0);

    double newX;
    while (true) {
      double xSqr = x * x;
      newX = (0.66666666666666666 * xSqr * x + y) / (xSqr + 1.0);
      if (newX >= x) break;
      x = newX;
    }
    q[i] = negate ? -newX : newX;
  }
}
#endif


void
squash(DeepVECPT::SquashFunction fn, double *q, unsigned len, float beta=1.0) {
  switch (fn) {
  case DeepVECPT::SOFTMAX:
    softmax(q, len);
    break;
  case DeepVECPT::LOGISTIC:
    logistic(q, len, beta);
    break;
  case DeepVECPT::TANH:
    hyptan(q, len);
    break;
  case DeepVECPT::ODDROOT:
    oddroot(q, len);
    break;
  case DeepVECPT::LINEAR:
    // q = q
    break;
  case DeepVECPT::RECTLIN:
    rectlin(q, len);
    break;
  default:
    assert(0); // impossible
  }
}


logpr
DeepVECPT::applyDeepModel(DiscRVType parentValue, DiscRV * drv) {
  register DiscRVType val = drv->val;

  logpr p((void*)NULL);

  unsigned frame = drv->frame();

  if (frame == cached_frame && obs->segmentNumber() == cached_segment) {
    p.valref() = cached_CPT[curParentValue];
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

  cached_frame = frame;
  cached_segment = obs->segmentNumber();

  unsigned num_inputs = nfs * ( 2 * window_radius + 1 ) + 1;
  double *input_vector = new double[num_inputs];
  input_vector[num_inputs-1] = 1.0; // homogeneous coordinates
  double *dest = input_vector;
  // guarantees [frame - window_radius, frame + window_radius] are in cache
  unsigned stride = obs->stride();
  float *src  = obs->floatVecAtFrame(frame) - window_radius * stride + obs_file_foffset; 
  for (unsigned i = 0; i < 1 + 2 * window_radius; i+=1, src += stride, dest += nfs) {
    for (unsigned j=0; j < nfs; j+=1) {
      dest[j] = (double)(src[j]);
    }
  }
#if 0
printf("%02u:", frame);
for(unsigned i=0; i < num_inputs; i+=1)
  printf(" %f", input_vector[i]);
#endif
  double *output_vector[2];
  output_vector[0] = new double[max_outputs+1]; // big enough to hold any layer's output (+1 for homogeneous coordinates)
  output_vector[1] = new double[max_outputs+1];

  mul_mdmd_md(layer_output_count[0], num_inputs, 1, 
	      layer_matrix[0]->values.ptr, input_vector, output_vector[0], 
	      num_inputs, 1, 1);
  delete[] input_vector;
  squash(layer_squash_func[0], output_vector[0], layer_output_count[0], layer_logistic_beta[0]);
  output_vector[0][layer_output_count[0]] = 1.0;

  unsigned cur_output_vector = 0;
  for (unsigned layer=1; layer < num_matrices; layer += 1) {
    input_vector = output_vector[cur_output_vector];
    cur_output_vector = (cur_output_vector + 1) % 2;
    mul_mdmd_md(layer_output_count[layer], layer_output_count[layer-1]+1, 1, 
		layer_matrix[layer]->values.ptr, input_vector, output_vector[cur_output_vector], 
		layer_output_count[layer-1]+1, 1, 1);
    squash(layer_squash_func[layer], output_vector[cur_output_vector], layer_output_count[layer], layer_logistic_beta[layer]);
    output_vector[cur_output_vector][layer_output_count[layer]] = 1.0;
  }
  memcpy(cached_CPT, output_vector[cur_output_vector], parentCardinality(0) * sizeof(double));
#if 0
printf(" ->");
for (unsigned i=0; i < parentCardinality(0); i+=1)
  printf(" %f", output_vector[cur_output_vector][i]);
printf("\n");
#endif
  p.valref() = output_vector[cur_output_vector][curParentValue];
  // The obseved value (1) is the one corresponding
  // to the value in the file. I.e., the score 
  // in the file corresponds to Pr(child = 1 | parent = j) = f_t(j), 
  // and f_t(j) is the value stored in the file.  
  if (val == 0) {
    // if the child RV has zero value, then we invert the probability.
    // see comment 'zero valued child case' elsewhere in this file.
    p = 1.0 - p;
  }
  delete[] output_vector[0];
  delete[] output_vector[1];
  return p;
}

logpr 
DeepVECPT::probGivenParents(vector <RV *>& parents, DiscRV * drv) {
  assert ( bitmask & bm_basicAllocated );
  curParentValue = RV2DRV(parents[0])->val;
  return applyDeepModel(curParentValue, drv);

#if 0
  register DiscRVType val = drv->val;

  logpr p((void*)NULL);

  unsigned num_inputs = nfs * ( 2 * window_radius + 1 ) + 1;
  float *input_vector = new float[num_inputs];
  input_vector[num_inputs-1] = 1.0; // homogeneous coordinates
  unsigned frame = drv->frame();
  float *dest = input_vector;
  for (unsigned t = frame - window_radius; t <= frame + window_radius; t+=1, dest += nfs) {
    memcpy(dest, obs->floatVecAtFrame(t), nfs * sizeof(float));
  }

  float *output_vector[2];
  output_vector[0] = new float[max_outputs]; // big enough to hold any layer's output (+1 for homogeneous coordinates)
  output_vector[1] = new float[max_outputs];

  mul_mfmf_mf(layer_output_count[0], num_inputs, 1, 
	      layer_matrix[0]->values.ptr, input_vector, output_vector[0], 
	      num_inputs, 1, 1);
  delete[] input_vector;
  squash(layer_squash_func[0], output_vector[0], layer_output_count[0]);
  output_vector[0][layer_output_count[0]] = 1.0;

  unsigned cur_output_vector = 0;
  for (unsigned layer=1; layer < num_matrices; layer += 1) {
    input_vector = output_vector[cur_output_vector];
    cur_output_vector = (cur_output_vector + 1) % 2;
    mul_mdmd_md(layer_output_count[layer], layer_output_count[layer-1]+1, 1, 
		layer_matrix[layer]->values.ptr, input_vector, output_vector[cur_output_vector], 
		layer_output_count[layer-1]+1, 1, 1);
    squash(layer_squash_func[layer], output_vector[cur_output_vector], layer_output_count[layer]);
    output_vector[cur_output_vector][layer_output_count[layer]] = 1.0;
  }
  p.valref() = output_vector[cur_output_vector][curParentValue];
  // The obseved value (1) is the one corresponding
  // to the value in the file. I.e., the score 
  // in the file corresponds to Pr(child = 1 | parent = j) = f_t(j), 
  // and f_t(j) is the value stored in the file.  
  if (val == 0) {
    // if the child RV has zero value, then we invert the probability.
    // see comment 'zero valued child case' elsewhere in this file.
    p = 1.0 - p;
  }
  delete[] output_vector[0];
  delete[] output_vector[1];
  return p;
#endif
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
  p = applyDeepModel(curParentValue, drv);
}

bool DeepVECPT::next(iterator &it,logpr& p)
{
  if (it.drv->val == 1)
    return false;
  // here, the RV child now gets value 1.
  it.drv->val = 1;
  p = applyDeepModel(curParentValue, it.drv);
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
