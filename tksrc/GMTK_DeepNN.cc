/*-
 *
 * GMTK_DeepNN.cc
 *     Deep Neural Network.
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes 
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 * An example of how you'd specify a DeepNN:
 * (in a masterfile:)


DOUBLE_MAT_IN_FILE inline <m>
0
matrixName0
<rows>              % # of outputs
<cols>              % # of inputs +1 for bias vector
<value>...<value>
...

DEEP_NN_IN_FILE inline 1
0
NNName
<input card> <output card>
matrices:m M_0 M_1 ... M_{m-1}  % m is # of weight layers (matrices)
                                % matrix0 must be M_0 rows by 1 + input_card columns
                                % matrixX must be M_X rows by 1 + M_{X-1} columns
                                % M_{m-1} = <output card> rows by 1 + M_{m-2} columns 
matrix0:matrixName0  
squash0:<squash fn>  % softmax, logistic(b), tanh, odd root, etc.
matrix1:matrixName1
squash1:<squash fn>
...
matrix<m-1>:matrixName<m-1>
squash<m-1>:<squash fn>
END

Note the "EOF" at the end of each DeepNN.  Also, parsing of the DeepNNs is
rudimentary and no spaces are allowed around the ":".

Also the order in which the arguments are written is not important.

The default values are:

matrices:0 % will cause an error if not updated

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

#include "GMTK_DeepNN.h"
#include "GMTK_GMParms.h"
#include "GMTK_DiscRV.h"
#include "GMTK_FileSource.h"

VCID(HGID)


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////

 

/*-
 *-----------------------------------------------------------------------
 * DeepNN::read(is)
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
DeepNN::read(iDataStreamFile& is)
{

  //Initialize defaults
  

  string str;
  string option_name;
  string option_value;
  unsigned lineNum=0;

  // read the name of the object.
  NamedObject::read(is);

  is.read(num_inputs, "Can't read DeepNN's num inputs");
  unsigned num_outputs;
  is.read(num_outputs, "Can't read DeepNN's num outputs");

  try {

    ///////////////////////////////////////////////////////////////////////
    // Next, we have a set of arguments that the user may give.
    // These arguments consist of lines of a "flag : value" syntax, where "flag"
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
      
      if(option_name == "matrices") {
	if (!strIsInt(option_value.c_str(),&num_matrices)) {
	  string error_message;
	  stringprintf(error_message,"Invalid value '%s' for number of matrices option 'matrices. Must be integer.", option_value.c_str());
	  throw(error_message);
	}
	if (num_matrices == 0) {
	  string error_message;
	  stringprintf(error_message,"DeepNN requires at least 1 matrix");
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
			 layer, str.c_str());
	    throw(error_message);
	  }
	  layer_matrix_name[layer] = option_value;
	  if (GM_Parms.doubleMatsMap.find(option_value) == GM_Parms.doubleMatsMap.end()) {
	    string error_message;
	    stringprintf(error_message,"DeepNN '%s' specifies double matrix name '%s' that does not exist",
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

    if (num_inputs == 0) {
      string error_message;
      stringprintf(error_message,"specifies %u inputs, but must have > 0 inputs", num_inputs);
      throw(error_message);
    }

    // matrix rows = # of outputs for the layer
    // matrix cols = # of inputs for the layer (+ 1 for bias vector)
    
    // check that layer 0 matches the input vector size (+1 for bias element)
    if (layer_matrix[0]->cols() != (int)num_inputs+1) {
      string error_message;
      stringprintf(error_message, "DeepNN '%s' requires matrix0 to have %u columns, but double matrix '%s' has %d columns", 
		   _name.c_str(), num_inputs + 1, layer_matrix_name[0].c_str(), layer_matrix[0]->cols());
      throw(error_message);
    }

    // check that the # rows for each layer's matrix matches the declared
    // output size in the matrix:... option   Also find the maximum # of outputs
    max_outputs = 0;
    for (unsigned i=0; i < num_matrices; i+=1) {
      if (layer_matrix[i]->rows() != (int)layer_output_count[i]) {
	string error_message;
	stringprintf(error_message,"DeepNN '%s' expects double matrix '%s' to have %u rows, but it has %d",
		     _name.c_str(), layer_matrix_name[i].c_str(), layer_output_count[i], layer_matrix[i]->rows());
	throw(error_message);
      }
      if (layer_output_count[i] > max_outputs)
	max_outputs = layer_output_count[i];
    }
    unsigned vector_length = max_outputs > num_inputs ? max_outputs : num_inputs;
    output_vector[0] = new double[vector_length+1]; // big enough to hold any layer's (in|out)put (+1 for homogeneous coordinates)
    output_vector[1] = new double[vector_length+1];

    // now check that the input vector size for layer i matches the
    // output vector size from layer i-1
    for (unsigned i=1; i < num_matrices; i+=1) {
      if (layer_matrix[i]->cols() != (int)layer_output_count[i-1] + 1) {
	string error_message;
	stringprintf(error_message,"DeepNN '%s' expects double matrix '%s' to have %u cols, but it has %d",
		     _name.c_str(), layer_matrix_name[i].c_str(), layer_output_count[i-1]+1, layer_matrix[i]->cols());
	throw(error_message);
      }
    }
    // check that the output layer's # outputs matches
    if (layer_matrix[num_matrices-1]->rows() != (int)num_outputs) {
      string error_message;
      stringprintf(error_message,"DeepNN '%s' final output matrix '%s' must have %u rows to match specified number of outputs, but it has %d",
		   _name.c_str(), layer_matrix_name[num_matrices-1].c_str(), num_outputs, layer_matrix[num_matrices-1]->rows());
      throw(error_message);
    }

  } catch ( string error_message ) {
    error("ERROR: reading file '%s' line %u, of DeepNN spec '%s': %s",
	  is.fileName(),is.lineNo(),name().c_str(),error_message.c_str());
  } catch( const char * const error_message ) {
    error("ERROR: reading file '%s' line %u, of DeepNN spec '%s': %s",
	  is.fileName(),is.lineNo(),name().c_str(),error_message);
  }
//  setBasicAllocatedBit();
}



/*-
 *-----------------------------------------------------------------------
 * DeepNN::write(os)
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
DeepNN::write(oDataStreamFile& os)
{
#if 0
#define DeepNN_TMP_OUTPUT_STR_LEN 500

  char* tmp=new char[DeepNN_TMP_OUTPUT_STR_LEN];

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
#endif
}

void static
rectlin(double *q, unsigned len) {
  for (unsigned i=0; i < len; i+=1)
    if (q[i] < 0.0) q[i] = 0.0;
}

void static
softmax(double *q, unsigned len) {
  double k  = *q;
#if 1
  // Unrolled max by 2 was fastest of a few implementations I tried
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
  for (unsigned i=0; i < len; i+=1) {
    q[i] = exp ( q[i] - k - sum );
  }
}

void static
logistic(double *q, unsigned len, float beta=1.0) {
  for (unsigned i=0; i < len; i+=1)
    if (q[i] < -30)
      q[i] = 0;
    else if (q[i] > 30)
      q[i] = 1;
    else
      q[i] = 1.0 / (1.0 + exp( -beta * q[i] ));
}

void static
hyptan(double *q, unsigned len) {
  for (unsigned i=0; i < len; i+=1) 
    q[i] = tanh(q[i]);
}

void static 
oddroot(double *q, unsigned len) {
#if 0
// Mathematica's inverse of $\frac{x^3}{3} + x$
#define CUBE_ROOT_OF_2 1.25992104989
  for (unsigned i=0; i < len; i+=1) {
    double x = pow( 3 * q[i] + sqrt( 4.0 + 9.0 * q[i] * q[i] ), 1.0/3.0 );
    q[i] = x / CUBE_ROOT_OF_2 - CUBE_ROOT_OF_2 / x;
  }
}
#else
// Galen's faster implementation based on Newton's method 
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
#endif
}


void static
squash(DeepNN::SquashFunction fn, double *q, unsigned len, float beta=1.0) {
  switch (fn) {
  case DeepNN::SOFTMAX:
    softmax(q, len);
    break;
  case DeepNN::LOGISTIC:
    logistic(q, len, beta);
    break;
  case DeepNN::TANH:
    hyptan(q, len);
    break;
  case DeepNN::ODDROOT:
    oddroot(q, len);
    break;
  case DeepNN::LINEAR:
    // q = q
    break;
  case DeepNN::RECTLIN:
    rectlin(q, len);
    break;
  default:
    assert(0); // impossible
  }
}


double *
DeepNN::applyDeepModel(float *inputs) {
  assert(output_vector[0] && output_vector[1]);
  double *input_vector = output_vector[1];
  input_vector[num_inputs] = 1.0;                  // homogeneous coordinates
  for (unsigned i = 0; i < num_inputs; i+=1) {
    input_vector[i] = (double)inputs[i];
  }
  memset(output_vector[0], 0, (max_outputs+1) * sizeof(double));
  mul_mdmd_md(layer_output_count[0], num_inputs+1, 1, 
	      layer_matrix[0]->values.ptr, input_vector, output_vector[0], 
	      num_inputs+1, 1, 1);
  squash(layer_squash_func[0], output_vector[0], layer_output_count[0], layer_logistic_beta[0]);
  output_vector[0][layer_output_count[0]] = 1.0;

  unsigned cur_output_vector = 0;
  for (unsigned layer=1; layer < num_matrices; layer += 1) {
    input_vector = output_vector[cur_output_vector];
    cur_output_vector = (cur_output_vector + 1) % 2;
    memset(output_vector[cur_output_vector], 0, (max_outputs+1) * sizeof(double));
    mul_mdmd_md(layer_output_count[layer], layer_output_count[layer-1]+1, 1, 
		layer_matrix[layer]->values.ptr, input_vector, output_vector[cur_output_vector], 
		layer_output_count[layer-1]+1, 1, 1);
    squash(layer_squash_func[layer], output_vector[cur_output_vector], layer_output_count[layer], layer_logistic_beta[layer]);
    output_vector[cur_output_vector][layer_output_count[layer]] = 1.0;
  }
  return output_vector[cur_output_vector];
}

////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#endif
