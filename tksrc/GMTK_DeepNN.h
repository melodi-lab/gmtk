/*-
 * GMTK_DeepNN.h
 *
 *  A "deep neural network" class. Can be used for deep VECPT or deep CPTs.
 *
 *  Written by Richard Rogers <rprogers@ee.washington.edu>
 * 
 * 
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_DEEPNN_H
#define GMTK_DEEPNN_H

#include "fileParser.h"
#include "logp.h"

#include "sArray.h"

#include "GMTK_NamedObject.h"
#include "GMTK_DoubleMatrix.h"


class DeepNN : public NamedObject {

 public:
  enum SquashFunction {
    SOFTMAX,
    LOGISTIC,
    TANH,
    ODDROOT,
    LINEAR,
    RECTLIN
  };

 private:
  unsigned max_outputs;
  unsigned num_inputs;
  unsigned num_matrices;

  vector<unsigned>       layer_output_count;
  vector<string>         layer_matrix_name;
  vector<string>         layer_squash_name;
  vector<SquashFunction> layer_squash_func;
  vector<float>          layer_logistic_beta;
  vector<DoubleMatrix *> layer_matrix;

public:

  ///////////////////////////////////////////////////////////  
  // General constructor, 

  DeepNN() : num_inputs(0), num_matrices(0) { }
 
  ~DeepNN() { }

  // Get NN outputs
  double *applyDeepModel(float *inputs);

  // Total number of inputs
  unsigned numInputs() { return num_inputs; } 

  // Number of outputs from the top layer
  unsigned numOutputs() { return layer_output_count[num_matrices-1]; }

  // Number of outputs from the specified layer
  unsigned layerOutputs(unsigned layer) { 
    assert(layer < layer_output_count.size());
    return layer_output_count[layer];
  }

  // Specified layer's squash function
  SquashFunction getSquashFn(unsigned layer) {
    assert(layer < layer_squash_func.size());
    return layer_squash_func[layer];
  }

  // Set the parameters of the specified layer. ld is the weights' stride
  void setParams(unsigned layer, double const *weights, unsigned ld, double const *bias) {
    assert(layer < layer_matrix.size());
    DoubleMatrix *w = layer_matrix[layer];
    unsigned rows = w->_rows;
    unsigned cols = w->_cols;
    for (unsigned r=0; r < rows; r+=1) {
      unsigned c;
      for (c=0; c < cols - 1; c+=1) {
	// weights come in column-major order but transposed; store them in row-major
	w->values[ r * cols + c ] = weights[ r * ld + c ];
      }
      w->values[ r * cols + c ] = bias[r];
    }
  }

  unsigned numLayers() { return num_matrices; }

  vector<DoubleMatrix *> &getMatrices() { return layer_matrix; }

  float getBeta(unsigned layer) {
    assert(layer < layer_logistic_beta.size());
    assert(layer_squash_func[layer] == LOGISTIC);
    return layer_logistic_beta[layer];
  }


  ////////////////////////////////////////////////////////////////////
  // How many parameters does this consume? We return 0 since the
  // virtual evidence does not constitute parameters in the normal
  // sense.
  unsigned totalNumberParameters() { 
    unsigned count = 0;
    for (unsigned layer=0; layer < num_matrices; layer+=1) {
      DoubleMatrix *w = layer_matrix[layer];
      count += w->_rows * w->_cols;
    }
    return count;
  }

  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  const string typeName() { return "DeepNN"; }
  //////////////////////////////////

};


#endif // defined DeepNN
