#pragma once

/*-
 * DBN.h
 *     Deep Belief Network training code
 *
 * Written by Galen Andrew gmandrew@uw.edu
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#if defined(HAVE_MKL)
#  include "mkl.h"
#  include "mkl_lapacke.h"
#  include "mkl_spblas.h"
#  include "mkl_trans.h"
#elif defined(HAVE_BLAS)
extern "C" {            /* Assume C declarations for C++ */
#  include <cblas.h>
}
#endif

#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>      // std::setprecision

#include "rand.h"
#include "fileParser.h"
#include "debug.h"

#include "FileBackedMatrix.h"
#include "MMapMatrix.h"
#include "StdioMatrix.h"
#include "Layer.h"
#include "MatrixFunc.h"
#include "BatchSource.h"

using namespace std;

#define MIN_STEP_SIZE (1.0e-20)

/*
The DBN class maintains all of the parameters of the DBN, and includes operations
for training and testing.
 */

class DBN {
public:

  static bool     resumeTraining;  // should this invocation resume where the last left off?
  static bool     checkSignal;     // if true, print status as if end of check interval
  static bool     randomInitLayer; // randomly initialize layers if true, else use existing weights from master file
  static bool     sparseInitLayer; // select sparse or dense layer initialization method
  static unsigned nnChunkSize;     // O(1) space to use for (non-minibatch) incremental processing

  enum PretrainType {
    NONE, // no pretraining
    AE, // denoising autoencoder
    CD  // contrastive divergence
  };

	// The type of training objective, either minimize squared error, or softmax regression
  enum ObjectiveType { SQ_ERR, SOFT_MAX };

	// structure for storing all hyperparameters of pretraining or backpropagation
  struct HyperParams {
    double initStepSize, // step size, constant unless no progress observed, then decreased
			maxMomentum, // momentum reaches this value by the end of training (before annealing)
			minMomentum, // initial momentum value
			maxUpdate, // maximum l2 norm of update (including momentum). Updates truncated if larger than this
			l2; // l2 regularization, aka weight decay
    double iDropP, // dropout probability for input layer
			hDropP; // dropout probability for hidden layers
    int numUpdates, // total number of updates before annealing phase
			numAnnealUpdates, // number of updates in annealing phase
		                    // during which momentum and step size both decrease to zero 
			miniBatchSize, // number of instances in each minibatch
			checkInterval; // number of updates between console outputs

    PretrainType pretrainType;

    HyperParams() :
      initStepSize(1e-2),
      maxMomentum(0.99),
      minMomentum(0.5),
      maxUpdate(0.1),
      l2(1e-3),
      iDropP(0),
      hDropP(0),
      numUpdates(10000),
      numAnnealUpdates(2000),
      miniBatchSize(10),
      checkInterval(2000),
      pretrainType(CD)
    { }

    void print() const {
      printf("initStepSize     %f\n"
	     "minMomentum      %f\n"
             "maxMomentum      %f\n"
             "maxUpdate        %f\n"
             "l2               %f\n"
             "iDropP           %f\n"
             "hDropP           %f\n"
             "numUpdates       %d\n"
             "numAnnealUpdates %d\n"
	     "miniBatchSize    %d\n"
	     "checkInterval    %d\n",
	     initStepSize, minMomentum, maxMomentum, maxUpdate, l2, iDropP, hDropP,
	     numUpdates, numAnnealUpdates, miniBatchSize, checkInterval);
    }
  };

private:
  // width of input, all hidden, and output layers
  int _iSize, _oSize;
  vector<int> _hSize;

	// activation function of input layer
  Layer::ActFunc _iActFunc;

	// activation function of hidden layers
  vector<Layer::ActFunc> _hActFunc;
  
  // number of layers (not including input)
  int _numLayers;

	// layers
  mutable vector<Layer> _layers;

  // all model parameters
  AllocatingVector _params, _deltaParams, _savedParams;

  vector<MutableMatrix> _W; // weight matrices for each layer
  vector<MutableVector> _B; // biases for each layer
  vector<MutableMatrix> _deltaW; // cumulative weight update for each layer (with momentum)
  vector<MutableVector> _deltaB; // cumulative bias update for each layer (with momentum)
  vector<MutableMatrix> _savedW; // saved weights
  vector<MutableVector> _savedB; // saved biases

	// vectors pointing to all parameters _W and _B for each layer (and deltas and saved)
  vector<MutableVector> _layerParams, _layerDeltaParams, _layerSavedParams;

  // space for temporary values during CD training
  mutable AllocatingMatrix _tempTopSample, _tempBottomSample;

  // temporary storage used in several places
  mutable AllocatingMatrix _tempMat, _tempDropoutInput;
  mutable AllocatingVector _tempVec;
  mutable Layer _tempTopLayer, _tempBottomLayer;

	// apply weight decay and return weighted squared L2 norm
  static double Decay(MutableVector & v, double alpha, double step) {
    double val = 0.5 * alpha * (v * v);
		// decay amount is scaled by step size
    if (step * alpha > 0) v *= (1.0 - step * alpha);
    return val;
  }

	// initialize the weights and biases
  void InitLayer(int layer) {
    if (sparseInitLayer) {
      // sparse initialization strategy from Martens, 2010
      _B[layer] *= 0;
      MutableMatrix W = _W[layer];
      W *= 0;
      for (int c = 0; c < W.NumC(); ++c) {
	MutableVector col = W.GetCol(c);
	// sampling with replacement, but it doesn't really matter
	for (int i = 0; i < 15; ++i) {
	  col[rnd.uniformOpen(col.Len())] = rnd.normal() * 0.01;
	}
      }
    } else {
      // dense initialization strategy from Glorot & Bengio 2010
      _B[layer] *= 0;
      double max = 1.0 / sqrt(_W[layer].NumR());
      _W[layer].Vec().Replace([&]() { return (2 * rnd.drand48() - 1) * max; });
    }
  }

	// base class of functions to be optimized with SGD
	// stores parameters and input and output layers
  class TrainingFunction {
  protected:
		// DBN whose parameters we are training
    DBN & _dbn;

    BatchSource *batchSrc;

		// pointers to parameters
    MutableVector _params, _deltaParams, _savedParams;

		// temporary input and output matrices for making edge-case minibatches
    AllocatingMatrix _tempInput, _tempOutput;

		// training hyperparameters
    const HyperParams & _hyperParams;


    TrainingFunction(DBN & dbn, BatchSource *batchSrc, const HyperParams & hyperParams, MutableVector & params, MutableVector & deltaParams, MutableVector & savedParams)
      :
    _dbn(dbn), batchSrc(batchSrc), _params(params), _deltaParams(deltaParams), _savedParams(savedParams), _hyperParams(hyperParams) 
    { }

		// to be overridden in subclasses, implementing the specific training function update
    virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) = 0;

  public:

    // evaluate the function on the entire dataset and return the value (performing no update)
    virtual double Eval(float epochFraction=1.0) {
      double sum = 0;
      unsigned numRows = batchSrc->numDataRows() + batchSrc->numLabelRows();
      unsigned miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / numRows;
      if (miniBatchSize < 1) {
	miniBatchSize = 1;
	warning("WARNING: -nnChunkSize %u MiB is too small. One training instances requires %u bytes\n",
		nnChunkSize, numRows * sizeof(double));
      }
      unsigned numBatches = (unsigned)(0.5 + batchSrc->batchesPerEpoch(miniBatchSize) * epochFraction);
      unsigned data_NumC = 0;
      for (unsigned i=0; i < numBatches; i+=1) {
	Matrix  inputMiniBatch, outputMiniBatch;
	batchSrc->getBatch(miniBatchSize, inputMiniBatch, outputMiniBatch);
	data_NumC += outputMiniBatch.NumC();
	sum += PrivateEval(inputMiniBatch, outputMiniBatch, 0);
      }
      return sum / data_NumC + Decay(_params, _hyperParams.l2, 0);
    }

    double DoUpdate(Matrix &inputMiniBatch, Matrix &outputMiniBatch, double stepSize, double maxUpdate, double momentum) {
      int miniBatchSize = outputMiniBatch.NumC();
      _deltaParams *= momentum;
      double val = PrivateEval(inputMiniBatch, outputMiniBatch, stepSize * (1 - momentum) / miniBatchSize);
      if (!IsNaN(maxUpdate)) Trunc(_deltaParams, maxUpdate);
      _params -= _deltaParams;

      return val / miniBatchSize + Decay(_params, _hyperParams.l2, stepSize);
    }

		// update the parameters using the ith minibatch
		// with the specified step size, max update and momentum
		// calls PrivateEval for actual evaluation and gradient computation
    virtual double Update(double stepSize, double maxUpdate, double momentum) {
      Matrix  inputMiniBatch, outputMiniBatch;
      batchSrc->getBatch(_hyperParams.miniBatchSize, inputMiniBatch, outputMiniBatch);
      return DoUpdate(inputMiniBatch, outputMiniBatch, stepSize, maxUpdate, momentum);
    }

		// save parameters
    void Save() {
      _savedParams.CopyFrom(_params);
    }

		// load saved parameters
    void Restore() {
      _params.CopyFrom(_savedParams);
    }

		// reset saved and delta parameters
    virtual void Init() {
      _savedParams.CopyFrom(_params);
      _deltaParams *= 0;
    }
  };

	// base class for training functions that operate on a single layer's parameters
  class LayerTrainingFunction : public TrainingFunction {
  protected:
    int _layer;
    AllocatingVector _inputBiases;

  public:

    LayerTrainingFunction(DBN & dbn, int layer, BatchSource *batchSrc, const HyperParams & hyperParams, Layer::ActFunc actFunc, float epochFraction=1.0)
      :
    TrainingFunction(dbn, batchSrc, hyperParams, dbn._layerParams[layer], dbn._layerDeltaParams[layer], dbn._layerSavedParams[layer]),
      _layer(layer)
    {
      if (resumeTraining) {
        _inputBiases.CopyFrom(dbn._B[layer]);
      } else {
	DBN::InitializeInputBiases(batchSrc, _inputBiases, actFunc, 1e-3, epochFraction);
      }
    }

    // Evaluate the function on the entire dataset and return the value (performing no update)
    // For pretraining TrainingFunctions, use the input as the target output.
    virtual double Eval(float epochFraction=1.0) {
      double sum = 0;
      unsigned numRows = batchSrc->numDataRows();
      unsigned miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / numRows;
      if (miniBatchSize < 1) {
	miniBatchSize = 1;
	warning("WARNING: -nnChunkSize %u MiB is too small. One training instances requires %u bytes\n",
		nnChunkSize, numRows * sizeof(double));
      }
      unsigned numBatches = (unsigned)(0.5 + batchSrc->batchesPerEpoch(miniBatchSize) * epochFraction);
      unsigned data_NumC = 0;
      for (unsigned i=0; i < numBatches; i+=1) {
	Matrix  inputMiniBatch = batchSrc->getData(miniBatchSize);
	data_NumC += inputMiniBatch.NumC();
	sum += PrivateEval(inputMiniBatch, inputMiniBatch, 0);
      }
      return sum / data_NumC + Decay(_params, _hyperParams.l2, 0);
    }

    // Update the parameters using the ith minibatch with the specified step size, max update and momentum.
    // Calls PrivateEval for actual evaluation and gradient computation.
    // For pretraining TrainingFunctions, use input as the target output
    virtual double Update(double stepSize, double maxUpdate, double momentum) {
      Matrix  inputMiniBatch =  batchSrc->getData(_hyperParams.miniBatchSize);
      return DoUpdate(inputMiniBatch, inputMiniBatch, stepSize, maxUpdate, momentum);
    }
  };



	// training function to implement the contrastive divergence update
  class CDTrainingFunction : public LayerTrainingFunction {
  protected:
    AllocatingMatrix _particles;
		
    virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
      return _dbn.UpdateCD(outputMiniBatch, _inputBiases, _layer, stepSize);
    }

  public:
  CDTrainingFunction(DBN & dbn, int layer, BatchSource *batchSrc, const HyperParams & hyperParams, Layer::ActFunc actFunc, float epochFraction=1.0)
    : LayerTrainingFunction(dbn, layer, batchSrc, hyperParams, actFunc, epochFraction)
    { }
  };

	// training function to implement denoising autoencoder update
  class AETrainingFunction : public LayerTrainingFunction {

    Layer::ActFunc lowerActFunc;
      
  protected:
    virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
      return _dbn.UpdateAE(inputMiniBatch, _inputBiases, outputMiniBatch, _layer, stepSize);
    }

  public:

  AETrainingFunction(DBN & dbn, int layer, BatchSource *batchSrc, const HyperParams & hyperParams, Layer::ActFunc lowerActFunc, float epochFraction=1.0)
    : LayerTrainingFunction(dbn, layer, batchSrc, hyperParams, lowerActFunc, epochFraction), lowerActFunc(lowerActFunc)
    { }

    // evaluate the function on the entire (distorted) dataset and return the value (performing no update)
    virtual double Eval(float epochFraction=1.0) {
      double sum = 0;
      unsigned numRows = batchSrc->numDataRows();
      unsigned miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / numRows;
      if (miniBatchSize < 1) {
	miniBatchSize = 1;
	warning("WARNING: -nnChunkSize %u MiB is too small. One training instances requires %u bytes\n",
		nnChunkSize, numRows * sizeof(double));
      }
      unsigned numBatches = (unsigned)(0.5 + batchSrc->batchesPerEpoch(miniBatchSize) * epochFraction);
      unsigned data_NumC = 0;
      for (unsigned i=0; i < numBatches; i+=1) {
	Matrix  inputMiniBatch = batchSrc->getData(miniBatchSize);
	AllocatingMatrix temp(inputMiniBatch.NumR(), inputMiniBatch.NumC());
	// distort the input, use undistorted input as the target output
	lowerActFunc.Sample(inputMiniBatch.Vec(), temp.Vec());
	data_NumC += inputMiniBatch.NumC();
	sum += PrivateEval(temp, inputMiniBatch, 0);
      }
      return sum / data_NumC + Decay(_params, _hyperParams.l2, 0);
    }


    // Update on distorted input
    virtual double Update(double stepSize, double maxUpdate, double momentum) {
      Matrix  inputMiniBatch = batchSrc->getData(_hyperParams.miniBatchSize);
      AllocatingMatrix temp(inputMiniBatch.NumR(), inputMiniBatch.NumC());
      lowerActFunc.Sample(inputMiniBatch.Vec(), temp.Vec());
      return DoUpdate(temp, inputMiniBatch, stepSize, maxUpdate, momentum);
    }
  };

	// training function for backpropagation training
  class BPTrainingFunction : public TrainingFunction {

   private:
		// training objective (squared error or softmax)
    ObjectiveType _objectiveType;

   protected:
    virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
      return _dbn.UpdateBackProp(inputMiniBatch, outputMiniBatch, _objectiveType, false, _hyperParams.iDropP, _hyperParams.hDropP, stepSize);
    }

   public:

  BPTrainingFunction(BatchSource *batchSrc, DBN & dbn, HyperParams hyperParams, ObjectiveType objectiveType)
    : TrainingFunction(dbn, batchSrc, hyperParams, dbn._params, dbn._deltaParams, dbn._savedParams), _objectiveType(objectiveType)
    { }
  };
	
	// training function for output layer before full backpropagation
  class OutputLayerTrainingFunction : public TrainingFunction {

   private:
		// training objective (squared error or softmax)
    ObjectiveType _objectiveType;

   protected:
    virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
      return _dbn.UpdateBackProp(inputMiniBatch, outputMiniBatch, _objectiveType, true, _hyperParams.iDropP, _hyperParams.hDropP, stepSize);
    }

   public:

    OutputLayerTrainingFunction(BatchSource *batchSrc, DBN & dbn, HyperParams hyperParams, ObjectiveType objectiveType)
      :
    TrainingFunction(dbn, batchSrc, hyperParams, dbn._layerParams.back(), dbn._layerDeltaParams.back(), dbn._layerSavedParams.back()),  _objectiveType(objectiveType)
    { }
  };

  // Implments stochastic gradient descent training with momentum.
  // This version always does the full # of requested [anneal] updates and
  // assumes the training state won't be saved or restored. This version is
  // used for pre-training.
  void TrainSGD(TrainingFunction & trainer, const HyperParams & hyperParams) {
    float epochFraction = 1.0, annealEpochFraction = 1.0;
    double prevStepSize = 0.0;
    int curUpdate = 0, curAnnealUpdate = 0;
    assert(!resumeTraining); // this version only supports pre-training, which cannot be resumed
    TrainSGD(trainer, hyperParams, epochFraction, annealEpochFraction, prevStepSize, curUpdate, curAnnealUpdate);
  }


  // Implments stochastic gradient descent training with momentum.
  // This version supports partial training, and should only be used by back prop.
  void TrainSGD(TrainingFunction & trainer, const HyperParams & hyperParams, 
		float epochFraction,             // fraction of total non-anneal updates to do in this program invocation
		float annealEpochFraction,       // fraction of total anneal updates ...
		double &prevStepSize,            // non-anneal step size from previous program invocation
		int &curUpdate,                  // non-anneal update # from previous ...
		int &curAnnealUpdate)            // anneal update # from previous ...
  {
    double stepSize = resumeTraining ? prevStepSize : hyperParams.initStepSize;
    unsigned t = 0;

    if (!resumeTraining) {
      trainer.Init();

      // get a baseline score by evaluating on the whole dataset
      double bestTrainScore = trainer.Eval();

      infoMsg(IM::Training, IM::Default, "Baseline score: %e \n", bestTrainScore);

      // ensure that step size gives improvement over one check period before adding momentum
      for (;;) {
	double totalScore = 0;
	for (t = 0; t < (unsigned)hyperParams.checkInterval; t += 1) {
	  totalScore += trainer.Update(stepSize, hyperParams.maxUpdate, 0);
	}

	double score = totalScore / hyperParams.checkInterval;

	infoMsg(IM::Training, IM::Default, "Improvement check score: %e\n", score);

	if (score < bestTrainScore) {
	  infoMsg(IM::Training, IM::Default, "Keeping step size of %e\n", stepSize);
	  break;
	}

	// step size did not yield improvement, so it is probably too large
	stepSize /= 2;
	if (stepSize < MIN_STEP_SIZE) {
	  error("ERROR: step size < %e does not improve check score; giving up\n", MIN_STEP_SIZE);
	}
	infoMsg(IM::Training, IM::Default, "Step size reduced to %e\n", stepSize);
	trainer.Restore(); // undo failed training with previous stepSize
      }
    }
 
    double totalScore = 0;
    if (resumeTraining) t = curUpdate; // non-anneal update number to resume at
    unsigned numUpdates = (unsigned) (epochFraction * hyperParams.numUpdates + 0.5);
    unsigned lastUpdate = t + numUpdates; // last update # for this invocation
    if (lastUpdate > (unsigned)hyperParams.numUpdates) lastUpdate = hyperParams.numUpdates;  // in case numUpdates doesn't divide evenly
    infoMsg(IM::Training, IM::Moderate, "Training updates %u to %u of %d\n", t, lastUpdate, hyperParams.numUpdates);
    unsigned iterations = 0; // iterations since last check score
    for (; t <= lastUpdate; t+=1, iterations += 1) {
      // momentum increases from 0.5 to max
      // using formula from Sutskever et al. 2013
      double tFrac = (double)t / hyperParams.numUpdates;
      double m = tFrac / (1 - hyperParams.maxMomentum) + (1 - tFrac) / (1 - hyperParams.minMomentum);
      double momentum = 1.0 - 1.0 / m;

      totalScore += trainer.Update(stepSize, hyperParams.maxUpdate, momentum);

      // output online training objective value regularly
      if (checkSignal || t % hyperParams.checkInterval == 0) {
        double score;
	if (iterations > 0) {
	  score = totalScore / iterations;
	  infoMsg(IM::Training, IM::Default, "Step: %d\nMomentum: %e\nNew score: %e\n\n", t, momentum, score);
	}
        totalScore = 0;
	iterations = 0;
	checkSignal = false;
      }
    }
    curUpdate = t; // remember in case we save training state later

    double step = stepSize; // remember this in case need to save state but we don't anneal
    totalScore = 0;
    t = resumeTraining ? (unsigned)curAnnealUpdate : 0; // anneal update #
    numUpdates = (unsigned) (annealEpochFraction * hyperParams.numAnnealUpdates + 0.5);
    lastUpdate = t + numUpdates; // last anneal update # for this invocation
    if (lastUpdate > (unsigned)hyperParams.numAnnealUpdates) lastUpdate = hyperParams.numAnnealUpdates;

    if (t < lastUpdate) {
      infoMsg(IM::Training, IM::Default, "Annealing\n");
      infoMsg(IM::Training, IM::Moderate, "Training anneal updates %u to %u of %d\n", 
	      t, lastUpdate, hyperParams.numAnnealUpdates);
    }
    for (; t < lastUpdate; t+=1) {
      // anneal momentum and step size to zero
      double tFrac = (double)t / hyperParams.numAnnealUpdates;
      double m = tFrac + (1 - tFrac) / (1 - hyperParams.maxMomentum);
      double momentum = 1.0 - 1.0 / m;
      step = (1 - tFrac) * stepSize;

      totalScore += trainer.Update(step, hyperParams.maxUpdate, momentum);
    }
    curAnnealUpdate = t;
    prevStepSize = step;

    if (numUpdates > 0) {
      double score = totalScore / hyperParams.numAnnealUpdates;
      infoMsg(IM::Training, IM::Default, "Anneal score: %e\n\n", score);
    }
  }

	// perform the backpropagation update to update the weights of all layers
	// if lastLayerOnly, then backpropagation is terminated after updating the output layer
  double UpdateBackProp(const Matrix &input, const Matrix& output, ObjectiveType _objType, bool lastLayerOnly, double iDropP, double hDropP, double stepSize) {
    int startLayer = lastLayerOnly ? _numLayers - 1 : 0;

		// the input matrix (which may be altered if dropout is used
    Matrix altInput;
    if (iDropP > 0) {
      _tempDropoutInput.CopyFrom(input);
      _tempDropoutInput.Vec().Apply([&] (double x) { return (rnd.drand48() > iDropP) ? x : 0; });
      altInput = _tempDropoutInput;
    } else {
      altInput = input;
    }

		// map the data up to the output layer
    Matrix mappedInput = MapUp(altInput, hDropP, startLayer);

    double loss = 0;

		// compute the loss, depending on objective type
    // _tempMat is used to store the gradient w.r.t. the outputs at each layer
    switch (_objType) {
    case SQ_ERR:
      _tempMat.CopyFrom(mappedInput);
      _tempMat -= output;

      loss += _tempMat.Vec() * _tempMat.Vec();
      loss /= 2;
      break;

    case SOFT_MAX:
      _tempMat.CopyFrom(mappedInput);
      for (int i = 0; i < input.NumC(); ++i) {
        MutableVector errCol = _tempMat.GetCol(i);
        double max = Max(errCol);
        errCol.Apply([max] (double x) { return x - max; });
        loss += max;
      }
#if HAVE_MKL
      _tempMat.Vec().ApplyVML(vdExp);
#else
      _tempMat.Vec().Apply([](double x)->double {return exp(x);});
#endif

      for (int i = 0; i < input.NumC(); ++i) {
        MutableVector errCol = _tempMat.GetCol(i);
        Vector labelCol = output.GetCol(i);
        Vector predictionCol = mappedInput.GetCol(i);
        double sum = Sum(errCol);
        errCol /= sum;
        errCol -= labelCol;
        loss += log(sum) - predictionCol * labelCol;
      }
      break;
    }

    if (stepSize > 0) {
			// backpropagate errors
      Matrix mappedError = _layers.back().ComputeErrors(_tempMat, Layer::ActFunc::LINEAR);
      for (int l = _numLayers - 1; ; --l) {
        for (int i = 0; i < input.NumC(); ++i) _deltaB[l] += mappedError.GetCol(i) * stepSize;

				// inputs to this layer
				const Matrix & lowerProbs = l > startLayer ? _layers[l - 1].Activations() : altInput;

				// perform actual update
        _deltaW[l] += stepSize * lowerProbs * mappedError.Trans();

        if (lastLayerOnly || l == 0) break;

				// backpropagate error
        _tempMat = _W[l] * mappedError;
        mappedError = _layers[l-1].ComputeErrors(_tempMat, _hActFunc[l-1]);
      }
    }

    return loss;
  }

	// determine biases used for input during pretraining (either autoencoder or contrastive divergence)
	// essentially minimizes the loss for a zero weight matrix
  static void 
  InitializeInputBiases(BatchSource *batchSrc, AllocatingVector & inputBiases, Layer::ActFunc actFunc, double eps, float epochFraction=1.0) 
  {
    unsigned n = 0;
    unsigned numRows = batchSrc->numDataRows();
    inputBiases.Assign(numRows, 0);
    unsigned miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / numRows;
    if (miniBatchSize < 1) {
      miniBatchSize = 1;
      warning("WARNING: -nnChunkSize %u MiB is too small. One training instances requires %u bytes\n",
	      nnChunkSize, numRows * sizeof(double));
    }
    unsigned numBatches = (unsigned)(0.5 + batchSrc->batchesPerEpoch(miniBatchSize) * epochFraction);
    for (unsigned t=0; t < numBatches; t+=1) {
      Matrix miniBatch = batchSrc->getData(miniBatchSize);
      unsigned numCols = miniBatch.NumC();
      n += numCols;
      int miniBatchLD = miniBatch.Ld();
      double const *miniBatchP = miniBatch.Start();
      for (unsigned c=0; c < numCols; c+=1, miniBatchP += miniBatchLD) {
	for (unsigned r=0; r < numRows; r+=1) {
	  inputBiases[r] += miniBatchP[r];
	}
      }
    }
    for (unsigned r=0; r < numRows; r+=1) {
      inputBiases[r] /= n;
    }
		
    for (unsigned r=0; r < numRows; r+=1) {
      double x = inputBiases[r];
      switch (actFunc.actType) {
      case Layer::ActFunc::LOG_SIG:
	x = eps + (1 - 2 * eps) * x;
	if (x <= 0.0 || 1.0 <= x) {
	  error("ERROR: input %e is out of the allowed range (0,1) for the sigmoid input activation function\n", x);
	}
	x = log(x / (1-x));
	break;
				
      case Layer::ActFunc::TANH:
	x = (1 - eps) * x;
	if (x <= -1.0 || 1.0 <= x) {
	  error("ERROR: input %e is out of the allowed range (-1,1) for the hyperbolic tangent input activation function\n", x);
	}
	x = atanh(x);
	break;
				
      case Layer::ActFunc::CUBIC:
	x = x * x * x / 3 + x;
	break;
				
      default:
	// LINEAR and RECT_LIN require no transformation
	break;
      }

      inputBiases[r] = x;
    }
  }


	// update the layer's parameters using contrastive divergence
  double UpdateCD(const Matrix & input, const Vector & inputBiases, int layer, double stepSize) {
    Layer::ActFunc lowerActFunc = (layer == 0) ? _iActFunc : _hActFunc[layer-1];

    Matrix weights = _W[layer];
    Vector topBiases = _B[layer];
    Layer & topLayer = _layers[layer];

    _tempBottomSample.Resize(input);

    topLayer.ActivateUp(weights, topBiases, input, _hActFunc[layer]);
    _tempTopSample.Resize(topLayer.Size(), input.NumC());
    topLayer.Sample(_tempTopSample, _hActFunc[layer]);

    _tempBottomLayer.ActivateDown(weights, inputBiases, _tempTopSample, lowerActFunc);
    _tempBottomLayer.Sample(_tempBottomSample, lowerActFunc);
    _tempTopLayer.ActivateUp(weights, topBiases, _tempBottomSample, _hActFunc[layer]);

    if (stepSize > 0) {
			// actually perform update
      Matrix topProbsP = topLayer.Activations(), topProbsN = _tempTopLayer.Activations();

      _deltaW[layer] -= stepSize * input * topProbsP.Trans();
      _deltaW[layer] += stepSize * _tempBottomSample * topProbsN.Trans();

      for (int i = 0; i < input.NumC(); ++i) {
        _deltaB[layer] -= stepSize * topProbsP.GetCol(i);
        _deltaB[layer] += stepSize * topProbsN.GetCol(i);
      }
    }

    // often squared error is used for monitoring, like this

    _tempBottomSample = input - _tempBottomLayer.Activations();
    return _tempBottomSample.Vec() * _tempBottomSample.Vec();
  }

	// update layer's parameters using the autoencoder loss
  double UpdateAE(const Matrix & input, const Vector & inputBiases, const Matrix & targetInput, int layer, double stepSize) {
    Layer::ActFunc lowerActFunc = (layer == 0) ? _iActFunc : _hActFunc[layer-1];

    Matrix weights = _W[layer];
    Vector topBiases = _B[layer];

    Layer & topLayer = _layers[layer];
    Matrix topProbs = topLayer.ActivateUp(weights, topBiases, input, _hActFunc[layer]);

    double negll = _tempBottomLayer.ActivateDownAndGetNegLL(weights, inputBiases, topProbs, targetInput, lowerActFunc);
    if (stepSize > 0) {
      // fill bottomProbs with errors = activations - input
			// it is independent of activation type
      MutableMatrix & bottomProbs = _tempBottomLayer.Activations();
      bottomProbs -= targetInput;

      _deltaW[layer] += stepSize * bottomProbs * topProbs.Trans();

      _tempMat = weights.Trans() * bottomProbs;
      Matrix topError = topLayer.ComputeErrors(_tempMat, _hActFunc[layer]);

      _deltaW[layer] += stepSize * input * topError.Trans();

      for (int i = 0; i < input.NumC(); ++i) {
        _deltaB[layer] += stepSize * topError.GetCol(i);
      }
    }

    return negll;
  }

public:
  DBN(int numLayers, int iSize, vector<int> &hSize, int oSize, Layer::ActFunc iActFunc, vector<Layer::ActFunc> &hActFunc) 
  {
    Initialize(numLayers, iSize, hSize, oSize, iActFunc, hActFunc);
  }

  // resume training with previously learned W and B
  DBN(int numLayers, int iSize, vector<int> &hSize, int oSize, Layer::ActFunc iActFunc, vector<Layer::ActFunc> &hActFunc, vector<AllocatingMatrix> &W, vector<AllocatingVector> &B)
  {
    Initialize(numLayers, iSize, hSize, oSize, iActFunc, hActFunc);
    for (int i=0; i < numLayers; i+=1) {
      _W[i].CopyFrom(W[i]);
      _B[i].CopyFrom(B[i]);
    }
  }

	// allocate memory for parameters and activations at all layers
  void Initialize(int numLayers, int iSize, vector<int> &hSize, int oSize, Layer::ActFunc iActFunc, vector<Layer::ActFunc> &hActFunc)
  {
    _numLayers = numLayers;
    _iSize = iSize;
    _hSize = hSize;
    _oSize = oSize;

    _iActFunc = iActFunc;
    _hActFunc = hActFunc;
    assert(_hActFunc[_numLayers-1].actType == Layer::ActFunc::LINEAR);

    _layers.resize(numLayers);
    _layerParams.resize(numLayers);
    _layerDeltaParams.resize(numLayers);
    _layerSavedParams.resize(numLayers);

    int numParams = NumParams();
    _params.Resize(numParams);
    _deltaParams.Resize(numParams);
    _savedParams.Resize(numParams);

    _W.resize(numLayers);
    _B.resize(numLayers);
    _deltaW.resize(numLayers);
    _deltaB.resize(numLayers);
    _savedW.resize(numLayers);
    _savedB.resize(numLayers);

    int bSize = _iSize;
    int start = 0;
    for (int i = 0; i < _numLayers; ++i) {
      int tSize = (i == _numLayers - 1) ? _oSize : _hSize[i];
      int endW = start + bSize * tSize;
      _W[i] = _params.SubVector(start, endW).AsMatrix(bSize, tSize);
      _deltaW[i] = _deltaParams.SubVector(start, endW).AsMatrix(bSize, tSize);
      _savedW[i] = _savedParams.SubVector(start, endW).AsMatrix(bSize, tSize);

      int endB = endW + tSize;
      _B[i] = _params.SubVector(endW, endB);
      _deltaB[i] = _deltaParams.SubVector(endW, endB);
      _savedB[i] = _savedParams.SubVector(endW, endB);

      _layerParams[i] = _params.SubVector(start, endB);
      _layerDeltaParams[i] = _deltaParams.SubVector(start, endB);
      _layerSavedParams[i] = _savedParams.SubVector(start, endB);

      start = endB;
      bSize = tSize;
    }

    assert (start == numParams);
  }

  Matrix const &getWeights(int layer) {
    assert(0 <= layer && layer < _numLayers);
    return _W[layer];
  }

  Vector const &getBias(int layer) {
    assert(0 <= layer && layer < _numLayers);
    return _B[layer];
  }

	// read all parameters from disk
  void Deserialize(istream & inStream) {
    inStream.read((char *) &_numLayers, sizeof(int));
    inStream.read((char *) &_iSize, sizeof(int));
    inStream.read((char *) &_hSize[0], sizeof(int) * _numLayers);
    inStream.read((char *) &_oSize, sizeof(int));
    inStream.read((char *) &_iActFunc, sizeof(Layer::ActFunc));
    inStream.read((char *) &_hActFunc[0], sizeof(Layer::ActFunc) * _numLayers);

    Initialize(_numLayers, _iSize, _hSize, _oSize, _iActFunc, _hActFunc);

    inStream.read((char *) _params.Start(), sizeof(double) * _params.Len());
  }

	// write all parameters to disk
  void Serialize(ostream & outStream) const {
    outStream.write((const char *) &_numLayers, sizeof(int));
    outStream.write((const char *) &_iSize, sizeof(int));
    outStream.write((const char *) &_hSize[0], sizeof(int) * _numLayers);
    outStream.write((const char *) &_oSize, sizeof(int));
    outStream.write((const char *) &_iActFunc, sizeof(Layer::ActFunc));
    outStream.write((const char *) &_hActFunc[0], sizeof(Layer::ActFunc) * _numLayers);
    outStream.write((const char *) _params.Start(), sizeof(double) * _params.Len());
  }

	// randomize parameters
  void Randomize(double stdDev) {
    _params.Apply([&] (double x) { return rnd.normal() * stdDev; });
  }

  int NumLayers() const { return _numLayers; }

	// total number of parameters
  int NumParams() const {
    int numParams = 0;
    for (int layer = 0; layer < _numLayers; ++layer) numParams += NumParamsInLayer(layer);

    return numParams;
  }

	// number of parameters in a given layer
  int NumParamsInLayer(int layer) const {
    assert (0 <= layer && layer < _numLayers);
    if (_numLayers == 1) return (_iSize + 1) * _oSize;
    else if (layer == 0) return (_iSize + 1) * _hSize[layer];
    else if (layer == _numLayers - 1) return (_hSize[layer-1] + 1) * _oSize;
    else return (_hSize[layer-1] + 1) * _hSize[layer];
  }

	// output size of a layer
  int LayerOutSize(int layer) const {
    assert (0 <= layer && layer < _numLayers);
    return (layer < _numLayers - 1) ? _hSize[layer] : _oSize;
  }

	// expected input size of a layer
  int LayerInSize(int layer) const {
    assert (0 <= layer && layer < _numLayers);
    return (layer == 0) ? _iSize : _hSize[layer-1];
  }

	// given an input matrix to a given layer, compute the output matrix
	// by multiplying by weights and adding biases
  StdioMatrix MapLayer(StdioMatrix &input, int layer) const {
    int miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / input.Ld();
    if (miniBatchSize < 1) {
      miniBatchSize = 1;
      warning("WARNING: -nnChunkSize %u MiB is too small. One training instances requires %u bytes\n",
	      nnChunkSize, input.Ld() * sizeof(double));
    }
    StdioMatrix output(_B[layer].Len(), input.NumC(), _B[layer].Len());
    int start=0, end, lastCol = input.NumC();
    for ( ; start < lastCol ; start += miniBatchSize) {
      end = (start + miniBatchSize <= lastCol) ? start + miniBatchSize : lastCol;
      Matrix const &m = input.GetCols(start, end);
      Matrix const &activated = _layers[layer].ActivateUp(_W[layer], _B[layer], m, _hActFunc[layer]);
      output.PutColsM(activated, (int)start);
    }
    return output;
  }

	// given input data to a given layer from observation files, compute the output matrix
	// by multiplying by weights and adding biases
  StdioMatrix MapLayer(BatchSource *batchSrc, int layer, StdioMatrix &batchSrcLabels, unsigned numCols) const {
    unsigned numRows = batchSrc->numDataRows() + batchSrc->numLabelRows();
    unsigned miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / numRows;
    if (miniBatchSize < 1) {
      miniBatchSize = 1;
      warning("WARNING: -nnChunkSize %u MiB is too small. One training instances requires %u bytes\n",
	      nnChunkSize, numRows * sizeof(double));
    }
    unsigned start, end, lastCol = min(numCols, batchSrc->epochSize());
    StdioMatrix output(_B[layer].Len(), lastCol, _B[layer].Len());
    //    assert(batchSrcLabels.NumC() == (int)lastCol);
    for (start=0, end=miniBatchSize; start < lastCol ; start += miniBatchSize, end += miniBatchSize) {
      unsigned batchSize;
      if (end <= lastCol) {
	// take the full minibatch
	batchSize = miniBatchSize;
      } else {
	end = lastCol; // partial batch
	batchSize = lastCol - start;
      }
      //      end = (start + miniBatchSize <= lastCol) ? start + miniBatchSize : lastCol;
      Matrix d, l; // instance data, labels
      batchSrc->getBatch(batchSize, d, l);
assert(d.NumC() == l.NumC());
      output.PutColsM( _layers[layer].ActivateUp(_W[layer], _B[layer], d, _hActFunc[layer]) , start);
      batchSrcLabels.PutColsM(l, start);
    }
    return output;
  }


	// given an input matrix to a given layer, compute the output matrix
	// by multiplying by weights and adding biases
  Matrix MapLayer(const Matrix & input, int layer) const {
    return _layers[layer].ActivateUp(_W[layer], _B[layer], input, _hActFunc[layer]);
  }


	// given an input matrix to a given layer, compute the output matrix of a higher layer
	// by consecutively multiplying by weights and adding biases
  Matrix MapUp(const Matrix &input, double dropoutRate = 0, int startLayer = -1, int endLayer = -1) const {
    if (startLayer == -1) startLayer = 0;
    if (endLayer == -1) endLayer = _numLayers;

    Matrix mappedInput = input;

    for (int layer = startLayer; layer < endLayer; ++layer) {
      mappedInput = MapLayer(mappedInput, layer);
      if (layer + 1 < _numLayers && dropoutRate > 0) {
        mappedInput = _layers[layer].Dropout(dropoutRate);
      }
    }
    return mappedInput;
  }

	// (Pretrain and) train parameters of entire network
  void Train(BatchSource *batchSrc, ObjectiveType objectiveType, 
	     const vector<HyperParams> &hyperParams_pt, const HyperParams &hyperParams_bp, 
	     float ptFraction = 1.0, float epochFraction = 1.0, float annealEpochFraction = 1.0,
	     char const *loadFilename=NULL, char const *saveFilename=NULL) 
  {
    if (IM::messageGlb(IM::Training, IM::Moderate)) {
      printf("\nPretrain hyperparams:\n\n");
      hyperParams_pt[0].print();
      printf("\nBackprop hyperparams:\n\n");
      hyperParams_bp.print();
      printf("\n");
    }

    // numCols is the # of columns in the temp file holding the mapped up input to the next layer
    unsigned numCols = (ptFraction >= 1.0) ? batchSrc->epochSize() : (ptFraction * batchSrc->epochSize());

    StdioMatrix   input,  // current layer's input data
                 output, // current layer's output data (will become next layer's input)
      // labels for the MapLayer()ed input instances
      batchSrcLabels(batchSrc->numLabelRows(), numCols, batchSrc->numLabelRows());

    BatchSource *bs;


    // unnormalize weights if dropout training was used
    if (resumeTraining) {
      if (hyperParams_bp.iDropP > 0) _W[0] *= 1.0 / (1 - hyperParams_bp.iDropP);
      if (hyperParams_bp.hDropP > 0) {
	for (unsigned l = 1; l < _W.size(); ++l) {
	  _W[l] *= 1.0 / (1 - hyperParams_bp.hDropP);
	}
      }
    }

    if (!resumeTraining) { // we only support resuming backprop, not pretraining

      // pretrain hidden layers
      unsigned layer;
      for (layer = 0; layer < _layers.size() - 1; ++layer) {
	const HyperParams & hyperParams = hyperParams_pt[layer];
	if (randomInitLayer) InitLayer(layer);
	
	Layer::ActFunc lowerActFunc = (layer == 0) ? _iActFunc : _hActFunc[layer-1];
	if (layer == 0) {
	  bs = batchSrc; // use observation file data
	} else {
	  bs = new MatrixBatchSource(input); // use previous layer's output
	}
	// For the first layer, we only want to use ptFraction of an epoch as
	// input to InitializeInputBiases() (and there's no point using more
	// than one epoch for it, hence the min()).

	// For the remaining layers, we already have the truncated data size,
	// so we want to use all available input.
	float curPtFraction = (layer > 0) ? 1.0 : min(1.0f, ptFraction);
	switch (hyperParams.pretrainType) {
	case CD:
	  {
	    infoMsg(IM::Training, IM::Default, "Pretraining layer %d of size %d x %d with CD\n", 
		    layer, LayerInSize(layer), LayerOutSize(layer));
	    CDTrainingFunction cdFunc(*this, layer, bs, hyperParams, lowerActFunc, curPtFraction);
	    TrainSGD(cdFunc, hyperParams);
	  }
	  break;
	  
	case AE:
	  {
	    infoMsg(IM::Training, IM::Default, "Pretraining layer %d of size %d x %d with AE\n", 
		    layer, LayerInSize(layer), LayerOutSize(layer));
	    AETrainingFunction aeFunc(*this, layer, bs, hyperParams, lowerActFunc, curPtFraction);
	    TrainSGD(aeFunc, hyperParams);
	  }
	  break;
	  
	case NONE:
	  break;
	  
	default: assert(false); // don't specify illegal pretraining types...
	}
	if (layer == 0) {
	  // apply input layer to observation file data
	  output = MapLayer(batchSrc, layer, batchSrcLabels, numCols);
	} else {
	  output = MapLayer(input, layer);    // apply hidden layer to previous layer's output
	}
	input = output; // current layer's output becomes next layer's input
	//	assert(input.NumC() == batchSrcLabels.NumC());  no longer valid with ptNumEpoch < 1
	if (layer > 0) _layers[layer - 1].Clear(); // save some memory
      }

      if (randomInitLayer) InitLayer(_layers.size() - 1);
      
      if (hyperParams_pt.back().pretrainType != NONE) {
	// pretrain output layer
	infoMsg(IM::Training, IM::Default, "Pretraining output layer\n");

	if (layer == 0) {
	  bs = batchSrc; // use observation file data
	} else {
	  bs = new MatrixBatchSource(input, batchSrcLabels); // use previous layer's output
	}
	OutputLayerTrainingFunction outputFunc(bs, *this, hyperParams_pt.back(), objectiveType);
	TrainSGD(outputFunc, hyperParams_pt.back());
      }
    }

    // and now, fine tune all layers with backpropagation
    
    double prevStepSize = hyperParams_bp.initStepSize;
    int curUpdate = 0, curAnnealUpdate = 0;
    iDataStreamFile *loadFile = NULL;
    oDataStreamFile *saveFile = NULL;
    if (resumeTraining) {
      loadFile = new iDataStreamFile(loadFilename, true);
      loadFile->read(prevStepSize, "previous backprop step size");
      loadFile->read(curUpdate, "previous backprop update number");
      loadFile->read(curAnnealUpdate, "previous backprop anneal update number");
      // read momentum & update schedule -- squak if inconsistant w/ args
      loadFile->read(_params.Start(), _params.Len(), "previous network parameters");
      loadFile->read(_deltaParams.Start(), _deltaParams.Len(), "previous network momentum");
      delete loadFile;
    }
    
    if (epochFraction > 0.0 || annealEpochFraction > 0.0) {
      infoMsg(IM::Training, IM::Default, "Fine tuning\n");
      BPTrainingFunction bpFunc(batchSrc, *this, hyperParams_bp, objectiveType);
      TrainSGD(bpFunc, hyperParams_bp, epochFraction, annealEpochFraction, prevStepSize, curUpdate, curAnnealUpdate);
    }

    // renormalize weights if dropout training was used
    if (hyperParams_bp.iDropP > 0) _W[0] *= (1 - hyperParams_bp.iDropP);
    if (hyperParams_bp.hDropP > 0) {
      for (unsigned l = 1; l < _W.size(); ++l) {
	_W[l] *= (1 - hyperParams_bp.hDropP);
      }
    }

    if (saveFilename) {
      saveFile = new oDataStreamFile(saveFilename, true);
      saveFile->write(prevStepSize, "previous backprop step size");
      saveFile->write(curUpdate, "previous backprop update number");
      saveFile->write(curAnnealUpdate, "previous backprop anneal update number");
      // write momentum & update schedule
      saveFile->write(_params.Start(), _params.Len(), "previous network parameters");
      saveFile->write(_deltaParams.Start(), _deltaParams.Len(), "previous network momentum");
      delete saveFile;
    }
  }
};
