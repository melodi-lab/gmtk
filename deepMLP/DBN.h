#pragma once

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

#include "MMapMatrix.h"
#include "rand.h"
#include "Layer.h"
#include "MatrixFunc.h"


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
  static bool     sparseInitLayer; // select layer initialization method
  static unsigned nnChunkSize;     // O(1) space to use for (non-minibatch) incremental processing

  enum PretrainType {
    NONE, // no pretraining
    AE, // denoising autoencoder
    CD, // contrastive divergence
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
      minMomentum(0.5),
      maxMomentum(0.99),
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
#if 0
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
#endif
  }

	// base class of functions to be optimized with SGD
	// stores parameters and input and output layers
  class TrainingFunction {
  protected:
		// DBN whose parameters we are training
    DBN & _dbn;

		// input and output data
    MMapMatrix _input, _output;

		// pointers to parameters
    MutableVector _params, _deltaParams, _savedParams;

		// temporary input and output matrices for making edge-case minibatches
    AllocatingMatrix _tempInput, _tempOutput;

		// training hyperparameters
    const HyperParams & _hyperParams;

    TrainingFunction(DBN & dbn, const MMapMatrix & input, const MMapMatrix & output, const HyperParams & hyperParams, MutableVector & params, MutableVector & deltaParams, MutableVector & savedParams)
      :
    _dbn(dbn), _hyperParams(hyperParams), _input(input), _output(output), _params(params), _deltaParams(deltaParams), _savedParams(savedParams) { }

    TrainingFunction(DBN & dbn, const MMapMatrix & output, const HyperParams & hyperParams, MutableVector & params, MutableVector & deltaParams, MutableVector & savedParams)
      :
    _dbn(dbn), _hyperParams(hyperParams), _input(), _output(output), _params(params), _deltaParams(deltaParams), _savedParams(savedParams) { }

		// to be overridden in subclasses, implementing the specific training function update
    virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) = 0;

  public:
		// evaluate the function on the entire dataset and return the value (performing no update)
    double Eval() {
#if 0
      return PrivateEval(_input, _output, 0) / _output.NumC() + Decay(_params, _hyperParams.l2, 0);
#else
      double sum = 0;
      int data_NumC = _output.NumC();
      int miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / _input.Ld();
      if (miniBatchSize < 1) miniBatchSize = 1;
      for (int startCol=0; startCol < data_NumC; startCol += miniBatchSize) {
	int endCol = (startCol + miniBatchSize <= data_NumC) ? startCol + miniBatchSize : data_NumC;
	Matrix  inputMiniBatch = (_input.NumC() > 0) ? _input.GetCols(startCol, endCol) : _input;
        Matrix outputMiniBatch = _output.GetCols(startCol, endCol);
	sum += PrivateEval(inputMiniBatch, outputMiniBatch, 0);
      }
      return sum / _output.NumC() + Decay(_params, _hyperParams.l2, 0);
#endif
    }

		// update the parameters using a minibatch beginning from instance start
		// with the specified step size, max update and momentum
		// calls PrivateEval for actual evaluation and gradient computation
    double Update(int start, double stepSize, double maxUpdate, double momentum) {
      int numIns = _output.NumC();
      start %= numIns;

      int miniBatchSize = _hyperParams.miniBatchSize;
      if (miniBatchSize <= 0) miniBatchSize = numIns;

      Matrix inputMiniBatch, outputMiniBatch;
      int end = start + miniBatchSize;
      if (end <= numIns) {
        inputMiniBatch = (_input.NumC() > 0) ? _input.GetCols(start, end) : _input;
        outputMiniBatch = _output.GetCols(start, end);
      } else {
        if (_input.NumC() > 0) {
          _tempInput.Resize(_input.NumR(), miniBatchSize);
          _tempInput.GetCols(0, numIns - start).CopyFrom(_input.GetCols(start, -1));
          _tempInput.GetCols(numIns - start, -1).CopyFrom(_input.GetCols(0, end - numIns));
          inputMiniBatch = _tempInput;
        } else inputMiniBatch = _input;

        _tempOutput.Resize(_output.NumR(), miniBatchSize);
        _tempOutput.GetCols(0, numIns - start).CopyFrom(_output.GetCols(start, -1));
        _tempOutput.GetCols(numIns - start, -1).CopyFrom(_output.GetCols(0, end - numIns));
        outputMiniBatch = _tempOutput;
      }

      _deltaParams *= momentum;
      double val = PrivateEval(inputMiniBatch, outputMiniBatch, stepSize * (1 - momentum) / miniBatchSize);
      if (!IsNaN(maxUpdate)) Trunc(_deltaParams, maxUpdate);
      _params -= _deltaParams;

      return val / miniBatchSize + Decay(_params, _hyperParams.l2, stepSize);
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
  LayerTrainingFunction(DBN & dbn, int layer, MMapMatrix & trainData, const HyperParams & hyperParams, Layer::ActFunc actFunc)
      :
    TrainingFunction(dbn, trainData, hyperParams, dbn._layerParams[layer], dbn._layerDeltaParams[layer], dbn._layerSavedParams[layer]),
      _layer(layer)
    {
      if (resumeTraining) {
        _inputBiases.CopyFrom(dbn._B[layer]);
      } else {
	DBN::InitializeInputBiases(trainData, _inputBiases, actFunc, 1e-3);
      }
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
    CDTrainingFunction(DBN & dbn, int layer, MMapMatrix & trainData, const HyperParams & hyperParams, Layer::ActFunc actFunc)
      :
    LayerTrainingFunction(dbn, layer, trainData, hyperParams, actFunc)
    {
    }
  };

	// training function to implement denoising autoencoder update
  class AETrainingFunction : public LayerTrainingFunction {
    //    AllocatingMatrix _distortedInput;

  protected:
    virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
      return _dbn.UpdateAE(inputMiniBatch, _inputBiases, outputMiniBatch, _layer, stepSize);
    }

  public:		
    AETrainingFunction(DBN & dbn, int layer, MMapMatrix & trainData, const HyperParams & hyperParams, Layer::ActFunc lowerActFunc, bool fixTrainDistortion)
      : LayerTrainingFunction(dbn, layer, trainData, hyperParams, lowerActFunc)
    {
      if (fixTrainDistortion) {
				// precompute the distorted input data on the whole dataset and save it
				// if fixTrainDistortion is false, new noisy input will be created for each update
#if 1
	MMapMatrix _distortedInput(trainData.NumR(), trainData.NumC(), trainData.NumR());
	int input_NumC = trainData.NumC();
	int miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / trainData.Ld();
	if (miniBatchSize < 1) miniBatchSize = 1;
	for (int startCol=0; startCol < input_NumC; startCol += miniBatchSize) {
	  int endCol = (startCol + miniBatchSize <= input_NumC) ? startCol + miniBatchSize : input_NumC;
	  Matrix miniBatch = trainData.GetCols(startCol, endCol);
	  AllocatingMatrix temp(miniBatch.NumR(), miniBatch.NumC());
	  lowerActFunc.Sample(miniBatch.Vec(), temp.Vec());
	  _distortedInput.PutCols(temp, startCol);
	}
	_input = _distortedInput;
#else
        _distortedInput.Resize(trainData);
        lowerActFunc.Sample(trainData.Vec(), _distortedInput.Vec());
        _input = _distortedInput;
#endif
      }
    }
  };

	// training function for backpropagation training
  class BPTrainingFunction : public TrainingFunction {
    AllocatingMatrix _errors; // I think this is not referenced

  private:
		// training objective (squared error or softmax)
    ObjectiveType _objectiveType;

  protected:
    virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
      return _dbn.UpdateBackProp(inputMiniBatch, outputMiniBatch, _objectiveType, false, _hyperParams.iDropP, _hyperParams.hDropP, stepSize);
    }

  public:
    BPTrainingFunction(const MMapMatrix & input, const MMapMatrix & output, DBN & dbn, HyperParams hyperParams, ObjectiveType objectiveType)
      :
    TrainingFunction(dbn, input, output, hyperParams, dbn._params, dbn._deltaParams, dbn._savedParams),
      _objectiveType(objectiveType)
    { }
  };
	
	// training function for output layer before full backpropagation
  class OutputLayerTrainingFunction : public TrainingFunction {
    AllocatingMatrix _errors; // I think this is not referenced

  private:
		// training objective (squared error or softmax)
    ObjectiveType _objectiveType;

  protected:
    virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
      return _dbn.UpdateBackProp(inputMiniBatch, outputMiniBatch, _objectiveType, true, _hyperParams.iDropP, _hyperParams.hDropP, stepSize);
    }

  public:
    OutputLayerTrainingFunction(const MMapMatrix & input, const MMapMatrix & output, DBN & dbn, HyperParams hyperParams, ObjectiveType objectiveType)
      :
    TrainingFunction(dbn, input, output, hyperParams, dbn._layerParams.back(), dbn._layerDeltaParams.back(), dbn._layerSavedParams.back()),
      _objectiveType(objectiveType)
    { }
  };

	// implments stochastic gradient descent training with momentum
  void TrainSGD(TrainingFunction & trainer, const HyperParams & hyperParams) {
    double stepSize = hyperParams.initStepSize;
    bool quiet = false;

    if (!quiet) cout << fixed << setprecision(4);

    trainer.Init();

    // get a baseline score by evaluating on the whole dataset
    double bestTrainScore = trainer.Eval();

    if (!quiet) cout << "Baseline score: " << bestTrainScore << endl;

    int startInstance;

    // ensure that step size gives improvement over one check period before adding momentum
    for (;;) {
      startInstance = 0;

      double totalScore = 0;
      for (int t = 0; t < hyperParams.checkInterval; ++t) {
        totalScore += trainer.Update(startInstance, stepSize, hyperParams.maxUpdate, 0);
        startInstance += hyperParams.miniBatchSize;
      }

      double score = totalScore / hyperParams.checkInterval;

      if (!quiet) cout << "Improvement check score: " << score << endl;
      if (score < bestTrainScore) {
        if (!quiet) cout << "Keeping step size of " << stepSize << endl;
        break;
      }

			// step size did not yield improvement, so it is probably too large
      stepSize /= 2;
      if (stepSize < MIN_STEP_SIZE) {
	error("ERROR: step size < %e does not improve check score; giving up\n", MIN_STEP_SIZE);
      }
      if (!quiet) cout << "Step size reduced to " << stepSize << endl;
      trainer.Restore();
    }

    double totalScore = 0;
    for (int t = 1; t <= hyperParams.numUpdates; ++t) {
      // momentum increases from 0.5 to max
      // using formula from Sutskever et al. 2013
      double tFrac = (double)t / hyperParams.numUpdates;
      double m = tFrac / (1 - hyperParams.maxMomentum) + (1 - tFrac) / (1 - hyperParams.minMomentum);
      double momentum = 1.0 - 1.0 / m;

      totalScore += trainer.Update(startInstance, stepSize, hyperParams.maxUpdate, momentum);
      startInstance += hyperParams.miniBatchSize;

			// output online training objective value regularly
      if (checkSignal || !quiet && t % hyperParams.checkInterval == 0) {
        double score = totalScore / hyperParams.checkInterval;
        cout << "Step: " << t << endl;
        cout << "Momentum: " << momentum << endl;
        cout << "New score: " << score << endl << endl;
        totalScore = 0;
	checkSignal = false;
      }
    }

    totalScore = 0;
    for (int t = 0; t < hyperParams.numAnnealUpdates; ++t) {
      // anneal momentum and step size to zero
      double tFrac = (double)t / hyperParams.numAnnealUpdates;
      double m = tFrac + (1 - tFrac) / (1 - hyperParams.maxMomentum);
      double momentum = 1.0 - 1.0 / m;
      double step = (1 - tFrac) * stepSize;

      totalScore += trainer.Update(startInstance, step, hyperParams.maxUpdate, momentum);
      startInstance += hyperParams.miniBatchSize;
    }

    if (!quiet) {
      double score = totalScore / hyperParams.numAnnealUpdates;
      cout << "Anneal score: " << score << endl << endl;
    }

  }

	// perform the backpropagation update to update the weights of all layers
	// if lastLayerOnly, then backpropagation is terminated after updating the output layer
  double UpdateBackProp(const Matrix & input, const Matrix & output, ObjectiveType _objType, bool lastLayerOnly, double iDropP, double hDropP, double stepSize) {
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
  static void InitializeInputBiases(MMapMatrix & input, AllocatingVector & inputBiases, Layer::ActFunc actFunc, double eps) {
    int n = input.NumC();
    inputBiases.Assign(input.NumR(), 0);
#if 0
    AllocatingVector transformedRow(input.NumC());
    for (int r = 0; r < input.NumR(); ++r) {
      Vector row = input.GetRow(r);
      switch (actFunc.actType) {
      case Layer::ActFunc::LOG_SIG:
	transformedRow.Replace([&] (double x) { 
	    x = eps + (1 - 2 * eps) * x;
	    return log(x / (1-x));
	  }, row);
	row = transformedRow;
	break;

      case Layer::ActFunc::TANH:
	transformedRow.Replace([&] (double x) { 
	    x = (1 - eps) * x;
	    return atanh(x);
	  }, row);
	row = transformedRow;
	break;
	
      default:
	break;
      }
      inputBiases[r] = Sum(row) / n;
    }
#else
    int miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / input.Ld();
    if (miniBatchSize < 1) miniBatchSize = 1;
    int input_NumC = input.NumC();
    for (int startCol=0; startCol < input_NumC; startCol += miniBatchSize) {
      int endCol = (startCol + miniBatchSize <= input_NumC) ? startCol + miniBatchSize : input_NumC;
      Matrix miniBatch = input.GetCols(startCol, endCol);
      int numCols = miniBatch.NumC();
      int numRows = miniBatch.NumR();
      int miniBatchLD = miniBatch.Ld();
      double const *miniBatchP = miniBatch.Start();
      for (int c=0; c < numCols; c+=1, miniBatchP += miniBatchLD) {
	for (int r=0; r < numRows; r+=1) {
	  double transformed = miniBatchP[r], x;
assert(!std::isnan(transformed));
	  switch (actFunc.actType) {
	  case Layer::ActFunc::LOG_SIG:
	    x = eps + (1 - 2 * eps) * transformed;
assert(0.0 < x && x < 1.0);
	    transformed = log(x / (1-x));
assert(!std::isnan(transformed));
	    break;
			      
	  case Layer::ActFunc::TANH:
	    x = (1 - eps) * transformed;
assert(-1.0 < x && x < 1.0);
	    transformed = atanh(x);
assert(!std::isnan(transformed));
	    break;
	    
	  default:
	    break;
	  }
assert(!std::isnan(inputBiases[r]));
	  inputBiases[r] += transformed;
	}
      }
    }
    int numRows = input.NumR();
    for (int r=0; r < numRows; r+=1) {
      assert(!std::isnan(inputBiases[r]));
      inputBiases[r] /= n;
    }
#endif
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
  MMapMatrix MapLayer(MMapMatrix &input, int layer) const {
    int miniBatchSize = nnChunkSize * ((1<<20) / sizeof(double)) / input.Ld();
    if (miniBatchSize < 1) miniBatchSize = 1;
    MMapMatrix output(_B[layer].Len(), input.NumC(), _B[layer].Len());
    int start=0, end, lastCol = input.NumC();
    for ( ; start < lastCol ; start += miniBatchSize) {
      end = (start + miniBatchSize <= lastCol) ? start + miniBatchSize : lastCol;
      Matrix const &m = input.GetCols(start, end);
      output.PutCols( _layers[layer].ActivateUp(_W[layer], _B[layer], m, _hActFunc[layer]) , start);
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
  Matrix MapUp(const Matrix & input, double dropoutRate = 0, int startLayer = -1, int endLayer = -1) const {
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
  void Train(MMapMatrix & input, const MMapMatrix & output, ObjectiveType objectiveType, bool quiet, const vector<HyperParams> & hyperParams_pt, const HyperParams & hyperParams_bp) {
    MMapMatrix mappedInput = input;
    if (!quiet) {
      printf("\nPretrain hyperparams:\n\n");
      hyperParams_pt[0].print();
      printf("\nBackprop hyperparams:\n\n");
      hyperParams_bp.print();
      printf("\n");
    }

    if (!resumeTraining) { // we only support resuming backprop, not pretraining

      // pretrain hidden layers
      for (int layer = 0; layer < _layers.size() - 1; ++layer) {
	const HyperParams & hyperParams = hyperParams_pt[layer];
	InitLayer(layer);
	
	Layer::ActFunc lowerActFunc = (layer == 0) ? _iActFunc : _hActFunc[layer-1];
	switch (hyperParams.pretrainType) {
	case CD:
	  {
	    if (!quiet) {
	      cout << "Pretraining layer " << layer << " of size " << LayerInSize(layer) << "x" << LayerOutSize(layer) << " with CD" << endl;
	    }
	    CDTrainingFunction cdFunc(*this, layer, mappedInput, hyperParams, lowerActFunc);
	    TrainSGD(cdFunc, hyperParams);
	  }
	  break;
	  
	case AE:
	  {
	    if (!quiet) {
	      cout << "Pretraining layer " << layer << " of size " << LayerInSize(layer) << "x" << LayerOutSize(layer) << " with AE" << endl;
	    }
	    AETrainingFunction aeFunc(*this, layer, mappedInput, hyperParams, lowerActFunc, true);
	    TrainSGD(aeFunc, hyperParams);
	  }
	  break;
	  
	case NONE:
	  break;
	  
	default:
	  abort();
	}
	int l = layer;
	mappedInput = MapLayer(mappedInput, l);
	if (layer > 0) _layers[layer - 1].Clear(); // save some memory
      }

      InitLayer(_layers.size() - 1);
      
      if (hyperParams_pt.back().pretrainType != NONE) {
	// pretrain output layer
	if (!quiet) {
	  cout << "Pretraining output layer" << endl;
	}

	OutputLayerTrainingFunction outputFunc(mappedInput, output, *this, hyperParams_pt.back(), objectiveType);
	TrainSGD(outputFunc, hyperParams_pt.back());
      }
    }

    // and now, fine tune all layers with backpropagation
    if (!quiet) {
      cout << "Fine tuning" << endl;
    }

    BPTrainingFunction bpFunc(input, output, *this, hyperParams_bp, objectiveType);
    TrainSGD(bpFunc, hyperParams_bp);

		// renormalize weights if dropout training was used
    if (hyperParams_bp.iDropP > 0) _W[0] *= (1 - hyperParams_bp.iDropP);
    if (hyperParams_bp.hDropP > 0) {
      for (int l = 1; l < _W.size(); ++l) {
        _W[l] *= (1 - hyperParams_bp.hDropP);
      }
    }
  }
};
