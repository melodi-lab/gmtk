#pragma once

#if defined(HAVE_MKL)
#  include "mkl.h"
#  include "mkl_lapacke.h"
#  include "mkl_spblas.h"
#  include "mkl_trans.h"
#elif defined(HAVE_BLAS)
#  ifdef __cplusplus
     extern "C" {            /* Assume C declarations for C++ */
#  endif /* __cplusplus */

#  include <cblas.h>
#  include <clapack.h>

#  ifdef __cplusplus
     }
#  endif    /* __cplusplus */
#endif

#include <vector>
#include <algorithm>
#include <boost/timer.hpp>
#include <fstream>

#include "Layer.h"
#include "MatrixFunc.h"

using namespace std;

class DBN {
public:
	enum PretrainType {
		NONE, // no pretraining
		AE, // denoising autoencoder
		CD, // contrastive divergence
	};

	enum ObjectiveType { SQ_ERR, SOFT_MAX };

	struct HyperParams {
		double initStepSize, maxMomentum, maxUpdate, l2;

		int numUpdates, numAnnealUpdates, miniBatchSize, checkInterval;

		PretrainType pretrainType;

		HyperParams() :
			initStepSize(1e-2),
			maxMomentum(0.99),
			maxUpdate(0.1),
			l2(1e-3),
			numUpdates(10000),
			numAnnealUpdates(2000),
			miniBatchSize(10),
			checkInterval(2000),
			pretrainType(CD)
			{ }
	};

private:
	// width of input, all hidden, and output layers
	int _iSize, _hSize, _oSize;
	Layer::ActFunc _iActFunc, _hActFunc, _oActFunc;

	// number of layers (not including input)
	int _numLayers;

	mutable vector<Layer> _layers;

	// all model parameters
	AllocatingVector _params, _deltaParams, _savedParams;

	vector<MutableMatrix> _W;
	vector<MutableVector> _B;
	vector<MutableMatrix> _deltaW;
	vector<MutableVector> _deltaB;
	vector<MutableMatrix> _savedW;
	vector<MutableVector> _savedB;
	vector<MutableVector> _layerParams, _layerDeltaParams, _layerSavedParams;

	// space for temporary values during CD training
	mutable AllocatingMatrix _tempTopSample, _tempBottomSample;

	// temporary storage used in several places
	mutable AllocatingMatrix _tempMat;
	mutable AllocatingVector _tempVec;
	mutable Layer _tempTopLayer, _tempBottomLayer;

	static double Decay(MutableVector & v, double alpha, double step) {
		double val = 0.5 * alpha * (v * v);
		if (step * alpha > 0) v *= (1.0 - step * alpha);
		return val;
	}

	void InitLayer(int layer, Random & rand) {
		// sparse initialization strategy from Martens, 2010
		_B[layer] *= 0;
		MutableMatrix W = _W[layer];
		W *= 0;
		for (int c = 0; c < W.NumC(); ++c) {
			MutableVector col = W.GetCol(c);
			// sampling with replacement, but it doesn't really matter
			for (int i = 0; i < 15; ++i) {
				col[rand.Rand(col.Len())] = rand.Normal() * 0.01;
			}
		}
	}

	class TrainingFunction {
	protected:
		DBN & _dbn;

		Matrix _input, _output;

		MutableVector _params, _deltaParams, _savedParams;

		AllocatingMatrix _tempInput, _tempOutput;

		const HyperParams & _hyperParams;

		TrainingFunction(DBN & dbn, const Matrix & input, const Matrix & output, const HyperParams & hyperParams, MutableVector & params, MutableVector & deltaParams, MutableVector & savedParams)
			:
		_dbn(dbn), _hyperParams(hyperParams), _input(input), _output(output), _params(params), _deltaParams(deltaParams), _savedParams(savedParams) { }

		TrainingFunction(DBN & dbn, const Matrix & output, const HyperParams & hyperParams, MutableVector & params, MutableVector & deltaParams, MutableVector & savedParams)
			:
		_dbn(dbn), _hyperParams(hyperParams), _output(output), _params(params), _deltaParams(deltaParams), _savedParams(savedParams) { }

		virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) = 0;

	public:
		double Eval() {
			return PrivateEval(_input, _output, 0) / _output.NumC() + Decay(_params, _hyperParams.l2, 0);
		}

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

		void Save() {
			_savedParams.CopyFrom(_params);
		}

		void Restore() {
			_params.CopyFrom(_savedParams);
		}

		virtual void Init(Random & rand) {
			_savedParams.CopyFrom(_params);
			_deltaParams *= 0;
		}
	};

	class LayerTrainingFunction : public TrainingFunction {
	protected:
		int _layer;
		AllocatingVector _inputBiases;

	public:
		LayerTrainingFunction(DBN & dbn, int layer, const Matrix & trainData, const HyperParams & hyperParams)
			:
		TrainingFunction(dbn, trainData, hyperParams, dbn._layerParams[layer], dbn._layerDeltaParams[layer], dbn._layerSavedParams[layer]),
			_layer(layer)
		{
			DBN::InitializeInputBiases(trainData, _inputBiases, 1e-3);
		}
	};

	class CDTrainingFunction : public LayerTrainingFunction {
	protected:
		Random & _rand;
		AllocatingMatrix _particles;

		virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
			return _dbn.UpdateCD(outputMiniBatch, _inputBiases, _layer, _rand, stepSize);
		}

	public:
		CDTrainingFunction(DBN & dbn, int layer, const Matrix & trainData, Random & rand, const HyperParams & hyperParams)
			:
		LayerTrainingFunction(dbn, layer, trainData, hyperParams), _rand(rand)
		{
		}
	};

	class AETrainingFunction : public LayerTrainingFunction {
		AllocatingMatrix _distortedInput;
		Random & _rand;

	protected:
		virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
			return _dbn.UpdateAE(inputMiniBatch, _inputBiases, outputMiniBatch, _layer, stepSize);
		}

	public:		
		AETrainingFunction(DBN & dbn, int layer, const Matrix & trainData, Random & rand, const HyperParams & hyperParams, Layer::ActFunc lowerActFunc, bool fixTrainDistortion)
			: LayerTrainingFunction(dbn, layer, trainData, hyperParams), _rand(rand)
		{
			if (fixTrainDistortion) {
				_distortedInput.Resize(trainData);
				lowerActFunc.Sample(trainData.Vec(), _distortedInput.Vec(), _rand);
				_input = _distortedInput;
			}
		}
	};

	class BPTrainingFunction : public TrainingFunction {
		AllocatingMatrix _errors;

	private:
		ObjectiveType _objectiveType;

	protected:
		virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
			return _dbn.UpdateBackProp(inputMiniBatch, outputMiniBatch, _objectiveType, false, stepSize);
		}

	public:
		BPTrainingFunction(const Matrix & input, const Matrix & output, DBN & dbn, HyperParams hyperParams, ObjectiveType objectiveType)
			:
		TrainingFunction(dbn, input, output, hyperParams, dbn._params, dbn._deltaParams, dbn._savedParams),
			_objectiveType(objectiveType)
		{ }
	};

	class OutputLayerTrainingFunction : public TrainingFunction {
		AllocatingMatrix _errors;

	private:
		ObjectiveType _objectiveType;

	protected:
		virtual double PrivateEval(Matrix inputMiniBatch, Matrix outputMiniBatch, double stepSize) {
			return _dbn.UpdateBackProp(inputMiniBatch, outputMiniBatch, _objectiveType, true, stepSize);
		}

	public:
		OutputLayerTrainingFunction(const Matrix & input, const Matrix & output, DBN & dbn, HyperParams hyperParams, ObjectiveType objectiveType)
			:
		TrainingFunction(dbn, input, output, hyperParams, dbn._layerParams.back(), dbn._layerDeltaParams.back(), dbn._layerSavedParams.back()),
			_objectiveType(objectiveType)
		{ }

		virtual void Init(Random & rand) {
			// sparse initialization strategy from Martens, 2010
			_dbn._B.back() *= 0;
			MutableMatrix W = _dbn._W.back();
			W *= 0;
			for (int c = 0; c < W.NumC(); ++c) {
				MutableVector col = W.GetCol(c);
				// sampling with replacement, but it doesn't really matter
				for (int i = 0; i < 15; ++i) {
					col[rand.Rand(col.Len())] = rand.Normal() * 0.01;
				}
			}

			TrainingFunction::Init(rand);
		}
	};

	void TrainSGD(TrainingFunction & trainer, Random & rand, const HyperParams & hyperParams) {
		double stepSize = hyperParams.initStepSize;
		bool quiet = false;

		if (!quiet) cout << fixed << setprecision(4);

		trainer.Init(rand);

		// get a baseline
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

			stepSize /= 2;
			if (!quiet) cout << "Step size reduced to " << stepSize << endl;
			trainer.Restore();
		}

		double totalScore = 0;
		for (int t = 1; t <= hyperParams.numUpdates; ++t) {
			// momentum increases from 0.5 to max
			// using formula from Sutskever et al. 2013
			double tFrac = (double)t / hyperParams.numUpdates;
			double momentum = 1.0 - 1.0 / (2.0 + tFrac * ( 1.0 / (1.0 - hyperParams.maxMomentum) - 2.0));

			totalScore += trainer.Update(startInstance, stepSize, hyperParams.maxUpdate, momentum);
			startInstance += hyperParams.miniBatchSize;

			if (!quiet && t % hyperParams.checkInterval == 0) {
				double score = totalScore / hyperParams.checkInterval;
				cout << "Step: " << t << endl;
				cout << "Momentum: " << momentum << endl;
				cout << "New score: " << score << endl << endl;
				totalScore = 0;
			}
		}

		totalScore = 0;
		for (int t = 0; t < hyperParams.numAnnealUpdates; ++t) {
			// momentum increases from 0.5 to max
			// using formula from Sutskever et al. 2013
			double tFrac = (double)t / hyperParams.numAnnealUpdates;
			double momentum = (1 - tFrac) * hyperParams.maxMomentum + tFrac * 0.9;
			double step = (1 - tFrac) * stepSize;

			totalScore += trainer.Update(startInstance, step, hyperParams.maxUpdate, momentum);
			startInstance += hyperParams.miniBatchSize;
		}

		if (!quiet) {
			double score = totalScore / hyperParams.numAnnealUpdates;
			cout << "Anneal score: " << score << endl << endl;
		}

	}

	double UpdateBackProp(const Matrix & input, const Matrix & output, ObjectiveType _objType, bool lastLayerOnly, double stepSize) {
		int startLayer = lastLayerOnly ? _numLayers - 1 : 0;

		Matrix mappedInput = MapUp(input, startLayer);
		double loss = 0;

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
			Matrix mappedError = _layers.back().ComputeErrors(_tempMat, _oActFunc);
			for (int l = _numLayers - 1; ; --l) {
				for (int i = 0; i < input.NumC(); ++i) _deltaB[l] += mappedError.GetCol(i) * stepSize;

				const Matrix & lowerProbs = l > startLayer ? _layers[l - 1].Activations() : input;

				_deltaW[l] += stepSize * lowerProbs * mappedError.Trans();

				if (lastLayerOnly || l == 0) break;

				_tempMat = _W[l] * mappedError;
				mappedError = _layers[l-1].ComputeErrors(_tempMat, _hActFunc);
			}
		}

		return loss;
	}

	static void InitializeInputBiases(const Matrix & input, AllocatingVector & inputBiases, double eps) {
		int n = input.NumC();
		inputBiases.Resize(input.NumR());
		//		inputBiases.MultInto(input, false, AllocatingVector(n, 1.0), 1.0/n);
		for (int r = 0; r < input.NumR(); ++r) {
			Vector row = input.GetRow(r);
			inputBiases[r] = Sum(row) / n;
		}
	}

	double UpdateCD(const Matrix & input, const Vector & inputBiases, int layer, Random & rand, double stepSize) {
		Layer::ActFunc lowerActFunc = (layer == 0) ? _iActFunc : _hActFunc;

		Matrix weights = _W[layer];
		Vector topBiases = _B[layer];
		Layer & topLayer = _layers[layer];

		_tempBottomSample.Resize(input);

		topLayer.ActivateUp(weights, topBiases, input, _hActFunc);
		_tempTopSample.Resize(topLayer.Size(), input.NumC());
		topLayer.Sample(_tempTopSample, rand, _hActFunc);

		_tempBottomLayer.ActivateDown(weights, inputBiases, _tempTopSample, lowerActFunc);
		_tempBottomLayer.Sample(_tempBottomSample, rand, lowerActFunc);
		_tempTopLayer.ActivateUp(weights, topBiases, _tempBottomSample, _hActFunc);

		if (stepSize > 0) {
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

	double UpdateAE(const Matrix & input, const Vector & inputBiases, const Matrix & targetInput, int layer, double stepSize) {
		Layer::ActFunc lowerActFunc = (layer == 0) ? _iActFunc : _hActFunc;

		Matrix weights = _W[layer];
		Vector topBiases = _B[layer];

		Layer & topLayer = _layers[layer];
		Matrix topProbs = topLayer.ActivateUp(weights, topBiases, input, _hActFunc);

		double negll = _tempBottomLayer.ActivateDownAndGetNegLL(weights, inputBiases, topProbs, targetInput, lowerActFunc);

		if (stepSize > 0) {
			// fill bottomProbs with errors
			MutableMatrix & bottomProbs = _tempBottomLayer.Activations();
			bottomProbs -= targetInput;

			_deltaW[layer] += stepSize * bottomProbs * topProbs.Trans();

			_tempMat = weights.Trans() * bottomProbs;
			Matrix topError = topLayer.ComputeErrors(_tempMat, _hActFunc);

			_deltaW[layer] += stepSize * input * topError.Trans();

			for (int i = 0; i < input.NumC(); ++i) {
				_deltaB[layer] += stepSize * topError.GetCol(i);
			}
		}

		return negll;
	}

public:
	DBN(int numLayers, int iSize, int hSize, int oSize, Layer::ActFunc iActFunc, Layer::ActFunc hActFunc, Layer::ActFunc oActFunc) 
	{
		Initialize(numLayers, iSize, hSize, oSize, iActFunc, hActFunc, oActFunc);
	}

	void Initialize(int numLayers, int iSize, int hSize, int oSize, Layer::ActFunc iActFunc, Layer::ActFunc hActFunc, Layer::ActFunc oActFunc)
	{
		_numLayers = numLayers;
		_iSize = iSize;
		_hSize = hSize;
		_oSize = oSize;

		_iActFunc = iActFunc;
		_hActFunc = hActFunc;
		_oActFunc = oActFunc;

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
			int tSize = (i == _numLayers - 1) ? _oSize : _hSize;
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

	void Deserialize(istream & inStream) {
		inStream.read((char *) &_numLayers, sizeof(int));
		inStream.read((char *) &_iSize, sizeof(int));
		inStream.read((char *) &_hSize, sizeof(int));
		inStream.read((char *) &_oSize, sizeof(int));
		inStream.read((char *) &_iActFunc, sizeof(Layer::ActFunc));
		inStream.read((char *) &_hActFunc, sizeof(Layer::ActFunc));
		inStream.read((char *) &_oActFunc, sizeof(Layer::ActFunc));

		Initialize(_numLayers, _iSize, _hSize, _oSize, _iActFunc, _hActFunc, _oActFunc);

		inStream.read((char *) _params.Start(), sizeof(double) * _params.Len());
	}

	void Serialize(ostream & outStream) const {
		outStream.write((const char *) &_numLayers, sizeof(int));
		outStream.write((const char *) &_iSize, sizeof(int));
		outStream.write((const char *) &_hSize, sizeof(int));
		outStream.write((const char *) &_oSize, sizeof(int));
		outStream.write((const char *) &_iActFunc, sizeof(Layer::ActFunc));
		outStream.write((const char *) &_hActFunc, sizeof(Layer::ActFunc));
		outStream.write((const char *) &_oActFunc, sizeof(Layer::ActFunc));
		outStream.write((const char *) _params.Start(), sizeof(double) * _params.Len());
	}

	void Randomize(Random & rand, double stdDev) {
		_params.Apply([&] (double x) { return rand.Normal() * stdDev; });
	}

	int NumLayers() const { return _numLayers; }

	int NumParams() const {
		int numParams = 0;
		for (int layer = 0; layer < _numLayers; ++layer) numParams += NumParamsInLayer(layer);

		return numParams;
	}

	int NumParamsInLayer(int layer) const {
		assert (layer < _numLayers);
		if (_numLayers == 1) return (_iSize + 1) * _oSize;
		else if (layer == 0) return (_iSize + 1) * _hSize;
		else if (layer == _numLayers - 1) return (_hSize + 1) * _oSize;
		else return (_hSize + 1) * _hSize;
	}

	int LayerOutSize(int layer) const {
		assert (layer < _numLayers);
		return (layer < _numLayers - 1) ? _hSize : _oSize;
	}

	int LayerInSize(int layer) const {
		assert (layer < _numLayers);
		return (layer == 0) ? _iSize : _hSize;
	}

	Matrix MapLayer(const Matrix & input, int layer) const {
		Layer::ActFunc actFunc = (layer == _numLayers - 1) ? _oActFunc : _hActFunc;
		return _layers[layer].ActivateUp(_W[layer], _B[layer], input, actFunc);
	}

	Matrix MapUp(const Matrix & input, int startLayer = -1, int endLayer = -1) const {
		if (startLayer == -1) startLayer = 0;
		if (endLayer == -1) endLayer = _numLayers;

		Matrix mappedInput = input;

		for (int layer = startLayer; layer < endLayer; ++layer) mappedInput = MapLayer(mappedInput, layer);

		return mappedInput;
	}

	void Train(const Matrix & input, const Matrix & output, ObjectiveType objectiveType, Random & rand, bool quiet, const vector<HyperParams> & hyperParams_pt, const HyperParams & hyperParams_bp) {
		Matrix mappedInput = input;

		// pretrain hidden layers
		for (int layer = 0; layer < _layers.size() - 1; ++layer) {
			const HyperParams & hyperParams = hyperParams_pt[layer];
			InitLayer(layer, rand);

			switch (hyperParams.pretrainType) {
			case CD:
				{
					if (!quiet) {
						cout << "Pretraining layer " << layer << " of size " << LayerInSize(layer) << "x" << LayerOutSize(layer) << " with CD" << endl;
					}
					CDTrainingFunction cdFunc(*this, layer, mappedInput, rand, hyperParams);
					TrainSGD(cdFunc, rand, hyperParams);
				}
				break;

			case AE:
				{
					if (!quiet) {
						cout << "Pretraining layer " << layer << " of size " << LayerInSize(layer) << "x" << LayerOutSize(layer) << " with AE" << endl;
					}
					Layer::ActFunc lowerActFunc = (layer == 0) ? _iActFunc : _hActFunc;
					AETrainingFunction aeFunc(*this, layer, mappedInput, rand, hyperParams, lowerActFunc, true);
					TrainSGD(aeFunc, rand, hyperParams);
				}
				break;

			case NONE:
				break;

			default:
				abort();
			}

			mappedInput = MapLayer(mappedInput, layer);
			if (layer > 0) _layers[layer - 1].Clear(); // save some memory
		}

		if (hyperParams_pt.back().pretrainType != NONE) {
			// pretrain output layer
			if (!quiet) {
				cout << "Pretraining output layer" << endl;
			}

			OutputLayerTrainingFunction outputFunc(mappedInput, output, *this, hyperParams_pt.back(), objectiveType);
			TrainSGD(outputFunc, rand, hyperParams_pt.back());
		}

		// and now, fine tune

		if (!quiet) {
			cout << "Fine tuning" << endl;
		}

		BPTrainingFunction bpFunc(input, output, *this, hyperParams_bp, objectiveType);
		TrainSGD(bpFunc, rand, hyperParams_bp);
	}
};
