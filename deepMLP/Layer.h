#pragma once

#include <assert.h>
#include <vector>

#include "Matrix.h"
#include "Random.h"

using namespace std;

class Layer {
public:
	struct ActFunc {
		enum ActType { LOG_SIG, TANH, CUBIC, LINEAR, RECT_LIN } actType;

		ActFunc() : actType(LINEAR) { }

		ActFunc(ActType actType) : actType(actType) { }

		static double CubicSigmoid(double y) {
			// I have a vectorized implementation of this using SSE but
			// it uses VS intrinsics
			bool negate = false;
			if (y < 0) { negate = true; y = -y; }
			double x = (y <= 20) ? y : pow(3 * y, 1.0/3);

			double newX;
			while (true) {
				double xSqr = x * x;
				newX = (0.66666666666666666 * xSqr * x + y) / (xSqr + 1.0);
				if (newX >= x) break;
				x = newX;
			}

			return negate ? -newX : newX;
		}

		void Apply(const MutableVector & inputs) {
			switch (actType) {
			case LOG_SIG:
				{
					auto func = [] (double x)->double {
						if (x < -30) return 0;
						else if (x > 30) return 1;
						else return 1.0 / (1 + exp(-x));
					};
					inputs.Apply(func);
				}
				break;

			case TANH:
#if HAVE_MKL
				inputs.ApplyVML(vdTanh);
#else
				inputs.Apply([](double x)->double {return tanh(x);});
#endif
				break;

			case CUBIC:
				inputs.Apply(CubicSigmoid);
				break;

			case LINEAR: break;

			case RECT_LIN:
				{
					auto func = [] (double x) { return max(x, 0.0); };
					inputs.Apply(func);
				}
				break;

			default: abort();
			}
		}

		double ApplyAndGetNegLL(const MutableVector & inputs, const Vector & labels) {
			double negll = 0;

			switch (actType) {
			case TANH:
				{
					auto func = [&](double aVal, double xVal)->double {
						if (aVal < -14) {
							double xP = (1 + xVal) / 2;
							negll -= 2 * xP * aVal;
							return -1;
						} else if (aVal > 14) {
							double xP = (1 + xVal) / 2;
							negll += 2 * (1 - xP) * aVal;
							return 1;
						} else {
							double a = tanh(aVal);
							double p = (1 + a) / 2;
							double xP = (1 + xVal) / 2;
							negll -= xP * log(p) + (1 - xP) * log (1.0 - p);
							return a;
						}
					};
					inputs.Apply(func, labels);
				}
				break;

			case LINEAR:
			case CUBIC:
				{
					auto func = [&](double aVal, double xVal) {
						negll += (aVal - xVal) * (aVal - xVal);
					};
					inputs.Visit(func, labels);
					negll /= 2;
				}
				break;

			case LOG_SIG:
				{
					auto func = [&](double aVal, double xVal)->double {
						if (aVal < -30) {
							negll -= xVal * aVal;
							return 0;
						} else if (aVal > 30) {
							negll += (1 - xVal) * aVal;
							return 1;
						} else {
							double p = 1.0 / (1.0 + exp(-aVal));
							negll -= xVal * log(p) + (1 - xVal) * log (1.0 - p);
							return p;
						}
					};
					inputs.Apply(func, labels);
				}
				break;

			default:
				abort();
			}

			return negll;
		}

		void ComputeErrors(const MutableVector & activations, const Vector & inError) {
			switch (actType) {

			case LOG_SIG:
				{
					auto func = [](double aVal, double eVal) { return (1.0 - aVal) * aVal * eVal; };
					activations.Apply(func, inError);
				}
				break;

			case TANH:
				{
					auto func = [](double aVal, double eVal) { return (1.0 - aVal * aVal) * eVal; };
					activations.Apply(func, inError);
				}
				break;

			case CUBIC:
				{
					auto func = [](double aVal, double eVal) { return eVal / (1.0 + aVal * aVal); };
					activations.Apply(func, inError);
				}
				break;

			case LINEAR:
				activations.CopyFrom(inError);
				break;

			case RECT_LIN:
				{
					auto func = [](double aVal, double eVal) { return (aVal <= 0) ? 0 : eVal; };
					activations.Apply(func, inError);
				}
				break;

			default: abort();
			}
		}

		void Sample(const Vector & activations, const MutableVector & sample, Random & rand) const {
			switch (actType) {
			case LOG_SIG:
				{
					auto func = [&](double aVal) { return (rand.Uniform() < aVal) ? 1.0 : 0.0; };
					sample.Replace(func, activations);
				}
				break;

			case TANH:
				{
					auto func = [&](double aVal) { return (2 * rand.Uniform() - 1) < aVal ? 1.0 : -1.0; };
					sample.Replace(func, activations);
				}
				break;

			case LINEAR:
			case CUBIC:
				{
					auto func = [&](double aVal) { return aVal + rand.Normal(); };
					sample.Replace(func, activations);
				}
				break;

			default:
				abort();
			}
		}
	};

private:
	AllocatingMatrix _a;
	
	void ComputeInputs(const Matrix & weights, const Vector & biases, const Matrix & values, bool trans) {
		int numIns = values.NumC();
		int size = trans ? weights.NumC() : weights.NumR();

		_a = (trans ? weights.Trans() : weights) * values;
		for (int i = 0; i < numIns; ++i) _a.GetCol(i) += biases;
	}

public:
	Layer() { }

	const Matrix & ActivateUp(const Matrix & weights, const Vector & biases, const Matrix & lowerValues, ActFunc actFunc) {
		ComputeInputs(weights, biases, lowerValues, true);
		actFunc.Apply(_a.Vec());
		return _a;
	}
	
	double ActivateDownAndGetNegLL(const Matrix & weights, const Vector & biases, const Matrix & upperValues, const Matrix & lowerValues, ActFunc actFunc) {
		ComputeInputs(weights, biases, upperValues, false);
		return actFunc.ApplyAndGetNegLL(_a.Vec(), lowerValues.Vec());
	}
	
	const Matrix & ActivateDown(const Matrix & weights, const Vector & biases, const Matrix & upperValues, ActFunc actFunc) {
		ComputeInputs(weights, biases, upperValues, false);
		actFunc.Apply(_a.Vec());
		return _a;
	}

	const Matrix & ComputeErrors(const Matrix & inError, ActFunc actFunc) {
		actFunc.ComputeErrors(_a.Vec(), inError.Vec());
		return _a;
	}

	MutableMatrix & Activations() { return _a; }

	int Size() const { return _a.NumR(); }

	int Count() const { return _a.NumC(); }
	
	void Sample(MutableMatrix & sample, Random & rand, ActFunc actFunc) const {
		actFunc.Sample(_a.Vec(), sample.Vec(), rand);
	}

	void Clear() {
		_a.Resize(0, 0);
	}
};
