#pragma once

/*
 * Written by Galen Andrew gmandrew@uw.edu
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include <assert.h>
#include <vector>

#include "Matrix.h"
#include "rand.h"

using namespace std;

// The Layer class holds the activations of a layer
// and is responsible for computing them and backpropagating error
class Layer {
public:

	// The ActFunc struct stores a type of activation function
	// and performs the operations pertinant to it
	// 
	// possible activation functions:
	// LOG_SIG: logistic sigmoid (1+exp(-x))^{-1}
	// TANH: hyperbolic tangent
	// CUBIC: modified cube root, described in Andrew et al., 2013
	// LINEAR: identity function
	// RECT_LIN: rectified linear activation
	// 
	// Operations include
	//   evaluating the function,
	//   multiplying by the derivative (for backpropagation)
	//   sampling a value, if it is interpreted as the parameter of a distribution
  //      (LINEAR and CUBIC are mean of Gaussian, LOG_SIG is p of Bernoulli)
  struct ActFunc {
    enum ActType { LOG_SIG, TANH, CUBIC, LINEAR, RECT_LIN } actType;

  ActFunc() : actType(LINEAR) {}

  ActFunc(ActType actType) : actType(actType) {}

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

		// compute function value on all inputs
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

      default: assert(false); // don't construct illegal NNs...
      }
    }

		// compute function on all inputs, and return negative log-likelihood
		// of single-parameter distribution evaluated on label vector
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
        {
          auto func = [&](double aVal, double xVal) {
            negll += (aVal - xVal) * (aVal - xVal);
          };
          inputs.Visit(func, labels);
          negll /= 2;
        }
        break;

      case CUBIC:
        {
          auto func = [&](double aVal, double xVal)->double {
	    double res = CubicSigmoid(aVal);
            negll += (res - xVal) * (res - xVal);
	    return res;
          };
          inputs.Apply(func, labels);
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

      case RECT_LIN:
	{
          auto func = [&](double aVal, double xVal)->double {
	    double res = max(aVal, 0.0);
            negll += (res - xVal) * (res - xVal);
	    return res;
          };
          inputs.Apply(func, labels);
          negll /= 2;
	}
	break;

      default: assert(false); // don't construct illegal NNs...
      }

      return negll;
    }

		// Given incoming error vector (derivative of loss wrt output)
		// multiply by derivative to produce derivative of loss wrt input
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
	  auto func = [](double aVal, double eVal) {
	    // if activation is numerically equal to zero, assume it is a dropped-out unit
	    return (aVal == 0.0) ? 0 : (1.0 - aVal * aVal) * eVal;
	  };
          activations.Apply(func, inError);
        }
        break;

      case CUBIC:
        {
	  auto func = [](double aVal, double eVal) {
	    // if activation is numerically equal to zero, assume it is a dropped-out unit
	    return (aVal == 0.0) ? 0 : eVal / (1.0 + aVal * aVal);
	  };
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

      default: assert(false); // don't construct illegal NNs...
      }
    }

    void Sample(const Vector & activations, const MutableVector & sample) const {
      switch (actType) {
      case LOG_SIG:
        {
          auto func = [&](double aVal) { return (rnd.drand48() < aVal) ? 1.0 : 0.0; };
          sample.Replace(func, activations);
        }
        break;

      case TANH:
        {
          auto func = [&](double aVal) { return (2 * rnd.drand48() - 1) < aVal ? 1.0 : -1.0; };
          sample.Replace(func, activations);
        }
        break;

      case RECT_LIN:
      case LINEAR:
      case CUBIC:
        {
          auto func = [&](double aVal) { return aVal + rnd.normal(); };
          sample.Replace(func, activations);
        }
        break;

      default: assert(false); // don't construct illegal NNs...
      }
    }
  };

private:
  AllocatingMatrix _a;

	// compute the input values by multiplying by weights and adding biases
  void ComputeInputs(const Matrix & weights, const Vector & biases, const Matrix & values, bool trans) {
    int numIns = values.NumC();
    //    int size = trans ? weights.NumC() : weights.NumR();

    _a = (trans ? weights.Trans() : weights) * values;
    for (int i = 0; i < numIns; ++i) _a.GetCol(i) += biases;
  }

public:
  Layer() { }

	// compute inputs and apply activation function
  Matrix ActivateUp(const Matrix & weights, const Vector & biases, const Matrix & lowerValues, ActFunc actFunc) {
    ComputeInputs(weights, biases, lowerValues, true);
    actFunc.Apply(_a.Vec());
    return _a;
  }

	// dropout values
  const Matrix & Dropout(double dropP) {
    _a.Vec().Apply([&] (double a) { return rnd.drand48() > dropP ? a : 0; });
    return _a;
  }

	// compute activation downward, and return the negative log-likelihood, for autoencoder training
  double ActivateDownAndGetNegLL(const Matrix & weights, const Vector & biases, const Matrix & upperValues, const Matrix & lowerValues, ActFunc actFunc) {
    ComputeInputs(weights, biases, upperValues, false);
    return actFunc.ApplyAndGetNegLL(_a.Vec(), lowerValues.Vec());
  }
	
	// compute activations downward
  const Matrix & ActivateDown(const Matrix & weights, const Vector & biases, const Matrix & upperValues, ActFunc actFunc) {
    ComputeInputs(weights, biases, upperValues, false);
    actFunc.Apply(_a.Vec());
    return _a;
  }

	// multiply errors by derivative, for backpropagation
  const Matrix & ComputeErrors(const Matrix & inError, ActFunc actFunc) {
    actFunc.ComputeErrors(_a.Vec(), inError.Vec());
    return _a;
  }

  MutableMatrix & Activations() { return _a; }

  int Size() const { return _a.NumR(); }

  int Count() const { return _a.NumC(); }

  void Sample(MutableMatrix & sample, ActFunc actFunc) const {
    actFunc.Sample(_a.Vec(), sample.Vec());
  }

  void Clear() {
    _a.Resize(0, 0);
  }
};
