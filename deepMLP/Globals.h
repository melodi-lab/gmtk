#pragma once

/*********
 * This file, Globals.h, contains some basic functions that are used all
 * over the NN training code.
 *
 * Written by Galen Andrew gmandrew@uw.edu
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *********/

#include <emmintrin.h>
#include <xmmintrin.h>
#include <limits>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>

#include "error.h"

using namespace std;

// swap the values of a and b
template<class T>
void Swap(T a, T b) {
  T t = a;
  a = b;
  b = t;
}

static const double INFTY = std::numeric_limits<double>::infinity();

static const double NaN = std::numeric_limits<double>::quiet_NaN();

static const double TOL = pow(std::numeric_limits<double>::epsilon(), (double)1.0 / 3);

#if 0
// unused
static bool IsClose(double a, double b) {
  return fabs(a - b) < TOL;
}
#endif

#pragma GCC diagnostic ignored "-Wunused-function"
static bool IsNaN(double x) { return std::isnan(x); }

#if 0
// unused
static bool IsInf(double x) { return std::isinf(x); }

static bool IsDangerous(double x) { return IsNaN(x) || IsInf(x); }

// Computes the logsum of x and y, that is
// log(exp(x) + exp(y))
// in an efficient and numerical stable way
static double LogSum(double x, double y) {
  double d = x - y;
  if (d < -30) return y;
  else if (d > 30) return x;
  else if (d > 0) return x + log(1.0 + exp(-d));
  else return y + log(1.0 + exp(d));
}

// Computes the logistic sigmoid of x
// 1/(1+exp(-x))
// in an efficient and numerical stable way
static double Logistic(double x) {
  if (x < -30) return 0;
  else if (x > 30) return 1;
  else return 1.0 / (1.0 + exp(-x));
}

// Computes the negative log of the logistic sigmoid
// log(1+exp(-x))
// in an efficient and numerical stable way
static double LogLoss(double x) {
  if (x < -30) return -x;
  else if (x > 30) return 0;
  else return log(1 + exp(-x));
}
#endif

// Assuming the class C implements a method
// void Serialize(std::ofstream)
// Serialize c to a file with the given path
template<class C>
void Serialize(const C & c, const string & filename) {
  ofstream outStream(filename.c_str(), ios::out|ios::binary);
  if (!outStream.is_open()) {
    error("Error: Couldn't open file '%s' to serialize to\n", filename.c_str());
  }

  c.Serialize(outStream);

  outStream.close();
}

// Assuming the class C implements a method
// void Deserialize(std::ofstream)
// Deserialize c from a file with the given path
template<class C>
void Deserialize(C & c, const string & filename) {
  ifstream inStream(filename.c_str(), ios::in|ios::binary);
  if (!inStream.is_open()) {
    error("Error: Couldn't open serialized file '%s'\n", filename.c_str());
  }

  c.Deserialize(inStream);

  inStream.close();
}
