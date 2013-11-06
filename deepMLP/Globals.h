#pragma once

/*********
This file, Globals.h, contains some basic functions that are used all
over the NN training code.
 *********/

#include <emmintrin.h>
#include <xmmintrin.h>
#include <limits>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

// swap the values of a and b
template<class T>
void Swap(T a, T b) {
  T t = a;
  a = b;
  b = t;
}

// macro defining hard assertion: aborts program if condition fails
#define ASSERT(TST) ( (TST) ? (void)0 : (std::cerr << __FILE__ "(" << __LINE__	<< "): Assertion failed " #TST << std::endl,abort()) )

static const double INFTY = std::numeric_limits<double>::infinity();

static const double NaN = std::numeric_limits<double>::quiet_NaN();

static const double TOL = pow(std::numeric_limits<double>::epsilon(), (double)1.0 / 3);

static bool IsClose(double a, double b) {
  return fabs(a - b) < TOL;
}

static bool IsNaN(double x) { return std::isnan(x); }

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

// Assuming the class C implements a method
// void Serialize(std::ofstream)
// Serialize c to a file with the given path
template<class C>
void Serialize(const C & c, const string & filename) {
  ofstream outStream(filename, ios::out|ios::binary);
  if (!outStream.is_open()) {
    cout << "Couldn't open serialized file " << filename.c_str() << endl;
    abort();
  }

  c.Serialize(outStream);

  outStream.close();
}

// Assuming the class C implements a method
// void Deserialize(std::ofstream)
// Deserialize c from a file with the given path
template<class C>
void Deserialize(C & c, const string & filename) {
  ifstream inStream(filename, ios::in|ios::binary);
  if (!inStream.is_open()) {
    cout << "Couldn't open serialized file " << filename.c_str() << endl;
    abort();
  }

  c.Deserialize(inStream);

  inStream.close();
}
