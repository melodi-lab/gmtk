#pragma once

/*
 * This file, MatrixFunc.h, defines several functions involving
 * the vector classes defined in Matrix.h
 *
 * Written by Galen Andrew gmandrew@uw.edu
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "Globals.h"
#include "Matrix.h"

// returns the maximally valued element
double Max(const Vector & vec);

// Computes the log sum of the entries of vec
// log sum_i exp(x_i)
// in a particularly fast and stable way that
// does not require temporary storage, but
// overwrites (destroys) the original data in vec
double LogSumFastDestroy(MutableVector & vec);

// Replaces each element x_i of the vector with
// log(exp(x_i) + exp(-x_i))
// in a particularly fast and stable way requiring
// temporary storage in the vector temp
void LogSumWithNeg(MutableVector & vec, MutableVector & temp);

// Projection onto L-infinity ball of radius maxVal.
// Replaces any element of vec that exceeds maxVal in absolute value
// with the signed value maxVal
// maxVal is assumed to be positive
void Trunc(MutableVector & vec, double maxVal);

// For each element x, if x is larger than d in magnitude, reduces
// the value by d. If x is smaller than d, it is set to zero.
// d is assumed to be positive.
void Shrink(MutableVector & vec, double d);

// Returns the L2 norm of the vector
double Norm(const Vector & vec);

// Print the values of vector to the console
void Print(const Vector & vec);

// Return the index of the maximal element
int ArgMax(const Vector & vec);

// Returns the L1 norm of the vector
double NormL1(const Vector & vec);

// Returns the sum of the vector elements
double Sum(const Vector & vec);

// L1 or L2 weight decay
// takes a step in the direction of the negative (sub-)gradient
// of the function 0.5 * decayRate ||vec||_2^2 in case decayType==2
// or decayRate ||vec||_1 in case decayType==1
// and returns the value of that function
// For L2, this is functionally equivalent to multiplication by
// (1 - stepSize * decayRate)
double Decay(MutableVector vec, double stepSize, double decayRate, int decayType = 2);
