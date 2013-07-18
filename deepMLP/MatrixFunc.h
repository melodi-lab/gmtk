#pragma once

#include "Globals.h"
#include "Matrix.h"

double Max(const Vector & vec);

double LogSumFastDestroy(MutableVector & vec);

void LogSumWithNeg(MutableVector & vec, MutableVector & temp);

void Trunc(MutableVector & vec, double maxVal);

void Shrink(MutableVector & vec, double d);

double Norm(const Vector & vec);

void Print(const Vector & vec);

int ArgMax(const Vector & vec);

double NormL1(const Vector & vec);

double Sum(const Vector & vec);

double Decay(MutableVector vec, double stepSize, double decayRate, int decayType = 2);