#pragma once

#include <iostream>

#include "Matrix.h"
#include "Random.h"

using namespace std;

void TrainSGD2(Random & rand) {
	cout << "TrainSGD2" << endl;
}

template <class Func>
void TrainSGD(Func & func, MutableVector & params, int miniBatchSize, int checkInterval, bool uniformSamples, Random & rand) {
	int maxSteps = 1000000;

	// TODO: make this an argument
	bool quiet = false;

	int count = func.NumIns();
	if (miniBatchSize == -1) miniBatchSize = count;
	if (checkInterval < 0) checkInterval *= -count / miniBatchSize;

	double maxUpdate = 1e-2;

	AllocatingVector savedParams(params), writeParams(params);
	vector<int> updateIndices(miniBatchSize);
	double stepSize = (double)miniBatchSize / count;
	const double minImprovement = 1e-4;
	const int minStepSizeChanges = 5;
	if (!quiet) printf("Setting initial step size to %8.4e\n", stepSize);

	bool _useOnlineEval = true;
	
RESTART:
	double bestValLastStepSizeChange = INFTY;
	int numStepSizeChanges = 0;
			
	double bestVal = func.Eval(params);
	if (!quiet) printf("\nbaseline score: %-8.4f\n", bestVal);

	int epoch = 0;
	int startMiniBatch = 0;
	double totalVal = 0;
	int numEvals = 0;

	for (int step = 1; step != maxSteps; ++step) {
		if (uniformSamples) {
			for (int i = 0; i < miniBatchSize; ++i) {
				updateIndices[i] = rand.Rand(count);
			}
			totalVal += func.Eval(params, writeParams, updateIndices);
		} else {
			totalVal += func.Eval(params, writeParams, startMiniBatch, miniBatchSize);
			startMiniBatch += miniBatchSize;
		}
		numEvals++;

		writeParams *= stepSize;
		writeParams.Trunc(maxUpdate);

		params -= writeParams;

		if (!params.AllSafe()) {
			params.CopyFrom(savedParams);
			stepSize /= 2;
			if (!quiet) printf("Got NaN. Reducing step size to %8.4f\n", stepSize);
			goto RESTART;
		}

		if (step % checkInterval == 0) {
			++epoch;

			double val;
			if (_useOnlineEval) {
				val = totalVal / numEvals;
				totalVal = 0;
				numEvals = 0;

				if (miniBatchSize > count / 2 && epoch == 1) continue;
			} else {
				val = func.Eval(params);
			}

			if (!quiet) {
				printf("%-5d", epoch);
				printf("norm: %8.4f", params.Norm());
				printf("  val: %8.4f\n", val);
			}
		
			if (val >= bestVal * (1.0 - minImprovement) || IsDangerous(val)) {
				if (numStepSizeChanges >= minStepSizeChanges && bestVal > (1.0 - minImprovement) * bestValLastStepSizeChange) break;
				bestValLastStepSizeChange = bestVal;
				numStepSizeChanges++;

				stepSize /= 2;
				if (epoch == 1) stepSize /= 5;
				if (val > bestVal) params.CopyFrom(savedParams);
				if (!quiet)
					printf("stepSize reduced to %8.4e\n", stepSize);

				if (epoch == 1 && (val >= bestVal || IsDangerous(val))) {
					if (!quiet) printf("Resetting.\n");
					goto RESTART;
				}
			}

			if (!IsDangerous(val) && val < bestVal) {
				savedParams.CopyFrom(params);
				bestVal = val;
			}
		}
	}
}
