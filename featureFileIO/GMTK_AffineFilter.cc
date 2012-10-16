
/*
 * GMTK_AffineFilter.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif

#include "error.h"
#include "general.h"
VCID(HGID)

#include <string.h>
#include <stdio.h>

#include "GMTK_AffineFilter.h"

AffineFilter::AffineFilter(char const *fileName, Filter *nextFilter) 
  : Filter(nextFilter)
{
  FILE *f = fopen(fileName, "r");
  if (!f) {
    error("AffineFilter: cannot open '%s' for reading\n", fileName);
  }
  if (fscanf(f, "%u %u", &rows, &cols) != 2) {
    error("AffineFilter: error reading affine filter size from '%s'\n", fileName);
  }

  c = new double[cols];
  if (!c) {
    error("AffineFilter: failed to allocate constant vector\n");
  }
  unsigned Bsize = rows * cols;
  B = new double[Bsize];
  if (!B) {
    error("AffineFilter: failed to allocate coefficient matrix\n");
  }
  for (unsigned i=0; i < Bsize; i+=1) {
    if (fscanf(f,"%lf", B+i) != 1) {
      error("AffineFilter: error reading coefficient matrix\n");
    }
  }
  for (unsigned i=0; i < cols; i+=1) {
    if (fscanf(f,"%lf", c+i) != 1) {
      error("AffineFilter: error reading constant vector\n");
    }
  }
  fclose(f);
  buffer=NULL; buffSize=0;
  workBuffer=NULL; workBuffSize=0;
}


subMatrixDescriptor *
AffineFilter::getRequiredInput(unsigned first, unsigned count, 
			    unsigned inputContinuous, unsigned inputDiscrete,
			    unsigned inputTotalFrames)
{
  subMatrixDescriptor *nextFilterInput = NULL;
  if (nextFilter) {
    unsigned outputContinuous = cols;
    unsigned outputDiscrete   = inputDiscrete;
    unsigned outputFrames     = inputTotalFrames;
    nextFilterInput = nextFilter->getRequiredInput(first, count, 
						   outputContinuous, outputDiscrete,
						   outputFrames);
    assert(nextFilterInput);
    first = nextFilterInput->requestedFirst;
    count = nextFilterInput->requestedCount;
  } 
  return subMatrixDescriptor::getSMD(first, count, 0, 0, 
				     inputContinuous, inputDiscrete, inputTotalFrames, 
				     first, count, nextFilterInput);
}


subMatrixDescriptor
AffineFilter::describeLocalOutput(subMatrixDescriptor const &inputDescription) {
  subMatrixDescriptor myOutput = inputDescription;
  myOutput.numContinuous = cols;
  myOutput.historyFrames = 0;
  myOutput.futureFrames = 0;
  myOutput.next = NULL;
  return myOutput;
}

extern "C" void mul_mfmf_mf(const int M, const int K, const int N, 
		       const float *const A, const float *const B, float *const C, 
                       const int Astride, const int Bstride, const int Cstride);

extern "C" void mul_mdmd_md(const int M, const int K, const int N, 
                       const double *const A, const double *const B, double *const C, 
                       const int Astride, const int Bstride, const int Cstride);

Data32 const *
AffineFilter::localTransform(Data32 const *inputSubMatrix, 
			  subMatrixDescriptor const &inputDescription,
			  subMatrixDescriptor *outputDescription) 
{

  // FIXME - error checking that B and c are compatible with X
  if (cols != inputDescription.numContinuous) {
    error("AffineFilter: filter expects %u features, but input has %u\n", 
	  cols, inputDescription.numContinuous);
  }

  if (B) {
    unsigned workNeeded = inputDescription.numFrames * cols;
    if (workBuffSize < workNeeded) {
      workBuffer = (double *)realloc(workBuffer, workNeeded * sizeof(double));
      assert(workBuffer);
      workBuffSize = workNeeded;
    }
    
    unsigned xNeeded = inputDescription.numFrames * inputDescription.numContinuous;
    if (xBuffSize < xNeeded) {
      xBuffer = (double *)realloc(xBuffer, xNeeded * sizeof(double));
      assert(xBuffer);
      xBuffSize = xNeeded;
    }

    // FIXME - PHiPAC optimize
    // make it a double
    float *xSource = (float *) inputSubMatrix;
    double *xDest = xBuffer;
    for (unsigned i=0; i < inputDescription.numFrames; i+=1) {
      for (unsigned j=0; j < inputDescription.numContinuous; j+=1) {
	*(xDest++) = (double) xSource[j];
      }
      xSource += inputDescription.numContinuous+inputDescription.numDiscrete;
    }

    // Y = XB
    mul_mdmd_md(inputDescription.numFrames,      // X rows
		inputDescription.numContinuous,  // X columns & B rows
		cols,                            // B columns
		xBuffer, B, workBuffer,          // X, B, Y
		inputDescription.numContinuous,  // X stride = # X columns
		cols,                            // B stride = # B columns
		cols);                           // Y stride = # Y columns
  } else {
    // FIXME - PHiPAC optimize
    unsigned workNeeded = inputDescription.numFrames * inputDescription.numContinuous;
    if (workBuffSize < workNeeded) {
      workBuffer = (double *)realloc(workBuffer, workNeeded * sizeof(double));
      assert(workBuffer);
      workBuffSize = workNeeded;
    }
    
    float *xSource = (float *) inputSubMatrix;
    double *xDest = workBuffer;
    for (unsigned i=0; i < inputDescription.numFrames; i+=1) {
      for (unsigned j=0; j < inputDescription.numContinuous; j+=1) {
	*(xDest++) = (double) xSource[j];
      }
      xSource += inputDescription.numContinuous+inputDescription.numDiscrete;
    }
  }

  subMatrixDescriptor myOutput = describeLocalOutput(inputDescription);
  assert(myOutput.numContinuous == cols);
  assert(myOutput.numDiscrete == inputDescription.numDiscrete);

  unsigned stride = myOutput.numContinuous + myOutput.numDiscrete;
  unsigned needed = myOutput.numFrames * stride;
  if (needed > buffSize) {
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    assert(buffer);
    buffSize = needed;
  }

  // FIXME - PHiPAC optimize
  for (unsigned i=0; i < myOutput.numFrames; i+=1) {
    float *outputCont = (float *)(buffer + i * stride);
    if (c) {
      for (unsigned j=0; j < myOutput.numContinuous; j+=1) {
	outputCont[j] = (float)(workBuffer[i * myOutput.numContinuous + j] + c[j]);
      }
    } else {
      for (unsigned j=0; j < myOutput.numContinuous; j+=1) {
	outputCont[j] = (float)(workBuffer[i * myOutput.numContinuous + j]);
      }
    }
    Uint32 *outputDisc = (Uint32 *)(buffer + i * stride + myOutput.numContinuous);
    Uint32 *inputDisc  = 
      (Uint32 *)(inputSubMatrix + 
		 i * (inputDescription.numContinuous + inputDescription.numDiscrete) +
		 inputDescription.numContinuous);
    memcpy(outputDisc, inputDisc, myOutput.numDiscrete * sizeof(Data32));
  }


  if (outputDescription) *outputDescription = myOutput;
  return buffer;
}
