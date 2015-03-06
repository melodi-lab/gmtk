#pragma once

/*-
 * BatchSource.h
 *     Mini-batch stream source
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_INTTYPES_H
   // The ISO C99 standard specifies that the macros in inttypes.h must
   //  only be defined if explicitly requested. 
#  ifndef __STDC_FORMAT_MACROS
#    define __STDC_FORMAT_MACROS 1
#  endif
#  include <inttypes.h>
#endif
#if HAVE_STDINT_H
#  include <stdint.h>
#endif

#include <assert.h>

#include "error.h"
#include "debug.h"
#include "GMTK_TrainingSchedule.h"

#include "Matrix.h"


// BatchSource subclasses provide a stream of [mini]batches for NN trainers.
// MatrixBatchSource produces the [mini]batch stream from instances of Galen's 
//   Matrix class (see Matrix.h). 
// The ScheduleBatchSource produces the stream from observation files permuted 
//   according to a TrainingSchedule (see featureFileIO/TrainingSchedule.h).

class BatchSource {

 public:

  // Get the next batchSize training instances' data. You won't be able to get
  // their labels -- use getBatch() if you need them. May return fewer instances
  // than you request. Instances are columns, Matrix in column-major order.
  virtual Matrix getData(unsigned batchSize) = 0;

  // Get the next batchSize labels. You won't be able to get the input data
  // corresponding to the labels -- use getBatch if you need it...
  virtual Matrix getLabels(unsigned batchSize) = 0;
  
  // Get the next [mini]batch's training data and labels.
  virtual void getBatch(unsigned batchSize, Matrix &data, Matrix &labels) = 0;

  // Smallest number of [mini]batches such that the total number of [mini]batch
  // columns will be >= the total number of training instances.
  virtual unsigned batchesPerEpoch(unsigned batchSize) = 0;

  // Number of instances in an epoch
  virtual unsigned epochSize() { return batchesPerEpoch(1); }

  // So you don't have to actually produce a Matrix to know the # of rows...
  // You can't know the number of columns in advance, since returned [mini]batches
  // may be smaller than batchSize.
  virtual unsigned numDataRows() = 0;

  virtual unsigned numLabelRows() = 0;
};


// Stream batches from instances of Galen's Matrix class
class MatrixBatchSource : public BatchSource {

  Matrix const     &dataSource;   // training data
  Matrix const     &labelSource;  // labels thereof

  unsigned nextColumn;            // index of start of next batch

  AllocatingMatrix  dataTemp;     // used to assemble "wrap-around" mini-batches
  AllocatingMatrix  labelTemp;
  

  // get the batch or labels starting at column i
  Matrix getMatrix(unsigned i, unsigned batchSize, Matrix const &source, AllocatingMatrix &temp, unsigned &actualSize) {
    temp.Resize(source.NumR(), batchSize);
    unsigned srcCols = source.NumC();

    Matrix batch;
    uint64_t c = i % srcCols;
    unsigned startCol = (unsigned) c;
    unsigned endCol = startCol + batchSize;
    if (endCol <= srcCols) {
      batch = source.GetCols(startCol, endCol);
      actualSize = batchSize;
      return batch;
    }

    // take the last few columns, then wrap around to the first columns
    int numLastCols = srcCols - startCol;
    temp.GetCols(0, numLastCols).CopyFrom(source.GetCols(startCol, -1)); // last few
    if (batchSize <= srcCols) { // batch is smaller than matrix
      endCol %= srcCols;
      actualSize = batchSize;
    } else {
      endCol = startCol; // use entire matrix
      actualSize = srcCols;
    }
    temp.GetCols(numLastCols, -1).CopyFrom(source.GetCols(0, endCol)); // wrap around
    batch = temp;
    return batch;
  }

 public:

  MatrixBatchSource(Matrix const &dataSource, Matrix const &labelSource) 
    : dataSource(dataSource), labelSource(labelSource), nextColumn(0)
  {
    assert(dataSource.NumC() == labelSource.NumC());
  }


  // Use the same matrix as training data & labels (e.g., for pre-training)
  MatrixBatchSource(Matrix const &source) 
    : dataSource(source), labelSource(source), nextColumn(0)
  { }


  Matrix getData(unsigned batchSize) {
    unsigned actualSize;
    Matrix batch = getMatrix(nextColumn, batchSize, dataSource, dataTemp, actualSize);
    nextColumn += actualSize;
    nextColumn %= dataSource.NumC();
    return batch;
  }


  Matrix getLabels(unsigned batchSize) {
    unsigned actualSize;
    Matrix batch = getMatrix(nextColumn, batchSize, labelSource, labelTemp, actualSize);
    nextColumn += actualSize;
    nextColumn %= dataSource.NumC();
    return batch;
  }


  void getBatch(unsigned batchSize, Matrix &data, Matrix &labels) {
    unsigned actualBatchSize, actualLabelSize;
    data   = getMatrix(nextColumn, batchSize, dataSource,  dataTemp, actualBatchSize);
    labels = getMatrix(nextColumn, batchSize, labelSource, labelTemp, actualLabelSize);
    assert(actualBatchSize == actualLabelSize);
    nextColumn += actualBatchSize;
    nextColumn %= dataSource.NumC();
  }


  unsigned batchesPerEpoch(unsigned batchSize) { 
    return (dataSource.NumC() / batchSize) + ((dataSource.NumC() % batchSize) ? 1 : 0);
  }


  unsigned numDataRows() { return dataSource.NumR(); }

  virtual unsigned numLabelRows() { return labelSource.NumR(); }

};


// Stream batches from a TrainingSchedule
class ScheduleBatchSource : public BatchSource {
  TrainingSchedule *schedule;
  float            *TSfeatures;   // current unit's instance data
  float            *TSlabels;     // current unit's label data
  unsigned          unitSize;     // # instances in current training unit
  unsigned          curInstance;  // index of current instance within current unit

  AllocatingMatrix  dataTemp;     // used to assemble "wrap-around" mini-batches
  AllocatingMatrix  labelTemp;

  unsigned dataRows, labelRows, labelStride;

 public:

  ScheduleBatchSource(TrainingSchedule *schedule) 
    : schedule(schedule), TSfeatures(NULL), TSlabels(NULL), unitSize(0), curInstance(0)
  {
    assert(schedule);
    unsigned dummy;
    schedule->describeFeatures(dataRows, dummy);
    schedule->describeLabels(labelRows, dummy, labelStride);
    infoMsg(IM::ObsFile, IM::Moderate, "Training schedule batch source: data %u x %u;  labels %u x %u + %u\n",
	    dataRows, epochSize(), labelRows, epochSize(), labelStride);
  }


  Matrix getData(unsigned batchSize) {
    Matrix result, dummy;
    getBatch(batchSize, result, dummy);
    return result;
  }


  Matrix getLabels(unsigned batchSize) {
    Matrix result, dummy;
    getBatch(batchSize, dummy, result);
    return result;
  }


  void getBatch(unsigned batchSize, Matrix &data, Matrix &labels) {
    infoMsg(IM::ObsFile, IM::High, "ScheduleBatchSource::getBatch(%u)\n", batchSize);
    dataTemp.Resize(dataRows, batchSize);
    labelTemp.Resize(labelRows, batchSize);
    double *dataP = dataTemp.Start();
    double *labelsP = labelTemp.Start();
    // copy columns from {data,label}Temp until they're all used, then load more according to the schedule
    for (unsigned i = 0; i < batchSize; i += 1, 
	   curInstance += 1, dataP += dataRows, labelsP += labelRows, TSfeatures += dataRows, TSlabels += labelStride) 
    {
      if (curInstance == unitSize) { // need to load a new training unit
	unsigned segment, frame;
	schedule->nextTrainingUnit(segment, frame);
	TSfeatures = schedule->getFeatures(segment, frame, unitSize);
	TSlabels = schedule->getLabels(segment, frame, unitSize);
	infoMsg(IM::ObsFile, IM::High, "ScheduleBatchSource: loading training unit %u instances @ segment %u frame %u\n",
		unitSize, segment, frame);
	curInstance = 0;
      }
      for (unsigned r=0; r < dataRows; r+=1)
	dataP[r] = (double) (TSfeatures[r]);
      for (unsigned r=0; r < labelRows; r+=1)
	labelsP[r] = (double) (TSlabels[r]);
    }
    data = dataTemp;
    labels = labelTemp;
  }


  unsigned batchesPerEpoch(unsigned batchSize) {
    return (schedule->numInstances() / batchSize) + 
           ((schedule->numInstances() % batchSize) ? 1 : 0);
  }


  unsigned numDataRows() {
    return dataRows;
  }


  unsigned numLabelRows() {
    return labelRows;
  }
};
