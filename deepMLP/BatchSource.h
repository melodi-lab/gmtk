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
class BatchSource {

 protected:

  unsigned batchSize; // number of training instances in a batch

 public:

  BatchSource(unsigned batchSize) : batchSize(batchSize) {}

  // Return the ith [mini]batch. Note that there's no guarantee that
  // you will get the same results from multiple calls for a particular i.
  // Note that the returned batch may be shorter than the batchSize.
  // Matrix is column-major, with a column corresponding to a training instance.
  virtual Matrix getBatch(unsigned i) = 0;

  // Return the label Matrix for the ith [mini]batch. The result of 
  // getLabels(i) will correspond to the result of getBatch(i) as long
  // as there aren't any intervening calls of getLabels(j) or getBatch(j) 
  // where j != i.
  virtual Matrix getLabels(unsigned i) = 0;

  // Smallest number of [mini]batches such that the total number of [mini]batch
  // columns will be >= the total number of training instances.
  virtual unsigned batchesPerEpoch() = 0;

  // So you don't have to actually produce a Matrix to know the # of rows...
  // You can't know the number of columns in advance, since returned [mini]batches
  // may be smaller than batchSize.
  virtual unsigned numBatchRows() = 0;

  virtual unsigned numLabelRows() = 0;
};


// Stream batches from instances of Galen's Matrix class
class MatrixBatchSource : public BatchSource {

  Matrix const     &batchSource;  // training data
  Matrix const     &labelSource;  // labels thereof

  bool              small;        // mini-batch size <= source matrix size

  AllocatingMatrix  batchTemp;    // used to assemble "wrap-around" mini-batches
  AllocatingMatrix  labelTemp;
  

  // get the ith batch or labels
  Matrix getMatrix(unsigned i, Matrix const &source, AllocatingMatrix &temp) {
    if (small) return source; // that's all there is...

    Matrix batch;
    uint64_t c = (i * batchSize) % source.NumC();
    unsigned startCol = (unsigned) c;
    unsigned endCol = startCol + batchSize;
    if (endCol <= source.NumC()) {
      batch = source.GetCols(startCol, endCol);
      return batch;
    }

    // take the last few columns, then wrap around to the first columns
    int numLastCols = source.NumC() - startCol;
    temp.GetCols(0, numLastCols).CopyFrom(source.GetCols(startCol, -1)); // last few
    endCol %= source.NumC();
    temp.GetCols(numLastCols, -1).CopyFrom(source.GetCols(0, endCol)); // wrap around
    batch = temp;
    return batch;
  }

 public:

  MatrixBatchSource(unsigned batchSize, Matrix const &batchSource, Matrix const &labelSource) 
    : BatchSource(batchSize), batchSource(batchSource), labelSource(labelSource),
      small(batchSource.NumC() <= batchSize)
  {
    assert(batchSource.NumC() == labelSource.NumC());
    batchTemp.Resize(batchSource.NumR(), batchSize);
    labelTemp.Resize(labelSource.NumR(), batchSize);
  }


  Matrix getBatch(unsigned i) {
    return getMatrix(i, batchSource, batchTemp);
  }


  Matrix getLabels(unsigned i) {
    return getMatrix(i, labelSource, labelTemp);
  }


  unsigned batchesPerEpoch() { 
    return (batchSource.NumC() / batchSize) + ((batchSource.NumC() % batchSize) ? 1 : 0);
  }


  unsigned numBatchRows() { return batchSource.NumR(); }

  virtual unsigned numLabelRows() { return labelSource.NumR(); }

};
