#pragma once

/*-
 * AsynchronousBatchSource.h
 *     A threaded mini-batch stream source adaptor. This
 *     class will spawn a producer thread to fill a queue 
 *     of mini-batches produced by another BatchSource 
 *     instance.
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
#if HAVE_PTHREAD
#  include <pthread.h>
#endif

#include <assert.h>

#include "debug.h"

#include "BatchSource.h"

#include "Matrix.h"


void *MinibatchProducer(void *asynchBatchSource);


class AsynchronousBatchSource : public BatchSource {

  BatchSource      *src;         // BatchSource run by the producer thread
  AllocatingMatrix *dataQueue;   // queue of training features
  AllocatingMatrix *labelQueue;  // queue of labels
  unsigned          queueSize;   // queue capacity
  unsigned          queueCount;  // # instances in the queues
  unsigned          queueStart;  // # index of head of queue

  // temporary batch feature/label storage
  AllocatingMatrix batchData;
  AllocatingMatrix batchLabels;

#if HAVE_PTHREAD
  pthread_t         producerThread;
  pthread_mutex_t   queueMutex;
  pthread_cond_t    queueNotFull;  // producerThread waits on this
  pthread_cond_t    queueNotEmpty; // consumer waits on this
#endif

 public:

  AsynchronousBatchSource(BatchSource *src, unsigned queueSize) 
    : src(src), queueSize(queueSize), queueCount(0), queueStart(0)
#if HAVE_PTHREAD
      , queueNotFull(PTHREAD_COND_INITIALIZER),
      queueNotEmpty(PTHREAD_COND_INITIALIZER)
#endif
  { 
    assert(src);
    dataQueue = new AllocatingMatrix[queueSize];
    labelQueue = new AllocatingMatrix[queueSize];
    for (unsigned i=0; i < queueSize; i+=1) {
       dataQueue[i].Resize(src->numDataRows(), 1);
      labelQueue[i].Resize(src->numLabelRows(), 1);
    }
#if HAVE_PTHREAD
    pthread_mutex_init(&queueMutex, NULL);
    if (pthread_create(&producerThread, NULL, MinibatchProducer, (void *)this)) {
      error("ERROR: unable to start asynchronous minibatch producer thread\n");
    }
#endif
  }

  ~AsynchronousBatchSource() {
#if HAVE_PTHREAD
    assert(pthread_cancel(producerThread) == 0);
#endif
    if (dataQueue)  delete[] dataQueue;
    if (labelQueue) delete[] labelQueue;
  }

  // Call getBatch to keep the two queues aligned
  Matrix getData(unsigned batchSize) { 
    Matrix data, labels;
    getBatch(batchSize, data, labels);
    return data;
  }

  // Call getBatch to keep the two queues aligned
  Matrix getLabels(unsigned batchSize) {
    Matrix data, labels;
    getBatch(batchSize, data, labels);
    return labels;
  }
  
  // consume a batch from the queue
  void getBatch(unsigned batchSize, Matrix &data, Matrix &labels) {
#if HAVE_PTHREAD
    if (batchSize >= queueSize) {
      error("ERROR: minibatch of size %u requested, but the batch queue can only hold %u\n", 
	    batchSize, queueSize);
    }
    pthread_mutex_lock(&queueMutex); // lock queue, wait for queue to contain at least batchSize entries
    while (queueCount < batchSize) {
      infoMsg(IM::ObsFile, IM::Mod, "AsynchBatchSource: waiting for %u batches to be enqueued\n", 
	      batchSize - queueCount);
      pthread_cond_wait(&queueNotEmpty, &queueMutex);
    }
      batchData.Resize(src->numDataRows(),  batchSize);
    batchLabels.Resize(src->numLabelRows(), batchSize);
    // pull the batch out of the queue
    for (unsigned n=0; n < batchSize; n+=1, queueCount-=1, queueStart = (queueStart+1) % queueSize) {
      unsigned numC = dataQueue[queueStart].NumC();
      infoMsg(IM::ObsFile, IM::Mod, "AsynchBatchSource: dequeuing %u from %u\n", numC, queueStart);
        batchData.GetCols(n, n + numC).CopyFrom(dataQueue[queueStart]);
      batchLabels.GetCols(n, n + numC).CopyFrom(labelQueue[queueStart]);
    }
    pthread_cond_signal(&queueNotFull); // now there's room for the producer to fill, unlock
    pthread_mutex_unlock(&queueMutex);
    data = batchData;
    labels = batchLabels;
#else
    src->getBatch(batchSize, data, labels);
#endif
  }


  // Ye Olde Producer
  void fill() {
#if HAVE_PTHREAD
    do {
      pthread_mutex_lock(&queueMutex); // lock queue & wait for room to produce
      while (queueCount == queueSize) {
	infoMsg(IM::ObsFile, IM::Mod, "AsynchBatchSource: waiting for queue to not be full\n");
	pthread_cond_wait(&queueNotFull, &queueMutex);
      }
      unsigned numSpaces = queueSize - queueCount; // try to fill the whole queue 
      infoMsg(IM::ObsFile, IM::Mod, "AsynchBatchSource: enqueuing %u -> %u .. %u\n", numSpaces, 
	      (queueStart+queueCount) % queueSize, 
	      (queueStart+queueCount + numSpaces) % queueSize);
      Matrix data, labels;
      src->getBatch(numSpaces, data, labels);
      for (int i=0, queueEnd = (queueStart + queueCount) % queueSize; 
	   i < data.NumC(); 
	   i+=1, queueCount+=1, queueEnd = (queueEnd + 1) % queueSize )
      {
#if 0
	// paranoia
	assert(data.Start());
	assert(labels.Start());
	assert(data.NumC() == labels.NumC());
	assert(i < data.NumC());
	assert(data.GetCols(i,i+1).NumC() > 0);
	assert(data.GetCols(i,i+1).Start());
#endif
	// enqueue the data
         dataQueue[queueEnd].CopyFrom(   data.GetCols(i, i+1) );
	labelQueue[queueEnd].CopyFrom( labels.GetCols(i, i+1) );
	infoMsg(IM::ObsFile, IM::Mod, "AsynchBatchSource: enqueuing %u\n", queueEnd);
      }
      pthread_cond_signal(&queueNotEmpty); // now there's data available for the consumer, unlock
      pthread_mutex_unlock(&queueMutex);
    } while (1);
#endif
  }


  unsigned batchesPerEpoch(unsigned batchSize) {
    return src->batchesPerEpoch(batchSize);
  }

  unsigned epochSize() { return src->batchesPerEpoch(1); }

  unsigned numDataRows() { return src->numDataRows(); }

  unsigned numLabelRows() { return src->numLabelRows(); }
};
