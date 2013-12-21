
/*-
 * AsynchronousBatchSource.cc
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


#include "AsynchronousBatchSource.h"

void *
MinibatchProducer(void *asynchBatchSource) {
  AsynchronousBatchSource *bs = (AsynchronousBatchSource *) asynchBatchSource;
  infoMsg(IM::ObsFile, IM::Mod, "launching producer thread\n");
  bs->fill();
#if HAVE_PTHREAD
  pthread_exit(NULL);
#else
  return NULL;
#endif
}
