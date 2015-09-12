/*
 * GMTK_Separator.h
 *   
 * Eventially this may be a full junction forest representation
 * of the separator. Initially, it may be a single separator
 * clique as in the current inference implementation. 
 *
 * SectionInferenceAlgorithm classes must know how to send/recieve
 * the necessary messages to/from a Separator.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2014 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_SEPARAOTR_H
#define GMTK_SEPARAOTR_H

class Separator {

  void collectEvidence(bool doViterbi=false, unsigned nBest=1);

  void distributeEvidence(bool doViterbi=false, unsigned nBest=1);
};

#endif
