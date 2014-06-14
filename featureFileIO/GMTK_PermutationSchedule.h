/*
 * GMTK_PermutationSchedule.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_PERMUTATIONSCHEDULE_H
#define GMTK_PERMUTATIONSCHEDULE_H

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

#include <string.h>

#include "GMTK_TrainingSchedule.h"
#include "GMTK_ObservationSource.h"
#include "rand.h"
#include "prime.h"
#include "error.h"
#include "debug.h"


// Return non-overlapping training units of requested size according to a permutation
// of the $T$ total training units determined by $\sigma(i) = (ai + b)^3 \bmod p$ where
//    $p$ is the smallest prime such that $p >= T$ and $p \equiv_3 2$
//    $a$ is a random integer in $[1, p)$
//    $b$ is a random integer in $[0, p)$


class PermutationSchedule : public TrainingSchedule {

  vector<unsigned> segmentScanSum;  // total # of training units preceeding segment i
  vector<unsigned> segmentIndex;    // segmentIndex[i] is the ith segment # in the trrng

  uint32_t  i;               // position in permutation

  uint32_t  a, b, p;         // p is the smallest prime \equiv_3 2 larger than num_units
                             // a is an integer in [1, p)  and b is an integer in [0, p)
                             // \sigma(i) = (ai+b)^3 mod p is a permutation of 0, 1, ..., p-1
 public:

  PermutationSchedule(unsigned feature_offset, unsigned features_per_frame,
		 unsigned label_offset,  unsigned label_domain_size,
		 bool one_hot, unsigned window_radius, unsigned unit_size, 
		 FileSource *obs_source, char const *trrng_str)
    : TrainingSchedule(feature_offset, features_per_frame, label_offset, 
		       label_domain_size, one_hot, window_radius, unit_size,
		       obs_source, trrng_str), i(0)
  {
    // count viable training units
    segmentScanSum.resize(trrng->length()+1, 0); // +1 to avoid writing past end
    segmentIndex.resize(trrng->length());
    num_viable_units = 0;
    Range::iterator* trrng_it = new Range::iterator(trrng->begin());
    unsigned ii = 0;
    while (!trrng_it->at_end()) {
      unsigned i = (unsigned)(*(*trrng_it));
      if (!obs_source->openSegment(i)) {
	error("ERROR: Unable to open observation file segment %u\n", i);
      }
      unsigned num_frames = obs_source->numFrames();
      unsigned units = num_frames / unit_size;
      if (num_frames % unit_size) units += 1;
      if (units == 0) {
	warning("WARNING: segment %u contains no frames\n", i);
      }
      num_viable_units += units;
      segmentIndex[ii] = i;
      segmentScanSum[ii+1] = segmentScanSum[ii] + units; // ii+1 as frames in this segment preceed the next segment
      ii += 1;
      (*trrng_it)++;
    }
    if (num_viable_units == 0) {
      error("ERROR: observation files contain no viable training instances\n");
    }

    // find smallest prime p >= num_viable_units that is \equiv_3 2
    for (p = num_viable_units + (2 - num_viable_units % 3); !prime32(p); p += 3)
      ;
    a = rnd.uniform(1, p-1);     // pick a random a in [1,p)
    b = rnd.uniform(p-1);        // pick a random b in [0,p)
    infoMsg(IM::ObsFile, IM::Moderate, "T=%u <= %u = %u mod 3;  sigma(i) = (%u i + %u)^3 mod %u\n",
	    num_viable_units, p, p%3, a, b, p);
  } 

    
  void nextTrainingUnit(unsigned &segment, unsigned &frame) { 
    uint32_t sigma; // sigma(i) = (ai + b)^3 mod p
    do {
      uint64_t t = (a * i) % p;
      t = (t + b) % p;
      uint64_t tt = (t * t) % p;
      sigma = (t * tt) % p;
      i = (i + 1) % p;
    } while (sigma >= num_viable_units); // skip any extras since p >= num_viable_units

    // binary search for the segment containing the selected training unit
    unsigned l=1, m, r=segmentScanSum.size();
    do {
      m = (l+r)/2;
      if (segmentScanSum[m] <= sigma) {
	l = m+1;
      } else {
	r = m;
      }
    } while ( ! (segmentScanSum[m-1] <= sigma && sigma < segmentScanSum[m]) ); // m >= 1

    segment = segmentIndex[m-1];
    frame = sigma - segmentScanSum[m-1];
    TrainingSchedule::nextTrainingUnit(segment, frame);
  }
};

#endif
