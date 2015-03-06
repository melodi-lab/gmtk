#ifndef GMTK_OBSNORM_H
#define GMTK_OBSNORM_H

/*
 * Copyright (C) 2004 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 */

#include "GMTK_ObsPrint.h"


void getSegMarkers(int* seg_markers, const unsigned size, char* segment_length_fname, unsigned segment_length);

void obsNorm(FILE* out_fp,
	     HDF5File *hdf5,
	     FileSource* obs_mat,
	     Range& srrng,
	     const double result_mean,
	     const double result_std,
	     char* segment_length_fname,
	     unsigned   segment_length,
	     const bool dontPrintFrameID,
	     const bool quiet,
	     unsigned ofmt,
	     int debug_level,
	     bool oswap);

#endif
