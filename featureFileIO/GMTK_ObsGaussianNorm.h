#ifndef GMTK_OBSGAUSSIANNORM_H
#define GMTK_OBSGAUSSIANNORM_H

/*
 * Copyright (C) 2004 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 */

#include "GMTK_ObsPrint.h"

void gaussianNorm(FILE* out_fp,
		  FileSource* obs_mat,
		  FILE *in_st_fp,
		  FILE *out_st_fp,
		  Range& srrng,
		  Range& frrng,
		  const char*pr_str,
		  const size_t hist_bins, 
		  const float num_stds,
		  const bool uniform_output,
		  const bool dontPrintFrameID,const bool quiet,unsigned ofmt,int debug_level,bool oswap);

#endif
