#ifndef GMTK_OBSGAUSSIANNORM_H
#define GMTK_OBSGAUSSIANNORM_H

#include "GMTK_ObsPrint.h"

void gaussianNorm(FILE* out_fp,
		  ObservationMatrix* obs_mat,
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
