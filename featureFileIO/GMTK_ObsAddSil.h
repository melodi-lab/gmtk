#ifndef GMTK_OBSADDSIL_H
#define GMTK_OBSADDSIL_H

#include "GMTK_ObsPrint.h"

void addSil(FILE* out_fp, 
	     ObservationMatrix* obs_mat,
	     Range& srrng,
	     const int nb,
	     const char *prb_str,
	     const int ne,
	     const char *pre_str,
	     const double mmf,
	     const double maf,
	     const double smf,
	     const double saf,
	     const bool dontPrintFrameID,
	     const bool quiet,
	     unsigned ofmt,
	     int debug_level,
	    bool oswap);


#endif
