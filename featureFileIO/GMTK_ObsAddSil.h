#ifndef GMTK_OBSADDSIL_H
#define GMTK_OBSADDSIL_H

/*
 * Copyright (C) 2004 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 */

#include "GMTK_ObsPrint.h"

void addSil(FILE* out_fp, 
	    HDF5File *hdf5,
	     FileSource* obs_mat,
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
