#ifndef GMTK_OBSKLT_H
#define GMTK_OBSKLT_H

/*
 * Copyright (C) 2004 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 */

#include "GMTK_ObsPrint.h"

void readStats(FILE*f, size_t N, bool ascii, double *cor, double *means, double *vecs, double *vals);
void writeStats(FILE*f, size_t N, bool ascii, double *cor, double *means, double *vecs, double *vals);
void obsKLT(FILE* out_fp, HDF5File *hdf5, FileSource* obs_mat, FILE *in_st_fp,FILE *out_st_fp, Range& ofrrng,const bool unity_variance,const bool ascii,const bool dontPrintFrameID,const bool quiet,unsigned ofmt,int debug_level,bool oswap);

#endif
