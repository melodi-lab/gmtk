#ifndef GMTK_OBSKLT_H
#define GMTK_OBSKLT_H

#include "GMTK_ObsPrint.h"

void readStats(FILE*f, size_t N, bool ascii, double *cor, double *means, double *vecs, double *vals);
void writeStats(FILE*f, size_t N, bool ascii, double *cor, double *means, double *vecs, double *vals);
void obsKLT(FILE* out_fp, ObservationMatrix* obs_mat, FILE *in_st_fp,FILE *out_st_fp, Range& srrng,Range& ofrrng,const bool unity_variance,const bool ascii,const bool dontPrintFrameID,const bool quiet,unsigned ofmt,int debug_level,bool oswap);

#endif
