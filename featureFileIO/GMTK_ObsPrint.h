#ifndef GMTK_OBSPRINT_H
#define GMTK_OBSPRINT_H

#include <cstdlib>
#include "range.h"
#include "GMTK_ObservationMatrix.h"

void printSegment(unsigned sent_no, FILE* out_fp, float* cont_buf, unsigned num_continuous, UInt32* disc_buf, unsigned num_discrete, unsigned num_frames, const bool dontPrintFrameID,const bool quiet,unsigned ofmt,int debug_level,bool oswap, OutFtrLabStream_PFile* out_stream);

void printHTKHeader(FILE* ofp, bool oswap, int numInts, int numFloats, int numSamples, short parameterKInd, Int32 samplePeriod);

void obsPrint(FILE* out_fp,Range& srrng,const char * pr_str,const bool dontPrintFrameID,const bool quiet,unsigned ofmt,int debug_level,bool oswap);

#endif
