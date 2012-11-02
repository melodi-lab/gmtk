#include "global-parameters.h"
#include "general.h"
#include "error.h"
#include "mixNormal.h"
#include "mixNormalCollection.h"
#include "readRange.h"
#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_FileSource.h"
#endif

void dumpDistribSampleData(FILE* ofp, 
			   FileSource *obsMat,
			   RangeSetCollection &tupleCol,
			   Range &lrrng,
			   Range &sentRange,
			   unsigned numMixtures,
			   unsigned maxIter,
			   int labpos,
			   unsigned tupleNum,
			   const bool quiet); 


