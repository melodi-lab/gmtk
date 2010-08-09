#include "global-parameters.h"
#include "general.h"
#include "error.h"
#include "mixNormal.h"
#include "mixNormalCollection.h"
#include "readRange.h"
#include "GMTK_ObservationMatrix.h"

void dumpDistribSampleData(FILE* ofp, 
			   ObservationMatrix * obsMat,
			   RangeSetCollection &tupleCol,
			   Range &lrrng,
			   Range &sentRange,
			   unsigned numMixtures,
			   unsigned maxIter,
			   int labpos,
			   unsigned tupleNum,
			   const bool quiet); 


