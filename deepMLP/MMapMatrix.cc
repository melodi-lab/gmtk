
#include "MMapMatrix.h"

map<string,unsigned> MMapMatrix::ref_count;

char const *MMapMatrix::dmlpTempDir = "/tmp";
unsigned MMapMatrix::fileNumber = 0;

