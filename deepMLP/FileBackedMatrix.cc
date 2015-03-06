
#include "FileBackedMatrix.h"

map<string,unsigned> FileBackedMatrix::ref_count;

char const *FileBackedMatrix::dmlpTempDir = NULL;
unsigned FileBackedMatrix::fileNumber = 0;

