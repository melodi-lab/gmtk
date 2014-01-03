
#include "FileBackedMatrix.h"

map<string,unsigned> FileBackedMatrix::ref_count;

char const *FileBackedMatrix::dmlpTempDir = "/tmp";
unsigned FileBackedMatrix::fileNumber = 0;

