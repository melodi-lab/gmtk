
#include <stdio.h>

#include <fstream>
#include <iostream>
using namespace std;

unsigned int SwapEndian(unsigned int x) {
  char * b = (char*) &x;
  unsigned int ret;
  char* br = (char*) &ret;
  br[0] = b[3];
  br[1] = b[2];
  br[2] = b[1];
  br[3] = b[0];
  return ret;
}

// mnist2lbl labels obsfile

int
main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "mnist2lbl labels obsfile\n";
    return 1;
  }

  ifstream labelStream(argv[1], ios::in|ios::binary);
  if (!labelStream.is_open()) {
    cout << "Couldn't open label file " << argv[1] << endl;
    return 1;
  }
  
  unsigned int header[2];
  labelStream.read((char*)header, 8);

  unsigned magic     = SwapEndian(header[0]);
  unsigned numLabels = SwapEndian(header[1]);
  unsigned numRows   = SwapEndian(header[2]);
  unsigned numCols   = SwapEndian(header[3]);

  if (magic != 2049) {
    cerr << "magic number " << magic << " != 2049\n";
    return 1;
  }

  FILE *f = fopen(argv[2], "w");
  if (!f) {
    perror(argv[2]);
    return 1;
  }

  unsigned char* lbl = new unsigned char[numLabels];
  labelStream.read((char*)lbl, numLabels);

  for (unsigned i=0; i < numLabels; i+=1) {
    fprintf(f, "0 %05u %u\n", i, lbl[i]);
  }
  if (fclose(f)) {
    perror(argv[2]);
    return 1;
  }

  return 0;
}
