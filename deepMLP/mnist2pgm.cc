
#include <stdio.h>

#include <fstream>
#include <iostream>
using namespace std;

const int MAX_PIXEL = 255;

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

// mnist2pgm images prefix

int
main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "mnist2pgm images prefix\n";
    return 1;
  }

  ifstream imageStream(argv[1], ios::in|ios::binary);
  if (!imageStream.is_open()) {
    cout << "Couldn't open feature file " << argv[1] << endl;
    return 1;
  }
  
  unsigned int header[4];
  imageStream.read((char*)header, 16);

  unsigned magic     = SwapEndian(header[0]);
  unsigned numImages = SwapEndian(header[1]);
  unsigned numRows   = SwapEndian(header[2]);
  unsigned numCols   = SwapEndian(header[3]);

  if (magic != 2051) {
    cerr << "magic number " << magic << " != 2051\n";
    return 1;
  }
  cout << numImages << " images of size " << numRows << "x" << numCols << endl;

  for (unsigned i=0; i < numImages; i+=1) {
    string prefix(argv[2]);
    char buf[16];
    sprintf(buf, "%05u", i);
    string fname = prefix.append(buf);
    FILE *f = fopen(fname.c_str(), "w");
    if (!f) {
      perror(fname.c_str());
      return 1;
    }
    fprintf(f, "P2\n# Image %u from %s\n", i, argv[1]);
    fprintf(f, "%u %u\n%u\n", numCols, numRows, MAX_PIXEL);
    for (unsigned r=0; r < numRows; r+=1) {
      unsigned char row[numCols];
      imageStream.read((char*)row, numCols);
      for (unsigned c=0; c < numCols; c+=1) {
	fprintf(f, "%u ", MAX_PIXEL - row[c]);
      }
      fprintf(f, "\n");
    }
    if (fclose(f)) {
      perror(fname.c_str());
      return 1;
    }
  }

  return 0;
}
