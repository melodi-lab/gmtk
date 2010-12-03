
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "pfile.h"
#include "GMTK_WordOrganization.h"

#define DEBUG 0

int
readFrame(OutFtrLabStream_PFile *pfile, char *fname, FILE *in, bool &eof, 
	  int &curSent, int numFloat, int numInt, float *floatBuf, int *intBuf) 
{
  int err = 0;
  int nread = fscanf(in, "%d", &curSent);
  eof = nread == EOF;
  if (eof) goto done;
  if (nread != 1) {
    err = errno;
    perror(fname);
    goto done;
  }
  
  int curFrame;
  nread = fscanf(in, "%d", &curFrame);
  eof = nread == EOF;
  if (eof) goto done;
  if (nread != EOF && nread != 1) {
    err = errno;
    perror(fname);
    goto done;
  }
  
#if DEBUG
  fprintf(stderr, "%3d %3d:", curSent, curFrame);
#endif
  for (int j=0; !eof && j < numFloat; j+=1) {
    int nread = fscanf(in, "%f", &(floatBuf[j]));
    eof = nread == EOF;
    if (nread != EOF && nread != 1) {
      err = errno;
      perror(fname);
      goto done;
    }
#if DEBUG
    fprintf(stderr, " %f", floatBuf[j]);
#endif
  }
  for (int j=0; !eof && j < numInt; j+=1) {
    int nread = fscanf(in, "%d", &(intBuf[j]));
    eof = nread == EOF;
    if (nread != EOF && nread != 1) {
      err = errno;
      perror(fname);
      goto done;
    }
#if DEBUG
    fprintf(stderr, " %d", intBuf[j]);
#endif
  }
  pfile->write_ftrslabs(1, floatBuf, (const UInt32*) intBuf);
#if DEBUG
  fprintf(stderr, "\n");
#endif
done:
  return err;
}

bool
littleEndian() {
  ByteEndian byteEndian = getWordOrganization();
  switch(byteEndian) {
  case BYTE_BIG_ENDIAN:
    return false;
  case BYTE_LITTLE_ENDIAN:
    return true;
  default:
    // We weren't able to figure the Endian out.  Leave the swap defaults as they are.
#ifdef INTV_WORDS_BIGENDIAN
    return true;
#else
    return false;
#endif
  }
}

int
main(int argc, char *argv[]) {
  if (argc < 5) {
    fprintf(stderr, "\n%s numFloat numInt pfile infile ...\n\n", argv[0]);
    exit(1);
  }
  int   numFloat  = atoi(argv[1]);
  int   numInt    = atoi(argv[2]);
  char *pfileName = argv[3];
  int   firstFile = 4;
  bool  doWeSwap  = littleEndian();

#if DEBUG
  fprintf(stderr, "writing to %s\n", pfileName);
#endif

  FILE *pfilef = fopen(pfileName, "w");
  if (!pfilef) {
    perror(pfileName);
    exit(1);
  }
  OutFtrLabStream_PFile *pfile = 
    new OutFtrLabStream_PFile(0, "", pfilef, (size_t) numFloat, (size_t) numInt, 1, doWeSwap);
  float floatBuf[numFloat];
  int32_t intBuf[numInt];

  int err = 0;
  for (int i=firstFile; i < argc; i+=1) {
#if DEBUG
    fprintf(stderr, "reading %s\n", argv[i]);
#endif
    FILE *in = fopen(argv[i], "r");
    if (!in) {
      err = errno;
      perror(argv[i]);
      goto done;
    }
    bool eof = false;
    int curSent;

    err = readFrame(pfile, argv[i], in, eof, curSent, numFloat, numInt, floatBuf, intBuf);
    if (err) goto done;
    int prevSent = curSent;

    while (!eof) {
      err = readFrame(pfile, argv[i], in, eof, curSent, numFloat, numInt, floatBuf, intBuf);
      if (err) goto done;
      if (curSent != prevSent) {
#if DEBUG
	fprintf(stderr, "finished sentence %d\n", prevSent);
#endif
	pfile->doneseg(SEGID_UNKNOWN);
	prevSent = curSent;
      }
    }
    pfile->doneseg(SEGID_UNKNOWN);

    if (fclose(in)) {
      err = errno;
      perror(argv[i]);
    }
  }
  delete pfile;

done:
  exit(err);
}
