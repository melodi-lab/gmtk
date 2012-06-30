
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#ifdef HAVE_CONFIG_H
#  include <config.h>
static const char * gmtk_version_id = PACKAGE_STRING;
#  ifdef HAVE_HG_H
#    include "hgstamp.h"
#  endif

#else 
// TODO: automate the process of updating this string.
static const char * gmtk_version_id = "GMTK Version 0.2b Tue Jan 20 22:59:41 2004";
#endif

#include "arguments.h"
#include "pfile.h"
#include "GMTK_WordOrganization.h"

#define DEBUG 0

int
readFrame(OutFtrLabStream_PFile *pfile, char *fname, FILE *in, bool &eof, 
	  int &curSent, int numFloat, int numInt, float *floatBuf, int *intBuf) 
{
  int err = 0;
  int prevSent = curSent;

  int nread = fscanf(in, "%d", &curSent);
  eof = nread == EOF;
  if (eof) goto done;
  if (nread != 1) {
    err = errno;
    perror(fname);
    goto done;
  }

  if (prevSent != curSent) {
    pfile->doneseg(SEGID_UNKNOWN);
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


char *pfileName = NULL;
#define MAX_OBJECTS 10

char *input_fname[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};  // Input file name.
unsigned int numInt;
unsigned int numFloat;
bool printVersion = false;
unsigned help=0;  // help=0...5 depending on the amount of info we want printed
Arg Arg::Args[] = {

  Arg("\n*** Input arguments ***\n"),

  Arg("i",    Arg::Req, input_fname,"input file. Replace X with the file number",Arg::ARRAY,MAX_OBJECTS),
  Arg("nf",   Arg::Req, numFloat,"number of floats in input file(s)"),
  Arg("ni",   Arg::Req, numInt,"number of ints (labels) in input file(s)"),

  Arg("\n*** Output arguments ***\n"),

  Arg("o",    Arg::Opt, pfileName,"output file"),

  Arg("\n*** Misc arguments ***\n"),

  Arg("help",  Arg::Help, help,  "Print this message. Add an argument from 1 to 5 for increasing help info."),
  Arg("version", Arg::Tog, printVersion, "Print GMTK version and exit."),
  // The argumentless argument marks the end of the above list.
  Arg()
};

int
main(int argc, char *argv[]) {

  bool parse_was_ok = Arg::parse(argc,(char**)argv);

  if(help) {
    Arg::usage();
    exit(0);
  }
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }
  if (printVersion) {
#ifdef HAVE_CONFIG_H
    printf("%s (Mercurial id: %s)\n",gmtk_version_id,HGID);
#else
    printf("%s\n", gmtk_version_id);
#endif
    exit(0);
  }

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
  float   *floatBuf = new float[numFloat];
  int32_t *intBuf   = new int32_t[numInt];

  int err = 0;
  for (int i=0; i < MAX_OBJECTS && input_fname[i] != NULL; i+=1) {
#if DEBUG
    fprintf(stderr, "reading %s\n", input_fname[i]);
#endif
    FILE *in = fopen(input_fname[i], "r");
    if (!in) {
      err = errno;
      perror(input_fname[i]);
      goto done;
    }
    bool eof = false;
    int curSent = 0;

    err = readFrame(pfile, input_fname[i], in, eof, curSent, numFloat, numInt, floatBuf, intBuf);
    if (err) goto done;

    while (!eof) {
      err = readFrame(pfile, input_fname[i], in, eof, curSent, numFloat, numInt, floatBuf, intBuf);
      if (err) goto done;
#if 0
      if (curSent != prevSent) {
#if DEBUG
	fprintf(stderr, "finished sentence %d\n", prevSent);
#endif
	pfile->doneseg(SEGID_UNKNOWN);
	prevSent = curSent;
      }
#endif
    }
    pfile->doneseg(SEGID_UNKNOWN);

    if (fclose(in)) {
      err = errno;
      perror(input_fname[i]);
    }
  }
  delete pfile;

done:
  exit(err);
}
