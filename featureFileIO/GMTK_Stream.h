
/* GMTK_FileDescription: stores info related to input stream 
 * 
 * Written by Katrin Kirchhoff <katrin@ee.washington.edu>
 *
 * Modified by Karim Filali <karim@cs.washington.edu> to add the
 * option to pipe the list of file names through CPP.  Made a few
 * other minor "bug" fixes.
 *
 * $Header$
 *
 * Copyright (c) 2001
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 */

#ifndef GMTK_Stream_h
#define GMTK_Stream_h

#include <ctype.h>
#include "bp_range.h"
#include "range.h"
#include "pfile.h"

// possible file formats

enum {
  RAWBIN, 
  RAWASC, 
  PFILE,
  HTK,
  FLATBIN,
  FLATASC
};


#define MAXSTRLEN 1024 // max length of input file name
#define MAXFRAMES 400  // initial data buffer size (in frames)

#define BASEMASK  077 // Mask to remove HTK par. kind qualifiers 

//HTK parameter kinds - needed to determine HTK file type and for sanity checking

enum {
  WAVEFORM,
  LPC,LPREFC,LPCEPSTRA, LPDELCEP,
  IREFC,
  MFCC,
  FBANK,
  MELSPEC,
  USER,
  DISCRETE,
  ANON
};


/* file description class: contains information about (list of) input files */

class StreamInfo {
  bool     cppIfAscii;
  char*    cppCommandOptions;
  size_t   numFileNames;               // number of filenames
  unsigned nFloatsUsed;                // number of floats actually used (of input)
  unsigned nIntsUsed;                  // number of ints actually used (of input)
  
  size_t   prrngCurNumFrames;          // size of current data file
  size_t   afterTransformCurNumFrames;

  size_t calcNumFileNames(FILE*&f);
public:

  unsigned nFloats;                    // number of floats (cont. features) in file
  unsigned nInts;                      // number of ints (disc. features) in file
  unsigned dataFormat;                 // file format

  

  char *   fofName;                    // this file's file name (name of list of file names)
  FILE *   fofFile;                    // this file (list of file names)
  size_t   fullFofSize;                // full size of list of file names, i.e., before the sentence range is applied
  Range *  srRng;                      // sentence range
  size_t   fofSize;                    // size of list of file names, after the sentence range is applied



  bool bswap;              // true if file needs to be byte-swapped

  InFtrLabStream_PFile *pfile_istr;  // pfile input stream

  FILE *curDataFile;

  char **dataNames;        // pointers to individual filenames (into fofBuf)
  size_t curNumFrames;      // size of current data file
 
  size_t curPos;           // index of current position in list of files

  BP_Range *cont_rng;      // range of cont. features used
  BP_Range *disc_rng;      // range of disc. features used

  StreamInfo(const char *,
	     const char *, 
	     const char *,
	     unsigned *, 
	     unsigned *, 
	     unsigned *,
	     bool,
	     unsigned,
	     bool cppIfAscii=false,
	     char* cppCommandOptions=NULL,
	     const char* sr_range_str=NULL);

  ~StreamInfo();

  size_t   readFof(FILE *);       // read file of file names

  size_t   getFofSize()            { return fofSize; } 
  void     setFofSize(size_t size) { fofSize = size; } 

  size_t   getNumFileNames()  { return numFileNames; }
  unsigned getNumFloatsUsed() { return nFloatsUsed;  }
  unsigned getNumIntsUsed()   { return nIntsUsed;    }
  unsigned getNumFloats()     { return nFloats;      }
  unsigned getNumInts()       { return nInts;        }

  size_t   getPrrngCurNumFrames()         { return prrngCurNumFrames; }
  void     setPrrngCurNumFrames(size_t n) { prrngCurNumFrames=n;      }
  size_t   getCurNumFrames()              { return curNumFrames;      }
  void     setAfterTransformCurNumFrames(size_t n) { afterTransformCurNumFrames=n;      }
  size_t   getAfterTransformCurNumFrames()         { return afterTransformCurNumFrames; }
  unsigned getDataFormat() { return dataFormat; }
  bool     swap()             { return bswap;   }
  void     setSwap(bool swap) { bswap = swap;   }

  unsigned mapToValueInRange(unsigned segno) { return (unsigned) srRng->index(segno); }
};

#endif
