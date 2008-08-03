
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
#include "range.h"
#include "pfile.h"
#include <string>

using namespace std;

// possible file formats

enum {
  RAWBIN, 
  RAWASC, 
  PFILE,
  HTK,
  FLATBIN,
  FLATASC
};


#define MAXSTRLEN (16*1024) // max length of input file name
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


///HTK parameter kind flags
enum {
	HAS_ENERGY=000100, 						// has energy
	ABS_ENERGY_SUPPRESSED=000200, 			//absolute energy suppressed
	HAS_DELTA_COEFS=000400,					//has delta coeﬃcients
	HAS_ACCEL_COEFS=001000,					//has acceleration coeﬃcients
	IS_COMPRESSED=002000, 					//is compressed
	HAS_ZERO_MEAN=004000,					//has zero mean static coef.
	HAS_CRC_CHECKSUM=010000,				//has CRC checksum
	HAS_ZEROTH_CEPSTRAL_COEF=020000			//has 0’th cepstral coef.
};

/** The description of an HTK observations file.
 * Note that this refers to the actual HTK file, and not to a row in the FoF file
 *  (the later can be a subset of the former).
 */
class HTKFileInfo{
public:
	HTKFileInfo(const string& fname, int samp_size, int n_samples, int startOfData,
				bool isCompressed, float* scale, float* offset);
	~HTKFileInfo();

	string fname;  ///the name of the file
	
	int samp_size; /// bytes per sample
	int n_samples; /// number of samples
	int startOfData; ///the offset from start of file to the start of actual data
	
	///HTK kind flags.  See sect. "5.10.1 HTK Format Parameter Files" of the HTK book
	bool isCompressed;
	 
	/**	if the file is compressed, these are the parameters need to reinflate the compressed shorts
	 *	into regular floats.  The memory for them is freed when the instance is destroyed
	 */
	float* scale;
	float* offset;
	  
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

  bool bswap;              // true if file needs to be byte-swapped

public:

  unsigned nFloats;                    // number of floats (cont. features) in file
  unsigned nInts;                      // number of ints (disc. features) in file
  unsigned dataFormat;                 // file format

  

  char *   fofName;                    // this file's file name (name of list of file names)
  FILE *   fofFile;                    // this file (list of file names)
  Range *  srRng;                      // sentence range
  size_t   fullFofSize;                // full size of list of file names, i.e., before the sentence range is applied
  size_t   fofSize;                    // size of list of file names, after the sentence range is applied





  InFtrLabStream_PFile *pfile_istr;  // pfile input stream

  FILE *curDataFile;
  HTKFileInfo* curHTKFileInfo;	//not null if and only if curDataFile is an HTK file 

  char **dataNames;        // pointers to individual filenames (into fofBuf)
  size_t curNumFrames;      // size of current data file
 
  size_t curPos;           // index of current position in list of files

  Range *cont_rng;      // range of cont. features used
  Range *disc_rng;      // range of disc. features used

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

  size_t   getFullFofSize()            { return fullFofSize; } 
  void     setFullFofSize(size_t size) { fullFofSize = size; } 

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
  //  void     setSwap(bool swap) { bswap = swap;   }

  unsigned mapToValueInRange(unsigned segno) { return (unsigned) srRng->index(segno); }
};

#endif
