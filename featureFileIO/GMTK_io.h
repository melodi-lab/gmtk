//static char *rcsid = "$Id$";

/* GMTK_io: defines basic input and output classes
 * Author: K. Kirchhoff, U Washington <katrin@ee.washington.edu>
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 */

#include <stdio.h>
#include "bp_range.h"
#include "pfile.h"
#include "GMTK_utils.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_FileDescription.h"


#define MAXSTRLEN 1024 // max length of input file name
#define MAXFRAMES 400  // initial data buffer size (in frames)

#define BASEMASK  077         /* Mask to remove HTK par. kind qualifiers */

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




class DataStream {

protected:

  int _numFiles;          // number of input files
  FileDescription **fd;	 // File descriptors
};
  
class DataInStream: public DataStream {

private:

  size_t _curPos;         // index of current position in list
  size_t _numSegs;        // length of input stream

  unsigned _totalFloats;  // number of floats per data record
  unsigned _totalInts;    // number of ints per data record

  ObservationMatrix *inp_buf; // input data buffer 

  // file opening routines

  size_t openBinaryFile(FileDescription *,size_t);
  size_t openAsciiFile(FileDescription *,size_t);
  size_t openHTKFile(FileDescription *,size_t);
  size_t openGMTKFile(FileDescription *,size_t);
  size_t openPFile(FileDescription *,size_t);
  size_t readFof(char *,char **);

  // close current open data files

  void closeDataFiles();
 
public:

   DataInStream(int n_files, 
		const char **fof_names,
		const char **cont_range_str, 
		const char **disc_range_str,
		unsigned *n_floats, 
		unsigned *n_ints, 
		unsigned *formats,
		bool *swapflags);

  ~DataInStream();


  // get names of current open files

  void getFileNames(size_t,char**);

  // get number of ints for all input streams

  void getNumInts(unsigned *); 

  // get number of floats for all input streams

  void getNumFloats(unsigned *);

  // get total number of features in stream

  unsigned getNumFea() { return _totalFloats + _totalInts;}

  // return number of segments (utterances) in stream

  size_t numSegs() { return _numSegs; }

  // get pointer to feature buffer

  ObservationMatrix *getObs();

  // read data for segment i

  size_t readFile(size_t i);

  
};

// output class


class DataOutStream: public DataStream {

  size_t _numFrames;
  size_t _bufSize;

  unsigned _maxFloats;
  unsigned _maxInts;

  FILE **outFiles;
  char **outNames;

  unsigned *dataFormats;

  bool *bswap;

  short *parKind;

  float *ftr_buf;  // temporary buffers for pfile output
  UInt32 *lab_buf;

  OutFtrLabStream_PFile **pfile_ostr; // pfile streams

  unsigned *sampRate;
  unsigned *nFloats;
  unsigned *nInts;

  // writes file header in HTK format

  void writeHTKHeader(FILE *, unsigned,unsigned,unsigned,size_t,short);

  // feature output routines

  void writeBinFeatures(FILE *,ObservationMatrix *,
                        BP_Range *, BP_Range *, BP_Range *);

  void writeAscFeatures(FILE *,ObservationMatrix *, 
                        BP_Range *, BP_Range *, BP_Range *);

  void writePFileFeatures(OutFtrLabStream_PFile *,ObservationMatrix *,
                          BP_Range *,BP_Range *, BP_Range *,size_t);

public:

  DataOutStream(int,
                unsigned *,unsigned *,
                unsigned *,bool *,int *,
                unsigned *,char **);

  ~DataOutStream();

  void writeData(size_t,char **,char **, char *,ObservationMatrix *);

  void initOutStream(unsigned *,unsigned *,unsigned *, bool *,int *,
		     unsigned *,char **,unsigned);

  // set new file names

  void initFiles(char **,char **,size_t);

  void resizeBuffers(size_t);

  // close current open files

  void closeFiles();

};






