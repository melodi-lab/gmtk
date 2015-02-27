
/*
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2014 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#if HAVE_HG_H
#  include "hgstamp.h"
#endif

#include "general.h"
#include "error.h"
#include "debug.h"
#include "arguments.h"
#include "vbyteswapping.h"
#include "rand.h"

#include "GMTK_WordOrganization.h"
#include "GMTK_ObservationArguments.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_ObservationFile.h"
#include "GMTK_FileSource.h"
#include "GMTK_CreateFileSource.h"
#include "GMTK_ASCIIFile.h"
#include "GMTK_FlatASCIIFile.h"
#include "GMTK_PFileFile.h"
#include "GMTK_HTKFile.h"
#include "GMTK_HDF5File.h"
#include "GMTK_BinaryFile.h"
#include "GMTK_Filter.h"
#include "GMTK_FilterFile.h"
#include "GMTK_FIRFilter.h"
#include "GMTK_AffineFilter.h"
#include "GMTK_Stream.h"

RAND rnd(false);
FileSource *gomFS;
ObservationSource *globalObservationMatrix;

#define GMTK_ARG_OBS_FILES
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_DCDRNG
#define GMTK_ARG_START_END_SKIP
#define GMTK_ARG_OBS_MATRIX_OPTIONS
#define GMTK_ARG_OBS_MATRIX_XFORMATION
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION
#define GMTK_ARGUMENTS_DEFINITION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DEFINITION

unsigned obsOffset    = 0;
unsigned numFeatures  = 0;
unsigned radius       = 0;
unsigned verb = IM::Default;

char *outputFileName = NULL;
char *outputListName = NULL;
char *outputNameSeparator = NULL;
char const *outputFileFmt = "pfile";
bool outputFileSwap = false;

Arg Arg::Args[] = {
  Arg("\n*** Window Options ***\n"),
  Arg("featureOffset", Arg::Opt, obsOffset, "Offset in observation file where input features start"),
  Arg("numFeatures", Arg::Req, numFeatures, "Number of input features (per frame)"),
  Arg("radius", Arg::Opt, radius, "output frames consist of windows of 1+2r input frames"),

  Arg("outputFileName", Arg::Opt, outputFileName, "Output filename for windowed observation file"),
  Arg("outputListFileName", Arg::Opt, outputListName, "Output list filename for windowed observation file (HTK, ASCII, Binary)"),
  Arg("outputNameSeparator", Arg::Opt, outputNameSeparator, "String to use as separator when outputting HTK, ASCII, or binary windowed observation file"),
  Arg("outputFileFormat", Arg::Opt, outputFileFmt, "Output windowed observation file format (htk,binary,ascii,flatascii,pfile)"),
  Arg("outputFileSwap", Arg::Opt, outputFileSwap, "Do byte swapping when outputting PFile, HTK, or binary windowed observation file"),


#define GMTK_ARGUMENTS_DOCUMENTATION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION
  Arg("verbosity",Arg::Opt, verb, "Level of debugging output, [0-100]"),
  // final one to signal the end of the list
  Arg()
};

int 
main(int argc, char *argv[]) {

  CODE_TO_COMPUTE_ENDIAN

  bool parse_was_ok = Arg::parse(argc,(char**)argv);

  if(help) {
    Arg::usage();
    exit(0);
  }
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }
  if (print_version_and_exit) {
#ifdef HAVE_CONFIG_H
    printf("%s (Mercurial id: %s)\n",gmtk_version_id,HGID);
#else
    printf("%s\n", gmtk_version_id);
#endif
    exit(0);
  }
  infoMsg(IM::Max,"Finished parsing arguments\n");
#define GMTK_ARGUMENTS_CHECK_ARGS
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS

  if (startSkip < (int)radius) {
    error("-startSkip must be >= -radius\n");
  }
  if (endSkip < (int)radius) {
    error("-endSkip must be >= -radius\n");
  }

  infoMsg(IM::Max,"Opening Files ...\n");
  gomFS = instantiateFileSource();
  globalObservationMatrix = gomFS;

  if (obsOffset >= gomFS->numContinuous()) {
    error("ERROR: -featureOffset must be < number of continuous features (%u)\n", gomFS->numContinuous());
  }
  if (obsOffset + numFeatures > gomFS->numContinuous()) {
    error("ERROR: -featureOffset + -numFeatures must be <= number of continuous features (%u)\n", gomFS->numContinuous());
  }

  infoMsg(IM::Max,"Finished opening files.\n");

  gomFS->setMinPastFrames( radius );
  gomFS->setMinFutureFrames( radius );

  if (gomFS->numSegments()==0)
    error("ERROR: no segments are available in observation file");
  Range* dcdrng = new Range(dcdrng_str,0,gomFS->numSegments());
  if (dcdrng->length() <= 0) {
    infoMsg(IM::Default,"Decoding range '%s' specifies empty set. Exiting...\n",
          dcdrng_str);
    exit_program_with_status(0);
  }

  unsigned diameter = 1 + 2 * radius;
  unsigned numWindowFeatures = diameter * numFeatures;
  float *output_vector = new float[numWindowFeatures];
  unsigned stride = gomFS->stride();

  ObservationFile *outputFile = instantiateWriteFile(outputListName, outputFileName, const_cast<char *>(outputNameSeparator), 
						     const_cast<char *>(outputFileFmt), numWindowFeatures, 0, outputFileSwap);

  Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());
  while (!dcdrng_it->at_end()) {
    const unsigned segment = (unsigned)(*(*dcdrng_it));
    if (gomFS->numSegments() < (segment+1)) 
      error("ERROR: only %d segments in file, segment must be in range [%d,%d]\n",
            gomFS->numSegments(),
            0,gomFS->numSegments()-1);

    infoMsg(IM::Max,"Processing segment %d ...\n",segment);
    if (!gomFS->openSegment(segment)) {
      error("ERROR: failed to open segment %u\n", segment);
    }
    const unsigned numFrames = gomFS->numFrames();

    //  Data32 const *frame;

    for (unsigned frame=0; frame < numFrames; frame+=1) {

      float *dest = output_vector;
      // guarantees [frame - radius, frame + radius] are in cache
      float *src  = gomFS->floatVecAtFrame(frame) - radius * stride + obsOffset; 
      for (unsigned i = 0; i < diameter; i+=1, src += stride, dest += numFeatures) {
	memcpy(dest, src, numFeatures * sizeof(float));
      }
      outputFile->writeFrame((Data32 *)output_vector);
    }
    outputFile->endOfSegment();
    (*dcdrng_it)++;
  }
  delete outputFile;

  exit_program_with_status(0);
}

