/*
 * GMTK_FileSource.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILESOURCE_H
#define GMTK_FILESOURCE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"
#include "GMTK_FilterFile.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_Filter.h"

// FileSource handles random access data sources like the various 
// binary file formats. The FileSource class itself implements 
// frame caching and prefetching. It acts as an adaptor between
// the ObservationSource interface (see GMTK_ObservationSource.h)
// used by the inference code and the ObservationFile interface 
// (see GMTK_ObsevationFile.h) implemented by data file formats
// that support random access. 

// There are a couple of helper classes - FilterFile and MergeFile - that
// also implement the ObservationFile interface. (See GMTK_FilterFile.h
// and GMTK_MergeFile.h) The responsibilities for implementing the
// GMTK program command line arguments related to the observation 
// data are split up amoung instances of the classes as shown below

/*
 *
 *                 FileSource             -justification -startSkip -endSkip -constantSpace
 *                     |                  -fileBufferSize -fileWindowSize -fileWindowDelta
 *                     |                                        
 *                     |                                        
 *                 FilterFile             -posttrans -gpr       
 *                     |                                        
 *                     |                                        
 *                 MergeFile              -comb -fdiffactX -sdiffactX
 *                  /     \                                     
 *                 /       \                                    
 *        FilterFile       FilterFile     -transX -frX -irX -postprX
 *            |       ...      |                                
 *            |                |                                
 *         HDF5File        BinaryFile     -ofX -fmtX -niX -nfX -iswpX -srX 
 *                                        -prefrX -preirX -preprX
 *
 *
 */

// At the lowest level are instances of subclasses of ObservationFile
// that implement the file format specified by -fmtX. For convenience and
// consistancy, the instantiateFileSource() function (see GMTK_CreateFileSource.h)
// instantiates the above object hierarchy according to the command
// line arguments given to a GMTK program.

class FileSource: public ObservationSource {

 protected:

  // FileSource prefetches and caches transformed frames into the "cooked buffer"
  // (as opposed to untransformed "raw frames" - though this distinction became
  // less relevant when -posttrans was refactored out of FileSource). The FileSource
  // doesn't know whether the inference code is doing a forward or backward pass,
  // so the first group of prefected data is actually stored in the middle of the
  // cookedBuffer. When the inference code gets within delta frames of either end
  // of the currently cached frame range, it will append another window's worth of 
  // frames to that end. This means that half the cookedBuffer will be wasted,
  // but that is deemed acceptable compared to the memory saved by not loading
  // the entire segment of observed data into memory at once.

  // Note that the above only applies in constant space mode - in the default
  // O(T) space mode, the cookedBuffer is just resized to hold at least the number
  // of frames in the new segment on each openSegment() call, which are all loaded
  // at once.

  Data32 *cookedBuffer; // buffer to hold transformed cahced/prefected observations
  unsigned bufStride;   // increment between frames in cookedBuffer

  unsigned firstBufferedFrame;        // frame # of first frame currently in the cookedBuffer
  unsigned firstBufferedFrameIndex;   // index (in frames) in cookedBuffer where the first buffered frame starts
  unsigned numBufferedFrames;         // # of frames currently in cookedBuffer
  unsigned bufferFrames;              // size of cookedBuffer in frames
  unsigned bufferSize;                // size of cookedBuffer in Data32

  unsigned window;                    // # of frames to load at once
  unsigned delta;                     // load more frames when within delta frames of the edge of the window

  ObservationFile *file;              // Where the observation data comes from.
 
  // Due to the boundary algorithm M and S parameters or other size mis-matches,
  // a model may not be able to use all of the frames in a segment. The frames
  // used by the model can be left-most, center, or right-most of the segment.

  int justificationMode;              // left, center, right
  unsigned justificationOffset;       // # of frames to skip to effect the justification

  unsigned _startSkip;                // -startSkip frames
  unsigned _endSkip;                  // -endSkip frames

  // dlinks, VE CPTs, and filters may require a number of frames from before or 
  // after the frame being computed. Thus, if the inference code requests frames
  // [i,j] the cookedBuffer must actually contain [i-_minPastFrames, j+_minFutureFrames]
  // in order to satisfy the request.

  unsigned _minPastFrames;
  unsigned _minFutureFrames;

  unsigned _numCacheableFrames;   // after considering -fdiffactX  & -gpr only
                                  // This is the # of frames that are eligable for caching
                                  // in the cookedBuffer

  unsigned _numFrames;            // after -fdiffactX, -gpr, -startSkip, -endSkip, -justification
                                  // This is the # of frames actually accessible to clients

  // _numFrames <= _numCacheableFrames.  _numFrames will be less when -startSkip
  // or -endSkip reserve some frames at the start or end of each segment. These 
  // reserved frames (for dlinks, VECPTs) must be present in the cache when the first
  // or last several frames are accessed, but are not directly accessable to
  // clients via loadFrames(), floatVecAtFrame(), etc.

  int      segment;               // currently open segment; -1 if none yet

  bool     constantSpace;         // if true,  load only O(1) frames at a time
                                  // if false, load the entire segment, resizing the 
                                  //           cookedBuffer, if necessary

  // number of continuous/discrete/total features in the combined observation data
  unsigned numContinuousFeatures; 
  unsigned numDiscreteFeatures;   
  unsigned _numFeatures;

  // Load requested frames into cookedBuffer, starting at the specified index.
  // index is in frames, so the frames will start at 
  // cookedBuffer[bufferIndex * buffStride]
  Data32 const *loadFrames(unsigned bufferIndex, unsigned first, unsigned count);

 public:

#define DEFAULT_BUFFER_SIZE       (16 * 1024 * 1024)
#ifndef DEFAULT_FILE_WINDOW_SIZE
#define DEFAULT_FILE_WINDOW_SIZE  (4)
#define DEFAULT_FILE_WINDOW_DELTA (100)
#endif
#define DEFAULT_FILE_WINDOW_BYTES (DEFAULT_FILE_WINDOW_SIZE * 1024 * 1024)

  FileSource(ObservationFile *file,
	     unsigned windowBytes = DEFAULT_FILE_WINDOW_BYTES, 
	     unsigned deltaFrames = DEFAULT_FILE_WINDOW_DELTA,
	     unsigned bufferSize  = DEFAULT_BUFFER_SIZE, 
	     unsigned startSkip=0, unsigned endSkip=0,
	     int justificationMode=0, bool constantSpace = false); 

  // Construct an invalid empty FileSource. Turn it into a valid one
  // with the initialize() method.
  FileSource() {
    cookedBuffer = NULL;
    bufferSize = 0;
    bufStride = 0;
    window = 0;
    delta = 0;
    numBufferedFrames = 0;
    file = NULL;
    _startSkip = 0;
    _endSkip = 0;
    justificationMode = 0;
    justificationOffset = 0;
    _minPastFrames = 0;
    _minFutureFrames = 0;
    segment = -1;
    constantSpace = false;
  }

  virtual ~FileSource() {
    if (cookedBuffer) delete [] cookedBuffer;
    if (file) delete file;
  }

  // Turn an invalid FileSource created by the no-arg ctor into a
  // valid FileSource.
  void initialize(ObservationFile *file,
		  unsigned windowBytes = DEFAULT_FILE_WINDOW_BYTES, 
		  unsigned deltaFrames = DEFAULT_FILE_WINDOW_DELTA,
		  unsigned bufferSize  = DEFAULT_BUFFER_SIZE,
		  unsigned startSkip=0, unsigned endSkip=0,
		  int justificationMode = 0, bool constantSpace = false);

  // The number of available segments.
  unsigned numSegments() { 
    assert(file);
    return file->numLogicalSegments(); 
  }

  // The number of ObservationFiles combined into the observation matrix
  unsigned numFiles() {
    assert(file);
    return file->numFiles();
  }

  // The current segment 
  unsigned segmentNumber() {
    assert(segment >= 0);
    return (unsigned) segment;
  }

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  bool openSegment(unsigned seg);

  // Apply left, center, or right justification to the usable frames
  // in the current segment. numUsableFrames must be <= _numFrames
  // prior to calling justifySegment. After the call, _numFrames = numUsableFrames

  // The number of usable frames in a segment depends on properties
  // of the model (dlinks, triangulation parameters, etc). JunctionTree::unroll()
  // returns the number of frames that can actually be used. So, the outer 
  // loops of the GMTK programs that process observation data go something like this:
  //
  // while (!segment_iterator->at_end()) {
  //   segment =  (unsigned)(*(*segment_iterator));
  //   fileSouce->openSegment(segment);
  //   numFrames = fileSource->numFrames();         // Returns # of frames in the segment after
  //                                                // -gpr -startSkip and -endSkip are applied
  //
  //   numUsableFrames = myjt.unroll(numFrames);    // compute # frames usable by the model
  //   fileSource->justifySegment(numUsableFrames); // perform requested justification
  // }
  void justifySegment(unsigned numUsableFrames);

  bool active(); // Returns true iff a segment has been activated by an
                 // openSegment() call.

  // The number of frames available in the currently open segment.
  unsigned numFrames() {
    assert(segment >= 0);
    return _numFrames;
  }

  // Returns a pointer to the requested frames in the cookedBuffer.
  // Any required frames before/after (see _minPastFrames and
  // _minFutureFrames above) are guaranteed to be in the cookedBuffer
  // as well. 
  virtual Data32 const *loadFrames(unsigned first, unsigned count);


  // The number of continuous, discrete, total features

  unsigned numContinuous() {
    assert(file);
    return numContinuousFeatures;
  }


  unsigned numDiscrete() {
    assert(file);
    return numDiscreteFeatures;
  }


  unsigned numFeatures() {
    assert(file);
    return _numFeatures;
  }


  // The number of Data32's between each frame in the cookedBuffer
  unsigned stride() {
    assert(file);
    return _numFeatures;   // with the new implementation, it's always the case that stride == numFeatures
  }


  // number of frames to skip at the beginning of the segment
  unsigned startSkip() {return _startSkip;};

  // number of frames to skip at the end of the segment
  unsigned endSkip() {return _endSkip;} 

  unsigned minPastFrames() {return _minPastFrames;}
  unsigned minFutureFrames() {return _minFutureFrames;}

  void setMinPastFrames(unsigned n) {
    if (_startSkip < n) {
      error("ERROR: the model requires a -startSkip of at least %u (currently %u)",
	    n, _startSkip);
    }
    _minPastFrames = n;
  }

  void setMinFutureFrames(unsigned n) {
    if (_endSkip < n) {
      error("ERROR: the model requires an -endSkip of at least %u (currently %u)",
	    n, _endSkip);
    }
    _minFutureFrames = n;
  }

  // Note that the following call loadFrames() internally

  float *const floatVecAtFrame(unsigned f);

  float *const floatVecAtFrame(unsigned f, const unsigned startFeature);

  unsigned *const unsignedVecAtFrame(unsigned f);

  unsigned &unsignedAtFrame(const unsigned frame, const unsigned feature);

  Data32 const * const baseAtFrame(unsigned f);

};

#endif
