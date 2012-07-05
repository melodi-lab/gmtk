#ifndef _GLOBAL_PARAMETERS_
#define _GLOBAL_PARAMETERS_

#include "range.h"

#define ALLOW_REDUNDANT_PAIRS 1
#define CHECK_TUPLE_OVERLAP 1
#define ALLOW_ENTROPY 1

#define ERR_FILE "err-log.mi"

//// Don't change the buffer data type to double
///  The ObservationMatrix library is using floats.
#define BUFFER_DATA_TYPE float
#define PARAM_DATA_TYPE double
#define MAX_POINTER_SET_DIM 20  // == MAX NUM OF VARIABLES
#define MAX_NUM_MIXTURES 20


//////////// Verbosity //////////////////
#define VERBOSE_NUM_SECONDS_PER_PRINT 30
#define VERBOSE_NUM_SECONDS_PER_SENT_PRINT 30
#define VERBOSE_PRINT_FREQUENCY 1
#define VERBOSE_SENT_PRINT_FREQUENCY 10
#define VERBOSE_ACTIVE_PRINT_FREQUENCY 1
#define VERBOSE_MINTIMEPERPRINTNUMACTIVE (-1) ///< always print 


#define NUM_SECONDS_PER_PRINT (240)
#define NUM_SECONDS_PER_SENT_PRINT (120)
#define PRINT_FREQUENCY 10  
#define SENT_PRINT_FREQUENCY 1000 ///< frequency at which sentence processing information is printed
#define ACTIVE_PRINT_FREQUENCY 10
#define MINTIMEPERPRINTNUMACTIVE (60) ///< always print 
/////////////////////////////////////////


#define MINTIMEPERPARMSAVE (60*50)  /// Minimum amount of time between parameter saves, in seconds.

#define MAX_LABEL_VAL 16384

#define MIN(a,b) ((a)<(b)?(a):(b))

#define NO_DATA 0
#define DATA_LEFT 1
#define DONE 2

int readFeatures(Range::iterator krit, size_t &n_frames,
		 size_t &n_samps,
		 Range &lrrng, int labpos,
		 unsigned &frameStart,unsigned &firstFrame);


#define DOUBLE_PROCESSING_DEFINED

#ifdef DEBUG
#define DBGFPRINTF(_x_) fprintf _x_
#else
#define DBGFPRINTF(_x_)
#endif

#endif  // _GLOBAL_PARAMETERS_
