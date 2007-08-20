#ifndef PFILE_H_INCLUDED
#define PFILE_H_INCLUDED

// This file contains some miscellaneous routines for handling PFiles

// IMPORTANT NOTE - some places in this file refer to "sentence" - this
// is synonymous with "segment", the later being a generalization of the
// former with no functional difference.

#include <stdio.h>
#include "error.h"

// This is the size of the PFile header
#define PFILE_HEADER_SIZE (32768)

typedef size_t UInt32;
typedef int Int32;

typedef int intv_int32_t;
typedef unsigned int intv_uint32_t;

typedef long SegID;
enum {
    SEGID_UNKNOWN = 0,	// Returned when we do not know the segid.
    SEGID_BAD = -1		// Used for passing back exceptions.
}; 

enum
{
    SIZET_BAD = 0xffffffffu, // FIXME - this should be configured.
    SIZET_MAX = 0xfffffffeu,	// The maximum size of anything - infinity!
    ALL = SIZET_BAD	// If you want as many as possible of somat.
};

enum
{
    OK = 0,
    BAD = -1,
    ENDOF = -1
};

// The standard PFile 32 bit data type

typedef union
{
    Int32 l; float f;
} PFile_Val;

// Adapted from the latest (2004-10-26) Quicknet distribution to
// support 64 bit pfiles. I am defining _FILE_OFFSET_BITS here but it
// should be set during config time.  -- Karim

#define _FILE_OFFSET_BITS 64

// Some stuff to deal with large file support
//  First make sure things still work if we do not have fseeko/ftello
#ifdef PF_HAVE_FSEEKO
#define pfile_fseek(a,b,c) fseeko(a,b,c)
#define pfile_ftell(a) ftello(a)
typedef off_t pfile_off_t;
#else
#define pfile_fseek(a,b,c) fseek(a,b,c)
#define pfile_ftell(a) ftell(a)
typedef long pfile_off_t;
#endif

// Set up long long types if we have them, along with approriate format
// string segments
#if defined(PF_HAVE_LONG_LONG)
typedef long long pfile_longlong_t;
typedef unsigned long long pfile_ulonglong_t;
#define PF_LLU "%llu"
#define PF_LLD "%lld"
#else
typedef long pfile_longlong_t;
typedef unsigned long pfile_ulonglong_t;
#define PF_LLU "%lu"
#define PF_LLD "%ld"
#endif

// Define PFILE_LARGE if we support large PFiles
#if (defined(PF_HAVE_FSEEKO) && defined(PF_HAVE_LONG_LONG) && (_FILE_OFFSET_BITS==64))
#define PFILE_LARGEFILES 1
#else
#define PFILE_LARGEFILES 0
#endif


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// InFtrLabStream_PFile - access a PFile for input
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// This is the lowest level access to a PFile, returning feature and
// label data simultaneously.  For general use, the "FtrStream" and
// "LabStream" based interfaces are preferable.

class InFtrLabStream_PFile {

 public:
    // Constructor
    // "a_debug" controls the level of status message output.
    // "a_dbgname" is appended to status messages.
    // "a_file" is the stream we read the features and labels from.
    // "a_indexed" is non-zero if we want random access to the PFile.

    InFtrLabStream_PFile(int a_debug, 
			   const char* a_dbgname,
			   FILE* a_file, 
			   int a_indexed,
                           short swap);

    ~InFtrLabStream_PFile();

    // Return the number of labels in the PFile
    size_t num_labs();

    // Return the number of features in the PFile
    size_t num_ftrs();

    // Return the number of segments in the PFile
    size_t num_segs();

    // Read the next set of frames, up until the end of segment.
    // Returns the number of frames read, 0 if already at end of sentence
    size_t read_ftrslabs(size_t frames, float* ftrs, UInt32* labs);

    // Read just features from the next set of frames.
    // Returns the number of frames read, 0 if already at end of segment.
    // NOTE - once read_ftrs has been called, the label values for that
    // frame have been lost.  Use read_ftrslabs() if you want both features
    // and labels
    size_t read_ftrs(size_t frames, float* ftrs);

    // Read just labels from the next set of frames.
    // Returns the number of frames read, 0 if already at end of segment.
    // NOTE - once read_labs has been called, the feature values for that
    // frame have been lost.  Use read_ftrslabs() if you want both features
    // and labels
    size_t read_labs(size_t frames, UInt32* labs);

    // Return details of where we are (next frame to be read)
    // Always returns OK.
    int get_pos(size_t* segnop, size_t* framenop);

    // Move on to the next segment.
    // Returns ID of next segment if succeeds, else SEGID_BAD if at end
    // of file.
    SegID nextseg();

    // Move back to the beginning of the PFile.  This will work for
    // unindexed PFiles, but not streams.  If this fails, it returns BAD.
    // After a rewind, nextseg() must be used to move to the start of the
    // first segment.
    int rewind();

// The following operations are only available if we selected indexing in
// the constructor

    // Return the number of frames in the given segment.
    // (or the whole file if segno==ALL)
    // Returns SIZET_BAD if info not known

    size_t num_frames(unsigned int); 

    size_t num_frames();

    SegID set_pos(size_t segno, size_t frameno);

private:
    // Some constants for pfiles
    enum
    {
	SENT_EOF = -1		// A suitable segment no. to mark EOF
    };

    FILE* const file;		// The stream for reading the PFile.
    const char *filename;       // the filename
    const int indexed;		// 1 if indexed.
    unsigned int num_cols;	// Columns in data section.
    unsigned int total_frames;	// Frames in data section.
    unsigned int total_sents;	// The number of segments in the data section.
    unsigned int first_ftr_col;	// First column of features.
    unsigned int num_ftr_cols;	// Number of feature columns.
    unsigned int first_lab_col;	// First column of labels.
    unsigned int num_lab_cols;	// Number of label columns.
// Do not support targets yet
//    unsigned int first_target_col; // First column of targets.
//    unsigned int num_target_cols; // Number of target columns.

    pfile_longlong_t data_offset;		// Offset of data section from pfile header.
    pfile_longlong_t sentind_offset;	// The offset of the segment index section
				// ..in the pfile.
    size_t bytes_in_row;		// Bytes in one row of the PFile.

    long current_sent;		// Current segment number.
    long current_frame;		// Current frame number within segment.
    long current_row;		// Current row within pfile.

    long pfile_sent;		// Segment number read from PFile.
    long pfile_frame;		// Frame number read from PFile.

    PFile_Val* buffer;	// A buffer for one frame.
    UInt32* sentind;		// The segment start index.

    short bswap;       // byteswapping flag




//// Private functions

    // Read the PFile header.
    void read_header();
    // Read one frame of PFile data.
    void read_frame();
    // Two diffent implementations for "build_index", depending on whether
    // there is already a segment index section in the pfile.
    void build_index_from_sentind_sect();
    void build_index_from_data_sect();
};

// Return the number of label columns in the PFile.

inline size_t
InFtrLabStream_PFile::num_labs()
{
    return num_lab_cols;
}

// Return the number of feature columns in the PFile.
inline size_t
InFtrLabStream_PFile::num_ftrs()
{
    return num_ftr_cols;
}

// Return the number of segments in the PFile.
inline size_t
InFtrLabStream_PFile::num_segs()
{
    return total_sents;
}

// Return just the feature values.
inline size_t
InFtrLabStream_PFile::read_ftrs(size_t frames, float* ftrs)
{
    return read_ftrslabs(frames, ftrs, NULL);
}

// Return just the label values.
inline size_t
InFtrLabStream_PFile::read_labs(size_t frames, UInt32* labs)
{
    return read_ftrslabs(frames, NULL, labs);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// OutFtrLabStream_PFile - access a PFile for output
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// This is the lowest level access to a PFile, requiring feature and
// label data simultaneously.  For general use, the "OutFtrStream" and
// "OutLabStream" based interfaces are preferable.

class OutFtrLabStream_PFile 
{
public:
    // Constructor.
    // "a_debug" controls the level of status message output.
    // "a_dbgname" is the debugging output tag.
    // "a_file" is the stream we write the output data to.
    // "a_ftrs" is the number of features in the resulting PFile.
    // "a_labels" is the number of labels in the resulting PFile.
    // "a_indexed" is non-zero if we want an indexed PFile.

    OutFtrLabStream_PFile(int a_debug, const char* a_dbgname,
			     FILE* a_file, size_t a_ftrs, size_t a_labs,
			     int a_indexed,short swap);
    ~OutFtrLabStream_PFile();

    // Return the number of labels in the PFile.
    size_t num_labs();

    // Return the number of features in the PFile.
    size_t num_ftrs();

    // Write the next set of frames, keeping them part of the current segment.
    void write_ftrslabs(size_t frames, const float* ftrs,
		       const UInt32* labs);

    // Wrt. just features - compatible with OutFtrStream abstract interface.
    // NOTE - cannot use write_labs then write_ftrslabs - must use
    // write_ftrslabs if we need none-zero values for both.
    void write_ftrs(size_t frames, const float* ftrs);

    // Write just labels - compatible with OutLabStream abstract interface.
    // NOTE - cannot use write_ftrs then write_labs - must use
    // write_ftrslabs if we need non-zero values for both.
    void write_labs(size_t frames, const UInt32* labs);


    // Finish writing the current segment, identify the current segment
    // then move on to the next segment.
    void doneseg(SegID segid);

private:
    FILE* const file;		// The stream for reading the PFile.
    const char *filename;             // the file name
    const int indexed;		// 1 if indexed.

    const unsigned int num_ftr_cols; // Number of feature columns.
    const unsigned int num_lab_cols; // Number of label columns.

    long current_sent;		// The number of the current segment.
    long current_frame;		// The number of the current frame.
    long current_row;		// The number of the current row.
  


    PFile_Val* buffer;	// A buffer for one frame of PFile data.
    char* header;		// Space for the header data.
    // The default length of the sentence index - it grows by doubling in
    // size.
    // Note that this is small so index enlargement can be tested with small
    // test files.
    enum { DEFAULT_INDEX_SIZE = 32 }; 
    UInt32* index;		// Sentence index (if building indexed file).
    size_t index_len;		// The length of the sentence index that has
				// been allocated (in words).

    short bswap;               // byte swapping on output - yes/no

// Private member functions.
    // Write the header information at the start of the file.
    void write_header();

    // Write the sentence index information at the end of the file.
    void write_index();
};

// Return the number of labels in the PFile
inline size_t
OutFtrLabStream_PFile::num_labs()
{
    return num_lab_cols;
}

// Return the number of features in the PFile
inline size_t
OutFtrLabStream_PFile::num_ftrs()
{
    return num_ftr_cols;
}

// Write just features to the file
inline void
OutFtrLabStream_PFile::write_ftrs(size_t frames, const float* ftrs)
{
    write_ftrslabs(frames, ftrs, NULL);
}

// Write just labels to the file

inline void
OutFtrLabStream_PFile::write_labs(size_t frames, const UInt32* labs)
{
    write_ftrslabs(frames, NULL, labs);
}

#endif // #ifndef PFILE_H_INCLUDED

