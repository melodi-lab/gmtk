// $Header$

#ifndef QN_types_h_INCLUDED
#define QN_types_h_INCLUDED

/* Must include the config.h file first */
#include <QN_config.h>
#include <stdio.h>
#ifdef QN_HAVE_FCNTL_H
#include <fcntl.h>
#endif

typedef unsigned int QNUInt32;
typedef int QNInt32;

// Use this to indicate in bad error info when returning a size_t

enum
{
    QN_SIZET_BAD = 0xffffffffu, // FIXME - this should be configured.
    QN_SIZET_MAX = 0xfffffffeu,	// The maximum size of anything - infinity!
    QN_ALL = QN_SIZET_BAD	// If you want as many as possible of somat.
};

// The type used to uniquely identify a segment (typically, one sentence = one
// segment) within QuickNet.  Note that this does not necessarily correspond
// to the segment/sentence number in any specific file.

typedef long QN_SegID;
enum {
    QN_SEGID_UNKNOWN = 0,	// Returned when we do not know the segid.
    QN_SEGID_BAD = -1		// Used for passing back exceptions.
}; 

#define QN_UINT32_MAX 0xffffffffu;

enum QN_LayerSelector
{
    QN_LAYER0 = 0,
    QN_LAYER1,
    QN_LAYER2,
    QN_LAYER3,
    QN_LAYER4,
    QN_MLP3_INPUT = QN_LAYER0,
    QN_MLP3_HIDDEN = QN_LAYER1,
    QN_MLP3_OUTPUT = QN_LAYER2,
    QN_MLP4_INPUT = QN_LAYER0,
    QN_MLP4_HIDDEN1 = QN_LAYER1,
    QN_MLP4_HIDDEN2 = QN_LAYER2,
    QN_MLP4_OUTPUT = QN_LAYER3,
    QN_LAYER_UNKNOWN
};

// An indicator for different weight sections
enum QN_SectionSelector
{
    QN_LAYER01_WEIGHTS = 0,
    QN_LAYER1_BIAS,
    QN_LAYER12_WEIGHTS,
    QN_LAYER2_BIAS,
    QN_LAYER23_WEIGHTS,
    QN_LAYER3_BIAS,
    QN_MLP3_INPUT2HIDDEN = QN_LAYER01_WEIGHTS,
    QN_MLP3_HIDDENBIAS = QN_LAYER1_BIAS,
    QN_MLP3_HIDDEN2OUTPUT = QN_LAYER12_WEIGHTS,
    QN_MLP3_OUTPUTBIAS = QN_LAYER2_BIAS,
    QN_WEIGHTS_UNKNOWN,
    QN_WEIGHTS_NONE = QN_WEIGHTS_UNKNOWN
};

// SectionSelectors used to be called WeightSelectors
typedef QN_SectionSelector QN_WeightSelector;

// The type of output layer for an MLP
enum QN_OutputLayerType {
    QN_OUTPUT_LINEAR,
    QN_OUTPUT_SIGMOID,
    QN_OUTPUT_SIGMOID_XENTROPY,
    QN_OUTPUT_SOFTMAX
};

// The type of files used for feature/label/activation input/output

enum QN_StreamType {
    QN_STREAM_UNKNOWN,		// Not yet known
    QN_STREAM_PFILE,		// A PFile
    QN_STREAM_BERPFTR,		// An online equivalent of a PFile used in BERP
    QN_STREAM_RAPACT_ASCII,	// Rap-stle output activitions - ASCII
    QN_STREAM_RAPACT_HEX,	// Rap-stle output activitions - hex
    QN_STREAM_RAPACT_BIN,	// Rap-stle output activitions - binary
    QN_STREAM_LNA8,		// Cambridge 8-bit output activations
    QN_STREAM_LNA16,		// Cambridge 16-bit output activations
    QN_STREAM_PREFILE		// Cambridge feature file format
};

// The types of files used for storing weights.

enum QN_WeightFileType {
    QN_WEIGHTFILE_RAP3		// Traditional RAP-style 3 layer weight files.
};


// An indicator for the order of weight matrices of files
enum QN_WeightMaj
{
    QN_INPUTMAJOR,
    QN_OUTPUTMAJOR
};

// Modes for files
enum QN_FileMode
{
    QN_READ = O_RDONLY,
    QN_WRITE = O_RDWR
};

// General error codes

enum
{
    QN_OK = 0,
    QN_BAD = -1,
    QN_EOF = -1
};

extern char* QN_progname;
extern FILE* QN_logfile;

#endif /* #ifndef QN_types_h_INCLUDED */
