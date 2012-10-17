#ifndef GMTK_FILE_UTILS
#define GMTK_FILE_UTILS

/*
 * file_utils.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

// Some stuff to deal with large file support
//  First make sure things still work if we do not have fseeko/ftello
#if HAVE_FSEEKO
#  define gmtk_fseek(a,b,c) fseeko(a,b,c)
#  define gmtk_ftell(a) ftello(a)
   typedef off_t gmtk_off_t;
#else
#  define gmtk_fseek(a,b,c) fseek(a,b,c)
#  define gmtk_ftell(a) ftell(a)
   typedef long gmtk_off_t;
#endif


// FIXME - move this to some common constant declaration place

#define MAXSTRLEN (16*1024) // max length of input file name
#define CPP_DIRECTIVE_CHAR '#'

#ifdef DEBUG
#define DBGFPRINTF(_x_) fprintf _x_
#else
#define DBGFPRINTF(_x_)
#endif


/**
 * openCPPableFile -- open an ASCII file that may need to be preprocessed by CPP
 *
 */
FILE *
openCPPableFile(char const *filename, bool cppIfAscii, 
		char const *cppCommandOptions);


/**
 * closeCPPableFile -- close a file opened by openCPPableFile
 *
 */
void
closeCPPableFile(FILE * &f, bool cppIfAscii);


/**
 *  calcNumFileNames -- calculate the number of file names in the 
 *                      file pointed to by the file handle f
 *
 *  pre-conditions: the file f and its file name, fofName, must be initialized 
 *
 *  side effects: if f is a CPP pipe, the file is closed and re-opened
 *  in order to achieve the same effect as a rewind
 *
 */
unsigned 
calcNumFileNames(FILE* &f, char const *fofName, 
		 bool cppIfAscii, char const *cppCommandOptions);


/**
 * readFof -- reads the list of filenames in f into dataNames,
 *            stores the filename count in numFileNames
 * 
 *  pre-conditions: the file f and its file name, fofName, must be initialized 
 *
 *  side effects: if f is a CPP pipe, the file is closed and re-opened
 *  in order to achieve the same effect as a rewind
 *
 */
unsigned
readFof(FILE * &f, char const *fofName, unsigned &numFileNames, char **&dataNames, 
	bool cppIfAscii, char const *cppCommandOptions);

#endif
