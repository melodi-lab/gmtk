/*
 * GMTK_CreateFileSource.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_CREATEFILESOURCE_H
#define GMTK_CREATEFILESOURCE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "GMTK_FileSource.h"

FileSource *instantiateFileSource();

void instantiateFileSource(FileSource *source);

#endif
