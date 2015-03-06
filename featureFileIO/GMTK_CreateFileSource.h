/*
 * GMTK_CreateFileSource.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2011 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 * 
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
