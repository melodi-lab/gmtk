/*
 * GMTK_ShowBuildOptions.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2016 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 * 
 *
 */

#ifndef GMTK_SHOWBUILDOPTIONS_H
#define GMTK_SHOWBUILDOPTIONS_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdio.h>
using namespace std;
bool showConfigOptions = false;

void
showBuildOptions() {
  if (showConfigOptions) {
    printf("\nBuild configuration options:\n\n");
#    if PIPE_ASCII_FILES_THROUGH_CPP
    printf("Preprocess ASCII files: YES\n");
    printf("  Current preprocessing command: %s\n", CPP_Command());
#    else
    printf("Preprocess ASCII files: NO\n");
#    endif
#    if ENABLE_GZIP
    printf("Allow GZIP compressed ASCII files: YES\n");
#    else
    printf("Allow GZIP compressed ASCII files: NO\n");
#    endif
#    if ENABLE_BZIP2
    printf("Allow BZIP2 compressed ASCII files: YES\n");
#    else
    printf("Allow BZIP2 compressed ASCII files: NO\n");
#    endif
#    if HAVE_HDF5
    printf("HDF5 support: YES\n");
#    else
    printf("HDF5 support: NO\n");
#    endif
    printf("\n");
  } // if (showConfigOptions)
}

#endif
