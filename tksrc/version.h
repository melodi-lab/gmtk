/*
 * version.h
 * Keeps track of the major and minor GMTK version number.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * $Header$
 *
 */

#if HAVE_CONFIG_H
#include <config.h>
static const char * gmtk_version_id = PACKAGE_VERSION;
#else
// TODO: automate the process of updating this string.
static const char * gmtk_version_id = "GMTK Version 0.2b Tue Jan 20 22:59:41 2004";
#endif

// TODO: figure out a better way to keep compiler from complaining
static const char *___gvi_tmp___ = gmtk_version_id;


