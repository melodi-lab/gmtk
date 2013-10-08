/*
 * version.h
 * Keeps track of the major and minor GMTK version number.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 *
 * $Header$
 *
 */

#if HAVE_CONFIG_H
#include <config.h>
static const char * gmtk_version_id = PACKAGE_STRING;
#else
// TODO: automate the process of updating this string.
static const char * gmtk_version_id = "GMTK Version 0.2b Tue Jan 20 22:59:41 2004";
#endif

// TODO: figure out a better way to keep compiler from complaining
static const char *___gvi_tmp___ = gmtk_version_id;


