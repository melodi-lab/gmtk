/*-
 * machine_dependent.h
 *        Definition of a variety of machine dependnet functions, types, etc.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#ifndef MACHINE_DEPENDENT_H
#define MACHINE_DEPENDENT_H


// 16 bit int
typedef short Int16;
// 16 bit unsigned
typedef unsigned short Uint16;
// generic 16 bit data value
typedef short Data16;


// 32 bit int
typedef int Int32;
// 32 bit unsigned
typedef unsigned int Uint32;
// generic 32 bit data value
typedef int Data32;
// 32 bit pointer
typedef Data32* Datap32;
// generic 64 bit data value
typedef double Data64;
// 64 bit pointer
typedef Data64* Datap64;


#endif






