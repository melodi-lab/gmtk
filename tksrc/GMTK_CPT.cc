/*-
 * GMTK_CPT.cc
 *     Trainable (with say EM) CPT
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ieeefp.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"


VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * CPT::CPT()
 *      Constructor
 *
 * Results:
 *      Constructs the object.
 *
 * Side Effects:
 *      None so far.
 *
 *-----------------------------------------------------------------------
 */
CPT::CPT() 
{

}



/*-
 *-----------------------------------------------------------------------
 * CPT::read(is)
 *      read in a distribution from file 'is'. 
 *      The data probs are stored as doubles, but when they are read in
 *      they are converted to the log domain.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the pmf member function in the object.
 *
 *-----------------------------------------------------------------------
 */
void
CPT::read(iDataStreamFile& is)
{
  is.read(rows,"Discrete1DPDF::read, rows");
  if (rows <= 0)
    error("Discrete1DPDF: read rows (%d) < 0 in input",rows);
  is.read(cols,"Discrete1DPDF::read, cols");
  if (rows <= 0)
    error("Discrete1DPDF: read cols (%d) < 0 in input",cols);

  pmf.resize(rows*cols);

  logpr * ptr = pmf.ptr;
  for (int r=0;i<rows;i++) {
    for (int c=0;c<cols;c++) {
      double val;
      is.readDouble(val,"CPT::read, reading value");
      *ptr++ = val;
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * CPT::write(os)
 *      write out distribution to file 'os'. 
 *      The data probs are stored as doubles not in log domain.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effectcs other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
CPT::write(oDataStreamFile& os)
{
  is.write(rows,"Discrete1DPDF::write, rows");
  is.write(cols,"Discrete1DPDF::write, cols");

  logpr * ptr = pmf.ptr;
  for (int r=0;i<rows;i++) {
    for (int c=0;c<cols;c++) {
      is.writeDouble((*ptr).unlog(),"CPT::read, reading value");
      ptr++;
    }
  }
}




////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

