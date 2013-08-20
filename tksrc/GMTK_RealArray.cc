/*-
 * GMTK_RealArray.cc
 *     General real array (for means, etc.)
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"

#include "GMTK_RealArray.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)



////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * RealArray::RealArray()
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
RealArray::RealArray() 
{


}


/*-
 *-----------------------------------------------------------------------
 * RealArray::read(is)
 *      read in the array from file 'is'. 
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
RealArray::read(iDataStreamFile& is)
{
  int length;

  is.read(length,"Can't read array length");
  if (length <= 0)
    error("RealArray: read length (%d) < 0 in input",length);

  resize(length);

  for (int i=0;i<length;i++) {
    is.read(operator[](i),"Can't read array value");
  }
}




/*-
 *-----------------------------------------------------------------------
 * RealArray::write(os)
 *      write out data to file 'os'. 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effects other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
RealArray::write(oDataStreamFile& os)
{

  os.write(len(),"RealArray::write, length");
  for (int i=0;i<len();i++) {
    os.write(operator[](i),"RealArray::write, values");
  }
  os.nl();
}



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////
