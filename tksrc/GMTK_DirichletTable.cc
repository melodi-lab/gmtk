/*-
 * GMTK_DirichletTable.cc
 *     General Dirichlet Table class (for storing Dirichlet priors for CPTs).
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2005, < fill in later >
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
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"

#include "GMTK_DirichletTable.h"


VCID("$Header$")


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * DirichletTable::DirichletTable()
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
DirichletTable::DirichletTable() 
{
}


/*-
 *-----------------------------------------------------------------------
 * DirichletTable::read(is)
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
DirichletTable::read(iDataStreamFile& is)
{
  NamedObject::read(is);

  is.read(_numDimensions,"Can't read Dirichlet Table number of dimensions");

  if (_numDimensions < 1) 
    error("ERROR: reading file '%s' line %d, DirichletTable '%s' trying to use non-positive (%d) num dimensions.",
	  is.fileName(),is.lineNo(),name().c_str(),_numDimensions);

  dimensionLengths.resize(_numDimensions);
  // read the dimension lengths, the last one corresponds to the "self cardinality"
  int numValues = 1;
  for (unsigned i=0;i<_numDimensions;i++) {
    is.read(dimensionLengths[i],"Can't read Dirichlet Table dimension length");
    if (dimensionLengths[i] <= 0)
      error("ERROR: reading file '%s' line %d, Dirichlet Table '%s' trying to use non-positive (%d) dimension length, position %d.",
	    is.fileName(),is.lineNo(),name().c_str(),dimensionLengths[i],i);
    numValues *= dimensionLengths[i];
  }

  table.resize(numValues);

  logpr* ptr = table.ptr;
  for (int i=0;i<numValues;i++) {
    double val;
    is.read(val,"DirichletTable::read, reading value");
    logpr pr(val);
    *ptr++ = pr;
  }
}




/*-
 *-----------------------------------------------------------------------
 * DirichletTable::write(os)
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
DirichletTable::write(oDataStreamFile& os)
{

  NamedObject::write(os);

  os.write(_numDimensions);
  for (unsigned i=0;i<_numDimensions;i++) {
    os.write(dimensionLengths[i]);
  }

  for (unsigned i=0;i<table.size();i++) {
    if (i % lastDimension() == 0)
      os.nl();
    os.write(table[i].unlog());
  }
  os.nl();
}



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////
