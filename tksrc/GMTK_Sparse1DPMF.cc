/*-
 * GMTK_Sparse1DPMF.cc
 *     Trainable (with say EM) sparse 1D discrete probability
 *     distributions.
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
 * Sparse1DPMF::Sparse1DPMF()
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
Sparse1DPMF::Sparse1DPMF() 
{

}



/*-
 *-----------------------------------------------------------------------
 * Sparse1DPMF::read(is)
 *      read in a distribution from file 'is'. 
 *      The data probs are stored as doubles, but when they are read in
 *      they are converted into the log domain (so discrete data on disk
 *      is NOT represented as log probabilities. This is because 
 *          1) discrete data typically doesn't need such a huge dynamic range
 *              (like Gaussian probabilties do).
 *          2) it is easier to examine data on disk when it is not in log domain.
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
Sparse1DPMF::read(iDataStreamFile& is)
{
  int len;
  is.read(card,"Sparse1DPMF::read, card");
  if (card <= 0)
    error("Sparse1DPMF: read length (%d) < 0 in input",card);

  is.read(len,"Sparse1DPMF::read, len");
  if (len <= 0)
    error("Sparse1DPMF: read length (%d) < 0 in input",len);
  if (len > card)
    error("Sparse1DPMF: read length (%d) > card (%d) in input",len,card);

  pmf.resize(len);

  for (int i=0;i<length;i++) {
    int val;
    double prob;

    is.read(val,"Sparse1DPMF::read, reading value");
    if (val < 0 || val > card-1)
      error("Sparse1DPMF::read, bad value = %d, must be in range [0:%d]",val,
	    card-1);

    is.readDouble(prob,"Sparse1DPMF::read, reading prob");
    if (prob < 0.0 || prob > 1.0)
      error("Sparse1DPMF: read, invalid pmf value (%g)",val);
    pmf[i].val = val;
    pmf[i].prob = prob;
  }
}




/*-
 *-----------------------------------------------------------------------
 * Sparse1DPMF::write(is)
 *      write out distribution to file 'os'. 
 *      the data probs are stored on disk as doubles,  NOT in log domain.
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
Sparse1DPMF::write(oDataStreamFile& os)
{
  os.write(card,"Sparse1DPMF::write, card");
  os.write(pmf.len(),"Sparse1DPMF::write, len");
  for (int i=0;i<pmf.len();i++) {
    os.write(pmf[i].val,"Sparse1DPMF::write, writing value");
    os.writeDouble(pmf[i].prob,"Sparse1DPMF::write, writing prob");
  }
  os.nl();
}





////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

