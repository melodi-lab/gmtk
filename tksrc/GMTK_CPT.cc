/*-
 * GMTK_CPT.cc
 *     Trainable (with say EM) CPT
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

#include "GMTK_CPT.h"


#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)


////////////////////////////////////////////////////////////////////
//        Static Data
////////////////////////////////////////////////////////////////////

unsigned CPT::warningNumParents = 50;

/////////////////////////////////////////////////////////////////
// threshold to ensure input data sums to approximately unity
double CPT::normalizationThreshold = 1e-2;

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * CPT::setNumParents()
 *      Just sets the number of parents and resizes arrays as appropriate.
 *
 * Results:
 *      no results.
 *
 * Side Effects:
 *      Will change internal arrays of this object.
 *
 *-----------------------------------------------------------------------
 */
void CPT::setNumParents(const unsigned _nParents)
{
  _numParents = _nParents;
  cardinalities.resize(_numParents);
}



/*-
 *-----------------------------------------------------------------------
 * CPT::setNumCardinality(var,card)
 *      sets the cardinality of var to card
 *      to set self, var == _numParents
 *
 * Results:
 *      no results.
 *
 * Side Effects:
 *      Will change internal array content of this object.
 *
 *-----------------------------------------------------------------------
 */
void CPT::setNumCardinality(const unsigned var, const int card)
{
  if (var > _numParents) 
    error("CPT: setNumCardinality, trying to use illegal (%d) var.",
	  var);
  if (card <= 0)
    error("CPT: setNumCardinality, trying to use illegal (%d) card.",
	  card);

  // assertion should be satisifed by the way that cardinalities
  // is allocated allong with setting num parents.
  if ( var == _numParents ) 
    _card = card;
  else { 
    assert ( var < cardinalities.size() );
    cardinalities[var] = card;
  }

}



/*-
 *-----------------------------------------------------------------------
 * Function
 *      compareCardinalities: compare the cardinalities of this CPT with that of an other. REturn
 *      true if they are equal false otherwise.
 *
 * Results:
 *      returns true if cards are equal.
 *
 * Side Effects:
 *      none
 *
 *-----------------------------------------------------------------------
 */
bool 
CPT::compareCardinalities(CPT& cpt)
{
  if (cardinalities.size() != cpt.cardinalities.size())
    return false;

  for (unsigned int i=0;i<cardinalities.size();i++) {
    if (cardinalities[i] != cpt.cardinalities[i])
      return false;
  }
  if (_card != cpt._card)
    return false;

  return true;
}


