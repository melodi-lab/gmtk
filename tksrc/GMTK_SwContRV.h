/*
 * GMTK_SwContRV.h
 *
 *  Switching functionality for continuous observed RV.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
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
 * The top level GMTK random variable object for the RV class hierarchy.
 *
 *
 *
 */

#ifndef GMTK_SW_CONT_RV_H
#define GMTK_SW_CONT_RV_H

#include <vector>
#include <string>
#include <set>

#include "logp.h"
#include "GMTK_RVInfo.h"
#include "GMTK_NamedObject.h"
#include "GMTK_RngDecisionTree.h"

#include "GMTK_SwRV.h"

/*
 * Note that this class is empty. It is here
 * only to maintain symmetry with the Discrete RV object
 * hierarchy.
 */
class SwContRV : public SwRV {

  friend class FileParser;

protected:

public:

  SwContRV() {}
  ~SwContRV() {}

};

#endif
