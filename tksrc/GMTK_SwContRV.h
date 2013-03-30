/*
 * GMTK_SwContRV.h
 *
 *  Switching functionality for continuous observed RV.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
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

  ////////////////////////////////////////////////////////////////////////
  // Ties the parameters of 'this' with whatever those of 'other' are. 
  void tieParametersWith(SwContRV* other) {
    SwRV::tieParametersWith(other);
  }

  bool iterableSw() const { 
    return dtMapper->iterable();
  }



};

#endif
