/*-
 * GMTK_DirichletPrior.h
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

#ifndef GMTK_DIRICHLETPRIOR_H
#define GMTK_DIRICHLETPRIOR_H

#include "sArray.h"
#include "GMTK_DirichletTable.h"

class DirichletPrior {

protected:

  enum PriorType { NoneVal, DirichletConstVal, DirichletTableVal };
  PriorType smoothingType;
  // in case smoothingType == DirichletConstVal, then the constant value
  // corresponding to the one hyperparameter of the Dirichlet prior.
  double dirichletAlpha;
  // in case smoothingType == DirichletTableVal, then the pointer to the
  // table of DirichletTable hyperparameters.
  DirichletTable* dirichletTable;

  // strings for the master files
  static char* DirichletConstStr;
  static char* DirichletTableStr;

  DirichletPrior() {
    smoothingType = NoneVal;
    dirichletAlpha = 0.0;
    dirichletTable = NULL;
  }

};





#endif 
