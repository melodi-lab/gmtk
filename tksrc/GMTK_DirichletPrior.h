/*-
 * GMTK_DirichletPrior.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
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
  static const char* DirichletConstStr;
  static const char* DirichletTableStr;

  DirichletPrior() {
    smoothingType = NoneVal;
    dirichletAlpha = 0.0;
    dirichletTable = NULL;
  }

};





#endif 
