/*
 * GMTK_ScPnShRV.cc
 *
 * Switching support functionality for random variables with switching parents.
 * 
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 *
 * The top level GMTK random variable object for the RV class hierarchy.
 *
 *
 *
 */

#include "GMTK_ScPnShRV.h"
#include "GMTK_RV.h"
#include "GMTK_RVInfo.h"



/*-
 *-----------------------------------------------------------------------
 * printSelf()
 *      prints a one-line summary of the detailed information about this RV.
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
void ScPnShRV::printSelf(RVInfo::WeightInfo& wi,FILE *f,bool nl)
{
    if (wi.penalty.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      fprintf(f,"penalty = %f, ", wi.penalty.weight_value);
    } else if (wi.penalty.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      fprintf(f,"penalty = %d:%d, ",
	      wi.penalty.firstFeatureElement,
	      wi.penalty.firstFeatureElement);
    }
    if (wi.scale.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      fprintf(f,"scale = %f, ", wi.scale.weight_value);
    } else if (wi.scale.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      fprintf(f,"scale = %d:%d, ",
	      wi.scale.firstFeatureElement,
	      wi.scale.firstFeatureElement);
    }
    if (wi.shift.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      fprintf(f,"shift = %f, ", wi.shift.weight_value);
    } else if (wi.shift.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      fprintf(f,"shift = %d:%d, ",
	      wi.shift.firstFeatureElement,
	      wi.shift.firstFeatureElement);
    }
    if (nl) fprintf(f,"\n");
}



/*-
 *-----------------------------------------------------------------------
 * printSelfVerbose()
 *      prints a multi-line verbose description of this RV.
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
void ScPnShRV::printSelfVerbose(RVInfo::WeightInfo& wi,FILE *f)
{
  printSelf(wi,f);
}

