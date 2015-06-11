/*
 * GMTK_ScPnShV.h
 *
 *  Scale/Penalty/Shift functionality for a RV.
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

/*
  Support shift, scale, and penalty. 
  given a probabilty p, we to modify it, we do: penalty*p^scale+shift
  Or logged:
         log(penalty*p^scale+shift)
      =  log(penalty*p^scale) ++ log(shift)
      =  (log(penalty) + scale*log(p)) ++ log(shift)
  
  Example syntax:
       weight: 
         scale 1.0 0:0 shift 0.5 ;
       | scale 0:0 penalty 1:1
       | scale 2:2
       | nil;

*/



#ifndef GMTK_SC_PN_SH_RV_H
#define GMTK_SC_PN_SH_RV_H

#include <vector>
#include <string>
#include <set>

#include "logp.h"
#include "GMTK_RVInfo.h"
#include "GMTK_NamedObject.h"
#include "GMTK_RngDecisionTree.h"
#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_FileSource.h"
#endif
#include "GMTK_RV.h"

class ScPnShRV {

  friend class FileParser;

protected:

public:

  ScPnShRV() {}
  ~ScPnShRV() {}

  // For ticket #6: returns true iff it's "safe" to call modifyProbability()
  // with the specified WeightInfo at this time. Safe means that all the
  // information needed to modify the probability is currently available.
  // In particular, if data from the globalObservationMatrix is required, the
  // globalObservationMatrix must be non-NULL and ready to provide data.

  // modifyProbability() can be called early in the course of GMTK setup to
  // prepare for cpbeam pruning, before the globalObservationMatrix has been
  // instantiated. That is why this check is necessary.

  bool safeToModifyProbability(RVInfo::WeightInfo &wi) {
    if (wi.penalty.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      // this is safe; the penalty is known
    } else if (wi.penalty.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      if (!globalObservationMatrix) return false; // need observations, but they're not here yet
      if (!globalObservationMatrix->active()) return false; // GOM exists, but not ready yet
    }
    if (wi.scale.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      // this is safe; the scale is known
    } else if (wi.scale.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      if (!globalObservationMatrix) return false; // need observations, but they're not here yet
      if (!globalObservationMatrix->active()) return false; // GOM exists, but not ready yet
    }
    if (wi.shift.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      // this is safe; the shift is known
    } else if (wi.shift.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      if (!globalObservationMatrix) return false; // need observations, but they're not here yet
      if (!globalObservationMatrix->active()) return false; // GOM exists, but not ready yet
    }
    return true; // everything we need is available
  }

  // Version with branches.
  // Given a p, modify p according to:
  //         penalty*p^scale+shift
  // or alternatively:
  //  (log(penalty) + scale*log(p)) ++ log(shift)
  // where ++ is the log_add operator.
  inline void modifyProbability(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    if (wi.scale.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      p.valref() *= wi.scale.weight_value;
    } else if (wi.scale.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      p.valref() *= 
	(*(globalObservationMatrix->floatVecAtFrame(rv->frame(), 
						    wi.scale.firstFeatureElement)));
    }
    if (wi.penalty.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      p.valref() += wi.penalty.weight_value;
    } else if (wi.penalty.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      // TODO: check here that global obsevation matrix is active.
      p.valref() += 
	(*(globalObservationMatrix->floatVecAtFrame(rv->frame(),
						    wi.penalty.firstFeatureElement)));
    }
    if (wi.shift.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      logpr shift((void*)NULL,wi.shift.weight_value);
      p += shift;
    } else if (wi.shift.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      logpr shift((void*)NULL,
		  (*(globalObservationMatrix->floatVecAtFrame(rv->frame(), 
							      wi.shift.firstFeatureElement))));
      p += shift;
    }
  }

  // Versions without branches.
  // CP = constant penalty
  // OP = observed penalty
  // CS = constant scale
  // OS = observed scale
  // CO = constant offset
  // OO = observed offset

  // 3 values, 3 options (nothing, const, obsered) = 27 routines (actually
  // 26 since we don't include the do nothing routine).

  inline void modifyProbabilityCP(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    p.valref() += wi.penalty.weight_value;
  }
  inline void modifyProbabilityOP(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
      p.valref() += 
	(*(globalObservationMatrix->floatVecAtFrame(rv->frame(),
						    wi.penalty.firstFeatureElement)));
  }
  inline void modifyProbabilityCS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
      p.valref() *= wi.scale.weight_value;
  }
  inline void modifyProbabilityOS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
      p.valref() *= 
	(*(globalObservationMatrix->floatVecAtFrame(rv->frame(), 
						  wi.scale.firstFeatureElement)));
  }
  inline void modifyProbabilityCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
      logpr shift((void*)NULL,wi.shift.weight_value);
      p += shift;
  }
  inline void modifyProbabilityOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) { 
      logpr shift((void*)NULL,
		  (*(globalObservationMatrix->floatVecAtFrame(rv->frame(), 
							      wi.shift.firstFeatureElement))));
      p += shift;
  }

  inline void modifyProbabilityCPCSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityCPCSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }
  inline void modifyProbabilityCPOSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityCPOSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }
  inline void modifyProbabilityOPCSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityOPCSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }
  inline void modifyProbabilityOPOSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityOPOSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }



  inline void modifyProbabilityCPCS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityCP(p,wi,rv);
  }
  inline void modifyProbabilityCPOS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityCP(p,wi,rv);
  }
  inline void modifyProbabilityOPCS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityOP(p,wi,rv);
  }
  inline void modifyProbabilityOPOS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityOP(p,wi,rv);
  }


  inline void modifyProbabilityCPCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityCPOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }
  inline void modifyProbabilityOPCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityOPOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }

  inline void modifyProbabilityCSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityCSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }
  inline void modifyProbabilityOSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityOSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }


  void printSelf(RVInfo::WeightInfo& wi,FILE *f,bool nl=true);
  void printSelfVerbose(RVInfo::WeightInfo& wi,FILE *f);

};



#endif
