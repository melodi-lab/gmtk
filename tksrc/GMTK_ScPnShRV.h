/*
 * GMTK_ScPnShV.h
 *
 *  Scale/Penalty/Shift functionality for a RV.
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

#ifndef GMTK_SC_PN_SH_RV_H
#define GMTK_SC_PN_SH_RV_H

#include <vector>
#include <string>
#include <set>

#include "logp.h"
#include "GMTK_RVInfo.h"
#include "GMTK_NamedObject.h"
#include "GMTK_RngDecisionTree.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_RV.h"

class ScPnShRV {

  friend class FileParser;

protected:

public:

  ScPnShRV() {}
  ~ScPnShRV() {}

  // Version with branches.
  // Given a p, modify p according to:
  //         scale*p^penalty+shift
  // or alternatively:
  //  (log(scale) + penalty*log(p)) ++ log(shift)
  // where ++ is the log_add operator.
  inline void modifyProbability(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    if (wi.penalty.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      p.valref() *= wi.penalty.weight_value;
    } else if (wi.penalty.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      p.valref() *= 
	(*globalObservationMatrix.floatVecAtFrame(rv->frame(),
						  wi.penalty.firstFeatureElement));
    }
    if (wi.scale.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      p.valref() += wi.scale.weight_value;
    } else if (wi.scale.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      p.valref() += 
	(*globalObservationMatrix.floatVecAtFrame(rv->frame(), 
						  wi.scale.firstFeatureElement));
    }
    if (wi.shift.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Constant) {
      logpr shift((void*)NULL,wi.shift.weight_value);
      p += shift;
    } else if (wi.shift.wt_Status == RVInfo::WeightInfo::WeightItem::wt_Observation) {
      logpr shift((void*)NULL,
		  (*globalObservationMatrix.floatVecAtFrame(rv->frame(), 
							    wi.shift.firstFeatureElement)));
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
    p.valref() *= wi.penalty.weight_value;
  }
  inline void modifyProbabilityOP(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
      p.valref() *= 
	(*globalObservationMatrix.floatVecAtFrame(rv->frame(),
						  wi.penalty.firstFeatureElement));
  }
  inline void modifyProbabilityCS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
      p.valref() += wi.scale.weight_value;
  }
  inline void modifyProbabilityOS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
      p.valref() += 
	(*globalObservationMatrix.floatVecAtFrame(rv->frame(), 
						  wi.scale.firstFeatureElement));
  }
  inline void modifyProbabilityCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
      logpr shift((void*)NULL,wi.shift.weight_value);
      p += shift;
  }
  inline void modifyProbabilityOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) { 
      logpr shift((void*)NULL,
		  (*globalObservationMatrix.floatVecAtFrame(rv->frame(), 
							    wi.shift.firstFeatureElement)));
      p += shift;
  }

  inline void modifyProbabilityCPCSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityCPCSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }
  inline void modifyProbabilityCPOSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityCPOSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }
  inline void modifyProbabilityOPCSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityOPCSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityCS(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }
  inline void modifyProbabilityOPOSCO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityCO(p,wi,rv);
  }
  inline void modifyProbabilityOPOSOO(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityOS(p,wi,rv);
    modifyProbabilityOO(p,wi,rv);
  }



  inline void modifyProbabilityCPCS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityCS(p,wi,rv);
  }
  inline void modifyProbabilityCPOS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityCP(p,wi,rv);
    modifyProbabilityOS(p,wi,rv);
  }
  inline void modifyProbabilityOPCS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityCS(p,wi,rv);
  }
  inline void modifyProbabilityOPOS(logpr& p,RVInfo::WeightInfo& wi,RV* rv) {
    modifyProbabilityOP(p,wi,rv);
    modifyProbabilityOS(p,wi,rv);
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
