/*
 * GMTK_NetworkFlowTriangulate.h
 *   The GMTK Network flow triangulation support routines.
 *
 * Written by Mukund Narasimhan <mukundn@ee.washington.edu> 
 *
 * Copyright (c) 2005, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 * $Header$
 */


#ifndef GMTK_NETWORKFLOWTRIANGULATE_H
#define GMTK_NETWORKFLOWTRIANGULATE_H

#include <string>
#include <vector>
#include <deque>
#include <set>
#include <map>
#include <numeric>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "GMTK_NetworkFlow.h"
#include "GMTK_RV.h"


//////////////////////////////////////////////////////////////////////
//  The GMTK2Network class provides routines that translate the     //
//  boundary triangulation problem (for some objective functions    //
//  such as the number of vertices in the interface or the sum of   //
//  weights of the vertices in the interface) into a flow problem.  //
//////////////////////////////////////////////////////////////////////
class GMTK2Network : public networkFlow::Network {
 public:
  // Argument Summary.
  // computeWeight    : a function that computes the weight of each  
  //                    set of vertices. It will only be called with
  //                    one element sets, but a generic calling 
  //                    interface is used for future changes.
  // cLeft            : The set of RVs that constitute the left interface.
  // cMiddle          : The set of all the RVs
  // cRight           : The set of RVs that constitute the right interface.
  GMTK2Network(float (*computeWeight)(const set<RV*>&), 
	       const std::set<RV*>& cLeft, 
	       const std::set<RV*>& cMiddle, 
	       const std::set<RV*>& cRight) : 
    _cLeft(cLeft), _cRight(cRight), _cMiddle(cMiddle) {
    init(computeWeight, cLeft, cMiddle, cRight);
  }
  

  
  void init(float (*computeWeight)(const set<RV*>&), 
	    const std::set<RV*>& cLeft, 
	    const std::set<RV*>& cMiddle, 
	    const std::set<RV*>& cRight);
	
  // This routine finds the best interface and this interface is returned
  // in the variable boundary. All vertices to the "left" of the bounary
  // are returned in preBounary.
  void findBoundary(std::set<RV*>& boundary, std::set<RV*>& preBoundary, bool forward=true);
    
 protected:
  const std::set<RV*>& _cLeft;
  const std::set<RV*>& _cRight;
  const std::set<RV*>& _cMiddle;
  std::map<RV*, int> _translation;
  std::map<int, RV*> _reverseTranslation;
  
  int inVertex(RV* rv) {
    if (_translation.count(rv) != 0)
      return _translation[rv]+1;
    else
      return -1;
  }

  int outVertex(RV* rv) {
    if (_translation.count(rv) != 0)
      return _translation[rv];
    else
      return -1;
  }

  bool isInVertex(int i) const {
    return ((i%2) == 1);
  }

  bool isOutVertex(int i) const {
    return ((i%2) == 0);
  }
};
#endif
