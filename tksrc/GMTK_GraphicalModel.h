/*
 * GMTK_GraphicalModel.h
 * General graphical model class
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 * $Header$
 *
 */

#ifndef GMTK_GRAPHICALMODEL_H
#define GMTK_GRAPHICALMODEL_H

#include <vector>
#include <string>
#include <map>
#include <set>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_RV.h"
#include "GMTK_CPT.h"
#include "GMTK_MixtureCommon.h"

class RV;

class GraphicalModel
{
private:


  static bool topologicalSortRecurse(vector<RV*>& outputVarList,
				     RV* node,
				     unsigned& position,
				     map<RV*,unsigned>&);

  static bool topologicalSortRecurse(const set<RV*>& sortSet,
				     vector<RV*>& outputVarList,
				     RV* node,
				     unsigned& position,
				     map<RV*,unsigned>&);

  static bool topologicalSortRecurseRandom(const set<RV*>& sortSet,
					   vector<RV*>& outputVarList,
					   RV* node,
					   unsigned& position,
					   map<RV*,unsigned>&);

  static bool topologicalSortRecurseWPriorityRecurse(const set<RV*>& sortSet,
						     vector<RV*>& outputVarList,
						     RV* node,
						     unsigned& position,
						     map<RV*,unsigned>&);


public:

  GraphicalModel() {}
  ~GraphicalModel() {}

  static bool topologicalSort(vector<RV*> &inputVarList,
			      vector<RV*> &outputVarList);



  // topoligical sort where 
  //  1) input is a set
  //  2) constrain sort to be only variables in sortSet (which
  //     might be a subset of the ancestral or 'descendal' set.
  static bool topologicalSort(const set<RV*> &inputVarList,
			      const set<RV*> &sortSet,
			      vector<RV*> &outputVarList);


  // topoligical sort where 
  //  1) input is a set
  //  2) constrain sort to be only variables in sortSet (which
  //     might be a subset of the ancestral or 'descendal' set.
  //  3) we do a random order
  static bool topologicalSortRandom(const set<RV*> &inputVarList,
				    const set<RV*> &sortSet,
				    vector<RV*> &outputVarList);

  // A version of topological sort that places the continuous
  // variables as *early* in the sort as possible, all other things
  // being equal.
  static bool topologicalSortWPriority(const set<RV*> &inputVarList,
				       const set<RV*> &sortSet,
				       vector<RV*> &outputVarList,
				       const string priorityStr = "COB");

  ///////////////////////////////////////////
  // print this graphical model
  void print();
};

#endif
