/*
 * GMTK_GraphicalModel.h
 * General graphical model class
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
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

#include "GMTK_RandomVariable.h"
#include "GMTK_GM.h" 
#include "GMTK_CPT.h"
#include "GMTK_MixtureCommon.h"

class RandomVariable;

class GraphicalModel
{
private:



  struct frame {
    // list of random variables for this frame
    vector <RandomVariable*> rvs;
  };

  // the set of frames for this graphical model
  vector < frame > frames;


  static bool topologicalSortRecurse(vector<RandomVariable*>& outputVarList,
				     RandomVariable* node,
				     unsigned& position);

  static bool topologicalSortRecurse(const set<RandomVariable*>& sortSet,
				     vector<RandomVariable*>& outputVarList,
				     RandomVariable* node,
				     unsigned& position);

  static bool topologicalSortRecurseRandom(const set<RandomVariable*>& sortSet,
					   vector<RandomVariable*>& outputVarList,
					   RandomVariable* node,
					   unsigned& position);

  static bool topologicalSortRecurseContFirst(const set<RandomVariable*>& sortSet,
					      vector<RandomVariable*>& outputVarList,
					      RandomVariable* node,
					      unsigned& position);


public:

  GraphicalModel() {}
  ~GraphicalModel() {}

  static bool topologicalSort(vector<RandomVariable*> &inputVarList,
			      vector<RandomVariable*> &outputVarList);



  // topoligical sort where 
  //  1) input is a set
  //  2) constrain sort to be only variables in sortSet (which
  //     might be a subset of the ancestral or 'descendal' set.
  static bool topologicalSort(const set<RandomVariable*> &inputVarList,
			      const set<RandomVariable*> &sortSet,
			      vector<RandomVariable*> &outputVarList);


  // topoligical sort where 
  //  1) input is a set
  //  2) constrain sort to be only variables in sortSet (which
  //     might be a subset of the ancestral or 'descendal' set.
  //  3) we do a random order
  static bool topologicalSortRandom(const set<RandomVariable*> &inputVarList,
				    const set<RandomVariable*> &sortSet,
				    vector<RandomVariable*> &outputVarList);

  // A version of topological sort that places the continuous
  // variables as *early* in the sort as possible, all other things
  // being equal.
  static bool topologicalSortContFirst(const set<RandomVariable*> &inputVarList,
				       const set<RandomVariable*> &sortSet,
				       vector<RandomVariable*> &outputVarList,
				       const string priorityStr = "COB");

  ///////////////////////////////////////////
  // print this graphical model
  void print();
};

#endif
