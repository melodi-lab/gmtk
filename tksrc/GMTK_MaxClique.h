/*
 * GMTK_MAXCLIQUE.h
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 *
 *   A clique class.
 *   Note: some texts define a 'clique' as any complete set
 *   while other texts define a 'clique' as a maximally
 *   complete set with respect to the subset operator (i.e., a
 *   clique is one such that no proper superset of the set
 *   of nodes is a clique). In order to avoid confusion,
 *   I adopt here the term 'maxclique' which corresponds
 *   to a maximally complete set. Note, however, that in this
 *   program, the concepts are such that 
 *
 *               'clique == maxclique != complete set'
 *
 *   meaning that cliques are taken to be max cliques.
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
 */

#ifndef GMTK_MAXCLIQUE_H
#define GMTK_MAXCLIQUE_H

#include "GMTK_RandomVariable.h"

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

class MaxClique
{
  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;

public:

  // the set of nodes which form a max clique
  set<RandomVariable*> nodes;

  MaxClique(set<RandomVariable*> arg) {
    nodes = arg;
  }
  ~MaxClique() {}

  float weight(const bool useDeterminism = true) const { 
    return computeWeight(nodes,NULL,useDeterminism); 
  }
  // static version of routine
  static float computeWeight(const set<RandomVariable*>& nodes,
			     const RandomVariable* node = NULL,
			     const bool useDeterminism = true);

  // complete the max clique
  void makeComplete() { makeComplete(nodes); }
  // static version of variable set completion, to
  // complete set of random variables passed in.
  static void makeComplete(set<RandomVariable*> &rvs);



};

#endif

