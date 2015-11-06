/*
 * GMTK_MAXCLIQUE.h
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
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
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */

// TODO: perhaps create a subclass or member of maxClique at some point, rather than
// adding everything for exact inference to the base class.


#ifndef GMTK_INFERENCEFACTORCLIQUE_H
#define GMTK_INFERENCEFACTORCLIQUE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "general.h"
#include "cArray.h"
#include "sArray.h"
#include "debug.h"

#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_FactorClique.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <map>

class InferenceFactorClique : public IM
{

  // Non-STL "f=fast" versions of arrays that are instantiated
  // only in the final unrolled versions of the separator cliques.
  sArray< RV*> fOrderedNodes;

  // the original factor clique from which this object has been cloned.
  FactorClique& origin;


public:

  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  InferenceFactorClique() : origin(*((FactorClique*)NULL)) {}
  // normal (or re-)constructor
  InferenceFactorClique(FactorClique& _origin,
			vector <RV*>& newRvs,
			map < RVInfo::rvParent, unsigned >& ppf,
			const unsigned int frameDelta);
  // version for VE separators
  InferenceFactorClique(FactorClique& _origin);

};





#endif
