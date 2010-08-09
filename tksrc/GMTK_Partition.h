/*
 * GMTK_GMTemplate.h
 *   Basic GM Template and Basic Triangulation Routines.
 *   This includes code that is common to both triangulation and inference,
 *   so does not contain the more elaborate triangulation methods so that
 *   they are not appart of inference code.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * $Header$
 *
 */

#ifndef GMTK_PARTITION_H
#define GMTK_PARTITION_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_RV.h"
#include "GMTK_FileParser.h"
#include "GMTK_MaxClique.h"

#include "debug.h"

// class mention for forward references.
class GraphicalModel;
class BoundaryTriangulate;
class Partition;
class GMTemplate;

class Partition : public IM {

  friend class GMTemplate;
  friend class BoundaryTriangulate;

public:

  // variables comprising this partition.
  set<RV*> nodes;  

  // The cliques themselves, used to store the current triangulation
  // of each of the partitions.
  vector<MaxClique> cliques;
  
  // a string with information about the method used to form the cliques
  string triMethod;

  Partition() {}

  // Clone constructor from another Partition, but that uses a new set
  // of random variables, and adjusts the frame of each new set of
  // random variable with offset

#if 0
  Partition(Partition& from_part,
	    vector <RV*>& newRvs,
	    map < RVInfo::rvParent, unsigned >& ppf,
	    const unsigned int frameDelta = 0);
#endif

  void clear() { 
    set<RV*>::iterator it;
    for (it = nodes.begin();it != nodes.end();it++) 
      delete (*it);
    nodes.clear(); 
    cliques.clear(); 
    triMethod.clear(); 
  }

  void clearCliques() { cliques.clear(); triMethod.clear(); }

  void writeMaxCliques(oDataStreamFile& os);  
  void readMaxCliques(iDataStreamFile& is);
  void triangulatePartitionsByCliqueCompletion();
  void setCliquesFromAnotherPartition(Partition& p);
  void reportScoreStats();

};


#endif

