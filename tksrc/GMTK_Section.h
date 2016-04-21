/*
 * GMTK_Section.h
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 *
 * $Header$
 *
 */

#ifndef GMTK_SECTION_H
#define GMTK_SECTION_H

#include <vector>
#include <string>
#include <map>
#include <utility>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_RV.h"
#include "GMTK_FileParser.h"
#include "GMTK_MaxClique.h"

#include "debug.h"

// class mention for forward references.
class GraphicalModel;
class BoundaryTriangulate;
class Section;
class GMTemplate;

class Section : public IM {

  friend class GMTemplate;
  friend class BoundaryTriangulate;

public:

  // variables comprising this section.
  set<RV*> nodes;  

  // The cliques themselves, used to store the current triangulation
  // of each of the sections.
  vector<MaxClique> cliques;

  // a string with information about the method used to form the cliques
  string triMethod;



  // FIXME - this is section inference algorithm specific
  vector<unsigned> section_ri, section_li;
  map<string, vector<pair<unsigned, unsigned > > > ia_message_order;


  //  map<string, string> ia_name_to_section_inf_alg;
  map<string, map<string, unsigned> > clique_name_dictionary;



  Section() {}

  // Clone constructor from another Section, but that uses a new set
  // of random variables, and adjusts the frame of each new set of
  // random variable with offset

#if 0
  Section(Section& from_part,
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
  void writeInferenceArchitecture(oDataStreamFile& os, char section_name);
  void readMaxCliques(iDataStreamFile& is,  
		      string const &ia_name,
		      char section_type,
		      string const &section_inf_alg,
		      map< RVInfo::rvParent, RV* > &model_namePos2Var);
  void readInferenceArchitectureDefinition(iDataStreamFile &is,
					   string const &ia_name,
					   char section_type,
					   string const &section_inf_alg);
  void triangulateSectionsByCliqueCompletion();
  void setCliquesFromAnotherSection(Section& p);
  void reportScoreStats();

};


#endif

