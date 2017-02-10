/*-
 * GMTK_Section.cc
 *    Basic Section for a given graph file.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2009 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */


#if HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_HG_H
#  include "hgstamp.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <ctype.h>

#include <iterator>
#include <map>
#include <set>
#include <algorithm>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_Section.h"

VCID(HGID)


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors/Destructors 
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

#if 0
// This constructor isn't being used currently. as of Wed Jul 13, 2005
// TODO: eventually remove
Section::Section(Section& from_part,
		     vector <RV*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta)
{

  triMethod = from_part.triMethod;

  set<RV*>::iterator it;

  // clone over nodes RVs.  
  // TODO: make this next code a routine
  //  nodesClone() since it is used in several places.
  for (it = from_part.nodes.begin();
       it != from_part.nodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
    }

    RV* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }
  cliques.reserve(from_part.cliques.size());
  // 
  // NOTE: It is Crucial for the cliques in the cloned section to be
  // inserted in the *SAME ORDER* as in the section being cloned.
  for (unsigned i=0;i<from_part.cliques.size();i++) {
    cliques.push_back(MaxClique(from_part.cliques[i],
				newRvs,ppf,frameDelta));
  }
}
#endif

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Section support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



void
writeRVSet(oDataStreamFile& os, set<RV*> const &rvs, char const *comment) {
  // assigned prob nodes n (name frame)^n
  os.write(rvs.size()); 
  for (set<RV*>::iterator j=rvs.begin(); j != rvs.end(); j++) {
    RV* rv = (*j);
    os.write(rv->name().c_str());
    os.write(rv->frame());
  }
  os.writeComment(comment);
  os.nl();
}

void
Section::
writeInferenceArchitecture(oDataStreamFile& os, char section_name) {
  // First write out the cliques in commented form for the user

  // TODO: convert tri method string to have no spaces/tabs.


  os.writeComment("---- %d total max-cliques \n",cliques.size());
  if (triMethod.size() > 0)
    os.writeComment("---- Triangulation came from method: %s\n",triMethod.c_str());


  double maxWeight = -1.0;
  double totalWeight = -1.0; // starting flag
  for (unsigned i=0;i<cliques.size();i++) {
    double curWeight = cliques[i].weight();
    if (curWeight > maxWeight) maxWeight = curWeight;
    if (totalWeight == -1.0)
      totalWeight = curWeight;
    else
      totalWeight = log10add(curWeight,totalWeight);
    os.writeComment("%d : %d  %f\n",
		    i,
		    cliques[i].nodes.size(),curWeight);
    for (set<RV*>::iterator j=cliques[i].nodes.begin();
	 j != cliques[i].nodes.end(); j++) {
      RV* rv = (*j);
      os.writeComment("   %s(%d)\n",rv->name().c_str(),rv->frame());
    }
  }
  os.writeComment("Maximum clique state space = 1e%f, total state space = 1e%f\n",maxWeight,totalWeight);
  // Then write out the same information in a less human-readable but more machine
  // readable format.


  os.write("default "); os.write(section_name); os.write("SPARSE_JOIN"); os.nl();
  if (triMethod.size() > 0) {
    // remove all white space to turn into string token.
    for (unsigned i=0;i<triMethod.size();i++) {
      if (isspace(triMethod[i]))
	triMethod[i] = '_';
    }
    os.write(triMethod.c_str());
  } else {
    os.write("UNKNOWN_TRIANGULATION_METHOD");
  }
  os.nl();

  os.write(cliques.size()); // number of cliques
  os.nl();
  for (unsigned i=0;i<cliques.size();i++) {
    os.write(i); // clique number i
    char buf[16];
    sprintf(buf, "C%03u ", i);
    os.write(buf);
    os.write(cliques[i].nodes.size());  // number of nodes in clique number i
    for (set<RV*>::iterator j=cliques[i].nodes.begin();
	 j != cliques[i].nodes.end(); j++) {
      RV* rv = (*j);
      os.write(rv->name().c_str());
      os.write(rv->frame());
    }
    os.nl();

    // sorted assigned nodes n (name frame disposition)^n
    os.write(cliques[i].sortedAssignedNodes.size());
    for (unsigned j=0; j < cliques[i].sortedAssignedNodes.size(); ++j) {
      os.write(cliques[i].sortedAssignedNodes[j]->name());
      os.write(cliques[i].sortedAssignedNodes[j]->frame());
      os.write((int)cliques[i].dispositionSortedAssignedNodes[j]);
    }
    os.writeComment("sorted assigned (name frame disposition)");
    os.nl();

    // assigned prob nodes n (name frame)^n
    writeRVSet(os, cliques[i].assignedProbNodes, "assigned prob nodes");
    writeRVSet(os, cliques[i].cumulativeAssignedNodes, "cumulative assigned nodes");
    writeRVSet(os, cliques[i].cumulativeAssignedProbNodes, "cumulative assigned prob nodes");
    writeRVSet(os, cliques[i].unionIncommingCESeps, "union incomming seps");
    writeRVSet(os, cliques[i].unassignedIteratedNodes, "unassigned iterated");
    writeRVSet(os, cliques[i].cumulativeUnassignedIteratedNodes, "cumulative unassigned iterated");
    os.writeComment(""); os.nl();
  }
}


/*-
 *-----------------------------------------------------------------------
 * Section::writeMaxCliques()
 *   Write out the max cliques of the given sections.
 *
 * Preconditions:
 *   The maxclique variable must be instantiated.
 *
 * Postconditions:
 *   Information about the maxclique variable is written out.
 *
 * Side Effects:
 *   Moves file pointer
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
Section::
writeMaxCliques(oDataStreamFile& os)
{
  // First write out the cliques in commented form for the user

  // TODO: convert tri method string to have no spaces/tabs.


  os.writeComment("---- %d total max-cliques \n",cliques.size());
  if (triMethod.size() > 0)
    os.writeComment("---- Triangulation came from method: %s\n",triMethod.c_str());


  double maxWeight = -1.0;
  double totalWeight = -1.0; // starting flag
  for (unsigned i=0;i<cliques.size();i++) {
    double curWeight = cliques[i].weight();
    if (curWeight > maxWeight) maxWeight = curWeight;
    if (totalWeight == -1.0)
      totalWeight = curWeight;
    else
      totalWeight = log10add(curWeight,totalWeight);
    os.writeComment("%d : %d  %f\n",
		    i,
		    cliques[i].nodes.size(),curWeight);
    for (set<RV*>::iterator j=cliques[i].nodes.begin();
	 j != cliques[i].nodes.end(); j++) {
      RV* rv = (*j);
      os.writeComment("   %s(%d)\n",rv->name().c_str(),rv->frame());
    }
  }
  os.writeComment("Maximum clique state space = 1e%f, total state space = 1e%f\n",maxWeight,totalWeight);
  // Then write out the same information in a less human-readable but more machine
  // readable format.

  if (triMethod.size() > 0) {
    // remove all white space to turn into string token.
    for (unsigned i=0;i<triMethod.size();i++) {
      if (isspace(triMethod[i]))
	triMethod[i] = '_';
    }
    os.write(triMethod.c_str());
  } else {
    os.write("UNKNOWN_TRIANGULATION_METHOD");
  }
  os.nl();

  os.write(cliques.size()); // number of cliques
  os.nl();
  for (unsigned i=0;i<cliques.size();i++) {
    os.write(i); // clique number i
    os.write(cliques[i].nodes.size());  // number of nodes in clique number i
    for (set<RV*>::iterator j=cliques[i].nodes.begin();
	 j != cliques[i].nodes.end(); j++) {
      RV* rv = (*j);
      os.write(rv->name().c_str());
      os.write(rv->frame());
    }
    os.nl();
  }
}



string
makeDictionaryKey(string const &ia_name, char section_type) {
  return string(ia_name + ":" + section_type);
}

RV *
getRVfromIS(iDataStreamFile &is, 
	    map < RVInfo::rvParent, RV* > &namePos2Var,
	    char const *member_name, 
	    unsigned clique_index, unsigned rv_index)
{
  RVInfo::rvParent par;
  is.read(par.first,"RV name");
  is.read(par.second,"RV frame");
  
  map < RVInfo::rvParent, RV* >::iterator loc;
  loc = namePos2Var.find(par);
  if (loc == namePos2Var.end())
    error("ERROR: reading file %s line %d, %s %d has %d'th variable %s(%d) that does not exist.\n",
	  is.fileName(), is.lineNo(), member_name, clique_index, rv_index, par.first.c_str(), par.second);
  return (*loc).second;
}


RV *
getRVfromIS(iDataStreamFile &is, 
      map < RVInfo::rvParent, RV* > &namePos2Var,
      char const *member_name, 
      unsigned clique_index, unsigned rv_index,
      Json::Value const var_frame)
{
  RVInfo::rvParent par;
  par.first = var_frame[0].asString();
  par.second = var_frame[1].asInt();
  // is.read(par.first,"RV name");
  // is.read(par.second,"RV frame");
  
  map < RVInfo::rvParent, RV* >::iterator loc;
  loc = namePos2Var.find(par);
  if (loc == namePos2Var.end())
    error("ERROR: reading file %s line %d, %s %d has %d'th variable %s(%d) that does not exist.\n",
    is.fileName(), is.lineNo(), member_name, clique_index, rv_index, par.first.c_str(), par.second);
  return (*loc).second;
}



void
readRVSetFromIS(iDataStreamFile &is,
		char const *err_msg, char const *member_name,
		unsigned clique_idx, 
		map < RVInfo::rvParent, RV* > &namePos2Var,
		set<RV*> &rv_set)
{
  unsigned num_rvs;
  is.read(num_rvs, err_msg);
  for (unsigned j=0; j < num_rvs; ++j) {
    RV *rv = getRVfromIS(is, namePos2Var, member_name, clique_idx, j);
    rv_set.insert(rv);
  }
}



void
readRVSetFromIS(iDataStreamFile &is,
    char const *err_msg, char const *member_name,
    unsigned clique_idx, 
    map < RVInfo::rvParent, RV* > &namePos2Var,
    set<RV*> &rv_set,
    Json::Value const var_frame)
{
  unsigned num_rvs = var_frame.size();
  // is.read(num_rvs, err_msg);
  for (unsigned j=0; j < num_rvs; ++j) {
    RV *rv = getRVfromIS(is, namePos2Var, member_name, clique_idx, j, var_frame[j]);
    rv_set.insert(rv);
  }
}


/*-
 *-----------------------------------------------------------------------
 * Section::readMaxCliques()
 *   Write out the max cliques of the given sections.
 *
 * Preconditions:
 *   The maxclique variable must be instantiated.
 *
 * Postconditions:
 *   Information about the maxclique variable is written out.
 *
 * Side Effects:
 *   Moves file pointer
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
Section::
readMaxCliques(iDataStreamFile& is,
	       string const &ia_name,
	       char section_type,
	       string const &section_inf_alg,
	       map< RVInfo::rvParent, RV* > &model_namePos2Var)
{

  string dictionary_key = makeDictionaryKey(ia_name, section_type);

  // read triangulation method used to produce these cliques.
  is.read(triMethod,"triangulation method string");

  // read number of cliques
  unsigned numCliques;
  is.read(numCliques,"number of cliques");
#if 0
  // remove check for numCliques being > 0 since we now allow for empty sections.
  if (numCliques == 0)
    error("ERROR: reading file '%s' line %d, numCliques must be >= 1\n",
	  is.fileName(),is.lineNo());
#endif

  // create a map for easy access to set of nodes in this section
  // model_namePos2Var includes all model RVs for cumulative sets
  map < RVInfo::rvParent, RV* > namePos2Var;
  for (set<RV*>::iterator i=nodes.begin();
       i != nodes.end(); i++) {
    RV* rv = (*i);
    RVInfo::rvParent par;
    par.first = rv->name();
    par.second = rv->frame();
    namePos2Var[par] = rv;
  }

  vector<unsigned> disposition_vector;
  cliques.reserve(numCliques); // required to avoid dtor of MaxCliques due to resizing, which breaks sArray members
  for (unsigned i=0;i<numCliques;i++) {
    set<RV*> clique;
    
    unsigned cliqueNo;
    is.read(cliqueNo,"clique number value");
    if (cliqueNo != i)
      error("ERROR: reading file %s, line %d, bad cliqueNo (= %d) when reading cliques, out of sequence, should be = %d instead.\n",
	    is.fileName(),is.lineNo(),cliqueNo,i);
    string clique_name;
    is.read(clique_name, "clique name");

    if (clique_name_dictionary[ dictionary_key ].find(clique_name) != clique_name_dictionary[ dictionary_key ].end()) {
      error("ERROR: clique name '%s' already defined in file '%s' line %d\n",
	    clique_name.c_str(), is.fileName(), is.lineNo());
    }
    clique_name_dictionary[ dictionary_key ][ clique_name ] = i;

    unsigned cliqueSize;
    is.read(cliqueSize,"clique size value");

#if 0
    // unsigned can never be less than 0
    if (cliqueSize <= 0)
      error("ERROR: reading file %s line %d, reading clique number %d, but clique size %d must be >= 1\n",
	    is.fileName(),is.lineNo(),i,cliqueSize);
#endif

    for (unsigned j=0;j<cliqueSize;j++) {
      RV *rv = getRVfromIS(is, namePos2Var, "clique RV node specification", i, j);
      clique.insert(rv);
    }
    assert(cliques.size() == i);

    cliques.push_back(MaxClique(clique));
assert(cliques[cliques.size()-1].dispositionSortedAssignedNodes.ptr == NULL);
    disposition_vector.clear();
    // read cliques[i].assigedNodes & sortedAssignedNodes & dispostitions
    unsigned num_sorted_assigned;
    is.read(num_sorted_assigned, "number of sorted assigned RVs");
    for (unsigned j=0; j < num_sorted_assigned; ++j) {
      RV *rv = getRVfromIS(is, namePos2Var, "clique assigned RV specification", i, j);
      unsigned disposition;
      is.read(disposition, "assigned RV disposition");
//printf("read RV: %s(%d) %u", rv->name().c_str(), rv->frame(), disposition);
      cliques[i].assignedNodes.insert(rv);
      if (disposition != MaxClique::AN_CONTINUE) {
//printf("  sass %lu/%lu = %u", cliques[i].sortedAssignedNodes.size(), disposition_vector.size(),disposition);
	cliques[i].sortedAssignedNodes.push_back(rv);
	disposition_vector.push_back(disposition);
      }
//printf("\n");
    }
    cliques[i].dispositionSortedAssignedNodes.resize(disposition_vector.size());
    for (unsigned j=0; j < disposition_vector.size(); ++j) {
      cliques[i].dispositionSortedAssignedNodes[j] = (MaxClique::AssignedNodeDisposition)disposition_vector[j];
    }
    // read cliques[i].assignedProbNodes
    readRVSetFromIS(is, "number of assigned probability RVs", "clique assigned probability RV specification",
		    i, namePos2Var, cliques[i].assignedProbNodes);
    // read cliques[i].cumulativeAssignedNodes
    readRVSetFromIS(is, "number of cumulative assigned RVs", "cumulative assigned RV specification",
		    i, model_namePos2Var, cliques[i].cumulativeAssignedNodes);

    // read cliques[i].cumulativeAssignedProbNodes
    readRVSetFromIS(is, "number of cumulative assigned probability RVs", 
		    "cumulative assigned probability RV specification",
		    i, model_namePos2Var, cliques[i].cumulativeAssignedProbNodes);

    // read cliques[i].unionIncomingCESeps
    readRVSetFromIS(is, "number of incomming separator RVs", 
		    "incomming separator RV specification",
		    i, namePos2Var, cliques[i].unionIncommingCESeps);

    // read cliques[i].unassignedIteratedNodes
    readRVSetFromIS(is, "number of unassigned iterated RVs", 
		    "unassigned iterated RV specification",
		    i, namePos2Var, cliques[i].unassignedIteratedNodes);

    // read cliques[i].cumulativteUnassignedIteratedNodes
    readRVSetFromIS(is, "number of cumulative unassigned iterated RVs", 
		    "cumulative unassigned iterated RV specification",
		    i, model_namePos2Var, cliques[i].cumulativeUnassignedIteratedNodes);
  }

}




void
Section::
readMaxCliques(iDataStreamFile& is,
         string const &ia_name,
         char section_type,
         string const &section_inf_alg,
         map< RVInfo::rvParent, RV* > &model_namePos2Var,
         Json::Value const & inference)
{

  string dictionary_key = makeDictionaryKey(ia_name, section_type);

  // read triangulation method used to produce these cliques.
  if (!inference.isMember("tri_method")) {
    error("ERROR: json file %s section %s %c inference_arch does not have tri_method field\n", 
      is.fileName(), ia_name.c_str(), section_type);
  }
  triMethod = inference["tri_method"].asString();
  // is.read(triMethod,"triangulation method string");


  // read number of cliques
  unsigned numCliques = inference.get("num", 0).asInt();

  // is.read(numCliques,"number of cliques");
#if 0
  // remove check for numCliques being > 0 since we now allow for empty sections.
  if (numCliques == 0)
    error("ERROR: reading file '%s' line %d, numCliques must be >= 1\n",
    is.fileName(),is.lineNo());
#endif

  // create a map for easy access to set of nodes in this section
  // model_namePos2Var includes all model RVs for cumulative sets
  map < RVInfo::rvParent, RV* > namePos2Var;
  for (set<RV*>::iterator i=nodes.begin();
       i != nodes.end(); i++) {
    RV* rv = (*i);
    RVInfo::rvParent par;
    par.first = rv->name();
    par.second = rv->frame();
    namePos2Var[par] = rv;
  }

  if (!inference.isMember("clique_inference")) {
    error("ERROR: json file %s section %s %c inference_arch does not have clique_inference field\n", 
      is.fileName(), ia_name.c_str(), section_type);
  }

  const Json::Value cinfer = inference["clique_inference"];
  if (numCliques != cinfer.size()) {
    error("ERROR: json file %s section %s %c inference_arch num of cliques %d does not match clique_inference list\n", 
      is.fileName(), ia_name.c_str(), section_type, numCliques, cinfer.size());
  }

  vector<unsigned> disposition_vector;
  cliques.reserve(numCliques); // required to avoid dtor of MaxCliques due to resizing, which breaks sArray members
  for (unsigned i=0;i<numCliques;i++) {
    set<RV*> clique;
    
    unsigned cliqueNo = cinfer[i].get("index", numCliques + 1).asInt();
    // is.read(cliqueNo,"clique number value");

    if (cliqueNo != i)
      error("ERROR: reading file %s, line %d, bad cliqueNo (= %d) when reading cliques, out of sequence, should be = %d instead.\n",
      is.fileName(),is.lineNo(),cliqueNo,i);
    string clique_name = cinfer[i].get("name", "nullGarbageuaidsofnqewnmzxcvString").asString();
    // is.read(clique_name, "clique name");

    if (clique_name_dictionary[ dictionary_key ].find(clique_name) != clique_name_dictionary[ dictionary_key ].end()) {
      error("ERROR: clique name '%s' already defined in file '%s' line %d\n",
      clique_name.c_str(), is.fileName(), is.lineNo());
    }
    clique_name_dictionary[ dictionary_key ][ clique_name ] = i;

    if (!cinfer[i].isMember("var_frame")) {
      error("ERROR: json file %s section %s %c inference_arch clique_inference %d does not have var_frame field\n", 
      is.fileName(), ia_name.c_str(), section_type, i);
    }

    unsigned cliqueSize = cinfer[i]["var_frame"].size();
    // is.read(cliqueSize,"clique size value");

#if 0
    // unsigned can never be less than 0
    if (cliqueSize <= 0)
      error("ERROR: reading file %s line %d, reading clique number %d, but clique size %d must be >= 1\n",
      is.fileName(),is.lineNo(),i,cliqueSize);
#endif

    for (unsigned j=0;j<cliqueSize;j++) {
      RV *rv = getRVfromIS(is, namePos2Var, "clique RV node specification", i, j, cinfer[i]["var_frame"][j]);
      clique.insert(rv);
    }
    assert(cliques.size() == i);

    cliques.push_back(MaxClique(clique));
assert(cliques[cliques.size()-1].dispositionSortedAssignedNodes.ptr == NULL);
    disposition_vector.clear();
    // read cliques[i].assigedNodes & sortedAssignedNodes & dispostitions

    if (!cinfer[i].isMember("sorted_assigned")) {
      error("ERROR: json file %s section %s %c inference_arch clique_inference %d does not have sorted_assigned field\n", 
      is.fileName(), ia_name.c_str(), section_type, i);
    }
    unsigned num_sorted_assigned = cinfer[i]["sorted_assigned"].size();
    // is.read(num_sorted_assigned, "number of sorted assigned RVs");

    for (unsigned j=0; j < num_sorted_assigned; ++j) {
      RV *rv = getRVfromIS(is, namePos2Var, "clique assigned RV specification", i, j, cinfer[i]["sorted_assigned"][j]);
      unsigned disposition = cinfer[i]["sorted_assigned"][j][2].asInt();
      // is.read(disposition, "assigned RV disposition");
//printf("read RV: %s(%d) %u", rv->name().c_str(), rv->frame(), disposition);
      cliques[i].assignedNodes.insert(rv);
      if (disposition != MaxClique::AN_CONTINUE) {
//printf("  sass %lu/%lu = %u", cliques[i].sortedAssignedNodes.size(), disposition_vector.size(),disposition);
  cliques[i].sortedAssignedNodes.push_back(rv);
  disposition_vector.push_back(disposition);
      }
//printf("\n");
    }
    cliques[i].dispositionSortedAssignedNodes.resize(disposition_vector.size());
    for (unsigned j=0; j < disposition_vector.size(); ++j) {
      cliques[i].dispositionSortedAssignedNodes[j] = (MaxClique::AssignedNodeDisposition)disposition_vector[j];
    }
    // read cliques[i].assignedProbNodes

    if (!cinfer[i].isMember("assigned_prob")) {
      error("ERROR: json file %s section %s %c inference_arch clique_inference %d does not have assigned_prob field\n", 
      is.fileName(), ia_name.c_str(), section_type, i); 
    }
    readRVSetFromIS(is, "number of assigned probability RVs", "clique assigned probability RV specification",
        i, namePos2Var, cliques[i].assignedProbNodes, cinfer[i]["assigned_prob"]);


    if (!cinfer[i].isMember("cumulative_assigned")) {
      error("ERROR: json file %s section %s %c inference_arch clique_inference %d does not have cumulative_assigned field\n", 
      is.fileName(), ia_name.c_str(), section_type, i); 
    }
    // read cliques[i].cumulativeAssignedNodes
    readRVSetFromIS(is, "number of cumulative assigned RVs", "cumulative assigned RV specification",
        i, model_namePos2Var, cliques[i].cumulativeAssignedNodes, cinfer[i]["cumulative_assigned"]);


    if (!cinfer[i].isMember("cumulative_assigned_prob")) {
      error("ERROR: json file %s section %s %c inference_arch clique_inference %d does not have cumulative_assigned_prob field\n", 
      is.fileName(), ia_name.c_str(), section_type, i); 
    }
    // read cliques[i].cumulativeAssignedProbNodes
    readRVSetFromIS(is, "number of cumulative assigned probability RVs", 
        "cumulative assigned probability RV specification",
        i, model_namePos2Var, cliques[i].cumulativeAssignedProbNodes, cinfer[i]["cumulative_assigned_prob"]);

    if (!cinfer[i].isMember("union_incomming_seps")) {
      error("ERROR: json file %s section %s %c inference_arch clique_inference %d does not have union_incomming_seps field\n", 
      is.fileName(), ia_name.c_str(), section_type, i); 
    }
    // read cliques[i].unionIncomingCESeps
    readRVSetFromIS(is, "number of incomming separator RVs", 
        "incomming separator RV specification",
        i, namePos2Var, cliques[i].unionIncommingCESeps, cinfer[i]["union_incomming_seps"]);


    if (!cinfer[i].isMember("unassigned_iterated")) {
      error("ERROR: json file %s section %s %c inference_arch clique_inference %d does not have unassigned_iterated field\n", 
      is.fileName(), ia_name.c_str(), section_type, i); 
    }
    // read cliques[i].unassignedIteratedNodes
    readRVSetFromIS(is, "number of unassigned iterated RVs", 
        "unassigned iterated RV specification",
        i, namePos2Var, cliques[i].unassignedIteratedNodes, cinfer[i]["unassigned_iterated"]);


    if (!cinfer[i].isMember("cumulative_unassigned_iterated")) {
      error("ERROR: json file %s section %s %c inference_arch clique_inference %d does not have cumulative_unassigned_iterated field\n", 
      is.fileName(), ia_name.c_str(), section_type, i); 
    }    
    // read cliques[i].cumulativteUnassignedIteratedNodes
    readRVSetFromIS(is, "number of cumulative unassigned iterated RVs", 
        "cumulative unassigned iterated RV specification",
        i, model_namePos2Var, cliques[i].cumulativeUnassignedIteratedNodes, cinfer[i]["cumulative_unassigned_iterated"]);
  }

}


/*-
 *-----------------------------------------------------------------------
 * Section::reportScoreStats()
 *   print out stats about the cliques to stdout
 *
 * Preconditions:
 *   The maxclique variable must be instantiated and have valid RVs with cpts/factors associated.
 *
 * Postconditions:
 *   Information about the maxclique variables are written out.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
Section::reportScoreStats()
{
  for (unsigned i=0;i<cliques.size();i++) {
    printf("Clique %d:\n",i);
    cliques[i].reportScoreStats();
  }
}


/*-
 *-----------------------------------------------------------------------
 * Section::triangulateSectionsByCliqueCompletion()
 *   Triangulate the sections by completing the cliques that have been read in.
 *
 * Preconditions:
 *   The corresponding section  must be instantiated.
 *   The maxclique variables cliques must be instantiated!!
 *
 * Postconditions:
 *   the variables pointed to by the cliques will be made complete.
 *
 * Side Effects:
 *   Variables pointed to by cliques will have their
 *   neighbors adjusted.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
Section::
triangulateSectionsByCliqueCompletion()
{
  for (unsigned i=0;i<cliques.size();i++)
    MaxClique::makeComplete(cliques[i].nodes);
}


void
readSectionInterface(iDataStreamFile &is, char const *side, map<string, unsigned> &dictionary, 
		     vector<unsigned> &interface) 
{
  unsigned i_size;
  string msg(string(side) + string(" interface size"));
  is.read(i_size, msg.c_str());
  for (unsigned j=0; j < i_size; ++j) {
    unsigned root_clique_idx;
    is.read(root_clique_idx, string(string(side)+string(" interface source clique index")).c_str());
    if (root_clique_idx != j) {
      error("ERROR: expected %s interface source clique number %u, but got %u at line %d of '%s'\n",
	    side, j, root_clique_idx, is.lineNo(), is.fileName());
    }
    string root_clique_name;
    is.read(root_clique_name, string(string(side) + string(" interface source clique name")).c_str());
    if (dictionary.find(root_clique_name) == dictionary.end()) {
      error("ERROR: unknown clique name '%s' at %s interface %u at line %d in file '%s'\n",
	    root_clique_name.c_str(), side, j, is.lineNo(), is.fileName());
    }
    interface.push_back( dictionary[root_clique_name] );
  }
}


void
readSectionInterface(iDataStreamFile &is, char const *side, map<string, unsigned> &dictionary, 
         vector<unsigned> &interface, Json::Value const jsonInter) 
{
  unsigned i_size;
  string msg(string(side) + string(" interface size"));
  i_size = jsonInter.size();
  // is.read(i_size, msg.c_str());
  for (unsigned j=0; j < i_size; ++j) {
    unsigned root_clique_idx;
    root_clique_idx = j;
    // is.read(root_clique_idx, string(string(side)+string(" interface source clique index")).c_str());
    // if (root_clique_idx != j) {
      // error("ERROR: expected %s interface source clique number %u, but got %u at line %d of '%s'\n",
      // side, j, root_clique_idx, is.lineNo(), is.fileName());
    // }
    string root_clique_name = jsonInter[j].asString();
    is.read(root_clique_name, string(string(side) + string(" interface source clique name")).c_str());
    if (dictionary.find(root_clique_name) == dictionary.end()) {
      error("ERROR: unknown clique name '%s' at %s interface %u at line %d in file '%s'\n",
      root_clique_name.c_str(), side, j, is.lineNo(), is.fileName());
    }
    interface.push_back( dictionary[root_clique_name] );
  }
}



void 
Section::readInferenceArchitectureDefinition(iDataStreamFile &is,
					     string const &ia_name,
					     char section_type,
					     string const &section_inf_alg)
{
  string dictionary_key = makeDictionaryKey(ia_name, section_type);
  map<string, unsigned> &dictionary = clique_name_dictionary[ dictionary_key ];

  vector< pair<unsigned, unsigned> > msg_order;

  unsigned num_msgs;
  is.read(num_msgs, "number of messages");
  for (unsigned i=0; i < num_msgs; ++i) {
    unsigned index;
    is.read(index, "message number");
    if (index != i) {
      error("ERROR: reading file '$s' line %d, bad message number (= %u) out of sequence, should be %u instead.\n",
	    is.fileName(), is.lineNo(), index, i);
    }
    string source_clique_name, dest_clique_name;
    is.read(source_clique_name, "message source clique name");
    if (dictionary.find(source_clique_name) == dictionary.end()) {
      error("ERROR: unknown source clique name '%s' in message %u at line %d in file '%s'\n", 
	    source_clique_name.c_str(), i, is.lineNo(), is.fileName());
    }
    is.read(dest_clique_name, "message destination clique name");
    if (dictionary.find(dest_clique_name) == dictionary.end()) {
      error("ERROR: unknown source clique name '%s' in message %u at line %d in file '%s'\n", 
	    dest_clique_name.c_str(), i, is.lineNo(), is.fileName());
    }
    pair<unsigned, unsigned> msg(dictionary[source_clique_name], dictionary[dest_clique_name]);
    msg_order.push_back(msg);
  }
  ia_message_order[dictionary_key] = msg_order;

  // forward pass backwards messages to make multiple right interface cliques consistent
  // before projecting down for next section

  msg_order.clear();
  is.read(num_msgs, "number of reverse messages");
  for (unsigned i=0; i < num_msgs; ++i) {
    unsigned index;
    is.read(index, "message number");
    if (index != i) {
      error("ERROR: reading file '$s' line %d, bad message number (= %u) out of sequence, should be %u instead.\n",
	    is.fileName(), is.lineNo(), index, i);
    }
    string source_clique_name, dest_clique_name;
    is.read(source_clique_name, "message source clique name");
    if (dictionary.find(source_clique_name) == dictionary.end()) {
      error("ERROR: unknown source clique name '%s' in message %u at line %d in file '%s'\n", 
	    source_clique_name.c_str(), i, is.lineNo(), is.fileName());
    }
    is.read(dest_clique_name, "message destination clique name");
    if (dictionary.find(dest_clique_name) == dictionary.end()) {
      error("ERROR: unknown source clique name '%s' in message %u at line %d in file '%s'\n", 
	    dest_clique_name.c_str(), i, is.lineNo(), is.fileName());
    }
    pair<unsigned, unsigned> msg(dictionary[source_clique_name], dictionary[dest_clique_name]);
    msg_order.push_back(msg);
  }
  ia_rev_msg_order[dictionary_key] = msg_order;

  // read section's left interface cliques
  readSectionInterface(is, "left", dictionary, section_li);
  readSectionInterface(is, "right", dictionary, section_ri);
}




void 
Section::readInferenceArchitectureDefinition(iDataStreamFile &is,
               string const &ia_name,
               char section_type,
               string const &section_inf_alg,
               Json::Value const & inference)
{

  if (!inference.isMember("message_order")) {
    error("ERROR: json file %s section %s %c inference_arch does not have message_order field\n", 
      is.fileName(), ia_name.c_str(), section_type);
  }
  const Json::Value mOrder = inference["message_order"];

  string dictionary_key = makeDictionaryKey(ia_name, section_type);
  map<string, unsigned> &dictionary = clique_name_dictionary[ dictionary_key ];

  vector< pair<unsigned, unsigned> > msg_order;

  unsigned num_msgs = mOrder.get("num_msg_forward", 0).asInt();
  if (!mOrder.isMember("forward_messages")) {
    error("ERROR: json file %s section %s %c message_order does not have forward_messages field\n", 
      is.fileName(), ia_name.c_str(), section_type); 
  }
  if (num_msgs != mOrder["forward_messages"].size()){
    error("ERROR: json file %s section %s %c message_order number of forward messages %d does not match num in the foward_messages list %d\n", 
      is.fileName(), ia_name.c_str(), section_type, num_msgs, mOrder["forward_messages"].size());
  }

  // is.read(num_msgs, "number of messages");
  for (unsigned i=0; i < num_msgs; ++i) {
    unsigned index = i;
    // is.read(index, "message number");
    // if (index != i) {
    //   error("ERROR: reading file '$s' line %d, bad message number (= %u) out of sequence, should be %u instead.\n",
    //   is.fileName(), is.lineNo(), index, i);
    // }
    string source_clique_name, dest_clique_name;
    source_clique_name = mOrder["forward_messages"][i][0].asString();
    // is.read(source_clique_name, "message source clique name");
    if (dictionary.find(source_clique_name) == dictionary.end()) {
      error("ERROR: unknown source clique name '%s' in message %u at line %d in file '%s'\n", 
      source_clique_name.c_str(), i, is.lineNo(), is.fileName());
    }
    dest_clique_name = mOrder["forward_messages"][i][1].asString();
    // is.read(dest_clique_name, "message destination clique name");
    if (dictionary.find(dest_clique_name) == dictionary.end()) {
      error("ERROR: unknown source clique name '%s' in message %u at line %d in file '%s'\n", 
      dest_clique_name.c_str(), i, is.lineNo(), is.fileName());
    }
    pair<unsigned, unsigned> msg(dictionary[source_clique_name], dictionary[dest_clique_name]);
    msg_order.push_back(msg);
  }
  ia_message_order[dictionary_key] = msg_order;

  // forward pass backwards messages to make multiple right interface cliques consistent
  // before projecting down for next section

  msg_order.clear();

  num_msgs = mOrder.get("num_msg_backward", 0).asInt();
  if (!mOrder.isMember("backward_messages")) {
    error("ERROR: json file %s section %s %c message_order does not have backward_messages field\n", 
      is.fileName(), ia_name.c_str(), section_type); 
  }
  if (num_msgs != mOrder["backward_messages"].size()){
    error("ERROR: json file %s section %s %c message_order number of backward messages %d does not match num in the backward_messages list %d\n", 
      is.fileName(), ia_name.c_str(), section_type, num_msgs, mOrder["backward_messages"].size());
  }
  // is.read(num_msgs, "number of reverse messages");
  for (unsigned i=0; i < num_msgs; ++i) {
    unsigned index = i;
    // is.read(index, "message number");
    // if (index != i) {
    //   error("ERROR: reading file '$s' line %d, bad message number (= %u) out of sequence, should be %u instead.\n",
    //   is.fileName(), is.lineNo(), index, i);
    // }
    string source_clique_name, dest_clique_name;
    source_clique_name = mOrder["backward_messages"][i][0].asString();
    // is.read(source_clique_name, "message source clique name");
    if (dictionary.find(source_clique_name) == dictionary.end()) {
      error("ERROR: unknown source clique name '%s' in message %u at line %d in file '%s'\n", 
      source_clique_name.c_str(), i, is.lineNo(), is.fileName());
    }
    dest_clique_name = mOrder["backward_messages"][i][1].asString();
    // is.read(dest_clique_name, "message destination clique name");
    if (dictionary.find(dest_clique_name) == dictionary.end()) {
      error("ERROR: unknown source clique name '%s' in message %u at line %d in file '%s'\n", 
      dest_clique_name.c_str(), i, is.lineNo(), is.fileName());
    }
    pair<unsigned, unsigned> msg(dictionary[source_clique_name], dictionary[dest_clique_name]);
    msg_order.push_back(msg);
  }
  ia_rev_msg_order[dictionary_key] = msg_order;

  // read section's left interface cliques
  if (!mOrder.isMember("li")) {
    error("ERROR: json file %s section %s %c inference_arch message_order does not have li field\n", 
      is.fileName(), ia_name.c_str(), section_type);
  }
  readSectionInterface(is, "left", dictionary, section_li, mOrder["li"]);
  if (!mOrder.isMember("li")) {
    error("ERROR: json file %s section %s %c inference_arch message_order does not have ri field\n", 
      is.fileName(), ia_name.c_str(), section_type);
  }
  readSectionInterface(is, "right", dictionary, section_ri, mOrder["ri"]);
}


/*-
 *-----------------------------------------------------------------------
 * Section::setCliquesFromAnotherSection()
 *   Set the cliques from anohter section. The other section
 *   must be from the same structure file, and if the current section is  a P (resp. C, E)
 *   than the other section must also be a P (resp. C, and E). Also, the sections
 *   are assumed to come from the same boundary for the .str file.
 *   The routine is used to merge together gm_templates for different triangulations
 *   of the same boundary but say a paralle triangulation of P, C, and E.
 *  
 *   Note that it is assumed that the different sections refer to different instantiations of
 *   the same set of random variables (so we can't use rv1 == rv2, but instead must use
 *   name and frame equality).
 *
 * Preconditions:
 *   The corresponding section must be instantiated with nodes.
 *   
 *
 * Postconditions:
 *   The clique variables now refer to the cliques in the other section.
 *
 * Side Effects:
 *   Current cliques will be destroyed and set to new versions.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
Section::setCliquesFromAnotherSection(Section& from_part)
{
  cliques.clear();

  // create a rv set of just the IDs for the dest
  map < RVInfo::rvParent, unsigned > ppf;
  vector <RV*> newRvs;
  set<RV*>::iterator it;

  newRvs.reserve(nodes.size());
  unsigned i;
  for (i=0,it = nodes.begin();
       it != nodes.end();
       i++,it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame();
    ppf[rvp] = i;
    newRvs.push_back(rv);
  }

  // make sure nodes refer to same section.
  for (it = from_part.nodes.begin();
       it != from_part.nodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame();

    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find rv %s(%d) in RV set of dest section\n",
	    rvp.first.c_str(),rvp.second);
      assert ( 0 );
    }

  }  
  cliques.reserve(from_part.cliques.size());
  
  for (unsigned i=0;i<from_part.cliques.size();i++) {
    cliques.push_back(MaxClique(from_part.cliques[i],
				newRvs,ppf,0));
  }
  // copy tri-method string as well.
  triMethod = from_part.triMethod;
}


