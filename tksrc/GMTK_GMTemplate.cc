/*-
 * GMTK_GMTemplate.cc
 *    Basic GM Template structure for a given graph file.
 *    This includes code that is common to both triangulation and inference,
 *    so does not contain the more elaborate triangulation methods so that
 *    they are not appart of inference code.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */



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
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_JunctionTree.h"

VCID("$Header$");

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


const string GMTemplate::P_partition_name("P_PARTITION");
const string GMTemplate::C_partition_name("C_PARTITION");
const string GMTemplate::E_partition_name("E_PARTITION");
const string GMTemplate::PC_interface_name("PC_PARTITION");
const string GMTemplate::CE_interface_name("CE_PARTITION");

const string GMTemplate::fileExtension(".trifile");


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors/Destructors 
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

Partition::Partition(Partition& from_part,
		     vector <RandomVariable*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta)
{

  triMethod = from_part.triMethod;

  set<RandomVariable*>::iterator it;

  // clone over nodes RVs.  
  // TODO: make this next code a routine
  //  nodesClone() since it is used in several places.
  for (it = from_part.nodes.begin();
       it != from_part.nodes.end();
       it++) {
    RandomVariable* rv = (*it);
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

    RandomVariable* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }
  cliques.reserve(from_part.cliques.size());
  // 
  // NOTE: It is Crucial for the cliques in the cloned partition to be
  // inserted in the *SAME ORDER* as in the partition being cloned.
  for (unsigned i=0;i<from_part.cliques.size();i++) {
    cliques.push_back(MaxClique(from_part.cliques[i],
				newRvs,ppf,frameDelta));
  }
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Partition support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * Partition::writeMaxCliques()
 *   Write out the max cliques of the given partitions.
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
Partition::
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
      totalWeight = totalWeight + log10(1+pow(10,curWeight-totalWeight));
    os.writeComment("%d : %d  %f\n",
		    i,
		    cliques[i].nodes.size(),curWeight);
    for (set<RandomVariable*>::iterator j=cliques[i].nodes.begin();
	 j != cliques[i].nodes.end(); j++) {
      RandomVariable* rv = (*j);
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
    for (set<RandomVariable*>::iterator j=cliques[i].nodes.begin();
	 j != cliques[i].nodes.end(); j++) {
      RandomVariable* rv = (*j);
      os.write(rv->name().c_str());
      os.write(rv->frame());
    }
    os.nl();
  }
}


/*-
 *-----------------------------------------------------------------------
 * Partition::readMaxCliques()
 *   Write out the max cliques of the given partitions.
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
Partition::
readMaxCliques(iDataStreamFile& is)
{
  // read triangulation method used to produce these cliques.
  is.read(triMethod,"triangulation method string");

  // read number of cliques
  unsigned numCliques;
  is.read(numCliques,"numCliques value");
  if (numCliques == 0)
    error("ERROR: reading file %s, numCliques must be >= 1\n",is.fileName());

  // create a map for easy access to set of nodes
  map < RVInfo::rvParent, RandomVariable* > namePos2Var;
  for (set<RandomVariable*>::iterator i=nodes.begin();
       i != nodes.end(); i++) {
    RandomVariable* rv = (*i);
    RVInfo::rvParent par;
    par.first = rv->name();
    par.second = rv->frame();
    namePos2Var[par] = rv;
  }

  for (unsigned i=0;i<numCliques;i++) {
    set<RandomVariable*> clique;
    
    unsigned cliqueNo;
    is.read(cliqueNo,"cliqueNo value");
    if (cliqueNo != i)
      error("ERROR: reading file %s, bad cliqueNo (= %d) when reading cliques, out of sequence, should be = %d, file %s\n",is.fileName(),cliqueNo,i);
    
    unsigned cliqueSize;
    is.read(cliqueSize,"cliqueSize value");
    if (cliqueSize <= 1)
      error("ERROR: reading file %s, cliqueSize %d must be >= 2\n",
	    is.fileName(),cliqueSize);
    
    for (unsigned j=0;j<cliqueSize;j++) {

      RVInfo::rvParent par;
      is.read(par.first,"parent name");
      is.read(par.second,"parent position");

      map < RVInfo::rvParent, RandomVariable* >::iterator loc;
      loc = namePos2Var.find(par);
      if (loc == namePos2Var.end())
	error("ERROR: reading file %s, clique specification %d has %d'th variable %s(%d) that does not exist in partition.\n",
	      is.fileName(),i,j,par.first.c_str(),par.second);
      RandomVariable* rv = (*loc).second;
      clique.insert(rv);
    }

    cliques.push_back(MaxClique(clique));
  }

}


/*-
 *-----------------------------------------------------------------------
 * Partition::triangulatePartitionsByCliqueCompletion()
 *   Triangulate the partitions by completing the cliques that have been read in.
 *
 * Preconditions:
 *   The corresponding partition  must be instantiated.
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
Partition::
triangulatePartitionsByCliqueCompletion()
{
  for (unsigned i=0;i<cliques.size();i++)
    MaxClique::makeComplete(cliques[i].nodes);
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Private Support Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::setUpClonedPartitionGraph()
 *   Given a P,C, and E that has been cloned without parents, this 
 *   will set up the partitions as follows.
 *   only the members that are in the corresponding 'in' set (i.e., if
 *   any parents, children, neighbors, in 'in' pointed to variables
 *   outside of 'in', then the corresponding variables in 'out' do not
 *   contain those parents,children,neighbors.
 *
 *   There are two versions of this routine, one which returns
 *   the in to out variable map in case that might be useful.
 *
 * Preconditions:
 *   'in' is a set of random variables to be cloned.
 *
 * Postconditions:
 *   'out' is a clone of 'in' but without parents,children,neighbors.
 *   
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     returns the in_to_out mapping
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
setUpClonedPartitionGraph(const set<RandomVariable*>& P,
			  const set<RandomVariable*>& C,
			  const set<RandomVariable*>& E,
			  // cloned variables
			  set<RandomVariable*>& Pc,
			  set<RandomVariable*>& Cc,
			  set<RandomVariable*>& Ec,
			  // next 3 should be const but ther eis no "op[] const"
			  map < RandomVariable*, RandomVariable* >& P_in_to_out,
			  map < RandomVariable*, RandomVariable* >& C_in_to_out,
			  map < RandomVariable*, RandomVariable* >& E_in_to_out)
{

  // just do neighbors for now, don't bother with parents, children,
  // and so on.
  // Set the neighbors of out to be the correctly associated
  // variables, but do not include neighbors that are not in the
  // current set (i.e., dissociate with any other possible portion of
  // the network)

  cloneWithoutParents(P,Pc,P_in_to_out);
  cloneWithoutParents(C,Cc,C_in_to_out);
  cloneWithoutParents(E,Ec,E_in_to_out);

  setPartitionParentsChildrenNeighbors(P,Pc,P_in_to_out,C_in_to_out,E_in_to_out);
  setPartitionParentsChildrenNeighbors(C,Cc,C_in_to_out,P_in_to_out,E_in_to_out);
  setPartitionParentsChildrenNeighbors(E,Ec,E_in_to_out,P_in_to_out,C_in_to_out);

}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::setPartitionChildrenNeighbors()
 *   This routine takes as input an original partition S from an unrolled
 *   graph, and its unfinished cloned version Sc (unfinished in that
 *   it has no parents, neighbors (and children) set up yet). It also
 *   takes a mapping from S to Sc, and the corresponding mappings from
 *   the other partitions whatever they are, called O1 and O1. 
 *   
 *   It then sets the neighbors structures for all of Sc from S. The neighbors
 *   are variables that are forced to be part of the partition S itself.
 *
 *   It then sets the parents variables for each variable in Sc. The parents
 *   (and also the children) might NOT be fully contained in Sc, so it needs
 *   to use the mappings for O1 and O2 to get the location of the correspondly
 *   cloned variables for those partitions.
 *
 * Preconditions:
 *   Sc is unfishined, i.e., variables have been cloned without parents
 *
 * Postconditions:
 *   Sc is finished, as described above.
 *   
 *
 * Side Effects:
 *   Changes parents, neighbors, and children of all variables in Sc
 *
 * Results:
 *     results returned via output variable Sc and its parents, neighbors, children
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
setPartitionParentsChildrenNeighbors(const set<RandomVariable*>& S,
				     set<RandomVariable*>& Sc,
				     // next 3 should be const but there is no "op[] const" in STL
				     map < RandomVariable*, RandomVariable* >& S_in_to_out,
				     map < RandomVariable*, RandomVariable* >& O1_in_to_out,
				     map < RandomVariable*, RandomVariable* >& O2_in_to_out)
{

  for (set<RandomVariable*>::iterator i=S.begin();i != S.end(); i++) {

    RandomVariable*rv = (*i);

    // first set up new neighbors for S_in_to_out[rv]
    // Note: Neighbors are defined to point ONLY TO OTHER VARIABLES
    // IN THE SET S_in_to_out (i.e., the members of the partition). Any
    // other parents or children are not contained in that set are not included.
    set<RandomVariable*> tmp;
    for (set<RandomVariable*>::iterator j = rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {    
      if (S_in_to_out.find((*j)) != S_in_to_out.end()) {
	// then it is included in this set.
	tmp.insert(S_in_to_out[(*j)]);
      }
    }
    S_in_to_out[rv]->neighbors = tmp;
    // assertion to make sure that no node has itself as neighbor.
    assert( S_in_to_out[rv]->neighbors.find(S_in_to_out[rv])
	    == S_in_to_out[rv]->neighbors.end() );

    // next, set new sparents for in_to_out[rv].
    // Note that the parents might be outside of the set S.
    vector<RandomVariable *> sParents;
    for (unsigned l=0;l<rv->switchingParents.size();l++) {
      // grab a copy for readability
      RandomVariable* const par = rv->switchingParents[l];
      if (S_in_to_out.find(par) != S_in_to_out.end())
	sParents.push_back(S_in_to_out[par]);
      else if (O1_in_to_out.find(par) != O1_in_to_out.end())
	sParents.push_back(O1_in_to_out[par]);
      else if (O2_in_to_out.find(par) != O2_in_to_out.end())
	sParents.push_back(O2_in_to_out[par]);
      else
	// this shouldn't happen since the parent should live
	// in one of the three partitions.
	assert ( 0 );
    }

    // next, set conditional parents
    vector< vector < RandomVariable* > > cParentsList;
    cParentsList.resize(rv->conditionalParentsList.size());
    for (unsigned l=0;l<rv->conditionalParentsList.size();l++) {
      for (unsigned m=0;m<rv->conditionalParentsList[l].size();m++) {
	// grab a copy for readability
	RandomVariable* const par = rv->conditionalParentsList[l][m];
	if (S_in_to_out.find(par) != S_in_to_out.end())
	  cParentsList[l].push_back(S_in_to_out[par]);
	else if (O1_in_to_out.find(par) != O1_in_to_out.end())
	  cParentsList[l].push_back(O1_in_to_out[par]);
	else if (O2_in_to_out.find(par) != O2_in_to_out.end())
	  cParentsList[l].push_back(O2_in_to_out[par]);
	else
	  // this shouldn't happen since the parent should live
	  // in one of the three partitions.
	  assert ( 0 );
      }
    }
    S_in_to_out[rv]->setParents(sParents,cParentsList);
  }
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::cloneWithtoutParents()
 *   Clone a set of random variables from 'in' to 'out'. The cloned
 *   variables have only empty parents, children, and neighbors structures.

 * Preconditions:
 *   'in' is a set of random variables to be cloned.
 *
 * Postconditions:
 *   'out' is a clone of 'in' but without parents,children,neighbors.
 *   
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     returns the in_to_out mapping
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
cloneWithoutParents(const set<RandomVariable*>& in, 
		    set<RandomVariable*>& out,
		    map < RandomVariable*, RandomVariable* >& in_to_out)
{
  in_to_out.clear();
  out.clear();

  for (set<RandomVariable*>::iterator i=in.begin();
       i != in.end(); i++) {

    // sanity check, to ensure a node is not its own neighbor
    assert ( (*i)->neighbors.find((*i)) == (*i)->neighbors.end() );

    RandomVariable*rv = (*i)->cloneWithoutParents();
    out.insert(rv);
    in_to_out[(*i)] = rv;

  }
}



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Reading And Writing Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::writePartitions()
 *
 *  Writes the given argument partitions into a file for later
 *  retreival. This routine writes out the information both in a more
 *  human readable format (as comments preceeded by a coment character
 *  which includes other useful information) and in machine readable
 *  form to be read in again (lines that do not begin with comment
 *  characters).  The information written includes;
 *
 *    P partition
 *    C partition
 *    E partition
 *    PC interface
 *    CE interface
 *
 * Preconditions:
 *
 *   Each partition must corresond to a valid and separte GM. Each
 *   variable in each GM must have a valid parent and neighbor members
 *   and the parents/neighbors must only point to other members of a
 *   given partition.
 * 
 * Postconditions:
 *   Each of the partitions have been printed
 *
 * Side Effects:
 *   none, other than changing the file pointer of os
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
writePartitions(oDataStreamFile& os, string& str)
{

  string buffer;
  char buff[2048];

  // Write out current time/date.
  os.nl();
  os.writeComment("---\n");
  {
    time_t tloc;
    struct tm*tms;
    time(&tloc);
    tms = localtime(&tloc);
    strftime(buff,2048,"%A %B %d %Y, %H:%M:%S %Z",tms);
  }
  os.writeComment("File Created: %s\n",buff);
  os.writeComment("Options-: %s\n",str.c_str());
  os.writeComment("---\n");
  os.nl();


  // number of chunks in which to find interface boundary
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- M, number of chunks in which to find interface boundary\n");
  os.write(M);
  os.nl();

  // chunk skip
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- S, chunk skip\n");
  os.write(S);
  os.nl();

  // interface method
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- interface method\n");
  os.write((leftInterface?"LEFT":"RIGHT"));
  os.nl();


  // write out information about method used to create current
  // boundary
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- boundary method\n");
  printf("about to write boundary method = %s\n",boundaryMethod.c_str());
  if (boundaryMethod.size() > 0) {
    // make sure string has no white space.
    for (unsigned i=0;i<boundaryMethod.size();i++) {
      if (isspace(boundaryMethod[i]))
	boundaryMethod[i] = '_';
    }
    os.write(boundaryMethod.c_str());
  } else {
    os.write("UNKNOWN_BOUNDARY_METHOD");
  }
  os.nl();  


  // next write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- P partition information: variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=P.nodes.begin();
       i != P.nodes.end(); i++) {
    RandomVariable* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      sprintf(buff," %s(%d),",
	      (*j)->name().c_str(),(*j)->frame());
      buffer += buff;
    }
    os.writeComment("%s\n",buffer.c_str());
  }
  // Then write it out in machine readable form not as a comment
  os.writeComment("--- P partition definition\n");
  os.write(P_partition_name);
  os.write(P.nodes.size());
  for (set<RandomVariable*>::iterator i = P.nodes.begin();
       i != P.nodes.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();

  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- C partition information: variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=C.nodes.begin();
       i != C.nodes.end(); i++) {
    RandomVariable* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      sprintf(buff," %s(%d),",
	      (*j)->name().c_str(),(*j)->frame());
      buffer += buff;
    }
    os.writeComment("%s\n",buffer.c_str());
  }
  // Then write it out in machine readable form not as a comment
  os.writeComment("--- C partition definition\n");
  os.write(C_partition_name);
  os.write(C.nodes.size());
  for (set<RandomVariable*>::iterator i = C.nodes.begin();
       i != C.nodes.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- E partition information: variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=E.nodes.begin();
       i != E.nodes.end(); i++) {
    RandomVariable* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      sprintf(buff," %s(%d),",
	      (*j)->name().c_str(),(*j)->frame());
      buffer += buff;
    }
    os.writeComment("%s\n",buffer.c_str());
  }
  // Then write it out in machine readable form not as a comment
  os.writeComment("--- E partition definition\n");
  os.write(E_partition_name);
  os.write(E.nodes.size());
  for (set<RandomVariable*>::iterator i = E.nodes.begin();
       i != E.nodes.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- PC information : variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=PCInterface_in_C.begin();
       i != PCInterface_in_C.end(); i++) {
    RandomVariable* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      sprintf(buff," %s(%d),",
	      (*j)->name().c_str(),(*j)->frame());
      buffer += buff;
    }
    os.writeComment("%s\n",buffer.c_str());
  }
  // Then write it out in machine readable form not as a comment
  os.writeComment("--- PC interface definition\n");
  os.write(PC_interface_name);
  os.write(PCInterface_in_C.size());
  for (set<RandomVariable*>::iterator i = PCInterface_in_C.begin();
       i != PCInterface_in_C.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- CE information : variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=CEInterface_in_C.begin();
       i != CEInterface_in_C.end(); i++) {
    RandomVariable* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      sprintf(buff," %s(%d),",
	      (*j)->name().c_str(),(*j)->frame());
      buffer += buff;
    }
    os.writeComment("%s\n",buffer.c_str());
  }
  // Then write it out in machine readable form not as a comment
  os.writeComment("--- CE interface definition\n");
  os.write(CE_interface_name);
  os.write(CEInterface_in_C.size());
  for (set<RandomVariable*>::iterator i = CEInterface_in_C.begin();
       i != CEInterface_in_C.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();

}

/*-
 *-----------------------------------------------------------------------
 * GMTemplate::readPartitions()
 *  Create the three partitions (P,C,E) of the template using the
 *  given information stored in the input file. See
 *  findPartitions() routine above for argument definitions.
 *
 * Preconditions:
 *   Object must be instantiated and have the use of the information
 *   in a valid FileParser object (which stores the parsed structure
 *   file) Arguments must indicate valid heuristics to use.
 *
 * Postconditions:
 *   Arguments Pc, Cc, and Ec now contain partitions in a separate
 *   graph from the FileParser (i.e., file parser information is not
 *   disturbed)
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   Pc, Cc, and Ec as arguments.
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
readPartitions(iDataStreamFile& is)
{
  unsigned loc_S,loc_M;
  string loc_I;

  is.read(loc_M,"M value");

  if (loc_M == 0)
    error("ERROR: reading file %s, M (number of chunks in which to find interface boundary) must be >= 1\n",
	  is.fileName());

  if (M == GMTEMPLATE_UNINITIALIZED_MS) {
    // This is a const cast hack to get around the fact that member M
    // is declared constant, in cases where we create an object with
    // an uninitialized M.  Perhaps the right thing to do is to
    // undeclare M constant.
    unsigned* Mp = &(unsigned)M;
    *Mp = loc_M;
  } else {
    if (loc_M != M)
      error("ERROR: reading file %s, M (=%d) given in tri-file does not equal %d\n",
	    is.fileName(),loc_M,M);
  }

  is.read(loc_S,"S value");
  if (loc_S == 0)
    error("ERROR: reading file %s, S (chunk skip) must be >= 1\n",
	  is.fileName());

  if (S == GMTEMPLATE_UNINITIALIZED_MS) {
    unsigned *Sp = &(unsigned)S;
    *Sp = loc_S;
  } else {
    if (loc_S != S)
      error("ERROR: reading file %s, S in file (%d) does not equal %d\n",
	    is.fileName(),loc_S,S);
  }

  // interface method
  is.read(loc_I,"interface method value");
  if (loc_I == "LEFT") {
    leftInterface = true;
  } else if  (loc_I == "RIGHT") {
    leftInterface = false;
  } else {
      error("ERROR: reading file %s, interface in file must be 'LEFT' or 'RIGHT' but got string '%s'\n",
	    is.fileName(),loc_I.c_str());
  }

  // read in information about method used to create current boundary
  is.read(boundaryMethod,"boundary method string");

  vector <RandomVariable*> unrolled_rvs;
  map < RVInfo::rvParent, unsigned > positions;
  fp.unroll(M+S-1,unrolled_rvs,positions);

  // need to moralize.
  for (unsigned i=0;i<unrolled_rvs.size();i++) {
    unrolled_rvs[i]->createNeighborsFromParentsChildren();
  }
  for (unsigned i=0;i<unrolled_rvs.size();i++) {
    unrolled_rvs[i]->moralize();    
  }

  // create temporary local variables.
  set<RandomVariable*> loc_P;
  set<RandomVariable*> loc_C;
  set<RandomVariable*> loc_E;

  unsigned setSize;
  string str_tmp;

  is.read(str_tmp,"name");
  if (str_tmp != P_partition_name)
    error("ERROR: P partition information in file %s is invalid for given graph structure\n",is.fileName());
  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: P partition information in file %s is invalid for given graph structure\n",
	    is.fileName());
    loc_P.insert(unrolled_rvs[(*loc).second]);
  }

  is.read(str_tmp,"name");
  if (str_tmp != C_partition_name)
    error("ERROR: C partition information in file %s is invalid for given graph structure\n",is.fileName());
  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: C partition information in file %s is invalid for given graph structure\n",
	    is.fileName());
    loc_C.insert(unrolled_rvs[(*loc).second]);
  }

  is.read(str_tmp,"name");
  if (str_tmp != E_partition_name)
    error("ERROR: E partition information in file %s is invalid for given graph structure\n",is.fileName());
  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: E partition information in file %s is invalid for given graph structure\n",
	    is.fileName());
    loc_E.insert(unrolled_rvs[(*loc).second]);
  }

  //////////////////////////////////////////////
  // next, read in the interface definitions. //
  //////////////////////////////////////////////
  set<RandomVariable*> loc_PCInterface;
  set<RandomVariable*> loc_CEInterface;

  // get PC interface
  is.read(str_tmp,"name");
  if (str_tmp != PC_interface_name)
    error("ERROR: PC interface information in file %s is invalid for given graph structure\n",is.fileName());
  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: PC interface information in file %s is invalid for given graph structure\n",
	    is.fileName());
    loc_PCInterface.insert(unrolled_rvs[(*loc).second]);
  }

  // get CE interface
  is.read(str_tmp,"name");
  if (str_tmp != CE_interface_name)
    error("ERROR: CE interface information in file %s is invalid for given graph structure\n",is.fileName());
  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: CE interface information in file %s is invalid for given graph structure\n",
	    is.fileName());
    loc_CEInterface.insert(unrolled_rvs[(*loc).second]);
  }

  /////////////////////////////////////////////////////////////////////
  // finally create a new variable set for each, make the interfaces
  // complete, and finish up.

  createPartitions(loc_P,loc_C,loc_E,loc_PCInterface,loc_CEInterface);

}

/*-
 *-----------------------------------------------------------------------
 * GMTemplate::createPartitions()
 *   Create the internal partitions for this GMTemplate. 
 *
 * Preconditions:
 *   The arguments P, C, and E must point to a set of random variables
 *   and the interface arguments must point to the interfaces in that variable set
 *   all from the same graph that have been unrolled M+S-1 times.
 *
 * Postconditions:
 *   the internal partition and interface variables are created. 
 *
 * Side Effects:
 *   Deletes any exiting partitions in the object.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
createPartitions(const set<RandomVariable*>& arg_P,
		 const set<RandomVariable*>& arg_C,
		 const set<RandomVariable*>& arg_E,
		 const set<RandomVariable*>& arg_PCInterface,
		 const set<RandomVariable*>& arg_CEInterface)
{

  // delete old stuff
  clear();

  /////////////////////////////////////////////////////////////////////
  // create a new variable set for each, make the interfaces
  // complete, and finish up.

  map < RandomVariable*, RandomVariable* > P_in_to_out;
  map < RandomVariable*, RandomVariable* > C_in_to_out;
  map < RandomVariable*, RandomVariable* > E_in_to_out;

  setUpClonedPartitionGraph(arg_P,arg_C,arg_E,P.nodes,C.nodes,E.nodes,P_in_to_out,C_in_to_out,E_in_to_out);

  // complete the PC interface in P
  PCInterface_in_P.clear();
  for (set<RandomVariable*>::iterator i=arg_PCInterface.begin();
       i != arg_PCInterface.end();i++) {
    PCInterface_in_P.insert(P_in_to_out[(*i)]);
  }
  MaxClique::makeComplete(PCInterface_in_P);

  // complete the PC interface in C
  PCInterface_in_C.clear();
  for (set<RandomVariable*>::iterator i=arg_PCInterface.begin();
       i != arg_PCInterface.end();i++) {
    PCInterface_in_C.insert(C_in_to_out[(*i)]);
  }
  MaxClique::makeComplete(PCInterface_in_C);

  // complete the CE interface in C
  CEInterface_in_C.clear();
  for (set<RandomVariable*>::iterator i=arg_CEInterface.begin();
       i != arg_CEInterface.end();i++) {
    CEInterface_in_C.insert(C_in_to_out[(*i)]);
  }
  MaxClique::makeComplete(CEInterface_in_C);  

  // complete the CE interface in E
  CEInterface_in_E.clear();
  for (set<RandomVariable*>::iterator i=arg_CEInterface.begin();
       i != arg_CEInterface.end();i++) {
    CEInterface_in_E.insert(E_in_to_out[(*i)]);
  }
  MaxClique::makeComplete(CEInterface_in_E);

}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::writeMaxCliques()
 *   Write out the max cliques of the three partitions.
 *
 *
 * Preconditions:
 *   The maxclique variables Pcliques, Ccliques, and Ecliques must be instantiated!!
 *
 * Postconditions:
 *   Information about the maxclique variables is written out.
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
GMTemplate::
writeMaxCliques(oDataStreamFile& os)
{
  os.nl();
  os.writeComment("---\n");
  os.writeComment("---- P Partitions Cliques and their weights\n");
  P.writeMaxCliques(os);

  os.nl();
  os.writeComment("---\n");
  os.writeComment("---- C Partitions Cliques and their weights\n");
  C.writeMaxCliques(os);

  os.nl();
  os.writeComment("---\n");
  os.writeComment("---- E Partitions Cliques and their weights\n");
  E.writeMaxCliques(os);

  // and write out all the clique information to the file as comments.
  {
    os.nl();
    os.nl();
    os.writeComment("----\n");
    os.writeComment("----\n");
    os.writeComment("----\n");
    os.writeComment("---- Printing final clique set and clique weights---\n");

    double p_maxWeight = -1.0;
    double p_totalWeight = -1.0; // starting flag
    os.writeComment("  --- Prologue summary, %d cliques\n",P.cliques.size());
    for (unsigned i=0;i<P.cliques.size();i++) {
      double curWeight = MaxClique::computeWeight(P.cliques[i].nodes);
      os.writeComment("   --- P curWeight = %f\n",curWeight);
      if (curWeight > p_maxWeight) p_maxWeight = curWeight;
      if (p_totalWeight == -1.0)
	p_totalWeight = curWeight;
      else
	p_totalWeight = p_totalWeight + log10(1+pow(10,curWeight-p_totalWeight));
    }
    os.writeComment("  --- Prologue max clique weight = %f, total weight = %f, jt weight = %f\n",
	   p_maxWeight,p_totalWeight,
		    JunctionTree::junctionTreeWeight(P.cliques,PCInterface_in_P));

    double c_maxWeight = -1.0;
    double c_totalWeight = -1.0; // starting flag
    os.writeComment("  --- Chunk summary, %d cliques\n",C.cliques.size());
    for (unsigned i=0;i<C.cliques.size();i++) {
      double curWeight = MaxClique::computeWeight(C.cliques[i].nodes);
      os.writeComment("   --- C curWeight = %f\n",curWeight);
      if (curWeight > c_maxWeight) c_maxWeight = curWeight;
      if (c_totalWeight == -1.0)
	c_totalWeight = curWeight;
      else
	c_totalWeight = c_totalWeight + log10(1+pow(10,curWeight-c_totalWeight));
    }
    os.writeComment("  --- Chunk max clique weight = %f, total Cx%d weight = %f, per-chunk total C weight = %f, jt weight = %f\n",
	   c_maxWeight,
	   S,
	   c_totalWeight,
	   c_totalWeight - log10((double)S),
           JunctionTree::junctionTreeWeight(C.cliques,CEInterface_in_C));


    double e_maxWeight = -1.0;
    double e_totalWeight = -1.0; // starting flag
    os.writeComment("  --- Epilogue summary, %d cliques\n",E.cliques.size());
    for (unsigned i=0;i<E.cliques.size();i++) {
      double curWeight = MaxClique::computeWeight(E.cliques[i].nodes);
      os.writeComment("   --- E curWeight = %f\n",curWeight);
      if (curWeight > e_maxWeight) e_maxWeight = curWeight;
      if (e_totalWeight == -1.0)
	e_totalWeight = curWeight;
      else
	e_totalWeight = e_totalWeight + log10(1+pow(10,curWeight-e_totalWeight));
    }
    const set <RandomVariable*> emptySet;
    os.writeComment("  --- Epilogue max clique weight = %f, total weight = %f, jt_weight = %f\n",
	   e_maxWeight,e_totalWeight,
           JunctionTree::junctionTreeWeight(E.cliques,emptySet));

    double maxWeight
      = (p_maxWeight>c_maxWeight?p_maxWeight:c_maxWeight);
    maxWeight =
      (maxWeight>e_maxWeight?maxWeight:e_maxWeight);
    double totalWeight = p_totalWeight;
    // log version of: totalWeight += c_totalWeight
    totalWeight += log10(1+pow(10,c_totalWeight-totalWeight));
    // log version of: totalWeight += e_totalWeight
    totalWeight += log10(1+pow(10,e_totalWeight-totalWeight));

    os.writeComment("--- Final set (P,Cx%d,E) has max clique weight = %f, total state space = %f ---\n",
	   S,
	   maxWeight,
	   totalWeight);

    // print out a couple of total state spaces for various unrollings
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+S-1,totalWeight);

    totalWeight += log10(1+pow(10,c_totalWeight-totalWeight));	
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+2*S-1,totalWeight);

    totalWeight += log10(1+pow(10,log10(3.0) + c_totalWeight-totalWeight));
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+5*S-1,totalWeight);

    totalWeight += log10(1+pow(10,log10(5.0) + c_totalWeight-totalWeight));
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+10*S-1,totalWeight);

    totalWeight += log10(1+pow(10,log10(10.0) + c_totalWeight-totalWeight));
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+20*S-1,totalWeight);

    totalWeight += log10(1+pow(10,log10(30.0) + c_totalWeight-totalWeight));
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+50*S-1,totalWeight);

    totalWeight += log10(1+pow(10,log10(50.0) + c_totalWeight-totalWeight));
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+100*S-1,totalWeight);

    totalWeight += log10(1+pow(10,log10(400.0) + c_totalWeight-totalWeight));
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+500*S-1,totalWeight);

    totalWeight += log10(1+pow(10,log10(500.0) + c_totalWeight-totalWeight));
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+1000*S-1,totalWeight);

  }

}




/*-
 *-----------------------------------------------------------------------
 * GMTemplate::readMaxCliques()
 *   Read in the max cliques of the three partitions. 
 *
 * Preconditions:
 *   The maxclique variables Pcliques, Ccliques, and Ecliques must be instantiated!!
 *
 * Postconditions:
 *   Information about the maxclique variables is written out.
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
GMTemplate::
readMaxCliques(iDataStreamFile& is)
{
  P.readMaxCliques(is);
  C.readMaxCliques(is);
  E.readMaxCliques(is);
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::triangulatePartitionsByCliqueCompletion()
 *   Triangulate the partitions by completing the cliques that have been read in.
 *
 * Preconditions:
 *   The partition variables P, C, and E must be instantiated.
 *   The maxclique variables Pcliques, Ccliques, and Ecliques must be instantiated!!
 *
 * Postconditions:
 *   the variables pointed to by the cliques will be made complete.
 *
 * Side Effects:
 *   Variables pointed to by Pcliques, Ccliques, and Ecliques will have their
 *   neighbors adjusted.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
triangulatePartitionsByCliqueCompletion()
{
  P.triangulatePartitionsByCliqueCompletion();
  C.triangulatePartitionsByCliqueCompletion();
  E.triangulatePartitionsByCliqueCompletion();
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::computeUnrollParamaters()
 *   compute the unrolling amount for the number of frames in an observation matrix.
 *
 *
 * Preconditions:
 *   None.
 *
 * Postconditions:
 *   None.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   computes the correct unrolling amount in non-const arguments.
 *   Returns 'true' if success and 'false' otherwise.
 *
 *-----------------------------------------------------------------------
 */
bool
GMTemplate::
computeUnrollParamaters(const unsigned numFrames,
			unsigned& basicTemplateUnrollAmount,
			unsigned& modifiedTemplateUnrollAmount,
			unsigned& numUsableFrames,
			unsigned& frameStart,
			const JustifyType justifyType)

{
  // useful facts:
  // T = number of frames.
  const unsigned T = numFrames;
  // p = number of frames in basic template prologue (as read by file parser)
  const int p = fp.numFramesInP();
  // c = number of frames in basic template chunk (as read by file parser)
  const int c = fp.numFramesInC();
  // e = number of frames in basic template epilogue (as read by file parser)
  const int e = fp.numFramesInE();

  // we have relationships:
  // from basic template:
  //     T = p + i*c + e  for some i = positive integer.
  // from modified template:
  //     T = p + e + (M+jS)*c for j = positive integer.
  // check both here:
  if ( T < p + e + (M+S)*c )
    return false;
  // Ok, it's possible to use this sentence.
  // The above puts constraints:
  //    (T - p - e) = i*c
  //    (T - p - e) = (M+jS)*c = l*c, with l = (M+jS)
  // which means that we must have:
  //    (T - p - e - M*c) = jS*c  for j = positive integer.

  // compute T-p-e-M*c % S*c, the remainder.
  const unsigned remainder = (T-p-e-M*c) % (S*c);
  // number of usable frames subtracts this off, assign to
  // T' = numUsableFrames
  numUsableFrames = numFrames - remainder;
  // since T' - p - e - M*c = k*S*c for positive k,
  // then T' - p - e = M*c +k*S*c = (M+kS)*c = i*c, i = (M+kS)
  // as asked for above, so T' is a perfect number to use.

  if (justifyType == leftJustify)
    frameStart = 0;
  else if (justifyType == rightJustify)
    frameStart = remainder;
  else // center
    frameStart = remainder/2; // (round to left if not even)
  
  // all else should work
  assert ( (numUsableFrames - p - e) % c == 0 );
  basicTemplateUnrollAmount = (numUsableFrames-p-e)/c - 1;
  assert ( ((numUsableFrames-p-e)/c - M) % S == 0 );
  modifiedTemplateUnrollAmount = ((numUsableFrames-p-e)/c - M)/S - 1;

  return true;
}
