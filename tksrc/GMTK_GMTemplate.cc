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
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#if 0
#include "GMTK_ObservationMatrix.h"
#endif
#include "GMTK_JunctionTree.h"
#include "GMTK_Partition.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)


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


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::setUpClonedPartitionGraph()
 *
 *   Given a P,C, and E that has been cloned without parents, this
 *   will set up each partition (P,C, and E) as follows: only the nodes that are
 *   in the corresponding current partition set are pointed to by
 *   the resuling neighbors member in each RV in the partition. 
 *   In other words, neighbors of variables in a parititon never
 *   point to nodes outside of the partitions. Parents and children,
 *   however, of nodes in the partition, might very well point
 *   to variables outside of the current partition. See also
 *   comment below in main body of this routine.
 *
 * Preconditions:
 *   P, C, and E are the partitions to be cloned. 
 *
 * Postconditions:
 *   Pc,Cc,Ec are the cloned partitions of P,C, and E respectively, but
 *   with the properties as mentioned above. Also, the P_in_to_out,
 *   C_in_to_out, E_in_to_out, provide the mappings from P to Pc,
 *   C to Cc, and E to Ec respectively.
 *   
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     returns the cloned partitions and the in_to_out mappings
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
setUpClonedPartitionGraph(const set<RV*>& P,
			  const set<RV*>& C,
			  const set<RV*>& E,
			  // cloned variables
			  set<RV*>& Pc,
			  set<RV*>& Cc,
			  set<RV*>& Ec,
			  // next 3 should be const but there is no "op[] const"
			  map < RV*, RV* >& P_in_to_out,
			  map < RV*, RV* >& C_in_to_out,
			  map < RV*, RV* >& E_in_to_out)
{

  // These next few routine calls set the neighbors of the variables
  // in the output arguments (i..e, PC, Cc, and Ec) to be the
  // correctly associated variables, but those neighbors will not
  // include neighbors that are not in the current partition (i.e.,
  // they will be dissociated with any other possible portion of the
  // network).  This means that the rvs themselves in the intersection
  // between partitions need to be unique. Meaning, the intersection I
  // = intersection(P,C) is contained both in P, and C, but since P[I]
  // should have neighbors only in P and C[I] should have neighbors
  // only in C, we must use different actuall C++ random variables for
  // P[I] and C[I]. The reason for this is the triangulation code
  // which uses the neighbors of random variables as isolated
  // partitions (i.e., each partition triangulation is within a
  // partition, and has nothing to do with the triangulations of other
  // partitions) --- specifically, we triangulate each partition
  // separately, and after triangulation, the neighbors member of
  // variables in P[I] will not be the same as the corresponding
  // neighbors member for variables in C[I] (since they are
  // triangulated separately).  Indeed, this may seem a little weird,
  // since a variable v in C with a parent p in P will not have p as
  // v's neighbor.

  // Note further that the variables PCInterface_in_P,
  // PCInterface_in_C, CEInterface_in_C, and CEInterface_in_E will
  // have the unique interface variables in each partition since they
  // will sometimes be quite useful. I.e., even though
  // PCInterface_in_P and PCInterface_in_C correspond to the same
  // actual variables, they are STL sets with pointers to different
  // C++ objects (so their STL set intersection will in fact be
  // empty).


  cloneRVShell(P,Pc,P_in_to_out);
  cloneRVShell(C,Cc,C_in_to_out);
  cloneRVShell(E,Ec,E_in_to_out);
  // Note that now, Pc, Cc, and Ec consist of completely separate C++
  // RV objects, meaning that from STL's point of view, there is no
  // intersection between any of them (i..e, C++intersection(Pc,Cc) =
  // empty, etc.). Of course, in terms of actual real RVs, there is an
  // intersection, namely the interface variables.

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
 *   the other partitions whatever they are, called O1 and O2. 
 *   
 *   It then sets the neighbors structures for all of Sc from S. The
 *   neighbors are variables that are forced to be part of the
 *   partition S itself. The reason for this is that we triangulate
 *   each partition separately, so there must not be any neighbors
 *   edges pointing into adjacent partitions. Note, however, that
 *   partitions will of course have overlap in nodes (i.e., the
 *   interface nodes).
 *   
 *   It then sets the parents variables for each variable in Sc. The
 *   parents (and also the children) might NOT be fully contained in
 *   Sc, so it needs to use the mappings for O1 and O2 to get the
 *   location of the correspondly cloned variables for those
 *   partitions. In other words, when done the parents and children
 *   very well might point into other partitions.
 *
 *   To summarize
 *        1) Sc will be a cloned graph of S (same RVs)
 *        2) Sc will have undirected neighbors set to point only within Sc (not outside)
 *        3) Sc will have parents/childrens possibly pointing outside of Sc. If
 *           the parent/child lives in Sc, then the parent/child poitner will point
 *           to within Sc, but if the parent/child does not live within Sc, it will
 *           point to the corresponding variable in O1 or O2.
 *
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
setPartitionParentsChildrenNeighbors(const set<RV*>& S,
				     set<RV*>& Sc,
				     // next 3 should be const but there is no "op[] const" in STL
				     map < RV*, RV* >& S_in_to_out,
				     map < RV*, RV* >& O1_in_to_out,
				     map < RV*, RV* >& O2_in_to_out)
{

  for (set<RV*>::iterator i=S.begin();i != S.end(); i++) {

    RV*rv = (*i);

    // first set up new neighbors for S_in_to_out[rv]
    // Note: Neighbors are defined to point ONLY TO OTHER VARIABLES
    // IN THE SET S_in_to_out (i.e., the members of the partition). Any
    // other neighbors are not contained in that set are not included.
    set<RV*> tmp;
    for (set<RV*>::iterator j = rv->neighbors.begin();
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
    vector<RV *> sParents;
    for (unsigned l=0;l<rv->rv_info.switchingParents.size(); l++) {
      // grab a copy for readability
      RV* const par = rv->switchingParentsVec()[l];
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
    vector< vector < RV* > > cParentsList;
    cParentsList.resize(rv->rv_info.conditionalParents.size());
    for (unsigned l=0;l<rv->rv_info.conditionalParents.size();l++) {
      for (unsigned m=0;m<rv->condParentsVec(l).size();m++) {
	// grab a copy for readability
	RV* const par = rv->condParentsVec(l)[m];
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
    // this next routine sets both parents and children.
    S_in_to_out[rv]->setParents(sParents,cParentsList);
  }
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::cloneWithtoutParents()
 *   Clone a set of random variables from 'in' to 'out'. The cloned
 *   variables have entirely empty parents, children, and neighbors structures.
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
 *     returns the new set of RVs, and the input-to-output mapping as a STL map.
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
cloneRVShell(const set<RV*>& in, 
	     set<RV*>& out,
	     map < RV*, RV* >& in_to_out)
{
  in_to_out.clear();
  out.clear();

  for (set<RV*>::iterator i=in.begin();
       i != in.end(); i++) {

    // sanity check, to ensure a node is not its own neighbor
    assert ( (*i)->neighbors.find((*i)) == (*i)->neighbors.end() );

    RV*rv = (*i)->cloneRVShell();
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
  os.writeComment("GMTK Triangulation File Created: %s\n",buff);
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
  // printf("about to write boundary method = %s\n",boundaryMethod.c_str());
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
  for (set<RV*>::iterator i=P.nodes.begin();
       i != P.nodes.end(); i++) {
    RV* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RV*>::iterator j=rv->neighbors.begin();
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
  for (set<RV*>::iterator i = P.nodes.begin();
       i != P.nodes.end(); i++) {
    RV* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();

  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- C partition information: variables and their neighbors\n");
  buffer.clear();
  for (set<RV*>::iterator i=C.nodes.begin();
       i != C.nodes.end(); i++) {
    RV* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RV*>::iterator j=rv->neighbors.begin();
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
  for (set<RV*>::iterator i = C.nodes.begin();
       i != C.nodes.end(); i++) {
    RV* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- E partition information: variables and their neighbors\n");
  buffer.clear();
  for (set<RV*>::iterator i=E.nodes.begin();
       i != E.nodes.end(); i++) {
    RV* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RV*>::iterator j=rv->neighbors.begin();
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
  for (set<RV*>::iterator i = E.nodes.begin();
       i != E.nodes.end(); i++) {
    RV* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- PC information : variables and their neighbors\n");
  buffer.clear();
  for (set<RV*>::iterator i=PCInterface_in_C.begin();
       i != PCInterface_in_C.end(); i++) {
    RV* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RV*>::iterator j=rv->neighbors.begin();
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
  for (set<RV*>::iterator i = PCInterface_in_C.begin();
       i != PCInterface_in_C.end(); i++) {
    RV* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- CE information : variables and their neighbors\n");
  buffer.clear();
  for (set<RV*>::iterator i=CEInterface_in_C.begin();
       i != CEInterface_in_C.end(); i++) {
    RV* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RV*>::iterator j=rv->neighbors.begin();
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
  for (set<RV*>::iterator i = CEInterface_in_C.begin();
       i != CEInterface_in_C.end(); i++) {
    RV* rv = (*i);
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
    error("ERROR: reading file '%s' line %d, M (number of chunks in which to find interface boundary) must be >= 1\n",
	  is.fileName(),is.lineNo());

  if (M == GMTEMPLATE_UNINITIALIZED_MS) {
    // This is a const cast hack to get around the fact that member M
    // is declared constant, in cases where we create an object with
    // an uninitialized M.  Perhaps the right thing to do is to
    // undeclare M constant.
    unsigned* Mp = (unsigned*)&M;
    *Mp = loc_M;
  } else {
    if (loc_M != M)
      error("ERROR: reading file '%s' line %d, M (=%d) given in tri-file does not equal %d\n",
	    is.fileName(),is.lineNo(),loc_M,M);
  }

  is.read(loc_S,"S value");
  if (loc_S == 0)
    error("ERROR: reading file '%s' line %d, S (chunk skip) must be >= 1\n",
	  is.fileName(),is.lineNo());

  if (S == GMTEMPLATE_UNINITIALIZED_MS) {
    unsigned *Sp = (unsigned*)&S;
    *Sp = loc_S;
  } else {
    if (loc_S != S)
      error("ERROR: reading file '%s' line %d, S in file (%d) does not equal %d\n",
	    is.fileName(),is.lineNo(),loc_S,S);
  }

  // interface method
  is.read(loc_I,"interface method value");
  if (loc_I == "LEFT") {
    leftInterface = true;
  } else if  (loc_I == "RIGHT") {
    leftInterface = false;
  } else {
      error("ERROR: reading file '%s' line %d, interface in file must be 'LEFT' or 'RIGHT' but got string '%s'\n",
	    is.fileName(),is.lineNo(),loc_I.c_str());
  }

  // read in information about method used to create current boundary
  is.read(boundaryMethod,"boundary method string");

  vector <RV*> unrolled_rvs;
  map < RVInfo::rvParent, unsigned > positions;
  fp.unroll(M+S-1,unrolled_rvs,positions);

  // need to moralize.
  for (unsigned i=0;i<unrolled_rvs.size();i++) {
    unrolled_rvs[i]->createNeighborsFromParentsChildren();
  }
  // add edges from any extra factors in .str file
  fp.addUndirectedFactorEdges(M+S-1,unrolled_rvs,positions);
  for (unsigned i=0;i<unrolled_rvs.size();i++) {
    unrolled_rvs[i]->moralize();    
  }

  // create temporary local variables.
  set<RV*> loc_P;
  set<RV*> loc_C;
  set<RV*> loc_E;

  unsigned setSize;
  string str_tmp;

  is.read(str_tmp,"P partition name");
  if (str_tmp != P_partition_name)
    error("ERROR: P partition information in file '%s' line %d is invalid for given graph structure\n",
	  is.fileName(),is.lineNo());
  is.read(setSize,"P partition set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: P partition information in file '%s' line %d is invalid for given graph structure\n",
	    is.fileName(),is.lineNo());
    loc_P.insert(unrolled_rvs[(*loc).second]);
  }

  is.read(str_tmp,"C partition name");
  if (str_tmp != C_partition_name)
    error("ERROR: C partition information in file '%s' line %d is invalid for given graph structure\n",
	  is.fileName(),is.lineNo());
  is.read(setSize,"C partition set size");
  if (setSize == 0)
    error("ERROR: C partition information in file '%s' line %d specifies no variables. C partition must not be empty.\n",
	  is.fileName(),is.lineNo());
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: C partition information in file '%s' line %d is invalid for given graph structure\n",
	    is.fileName(),is.lineNo());
    loc_C.insert(unrolled_rvs[(*loc).second]);
  }



  is.read(str_tmp,"E partition name");
  if (str_tmp != E_partition_name)
    error("ERROR: E partition information in file '%s' line %d is invalid for given graph structure\n",
	  is.fileName(),is.lineNo());
  is.read(setSize,"E partition set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: E partition information in file '%s' line %d is invalid for given graph structure\n",
	    is.fileName(),is.lineNo());
    loc_E.insert(unrolled_rvs[(*loc).second]);
  }

  //////////////////////////////////////////////
  // next, read in the interface definitions. //
  //////////////////////////////////////////////
  set<RV*> loc_PCInterface;
  set<RV*> loc_CEInterface;

  // get PC interface
  is.read(str_tmp,"PC interface name");
  if (str_tmp != PC_interface_name)
    error("ERROR: PC interface information in file '%s' line %d is invalid for given graph structure\n",
	  is.fileName(),is.lineNo());
  is.read(setSize,"PC interface set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: PC interface information in file '%s' line %d is invalid for given graph structure\n",
	    is.fileName(),is.lineNo());
    loc_PCInterface.insert(unrolled_rvs[(*loc).second]);
  }

  // get CE interface
  is.read(str_tmp,"CE interface name");
  if (str_tmp != CE_interface_name)
    error("ERROR: CE interface information in file '%s' line %d is invalid for given graph structure\n",
	  is.fileName(),is.lineNo());
  is.read(setSize,"CE interface set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: CE interface information in file '%s' line %d is invalid for given graph structure\n",
	    is.fileName(),is.lineNo());
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
 *   The internal partition and interface variables are created. 
 *   Properties of the internal partitions are as follows:
 *     1) Each set of variasbles in each parititon are different C++ RVs, even the
 *        common interfaces are distinct C++ variables (but they correspond to the same
 *        actual RVs).
 *     2) All sets of interface variables in all partitions are completed.
 *     3) In each partition, the RV's 'neighbors strucures will point only to
 *        other RVs in its partition, and will not point outside the partition.
 *     4) parents/children will point to current partition if the parent/child
 *        lives in the current partiiton, but otherwise, those pointers will
 *        point to a corresponding variable in another partitiion.
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
createPartitions(const set<RV*>& arg_P,
		 const set<RV*>& arg_C,
		 const set<RV*>& arg_E,
		 const set<RV*>& arg_PCInterface,
		 const set<RV*>& arg_CEInterface)
{

  // delete old stuff
  clear();

  /////////////////////////////////////////////////////////////////////
  // create a new variable set for each, make the interfaces
  // complete, and finish up.

  map < RV*, RV* > P_in_to_out;
  map < RV*, RV* > C_in_to_out;
  map < RV*, RV* > E_in_to_out;

  setUpClonedPartitionGraph(arg_P,arg_C,arg_E,P.nodes,C.nodes,E.nodes,P_in_to_out,C_in_to_out,E_in_to_out);

  // complete the PC interface in P
  PCInterface_in_P.clear();
  for (set<RV*>::iterator i=arg_PCInterface.begin();
       i != arg_PCInterface.end();i++) {
    PCInterface_in_P.insert(P_in_to_out[(*i)]);
  }
  MaxClique::makeComplete(PCInterface_in_P);

  // complete the PC interface in C
  PCInterface_in_C.clear();
  for (set<RV*>::iterator i=arg_PCInterface.begin();
       i != arg_PCInterface.end();i++) {
    PCInterface_in_C.insert(C_in_to_out[(*i)]);
  }
  MaxClique::makeComplete(PCInterface_in_C);

  // complete the CE interface in C
  CEInterface_in_C.clear();
  for (set<RV*>::iterator i=arg_CEInterface.begin();
       i != arg_CEInterface.end();i++) {
    CEInterface_in_C.insert(C_in_to_out[(*i)]);
  }
  MaxClique::makeComplete(CEInterface_in_C);  

  // complete the CE interface in E
  CEInterface_in_E.clear();
  for (set<RV*>::iterator i=arg_CEInterface.begin();
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
  writePMaxCliques(os);
  writeCMaxCliques(os);
  writeEMaxCliques(os);
  writeCliqueInformation(os);
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::writePMaxCliques()
 *   Write out the max cliques of the P partition.
 *
 *
 * Preconditions:
 *   The maxclique variable Pcliques must be instantiated!!
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
writePMaxCliques(oDataStreamFile& os)
{
  os.nl();
  os.writeComment("---\n");
  os.writeComment("---- P Partitions Cliques and their weights\n");
  P.writeMaxCliques(os);
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::writeCMaxCliques()
 *   Write out the max cliques of the C partition.
 *
 *
 * Preconditions:
 *   The maxclique variable Ccliques must be instantiated!!
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
writeCMaxCliques(oDataStreamFile& os)
{
  os.nl();
  os.writeComment("---\n");
  os.writeComment("---- C Partitions Cliques and their weights\n");
  C.writeMaxCliques(os);
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::writeEMaxCliques()
 *   Write out the max cliques of the E partition.
 *
 *
 * Preconditions:
 *   The maxclique variable Ecliques must be instantiated!!
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
writeEMaxCliques(oDataStreamFile& os)
{
  os.nl();
  os.writeComment("---\n");
  os.writeComment("---- E Partitions Cliques and their weights\n");
  E.writeMaxCliques(os);
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::writeCliqueInformation()
 *   Write out information in comment form for the max cliques of the three partitions.
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
writeCliqueInformation(oDataStreamFile& os)
{
  // write out all the clique information to the file as comments.
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
	p_totalWeight = log10add(curWeight,p_totalWeight);
    }
    os.writeComment("  --- Prologue max clique weight = %f, total weight = %f, jt_weight = %f\n",
	   p_maxWeight,p_totalWeight,
		    JunctionTree::junctionTreeWeight(P.cliques,
						     PCInterface_in_P,
						     NULL,&PCInterface_in_P));


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
	c_totalWeight = log10add(curWeight,c_totalWeight);
    }
    os.writeComment("  --- Chunk max clique weight = %f, total Cx%d weight = %f, per-chunk total C weight = %f, jt_weight = %f\n",
	   c_maxWeight,
	   S,
	   c_totalWeight,
	   c_totalWeight - log10((double)S),
           JunctionTree::junctionTreeWeight(C.cliques,
					    CEInterface_in_C,
					    &PCInterface_in_C,&CEInterface_in_C));


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
	e_totalWeight = log10add(curWeight,e_totalWeight);
    }
    const set <RV*> emptySet;
    os.writeComment("  --- Epilogue max clique weight = %f, total weight = %f, jt_weight = %f\n",
	   e_maxWeight,e_totalWeight,
           JunctionTree::junctionTreeWeight(E.cliques,
					    emptySet,
					    &CEInterface_in_E,NULL));

    double maxWeight
      = (p_maxWeight>c_maxWeight?p_maxWeight:c_maxWeight);
    maxWeight =
      (maxWeight>e_maxWeight?maxWeight:e_maxWeight);
    double totalWeight = p_totalWeight;
    // log version of: totalWeight += c_totalWeight
    totalWeight = log10add(c_totalWeight,totalWeight);
    // log version of: totalWeight += e_totalWeight
    totalWeight = log10add(e_totalWeight,totalWeight);

    os.writeComment("--- Final set (P,Cx%d,E) has max clique weight = %f, total state space = %f ---\n",
	   S,
	   maxWeight,
	   totalWeight);

    // print out a couple of total state spaces for various unrollings
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+S-1,totalWeight);

    totalWeight = log10add(c_totalWeight,totalWeight);	
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+2*S-1,totalWeight);

    totalWeight = log10add(log10(3.0) + c_totalWeight,totalWeight);
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+5*S-1,totalWeight);

    totalWeight = log10add(log10(5.0) + c_totalWeight,totalWeight);
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+10*S-1,totalWeight);

    totalWeight = log10add(log10(10.0) + c_totalWeight,totalWeight);
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+20*S-1,totalWeight);

    totalWeight = log10add(log10(30.0) + c_totalWeight,totalWeight);
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+50*S-1,totalWeight);

    totalWeight = log10add(log10(50.0) + c_totalWeight,totalWeight);
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+100*S-1,totalWeight);

    totalWeight = log10add(log10(400.0) + c_totalWeight,totalWeight);
    os.writeComment("--- Total weight when unrolling %dx = %f ---\n",M+500*S-1,totalWeight);

    totalWeight = log10add(log10(500.0) + c_totalWeight,totalWeight);
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
  // C can't be empty.
  if (C.cliques.size() == 0)
    error("ERROR: reading file '%s' near line %d. Number of cliques in the C partition must be >= 1\n",
	  is.fileName(),is.lineNo());
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
 * GMTemplate::computeUnrollParameters()
 *
 *   compute the unrolling amount in various forms for the number of
 *   frames in an observation matrix.
 *
 *
 * Preconditions:
 *   None.
 *
 * Postconditions:
 *   None. Note that modifiedTemplateUnrollAmount might be (-1) which means
 *   that this corresponds to the basic template P C E  (or equivalently, P' E').
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
computeUnrollParameters(const unsigned numFrames,
			// for numFrames, the max basic template
			// unroll amount that corresponds to this numFrames.
			unsigned& basicTemplateMaxUnrollAmount,
			// for numFrames, the min amount by which we
			// can unroll and still re-use the rvs accross
			// time.
			unsigned& basicTemplateMinUnrollAmount,
			// This next two arguments are ints, since its
			// value can go as low as -1. A return value
			// of -1 means that we just use [P' E'] for
			// the modified template.
			int& modifiedTemplateMaxUnrollAmount,
			int& modifiedTemplateMinUnrollAmount,
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
  // 
  // From basic template, we need to find the largest positive i such that:
  //     T >= p + i*c + e  
  // where i = *positive* integer.
  // Then i is the basic template unrolling amount.
  // 
  // From modified template, we need to find the largest non-negative j such that
  //     T >= p + e + (M+jS)*c 
  // where j = non-negative integer (meaning minimum length is p+e+M*c
  // with M=1 this is the basic template).
  // Then j is the modified template unrolling amount.

  // First, we check to make sure that there exists such an i and j.
  // Since M>=1 and S>=1, we can check both here by making sure that T is at
  // least p+e+M*c
  if ( T < p + e + M*c ) {
    infoMsg(Info,"Unrolling Incompatibility: Trying to unroll segment with only %d frames using template of [P=%d,C=%d,E=%d] frames, but with M=%d,S=%d, minimum segment frame length is P+M*C+E=%d\n",T,p,c,e,M,S,p+e+M*c );
    return false;
  }

  // Ok, it's possible to use this segment.
  // The above puts constraints:
  //    (T - p - e) >= i*c,  i >= 1 
  //    (T - p - e) >= (M+jS)*c = l*c, with l = (M+jS), j>=0
  // which means that we must have:
  //    (T - p - e - M*c) >= jS*c  for j = non-negative integer.

  // compute (T-p-e-M*c) % (S*c), the remainder to find the
  // lenth T' = numUsableFrames that we can actually use.
  // Compute the number of frames to subtract off from T.
  const unsigned remainder = (T-p-e-M*c) % (S*c);
  numUsableFrames = numFrames - remainder;

  // Now, we have that:
  // since T' - p - e - M*c = k*S*c for some positive k,
  // then T' - p - e = M*c +k*S*c = (M+kS)*c = i*c, i = (M+kS)
  // as asked for above, so T' is the right number to use.

  // Since T might not be equal to T', we also return
  // a frameStart depending on the justify argument. This can be
  // usd as an observation offset.
  if (justifyType == leftJustify)
    frameStart = 0;
  else if (justifyType == rightJustify)
    frameStart = remainder;
  else // center
    frameStart = remainder/2; // (round to left if not even)
  
  // Finally, compute i = basicTemplateMaxUnrollAmount, and
  // j = modifiedTemplateMaxUnrollAmount to return.
  assert ( (numUsableFrames - p - e) % c == 0 );
  basicTemplateMaxUnrollAmount = (numUsableFrames-p-e)/c - 1;
  assert ( ((numUsableFrames-p-e)/c - M) % S == 0 );
  // solve for j in (T' - p - e) = (M+jS)*c
  modifiedTemplateMaxUnrollAmount = ((numUsableFrames-p-e)/c - M)/S - 1;

  // The number of actual C' partitions that will be used in this segment.
  const int numCoPartitions = modifiedTemplateMaxUnrollAmount+1;
  // The amount by which we should unroll the basic template.
  if (numCoPartitions == 0) {
    // then we use the non-unrolled basic template, i.e., graph
    // exactly in the str file. I.e., we should have random variables
    // for no replications of C', giving [P' E'].
    basicTemplateMinUnrollAmount = 0;
    modifiedTemplateMinUnrollAmount = -1;
  } else if  (numCoPartitions == 1) {
    // then we use the non-unrolled modified template. I.e., we should have random
    // variables for one replication of C', giving [P' C' E'].
    basicTemplateMinUnrollAmount = M + S - 1;
    modifiedTemplateMinUnrollAmount = 0;
  } else {
    // unroll the basic template by an amount that is equivalent to unrolling
    // the modified template by 1. I.e., we should have random
    // variables for two replications of C', giving [P' C' C' E'].
    basicTemplateMinUnrollAmount = M + 2*S - 1;
    modifiedTemplateMinUnrollAmount = 1;
  }

  return true;
}
