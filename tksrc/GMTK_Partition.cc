/*-
 * GMTK_Partition.cc
 *    Basic Partition for a given graph file.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2009, < fill in later >
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
#include "GMTK_ObservationMatrix.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_Partition.h"

VCID("$Header$")

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
Partition::Partition(Partition& from_part,
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
  // NOTE: It is Crucial for the cliques in the cloned partition to be
  // inserted in the *SAME ORDER* as in the partition being cloned.
  for (unsigned i=0;i<from_part.cliques.size();i++) {
    cliques.push_back(MaxClique(from_part.cliques[i],
				newRvs,ppf,frameDelta));
  }
}
#endif

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
  is.read(numCliques,"number of cliques");

#if 0
  // remove check for numCliques being > 0 since we now allow for empty partitions.
  if (numCliques == 0)
    error("ERROR: reading file '%s' line %d, numCliques must be >= 1\n",
	  is.fileName(),is.lineNo());
#endif

  // create a map for easy access to set of nodes
  map < RVInfo::rvParent, RV* > namePos2Var;
  for (set<RV*>::iterator i=nodes.begin();
       i != nodes.end(); i++) {
    RV* rv = (*i);
    RVInfo::rvParent par;
    par.first = rv->name();
    par.second = rv->frame();
    namePos2Var[par] = rv;
  }

  for (unsigned i=0;i<numCliques;i++) {
    set<RV*> clique;
    
    unsigned cliqueNo;
    is.read(cliqueNo,"clique number value");
    if (cliqueNo != i)
      error("ERROR: reading file %s, line %d, bad cliqueNo (= %d) when reading cliques, out of sequence, should be = %d instead.\n",
	    is.fileName(),is.lineNo(),cliqueNo,i);
    
    unsigned cliqueSize;
    is.read(cliqueSize,"clique size value");

#if 0
    // remove check for min clique size of 1.
    if (cliqueSize <= 1)
      error("ERROR: reading file %s line %d, cliqueSize %d must be >= 2\n",
	    is.fileName(),is.lineNo(),cliqueSize);
#endif    
    if (cliqueSize <= 0)
      error("ERROR: reading file %s line %d, reading clique number %d, but clique size %d must be >= 1\n",
	    is.fileName(),is.lineNo(),i,cliqueSize);


    for (unsigned j=0;j<cliqueSize;j++) {

      RVInfo::rvParent par;
      is.read(par.first,"parent name");
      is.read(par.second,"parent position");

      map < RVInfo::rvParent, RV* >::iterator loc;
      loc = namePos2Var.find(par);
      if (loc == namePos2Var.end())
	error("ERROR: reading file %s line %d, clique specification %d has %d'th variable %s(%d) that does not exist in partition.\n",
	      is.fileName(),is.lineNo(),i,j,par.first.c_str(),par.second);
      RV* rv = (*loc).second;
      clique.insert(rv);
    }
    cliques.push_back(MaxClique(clique));
  }

}


/*-
 *-----------------------------------------------------------------------
 * Partition::reportScoreStats()
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
Partition::reportScoreStats()
{
  for (unsigned i=0;i<cliques.size();i++) {
    printf("Clique %d:\n",i);
    cliques[i].reportScoreStats();
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



/*-
 *-----------------------------------------------------------------------
 * Partition::setCliquesFromAnotherPartition()
 *   Set the cliques from anohter partition. The other partition
 *   must be from the same structure file, and if the current partition is  a P (resp. C, E)
 *   than the other partition must also be a P (resp. C, and E). Also, the partitions
 *   are assumed to come from the same boundary for the .str file.
 *   The routine is used to merge together gm_templates for different triangulations
 *   of the same boundary but say a paralle triangulation of P, C, and E.
 *  
 *   Note that it is assumed that the different partitions refer to different instantiations of
 *   the same set of random variables (so we can't use rv1 == rv2, but instead must use
 *   name and frame equality).
 *
 * Preconditions:
 *   The corresponding partition must be instantiated with nodes.
 *   
 *
 * Postconditions:
 *   The clique variables now refer to the cliques in the other partition.
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
Partition::setCliquesFromAnotherPartition(Partition& from_part)
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

  // make sure nodes refer to same partition.
  for (it = from_part.nodes.begin();
       it != from_part.nodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame();

    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find rv %s(%d) in RV set of dest partition\n",
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


