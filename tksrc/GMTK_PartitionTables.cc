/*-
 * GMTK_JunctionTree.cc
 *     Junction Tree
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
#include <new>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_GMParms.h"


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

VCID("$Header$")


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        PartitionTables support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

PartitionTables::PartitionTables(JT_Partition& origin)
{

  // first allocate space with empty (and unusable) entries
  maxCliques = new MaxCliqueTable[origin.cliques.size()];
  separatorCliques = new ConditionalSeparatorTable[origin.separators.size()];

  // then actually re-construct the objects in the array appropriately.
  for (unsigned i=0;i<origin.cliques.size();i++) {
    new (&maxCliques[i]) 
      MaxCliqueTable(origin.cliques[i]);
  }
  for (unsigned i=0;i<origin.separators.size();i++) {
    new (&separatorCliques[i]) 
      ConditionalSeparatorTable(origin.separators[i]);
  }

}


/*-
 *-----------------------------------------------------------------------
 * PartitionTables::emIncrement()
 *
 *    Go through each clique in the partition and update the assigned probabiltiy
 *    variables for all entries in each clique, based on global probability of evidence
 *    given as the argument.
 *    If 'localNormalization' == true, then we ignore the evidence provided by probE
 *    and sum the clique first to get the local normalization constant.
 *
 * See Also:
 *    0) JunctionTree::collectEvidence()
 *    1) JunctionTree::distributeEvidence()
 *    2) JunctionTree::emIncrement()
 *       
 *
 * Preconditions:
 *    The cliques must be an instantiated table. All data structures must be set up.
 *
 * Postconditions:
 *   The accumulators are increment accordingly.
 *
 * Side Effects:
 *   This will update the accumulators of all trainable parameter objects.
 *
 * Results:
 *   None
 *
 *-----------------------------------------------------------------------
 */
void
PartitionTables::emIncrement(PartitionStructures& ps,
			     const logpr probE,
			     const bool localCliqueNormalization,
			     const double emTrainingBeam)
{
  for (unsigned cliqueNo=0;cliqueNo < ps.maxCliquesSharedStructure.size(); cliqueNo++ ) {
    maxCliques[cliqueNo].emIncrement(ps.maxCliquesSharedStructure[cliqueNo],
				     probE,localCliqueNormalization,emTrainingBeam);
  }
}


void PartitionTables::init(PartitionStructures& ps)
{
  for (unsigned cliqueNo=0;cliqueNo<ps.maxCliquesSharedStructure.size();cliqueNo++)
    maxCliques[cliqueNo].init(*(ps.maxCliquesSharedStructure[cliqueNo].origin));
  for (unsigned sepNo=0;sepNo<ps.separatorCliquesSharedStructure.size();sepNo++)
    separatorCliques[sepNo].init(*(ps.separatorCliquesSharedStructure[sepNo].origin));
}


/*-
 *-----------------------------------------------------------------------
 * PartitionTables::reportMemoryUsageTo()
 *   Reports memory usage of the template in ASCII format (suitable for stdout or an ascii file)
 *
 * Preconditions:
 *   The partitions must be validly instantiated with clique & separator structures.
 *
 * Postconditions:
 *   current memory usage is reported.
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
PartitionTables::reportMemoryUsageTo(PartitionStructures& ps,FILE *f)
{
  for (unsigned cliqueNo=0;cliqueNo < ps.maxCliquesSharedStructure.size(); cliqueNo++ ) {
    maxCliques[cliqueNo].reportMemoryUsageTo(*(ps.maxCliquesSharedStructure[cliqueNo].origin),f);
  }
  for (unsigned i=0;i<ps.separatorCliquesSharedStructure.size();i++) {
    separatorCliques[i].reportMemoryUsageTo(*(ps.separatorCliquesSharedStructure[i].origin),f);
  }
}




/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
