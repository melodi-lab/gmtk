/*-
 * GMTK_SparseJoinSectionTables.cc
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
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
#include <new>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_SectionIterator.h"
#include "GMTK_SparseJoinSectionTables.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
VCID(HGID)



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        SparseJoinSectionTables support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

SparseJoinSectionTables::SparseJoinSectionTables(JT_Partition& origin)
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


void 
SparseJoinSectionTables::projectToOutgoingSeparators(SectionIterator &stss_it,
						     PartitionStructures &sourceSectionStructures, 
						     ConditionalSeparatorTable *separatorTableArray,
						     ConditionalSeparatorTable::SharedLocalStructure &sepSharedStructure)
{

  // Here the section sending the message (projecting to outgoing separators) knows
  // it's using sparse join inference. So it's using SparseJoinSectionTables,
  // SparseJoinMaxCliques, etc. In particular, the section knows how to find its
  // cliques with outgoing separators (formerly just the single right interface clique),
  // how many MaxCliques are in the section, and the MaxCliqueTable::SharedLocalStructure
  // of each of those cliques...

  // for each maxCliques[i] with an outgoing separator
  unsigned i = stss_it.prev_ri();
  //   for each of maxClique[i]'s outgoing separators j
          maxCliques[i].ceSendToOutgoingSeparator(sourceSectionStructures.maxCliquesSharedStructure[i],
						  *separatorTableArray /* [j] */, sepSharedStructure);
}


void 
SparseJoinSectionTables::receiveBackwardsSeparators(SectionIterator &stss_it,
						    PartitionStructures &sourceSectionStructures, 
						    ConditionalSeparatorTable *separatorTableArray,
						    ConditionalSeparatorTable::SharedLocalStructure &sepSharedStructure)
{

  // Here the section receiving the message (separators) knows
  // it's using sparse join inference. So it's using SparseJoinSectionTables,
  // SparseJoinMaxCliques, etc. In particular, the section knows how to find its
  // cliques with incoming messages (formerly just the single right interface clique),
  // how many MaxCliques are in the section, and the MaxCliqueTable::SharedLocalStructure
  // of each of those cliques...

  // for each maxCliques[i] with an incoming message
  unsigned i = stss_it.prev_ri();
  //   for each of maxClique[i]'s cross-section separators j
          maxCliques[i].deReceiveFromIncommingSeparator(sourceSectionStructures.maxCliquesSharedStructure[i],
							*separatorTableArray /* [j] */, sepSharedStructure);
}


void 
SparseJoinSectionTables::printAllCliques(FILE *f, BP_Range *clique_print_range,
					 SectionIterator &stss_it, PartitionStructures &ss,
					 const bool normalize, const bool unlog,
					 const bool justPrintEntropy,
					 ObservationFile *obs_file)
{
  char buff[2048];
  BP_Range::iterator it = clique_print_range->begin();
  if (obs_file)
    obs_file->setFrame(stss_it.cur_st());
  while (!it.at_end()) {
    const unsigned cliqueNum = (unsigned)(*it);
    if (cliqueNum < ss.maxCliquesSharedStructure.size()) {
      if (obs_file) {
	maxCliques[cliqueNum].printCliqueEntries(ss.maxCliquesSharedStructure[cliqueNum],
						 obs_file,normalize, unlog);
      } else {
	sprintf(buff,"Section %d (%s), Clique %d:", stss_it.cur_st(), stss_it.cur_nm(), cliqueNum); 
	maxCliques[cliqueNum].printCliqueEntries(ss.maxCliquesSharedStructure[cliqueNum],
						 f,buff,normalize,unlog,justPrintEntropy);
      }
    }
    it++;
  }
}


/*-
 *-----------------------------------------------------------------------
 * SparseJoinSectionTables::emIncrement()
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
SparseJoinSectionTables::emIncrement(PartitionStructures& ps,
			     const logpr probE,
			     const bool localCliqueNormalization,
			     const double emTrainingBeam)
{
  for (unsigned cliqueNo=0;cliqueNo < ps.maxCliquesSharedStructure.size(); cliqueNo++ ) {
    maxCliques[cliqueNo].emIncrement(ps.maxCliquesSharedStructure[cliqueNo],
				     probE,localCliqueNormalization,emTrainingBeam);
  }
}


void SparseJoinSectionTables::init(PartitionStructures& ps)
{
  for (unsigned cliqueNo=0;cliqueNo<ps.maxCliquesSharedStructure.size();cliqueNo++)
    maxCliques[cliqueNo].init(*(ps.maxCliquesSharedStructure[cliqueNo].origin));
  for (unsigned sepNo=0;sepNo<ps.separatorCliquesSharedStructure.size();sepNo++)
    separatorCliques[sepNo].init(*(ps.separatorCliquesSharedStructure[sepNo].origin));
}


/*-
 *-----------------------------------------------------------------------
 * SparseJoinSectionTables::reportMemoryUsageTo()
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
SparseJoinSectionTables::reportMemoryUsageTo(PartitionStructures& ps,FILE *f)
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
