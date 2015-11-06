/*-
 * GMTK_MaxClique.cc
 *
 *     maxClique support and JT probabilistic inference. Includes the
 *     implementation of a message in the message passing algorithm.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */




/*
 * TODO: turn this into multiple files
 *   mc, mctable CE, mctable DE, mctable prune, csctable
 *
 */



/*
 * TODO:
 *   1) Maybe: In some of the structures the choice of when
 *   to use a hash table vs. an inlined value is done
 *   by either using an unsigned* or a unsigned. Create
 *   options where we do 'unsigned val[LEN]' and where
 *   only if packed value is > LEN words do we resort to
 *   the hash table.
 *
 *   2) Figure out way to remove the (unsigned**) type casts
 *   in front of call to packer object (see below).   
 *
 */

#ifdef HAVE_CONFIG_H
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
#include "psp.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MaxClique.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_JunctionTree.h"

VCID(HGID)


// for sorting an array of RVs ascending based on increasing cardinality
struct ParentCardinalityCompare 
{  
  bool operator() (const RV* rv1,
		   const RV* rv2) 
  {
    // place observed ones first.
    if (RV2DRV(rv1)->discreteObservedImmediate())
      return false;
    else if (RV2DRV(rv2)->discreteObservedImmediate())
      return true;
    else return (RV2DRV(rv1)->cardinality < RV2DRV(rv2)->cardinality);
  }
};



///////////////////////////////////////////////
// VE separator files information.
///////////////////////////////////////////////
bool SeparatorClique::recomputeVESeparatorTables = false;
const char* SeparatorClique::veSeparatorFileName = "veSeparatorFile.dat";
bool SeparatorClique::generatingVESeparatorTables = "true";
FILE* SeparatorClique::veSeparatorFile = NULL;
float SeparatorClique::veSeparatorLogProdCardLimit = 7.0; // i.e., 1e7=10M is default max.


/*
 * Memory management options
 */

// Set these to reasonable defaults (see https://j.ee.washington.edu/trac/gmtk/ticket/384)
// These defaults may be changed by the -memoryGrowth command line argument (see
// https://j.ee.washington.edu/trac/gmtk/ticket/342)
unsigned SeparatorClique::aiStartingSize = 1;
float    SeparatorClique::aiGrowthFactor = 1.05;

unsigned SeparatorClique::remStartingSize = 1;
float    SeparatorClique::remGrowthFactor = 1.05;

unsigned SeparatorClique::sepSpaceMgrStartingSize = 1;
float    SeparatorClique::sepSpaceMgrGrowthRate   = 1.05;
float    SeparatorClique::sepSpaceMgrDecayRate    = 0.0;

unsigned SeparatorClique::remSpaceMgrStartingSize = 1;
float    SeparatorClique::remSpaceMgrGrowthRate   = 1.05;
float    SeparatorClique::remSpaceMgrDecayRate    = 0.0;


/*
 *
 * separator beam width, for separator-based beam pruning.  Default value is
 * very large (1.0/0.0 = est. of infty) meaning that we do no beam
 * pruning.
 *
 */
double
SeparatorClique::separatorBeam=(-LZERO);



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        SeparatorClique support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

SeparatorClique::SeparatorClique(MaxClique& c1, MaxClique& c2)
  :  veSeparator(false),
     separatorValueSpaceManager(1,     // starting size
				sepSpaceMgrGrowthRate,   // growth rate
				1,     // growth addition
				sepSpaceMgrDecayRate),   // decay rate 
     remainderValueSpaceManager(1,     // starting size
				remSpaceMgrGrowthRate,   // growth rate
				1,     // growth addition
				remSpaceMgrDecayRate)    // decay rate
     
{
  nodes.clear();

  // not a ve sep clique.
  veSepClique = NULL;

  skipMe = false;

  // create a set of nodes that is the intersection of the two
  set_intersection(c1.nodes.begin(),c1.nodes.end(),
		   c2.nodes.begin(),c2.nodes.end(),
		   inserter(nodes,nodes.end()));
  // assert (nodes.size() > 0);
  
}

SeparatorClique::~SeparatorClique()
{
  if (veSeparator) {
    delete veSepClique;
  }
}


#if 0
// This constructor is not used for now.
SeparatorClique::SeparatorClique(SeparatorClique& from_sep,
				 vector <RV*>& newRvs,
				 map < RVInfo::rvParent, unsigned >& ppf,
				 const unsigned int frameDelta)
  :  separatorValueSpaceManager(1,     // starting size
				2.0,   // growth rate
				1,     // growth addition
				0.90), // decay rate 
     remainderValueSpaceManager(1,     // starting size
				2.0,   // growth rate
				1,     // growth addition
				0.90)  // decay rate 
{
  set<RV*>::iterator it;

  skipMe = false;

  // clone over nodes RVs and accumulated intersection.
  for (it = from_sep.nodes.begin();
       it != from_sep.nodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }

  // copy over accumulated intersection
  for (it=from_sep.accumulatedIntersection.begin();
       it != from_sep.accumulatedIntersection.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    


    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    accumulatedIntersection.insert(nrv);
  }

  // and 'remainder'
  for (it=from_sep.remainder.begin();
       it != from_sep.remainder.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    remainder.insert(nrv);
  }


}
#endif








/*-
 *-----------------------------------------------------------------------
 * SeparatorClique::clearSeparatorValueCache()
 *   
 *   Clear out hash table and other memories for use between segments.
 *
 * Preconditions:
 *   Data structures must have been set up.
 *
 * Postconditions:
 *   All shared hash tables and value holders for this clique are cleared if 'force'
 *   is true.
 *
 * Side Effects:
 *     internal data structures will change.
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
SeparatorClique::clearSeparatorValueCache(bool force) 
{
  if (force && accPacker.packedLen() > ISC_NWWOH_AI) {
    accValueHolder.prepare();
    accSepValHashSet.clear(aiStartingSize);
  }
  if (force && remPacker.packedLen() > ISC_NWWOH_RM) { 
    remValueHolder.prepare();
    remSepValHashSet.clear(remStartingSize);
  }
  // shrink space asked for by clique values. 
  separatorValueSpaceManager.decay();
  remainderValueSpaceManager.decay(); 
}



/*-
 *-----------------------------------------------------------------------
 * SeparatorClique::prepareForUnrolling()
 *   
 *   Sets a few last data structures that are needed for both unrolling
 *   and inference to work. 
 *
 * Preconditions:
 *   Set of nodes must be valid meaning that it has valid neighbors,
 *   parents, and children member variables.
 *
 * Postconditions:
 *   computed weight is provided.
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
SeparatorClique::prepareForUnrolling()
{

  // set up number of hidden ndoes within accumulated intersection.

  // create a vector form of the variables.
  set<RV*>::iterator it;
  for (it = accumulatedIntersection.begin();
       it != accumulatedIntersection.end();
       it++) {
    RV* rv = (*it);
    // VE separators have a different notion of hidden (see below).

    if ((!veSeparator && rv->hidden())
	||
	(veSeparator && !rv->discreteObservedImmediate()))
      hAccumulatedIntersection.push_back(rv);
  }

  if (hAccumulatedIntersection.size() > 0) {
    // we only use these structures if there is an intersection.
    new (&accPacker) PackCliqueValue(hAccumulatedIntersection);
    assert( accPacker.packedLen() > 0);
    if (accPacker.packedLen() > ISC_NWWOH_AI) {
      // only setup hash table if the packed accumulated insersection
      // set is larger than one machine word (unsigned).
      new (&accValueHolder) CliqueValueHolder(accPacker.packedLen());

      // TODO: optimize starting size.
      new (&accSepValHashSet) vhash_set< unsigned > (accPacker.packedLen(),aiStartingSize);
    }
  }

  // create a vector form of the variables.
  for (it = remainder.begin();
       it != remainder.end();
       it++) {
    RV* rv = (*it);

    if ((!veSeparator && rv->hidden())
	||
	(veSeparator && !rv->discreteObservedImmediate()))
      hRemainder.push_back(rv);
  }

  if (hRemainder.size() > 0) {
    new (&remPacker) PackCliqueValue(hRemainder);
    assert (remPacker.packedLen() > 0);
    if (remPacker.packedLen() > ISC_NWWOH_RM) { 
      // Only setup hash table if the packed remainder set is larger
      // than one machine word (unsigned).
      new (&remValueHolder) CliqueValueHolder(remPacker.packedLen());
      new (&remSepValHashSet) vhash_set< unsigned > (remPacker.packedLen(),remStartingSize);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // VE Separator support
  ////////////////////////////////////////////////////////////////////////

  // create any VE separators
  if (veSeparator) {

    ConditionalSeparatorTable::SharedLocalStructure veSepCliqueSharedStructure;

    // Need also to create the shared ConditionalSeparatorTable that
    // will be used for all instances of this SeparatorClique, this is
    // done below.

    // sanity check.
    assert ( veSepInfo.child != NULL && veSepInfo.parents.size() == veSepInfo.child->allParents.size() );

    // sort parents by increasing cardinality, for better branch
    // prediction.
    // printf("before sort");printRVSetAndCards(stdout,veSepInfo.parents);
    sort(veSepInfo.parents.begin(),veSepInfo.parents.end(),ParentCardinalityCompare());
    // printf("after sort");printRVSetAndCards(stdout,veSepInfo.parents);

    if (veSepInfo.grandChild == NULL) {
      ///////////////////////////
      // the PC case.

      ObsDiscRV* odc = (ObsDiscRV*)veSepInfo.child;

      // set any observed variables to their observed values
      // the child is guaranteed to be immediate const observed
      odc->setToObservedValue();
      // set any observed parents to their observed values.
      for (unsigned i=0;i<odc->allParents.size();i++) {
	// only discrete immediate observed RV parents are not stored,
	// the rest might vary per frame so we need to compute
	// the table for them with all possible values [0,card-1].
	if (odc->allParents[i]->discreteObservedImmediate()) {
	  ObsDiscRV*odrv = (ObsDiscRV*)odc->allParents[i];
	  odrv->setToObservedValue();
	}
      }

      // create a vector of just the hidden parents
      vector<RV*> hiddenParents;
      for (unsigned i = 0; i < odc->allParents.size() ; i++ ) {
	RV* rv = odc->allParents[i];
	// code will need to change when using continuous hidden
	// variables. Only consider a variable hidden if it is not a
	// discrete immediate observed, since any such varible might
	// still have varying values over time.
	if (!rv->discreteObservedImmediate())
	  hiddenParents.push_back(rv);
      }

      assert ( hiddenParents.size() == hAccumulatedIntersection.size() + hRemainder.size() );

      // Generate the table of parent values satisfying the observed
      // child, packing them into an array.
      PackCliqueValue parentPacker(hiddenParents);
      
      // create an array of direct pointers to the discrete hidden RV
      // values within the hidden nodes just computed above.
      sArray < DiscRVType*> hiddenNodeValPtrs;
      hiddenNodeValPtrs.resize(hiddenParents.size());
      for (unsigned i=0;i<hiddenParents.size();i++) {
	DiscRV* drv = RV2DRV(hiddenParents[i]);
	hiddenNodeValPtrs[i] = &(drv->val);
      }

      // create an array where we place the packed parent values that
      // satisfy the child.
      sArray < unsigned > packedParentVals(parentPacker.packedLen());

      unsigned num = 0;

      if (generatingVESeparatorTables == true) {

	if (message(Inference, Low)) {
	  float logProdCard = 
	    log10((double)RV2DRV(odc->allParents[0])->cardinality);
	  for (unsigned i=1;i<odc->allParents.size();i++) {
	    if (!odc->allParents[i]->discreteObservedImmediate())
	      logProdCard += log10((double)RV2DRV(odc->allParents[i])->cardinality);
	  }
	  printf("VE separator PC computation: Iterating 10^%f possible parent vals.\n",logProdCard);
	}
	odc->computeParentsSatisfyingChild(0,
					   veSepInfo.parents,
					   hiddenParents,
					   parentPacker,
					   hiddenNodeValPtrs,
					   odc,
					   packedParentVals,
					   num);
	if (num == 0) {
	  error("ERROR: found VE separator with 0 entries. No way for any parents to satisfy %s(%d)=const with non-zero probability. Makes all observations have probability 0.\n",
		odc->name().c_str(),odc->frame());
	}

	if (veSeparatorFile != NULL) {
	  // then we need to save this information to a file for next time.
	  unsigned tmp;
	  bool writeError = false;

	  // write a bunch of ID information.
	  // write a zero for case PC
	  tmp = 0;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the cardinality of the child.
	  tmp = odc->cardinality;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the number and cardinalities of the parents
	  tmp = odc->allParents.size();
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  for (unsigned i=0; i < odc->allParents.size(); i++) {
	    tmp = RV2DRV(odc->allParents[i])->cardinality;
	    if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	      writeError = true;
	  }

	  // write the number of elements 
	  if (!fwrite(&num,sizeof(num),1,veSeparatorFile))
	    writeError = true;

	  // write the size of the elements
	  tmp = parentPacker.packedLen();
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  
	  // finally write the set of parent values
	  if (!fwrite(packedParentVals.ptr,parentPacker.packedLen()*sizeof(unsigned),num,veSeparatorFile))
	    writeError = true;
	  
	  if (writeError)
	    error("ERROR: writing to PC VE separator file (%s)\n",veSeparatorFileName);


	}
      } else {
	// we must have a file to read from here.
	assert (veSeparatorFile != NULL);
	unsigned tmp;
	unsigned corrupt = 0;

	// read in and check the ID information for this separator information.
	if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	  corrupt = 1;
	if (!corrupt && tmp != 0)
	  corrupt = 2;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 3;
	}
	if (!corrupt && tmp != odc->cardinality)
	  corrupt = 4;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 5;
	}
	if (!corrupt && tmp != odc->allParents.size())
	  corrupt = 6;
	if (!corrupt) {
	  for (unsigned i=0; i < odc->allParents.size(); i++) {
	    if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1) {
	      corrupt = 7;
	      break;
	    }
	    if (!corrupt && tmp != RV2DRV(odc->allParents[i])->cardinality) {
	      corrupt = 8;
	      break;
	    }
	  }
	}
	if (!corrupt) {
	  if (fread(&num,sizeof(num),1,veSeparatorFile) != 1)
	    corrupt = 9; 
	  // we can't check num.
	}
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 10;
	}
	if (!corrupt && tmp != parentPacker.packedLen())
	  corrupt = 11;
	packedParentVals.resize(num*parentPacker.packedLen());
	if (!corrupt) {
	  if (fread(packedParentVals.ptr,parentPacker.packedLen()*sizeof(unsigned),num,veSeparatorFile) != num)
	    corrupt = 12;
	}

	if (corrupt)
	  error("ERROR: corrupt/wrong PC VE separator file (%s) with respect to current structure/triangulation/command options. Reason %d.\n",veSeparatorFileName,corrupt);
      }

      infoMsg(IM::Inference, Low,"VE separator PC generation: %d parent vals satisfying this case.\n",num);

      // Now create an ConditionalSeparatorTable and insert all
      // of the packed parent values we generated or got above.

      veSepClique = new ConditionalSeparatorTable(*this,
						  veSepCliqueSharedStructure);

#ifdef CHECK_VE_SEP_UNIQUE
      set < vector < unsigned > > tmpset;
#endif      

      for (unsigned i=0;i<num;i++) {
	parentPacker.unpack(&(packedParentVals.ptr[i*parentPacker.packedLen()]),
			    hiddenNodeValPtrs.ptr);
#ifdef CHECK_VE_SEP_UNIQUE
	{
	  vector < unsigned > foo(hiddenNodeValPtrs.size());
	  for (unsigned i=0;i<foo.size();i++) {
	    foo[i] = *(hiddenNodeValPtrs[i]);
	  }
	  if (tmpset.find(foo) != tmpset.end())
	    printf("Entry alreay in set\n");
	  else
	    tmpset.insert(foo);
	}
#endif
	if (message(Inference, Max+5)) {
	  printf("Inserting into PC VE separator:"); printRVSetAndValues(stdout,hiddenParents);
	}
	// insert current RV values into the separator. 
	veSepClique->insert(*this,
			    veSepCliqueSharedStructure);
      }

    } else {
      //////////////////////////////////////////////
      // the PCG case, grand-child is not NULL.

      // We know grandChild is observed.

      // set any observed parents to their observed values
      // Child and/or an parents of child might be observed. Need to check
      // those cases.
      // observed discrete grand child
      ObsDiscRV* odgc = (ObsDiscRV*)veSepInfo.grandChild;
      // dicrete child
      DiscRV* dc = (DiscRV*)veSepInfo.child;

      // the grandchild is guaranteed to be immediate const observed
      odgc->setToObservedValue();

      // child is going to be hidden, otherwise not PCG case.
      // if (dc->discreteObservedImmediate()) {
      //   ObsDiscRV* odc = (ObsDiscRV*)dc;
      //   odc->setToObservedValue();
      // }

      // set any observed parents to their observed values.
      for (unsigned i=0;i<dc->allParents.size();i++) {
	// only discrete immediate observed RV parents are not stored,
	// the rest might vary per frame so we need to compute
	// the table for them with all possible values [0,card-1].
	if (dc->allParents[i]->discreteObservedImmediate()) {
	  ObsDiscRV*odrv = (ObsDiscRV*)dc->allParents[i];
	  odrv->setToObservedValue();
	}
      }

      // create a vector of just the hidden (parents and child)
      vector<RV*> hiddenParents;
      for (unsigned i = 0; i < dc->allParents.size() ; i++ ) {
	RV* rv = dc->allParents[i];
	// code will need to change when using continuous hidden
	// variables.  only consider a variable hidden if it is not a
	// discrete immediate observed, since any such varible might
	// still have varying values over time.
	if (!rv->discreteObservedImmediate())
	  hiddenParents.push_back(rv);
      }
      // guaranteed not to be observed
      assert ( dc->hidden() );
      hiddenParents.push_back(dc);

      assert ( hiddenParents.size() == hAccumulatedIntersection.size() + hRemainder.size() );

      // Generate the table of parent values and child satisfying the
      // observed grandchild, packing them into an array.
      PackCliqueValue parentPacker(hiddenParents);
      
      // create an array of direct pointers to the discrete hidden
      // RV values within the hidden nodes just computed above.n 
      sArray < DiscRVType*> hiddenNodeValPtrs;
      hiddenNodeValPtrs.resize(hiddenParents.size());
      for (unsigned i=0;i<hiddenParents.size();i++) {
	DiscRV* drv = RV2DRV(hiddenParents[i]);
	hiddenNodeValPtrs[i] = &(drv->val);
      }

      // create an array where we place the packed parent values that
      // satisfy the child.
      sArray < unsigned > packedParentVals(parentPacker.packedLen());

      unsigned num = 0;

      if (generatingVESeparatorTables == true) {
	if (message(Inference, Low)) {
	  float logProdCard = 
	    log10((double)RV2DRV(dc->allParents[0])->cardinality);
	  for (unsigned i=1;i<dc->allParents.size();i++) {
	    if (!dc->allParents[i]->discreteObservedImmediate())
	      logProdCard += log10((double)RV2DRV(dc->allParents[i])->cardinality);
	  }
	  printf("VE separator PCG computation: Iterating 10^%f possible parent vals.\n",logProdCard);
	}

	dc->computeParentsChildSatisfyingGrandChild(0,
						    veSepInfo.parents,
						    hiddenParents,
						    parentPacker,
						    hiddenNodeValPtrs,
						    dc,
						    odgc,
						    packedParentVals,
						    num);

	if (num == 0) {
	  error("ERROR: found VE separator with 0 entries. No way for any parents to satisfy %s(%d)=const with non-zero probability. Makes all observations have probability 0.\n",
		odgc->name().c_str(),odgc->frame());
	}

	if (veSeparatorFile != NULL) {
	  // then we need to save this information to a file for next time.
	  unsigned tmp;
	  bool writeError = false;

	  // write a bunch of ID information.
	  // write a zero for case PC
	  tmp = 1;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the cardinality of the grandchild.
	  tmp = odgc->cardinality;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the cardinality of the child.
	  tmp = dc->cardinality;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the number and cardinalities of the parents
	  tmp = dc->allParents.size();
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  for (unsigned i=0; i < dc->allParents.size(); i++) {
	    tmp = RV2DRV(dc->allParents[i])->cardinality;
	    if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	      writeError = true;
	  }

	  // write the number of elements 
	  if (!fwrite(&num,sizeof(num),1,veSeparatorFile))
	    writeError = true;

	  // write the size of the elements
	  tmp = parentPacker.packedLen();
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  
	  // finally write the set of parent values
	  if (!fwrite(packedParentVals.ptr,parentPacker.packedLen()*sizeof(unsigned),num,veSeparatorFile))
	    writeError = true;
	  
	  if (writeError)
	    error("ERROR: writing to PCG VE separator file (%s)\n",veSeparatorFileName);


	}
      } else {
	// we must have a file to read from here.
	assert (veSeparatorFile != NULL);
	unsigned tmp;
	unsigned corrupt = 0;

	// read in and check the ID information for this separator information.
	if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	  corrupt = 1;
	if (!corrupt && tmp != 1)
	  corrupt = 2;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 3;
	}
	if (!corrupt && tmp != odgc->cardinality)
	  corrupt = 4;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 5;
	}
	if (!corrupt && tmp != dc->cardinality)
	  corrupt = 6;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 7;
	}
	if (!corrupt && tmp != dc->allParents.size())
	  corrupt = 8;
	if (!corrupt) {
	  for (unsigned i=0; i < dc->allParents.size(); i++) {
	    if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1) {
	      corrupt = 9;
	      break;
	    }
	    if (!corrupt && tmp != RV2DRV(dc->allParents[i])->cardinality) {
	      corrupt = 10;
	      break;
	    }
	  }
	}
	if (!corrupt) {
	  if (fread(&num,sizeof(num),1,veSeparatorFile) != 1)
	    corrupt = 11; 
	  // we can't check num.
	}
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 12;
	}
	if (!corrupt && tmp != parentPacker.packedLen())
	  corrupt = 13;
	packedParentVals.resize(num*parentPacker.packedLen());
	if (!corrupt) {
	  if (fread(packedParentVals.ptr,parentPacker.packedLen()*sizeof(unsigned),num,veSeparatorFile) != num)
	    corrupt = 14;
	}

	if (corrupt)
	  error("ERROR: corrupt/wrong PCG VE separator file (%s) with respect to current structure/triangulation/command options. Reason %d.\n",veSeparatorFileName,corrupt);

      }

      infoMsg(IM::Inference, Low,"VE separator PCG generation: %d (parent,child) combinations satisfying this case.\n",num);

      // now create an ConditionalSeparatorTable and insert all
      // of the packed parent values we generated or got above.

      veSepClique = new ConditionalSeparatorTable(*this,veSepCliqueSharedStructure);
      for (unsigned i=0;i<num;i++) {
	parentPacker.unpack(&(packedParentVals.ptr[i*parentPacker.packedLen()]),
			    hiddenNodeValPtrs.ptr);
	if (message(Inference, Max+5)) {
	  printf("Inserting into PCG VE separator:"); printRVSetAndValues(stdout,hiddenParents);
	}
	// insert current RV values into the separator. 
	// TODO: make this return bool to indicate if found or not.
	veSepClique->insert(*this,
			    veSepCliqueSharedStructure);
      }
    }

  }

}






/*-
 *-----------------------------------------------------------------------
 * SeparatorClique::printAllJTInfo()
 *   
 *   prints everything JT (not not inference) about this clique to file.
 *
 * Preconditions:
 *   all variables must have been set up. prepareForUnrolling() must have
 *   been called.
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
SeparatorClique::printAllJTInfo(FILE*f)
{
  // TODO: also print out nubmer of bits for acc and rem.
  fprintf(f,"%sSeparator information: %d acc packed bits (%d words, %d splits), %d rem packed bits (%d words, %d splits)\n",
	  (veSeparator?"VE ":""),
	  accPacker.packedLenBits(),accPacker.packedLen(),accPacker.numSplits(),
	  remPacker.packedLenBits(),remPacker.packedLen(),remPacker.numSplits());

  fprintf(f,"%ld Nodes: ",(unsigned long)nodes.size()); printRVSetAndCards(f,nodes);
  fprintf(f,"%ld Acc Inter: ",(unsigned long)accumulatedIntersection.size()); printRVSet(f,accumulatedIntersection);  
  fprintf(f,"%ld Hid Acc Inter: ",(unsigned long)hAccumulatedIntersection.size()); printRVSet(f,hAccumulatedIntersection);  
  fprintf(f,"%ld remainder: ",(unsigned long)remainder.size()); printRVSet(f,remainder);  
  fprintf(f,"%ld hRemainder: ",(unsigned long)hRemainder.size()); printRVSet(f,hRemainder);

  if (veSeparator) {
    fprintf(f,"%ld VE sep parents: ",(unsigned long)veSepInfo.parents.size()); printRVSet(f,veSepInfo.parents);
    fprintf(f,"VE sep child: "); veSepInfo.child->printNameFrame(f);
    if (veSepInfo.grandChild != NULL) {
      fprintf(f,"VE sep grand child: "); veSepInfo.grandChild->printNameFrame(f);
    }
  }

}



/*-
 *-----------------------------------------------------------------------
 * SeparatorClique::reportMemoryUsageTo()
 *
 *    Report current memory usage of this inference sep clique *and*
 *    the origin sep clique to the file in units of MBs.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      results printed.
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void
SeparatorClique::
reportMemoryUsageTo(FILE *f)
{
  // Memory: OSC=Origin Separator Clique, AI = accumulated intersection, RM = remainder
  if (accPacker.packedLen() > ISC_NWWOH_AI) {
    fprintf(f,"OSC(AIVH=%luMB,AIHS=%luMB,",
	    1+accValueHolder.bytesRequested()/(1024*1024),
	    1+accSepValHashSet.bytesRequested()/(1024*1024));
  } else {
    fprintf(f,"OSC(AIVH=0MB,AIHS=0MB,");
  }
  if (remPacker.packedLen() > ISC_NWWOH_RM) { 
    fprintf(f,"RMVH=%luMB,RMHS=%luMB)",
	    1+remValueHolder.bytesRequested()/(1024*1024),
	    1+remSepValHashSet.bytesRequested()/(1024*1024));
  } else {
    fprintf(f,"RMVH=0MB,RMHS=0MB)");
  }
}




/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
