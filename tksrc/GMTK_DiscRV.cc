/*
 * GMTK_RV.cc
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * The top level GMTK random variable object for the RV class hierarchy.
 *
 *
 *
 */



#include "general.h"
VCID("$Header$");

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <string.h>

#include "GMTK_DiscRV.h"
#include "GMTK_MTCPT.h"


/*-
 *-----------------------------------------------------------------------
 * printSelf()
 *      prints a one-line summary of the detailed information about this RV.
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
void DiscRV::printSelf(FILE *f,bool nl)
{
  printNameFrameValue(f,false);
  fprintf(f,"discrete cardinality = %d%s",cardinality,nls(nl));
}



/*-
 *-----------------------------------------------------------------------
 * printSelfVerbose()
 *      prints a multi-line verbose description of this RV.
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
void DiscRV::printSelfVerbose(FILE *f)
{
  fprintf(f,"Discrete Random variable:\n");
  printNameFrameValue(f,true);
  fprintf(f,"From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
  fprintf(f,"RV has cardinality = %d\n",cardinality);
}


#if 0
/*-
 *-----------------------------------------------------------------------
 * identicalStructureWith.
 *      Returns true if this rv has identical structure with that of other.
 *      "identical structure" means that the r.v. have the same
 *      number, type, and cardinality parents. If this returns
 *      true, then it will be valid to tie parameters between
 *      these two random variables. Note that this routine
 *      might need to change for each subclass of this class.
 * 
 * Preconditions:
 *      Both rvs must have parents filled in.
 *
 * Postconditions:
 *      If function returns true, then the variables have
 *      identical structure, otherwise not.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      indicating boolean.
 *
 *-----------------------------------------------------------------------
 */
bool
DiscRV::identicalStructureWith(RV& other)
{
  if (!other.discrete())
    return false;
  if (other.switching())
    return false;

  DiscRV* dother = (DiscRV*)&other;

  // Note: there are no switching parents for this type of RV and its
  // default subclasses.

  // Now check the set of conditional parents.
  if (allParents.size() != dother->allParents.size())
    return false;

  for (unsigned i=0;i<allParents.size();i++) {
    if (allParents[i]->discrete() != dother->allParents[i]->discrete())
      return false;
    if (allParents[i]->discrete()) {
	if (((DiscRV*)allParents[i])->cardinality 
	    != 
	    ((DiscRV*) ((DiscRV*)&other)->allParents[i])->cardinality)
	  return false;
    }
  }
  return true;
}



/*-
 *-----------------------------------------------------------------------
 * tieParametersWith()
 *      Ties the parameters of 'this' with whatever those of 'other' are. 
 *      'other' and 'this' must be identical structuraly, if the
 *       'checkStructure' option is true.
 * 
 * Preconditions:
 *      other must be a fully instantiated RV with parameters, and 'this'
 *      and 'other' must be structurally identical (if arg is true)
 *
 * Postconditions:
 *      'this' has the identical _tied_ parameters with 'other'
 *
 * Side Effects:
 *      Changes the internal parameter data structures, but does not delete anything.
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
void
DiscRV::tieParametersWith(RV*const _other,
			  bool checkStructure)
{
  assert ( _other -> discrete() );
  DiscRV * dother = (DiscRV*) _other;

  if (checkStructure && !identicalStructureWith(*dother))
    error("ERROR: trying to tie parameters of Discrete RV '%s' with Discrete RV '%s' but they have different structure.",
	  name().c_str(),dother->name().c_str());

  curCPT = dother->curCPT;
}


#endif

/*-
 *-----------------------------------------------------------------------
 * log10ProductCardOfParentsNotContainedInSet()
 *      Returns the log_10 product card of parents not contained in given set.
 *
 * Preconditions:
 *      allParents member must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      bool
 *
 *-----------------------------------------------------------------------
 */
double DiscRV::log10ProductCardOfParentsNotContainedInSet(const set <RV*> givenSet)
{
  // first get parents not contained in set.
  set <RV*> res;
  set_difference(allParents.begin(),allParents.end(),
		 givenSet.begin(),givenSet.end(),
		 inserter(res,res.end()));

  if (res.size() == 0)
    // all parents are contained in the given set
    return 0;

  set<RV*>::iterator it;
  double weight = -1;
  for (it = res.begin(); it != res.end(); it++) {
    // if it is a parent, it must be discrete
    RV* rv = (*it);
    assert ( rv->discrete() );
    DiscRV* drv = (DiscRV*)rv;
    double log_card = ::log10(drv->cardinality);
    if (weight == -1)
      weight = log_card;
    else
      weight = log10add(weight,log_card);
  } 

  return weight;
}



/*-
 *-----------------------------------------------------------------------
 * cloneRVShell()
 *      clones a shell of the current random variable (see GMTK_RV.h for docs)
 *
 * Preconditions:
 *      RV must be instantiated and with parameters (i.e., what lives in the template RVs).
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
DiscRV* DiscRV::cloneRVShell()
{
  DiscRV*rv = (DiscRV*)RV::cloneRVShell();
  rv->cardinality = cardinality;
  rv->curCPT = curCPT;
  return rv;
}



/*-
 *-----------------------------------------------------------------------
 * computeParentsChildSatisfyingGrandChild()
 *      pass arguments down to MTCPT to count cases for VE seps.
 *      This routine is here to avoid a cpp circular dependency
 *
 * Preconditions:
 *      see caller
 *
 * Postconditions:
 *      see caller
 *
 * Side Effects:
 *      see caller
 *
 * Results:
 *      see caller
 *
 *-----------------------------------------------------------------------
 */
void DiscRV::computeParentsChildSatisfyingGrandChild(
	    // input arguments
	    unsigned par, // parent number
	    vector <RV*> & parents, 
	    vector <RV*> & hiddenParents,
	    PackCliqueValue& hiddenParentPacker,
	    sArray < DiscRVType*>& hiddenNodeValPtrs,
	    RV* child,
	    RV* grandChild,
	    // output arguments
	    sArray < unsigned >& packedParentVals,
	    unsigned& num)
{
  assert ( !switching() && deterministic() && curCPT->cptType == CPT::di_MTCPT );
  MTCPT* mtcpt = (MTCPT*) curCPT;
  return mtcpt->computeParentsChildSatisfyingGrandChild(par,parents,hiddenParents,hiddenParentPacker,
							hiddenNodeValPtrs,child,grandChild,
							packedParentVals,num);
}





/*-
 *-----------------------------------------------------------------------
 * printRVSetAndCards}()
 *      Prints out the set of random variables and their cardinalities as well when discrete.
 *
 * Preconditions:
 *      f must be open, locset a set of RVs.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      void
 *
 *-----------------------------------------------------------------------
 */
void printRVSetAndCards(FILE*f,vector<RV*>& locset,const bool nl) 
{
  bool first = true;
  for (unsigned i=0;i<locset.size();i++) {
    RV* rv = locset[i];
    if (!first)
      fprintf(f,",");
    if (rv->discrete())
      RV2DRV(rv)->printNameFrameCard(f,false);
    else 
      rv->printNameFrame(f,false);
    first = false;
  }
  if (nl) fprintf(f,"\n");
}

void printRVSetAndCards(FILE*f,set<RV*>& locset,bool nl)
{
  bool first = true;
  set<RV*>::iterator it;
  for (it = locset.begin();
       it != locset.end();it++) {
    RV* rv = (*it);
    if (!first)
      fprintf(f,",");
    if (rv->discrete())
      RV2DRV(rv)->printNameFrameCard(f,false);
    else 
      rv->printNameFrame(f,false);
    first = false;
  }
  if (nl) fprintf(f,"\n");
}

