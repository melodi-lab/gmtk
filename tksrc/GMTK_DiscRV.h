/*
 * GMTK_DiscRV.h
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
 * The discrete random variable type.
 *
 *
 *
 */

#ifndef GMTK_DISC_RV_H
#define GMTK_DISC_RV_H

#include "GMTK_DiscRVType.h"

#include <vector>

#include "GMTK_RV.h"
#include "GMTK_CPT.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_PackCliqueValue.h"

class FileParser;

class DiscRV : public RV {
  friend class FileParser;
  friend class CPT;
  friend class MDCPT;
  friend class MSCPT;
  friend class MTCPT;


  /////////////////////////////////////////
  // Support for Setting members of a RV or members of a subclass of
  // RVs. It is assumed that the appropriate thing is done for each
  // subclass (including checking for errors by FileParser).
  /////////////////////////////////////////


  // set the cpts to what is given. Sub-classes may re-define as
  // appropriate.
  virtual void setCpts(vector<CPT*> &cpts) {
    // for a non-switching RV, we only have one CPT.
    assert (cpts.size() == 1);
    curCPT = cpts[0];
  }


public:

  // The "cardinality" or "state space" of this random variable,
  // meaning the number of values that this random can possibly take
  // on (some perhaps even with zero probability). Values are always
  // in the range [0:cardinality-1] (so they are never
  // negative). Making the cardinality as small as possible will save
  // memory and run-time.
  // TODO: take this from the rv_info object.
  unsigned cardinality;

  // The current value of this random variable. Note that there is no
  // real reason to make this smaller (e.g., short, uchar) as clique
  // values are stored in packed form in GMTK.
  DiscRVType val;

protected:

  ////////////////////////////////////////////////////////////////
  // The current CPT to use for this DiscRV. It is "current" in the
  // sence that it valid for the set of parent values that are
  // currently assigned. It is assigned within the FileParser class.
  CPT* curCPT;

  // CPT iterator to use in begin()/next() iterators of this RV.
  CPT::iterator it;



public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  DiscRV(RVInfo& _rv_info,
	unsigned _timeFrame = ~0x0,
	unsigned _cardinality = 0)
    : RV(_rv_info,_timeFrame),cardinality(_cardinality)
  {
    curCPT = NULL;
  }

  virtual ~DiscRV() {;}

  // printing routines.
  virtual void printNameFrameValue(FILE *f,bool nl=true) {
    RV::printNameFrame(f,false);
    fprintf(f,"=%d%s",val,nls(nl));
  }
  virtual void printSelf(FILE *f,bool nl=true);
  virtual void printSelfVerbose(FILE *f);
  
  void printNameFrameCard(FILE *f,const bool nl=true) {
    fprintf(f,"%s(%d)[%d]%s",name().c_str(),frame(),cardinality,nls(nl));
  }

  //////////////////////////////////////////////////////////
  // computing with probabilities
  //////////////////////////////////////////////////////////

  // Compute the probability of this RV given the current values of
  // the parents, usign whatever their values are currently set to.
  // See also begin() and next() iterators below, to iterate
  // through all values of this rv given current parent values.
  inline virtual void probGivenParents(logpr& p) {
    p = curCPT->probGivenParents(allParents,this);
  }

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Parameter and Structure adjustment/checking       /////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  // Sets the parameters determining probGivenParents() to "random" values.
  virtual void makeRandom() {
    curCPT->makeRandom();
  }

  ////////////////////////////////////////////////////////////////////////
  // Sets the parameters determining probGivenParents() to "uniform" values.
  virtual void makeUniform() {
    curCPT->makeUniform();
  }

#if 0
  ////////////////////////////////////////////////////////////////////////
  // Ties the parameters of 'this' with whatever those of 'other' are. 
  // 'other' and 'this' must be identical structuraly.
  virtual void tieParametersWith(RV*const other,
				 bool checkStructure=true);

  ////////////////////////////////////////////////////////////////////////
  // Returns true of 'this' and 'other' are structurally identical.
  virtual bool identicalStructureWith(RV& other);
#endif

  //////////////////////////////////////////////////////////
  // Iterate through values of this RV                   ///
  //////////////////////////////////////////////////////////

  // Iterators with probabilities for this RV. 
  // Begin starts iteration, placing the probability in p 
  virtual void begin(logpr& p) = 0;

  // next, moves to the next value, and if it returns true, then the
  // next value is set in the RV and the probability is returned in p.
  // Otherwise, if the function returns false, p is undefined.
  virtual bool next(logpr& p) = 0;


  // If this variable is deterministic (i.e., deterministic() returns
  // true), then this routine assumes that 1) the child is hidden, 2)
  // the parents are set sufficiently so that the child can be
  // determined with probability one, and 3) sets this child
  // accordingly.
  virtual void assignDeterministicChild() { error("INTERNAL ERROR: invalidly called RV::assignDeterministicChild()"); }

  ////////////////////////////////////////////////////////////////////////
  // Sample the current distribution setting the random variable to a set
  // of values depending on the current parent values.
  virtual void randomSample() {
    curCPT->becomeAwareOfParentValues( allParents, this );
    curCPT->randomSample(this); 
  }

  /////////////////////////////////////////
  // EM Learning                         //
  /////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // Increment the statistics with probabilty 'posterior' for the
  // current random variable's parameters, for the case where
  // the random variable and its parents are set to their current
  // values (i.e., the increment corresponds to the currently set
  // parent/child values). 
  inline virtual void emIncrement(logpr posterior) {
    curCPT->emIncrement(posterior,allParents,this);
  }

  /////////////////////////////////////////
  // Support for unrolling               //
  /////////////////////////////////////////

  /////////////////////////////////////////
  // Support for Triangulation           //
  /////////////////////////////////////////

  // Return the log10 of the product of the cardinalities of all
  // parents of the current random variable, but not in the given Set.
  double log10ProductCardOfParentsNotContainedInSet(const set <RV*> givenSet);

  // returns true if the RV is determinstic.
  bool deterministic() {  return rv_info.deterministic();  }

  // returns true if all cpts are sparse or determinstic.
  bool sparse() { return rv_info.sparse();  }

  // This routine returns the average cardinality (average number of
  // possible child random variable values with non-zero probability)
  // for the current random variable.  This routines assume that the
  // RV's parents are in the same clique (or an earlier clique
  // relative to the JT root) so that we only consider the child
  // values with non-zero probability.
  virtual unsigned averageCardinality() {
    assert ( rv_info.discImplementations.size() == 1 );
    if (curCPT != NULL)
      return curCPT->averageCardinality();
    // still here, we need to use rv_info.
    // TODO: make this a member of RVInfo.
    if (rv_info.discImplementations[0] == CPT::di_MTCPT)
      return 1;
    else 
      return cardinality;
  }

  // This routine returns the maximum possible cardinality (max number
  // of possible child random variable values with non-zero
  // probability) for the current random variable.  This routines
  // assume that the RV's parents are in the same clique (or an
  // earlier clique relative to the JT root) so that we only consider
  // the child values with non-zero probability.
  virtual unsigned maxCardinality() {
    assert ( rv_info.discImplementations.size() == 1 );
    if (curCPT != NULL)
      return curCPT->maxCardinality();
    // still here, we need to use rv_info.
    // TODO: make this a member of RVInfo.
    if (rv_info.discImplementations[0] == CPT::di_MTCPT)
      return 1;
    else 
      return cardinality;
  }

  bool iterable() const {
    if (curCPT == NULL)
      return false;
    else {
      return curCPT->iterable();
    }
  }


  // various other routines want a cardinality to 'use'.  they use the
  // routine useCardinality() which can be defined either as the avg
  // or max cardinality above.
  unsigned useCardinality() { return maxCardinality(); }  


  // only valid when this var is non-switching, observed, with a deterministic implementation.
  void computeParentsChildSatisfyingGrandChild(
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
	    unsigned& num);
    // body placed in .cc file to avoid circularity;


  /////////////////////////////////////////
  // Cloning Support                     //
  /////////////////////////////////////////

  virtual DiscRV* cloneRVShell();

};

// Define a few down-cast routine here for when we have a RV that is
// known to be a DiscRV. Note that this is not entirely type-safe, but
// speed is more important than type safety here, so we need to be
// able to do static down casts.
inline DiscRV* RV2DRV(RV* rv) { 
  return ((DiscRV*)rv);
}
inline DiscRV& RV2DRV(RV& rv) { 
  return *((DiscRV*)&rv);
}
inline const DiscRV* const RV2DRV(const RV* const rv) { 
  return ((const DiscRV* const)rv);
}

void printRVSetAndCards(FILE*f,set<RV*>& locvec,const bool nl=true);
void printRVSetAndCards(FILE*f,vector<RV*>& locvec,const bool nl=true);

// TODO: make a special subclass DummyDiscRV where RVs are easily
// creatable on the fly for certain test code (CPTs, DTs, etc).

#endif
