/*
 * GMTK_RV.h
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
 */

/***********************************************************************
 ***********************************************************************
 
 This is the top level GMTK random variable object for the GMTK RV
 class hierarchy.

 Below is the complete RV hierarchy. As can be seen, it is fairly
 deep, and many of the classes have only small differences in a few of
 their virtual functions. The class hierarchy was designed with the
 following goals in mind. C++ virtual functions are being used in two
 places, first to dispatch from an abstract RV to an appropriate
 routine at runtime, and then within a RV to dispatch to members of an
 appropriate CPT. Since a commitment is made to use virtual functions
 (and we are taking a speed hit for this decision), the RV hierarchy
 below is designed so that once a virtual dispatch has been made and
 we are at the callee location, a basic block of code (i.e., no
 conditionals) should be found whenever possible. Basically, all of
 the decisions about what should be done at each RV is done at
 STR-file parse time, which sets up the appropriate RV values thereby
 making these decisions. An alternative approach would have been to
 use a much smaller hierarchy and use conditionals at the callee side
 of virtual dispatch, but this will hurt branch-prediction hardware
 logic leading to much slower performance in some cases (since many of
 these dispatches lie within the inner most loops of inference).

 Regarding the hierarhcy below, class names preceded with a '-' are
 abstract and are never instantiated, class names preceded with a
 '+' may be instantiated.


 Switching Support Hierarchy (used with multiple inheritance below):

    - SwRV
      - SwDiscRV
      - SwContRV

Other classes to support implementation of scale, penalty, and shift
(also used with multiple inheritance):

    - ScPnShRV

Main RV hierarchy:

- RV
   - DiscRV
        + HidDiscRV
           + ScPnSh_HidDiscRV (also inherits from ScPnShRV)
           + Sw_HidDiscRV (also inherits from SwDiscRV)
              + ScPnSh_Sw_HidDiscRV (also inherits from ScPnShRV)
        + ObsDiscRV
           + ScPnSh_ObsDiscRV (also inherits from ScPnShRV)
           + Sw_ObsDiscRV  (also inherits from SwDiscRV)
              + ScPnSh_Sw_ObsDiscRV  (also inherits from ScPnShRV)
   - ContRV
        - // HidContRV (doesn't exist yet)
        + ObsContRV
           + ScPnSh_ObsContRV (also inherits from ScPnShRV)
           + Sw_ObsContRV  (also inherits from SwContRV)
              + ScPnSh_Sw_ObsContRV  (also inherits from ScPnShRV)

Please finish this if anyone cares to.

                               RV______________
                             /                 \
                       ContRV                 DiscRV
                     /    \
                          |
                         ObsContRV
                          /     \
              Sw_ObsContRV       ScPnSh_ObsContRV
                    /
        ScPnSh_Sw_ObsContRV  


 

Each routine has:

  create() - 
    Creates a variable of the same dynamic type as the object, even if called
    via abstract base class. Only needs to be defined for subclass that
    1) are not abstract and 2) that add member variables relatve to their
    parents. Current, create defined in:
       ObsDiscRV, Sw_ObsDiscRv.
 
  cloneRVShell()
    Creates a new random variable of the same dynamic type of the object. 
    Returned object has everything except for parents/children filled in (so
    this will have parameters). Only needs to be defined in each RV subclass that 
    adds a new member variable.

   (remove tieParametersWith, make new routine, setParametersFrom)

***********************************************************************
***********************************************************************/


#ifndef GMTK_RV_H
#define GMTK_RV_H

#include <vector>
#include <string>
#include <set>

#include "logp.h"
#include "GMTK_RVInfo.h"
#include "GMTK_NamedObject.h"

// There are many types of RVs. This structure encodes the various
// options in a current RV. Note that each different set of objects
// might or might not correspond to a different sub-class. In general,
// a separate sub-class is created if it otherwise would require a
// branch to decide what to do even after virutal dispatch to a
// routine that lives within the inner loops of inference.
struct RVType {
  unsigned discrete:1;    // discrete == 1 or continuous 
  unsigned hidden:1;      // hidden == 1, or continuous 
  unsigned switching:1;   // switching == 1, or no switching (one set of parents)
  // Probability modification: given a probability p, we can optionally
  // change the value using the formula: penalty*p^scale+offset.
  unsigned scale:1;       // scale == 1, or no scale (scale is prob. exponent)
  unsigned penalty:1;     // penalty == 1, or no current penalty (penalty is prob multiplication)
  unsigned shift:1;       // shift == 1, or no offset (offset is additive value in prob. domain)
};


class FileParser;
class RngDecisionTree;

class RV  {
  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;
  friend class JunctionTree;
  friend class MaxClique;
  friend class InferenceMaxClique;
  friend class BoundaryTriangulate;
  friend class SwRV;
  friend class SeparatorClique;
  friend class RngDecisionTree;


  /////////////////////////////////////////
  // Support for Setting members of a RV or members of a subclass of
  // RVs. It is assumed that the appropriate thing is done for each
  // subclass (including checking for errors by FileParser).
  /////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // Set the parents of this random variable to the parents given by
  // the following arguments. This routine allows to set both
  // switching parents and sets of conditonal parents. Of course, for
  // many RV types, there will be no switching parents and only one
  // set of conditional parents, and subclasses should define this
  // routine to accordingly signal an error if called with the
  // inappropriate parent configuration. Also, in the current GMTK,
  // all parents are assumed discrete, but they are passed as RVs for
  // simplicity. Lastly, this routine will also add 'this' as a child
  // to all given parents.
  virtual void setParents(vector<RV *> &sparents,vector<vector<RV *> > &cpl);

  // set the dtMapper, if this is a switching RV. Subclasses should
  // re-define accordingly.
  virtual void setDTMapper(RngDecisionTree *) { assert( 0 );  }

  virtual vector< RV* >& switchingParentsVec() {
    // a general RV has no switching parents. Subclasses must re-define as desired.
    assert ( 0 );
    // return something to satisfy compiler, but this will/should never happen.
    return allParents;
  }
  
  // get the parents (for non-switching) and conditional parents
  // for switching.
  virtual vector< RV* >& condParentsVec(unsigned j) {
    assert ( j == 0 );
    return allParents;
  }


protected:

  /////////////////////////////////////////////////////////////////////////
  // Pointer back to original R.V. file parser information structure
  // where much common info about this RV is stored.
  RVInfo& rv_info;

  /////////////////////////////////////////////////////////////////////////
  // The time frame (time slice) of the random variable.
  // Counting starts from 0.
  unsigned timeFrame;

  //////////////////////////////////////////////////////////
  // Support for the Graph of Random Variables
  //////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // The set of neighbors, for undirected graph formulation. This
  // consists of the union of the parents and the children. Note that
  // this variable is set only during triangulation. During inference,
  // this variable is likley to be empty, and the undirected graph is
  // defined only implicitly via the maxCliques.
  set<RV *> neighbors;

  ////////////////////////////////////////////////////////////////////////
  // all parents consists of all parents (including conditional and
  // switching parents).
  vector< RV* > allParents;

  ////////////////////////////////////////////////////////////////////////
  // all children are all children of this RV (via switching or
  // conditional parent of those children).
  vector< RV* > allChildren;

public:

  // return the name
  const string& name() const { return rv_info.name; }

  // return the current time frame
  unsigned frame() const { return timeFrame; }

  // return various aspects of the RV, based on the RVInfo of this RV.
  bool discrete() const { return (rv_info.rvType == RVInfo::t_discrete); }
  bool continuous() const { return !discrete(); }
  bool hidden() const { return (rv_info.rvDisp == RVInfo::d_hidden); }
  bool observed() const { return !hidden(); }
  bool discreteObservedImmediate() const {
    // returns true if the variable is discrete, observed, and the
    // observed value is a fixed constant (i.e., not variable like a
    // frame number, etc.) in the .str file.
    return (observed() && discrete() && 
	    (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_FirstIsValue));
  }
  bool switching() const { return (rv_info.switchingParents.size() > 0); } 
  // Right now, we don't distinguish between scale/penalty/shift
  // (i.e., either they're all on or all off from the class
  // hierarchy's point of view).  If we add subclasses for each
  // combination of cases for each switching condition, this might
  // change.
  bool scale() const { return (rv_info.rvWeightInfo.size() > 0); }
  bool penalty() const { return (rv_info.rvWeightInfo.size() > 0); }
  bool shift() const { return (rv_info.rvWeightInfo.size() > 0); }
  // returns true if the implementation of this RV changes with each segment.
  virtual bool iterable() const { return false; } 

  /////////////////////////////////////////////////////////////////////////
  // Initialize with the variable type.
  // The default timeFrame value of -1 indicates a static network.
  // The default value of "hidden" is true.
  // Discrete nodes must be specified with their cardinalities.
  RV(RVInfo& _rv_info,
     unsigned _timeFrame = ~0x0 )
    : rv_info(_rv_info),timeFrame(_timeFrame)
  {}

  virtual ~RV() {;}


  // set up the 'neighbors' member.
  void createNeighborsFromParentsChildren();

  // connect all neighbors other than the ones in the exclude set.
  void connectNeighbors(set<RV*> exclude);

  // moralize the node's parents
  void moralize();

  //////////////////////////////////////////////////////////
  // RV Printing Routines
  //////////////////////////////////////////////////////////

  // optionally print a new line to file.
  void pnl(FILE *f,const bool nl=true) { if (nl) fprintf(f,"\n"); }
  // return either a blank or a NL string.
  const char *const nls(const bool nl=true) { if (nl) return "\n"; else return ""; }

  // print all parent info to given file
  virtual void printParentInfo(FILE* f,bool nl=true);

  // print all children info to given file
  virtual void printChildrenInfo(FILE* f,bool nl=true);

  // print just the name and frame of this RV, optionally a new line.
  void printNameFrame(FILE *f,const bool nl=true) {
    fprintf(f,"%s(%d)%s",name().c_str(),frame(),nls(nl));
  }
  // print just the name frame and value of this RV.
  virtual void printNameFrameValue(FILE *f,bool nl=true) = 0;

  // Printing: show the information about the the current RV, in a
  // relatively concise single-line way.
  virtual void printSelf(FILE *f,bool nl=true) = 0;

  // Printing: show ALL the information about the the current RV in a
  // verbose way.
  virtual void printSelfVerbose(FILE *f) = 0;

  //////////////////////////////////////////////////////////
  // Abstract computing with probabilities
  //////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////

  // Compute the probability of this RV given the current values of
  // the parents, using whatever their values are currently set to.
  // See also begin() and next() iterators below, to iterate
  // through all values of this rv given current parent values.
  virtual void probGivenParents(logpr& p) = 0;

  // Versions that return the probability so require 
  // more temporary creation.
  logpr probGivenParents() { 
    logpr p((void*)0);
    probGivenParents(p);
    return p;
  }


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Parameter and Structure adjustment/checking       /////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  // Sets the parameters determining probGivenParents() to "random" values.
  virtual void makeRandom() = 0;

  ////////////////////////////////////////////////////////////////////////
  // Sets the parameters determining probGivenParents() to "uniform" values.
  virtual void makeUniform() = 0;

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

  ////////////////////////////////////////////////////////////////////////
  // Sample the current distribution setting the random variable to a set
  // of values depending on the current parent values.
  virtual void randomSample() = 0;

  /////////////////////////////////////////
  // EM Learning                         //
  /////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // Increment the statistics with probabilty 'posterior' for the
  // current random variable's parameters, for the case where
  // the random variable and its parents are set to their current
  // values (i.e., the increment corresponds to the currently set
  // parent/child values). 
  virtual void emIncrement(logpr posterior) = 0;

  /////////////////////////////////////////
  // Support for unrolling               //
  /////////////////////////////////////////

  /////////////////////////////////////////  
  // A routine that clones the "shell" of the current random
  // variable. The shell of a random variable is one without any of
  // the the parents and children. This is used in unrolling, where
  // the new (and appropriate) parents/children are set externally to
  // this class. Note that attributes other than parents/children are
  // retained in the clone, and that includes the parameters, the
  // cardinality (if discrete), disposition (hidden),
  // penalty/scale/offset, etc. Only the parents/children are not
  // included.
  virtual RV *cloneRVShell() { 
    RV* rv = create();
    return rv;
  }

  ////////////////////////////////////////////////////////
  // create an empty object of the same type as this (or for whatever
  // derived class this routine lives in). This routine need not fill
  // in any of the member parameters of the object (other than t),
  // rather it just needs to create a RV of the appropriate type.
  virtual RV *create() = 0;

  /////////////////////////////////////////
  // Support for Triangulation           //
  /////////////////////////////////////////

  ////////////////////////////////////////////////////////
  // Return true if all parents of this rv are contained
  // within given set.
  bool allParentsContainedInSet(const set <RV*> givenSet);

  // set the value of the RV, if it happens to be observed, to the
  // appropriate observed value, otherwise a noop.
  virtual void setToObservedValue() {}




};

void printRVSetAndValues(FILE*f,vector<RV*>& locset,const bool nl=true); 
void printRVSetAndValues(FILE*f,sArray<RV*>& locset,const bool nl=true);
void printRVSetAndValues(FILE*f,set<RV*>& locset,const bool nl=true);
void printRVSet(FILE*f,vector<RV*>& locvec,const bool nl=true);
void printRVSet(FILE*f,sArray<RV*>& locset,const bool nl=true);
void printRVSet(FILE*f,const set<RV*>& locset,const bool nl=true);
void printRVSetPtr(FILE*f,set<RV*>& locset,const bool nl=true);




#endif
