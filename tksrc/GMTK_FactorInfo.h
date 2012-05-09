/*
 * GMTK_FactorInfo.h
 *
 * A class that keeps all the information about a factor in a
 * structure file, unique to its location in the structure.  (i.e.,
 * this is a factor over a set of random variables, and the
 * information kept here is common to *all* instantiations of such a factor
 * over random variables that share this factor).


 *
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2006, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * $Header$
 *
 */

#ifndef GMTK_FACTORINFO_H
#define GMTK_FACTORINFO_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

class RV;
class DiscRV;
class ContRV;
class FileParser;
class GMTemplate;
class Partition;
class BoundaryTriangulate;

class RV;
class DiscRV;

#include "GMTK_CPT.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_RVInfo.h"

/*
 * the next structure is a set of undiredted edge 'factors', 'cliques'
 * (not max cliques), or completed sets of variables that are used for
 * adding edges to the graph and/or producing hybrid
 * directed-undirected graphical models. They are also used for
 * undirected style links between variables in the graph, where the
 * constraint among the variables involved can either be fixed
 * unlearned of the form:
 *     1) a constant form of constraint, specified in the .str file that
 *        inference knows how to deal with
 *
 *     2) a variable form of constraint, specified using a decision
 *        tree boolean evaluator or a learned form consisting of a
 *        log-linear model:
 *
 *     3) a set of feature functions of the variables involved in the clique along with
 *        weights that are learned along with EM using an iterative scaling form of algorithm.
 * 
 */
class FactorInfo {

  friend class RVInfo;
  friend class FileParser;
  friend class GMTemplate;
  friend class Partition;
  friend class BoundaryTriangulate;
  friend class JunctionTree;
  friend class StructPage;
  friend class VizNode;
  friend class GraphicalModel;
  friend class FactorClique;

  friend class RV;
  friend class DiscRV;
  friend class ObsDiscRV;
  friend class Sw_ObsDiscRV;
  friend class ScPnSh_ObsDiscRV;
  friend class ScPnSh_Sw_ObsDiscRV;
  friend class HidDiscRV;
  friend class Sw_HidDiscRV;
  friend class ScPnSh_HidDiscRV;
  friend class ScPnSh_Sw_HidDiscRV;
  friend class ObsContRV;
  friend class Sw_ObsContRV;
  friend class ScPnSh_ObsContRV;
  friend class ScPnSh_Sw_ObsContRV;
  friend class ScPnShRV;
  friend class ContRV;
  friend class SwRV;
  friend class SwDiscRV;

public:

  ////////////////////////////////////////////////////////////
  // types of factors
  enum FactorType { ft_unknown, // used to add undirected edges and for error checking
		    ft_symmetricConstraint,
		    ft_directionalConstraint,
		    ft_softConstraint };

  // types of symmetricConstraints
  enum SymmetricConstraintType {
    sct_allVarsEqual,     // all vars must be equal
    sct_allVarsUnequal,   // all vars must be unequal to each other (no pairwise equality allowed)
    sct_varsNotEqual,     // not all vars must be equal (i.e., at least two that are different)
    sct_varsSumTo,        // all vars must sum to n
    sct_varsMultiplyTo,   // all vars must multiply to n
    sct_varsSumMod,       // sum of all vars mod m = n
    sct_varsSatisfy       // vars must satisfy a DT constraint
  };

  // types of directionalConstraints
  enum DirectionalConstraintType {
    dct_functionOf        // one var is deterministic function of the others.
  };

  // types of softConstraints
  enum SoftConstraintType {
    fct_table,           // vars must live in a table, which provides a soft value
    fct_logLinear        // vars collectively index into a logLinear model with features.
  };

private:


  ///////////////////////////////////////////////////////////
  // data associated with a factor

  // the frame where it was defined (so is part of name really)
  unsigned frame;
  // line number of file where this factor was first declared
  unsigned fileLineNumber;
  // file name where this r.v. is first defined.
  string fileName;
  // the factor's name, used to identify it later.
  string name;

  // the list of variables (in file order) defined in this factor.
  vector < RVInfo::rvParent > variables;

  // pointers to the adam & eve actual RVs when they exist
  // (see RVInfo's rv for documentation on where these rvs come from)
  vector < RV* > rvs;

  // the type of the factor
  FactorType fType;

  // various structures active depending on fType.

  // if this is a SymmetricConstraint 
  struct SymmetricConstraintInfo {
    SymmetricConstraintType symmetricConstraintType;
    // use variable n, for 
    //     varsSumTo(n) --- all vars must sum to n 
    //     varsMultiplyTo(n) --- all vars must multiply to n
    unsigned n; 
    // use variable m,n for
    // varsSumMod(m,n) meaning, (sum(vars) mod m) = n to satisfy
    unsigned m;
    // name of a DT, must evaluate to non-zero to satisfy.
    string mappingDT;
  } symmetricConstraintInfo;

  // if this is a directionalConstraint
  struct DirectionalConstraintInfo {
    // Here, we satisfy if child = functionOf(parent variable list) using mapping
    vector < RVInfo::rvParent > parents;
    RVInfo::rvParent child;
    // name of a DT, to implement the function.
    string mappingDT;
  } directionalConstraintInfo;

  // if this is a softConstraint
  struct SoftConstraintInfo {
    SoftConstraintType softConstraintType;
    // name of either table or log-linear object model.
    string name;
  } softConstraintInfo;



public:

  /////////////////////////////////////////////////////////
  // constructor
  FactorInfo() { clear(); }

  // copy constructor
  FactorInfo(const FactorInfo&v) {
    frame = v.frame;
    fileLineNumber = v.fileLineNumber;
    fileName = v.fileName;
    name = v.name;
    variables = v.variables;
    fType = v.fType;
    if (fType == ft_symmetricConstraint) {
      symmetricConstraintInfo = v.symmetricConstraintInfo;
    } else if (fType == ft_directionalConstraint) {
      directionalConstraintInfo = v.directionalConstraintInfo;
    } else {
      softConstraintInfo = v.softConstraintInfo;
    }
  }

  // clear out the current RV structure when we
  // are parsing and encounter a new RV.
  void clear() {
    fType = ft_unknown;
    frame = ~0x0;
    fileLineNumber = ~0x0;
    fileName.clear();
    name.clear();
    variables.clear();
  }


};




#endif
