/* 
 * tieStructures.h
 * the myriad structures for storing questions, features, decision tree, etc.
 *
 * Written by Simon King <Simon.King@ed.ac.uk>
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
 */

#ifndef TIESTRUCTURES_H
#define TIESTRUCTURES_H


#include <stdlib.h>
#include <string.h>
#include <list>
#include <set>

////////////////////////////////////////////////////////////////////////
//
// defining the features being used
//
////////////////////////////////////////////////////////////////////////
// a feature set definition consists of a number of named features;
// for each named feature, there is a set of possible values
//

// example:
// 1 phonePos       (0 1 2 3 4)
//
// name=phonePos
// string_values=["0", "1", ..."4"]
typedef struct {
  std::string name; // the name of this feature
  std::vector<std::string> string_values;
  std::map<std::string,unsigned> string_value_to_index;
} FeatureDefinitionType;

// example
// my_features1
// 5
// 0 word           (and i oh okay really right sil so the well yes)
// 1 phonePos       (0 1 2 3 4)
// 2 phoneRelPos    (B M E A)
// 3 wordPhoneCount (1 2 3 4 5)
// 4 phoneState     (0 1 2)
//
// name=my_features1
// features[0] contains the FeatureDefinitionType "0 word (and i oh okay really right sil so the well yes)"
typedef struct {
  std::string name; // the name of this feature set
  std::vector<FeatureDefinitionType> features;
  std::map<std::string,unsigned> feature_name_to_index;
} FeatureDefinitionSetType;




////////////////////////////////////////////////////////////////////////
//
// storing the features for the parameters being clustered
//
////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////
// internally, features are stored as unsigned and each parameter has
// its own vector of features
// order of features here must match an associated FeatureDefinitionSetType
//
typedef std::vector<unsigned> FeatureVectorType;
//
// store the features indexed by parameter name
typedef struct {
  std::string name; // the name of this set of feature values
  FeatureDefinitionSetType *feature_definitions; // where to find the definition of feature "name"
  std::map<std::string,FeatureVectorType> values;
} FeatureValueSetType;


void print_feature_maps();

////////////////////////////////////////////////////////////////////////
//
// questions about the values of certain features
//
////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////
typedef struct {
  std::string name;                    // the name of the question
  std::string feature_name;            // what feature is this question about?
  unsigned feature_index;              // ditto, as an index into the FeatureDefinitionSetType named by feature_definition_name
  std::set<unsigned> valueSet;         // what values of feature_name result in the answer "yes"
} QuestionType;

typedef struct {
  std::string name; // the name of this question set
  FeatureDefinitionSetType *feature_definitions; // where to find the feature definitions for this question set
  FeatureValueSetType *feature_values;       // which loaded values should we query?
  std::vector<QuestionType> questions; // the questions
} QuestionSetType;

////////////////////////////////////////////////////////////////////////
//
// decision trees
//
////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////
// a binary branching tree stored as a vector. element 0 is
// unused. element 1 is the root and its children are elements 2
// and 3; element 2's children are 4 and 5, element n's children
// are 2n and 2n+1; element n's parent is n/2 (rounded down).
//
// the elements store the index of the question being asked at
// that node
//
// special entries are:
// -1 leaf node (like other GMTK trees)
// -2 unused node
// -3 node under construction


typedef struct {
  std::string name;
  FeatureDefinitionSetType *feature_definitions; // where to find the feature definitions for this question
  FeatureValueSetType *feature_values;       // which loaded values should we query?
  QuestionSetType *questions; // what questions are we asking
  std::vector<int> tree; // each int refers to one of the questions in the QuestionSetType object

  // after clustering and before tying, we do not know the name of the
  // centroid parameter of each cluster, so we store the cluster index
  // instead
  map<unsigned,unsigned> nodes_to_clusters;


  // which parameter lives in each leaf?  these will be the tied
  // parameter names, to which all items arriving at this leaf should
  // be tied
  map<unsigned,std::string> nodes_to_param_names;


} DecisionTreeType;


#endif
