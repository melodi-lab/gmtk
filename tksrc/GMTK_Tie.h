/* 
 * GMTK_Tie.h
 *
 * This class contains the main functionality for clustering, tying and untying parameters
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

#ifndef GMTK_TIE_H
#define GMTK_TIE_H


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"

#include "GMTK_GMParms.h"

#include <list>
#include <regex.h>


class GMTK_Tie {

public:
  ////////////////////////////////////////////////////////////////////////
  // constructor and destructor
  GMTK_Tie(GMParms *GM_Parms_ptr);
  ~GMTK_Tie();


  ////////////////////////////////////////////////////////////////////////
  // the following enumerated types are all the available commands and
  // options that the user can specify in the tying command file
  //
  // conversion from strings to enums happens as the file is parsed
  // using the maps
  ////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  // commands
  enum CommandType {
    CT_Unknown=0,
    CT_Tie=1,       // just tie, no questions asked
    CT_Cluster=2,   // bottom-up cluster, then tie
    CT_DTcluster=3, // decision-tree cluster, then tie
    CT_Untie=4      // does what it says on the tin
  };
  map<std::string, CommandType> string_to_CommandType;
  map<CommandType, std::string> CommandType_to_string;

  ////////////////////////////////////////////////////////////////////////
  // the type of parameter a command refers to
  enum ParameterType {
    PT_Unknown=0,
    PT_Mixture=1,
    PT_Mean=2,
    PT_Collection=3
  };
  map<std::string, ParameterType> string_to_ParameterType;
  map<ParameterType, std::string> ParameterType_to_string;


  ////////////////////////////////////////////////////////////////////////
  // pairwise dissimilarity measures between parameters (note: not
  // necessarily distance metrics)
  enum DissimilarityMeasureType{
    DMT_Unknown=0,

    // just uses the distance between two mean vectors; only works for
    // single Component Mixtures
    DMT_Euclidean=1,

    // for Mixtures with more than one Component, this simply
    // averages the DMT_Euclidean distance (between pairs of
    // corresponding means) across the Mixture components
    // 
    // for single Component Mixtures, it is identical to DMT_Euclidean
    //
    // note: no special care is taken to permute the order of the
    // Mixture components to get the most favourable pairings of mean
    // vectors
    DMT_AverageEuclidean=2,

    // use covariance of first Gaussian to normalize the Euclidean
    // distance; only works for single Component Mixtures; not
    // symmetric
    DMT_Mahalanobis=3,

    // like DMT_Mahalanobis, but made symmetric by summing across each
    // direction
    DMT_SymmetricMahalanobis=4,

    // K-L divergence; see, for example,
    // http://en.wikipedia.org/wiki/Kullback-Leibler_divergence
    DMT_KullbackLeiblerDivergence=5,

    // the K-L divergence summed over each direction (see under
    // "Symmetrised divergence" on Wikipedia!)
    DMT_SymmetricKullbackLeiblerDivergence=6,

    // compute the log likelihood of each mean vector under the
    // opposite PDF, and sum over all means
    DMT_CrossLogLikelihoodOfMeans=7,

    // emulate HTK's measures, as specified at
    // http://htk.eng.cam.ac.uk/prot-docs/HTKBook/node283_mn.html
    //
    // this is a Mahalanobis-like measure but uses the product of the
    // two Gaussians' standard deviations in place of the variance of
    // just one of them, so it is symmetric; only works for single
    // Component Mixtures of diagonal Gaussians
    //
    // HTK does something special for "fully tied mixture" systems,
    // but we are not emulating those here
    //
    // in all other cases, HTK uses DMT_CrossLikelihoodOfMeans (it
    // additionally restricts this to pairs of Mixtures with the same
    // number of components, which we currently do to - but may
    // change this)
    DMT_EmulateHTK=8


    // could add many more, e.g.:
    // * Goldberger & ? , Interspeech 2005
    // * squared Euclidean
    // * average squared Euclidean
    // * Bhattacharyya
  };
  map<std::string, DissimilarityMeasureType> string_to_DissimilarityMeasureType;
  map<DissimilarityMeasureType, std::string> DissimilarityMeasureType_to_string;


  ////////////////////////////////////////////////////////////////////////
  // ways of choosing the centroid of a cluster
  enum CentroidType{
    // use any existing cluster member
    CNT_Unknown=0, 

    // use any existing cluster member
    CNT_Arbitrary=1,

    // use the existing cluster member with the lowest total
    // dissimilarity to all other cluster members
    CNT_UseExistingCentroid=2,

    // create a new member at the exact centre of the cluster
    CNT_CreateCentroid_averageSingleComponentMixtures=3, 

    // not yet implemented
    CNT_CreateCentroid_permuteThenAverage=4,
    CNT_CreateCentroid_smartMerge=5,

    // create a new member at the exact centre of the cluster
    CNT_CreateCentroid_averageMeanVector=6,

    // when tying states, HTK does this:
    //
    // "the state with the largest total value of gConst in stream 1
    // (indicating broad variances) and the minimum number of defunct
    // mixture weights (see MU command) is selected from the item list
    // and all states are tied to this typical state." (HTK manual /
    // HHEd / TI command)
    //
    // http://htk.eng.cam.ac.uk/prot-docs/HTKBook/node299_mn.html
    //
    // and we will allow this to be selected for tying Mixtures of
    // Gaussians or single Gaussian Components
    //
    // for tying anything other than a state, HTK uses the last item
    // in the user-specified list as the centroid; we will use the
    // first item in the list (it's arbitrary after all)
    CNT_EmulateHTK=7,

    // these are non-user options which CNT_EmulateHTK will be changed
    // to, depending on parameter type (this avoids various functions
    // needing to know the parameter type)
    CNT_EmulateHTKMixturesOfGaussians=8,
    CNT_EmulateHTKOther=9// will be same as CNT_Arbitrary

  };
  map<std::string, CentroidType> string_to_CentroidType;
  map<CentroidType, std::string> CentroidType_to_string;


  ////////////////////////////////////////////////////////////////////////
  // ways of measuring the size of a cluster
  enum ClusterSizeMethodType {

    // find the farthest apart pair of cluster members (the "diameter"
    // of the cluster), as used in HTK
    CSM_MostDissimilarPair=0,
    CSM_EmulateHTK=1, // will be the same as CSM_MostDissimilarPair

    // take average over all members of dissimilarity to centroid (may
    // require recomputation of centroid)
    CSM_AverageDissimilarityToCentroid=2
  };
  map<std::string, ClusterSizeMethodType> string_to_ClusterSizeMethodType;
  map<ClusterSizeMethodType, std::string> ClusterSizeMethodType_to_string;


  ////////////////////////////////////////////////////////////////////////
  // all possible options that can be given to commands (not all of
  // them apply to all commands)
  typedef struct {

    // the stopping criterion for data driven agglomerative clustering
    // or decision tree clustering
    //
    // in the former case, clusters grow, starting from a size of 0.0
    //
    // in the latter case, clusters shrink (they are repeatedly split
    // in two); this criterion terminates tree growth only for one
    // branch of the tree, but other branches will continue to grow
    float MaxClusterSize;

    // the target number of parameters to aim for in decision tree
    // clustering (i.e. the number of leaves of the tree)
    //
    // normally used *instead* of MaxClusterSize (if both are given,
    // then tree growing is terminated for that branch if *either* of
    // them is satisfied - this may not be helpful behaviour?)
    unsigned MaxNumberParameters;

    // the first method for outlier removal
    double MinOccupancyCount;

    // the second method for outlier removal
    int MinClusterMembers; 

    // when parameters get tied, they will be renamed using this
    // prefix (not yet implemented)
    std::string NewParameterNamePrefix;

    // the measure used between pairs of cluster members
    DissimilarityMeasureType DissimilarityMeasure;

    // how the centroid should be chosen
    CentroidType Centroid;

    // whether weighting by occupancy should be used when computing
    // the centroid (not yet implemented)
    bool CentroidOccupancyWeighting;

    // how the size of the cluster should be measured
    ClusterSizeMethodType ClusterSizeMethod;

    // whether weighting by occupancy should be used when computing
    // cluster size (not yet implemented)
    bool ClusterSizeOccupancyWeighting;

  } Options;

  

  ////////////////////////////////////////////////////////////////////////
  // a struct to hold a single command
  typedef struct {

    // the command itself
    CommandType command;

    // the type of parameter it operates on (currently restricted to a
    // single type - in future we may want to be able to operate on a
    // heterogenous set of parameters, e.g. Diag and Full covar
    // Gaussians)
    ParameterType param_type;

    // the options for this command
    Options options;

    // a regular expression for specifing which parameters from
    // GM_Parms we apply this command to
    regex_t* param_expression;

    // the parameter names that match param_expression
    std::vector<std::string> params;

  } Command;

  ////////////////////////////////////////////////////////////////////////
  // an ordered list of commands - they will be executed in this order
  typedef std::vector<Command> CommandListType;
  CommandListType commands;

  ////////////////////////////////////////////////////////////////////////
  // execute a command
  bool execute_command(unsigned command_index);



private:

  ////////////////////////////////////////////////////////////////////////
  // we need to know about the global parameter object
  GMParms *GM_Parms;


  ////////////////////////////////////////////////////////////////////////
  // parsing the command file
  ////////////////////////////////////////////////////////////////////////
  void parse_option(CommandType ctype, Options &o, std::string opt);
  void fill_in_default_values(CommandType ctype, Options &opt);
  void enter_tie_option(GMTK_Tie::Options &opt, std::string &key, std::string &value);
  void enter_cluster_option(GMTK_Tie::Options &opt, std::string &key, std::string &value);
  void enter_DTcluster_option(GMTK_Tie::Options &opt, std::string &key, std::string &value);


public:


  ////////////////////////////////////////////////////////////////////////
  // load commands, create 'commands' list, fill in default values for
  // options, detect basic errors
  bool read_commands(iDataStreamFile& is);
  bool need_occupancy_counts();
  
  ////////////////////////////////////////////////////////////////////////
  // for every command, expand the regex param_expression and store
  // the results
  void find_matching_command_parameters();
  // validate the parameters (only simple error checking is done)
  bool validate_command_parameters();
  
  
  ////////////////////////////////////////////////////////////////////////
  // saving
  // - not yet implemented


  ////////////////////////////////////////////////////////////////////////
  // printing
  void print_Command(Command &cmd, int command_index, IM::VerbosityLevels v=IM::Tiny);
  void print_Options(Options &opt, IM::VerbosityLevels v=IM::Tiny);


  ////////////////////////////////////////////////////////////////////////
  // tying functions - either called as the direct result of a user
  // command, or after clustering has been used to determine groups of
  // parameters to tie
  bool tie_Mixtures(std::vector<std::string> param_expressions, CentroidType method, bool expand_expressions=true);
  bool untie_Mixtures(std::vector<std::string> param_expressions, bool expand_expressions=true);
  bool tie_Means(std::vector<std::string> param_expressions, CentroidType method, bool expand_expressions=true);

  ////////////////////////////////////////////////////////////////////////
  // to do - add more types:

  // * Covariance vectors/matrices (question: can we tie the diagonal
  //   and off-diagonal elements spearately?)

  // * dPmfs (problem: will need to tie sub-matrices, e.g. for tyng
  //   the transition matrix across all triphone models of a phone)

  // the main clustering functions - these are as generic as
  // possible (i.e. independent of parameter type and dissimilarity
  // measure)


  ////////////////////////////////////////////////////////////////////////
  // clustering
  bool data_driven_cluster(unsigned command_index, std::vector<std::vector<std::string> > &clustered_params);
  bool decision_tree_cluster(unsigned command_index, std::vector<std::vector<std::string> > &clustered_params); // not implemented


};


ostream& operator<<(ostream&s, GMTK_Tie::Command &cmd);






#endif
