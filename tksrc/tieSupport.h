/* 
 * tieSupport.h
 * support and convenience functions to support gmtkTie
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

#ifndef TIESUPPORT_H
#define TIESUPPORT_H


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

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"


#include "GMTK_Mixture.h"
#include "GMTK_Tie.h"


#include<list>
#include <regex.h>


#define MAX_REGEX_GROUPS 20


// ----------------------------------------------------------------------
// Clusterable
// a class to contain a parameter being clustered
//
// this is a virtual class and objects of this type cannot be instantiated


// this is too simple - cannot handle cases where a mean is already in
// use by more than one Component or Mixture...??

class Clusterable
{

protected:

  ////////////////////////////////////////////////////////////////////////
  // the name of the parameter associated with this clusterable object
  std::string _name;

  ////////////////////////////////////////////////////////////////////////
  // occupancy "count", will be loaded from saved accumulator file(s)
  logpr _occupancy;

  ////////////////////////////////////////////////////////////////////////
  // unlogged version of _occupancy (to do: take care about underflow)
  double _unlog_occupancy;

  ////////////////////////////////////////////////////////////////////////
  // a pointer to the dissimilarity function that will be used when
  // this clusterable object is compared to another one
  float (Clusterable::*dissimilarityFunction)(Clusterable *c);

public:

  Clusterable();
  virtual ~Clusterable();

  inline float dissimilarity(Clusterable *c){
    return (this->*dissimilarityFunction)((Clusterable*)c);
  }

  inline const std::string& name(){return _name;};
  inline const double occupancy(){return _unlog_occupancy;};
};


// ----------------------------------------------------------------------
// ClusterableObject<T>
// a class inherited from Clusterable to hold a parameter of a
// specific type
//
// this class must be instantiated with an appropriate type for T (see
// below)

template<typename T>
class ClusterableObject: public Clusterable
{
private:

protected:

  ////////////////////////////////////////////////////////////////////////
  // a pointer to the parameter object that this class contains
  T *param_ptr;

public:

  ClusterableObject<T>();
  ClusterableObject<T>(const std::string &n, T *m, GMTK_Tie::DissimilarityMeasureType d);
  ~ClusterableObject<T>();
  
};

// ----------------------------------------------------------------------
// Now some specific instantiations of the templated class
// ClusterableObject<T>
//
// these classes each implement one or more dissimilarity functions
// appropriate for the type of T
//
// the dissimilarityFunction will be set to point at one of them by
// the constructor

class ClusterableMean: public ClusterableObject<MeanVector>
{
private:
  void set_dissimilarity_function(GMTK_Tie::DissimilarityMeasureType d);
public:

  ////////////////////////////////////////////////////////////////////////
  // it makes no sense to construct an empty object
  ClusterableMean() {error("Cannot construct a ClusterableMean this way!");};
  // this is the constructor to use
  ClusterableMean(const std::string &n, MeanVector *m, GMTK_Tie::DissimilarityMeasureType d);

  ////////////////////////////////////////////////////////////////////////
  // available dissimilarity measures for ClusterableMean are:
  float Euclidean_distance(Clusterable *c);
  // add more here...
};


class ClusterableMixture: public ClusterableObject<Mixture>
{
private:
  void set_dissimilarity_function(GMTK_Tie::DissimilarityMeasureType d);
public:

  ////////////////////////////////////////////////////////////////////////
  // it makes no sense to construct an empty object
  ClusterableMixture() {error("Cannot construct a ClusterableMixture this way!");};
  // this is the constructor to use
  ClusterableMixture(const std::string &n, Mixture *m, GMTK_Tie::DissimilarityMeasureType d);

  ////////////////////////////////////////////////////////////////////////
  // available dissimilarity measures for ClusterableMixture are:
  float CrossLogLikelihoodOfMeans_distance(Clusterable *c);
  // add more here...

};


// ----------------------------------------------------------------------
// Cluster
// a structure for holding a list of Clusterable objects, i.e. a cluster

typedef struct {
  // a unique identifier
  unsigned ident;

  ////////////////////////////////////////////////////////////////////////
  // the first item in this list will be used as the centroid, in
  // cases where the centroid is an existing member of the cluster
  std::list<Clusterable*> items;

  ////////////////////////////////////////////////////////////////////////
  // to do: we could also store any constructed centroid (currently,
  // during clustering, we don't keep recomputing this and just use an
  // existing item, even if the centroid selection method indicates
  // that a centroid should be constructed)
  ////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  // a measure of the cluster's "impurity" (for decision tree
  // clustering) or "size" (for agglomerative clustering)
  float size;

  ////////////////////////////////////////////////////////////////////////
  // occupancy count summed over all items
  double occupancy;

} Cluster;



// ----------------------------------------------------------------------
// helper functions for operating on clusters

// re-arrange the items in the list above so that the centroid is the first item in the list
bool set_cluster_centroid(Cluster *c, GMTK_Tie::CentroidType method, bool force=false);

// calculate what the size of the cluster would be if two clusters were merged
float merged_cluster_size(std::list<Cluster>::iterator ci, std::list<Cluster>::iterator ci2,
			  GMTK_Tie::ClusterSizeMethodType ClusterSizeMethod,
			  GMTK_Tie::CentroidType Centroid);


// initial construction of a table containing all pairwise cluster merged sizes
void build_merged_sizes_table(sArray< sArray<float> > *merged_sizes,
			      std::list<Cluster> *clusters,
			      GMTK_Tie::ClusterSizeMethodType ClusterSizeMethod,
			      GMTK_Tie::CentroidType Centroid);



// update this table (to call after two clusters have been merged,
// with ci pointing to the newly constructed cluster)
void update_merged_sizes_table(sArray< sArray<float> > *merged_sizes,
			       std::list<Cluster> *clusters,
			       GMTK_Tie::ClusterSizeMethodType ClusterSizeMethod,
			       GMTK_Tie::CentroidType Centroid,
			       std::list<Cluster>::iterator ci);





// ----------------------------------------------------------------------
// expanding and matching user-supplied parameter names or regexs

void expand_param_names(GMParms::ObjectMapType* the_map, std::string &param_expression, 
			std::vector<string> *params);
void expand_param_names(GMParms::ObjectMapType* the_map, regex_t *param_expression, 
			std::vector<string> *params);
bool check_parameter_exists(GMParms::ObjectMapType* the_map, std::string param_name);


// ----------------------------------------------------------------------
// lookup from the maps in GM_Parms, with error checking

Mixture*    find_Mixture(std::string &name);
Component*  find_Component(std::string &name);
Dense1DPMF* find_Dense1DPMF(std::string &name);
MeanVector* find_MeanVector(std::string &name);
DiagGaussian* find_DiagGaussian(std::string &name);

// ----------------------------------------------------------------------
// checking properties of parameter objects

// are all Components in this Mixture of type DiagGaussian?
bool all_DiagGaussian(Mixture* mixture);


// ----------------------------------------------------------------------
// finding shared usage of parameter objects

// which Components are using this MeanVector?
unsigned find_Components_using_MeanVector(const MeanVector* const mean, const std::string mean_name, 
					  std::list<Component*> *components, 
					  std::list<std::string> *component_names);

// ----------------------------------------------------------------------
// finding sub-components, with sanity and type checking

// find the Components of a Mixture, with optional check that nothing is currently shared
std::vector<Component*> find_Components_of_Mixture(Mixture *mixture, bool forbid_sharing=true);

// these include a check that the Mixture has only Gaussian components
std::vector<MeanVector*> find_MeanVectors_of_Mixture(Mixture *mixture);
std::vector<MeanVector*> find_MeanVectors_of_Mixture(std::string &name);

// these include check that Mixture has exactly one Gaussian component
MeanVector* find_MeanVector_of_Mixture(Mixture *mixture);
MeanVector* find_MeanVector_of_Mixture(std::string &name);

// these include a check that the Component is DiagGaussian (may extend later)
MeanVector* find_MeanVector_of_Component(Component *component);
MeanVector* find_MeanVector_of_Component(std::string &name);

// these include a check that the Component is DiagGaussian 
MeanVector* find_MeanVector_of_DiagGaussian(DiagGaussian *diag_gaussian);
MeanVector* find_MeanVector_of_DiagGaussian(std::string &name);


// ----------------------------------------------------------------------
// regular expression helpers
//
// written because the interface to posix regex is ugly and needs to
// be hidden

regex_t* compile_regex(const std::string &str);
std::string get_regerror (int errcode, regex_t *compiled);
bool match(regex_t *compiled_regex, const std::string &str);


// ----------------------------------------------------------------------
// general helper functions
//
// some could be removed - were originally implemented using run-time
// type checking, before I realised that objects do know what type
// they are!

inline std::string get_Component_type(Component *component)
{
  return component->typeName();
  //return std::string(typeid(*(component)).name());
}

inline bool is_DiagGaussian(Component *component)
{
  return (component->typeName() == "Diag Gaussian");
  //return (get_Component_type(component).find("DiagGaussian") != std::string::npos);
}


// create a unique new name (one that doesn't yet exist in the given map)
std::string new_name(std::string basename, map< string, unsigned > *map);


#endif
