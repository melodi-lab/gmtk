/* 
 * tieSupport.cc
 * support and convenience functions to support gmtkTie
 *
 * Written by Simon King <Simon.King@ed.ac.uk>
 *
 * Copyright (C) 2006 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 */

#include "tieSupport.h"



/*-
 *-----------------------------------------------------------------------
 * Clusterable::Clusterable
 *      constructor - never call this directly
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *     none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
Clusterable::Clusterable()
{
  _name="";
}


/*-
 *-----------------------------------------------------------------------
 * Clusterable::~Clusterable
 *      destructor
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
 *      none
 *
 *-----------------------------------------------------------------------
 */
Clusterable::~Clusterable()
{
  // nothing currently needs freeing
}




/*-
 *-----------------------------------------------------------------------
 * ClusterableObject<T>::~ClusterableObject
 *      template destructor
 *
 * Preconditions:
 *      the template class should be instantiated with a specific type
 *      for T
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
template<typename T>
ClusterableObject<T>::~ClusterableObject()
{
  // nothing currently needs freeing
}



/*-
 *-----------------------------------------------------------------------
 * ClusterableObject<T>::ClusterableObject
 *      template constructor
 *
 * Preconditions:

 *      1) the template class should be instantiated with a specific
 *      type for T 
 *      2) m should be pointer to a valid T object, and its
 *      accumulators shoukd have been loaded
 *
 * Postconditions:
 *      1) the object is constructed
 *      2) its name is set
 *      3) the occupancy stats are set from the accumulator associated with m
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
template<typename T>
ClusterableObject<T>::ClusterableObject(const std::string &n, T *m, GMTK_Tie::DissimilarityMeasureType d)
{
  _name=n;

  _occupancy=m->get_accumulatedProbability();

  if(_occupancy.essentially_zero())
    infoMsg(IM::Max,"Essentially zero occupancy for %s\n",m->name().c_str());

  _unlog_occupancy = _occupancy.unlog();
  param_ptr=m;

  //infoMsg(IM::Huge,"Created a ClusterableObject<T> with occ %f\n",_unlog_occupancy);

}


/*-
 *-----------------------------------------------------------------------
 * ClusterableMean::ClusterableMean
 *      constructor
 *
 * Preconditions:
 *      1) d should be a valid GMTK_Tie::DissimilarityMeasureType
 *      2) plus preconditions for
           ClusterableObject<T>::ClusterableObject(const std::string
           &n, T *m, GMTK_Tie::DissimilarityMeasureType d)
 *
 * Postconditions:
 *      1) dissimilarityFunction points to the supplied function
 *      2) plus postconditions for 
 *         ClusterableObject<T>::ClusterableObject(const std::string
 *         &n, T *m, GMTK_Tie::DissimilarityMeasureType d)
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
ClusterableMean::ClusterableMean(const std::string &n, MeanVector *m, GMTK_Tie::DissimilarityMeasureType d) : ClusterableObject<MeanVector>(n,m,d)
{
  set_dissimilarity_function(d);
}



/*-
 *-----------------------------------------------------------------------
 * ClusterableMean::set_dissimilarity_function
 *      sets dissimilarityFunction to point at the appropriate function
 *
 * Preconditions:
 *      d should be a valid GMTK_Tie::DissimilarityMeasureType
 *
 * Postconditions:
 *      dissimilarityFunction points to the supplied function
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
ClusterableMean::set_dissimilarity_function(GMTK_Tie::DissimilarityMeasureType d)
{
  // set the pointer to the dissimilarity function
    switch(d){

    case GMTK_Tie::DMT_Unknown:
      error("Unknown DissimilarityMeasure in ClusterableMean constructor");
      break;

    case GMTK_Tie::DMT_Euclidean:
      dissimilarityFunction=(float (Clusterable::*)(Clusterable*))(&ClusterableMean::Euclidean_distance);
      break;

    case GMTK_Tie::DMT_EmulateHTK:
      error("EmulateHTK is not a valid DissimilarityMeasure for clustering MeanVectors");
      break;

    default:
      error("Inappropriate DissimilarityMeasure specified for clustering MeanVectors");
      break;
    }

}



/*-
 *-----------------------------------------------------------------------
 * ClusterableMean::Euclidean_distance
 *      computes the Euclidean distance from this ClusterableMean to
 *      another supplied ClusterableMean
 *
 * Preconditions:
 *      1) set_dissimilarity_function must have been called
 *      2) cc must point at a ClusterableMean object
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      float
 *
 *-----------------------------------------------------------------------
 */
float
ClusterableMean::Euclidean_distance(Clusterable *cc)
{

  // to do : assert type checking

  MeanVector* my_mean_vector=(MeanVector*)param_ptr;
  MeanVector* other_mean_vector=(MeanVector*)( ((ClusterableMean*)cc)->param_ptr);

  assert(my_mean_vector->means.len() == other_mean_vector->means.len());

  float e=0.0;

  for(signed i=0;i<my_mean_vector->means.len();i++){
    //cerr << m1->means[i] << "," m2->means[i]<<endl;
    e=e+ pow((double)(my_mean_vector->means[i] - other_mean_vector->means[i]),2.0);
  }
  return sqrt(e);
}











/*-
 *-----------------------------------------------------------------------
 * ClusterableDiagGaussian::ClusterableDiagGaussian
 *      constructor
 *
 * Preconditions:
 *      1) d should be a valid GMTK_Tie::DissimilarityMeasureType 
 *      2) plus preconditions for
           ClusterableObject<T>::ClusterableObject(const std::string
           &n, T *m, GMTK_Tie::DissimilarityMeasureType d)
 *
 * Postconditions:
 *      1) dissimilarityFunction points to the supplied function
 *      2) plus postconditions for 
 *         ClusterableObject<T>::ClusterableObject(const std::string
 *         &n, T *m, GMTK_Tie::DissimilarityMeasureType d)
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
ClusterableDiagGaussian::ClusterableDiagGaussian(const std::string &n, DiagGaussian *m, GMTK_Tie::DissimilarityMeasureType d) : ClusterableObject<DiagGaussian>(n,m,d)
{
  set_dissimilarity_function(d);
}



/*-
 *-----------------------------------------------------------------------
 * ClusterableDiagGaussian::set_dissimilarity_function
 *      sets dissimilarityFunction to point at the appropriate function
 *
 * Preconditions:
 *      d should be a valid GMTK_Tie::DissimilarityMeasureType
 *
 * Postconditions:
 *      dissimilarityFunction points to the supplied function
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
ClusterableDiagGaussian::set_dissimilarity_function(GMTK_Tie::DissimilarityMeasureType d)
{
  // set the pointer to the dissimilarity function
    switch(d){

    case GMTK_Tie::DMT_Unknown:
      error("Unknown DissimilarityMeasure in ClusterableMixture constructor");
      break;
    case GMTK_Tie::DMT_EmulateHTK:
      dissimilarityFunction=NULL; // not needed
      break;
    default:
      error("Inappropriate DissimilarityMeasure specified for a ClusterableDiagGaussian");
      break;
    }
}










/*-
 *-----------------------------------------------------------------------
 * ClusterableMixture::ClusterableMixture
 *      constructor
 *
 * Preconditions:
 *      1) d should be a valid GMTK_Tie::DissimilarityMeasureType
 *      2) plus preconditions for
           ClusterableObject<T>::ClusterableObject(const std::string
           &n, T *m, GMTK_Tie::DissimilarityMeasureType d)
 *
 * Postconditions:
 *      1) dissimilarityFunction points to the supplied function
 *      2) plus postconditions for 
 *         ClusterableObject<T>::ClusterableObject(const std::string
 *         &n, T *m, GMTK_Tie::DissimilarityMeasureType d)
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
ClusterableMixture::ClusterableMixture(const std::string &n, Mixture *m, GMTK_Tie::DissimilarityMeasureType d) : ClusterableObject<Mixture>(n,m,d)
{
  set_dissimilarity_function(d);
}



/*-
 *-----------------------------------------------------------------------
 * ClusterableMixture::set_dissimilarity_function
 *      sets dissimilarityFunction to point at the appropriate function
 *
 * Preconditions:
 *      d should be a valid GMTK_Tie::DissimilarityMeasureType
 *
 * Postconditions:
 *      dissimilarityFunction points to the supplied function
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
ClusterableMixture::set_dissimilarity_function(GMTK_Tie::DissimilarityMeasureType d)
{
  // set the pointer to the dissimilarity function
    switch(d){

    case GMTK_Tie::DMT_Unknown:
      error("Unknown DissimilarityMeasure in ClusterableMixture constructor");
      break;
    case GMTK_Tie::DMT_CrossLogLikelihoodOfMeans:
    case GMTK_Tie::DMT_EmulateHTK:
      dissimilarityFunction=(float (Clusterable::*)(Clusterable*))(&ClusterableMixture::CrossLogLikelihoodOfMeans_distance);
      break;
    default:
      error("Inappropriate DissimilarityMeasure specified for a ClusterableMixture");
      break;
    }

}



/*-
 *-----------------------------------------------------------------------
 * ClusterableMixture::CrossLogLikelihoodOfMeans_distance
 *      returns the total log likelihood of all the means in both
 *      Mixtures, under the PDF of the other Mixture
 *
 * Preconditions:
 *      c must point to a ClusterableMixture object, which in turn
 *      must contain a Mixture of Gaussian components
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      float
 *
 *-----------------------------------------------------------------------
 */
float
ClusterableMixture::CrossLogLikelihoodOfMeans_distance(Clusterable *c)
{
  // only for Mixtures of Gaussians
  // to do: assert some things here

  // for the mean of in each Component in this Mixture, compute its log
  // likelihood under the other Mixture

  // take sum over all Components

  // then sum over reverse direction

  std::vector<MeanVector*> means1 = find_MeanVectors_of_Mixture((Mixture*)param_ptr);
  std::vector<MeanVector*> means2 = find_MeanVectors_of_Mixture((Mixture*)( ((ClusterableMixture*)(c))->param_ptr));

  logpr rval;
  rval.set_to_zero();
  std::vector<MeanVector*>::iterator i;
	
  // values not needed !?
  Data32* const base = NULL;
  int stride=0;

  for (i=means1.begin();i!=means1.end();i++)
    rval += ((Mixture*)( ((ClusterableMixture*)(c))->param_ptr))->log_p((*i)->means.ptr,base,stride);

  for (i=means2.begin();i!=means2.end();i++)
    rval += ((Mixture*)param_ptr)->log_p((*i)->means.ptr,base,stride);
  
  // this converts double to float - possible precision problems??
  return rval.val();
}





/*-
 *-----------------------------------------------------------------------
 * check_parameter_exists
 *      checks that the supplied regex matches at least one parameter
 *      in the supplied map
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      if return val is true, we know that a matching parameter exists
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      bool
 *
 *-----------------------------------------------------------------------
 */
bool 
check_parameter_exists(GMParms::ObjectMapType* the_map, std::string param_expression)
{
  //  just look for the first matching param at this point

  regex_t* param_expression_compiled=compile_regex(param_expression);

  for (GMParms::ObjectMapType::iterator i=the_map->begin();i!=the_map->end();i++)
    if ( match(param_expression_compiled,i->first) )
      return true;
	
  regfree(param_expression_compiled);
  return false;
}



/*-
 *-----------------------------------------------------------------------
 * expand_param_names
 *      finds all matching parameter names from the supplied map using
 *      the regex given as a string, or as a pre-compiled regex
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      params is appended with all matching parameter names
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void 
expand_param_names(GMParms::ObjectMapType* the_map, std::string &param_expression, std::vector<string> *params)
{
  regex_t* param_expression_compiled=compile_regex(param_expression);
  expand_param_names(the_map,param_expression_compiled,params);
  regfree(param_expression_compiled);
}

void 
expand_param_names(GMParms::ObjectMapType* the_map, regex_t *param_expression_compiled, std::vector<string> *params)
{
  for (GMParms::ObjectMapType::iterator i=the_map->begin();i!=the_map->end();i++){
    if ( match(param_expression_compiled,i->first) )
      params->push_back(i->first);
  }

}


/*-
 *-----------------------------------------------------------------------
 * find_Mixture
 *      finds a Mixture in GM_Parms, given its name, with error checking
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      we know that the named Mixture exists
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a pointer to the Mixture
 *
 *-----------------------------------------------------------------------
 */
Mixture* 
find_Mixture(std::string &name)
{
  GMParms::ObjectMapType::iterator i=GM_Parms.mixturesMap.find(name);
  if (i == GM_Parms.mixturesMap.end())
    error("Cannot find Mixture called %s",name.c_str());
  return GM_Parms.mixtures[i->second];
}



/*-
 *-----------------------------------------------------------------------
 * find_Component
 *      finds a Component in GM_Parms, given its name, with error checking
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      we know that the named Component exists
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a pointer to the Component
 *
 *-----------------------------------------------------------------------
 */
Component* 
find_Component(std::string &name)
{
  GMParms::ObjectMapType::iterator i=GM_Parms.componentsMap.find(name);
  if (i == GM_Parms.componentsMap.end())
    error("Cannot find Component called %s",name.c_str());
  return GM_Parms.components[i->second];
}




/*-
 *-----------------------------------------------------------------------
 * find_Dense1DPMF
 *      finds a Dense1DPMF in GM_Parms, given its name, with error checking
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      we know that the named Dense1DPMF exists
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a pointer to the Dense1DPMF
 *
 *-----------------------------------------------------------------------
 */
Dense1DPMF* 
find_Dense1DPMF(std::string &name)
{
  GMParms::ObjectMapType::iterator i=GM_Parms.dPmfsMap.find(name);
  if (i == GM_Parms.dPmfsMap.end())
    error("Cannot find Dense1DPMF called %s",name.c_str());
  return GM_Parms.dPmfs[i->second];
}


/*-
 *-----------------------------------------------------------------------
 * find_MeanVector
 *      finds a MeanVector in GM_Parms, given its name, with error checking
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      we know that the named MeanVector exists
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a pointer to the MeanVector
 *
 *-----------------------------------------------------------------------
 */
MeanVector* 
find_MeanVector(std::string &name)
{
  GMParms::ObjectMapType::iterator i=GM_Parms.meansMap.find(name);
  if (i == GM_Parms.meansMap.end())
    error("Cannot find MeanVector called %s",name.c_str());
  return GM_Parms.means[i->second];
}



/*-
 *-----------------------------------------------------------------------
 * find_DiagGaussian
 *      finds a DiagGaussian in GM_Parms, given its name, with error checking
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      we know that the named DiagGaussian exists
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a pointer to the DiagGaussian
 *
 *-----------------------------------------------------------------------
 */
DiagGaussian*
find_DiagGaussian(std::string &name)
{
  // these are Components, so look in that map
  GMParms::ObjectMapType::iterator i=GM_Parms.componentsMap.find(name);
  if (i == GM_Parms.componentsMap.end())
    error("Cannot find DiagGaussian called %s",name.c_str());

  // check the type
  Component* component=GM_Parms.components[i->second];
  if( !is_DiagGaussian(component) )
    error("Found a Component called %s but it is not of type DiagGaussian");

  return (DiagGaussian*)GM_Parms.components[i->second];

}



/*-
 *-----------------------------------------------------------------------
 * all_DiagGaussian
 *      checks that all the Components of a Mixture are DiagGaussian
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      we know that all the Components of a Mixture are DiagGaussian
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      bool
 *
 *-----------------------------------------------------------------------
 */
bool 
all_DiagGaussian(Mixture* mixture)
{
  for(vector < Component* >::iterator i=mixture->components.begin(); i!=mixture->components.end(); i++)
    if( !is_DiagGaussian(*i) )
      return false;
  return true;
}


/*-
 *-----------------------------------------------------------------------
 * find_MeanVector_of_Mixture
 *      finds the only MeanVector of a Mixture
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      we know that Mixture only has one Component and that Component is Gaussian
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a MeanVector pointer 
 *
 *-----------------------------------------------------------------------
 */
MeanVector* 
find_MeanVector_of_Mixture(std::string &name)
{
  return find_MeanVector_of_Mixture(find_Mixture(name));
}


MeanVector* 
find_MeanVector_of_Mixture(Mixture *mixture)
{
  if (mixture->components.size() != 1)
    error("find_MeanVector_of_Mixture expected only one Component in this Mixture");

  return find_MeanVector_of_Component(mixture->components[0]);
}

/*-
 *-----------------------------------------------------------------------
 * find_Component_of_Mixture
 *      finds the only Component of a Mixture
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      we know that Mixture only has one Component
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a MeanVector pointer 
 *
 *-----------------------------------------------------------------------
 */
Component* 
find_Component_of_Mixture(std::string &name)
{
  return find_Component_of_Mixture(find_Mixture(name));
}


Component* 
find_Component_of_Mixture(Mixture *mixture)
{
  if (mixture->components.size() != 1)
    error("find_Component_of_Mixture expected only one Component in this Mixture");

  return mixture->components[0];
}


/*-
 *-----------------------------------------------------------------------
 * find_MeanVectors_of_Mixture
 *      finds all the MeanVectors of a Mixture
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      we know that Mixture only contains Gaussian Components
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a vector of MeanVector pointers
 *
 *-----------------------------------------------------------------------
 */
std::vector<MeanVector*> 
find_MeanVectors_of_Mixture(std::string &name)
{
  return find_MeanVectors_of_Mixture(find_Mixture(name));
}

std::vector<MeanVector*> 
find_MeanVectors_of_Mixture(Mixture *mixture)
{
  std::vector<MeanVector*> rval;
  rval.resize(mixture->components.size());
  for (unsigned i=0;i<mixture->components.size();i++)
    rval[i]=find_MeanVector_of_Component(mixture->components[i]);

  return rval;
}



/*-
 *-----------------------------------------------------------------------
 * find_MeanVector_of_Component
 *      finds the MeanVector of a Component
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      we know that the Component is Gaussian (currently in fact
 *      DiagGaussian, but this will be extended in the future)
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a MeanVector pointer
 *
 *-----------------------------------------------------------------------
 */
MeanVector* 
find_MeanVector_of_Component(std::string &name)
{
  return find_MeanVector_of_Component(find_Component(name));
}

MeanVector* 
find_MeanVector_of_Component(Component *component)
{
  if( !is_DiagGaussian(component) )
    error("find_MeanVector_of_Component expected a DiagGaussian Component but got %s",get_Component_type(component).c_str());

  return find_MeanVector_of_DiagGaussian((DiagGaussian*)component);
}

/*-
 *-----------------------------------------------------------------------
 * find_MeanVector_of_DiagGaussian
 *      finds the MeanVector of a DiagGaussian
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      for version taking a string arg, we know that the named
 *      DiagGaussian exists
 *
 * Side Effects:
 *      none
 *
 * Results:
 *       a pointer to a MeanVector
 *
 *-----------------------------------------------------------------------
 */
MeanVector* 
find_MeanVector_of_DiagGaussian(std::string &name)
{
  return find_MeanVector_of_DiagGaussian(find_DiagGaussian(name));
}

MeanVector* 
find_MeanVector_of_DiagGaussian(DiagGaussian *diag_gaussian)
{
  return diag_gaussian->mean;
}



/*-
 *-----------------------------------------------------------------------
 * find_Components_using_MeanVector
 *      finds the number of Components that are using (i.e. sharing)
 *      this MeanVector
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      unsigned
 *
 *-----------------------------------------------------------------------
 */
unsigned
find_Components_using_MeanVector(const MeanVector* const mean, const std::string mean_name, 
				      std::list<Component*> *components,
				      std::list<std::string> *component_names)
{
  unsigned n=0;

  // objects that might use a MeanVector are all sub-types of Component:
  // GMTK_DiagGaussian
  // GMTK_LinMeanCondDiagGaussian (not yet supported)

  // look at all the Components of valid types and return a list of those that use this partic

  for(GMParms::ObjectMapType::iterator i=GM_Parms.componentsMap.begin();i!=GM_Parms.componentsMap.end();i++){
    if( find_MeanVector_of_Component(GM_Parms.components[i->second]) == mean){
      components->push_front(GM_Parms.components[i->second]);
      component_names->push_front(i->first);
      n++;
    }
  }
  return n;
}




/*-
 *-----------------------------------------------------------------------
 * compile_regex
 *      makes a compiled posix regex from its string representation
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      we know that the given string is a valid regex
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      a pointer to a newly-created regex_t object, which the user
 *      will be responsible for deleting
 *
 *-----------------------------------------------------------------------
 */
regex_t* 
compile_regex(const std::string &str)
{
  regex_t *compiled = new regex_t;

  int errcode = regcomp(compiled, str.c_str(), REG_EXTENDED); 
  if (errcode != 0)
    error("regex compilation error for '%s': %s",str.c_str(),get_regerror(errcode,compiled).c_str());

  return compiled;
}


/*-
 *-----------------------------------------------------------------------
 * get_regerror
 *      a wrapper for the posix regerror function
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
 *      the error message as a string
 *
 *-----------------------------------------------------------------------
 */
std::string 
get_regerror (int errcode, regex_t *compiled)
{       

  // to do: assert compiled is valid

  size_t length = regerror (errcode, compiled, NULL, 0);
  char *buffer = (char*)malloc((size_t)(length*sizeof(char)));
  (void)regerror(errcode, compiled, buffer, length);

  //cerr << "REGERR: " << std::string(buffer) << endl;

  std::string str(buffer);
  free((void*)buffer);
  return str;
}

/*-
 *-----------------------------------------------------------------------
 * match
 *      a wrapper for the posix regexec functions
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
 *      true if compiled_regex matches str, false otherwise
 *
 *-----------------------------------------------------------------------
 */
bool 
match(regex_t *compiled_regex, const std::string &str)
{
  // a clean interface to an ugly Posix regex function

  size_t nmatch=MAX_REGEX_GROUPS;
  regmatch_t matchptr[MAX_REGEX_GROUPS];

  int errcode=regexec(compiled_regex,str.c_str(), nmatch, matchptr, 0);
  if (errcode != REG_NOMATCH){
    // make sure the whole string matched
    if ( (matchptr[0].rm_so == 0) and ((unsigned)(matchptr[0].rm_eo) == str.size() ) ){
      return true;
    }      
  }
  return false;
}



/*-
 *-----------------------------------------------------------------------
 * merged_cluster_size
 *      computes how big (according to the supplied ClusterSizeMethod)
 *      the resulting cluster would be if two clusters were merged;
 *      the clusters are NOT merged here though;
 *
 * Preconditions:
 *      the clusters must have been created and must contain valid
 *      Clusterable objects (all of the same type - since this is not
 *      currently checked) with appropriate dissimilarity measures set
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      float
 *
 *-----------------------------------------------------------------------
 */
float 
merged_cluster_size(std::list<Cluster>::iterator ci, std::list<Cluster>::iterator ci2,
		    GMTK_Tie::ClusterSizeMethodType ClusterSizeMethod,
		    GMTK_Tie::CentroidType Centroid)
{

  // assumes the two incoming clusters have at least one item each

  switch(ClusterSizeMethod){

  case GMTK_Tie::CSM_MostDissimilarPair:
  case GMTK_Tie::CSM_EmulateHTK:
    {
      // find farthest apart pair of items in pooled means from both
      // clusters (i.e. the two items that are farthest apart might both
      // be already in the same input cluster; if that *is* the case, we
      // already know the size: it's the size of that cluster
      
      float cross_size=0.0,s;
      std::list<Clusterable*>::iterator ii,ii2;
      for(ii=ci->items.begin();ii!=ci->items.end();ii++)
	for(ii2=ci2->items.begin();ii2!=ci2->items.end();ii2++){
	  
	  // this computes the distance, including symmetric versions
	  // (i.e. no need to also call ii2->dissimilarity(ii);
	  s=(*ii)->dissimilarity(*ii2);
	  
	  //cerr << "merged cross size: "<< s << endl;
	  if ( s > cross_size) cross_size=s;
	}
      
      //cerr << "merged sizes: "<<ci->size<<" "<<ci2->size<<" "<<cross_size<<endl;
      
      return max(max(ci->size,ci2->size),cross_size);
    }
    break;

  case GMTK_Tie::CSM_AverageDissimilarityToCentroid:
    {
      
      // have to actually compute the centroid first - this may get
      // computationally expensive, especially with this (lazy) method
      // involving copying...
      Cluster new_cluster;
      new_cluster.items.insert(new_cluster.items.begin(),ci->items.begin(),ci->items.end());
      new_cluster.items.insert(new_cluster.items.begin(),ci2->items.begin(),ci2->items.end());
      new_cluster.occupancy = ci->occupancy + ci2->occupancy;
      
      set_cluster_centroid(&new_cluster, Centroid, true);

      // now compute average dissimilarity of first item (centroid) to
      // all others
      float s=0;
      std::list<Clusterable*>::iterator ii=new_cluster.items.begin();
      ii++; // start with second item
      for(;ii!=new_cluster.items.end();ii++)
	s += (*(new_cluster.items.begin()))->dissimilarity(*ii);
      return s / (float)(new_cluster.items.size() -1);
	  

    }
    break;
    
  default:
    error("Invalid or unknown ClusterSizeMethod in merged_cluster_size - this is an internal error");

  }

  return -1.0;
}




/*-
 *-----------------------------------------------------------------------
 * build_merged_sizes_table
 *      computes all entries of a table containing pairwise merged
 *      cluster sizes
 *
 * Preconditions:
 *      1) the clusters must have been initialised
 *      2) the memory for merged_sizes must have been allocated
 *
 * Postconditions:
 *      table entries are filled in
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void 
build_merged_sizes_table(sArray< sArray<float> > *merged_sizes,
			 std::list<Cluster> *clusters,
			 GMTK_Tie::ClusterSizeMethodType ClusterSizeMethod,
			 GMTK_Tie::CentroidType Centroid)
{
  std::list<Cluster>::iterator ci,ci2;

  for(ci=clusters->begin();ci!=clusters->end();ci++){
    for(ci2=ci;ci2!=clusters->end();ci2++){
      if(ci==ci2) continue;
      (*merged_sizes)[ci->ident][ci2->ident]=merged_cluster_size(ci,ci2,ClusterSizeMethod,Centroid);
      (*merged_sizes)[ci2->ident][ci->ident]=(*merged_sizes)[ci->ident][ci2->ident];

    }
  }

}

/*-
 *-----------------------------------------------------------------------
 * update_merged_sizes_table
 *      recomputes one row and one column of the merged cluster sizes
 *      table
 *
 * Preconditions:
 *      1) the table must have been built with build_merged_sizes_table
 *      2) ci must point at the newly-merged cluster
 *
 * Postconditions:
 *      table is updated
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void 
update_merged_sizes_table(sArray< sArray<float> > *merged_sizes,
			  std::list<Cluster> *clusters,
			  GMTK_Tie::ClusterSizeMethodType ClusterSizeMethod,
			  GMTK_Tie::CentroidType Centroid,
			  std::list<Cluster>::iterator ci)
{
  // update the row and the column of the table containing entries
  // relating to cluster ci
  std::list<Cluster>::iterator ci2;
  for(ci2=clusters->begin();ci2!=clusters->end();ci2++){
    if(ci==ci2) continue;
    (*merged_sizes)[ci->ident][ci2->ident]=merged_cluster_size(ci,ci2,ClusterSizeMethod,Centroid);
    (*merged_sizes)[ci2->ident][ci->ident]=(*merged_sizes)[ci->ident][ci2->ident];
  }
}


/*-
 *-----------------------------------------------------------------------
 * set_cluster_centroid
 *      for those CentroidType methods which pick a centroid from
 *      amongst the existing members of a cluster, this function moves
 *      that member to the head of the list of items
 *
 * Preconditions:
 *      cluster must be valid
 *
 * Postconditions:
 *      centroid is now at head of list
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      true if centroid could be computed, false otherwise
 *      (indicating a bug, most likely)
 *
 *-----------------------------------------------------------------------
 */
bool 
set_cluster_centroid(Cluster *c, GMTK_Tie::CentroidType method, bool force)
{
  float min_dist=LZERO;

  // put the cluster centroid item at the beginning of the list of
  // items in this cluster (only applies to some CentroidType methods;
  // for others a more complicated calculation must be made later)

  // if force==true, always compute centroid using
  // CNT_UseExistingCentroid method, regardless of method specified

  // problem: here we only ever compute the "median" - should we use
  // that as the basis for inter-cluster distance calculations, even
  // when the users has specified, for example,
  // Centroid=CNT_CreateCentroid_averageSingleComponentMixtures)

  switch(method){
  case GMTK_Tie::CNT_Unknown:
    error("Unknown Centroid type in set_cluster_centroid");
    break;

  case GMTK_Tie::CNT_Arbitrary:
  case GMTK_Tie::CNT_CreateCentroid_averageSingleComponentMixtures:
  case GMTK_Tie::CNT_CreateCentroid_permuteThenAverage:
  case GMTK_Tie::CNT_CreateCentroid_smartMerge:
  case GMTK_Tie::CNT_CreateCentroid_averageMeanVector:
  case GMTK_Tie::CNT_EmulateHTKOther:
    if(!force)
      return true;
    // don't break here!

  case GMTK_Tie::CNT_UseExistingCentroid:
    {
      // compute "median" and re-order the items
      
      // cannot take median of two things so do nothing
      if(c->items.size() <= 2)
	return true;
      
      unsigned j,k;
      unsigned n=c->items.size();
      std::list<Clusterable*>::iterator ci,ci2;
      // build a (symmetric) table of pairwise distances
      std::vector< std::vector<float> > d;
      d.resize(n);
      for(j=0;j<n;j++){
	d[j].resize(n+1);  // the last entry (index n) holds the sum of that row
      }
      
      for(ci=c->items.begin(),j=0;ci!=c->items.end();ci++,j++){
	for(ci2=ci,k=j;ci2!=c->items.end();ci2++,k++){
	  
	  if(j==k) // on the diagonal
	    d[j][k]=0.0; // by definition
	  else{
	    d[j][k]=(*ci)->dissimilarity(*ci2);
	    d[k][j]=d[j][k];

	  }
	  //cerr << j << "," << k << "(" << d[j][k] << ") ";
	}
	//cerr << endl;
      }
      
      for(j=0;j<n;j++){
	d[j][n]=0.0;
	for(k=0;k<n;k++){
	  d[j][n] +=d[j][k];
	  //cerr << d[j][k] << " ";
	}
	//cerr << endl;
      }
      
      // find the centroid - just the min over j of d[j][n] 
      std::list<Clusterable*>::iterator current_centroid_index=c->items.begin();
      for(ci=c->items.begin(),j=0;ci!=c->items.end();ci++,j++){
	if ( (d[j][n] < min_dist) or (min_dist <= LZERO) ){
	  min_dist=d[j][n];
	  current_centroid_index=ci;
	}
      }
      
      // re-order the vector: swap current first position name with
      // centroid
      Clusterable* tmp=*current_centroid_index;
      c->items.erase(current_centroid_index);
      c->items.push_back(tmp);
      
    }
    break;
    
  case GMTK_Tie::CNT_EmulateHTKMixturesOfGaussians:
    {
      error("Centroid method EmulateHTKMixturesOfGaussians not yet implemented");
      break;
    }

  case GMTK_Tie::CNT_EmulateHTK:
    {
      error("Centroid method EmulateHTK should have been mapped internally to one of EmulateHTKMixturesOfGaussians or EmulateHTKOther - this is a bug!");
      break;
    }
    
    
  }

  return (min_dist > LZERO);

}








/*-
 *-----------------------------------------------------------------------
 * cluster_log_likelihood
 *      calculate log liklihood of data being generated by a clustered
 *      Gaussian using the method from Odell's thesis (equations 3.7
 *      and 3.8)
 *
 * Preconditions:
 *      1. clusters must be valid and have occupancy stats
 *      2. parameter type must be single Gaussians with diag covar
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      sets *tot_Occupancy to total occupancy if tot_Occupancy is not null
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */


double
cluster_scaled_log_likelihood(std::list<Clusterable*> &items, double* tot_Occupancy)
{
  infoMsg(IM::Max,"cluster_scaled_log_likelihood started\n");

  // method used is from J. Odell's thesis "THE USE OF CONTEXT IN
  // LARGE VOCABULARY SPEECH RECOGNITION", Cambridge, 1995

  // equation 3.8 approximates the variance of a Gaussian that fits to
  // all the data frames aligned with the Guassians being clustered

  // Sigma = E[o^2] - E[o]^2
  //            term1              (      term2     )^2
  //       =  ---------------   -  (----------------)
  //          total_occupancy      ( total_occupancy)
  

  // Arthur Kantor 12/10/08: I have put in the constants so the function now 
  // returns the actual log-likelihood approximation, exactly as in equation
  // 3.7 of Odell's thesis.  Perhaps this function should be renamed to cluster_log_likelihood.
  // The reason is that without the constants, log-likelihood can be positive, and when this 
  // function is used to compute relative improvement (Percent log likelihood improvement)
  // an improvement appears as worsening.
  // The following is no longer true.
  // NOT_TRUE:returns a scaled log likelihood because constant terms are not
  // NOT_TRUE:included - the return value is intended only for computing the
  // NOT_TRUE:difference in cluster_scaled_log_likelihood values between
  // NOT_TRUE:clusters
  // 

  // maybe better to work in double precision throughout this
  // function?
  
  if ( items.empty() )
    error("Cannot compute cluster_scaled_log_likelihood for an empty cluster");

  infoMsg(IM::Max,"cluster_scaled_log_likelihood for a cluster of %d items\n",items.size());

  ClusterableDiagGaussian *cdg = (ClusterableDiagGaussian*)(*(items.begin()));
  unsigned dim=cdg->param_ptr->mean->means.len();


  sArray<double> term1, term2;
  double temp;
  term1.resizeAndZero(dim);
  term2.resizeAndZero(dim);
  temp=0.0;
  double total_occupancy=0.0;

  //cerr << "starting loop" << endl;

  //infoMsg(IM::Max,"cluster_scaled_log_likelihood starting loop\n");

  for(std::list<Clusterable*>::iterator ci=items.begin();ci!=items.end();ci++){

    //cerr << " doing "<< (*ci)->name() << endl;

    ClusterableDiagGaussian* cim=(ClusterableDiagGaussian*)(*ci);
    sArray<float> *this_mean=&(cim->param_ptr->mean->means);
    sArray<float> *this_covar=&(cim->param_ptr->covar->covariances);

    total_occupancy += (*ci)->occupancy();

    infoMsg(IM::Max,"cluster_scaled_log_likelihood: processing item %s with occupancy %f\n",
	    (*ci)->name().c_str(),(*ci)->occupancy());


    for(unsigned i=0;i<dim;i++){
      temp = (double)((*this_mean)[i]);

      term2[i] += temp * (*ci)->occupancy();

      temp = temp * temp;
      temp += (double)((*this_covar)[i]);
      temp *= (double)((*ci)->occupancy());

      term1[i] += temp;
    }

  }

  
  infoMsg(IM::Max,"cluster_scaled_log_likelihood: total occupancy is %f\n",total_occupancy);
  if (tot_Occupancy){
    *tot_Occupancy=total_occupancy;
  }
  if(total_occupancy==0.0){
    infoMsg(IM::Max,"Total occupancy for a cluster is zero");
    return LZERO;

  }

  for(unsigned i=0;i<dim;i++){
    term1[i] = term1[i] / total_occupancy;

    term2[i] = term2[i] / total_occupancy;
    term2[i] = term2[i] * term2[i];
  }

  //cerr << "total_occupancy="<<total_occupancy<<endl;

  term1 -= term2; // now term1 is the desired covar

  // determinant of diagonal covar matrix
  logpr logdet;
  logdet.set_to_one();

  for(unsigned i=0;i<dim ;i++){
    //cerr << "covar = " << term1[i];
    //cerr << " log covar = " << log(term1[i]) << endl;

    // operators for logpr are defined as if we are working with
    // probs, not log probs, so (very confusingly) we do *= here and
    // not +=
    logdet *= logpr(term1[i]);
  }
  

  // equation 3.7 approximates the log likelihood of generating those
  // aligned frames using the new single Gaussian for this cluster

  // 

  //cerr << "logdet="<<logdet.val() << " occ=" << total_occupancy <<endl;
  //cerr << "logl="<< logdet.val() * total_occupancy <<endl;
  const double LOG_2PI = 1.83787706640934548356;
  return -0.5 * (dim * (1 + LOG_2PI) + logdet.val()) * total_occupancy;
}






/*-
 *-----------------------------------------------------------------------
 * new_name
 *      devises a unique new name for a RV
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      none 
 *
 * Side Effects:
 *      none (e.g., GM_Parms is unchanged)
 *
 * Results:
 *      the new name as a string
 *
 *-----------------------------------------------------------------------
 */
std::string
new_name(const std::string basename, map< string, unsigned > *map)
{
  std::string name;

  static unsigned cloneNo; 

  do {
    char buff[256];
    sprintf(buff,"%d",cloneNo);
    name = basename + string("_") + buff;
    cloneNo++;
  } while (map->find(name) != map->end());

  return name;
}






/*-
 *-----------------------------------------------------------------------
 * new_name
 *      devises a unique new name for a RV
 *
 * Preconditions:
 *      GM_Parms must exist
 *
 * Postconditions:
 *      none 
 *
 * Side Effects:
 *      sort order of sorted_table may become invalid, so the name
 *      collection is unsorted
 *
 * Results:
 *      nc->table is modified
 *
 *-----------------------------------------------------------------------
 */

/*
void search_and_replace_in_name_collection(NameCollection *nc,std::vector<std::string> old_names, std::string new_name)
{

  //cerr << "search_and_replace_in_name_collection for this many old names:" << old_names.size() << endl;
  infoMsg(IM::Huge,"search_and_replace_in_name_collection\n");
  unsigned nchanged=0;

  //cerr << "nc name is " << nc->name() << endl;

  if(! nc->is_sorted())
    nc->sort();

  //cerr << "sorting" << endl;
  std::sort(old_names.begin(),old_names.end());

  vector<std::string>::iterator nci,i,ncie,ie;

  i=old_names.begin();
  ie=old_names.end();

  nci=nc->sorted_table.begin();
  ncie=nc->sorted_table.end();

  //cerr << "running" << endl; 

  // jump to first match
  nci=std::find(nci,ncie,*i);

  while( (i!=ie) and (nci!=ncie) ){
  
    while ( (nci!=ncie) and (*nci == *i) ){
      //cerr << "replacing " << *nci << " with " << new_name << endl;
      *nci = new_name;
      nchanged++;
      nci++;
    }

    if ( (nci!=ncie)  and (*nci < *i) )
      nci=std::find(nci,ncie,*i);
    
    if (nci==ncie)
      break;
    
    if ( (i!=ie) and (*i < *nci) )
      i++;
    
  }

  nc->resort();

  infoMsg(IM::Huge,"search_and_replace_in_name_collection changed %d entries\n",nchanged);


}
*/
