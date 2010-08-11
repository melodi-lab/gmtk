/* 
 * GMTK_Tie.cc
 * 
 * This class contains the main functionality for clustering, tying
 * and untying parameters
 * 
 * 
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

#include "GMTK_Tie.h"
#include "tieSupport.h"
#include <limits>
#include <algorithm>
VCID("$Header$")



std::list<Cluster> null_list;
std::list<Cluster>::iterator NULL_ITERATOR=null_list.begin();
// temporary
#include <cstdlib>

// max line length seems a bit small?
//Yep.  -arthur
//#define MAX_LINE_LENGTH 1024
//#define MAX_FEATURE_VALUES 100
#define MAX_LINE_LENGTH 102400
#define MAX_FEATURE_VALUES 5000


// this will be divided by at some point, so must not result in
// overflow
// value is a wild guess!?
#define MIN_OCCUPANCY std::numeric_limits<float>::min() * 100



void print_map(std::map<std::string,unsigned> &map)
{
  cerr << "Map with " << map.size() << " entries: " << endl;
  for(std::map<std::string,unsigned>::iterator it=map.begin();it!=map.end();it++)
    cerr << it->first << " -> " << it->second << endl;
}


void print_map(std::map<unsigned,std::string> &map)
{
  cerr << "Map with " << map.size() << " entries: " << endl;
  for(std::map<unsigned,std::string>::iterator it=map.begin();it!=map.end();it++)
    cerr << it->first << " -> " << it->second << endl;
}





// commit changes to a name collection
void commit_nc_changes(string collection_name){
  if (collection_name =="")
    return;

  GMParms::ObjectMapType::iterator i=GM_Parms.nclsMap.find(collection_name);
  if (i == GM_Parms.nclsMap.end())
    error("Cannot find named collection called %s",collection_name.c_str());
  
  NameCollection* nc=GM_Parms.ncls[i->second];
  nc->commit_all_searches_and_replacements();

}



/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::GMTK_Tie
 *      constructor
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      1) GMTK_Tie object is constructed
 *      2) the maps between strings and enums are constructed
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *      
 *-----------------------------------------------------------------------
 */
GMTK_Tie::GMTK_Tie(GMParms *GM_Parms_ptr, char* cppCommandOptions)
{
  GM_Parms=GM_Parms_ptr;
  _cppCommandOptions=cppCommandOptions;
  ////////////////////////////////////////////////////////////////////////
  // make maps between tying command file keywords and the enumerated
  // types

  CommandType_to_string[CT_Unknown]="unknown";
  CommandType_to_string[CT_Tie]="tie";
  CommandType_to_string[CT_Cluster]="cluster";
  CommandType_to_string[CT_DTcluster]="DTcluster";
  CommandType_to_string[CT_Untie]="untie";
  CommandType_to_string[CT_DTsynthesise]="DTsynthesise";
  CommandType_to_string[CT_loadFeatureDefinitions]="loadFeatureDefinitionSet";
  CommandType_to_string[CT_loadQuestions]="loadQuestionSet";
  CommandType_to_string[CT_loadFeatureValues]="loadFeatureValueSet";
  CommandType_to_string[CT_saveTree]="saveTree";
  CommandType_to_string[CT_loadTree]="loadTree";

  ParameterType_to_string[PT_Unknown]="Unknown";
  ParameterType_to_string[PT_Mixture]="Mixture";
  ParameterType_to_string[PT_Component]="Component";
  ParameterType_to_string[PT_Mean]="Mean";
  ParameterType_to_string[PT_Collection]="Collection";
  ParameterType_to_string[PT_NotRequired]="NotRequired";

  DissimilarityMeasureType_to_string[DMT_Unknown]="Unknown";
  DissimilarityMeasureType_to_string[DMT_Euclidean]="Euclidean";
  DissimilarityMeasureType_to_string[DMT_AverageEuclidean]="AverageEuclidean";
  DissimilarityMeasureType_to_string[DMT_Mahalanobis]="Mahalanobis";
  DissimilarityMeasureType_to_string[DMT_SymmetricMahalanobis]="SymmetricMahalanobis";
  DissimilarityMeasureType_to_string[DMT_KullbackLeiblerDivergence]="KullbackLeiblerDivergence";
  DissimilarityMeasureType_to_string[DMT_SymmetricKullbackLeiblerDivergence]="SymmetricKullbackLeiblerDivergence";
  DissimilarityMeasureType_to_string[DMT_CrossLogLikelihoodOfMeans]="CrossLogLikelihoodOfMeans";
  DissimilarityMeasureType_to_string[DMT_EmulateHTK]="EmulateHTK";

  CentroidType_to_string[CNT_Unknown]="Unknown";
  CentroidType_to_string[CNT_Arbitrary]="Arbitrary";
  CentroidType_to_string[CNT_UseExistingCentroid]="UseExistingCentroid";
  CentroidType_to_string[CNT_CreateCentroid_averageSingleComponentMixtures]="CreateCentroid_averageSingleComponentMixtures";
  CentroidType_to_string[CNT_CreateCentroid_permuteThenAverage]="CreateCentroid_permuteThenAverage";
  CentroidType_to_string[CNT_CreateCentroid_smartMerge]="CreateCentroid_smartMerge";
  CentroidType_to_string[CNT_CreateCentroid_averageMeanVector]="CreateCentroid_averageMeanVector";
  CentroidType_to_string[CNT_EmulateHTK]="EmulateHTK";

  // non-user values
  CentroidType_to_string[CNT_EmulateHTKMixturesOfGaussians]="EmulateHTKMixturesOfGaussians";
  CentroidType_to_string[CNT_EmulateHTKOther]="EmulateHTKOther";

  ClusterSizeMethodType_to_string[CSM_Unknown]="Unknown";
  ClusterSizeMethodType_to_string[CSM_MostDissimilarPair]="MostDissimilarPair";
  ClusterSizeMethodType_to_string[CSM_EmulateHTK]="EmulateHTK";
  ClusterSizeMethodType_to_string[CSM_AverageDissimilarityToCentroid]="AverageDissimilarityToCentroid";

  // the reverse maps
  for(map<ParameterType, std::string>::iterator i=ParameterType_to_string.begin();
      i!= ParameterType_to_string.end(); i++)
    string_to_ParameterType[i->second] = i->first;

  for(map<CommandType, std::string>::iterator i=CommandType_to_string.begin();
      i!= CommandType_to_string.end(); i++)
    string_to_CommandType[i->second] = i->first;

  for(map<DissimilarityMeasureType, std::string>::iterator i=DissimilarityMeasureType_to_string.begin();
      i!= DissimilarityMeasureType_to_string.end(); i++)
    string_to_DissimilarityMeasureType[i->second] = i->first;

  for(map<CentroidType, std::string>::iterator i=CentroidType_to_string.begin();
      i!= CentroidType_to_string.end(); i++)
    string_to_CentroidType[i->second] = i->first;

  for(map<ClusterSizeMethodType, std::string>::iterator i=ClusterSizeMethodType_to_string.begin();
      i!= ClusterSizeMethodType_to_string.end(); i++)
    string_to_ClusterSizeMethodType[i->second] = i->first;

}





/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::~GMTK_Tie() 
 *      destructor
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      GMTK_Tie object is deleted
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
GMTK_Tie::~GMTK_Tie() 
{ 
  for (CommandListType::iterator i=commands.begin();i!=commands.end();i++)
    regfree(i->param_expression);

}






/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::read_commands(iDataStreamFile& is)
 *      loads a tying command file
 * 
 * Preconditions:
 *      s must be a valid iDataStreamFile
 *
 * Postconditions:
 *      'commands' is filled in 
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      true if file is successfully read and parsed
 *
 *-----------------------------------------------------------------------
 */
bool
GMTK_Tie::read_commands(iDataStreamFile& is)
{

  // stick with GMTK-style file format for now, viz:
  //  <int, giving number of commands>
  // then that many commands, each of which looks like:
  //  <int, giving command index>
  //  <string, givingcommand> <string, giving parameter type> <string, giving param name regex>
  //  <int, giving num options>
  //  <option> 
  //  <option>
  //  ...
  //
  // where each <option> looks like
  // <string, giving option name>=<string, giving option value>, e.g. MaxClusterSize=20
  //
  
  /* an example command
     0
     cluster Mixture gm1[234]
     7
     MaxClusterSize=4e-27
     MinClusterMembers=3
     MinOccupancyCount=1000.0
     DissimilarityMeasure=CrossLikelihoodOfMeans
     ClusterSizeMethod=AverageDissimilarityToCentroid
     NewParameterNamePrefix=tied_mean
     Centroid=UseExistingCentroid
  */
  unsigned ncmds,index,noptions;

  if(!is.readUnsigned(ncmds)) return false;
  commands.resize(ncmds);

  for (unsigned i=0;i<ncmds;i++){
    if (!is.readUnsigned(index)) return false;

    if (i != index)
      error("Command index is out of order");

    infoMsg(IM::High,"loading command %d\n", i);


    // the command name
    char* cmd;
    if (!is.readStr(cmd)) return false;
    if( string_to_CommandType.find(cmd) != string_to_CommandType.end() )
      commands[i].command=string_to_CommandType[cmd];
    else
      error("Unknown command type %s in command %d",cmd,i);

    // the parameter on which it operates - only applies to some commands
    switch(commands[i].command)
      {
      case CT_Tie:
      case CT_Cluster:
      case CT_DTcluster:
      case CT_Untie:
      case CT_DTsynthesise:
	{
	  char* param_type;
	  if (!is.readStr(param_type)) return false;
	  if( string_to_ParameterType.find(param_type) != string_to_ParameterType.end() )
	    commands[i].param_type=string_to_ParameterType[param_type];
	  else
	    error("Unknown parameter type %s in command %d",param_type,i);

	  // the regular expression specifying the parameter names
	  char* param_expression;
	  if (!is.readStr(param_expression)) return false;
	  commands[i].param_expression=compile_regex(param_expression);
	}
	break;

      default: 
	commands[i].param_type=PT_NotRequired;
	commands[i].param_expression=NULL;
	break;
      }
    


    // the list of optional arguments, each of the form key=value
    if (!is.readUnsigned(noptions)) return false;

    fill_in_default_values(commands[i].command, commands[i].options);

    for (unsigned j=0;j<noptions;j++){

      char* opt;
      if (!is.readStr(opt)) return false;
      parse_option(commands[i].command, commands[i].options, std::string(opt));

    }

    // now check whether the options make sense for this command 
    /*
      switch(ctype){
      case CT_Tie:
      enter_tie_option(o, key, value);
      break;
      case CT_Cluster:
      enter_cluster_option(o, key, value);
      break;
      case CT_Untie:
      error("Command %s does not take any parameters",CommandType_to_string[ctype].c_str());
      break;
      case CT_DTcluster:
      enter_DTcluster_option(o, key, value);
      break;
      default:
      error("Command type %s is not yet supported", CommandType_to_string[ctype].c_str());
      }
    */


    // depending on the parameter type, CentroidType CNT_EmulateHTK
    // needs to be mapped on to either
    // CNT_EmulateHTKMixturesOfGaussians or CNT_EmulateHTKOther

    if ( commands[i].options.Centroid ==  GMTK_Tie::CNT_EmulateHTK ){

      if ( commands[i].param_type == PT_Mixture )
	commands[i].options.Centroid = GMTK_Tie::CNT_EmulateHTKMixturesOfGaussians;
      else
	commands[i].options.Centroid = GMTK_Tie::CNT_EmulateHTKOther;

    }

  }

  return true;
}


/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::need_occupancy_counts
 *      determines whether any of the commands in 'commands' require
 *      occupancy stastics
 * 
 * Preconditions:
 *      'commands' should have been loaded
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      true, if occupancy stats are needed by one or more commands;
 *      false otherwise
 *
 *-----------------------------------------------------------------------
 */
bool
GMTK_Tie::need_occupancy_counts()
{
  for (CommandListType::iterator i=commands.begin();i!=commands.end();i++){
    if ( ! (i->options.MinOccupancyCount == 0.0) ) // 0.0 is special flag value
      return true;
    if ( ! (i->options.ThresholdOccupancyCount == 0.0) ) // 0.0 is special flag value
      return true;
    if ( i->command == CT_DTcluster )
      return true;
    if ( i->options.CentroidOccupancyWeighting or i->options.ClusterSizeOccupancyWeighting)
      return true;
  }

  return false;
}

/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::read_feature_definition_set(iDataStreamFile& is)
 *      loads the features for all parameters to be clustered
 * 
 * Preconditions:
 *      s must be a valid iDataStreamFile
 *
 * Postconditions:
 *      a named item is added to feature_definition_sets
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
GMTK_Tie::read_feature_definition_set(iDataStreamFile& is)
{

  unsigned n,i,index;
  char *line = new char [MAX_LINE_LENGTH];
  std::string s,val;

  regex_t *compiled = compile_regex("[[:space:]]*([^[:space:]]+)[[:space:]]+\\(([[:print:]]+)\\)"); 
  regmatch_t matchptr[MAX_FEATURE_VALUES];

  if(!is.readString(s))
    error("Failed to read name of this set of feature definitions");

  map<std::string,FeatureDefinitionSetType>::iterator it;
  it=feature_definition_sets.find(s);
  if (it != feature_definition_sets.end()){
    warning("Loading a feature definition set called '%s' will overwrite the one already in memory\n",s.c_str());
    feature_definition_sets.erase(s);
  }

  FeatureDefinitionSetType new_fds;
  new_fds.name=s;

  if(!is.readUnsigned(n))
    error("Failed to read number of features in this definition set");

  new_fds.features.resize(n);

  for (i=0;i<n;i++){
    if (!is.readUnsigned(index))
      error("Failed to read next feature definition index");

    if (i != index)
      error("Feature definitions index is out of order");

    infoMsg(IM::High,"loading feature definition %d\n", i);

    if (!is.readLine(line,MAX_LINE_LENGTH))
      error("Failed to read next line");

    s=string(line);
    int errcode=regexec(compiled, line, 3, matchptr, 0);

    if (errcode == REG_NOMATCH)
      error("Failed to parse feature definition");
    else if (errcode == REG_ESPACE)
      error("Regex ran out of memory - this is a bug");
    
    if (matchptr[1].rm_so < 0)
      error("Failed to parse feature definition - missing feature name");

    new_fds.features[i].name=s.substr(matchptr[1].rm_so,matchptr[1].rm_eo-matchptr[1].rm_so);
    new_fds.feature_name_to_index[new_fds.features[i].name]=i;

    infoMsg(IM::High," feature name %s\n", new_fds.features[i].name.c_str());

    //cerr << "feat name: " << feature_names[i] << endl;

    if (matchptr[2].rm_so < 0)
      error("Failed to parse feature definition - missing feature values");

    s = s.substr(matchptr[2].rm_so,matchptr[2].rm_eo-matchptr[2].rm_so);
    //cerr << "feat vals: " << s << endl;

    // now parse feature values - would be nicer to use regex groups,
    // but posix sucks, so we can't do that easily
    unsigned ws_index;
    unsigned j=0;
    while ( (ws_index=s.find_first_of(" ")) != std::string::npos){
      val = s.substr(0,ws_index);
      s = s.substr(ws_index+1,std::string::npos);

      new_fds.features[i].string_values.insert(new_fds.features[i].string_values.end(), val);
      new_fds.features[i].string_value_to_index[val]=j;

      j++;
    }
    val = s;
    new_fds.features[i].string_values.insert(new_fds.features[i].string_values.end(), val);
    new_fds.features[i].string_value_to_index[val]=j;


  }

  feature_definition_sets[new_fds.name]=new_fds;

  regfree(compiled);
  delete [] line;

}


/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::read_feature_value_set(iDataStreamFile& is)
 *      loads the feature values for some parameters to be clustered
 * 
 * Preconditions:
 *      s must be a valid iDataStreamFile
 *
 * Postconditions:
 *      a named item is added to feature_value_sets
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
GMTK_Tie::read_feature_value_set(iDataStreamFile& is)
{
  unsigned n,index,ival;
  char *line = new char [MAX_LINE_LENGTH];
  std::string s,val,pname;

  infoMsg(IM::High,"Loading feature values\n");

  // now load the feature values
  // example line:      0  gmMx0 and 0 B 3 0

  if(!is.readString(s))
    error("Failed to read name of this set of feature values");
  infoMsg(IM::High,"Found a feature value set called '%s'\n",s.c_str());

  map<std::string,FeatureValueSetType>::iterator it;
  it=feature_value_sets.find(s);
  if (it != feature_value_sets.end()){
    warning("Loading a feature value set called '%s' will overwrite the one already in memory\n",s.c_str());
    feature_value_sets.erase(s);
  }

  FeatureValueSetType new_fvs;
  new_fvs.name=s;

  if(!is.readString(s))
    error("Failed to read name of the feature definitions to use for this set of feature values");
  infoMsg(IM::High," which uses feature definition set '%s'\n",s.c_str());

  if (feature_definition_sets.find(s) == feature_definition_sets.end())
    error("There is no feature definition set called '%s' loaded into memory",s.c_str());
  
  new_fvs.feature_definitions=&(feature_definition_sets[s]);

  unsigned num_features=new_fvs.feature_definitions->features.size();

  if(!is.readUnsigned(n))
    error("Failed to read number of feature values (i.e number of parameters)");

  for (unsigned i=0;i<n;i++){
    if (!is.readUnsigned(index))
      error("Failed to read next feature value index");

    if (i != index)
      error("Feature values index is out of order");

    if (!is.readString(s))
      error("Failed to read parameter name");

    pname=s;

    //infoMsg(IM::Max,"loading feature values for parameter index %d called %s\n", i,pname.c_str());


    if (new_fvs.values.find(pname) != new_fvs.values.end())
      error("Repeated feature values provided for parameter %s\n",pname.c_str());

    new_fvs.values[pname].resize(num_features);
    //cerr << " feature " << val;

    for(unsigned j=0;j<num_features;j++){
      if (!is.readStr(line))
	error("Failed to read next feature value");
      val=string(line);

      if(new_fvs.feature_definitions->features[j].string_value_to_index.find(val) == new_fvs.feature_definitions->features[j].string_value_to_index.end()){
	print_map(new_fvs.feature_definitions->features[j].string_value_to_index);
	error("Invalid feature value: %s",line);
      }

      ival=new_fvs.feature_definitions->features[j].string_value_to_index[val];
      new_fvs.values[pname][j]=ival;

    }

  }
  
  feature_value_sets[new_fvs.name]=new_fvs;
  infoMsg(IM::High,"Finished loading feature values\n");

  delete [] line;
}

/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::read_question_set(iDataStreamFile& is)
 *      loads the questions for DT clustering
 * 
 * Preconditions:
 *      s must be a valid iDataStreamFile;
 *      the feature definitions and values used by this question set
 *      should already have been loaded
 *
 * Postconditions:
 *      a named item is added to question_sets
 *
 * Side Effects:
 *      all decision trees in memory are forgotten since they may not
 *      be valid with this question set
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
GMTK_Tie::read_question_set(iDataStreamFile& is)
{
  decision_trees.clear();

  unsigned nqs,index;
  //char *line = new char [MAX_LINE_LENGTH];
  std::string val,val2,s,fname;

  if(!is.readString(s))
    error("Failed to read name of this question set");

  infoMsg(IM::High,"Found a question set called '%s'\n",s.c_str());

  map<std::string,QuestionSetType>::iterator it;
  it=question_sets.find(s);
  if (it != question_sets.end()){
    warning("Loading a question set called '%s' will overwrite the one already in memory\n",s.c_str());
    question_sets.erase(s);
  }

  QuestionSetType new_qs;
  new_qs.name=s;

  if(!is.readString(s))
    error("Failed to read name of the feature definitions to use for this set of questions");
  infoMsg(IM::High," which uses feature definition set '%s'\n",s.c_str());

  if (feature_definition_sets.find(s) == feature_definition_sets.end())
    error("There is no feature definition set called '%s' loaded into memory",s.c_str());
  new_qs.feature_definitions=&(feature_definition_sets[s]);



  if(!is.readString(s))
    error("Failed to read name of the feature value set to use for this set of questions");
  infoMsg(IM::High," and uses feature value set '%s'\n",s.c_str());


  if (feature_value_sets.find(s) == feature_value_sets.end())
    error("There is no feature value set called '%s' loaded into memory",s.c_str());
  new_qs.feature_values=&(feature_value_sets[s]);





  if(!is.readUnsigned(nqs,"number of questions")) 
    error("Failed to read number of questions in this set");
  new_qs.questions.resize(nqs);
  infoMsg(IM::High,"There are %d questions in this set\n",nqs);

  for (unsigned i=0;i<nqs;i++){
    if (!is.readUnsigned(index,"question index")) 
      error("Failed to read question index");

    if (i != index)
      error("Question index is out of order");

    infoMsg(IM::High,"loading question %d\n", i);


    // the question's name
    std::string qname;
    if (!is.readString(qname,"question's name")) 
      error("Failed to read question's name");
    new_qs.questions[i].name=qname;


    // the question's feature name
    std::string fname;
    if (!is.readString(fname,"question's feature name")) 
      error("Failed to read question's feature name");

    if( new_qs.feature_definitions->feature_name_to_index.find(fname) != new_qs.feature_definitions->feature_name_to_index.end() )
      {
	new_qs.questions[i].feature_name=fname;
	new_qs.questions[i].feature_index=new_qs.feature_definitions->feature_name_to_index[fname];
      }
    else
      error("There is no feature called '%s' in the feature definition set '%s' which is in use for this question set (question %d)",fname.c_str(),new_qs.feature_definitions->name.c_str(),i);


    unsigned nvals;
    if (!is.readUnsigned(nvals,"number of feature values"))
      error("Failed to read number of feature values for this question");

    for(unsigned j=0;j<nvals;j++){
      if (!is.readString(val,"feature value"))
	error("Failed to read next feature value");



      if(new_qs.feature_definitions->features[new_qs.questions[i].feature_index].string_value_to_index.find(val)
	 == new_qs.feature_definitions->features[new_qs.questions[i].feature_index].string_value_to_index.end()){
	print_map(new_qs.feature_definitions->features[new_qs.questions[i].feature_index].string_value_to_index);
	error("Invalid value '%s' for feature '%s'",val.c_str(),fname.c_str());
      }
      new_qs.questions[i].valueSet.insert(new_qs.feature_definitions->features[new_qs.questions[i].feature_index].string_value_to_index[val]);
    }
  }

  question_sets[new_qs.name]=new_qs;

  //delete [] line;
}

/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::purge_features_and_questions()
 *      empty data structures that hold features and questions
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *     these members are emptied: 
 *      num_features, feature_names, feature_name_to_index,
 *      string_to_FeatureValueType_maps, FeatureValueType_to_string_maps
 *      features, questions
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

/*
void 
GMTK_Tie::purge_features_and_questions()
{
  feature_definition_sets.clear();
  feature_value_sets.clear();
  question_sets.clear();
}
*/

/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::parse_option
 *      parses a command options that has been read from the command files
 * 
 * Preconditions:
 *      command must have been loaded and its type determined
 *
 * Postconditions:
 *      'o' has the option added to it
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
GMTK_Tie::parse_option(CommandType ctype, Options &o, std::string opt)
{
  // parsing is dependent on ctype - some options do not apply to some
  // types

  // 'opt' should look like "key=value" with no spaces allowed 
  // to do: should change this

  // split opt into key and value - very simple strategy (too much
  // bother to use a proper regex with grouping - will do that later)
  std::string key,value;
  std::string::size_type eq_index;

  if ( (eq_index=opt.find_first_of("=")) == std::string::npos)
    error("Option '%s' does not contain a '='",opt.c_str());

  if (eq_index == opt.length()-1)
    error("Option '%s' has nothing after the '='",opt.c_str());

  if ( opt.find_first_of("=") != (eq_index=opt.find_last_of("=")) )
    error("Option '%s' contains more than one '='",opt.c_str());

  if (eq_index == 0)
    error("Option '%s' has nothing before the '='",opt.c_str());


  key = opt.substr(0,eq_index);
  value = opt.substr(eq_index+1,std::string::npos);

  enter_option(o, key, value);



}





/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::enter_option
 *      parses an option and enter values into the opt structure
 *      
 * 
 * Preconditions:
 *       an option must have been read in and parsed into key and value
 *
 * Postconditions:
 *      opt has the option added to it
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
GMTK_Tie::enter_option(GMTK_Tie::Options &opt, std::string &key, std::string &value)
{

  // the strings are hardwired here - this is a bad idea
  // to do: fix it

  if  (key == "MaxClusterSize"){
    opt.MaxClusterSize=atof(value.c_str());

  } else if  (key == "MaxNumberParameters"){
    opt.MaxNumberParameters=atoi(value.c_str());

  } else if  (key == "MinOccupancyCount"){
    opt.MinOccupancyCount=(double)atof(value.c_str());

  } else if  (key == "ThresholdOccupancyCount"){
    opt.ThresholdOccupancyCount=(double)atof(value.c_str());

  } else if  (key == "MinClusterMembers"){
    opt.MinClusterMembers=atoi(value.c_str());

  } else if  (key == "MinImprovementPercent"){
    opt.MinImprovementPercent=atof(value.c_str());
    if ( (opt.MinImprovementPercent<0.0) or (opt.MinImprovementPercent > 100.0))
      warning("Value specified for MinImprovementPercent is unlikely to work: %f\n",opt.MinImprovementPercent);

  } else if  (key == "MinImprovementAbsolute"){
    opt.MinImprovementAbsolute=(double)atof(value.c_str());

  } else if  (key == "NewParameterNamePrefix"){
    opt.NewParameterNamePrefix=value;

  } else if  (key == "DissimilarityMeasure"){
    map<std::string, DissimilarityMeasureType>::iterator i = string_to_DissimilarityMeasureType.find(value);
    if(i != string_to_DissimilarityMeasureType.end())
      opt.DissimilarityMeasure=i->second;
    else
      error("Invalid value of '%s' for option %s",value.c_str(),key.c_str());

  } else if  (key == "Centroid"){
    map<std::string, CentroidType>::iterator i = string_to_CentroidType.find(value);
    if(i != string_to_CentroidType.end())
      opt.Centroid=i->second;
    else
      error("Invalid value of '%s' for option %s",value.c_str(),key.c_str());

  } else if  (key == "CentroidOccupancyWeighting"){
    if(value == "true" or value == "True" or value == "TRUE")
      opt.CentroidOccupancyWeighting=true;
    else if (value == "false" or value == "False" or value == "FALSE")
      opt.CentroidOccupancyWeighting=true;
    else 
      error("Invalid value of '%s' for option %s",value.c_str(),key.c_str());

  } else if  (key == "ClusterSizeMethod"){
    map<std::string, ClusterSizeMethodType>::iterator i = string_to_ClusterSizeMethodType.find(value);
    if(i != string_to_ClusterSizeMethodType.end())
      opt.ClusterSizeMethod=i->second;
    else
      error("Invalid value of '%s' for option %s",value.c_str(),key.c_str());

  } else if  (key == "Filename")
    opt.Filename=value;
  
  else if (key == "TreeName")
    opt.TreeName=value;
  else if (key == "FeatureSetName")
    opt.FeatureSetName=value;
  else if (key == "QuestionSetName")
    opt.QuestionSetName=value;
  else if (key == "FeatureValuesName")
    opt.FeatureValuesName=value;
  else if (key == "CollectionName")
    opt.CollectionName=value;

  else
    error("Unrecognised option '%s'",key.c_str());

}




/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::fill_in_default_values
 *      fills in default values for some options for a particular command type
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      opt has some options added to it
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
GMTK_Tie::fill_in_default_values(CommandType ctype, GMTK_Tie::Options &opt)
{
  opt.MaxClusterSize=std::numeric_limits<float>::max();
  opt.MaxNumberParameters=std::numeric_limits<unsigned>::max();
  opt.MinOccupancyCount=std::numeric_limits<unsigned>::min();
  opt.ThresholdOccupancyCount=0.0; // ??? I think this is correct
  opt.MinClusterMembers=1;
  opt.MinImprovementPercent=0.0;
  opt.MinImprovementAbsolute=std::numeric_limits<double>::min();
  opt.NewParameterNamePrefix="";
  opt.DissimilarityMeasure=DMT_Unknown;
  opt.Centroid=CNT_Unknown;
  opt.CentroidOccupancyWeighting=false;
  opt.ClusterSizeMethod=CSM_Unknown;
  opt.ClusterSizeOccupancyWeighting=false;
  opt.Filename="";
  opt.TreeName="";
  opt.CollectionName="";

  switch(ctype){
  case CT_Tie:
    opt.NewParameterNamePrefix="Tied_";
    opt.DissimilarityMeasure=DMT_EmulateHTK;
    opt.Centroid=CNT_UseExistingCentroid;
    break;

  case CT_Untie:
    opt.NewParameterNamePrefix="Untied_";
    break;

case CT_Cluster:
    opt.NewParameterNamePrefix="Clustered_";
    opt.DissimilarityMeasure=DMT_EmulateHTK;
    opt.Centroid=CNT_UseExistingCentroid;
    opt.ClusterSizeMethod=CSM_EmulateHTK;
    break;

  case CT_DTcluster:
    opt.NewParameterNamePrefix="DT_clustered_";
    opt.DissimilarityMeasure=DMT_EmulateHTK;
    opt.Centroid=CNT_CreateCentroid_averageSingleComponentMixtures; // check if this matches HTK method
    opt.CentroidOccupancyWeighting=true;
    opt.ClusterSizeMethod=CSM_EmulateHTK;
    opt.ClusterSizeOccupancyWeighting=true;
    break;

  case CT_DTsynthesise:
    opt.NewParameterNamePrefix="DT_synthesised_";
    opt.Centroid=CNT_UseExistingCentroid;
    break;

  case CT_loadFeatureDefinitions:
  case CT_loadQuestions:
  case CT_loadFeatureValues:
  case CT_saveTree:
  case CT_loadTree:
    opt.Filename="";
    break;

  default:
    error("Command type %s is not yet supported", CommandType_to_string[ctype].c_str());
  }
}



/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::find_matching_command_parameters
 *      expands the regular expressions in each command into actual
 *      parameter names
 * 
 * Preconditions:
 *      commands should have been loaded; GM_Parms must contain the
 *      model parameters
 *
 * Postconditions:
 *      each command has 'params' filled in
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
GMTK_Tie::find_matching_command_parameters()
{
  for (CommandListType::iterator i=commands.begin();i!=commands.end();i++){

    switch(i->param_type){
      
    case PT_Mixture:
      expand_param_names(&(GM_Parms->mixturesMap), i->param_expression, &(i->params));
      break;
    case PT_Component:
      expand_param_names(&(GM_Parms->componentsMap), i->param_expression, &(i->params));
      break;
    case PT_Mean:
      expand_param_names(&(GM_Parms->meansMap), i->param_expression, &(i->params));
      break;
    case PT_Collection:
	error("Collection is not yet a supported parameter type");
	break;
    case PT_NotRequired:
      break;
    default:
	error("Unsupported parameter type");

    }
  }
}


/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::validate_command_parameters
 *      does some simple sanity checking on the options for a command;
 *      currently, this means just checking that the specified
 *      parameters exist and are of the correct type
 * 
 * Preconditions:
 *      'commands' must be filled in and GM_Parms must contain the
 *      model paramters
 *
 * Postconditions:
 *      we can be sure the commands are at least somewhat correct
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      true if no problems are found; false otherwise
 *
 *-----------------------------------------------------------------------
 */
bool 
GMTK_Tie::validate_command_parameters()
{
  bool rval=false;

  for (CommandListType::iterator i=commands.begin();i!=commands.end();i++){

    // check the parameter(s)
    for (std::vector<std::string>::iterator j=i->params.begin();j!=i->params.end();j++){
      
      switch(i->param_type) {
      case PT_Mixture:
	rval = check_parameter_exists(&(GM_Parms->mixturesMap),*j);
	break;
      case PT_Component:
	rval = check_parameter_exists(&(GM_Parms->componentsMap),*j);
	break;
      case PT_Mean:
	rval = check_parameter_exists(&(GM_Parms->meansMap),*j);
	break;
      case PT_Collection:
	rval = check_parameter_exists(&(GM_Parms->nclsMap),*j);
	break;
      default:
	error("Unsupported parameter type");
	rval=false;
      }
      
      if(!rval){
	error("Parameter '%s' does not exist or is listed as wrong type in command file",j->c_str());
	return false;
      }
    }

    /*
    // check the variables, if any
    for (std::vector<std::string>::iterator j=i->variables.begin();j!=i->variables.end();j++){
      
      // look in the graph to find the matching RV
      std::string rv=j->substr(0,j->find_first_of("="));

      fprintf(stderr,"Looking for RV called %s\n",rv.c_str());

      if(!rval){
	error("Variable '%s' does not exist in the graphs",j->c_str());
	return false;
      }
    }
    */
  }

  return true;
}


/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::print_Command
 *      pretty-prints a command
 * 
 * Preconditions:
 *      commands must be filled in
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      nonw
 *
 *-----------------------------------------------------------------------
 */
void
GMTK_Tie::print_Command(int index, IM::VerbosityLevels v)
{

  Command cmd;
  cmd=commands[index];

  infoMsg(v,"Command %d: %s\n",index,CommandType_to_string[cmd.command].c_str());

  switch(cmd.command){

  case CT_Unknown:
    break;


  case CT_Tie:
    infoMsg(v," ParameterType: %s\n",ParameterType_to_string[cmd.param_type].c_str());
    infoMsg(v," NewParameterNamePrefix: %s\n",cmd.options.NewParameterNamePrefix.c_str());
    infoMsg(v," %d Matching parameters",cmd.params.size());
    //infoMsg(IM::Max,": ");
    //for(std::vector<std::string>::iterator i=cmd.params.begin();i!=cmd.params.end();i++)
    //  infoMsg(IM::Max," %s", i->c_str());
    infoMsg(v,"\n");
    break;


  case CT_Cluster:
    infoMsg(v," ParameterType: %s\n",ParameterType_to_string[cmd.param_type].c_str());
    infoMsg(v," MaxClusterSize: %1.9e\n",cmd.options.MaxClusterSize);
    infoMsg(v," MinClusterMembers: %d\n",cmd.options.MinClusterMembers);
    infoMsg(v," DissimilarityMeasure: %s\n",DissimilarityMeasureType_to_string[cmd.options.DissimilarityMeasure].c_str());
    infoMsg(v," ClusterSizeMethod: %s\n",ClusterSizeMethodType_to_string[cmd.options.ClusterSizeMethod].c_str());
    infoMsg(v," Centroid: %s\n",CentroidType_to_string[cmd.options.Centroid].c_str());

    if(cmd.options.CentroidOccupancyWeighting)
      infoMsg(v," CentroidOccupancyWeighting is on\n");

    if(cmd.options.ClusterSizeOccupancyWeighting)
      infoMsg(v," ClusterSizeOccupancyWeighting is on\n");

    infoMsg(v," NewParameterNamePrefix: %s\n",cmd.options.NewParameterNamePrefix.c_str());
    infoMsg(v," %d Matching parameters: ",cmd.params.size());
    //infoMsg(v," Matching parameters: ");
    //for(std::vector<std::string>::iterator i=cmd.params.begin();i!=cmd.params.end();i++)
    //  infoMsg(v+10," %s", i->c_str());
    infoMsg(v,"\n");
    break;


  case CT_DTcluster:
    infoMsg(v," TreeName: %s\n",cmd.options.TreeName.c_str());
    infoMsg(v," ParameterType: %s\n",ParameterType_to_string[cmd.param_type].c_str());
    infoMsg(v," MaxClusterSize: %1.9e\n",cmd.options.MaxClusterSize);
    infoMsg(v," MaxNumberParameters: %d\n",cmd.options.MaxNumberParameters);
    infoMsg(v," MinOccupancyCount: %e\n",cmd.options.MinOccupancyCount);
    infoMsg(v," ThresholdOccupancyCount: %e\n",cmd.options.ThresholdOccupancyCount);
    infoMsg(v," MinClusterMembers: %d\n",cmd.options.MinClusterMembers);
    infoMsg(v," MinImprovementPercent: %3.2e\n",cmd.options.MinImprovementPercent);
    infoMsg(v," MinImprovementAbsolute: %e\n",cmd.options.MinImprovementAbsolute);

    if(cmd.options.CentroidOccupancyWeighting)
      infoMsg(v," CentroidOccupancyWeighting is on\n");

    if(cmd.options.ClusterSizeOccupancyWeighting)
      infoMsg(v," ClusterSizeOccupancyWeighting is on\n");

    infoMsg(v," NewParameterNamePrefix: %s\n",cmd.options.NewParameterNamePrefix.c_str());
    infoMsg(v," %d Matching parameters: ",cmd.params.size());
    //infoMsg(v," Matching parameters: ");
    //for(std::vector<std::string>::iterator i=cmd.params.begin();i!=cmd.params.end();i++)
    //  infoMsg(v+10," %s", i->c_str());
    infoMsg(v,"\n");
    break;

  case CT_Untie:
    infoMsg(v," NewParameterNamePrefix: %s\n",cmd.options.NewParameterNamePrefix.c_str());
    infoMsg(v," %d Matching parameters: ",cmd.params.size());
    //infoMsg(v," Matching parameters: ");
    //for(std::vector<std::string>::iterator i=cmd.params.begin();i!=cmd.params.end();i++)
    //  infoMsg(v+10," %s", i->c_str());
    //infoMsg(v,"\n");
    break;


  case CT_DTsynthesise:
    infoMsg(v," TreeName: %s\n",cmd.options.TreeName.c_str());
    infoMsg(v," ParameterType: %s\n",ParameterType_to_string[cmd.param_type].c_str());
    infoMsg(v," NewParameterNamePrefix: %s\n",cmd.options.NewParameterNamePrefix.c_str());
    infoMsg(v," %d Matching parameters: ",cmd.params.size());
    //infoMsg(v," Matching parameters: ");
    //for(std::vector<std::string>::iterator i=cmd.params.begin();i!=cmd.params.end();i++)
    //  infoMsg(v+10," %s", i->c_str());
    infoMsg(v,"\n");
    break;



  case CT_loadFeatureDefinitions:
    infoMsg(v,"  FeatureSetName: %s\n",cmd.options.FeatureSetName.c_str());
    infoMsg(v,"  Filename: %s\n",cmd.options.Filename.c_str());
    break;


  case CT_loadQuestions:
    infoMsg(v,"  QuestionSetName: %s\n",cmd.options.QuestionSetName.c_str());
    infoMsg(v,"  FeatureSetName: %s\n",cmd.options.FeatureSetName.c_str());
    infoMsg(v,"  FeatureValuesName: %s\n",cmd.options.FeatureValuesName.c_str());
    infoMsg(v,"  Filename: %s\n",cmd.options.Filename.c_str());
    break;


  case CT_loadFeatureValues:
    infoMsg(v,"  FeatureValuesName: %s\n",cmd.options.FeatureValuesName.c_str());
    infoMsg(v,"  FeatureSetName: %s\n",cmd.options.FeatureSetName.c_str());
    infoMsg(v,"  Filename: %s\n",cmd.options.Filename.c_str());
    break;


  case CT_saveTree:
    infoMsg(v,"  TreeName: %s\n",cmd.options.TreeName.c_str());
    infoMsg(v,"  Filename: %s\n",cmd.options.Filename.c_str());
    break;


  case CT_loadTree:
    infoMsg(v,"  TreeName: %s\n",cmd.options.TreeName.c_str());
    infoMsg(v,"  Filename: %s\n",cmd.options.Filename.c_str());
    break;


  }


}




/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::print_Options
 *      pretty-prints all options for a single command
 * 
 * Preconditions:
 *      commands must be filled in 
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
void
GMTK_Tie::print_Options(Options &opt, IM::VerbosityLevels v)
{
  
  // to do: detect special values (e.g. -1) and adjust what is printed


}






/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::execute_command
 *      executes a single command
 * 
 * Preconditions:
 *      commands must be filled in
 *
 * Postconditions:
 *      the specified clustering/tying/untying has been performed
 *
 * Side Effects:
 *      global GM_Parms may have changed
 *
 * Results:
 *      true, if successful; false otherwise
 *
 *-----------------------------------------------------------------------
 */
bool 
GMTK_Tie::execute_command(unsigned command_index)
{
  infoMsg(IM::Default,"\n----Executing-----\n");
  print_Command(command_index,IM::Default);

  if ( (commands[command_index].params.size() < 1) and (commands[command_index].param_type != PT_NotRequired) ){
    warning("There are no parameters matching this command - doing nothing!\n");
    return true;
  }

  bool rval=true;

  switch(commands[command_index].command){

  case CT_Unknown:
    error("Command is unknown");
    break;

  case CT_Tie:
    
    switch(commands[command_index].param_type) {


  case PT_Mixture:
    if (tie_Mixtures(commands[command_index].params, 
		       commands[command_index].options.Centroid, 
		       commands[command_index].options.CollectionName,
		       commands[command_index].options.NewParameterNamePrefix,
		       commands[command_index].options.CentroidOccupancyWeighting,
		       false) == "") // never expand expression - testing?
	error(" Failed to tie Mixtures");
      //commit_nc_changes(commands[command_index].options.CollectionName);
      break;

    case PT_Component:
      warning("Tying Components is not yet implemented; model not changed\n");
      break;

    case PT_Mean:
      warning("Tying Means is not yet implemented; model not changed\n");
      break;

    case PT_Collection:
      warning("Tying Collections is not yet implemented; model not changed\n");
      break;

    default:
      warning("Unsupported parameter type for command 'tie'. Model unchanged\n");
    }
    break;



  case CT_Cluster:
  case CT_DTcluster:
    {
      std::vector<std::vector<std::string> > clustered_params;
      std::vector<std::string> outliers;

      // clustering function is specific to the method
      switch(commands[command_index].command) {
      case CT_Cluster:
	{
	  if(!data_driven_cluster(command_index, clustered_params))
	    error(" Failed to do data_driven_cluster");
	}
	break;
	
      case CT_DTcluster:
	{
	  if(!decision_tree_cluster(command_index, clustered_params,outliers))
	    error(" Failed to do tree_based_cluster");
	}
	break;

      default:
	error("This line is here only to prevent a compiler warning about unhandled enumeration values in the switch statement\n");
      }

      // but once items are in clusters, the tying functions are the
      // same, regardless of whether we previously ran data-driven or
      // tree-based clustering
      switch(commands[command_index].param_type) {
	
      case PT_Mixture:
	{
	  for(unsigned j=0;j<clustered_params.size();j++){ // for each cluster

	    if (clustered_params[j].size() > 0 ){
	      infoMsg(IM::Mega,"Tying %d Mixtures within cluster %d\n",clustered_params[j].size(), j);
	      if (tie_Mixtures(clustered_params[j], 
			       commands[command_index].options.Centroid,
			       commands[command_index].options.CollectionName,
			       commands[command_index].options.NewParameterNamePrefix,
			       commands[command_index].options.CentroidOccupancyWeighting,false) == "") 
		error(" Failed to tie Mixtures after clustering");
	      //commit_nc_changes(commands[command_index].options.CollectionName);
	    } else
	      warning("There was an empty cluster during post-clustering tying\n");

	  }
	}
	break;
	
      case PT_Component:
	{
	  for(unsigned j=0;j<clustered_params.size();j++){ // for each cluster
	    infoMsg(IM::Mega,"Tying %d Components within cluster %d\n",clustered_params[j].size(), j);
	    //if (!(tie_Components(clustered_params[j], commands[command_index].options.Centroid,false))) 
	    //error(" Failed to tie Components after clustering");
	    error("tie_Components function is not yet implemented - cannot call it");
	  }
	}
	break;
	
      case PT_Mean:
	{
	  for(unsigned j=0;j<clustered_params.size();j++){ // for each cluster
	    infoMsg(IM::Mega,"Tying %d Means within cluster %d\n",clustered_params[j].size(), j);
	    if (tie_Means(clustered_params[j], commands[command_index].options.Centroid,false) == "") 
	      error(" Failed to tie Means after clustering");
	  }
	}
	break;


      default:
	error("Unsupported parameter type for tying after clustering");
      }

      if(commands[command_index].command==CT_DTcluster){

	// now that tying has happened, we know the name of the centroid
	// of each cluster, so we can fill in this information in the
	// decision tree; this makes it self-contained so it can be
	// saved later, or re-used (i.e. we don't need to also keep
	// clusters_params around)
	DecisionTreeType *decision_tree = &(decision_trees[commands[command_index].options.TreeName]);


	bool tree_is_valid=false;
	for(unsigned i=1;i<decision_tree->tree.size();i++)
	  if( decision_tree->tree[i]==-1 ){ // terminal node
	    decision_tree->nodes_to_param_names[i]=clustered_params[decision_tree->nodes_to_clusters[i]][0];
	    tree_is_valid=true;
	  }

	// an empty decision tree (only a root node, with no
	// associated parameters) may arise if we previously tried to
	// DT cluster a set of params, all of which were removed as
	// low-occupancy outliers prior to tree building
	if (!tree_is_valid){
	  warning("An invalid decision tree was found during command %d - probably because DT clustering failed earlier. Instead of synthesising, the parameters will simply all be tied together\n",command_index);
	  
	  // hardwire to Mixtures for now
	  if (commands[command_index].param_type != PT_Mixture)
	    error("Parameter type was not Mixtures - bailing out!");

	  if (tie_Mixtures(outliers, 
			   CNT_UseExistingCentroid,
			   commands[command_index].options.CollectionName,
			   commands[command_index].options.NewParameterNamePrefix,
			   false,false) == "") 
	    error(" Failed to tie Mixtures (as a fallback from synthesising them with a decision tree");



	} else {
	
	  // so now that the tree is complete, we can synthesise all the
	  // outliers that were set aside before
	  infoMsg(IM::Mod,"Synthesising %d parameters using tree '%s'\n",outliers.size(),decision_tree->name.c_str());
	  
	  decision_tree_synthesise(outliers.begin(),outliers.end(),decision_tree,commands[command_index].options.CollectionName);
	}


      }    
    }
    break;
    
    



  case CT_Untie:
    
    switch(commands[command_index].param_type) {
    case PT_Mixture:
      // two cases to deal with:
      //
      // 1) there are already two Mixture objects and we need to give
      //    them independent copies of their Components and weights
      //
      // 2) there is only one Mixture object, so we need to duplicate
      //    it and also make uniques names in the appropriate
      //    NameCollection
      //
      // the user tells us which case this is, by either providing a
      // name collection, or not
      if (!(untie_Mixtures(commands[command_index].params,
			   commands[command_index].options.CollectionName,
			   commands[command_index].options.NewParameterNamePrefix,
			   false))) // never expand expressions: for testing only??
			error(" Failed to untie Mixtures");
			break;
      
    case PT_Component:
      warning("Untying Components is not yet implemented; model not changed\n");
      break;
    case PT_Mean:
      warning("Untying Means is not yet implemented; model not changed\n");
      break;
    case PT_Collection:
      warning("Untying Collections is not yet implemented; model not changed\n");
      break;
    default:
      warning("Unsupported parameter type for command 'untie'. Model unchanged\n");
    }
    break;


  case CT_DTsynthesise:
    {

      if (commands[command_index].options.TreeName == "")
	error("No tree name was given in DTsynthesise command - this is a required parameter");

      if (decision_trees.find(commands[command_index].options.TreeName) == decision_trees.end())
	error("Requested a tree called '%s' but no such tree has been loaded or created",commands[command_index].options.TreeName.c_str());
      else{
	// this locates the tree leaf for each param to be
	// synthesised, then ties the parameter to the parameter found
	// at that leaf
	decision_tree_synthesise(command_index);
      }
    }
    break;

  case CT_loadFeatureDefinitions:
    {
      infoMsg(IM::Default,"Loading feature definitions from '%s'\n",commands[command_index].options.Filename.c_str());
      const char *f=commands[command_index].options.Filename.c_str();
      iDataStreamFile cf(f,false,true,_cppCommandOptions);
      read_feature_definition_set(cf);
    }
    break;

  case CT_loadQuestions:
    {
      infoMsg(IM::Default,"Loading questions from '%s'\n",commands[command_index].options.Filename.c_str());
      const char *f=commands[command_index].options.Filename.c_str();
      iDataStreamFile cf(f,false,true,_cppCommandOptions);
      read_question_set(cf);
    }
    break;

  case CT_loadFeatureValues:
    {
      infoMsg(IM::Default,"Loading feature values from '%s'\n",commands[command_index].options.Filename.c_str());
      const char *f=commands[command_index].options.Filename.c_str();
      iDataStreamFile cf(f,false,true,_cppCommandOptions);
      read_feature_value_set(cf);
    }
    break;

  case CT_saveTree:
    {
      save_decision_tree(command_index);
    }
    break;

  case CT_loadTree:
    {
      load_decision_tree(command_index);
    }
    break;

  }
  infoMsg(IM::Default,"----Finished-----\n");

  return rval;
}





/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::tie_Means
 *      ties the Means of the Mixtures specified
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      the MeanVectors are tied and GM_Parms is updated accordingly
 *
 * Side Effects:
 *      global GM_Parms may have changed
 *
 * Results:
 *      true if successful; false otherwise
 *
 *-----------------------------------------------------------------------
 */
std::string 
GMTK_Tie::tie_Means(std::vector<std::string> param_expressions, CentroidType method, bool expand_expressions)
{
  std::vector<std::string>::iterator i;
  std::vector<string> params;

  // to do: share this bit with other functions
  if(expand_expressions)
    // expand out the list of regex, param_expressions, into actual
    // parameter names. 
    for(i=param_expressions.begin();i!=param_expressions.end();i++)
      expand_param_names(&GM_Parms->mixturesMap,*i,&params);
  else
    params.swap(param_expressions);

  if ( params.size() == 0){
    warning("Tying with no parameters makes no sense - doing nothing for this command!\n");
    return "";
  }

  if ( params.size() == 1){
    warning("Tying a single parameter has no effect; model is unchanged.\n");
    return params[0];
  }


  switch(method){

  case CNT_Unknown:
    error("Unknown Centroid type in tie_Means");
    break;


  case CNT_CreateCentroid_averageMeanVector:
    {
      // first create the average, overwriting the first mean
      i=params.begin();
      MeanVector* average_mean = find_MeanVector(*i);

      i++; // now pointing at second item
      for(;i!=params.end();i++)
	average_mean->means += find_MeanVector(*i)->means;

      average_mean->means /= (float)(params.size());

      // then tie all others to this one, by continuing below...
    }
    // do not break here


  case CNT_Arbitrary:
  case CNT_UseExistingCentroid:
  case CNT_EmulateHTKOther:
    {
      // after clustering (if any has been done) the first item in the
      // parameter list will be the centroid (applies to
      // CNT_UseExistingCentroid method only)
            
      MeanVector* tied_mean = NULL;

      std::list<Component*> components_that_use_these_means;
      std::list<std::string> component_names;
      unsigned n_users, total_users=0;

      i++; // now pointing to second param
      for(i=params.begin();i!=params.end();i++){

	  MeanVector* this_mean = find_MeanVector(*i);
	  
	  // now, find out who is using this MeanVector, so that we can
	  // tell them to point to the tied one instead
	  n_users = find_Components_using_MeanVector(this_mean,(*i),&components_that_use_these_means,&component_names);
	  total_users += n_users;

	  if (n_users == 0)
	    warning("MeanVector '%s' will be tied, although could not find any Components that use it\n",i->c_str());

	  if (i==params.begin()){
	    // keep the first mean and tie the rest to it
	    tied_mean = this_mean;

	    // to do: rename the mean (and update the map) using
	    // user-supplied prefix
	    
	  } else {
	  
	    infoMsg(IM::Low,"Tying %s (used %d times) to %s\n",i->c_str(),n_users,params.begin()->c_str());
	  
	  }

      }

      // adjust all the found Components to use the tied mean
      for(std::list<Component*>::iterator ci=components_that_use_these_means.begin();
	  ci!=components_that_use_these_means.end(); ci++){

	if( !is_DiagGaussian(*ci) )
	  error("Found a non-DiagGaussian Component, so cannot adjust its MeanVector");

	((DiagGaussian*)(*ci))->mean = tied_mean;

      }

      
    }
    break;



  default:
    error("'%s' is not a valid Centroid selection method for tying Means",CentroidType_to_string[method].c_str());
    break;

  }

  // Clean up
  // --------

  //  delete any now-unused  Components from GM_Params
  // THIS IS TOO SLOW TO CALL HERE GM_Parms->markUsedMixtureComponents(); // does this actually delete anything?
  // and any now-unused  Dense1DPMF from GM_Params
  // TO DO - no function currently exists for this in class GMParms
  // TO DO - other param types too?
    

  return params[0];
}




/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::tie_Mixtures
 *      tie the underlying Components and weights of the specified
 *      Mixtures, OR tie within the named collection
 * 
 * Preconditions:
 *      all the Mixtures must have a single type of Component
 *
 * Postconditions:
 *      the Components and weights are tied and GM_Parms is updated
 *      accordingly
 *
 * Side Effects:
 *      1) global GM_Parms may have changed
 *      2) unused Components are deleted
 *
 * Results:
 *      true is successful; false otherwise
 *
 *-----------------------------------------------------------------------
 */
std::string 
GMTK_Tie::tie_Mixtures(std::vector<std::string> param_expressions, 
		       CentroidType method,
		       std::string collection_name,
		       std::string name_prefix,
		       bool occupancy_weighted, bool expand_expressions)
{
  infoMsg(IM::Huge,"tie_Mixtures\n");

  bool tie_at_collection_level=false;
  NameCollection* nc=NULL;

  if (collection_name !=""){
    // check that this collection exists, and if it does we will tie
    // Mixtures at the collection level (i.e. by manipulating entries
    // in that named collection table, rather than messing about with
    // components etc)
    GMParms::ObjectMapType::iterator i=GM_Parms->nclsMap.find(collection_name);
    if (i == GM_Parms->nclsMap.end())
      error("Cannot find named collection called %s",collection_name.c_str());
    else
      nc=GM_Parms->ncls[i->second];

    tie_at_collection_level=true;
  }

  infoMsg(IM::Huge,"tie_Mixtures - collection name is %s\n",collection_name.c_str());

  std::vector<std::string>::iterator i;
  std::vector<string> params;

  if(expand_expressions){
    // expand out the list of regex, param_expressions, into actual
    // parameter names. 
    infoMsg(IM::Huge,"tie_Mixtures - expanding regular expressions\n");
    for(i=param_expressions.begin();i!=param_expressions.end();i++)
      expand_param_names(&GM_Parms->mixturesMap,*i,&params);
  } else
    params=param_expressions;

  if ( params.size() == 0){
    warning("Tying with no parameters makes no sense - doing nothing for this command!\n");
    return "";
  }

  if ( params.size() == 1){
    if(expand_expressions)
      warning("Tying a single parameter has no effect; model is unchanged.\n");
    return params[0];
  }

  // first, locate the Mixtures concerned in the list of Mixtures and
  // find all their mixture weights and sets of Gaussian components

  /*
  // a list of the weight names
  std::vector<std::string> weights;
  weights.resize(params.size());

  // a list of the Gaussian component
  std::vector< std::vector<std::string> > components;
  components.resize(params.size());

  int p_counter=0;
  for(i=params.begin();i!=params.end();i++,p_counter++){

    Mixture* this_mixture = find_Mixture(*i); // GM_Parms.mixtures[GM_Parms.mixturesMap[*i]];
    components[p_counter].resize(this_mixture->components.size());

    int c_counter=0;
    for(vector< Component* >::iterator j=this_mixture->components.begin();j!=this_mixture->components.end();j++,c_counter++){
      std::string n=(*j)->name();
      components[p_counter][c_counter]=n;
    }

    weights[p_counter] = this_mixture->dense1DPMF->name();

    if ( weights[p_counter] == ""){
      error("Failed to locate the Dense1DPMF for Mixture %s",i->c_str());
      return "";
    }

  }

  */

  // create a replacement dense PMF weight and a set of replacement
  // Gaussian components 

  // first version: just use the first Mixture with
  // unchanged parameters

  // second version: merging the mixture components. User chooses one
  // of these methods:

  // * just take the current values of params for the first of the
  // * list of Mixtures

  // * or just take the current values of params for the centroid
  // * Mixture in each cluster

  // * for a set of single-component mixtures, a simple averaging of
  // * means, variances, mix weights

  // * for a set of N mixtures with the same number of components (M),
  // it's a permutation (to get the M components in the "same" order)
  // and then averaging of the set of N Gaussians and corresponding
  // weights

  // * for Mixtures with differing numbers of components ... to be
  // decided!


  // to do: try other methods of merging/averaging of the parameters


  infoMsg(IM::Huge,"tie_Mixtures: tying %d parameters to a Centroid selected by method '%s'\n",params.size(),CentroidType_to_string[method].c_str());

  switch(method){

  case CNT_Unknown:
    error("Unknown Centroid type in tie_Mixtures");
    break;

  case CNT_CreateCentroid_averageMeanVector:
    error("'%s' is not a valid Centroid selection method for tying Mixtures",CentroidType_to_string[method].c_str());
    break;

  case CNT_Arbitrary:
  case CNT_UseExistingCentroid:
    {
      // after clustering (if any) the first item in the parameter list
      // will be the centroid (CNT_UseExistingCentroid method only)
      
      // keep the first mixture and tie the rest to it
      i=params.begin();
      std::string tied_mixture_name=*i;
      Mixture* tied_mixture = find_Mixture(*i); // GM_Parms.mixtures[GM_Parms.mixturesMap[*i]];
      tied_mixture_name=*i;

      //rename by adding a user-supplied prefix
      //infoMsg(IM::Huge,"tie_Mixtures UseExistingCentroid - renamed %s",tied_mixture_name.c_str());
      //tied_mixture->setName(tied_mixture_name);
      //infoMsg(IM::Huge,"to %s\n",tied_mixture_name.c_str());

      if(tie_at_collection_level){
	infoMsg(IM::Huge,"tie_Mixtures UseExistingCentroid - tying at collection level with centroid %s\n",tied_mixture_name.c_str());
	
	std::vector<string> remaining_params = params;
	remaining_params.erase(remaining_params.begin());
	nc->queue_search_and_replace(remaining_params,tied_mixture_name);



      } else {
	// tie below the Mixture level by tying components & mixture weights

	i++; // now pointing to second param
	for(;i!=params.end();i++){
	  Mixture* this_mixture = find_Mixture(*i); // GM_Parms.mixtures[GM_Parms.mixturesMap[*i]];
	  
	  // delete existing Components from the Mixture object
	  this_mixture->components.clear();
	  
	  // insert tied components (works for multi-component mixtures too)
	  this_mixture->components.assign(tied_mixture->components.begin(),tied_mixture->components.end());
	  this_mixture->dense1DPMF = tied_mixture->dense1DPMF ;
	}
      }
    }
    break;



  case CNT_CreateCentroid_averageSingleComponentMixtures:
    {
      // only for single-component Mixtures of diagonal Gaussians (at
      // the moment, at least)

      // to avoid overflow, we need to calculate the total occupancy
      // up front, so we can be dividing by it inside the loop


      double this_occupancy=0.0,total_occupancy=0.0;
      i=params.begin();
      Mixture *this_mixture;
      DiagGaussian *this_gaussian;

      for(;i!=params.end();i++){
	
	this_mixture = find_Mixture(*i); // GM_Parms.mixtures[GM_Parms.mixturesMap[*i]];
	
	if (this_mixture->components.size() != 1)
	  error ("CreateCentroid_averageSingleComponentMixtures only available for Mixtures with exactly one Component");
	
	if( !is_DiagGaussian(this_mixture->components[0]))
	  error("CreateCentroid_averageSingleComponentMixtures only available for DiagGaussian Components, not %s",this_mixture->components[0]->typeName().c_str());
	if(occupancy_weighted){
	  this_gaussian = (DiagGaussian*)this_mixture->components[0];
	  total_occupancy +=this_gaussian->get_accumulatedProbability().unlog();
	}
      }
      
      
      
      // keep first Mixture's mean/var/weights, sum all other components
      // into them, then divide by total occupancy, or num parameters
      
      // CANNOT use same obj though!!!

      // this will currently do the wrong thing with occupancies of
      // shared parameters: they will be counted too many times; to
      // fix this we could zero each occupancy value as it is added
      // it, so that future attempts to add it again will simply add
      // zero - to do!


      if(occupancy_weighted){
	if(total_occupancy < MIN_OCCUPANCY)
	  warning("There is insufficient total occupancy (%e) in this cluster to do occupancy weighting when computing the centroid - falling back to unweighted version\n",total_occupancy);
	occupancy_weighted = false;
      } 

      //cerr << "total occ=" << total_occupancy << endl;

      i=params.begin();


      
      this_mixture = find_Mixture(*i);
      this_gaussian = (DiagGaussian*)this_mixture->components[0];

      sArray<float> tied_mean,tied_covar;
      unsigned k,dim=this_gaussian->mean->means.len();
      double weight=1/((double)(params.size())); // default weight when no occupancy weighting is done
      tied_mean.resizeAndZero(dim);
      tied_covar.resizeAndZero(dim);

      //for (k=0; k<dim; k++)
      //cerr << "at start: " <<tied_mean[k] << endl;


      for(i=params.begin();i!=params.end();i++){
	
	this_mixture = find_Mixture(*i);
	this_gaussian = (DiagGaussian*)this_mixture->components[0];

	if(occupancy_weighted){
	  this_occupancy=this_gaussian->get_accumulatedProbability().unlog();
	  weight=(this_occupancy / total_occupancy);
	  //cerr << "weight=" << weight << endl;
	}

	for (k=0; k<dim; k++){
	  tied_mean[k] += this_gaussian->mean->means[k] * weight;
	  tied_covar[k] += this_gaussian->covar->covariances[k] * weight;
	}
	
      }

      // now doing the tying
      //
      // first, the first parameter will the one to be kept, so set
      // its mean and covar to the computed average mean and covar


      i=params.begin();
      Mixture* tied_mixture = find_Mixture(*i);


      //rename by adding a user-supplied prefix
      std::string tied_mixture_name=*i;
      //tied_mixture->setName(tied_mixture_name);

      cerr << "Keeping Mixture " << tied_mixture_name << endl;

      for (k=0; k<dim; k++){
	this_gaussian = (DiagGaussian*)tied_mixture->components[0];
	this_gaussian->mean->means[k] = tied_mean[k];
	this_gaussian->covar->covariances[k] = tied_covar[k];
      }

      if(tie_at_collection_level){
	
	std::vector<string> remaining_params = params;
	remaining_params.erase(remaining_params.begin());
	nc->queue_search_and_replace(remaining_params,tied_mixture_name);


      } else {

	i++;
	for(;i!=params.end();i++){
	  
	  // tie below the Mixture level by tying components & mixture weights
	  this_mixture = find_Mixture(*i);
	  this_mixture->components.clear();
	  this_mixture->components.assign(tied_mixture->components.begin(),tied_mixture->components.end());
	  this_mixture->dense1DPMF = tied_mixture->dense1DPMF;
	}
      }


      //cerr << "total_occupancy=" << total_occupancy << endl;

      //float M=(float)params.size() * total_occupancy;

      
    }
    break;



  case CNT_CreateCentroid_permuteThenAverage:
    // number of components per Mixture must be the same; permute to
    // get components in "same" order (based on distances between the
    // means?), then average the matched-up parameters
    error("CreateCentroid_permuteThenAverage not yet implemented");
    break;

  case CNT_CreateCentroid_smartMerge:
    error("smartMerge not yet implemented");
    break;

  case CNT_EmulateHTKMixturesOfGaussians:
    error("CNT_EmulateHTKMixturesOfGaussians not yet implemented for tying Mixtures");
    break;

  default:
    error("'%s' is not a valid Centroid selection method for tying Mixtures",CentroidType_to_string[method].c_str());
    break;


  }

  // Clean up
  // --------

  //  delete any now-unused  Components from GM_Params
  // THIS IS TOO SLOW TO CALL HERE  GM_Parms->markUsedMixtureComponents(); // does this actually delete anything?
  // and any now-unused  Dense1DPMF from GM_Params
  // TO DO - no function currently exists for this in class GMParms
    

  return params[0];
}





/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::untie_Mixtures
 *      unties the specified Mixtures
 * 
 * Preconditions:
 *      none; no check is made that the Mixtures being untied actually
 *      share anything in the first place
 *
 * Postconditions:
 *      all Mixtures end up with their own unique copies of the
 *      Components and weights they used to share
 *
 * Side Effects:
 *      1) global GM_Parms may have changed
 *      2) unused Components are deleted
 *
 * Results:
 *      true if successful; false otherwise
 *
 *-----------------------------------------------------------------------
 */
bool 
GMTK_Tie::untie_Mixtures(std::vector<std::string> param_expressions, 
			 std::string collection_name,
			 std::string name_prefix,
			 bool expand_expressions)
{
  infoMsg(IM::Mod,"untie_Mixtures with %d parameter expressions\n",param_expressions.size());

  // this should work for any type of component - to do: remove any
  // previous checks for Gaussian/DiagGaussian

  // to do - share first part of this function with tie_Mixtures

  std::vector<std::string>::iterator i;
  std::vector<string> params;
  bool untie_at_collection_level=false;
  NameCollection* nc=NULL;

  if(expand_expressions)
    // expand out the list of regex, param_expressions, into actual
    // parameter names. 
    for(i=param_expressions.begin();i!=param_expressions.end();i++)
      expand_param_names(&GM_Parms->mixturesMap,*i,&params);
  else
    params=param_expressions;

  infoMsg(IM::Mod,"untie_Mixtures with %d actual parameters\n",params.size());

  if ( params.size() == 0){
    warning("Untying with no parameters makes no sense - doing nothing for this command!\n");
    return true;
  }

  if (collection_name !=""){
    // check that this collection exists, and if it does we will tie
    // Mixtures at the collection level (i.e. by manipulating entries
    // in that named collection table, rather than messing about with
    // components etc)
    GMParms::ObjectMapType::iterator i=GM_Parms->nclsMap.find(collection_name);
    if (i == GM_Parms->nclsMap.end())
      error("Cannot find named collection called %s",collection_name.c_str());
    else
      nc=GM_Parms->ncls[i->second];
    untie_at_collection_level=true;
  }

  // locate the Mixtures concerned in the list of Mixtures and find
  // all their mixture weights and sets of components



  // we don't care whether this leaves orphaned (unused) parameters -
  // they will get cleaned up somewhere else later; they will remain
  // in GM_Parms for now, so we can find them (no memory leak)
  
  int ii=0;
  int last_ii=-1001;
  int s=params.size();
  //int perc=-1,new_perc;
  
  if(untie_at_collection_level){
    // we actually need to make another Mixture object for each time
    // this name appears in the name collection
    
    vector<string>::iterator vsi;
    
    for(i=params.begin();i!=params.end();i++,ii++){
      
      if(ii > last_ii + 1000){
	last_ii = ii;
	infoMsg(IM::Mod,"untie_Mixtures: %d of %d\r",ii,s);
      }

      // could probably replace this with a fast version that looks in
      // nc->mxTable, but since we are messing about with the contents
      // of nc->sorted_table, that might be tricky to make reliable;
      // let's do it the safe way, by looking in GM_Parms


      // to do: use the queuing mechanism to avoid resorting the table
      // so many times

      Mixture* this_mixture = find_Mixture(*i);
      
      // find the start and end of the subrange 
      nc->sort();
      // the mxTable will be invalid now, so forget it, to be on the
      // safe side
      nc->mxTable.clear();

      pair<vector<string>::iterator,vector<string>::iterator> result = equal_range(nc->sorted_table.begin(),nc->sorted_table.end(),*i);
      
      // st now points at the first item in nc->sorted_table that
      // is equal to *i
      
      for(vsi=result.first;vsi!=result.second;vsi++){
	
	cerr << "Found " << *vsi << endl;
	
	// clone it (this also enters the new mixture into the
	// GM_Parms mixture table)
	this_mixture=this_mixture->identicalIndependentClone();
	
	// update the entry in nc->sorted_table with the new name (no
	// need to worry about fixing up any entries in nc->table or
	// nc->mxTable because these will be rebuilt later

	*vsi=this_mixture->name();

      }

      // finally, unsort
      nc->unsort();
      
    }

  } else {
    
    for(i=params.begin();i!=params.end();i++,ii++){
      
      if(ii > last_ii + 1000){
	last_ii = ii;
	infoMsg(IM::Mod,"untie_Mixtures: %d of %d\r",ii,s);
      }
      
      
      Mixture* this_mixture = find_Mixture(*i);
      
      // make an independent copy of Components and weights for each Mixture
      
      for(vector< Component* >::iterator j=this_mixture->components.begin();j!=this_mixture->components.end();j++){
	
	// clone this component and use the clone to replace the
	// (possibly shared) original; the identicalIndependentClone
	// function takes care of generating a new name and putting the
	// new object into the GM_Parms table
	
	// the old Component (which we are about to lose the pointer to)
	// needs to decrease its usage by 1
	(*j)->adjustNumTimesShared(-1);
	
	*j = (*j)->identicalIndependentClone();
	
	// the new Component we just made now gets a usage count of 1
	(*j)->adjustNumTimesShared(1);
	
	// now adjust the usage counts of any subobjects that have them
	// to do ....?
	
   
      }

      // Dense1DPMFs do not have numTimesShared
      this_mixture->dense1DPMF = this_mixture->dense1DPMF->identicalIndependentClone();

    }

  }


  infoMsg(IM::Mod,"\n",ii,params.size());

  return true;
}



/*-
 *-----------------------------------------------------------------------
 * print_cluster_stats
 *      reports some statistics about the current set of clusters
 * 
 * Preconditions:
 *      there should be some clusters
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
void
print_cluster_stats(std::list<Cluster> *clusters, IM::VerbosityLevels v=IM::Tiny)
{
  std::list<Cluster>::iterator ci=clusters->begin();
  double low_occ=ci->occupancy,high_occ=ci->occupancy;
  unsigned low_count=ci->items.size(), high_count=ci->items.size();

  for(;ci!=clusters->end();ci++){
    if ( ci->occupancy < low_occ)
      low_occ = ci->occupancy;
    if ( ci->occupancy > high_occ)
      high_occ = ci->occupancy;
    if ( ci->items.size() < low_count)
      low_count = ci->items.size();
    if ( ci->items.size() > high_count)
      high_count = ci->items.size();
  }

  infoMsg(v,"Cluster stats: %d clusters with occupancies ranging from %1.9e to %1.9e and # members from %d to %d\n",
	  clusters->size(),low_occ,high_occ,low_count,high_count);

}


/*-
 *-----------------------------------------------------------------------
 * merge_clusters
 *      merges a pair of clusters
 * 
 * Preconditions:
 *      1) the merged_size should already have been computed
 *      2) the clusters should be mergeable (i.e. have the same type
 *      of members)
 *
 * Postconditions:
 *      1) the clusters are merged
 *      2) clusters has one fewer member
 *      4) ci2 is no longer a valid iterator
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
merge_clusters(std::list<Cluster> *clusters, std::list<Cluster>::iterator ci1, 
			 std::list<Cluster>::iterator ci2, float merged_size)
{
  // set new cluster's size
  ci1->size=merged_size;
  ci1->occupancy += ci2->occupancy; 
  
  // merge items from cluster ci1 into cluster save
  for(std::list<Clusterable*>::iterator ii=ci2->items.begin();ii!=ci2->items.end();ii++)
    ci1->items.push_back(*ii);
  
  // delete ci
  clusters->erase(ci2);
  
}






/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::data_driven_cluster
 *      agglomerative clustering followed by outlier removal
 * 
 * Preconditions:
 *      GM_Parms must be valid
 *
 * Postconditions:
 *      clustered_params contains the clusters (e.g. for passing to a
 *      tying function)
 *
 * Side Effects:
 *      none (no actual tying takes place here)
 *
 * Results:
 *      true if successful; false otherwise
 *
 *-----------------------------------------------------------------------
 */
bool 
GMTK_Tie::data_driven_cluster(unsigned command_index, std::vector<std::vector<std::string> > &clustered_params)
{
  if (commands[command_index].command != CT_Cluster)
    error("data_driven_cluster was called, but the command is not cluster!");

  // simple bottom-up clustering using only the values of the parameters themselves
  
  // TO DO: insist on homogeneous parameter type within the set of
  // parameters being clustered


  std::list<Cluster> clusters;
  std::list<Cluster>::iterator ci,ci2,save1=NULL_ITERATOR,save2=NULL_ITERATOR;
  float this_size, saved_size=LZERO;
  int p_counter=0;
  // make the initial set of clusters - one per parameter
  for(std::vector<std::string>::iterator i=commands[command_index].params.begin();i!=commands[command_index].params.end();i++,p_counter++){
    
    Clusterable *c=NULL;
    
    switch(commands[command_index].param_type){
      
    case GMTK_Tie::PT_Unknown:
      error("Cannot call data_driven_cluster with parameter type PT_Unknown: this is a bug");
      break;

    case PT_Mixture:
      {
	Mixture *this_mixture = find_Mixture(*i);
	c = new ClusterableMixture(*i, this_mixture, commands[command_index].options.DissimilarityMeasure);
      }
      break;

    case PT_Component:
      error("data_driven_cluster PT_Component not yet implemented");
      break;

    case PT_Mean:
      {
	MeanVector *this_mean = find_MeanVector(*i);
	c = new ClusterableMean(*i, this_mean, commands[command_index].options.DissimilarityMeasure);
      }
      break;

    case PT_Collection:
      error("data_driven_cluster PT_Collection not yet implemented");
      break;

    case PT_NotRequired:
      error("Cannot call data_driven_cluster with parameter type PT_NotRequired: this is a bug!");
      break;
    }
    
    // each object starts off in its own cluster
    // so is automatically the centroid of that cluster
    Cluster new_cluster;
    new_cluster.size = 0.0;
    new_cluster.occupancy = c->occupancy();
    new_cluster.items.push_back(c);
    clusters.push_back(new_cluster);
  }
  unsigned n=clusters.size();


  if (n < 2){
    warning("Clustering with fewer than 2 parameters makes no sense - doing nothing for this command!\n");
    return true;
  }


  // -----------------------------------------------------------------------------------
  // Agglomerative clustering
  // -----------------------------------------------------------------------------------


  // somewhere to cache the resulting size from merging a pair of
  // clusters will be indexed by cluster "ident"s, so first set those
  // idents to be useful

  sArray< sArray<float> > merged_sizes;
  merged_sizes.resize(n);
  unsigned this_ident=0;
  for(ci=clusters.begin();ci!=clusters.end();ci++){
    ci->ident=this_ident;
    merged_sizes[ci->ident].resize(n);
    this_ident++;
  }

  build_merged_sizes_table(&merged_sizes,&clusters,
			   commands[command_index].options.ClusterSizeMethod,
			   commands[command_index].options.Centroid);


  bool finished=false;
  while (!finished){
    
    // find the two closest clusters and merge them, i.e. those two
    // clusters that, if merged, would have the smallest size
    
    // iterate over every pair of clusters - this is inefficient: we
    // could cache the computations for most pairs, although it would
    // be tedious to implement since clusters keep getting merged

    if (clusters.size() < 2){
      finished=true;
      continue;
    }

    saved_size=LZERO;
    for(ci=clusters.begin();ci!=clusters.end();ci++)
      for(ci2=ci;ci2!=clusters.end();ci2++){
	if(ci==ci2) continue;
	
	this_size=merged_sizes[ci->ident][ci2->ident];
	
	if ( (this_size < saved_size) or (saved_size <= LZERO) ){
	  save1=ci; save2=ci2;
	  saved_size=this_size;
	}
      }
    
    if (saved_size <= LZERO)
      error("Something went wrong during agglomerative clustering: couldn't find pair of clusters to merge");
      
    infoMsg(IM::Mod,"\nMerging clusters with sizes of %1.9e,%1.9e and new size %1.9e\n",
	    save1->size,save2->size,saved_size);
    
    if (saved_size >commands[command_index].options.MaxClusterSize){
      finished=true;
      continue;
    }

    
    // in the merging, save1 will consume save2; save2 will no longer exist
    merge_clusters(&clusters,save1,save2,saved_size);

    // recompute centroid
    if(!set_cluster_centroid(&(*save1),commands[command_index].options.Centroid))
      error("Failed to set_cluster_centroid of newly merged cluster");


    update_merged_sizes_table(&merged_sizes,&clusters,
			      commands[command_index].options.ClusterSizeMethod,
			      commands[command_index].options.Centroid, save1);


    infoMsg(IM::Default,"Size of biggest cluster out of %d clusters is currently %1.9e and occupancy is %1.9e       ",
	    clusters.size(),saved_size,save1->occupancy);
  }
  infoMsg(IM::Default,"\nAfter initial clustering:            ");
  print_cluster_stats(&clusters,IM::Default);

  // end of pure agglomerative phase
  





  // -----------------------------------------------------------------------------------
  // Outlier removal part 1: remove low occupancy clusters, only if occupancy stats are avaiable
  // -----------------------------------------------------------------------------------

  if (commands[command_index].options.MinOccupancyCount > 0.0){
    infoMsg(IM::Tiny,"Starting low occupancy cluster removal\n");

    bool no_more_low_occupancy=false;
    while (!no_more_low_occupancy){
    
      infoMsg(IM::Tiny,"Low occupancy cluster removal: # clusters is %d     ",clusters.size());

      if (clusters.size() < 2){
	infoMsg(IM::Tiny,"\n");
	warning("All items ended up in one cluster; MinOccupancyCount may not have been satisfied\n");
	no_more_low_occupancy=true;
	continue;
      }

      // find the lowest occupancy cluster
      save1=NULL_ITERATOR;
    
      ci=clusters.begin();
      double lowest_occupancy=ci->occupancy;
      save1=ci;
      for(;ci!=clusters.end();ci++){
	if (ci->occupancy < lowest_occupancy){
	  lowest_occupancy = ci->occupancy;
	  save1=ci;
	}
      }
    
      // if lowest occuapncy cluster exceeds the minimum, we are done
      if(lowest_occupancy > commands[command_index].options.MinOccupancyCount)
	save1=NULL_ITERATOR;
    
      if(save1==NULL_ITERATOR){
	no_more_low_occupancy=true;
	infoMsg(IM::Tiny,"\nMinOccupancyCount threshold reached:  ");
	print_cluster_stats(&clusters,IM::Tiny);
	continue;
      } 

      // now merge save1 into whichever other cluster results in the smallest cluster
      saved_size=LZERO;
      for(ci=clusters.begin();ci!=clusters.end();ci++){
	if(ci==save1) continue;

	this_size=merged_sizes[save1->ident][ci->ident];
      
	if ( (this_size < saved_size) or (saved_size <= LZERO) ){
	  save2=ci;
	  saved_size=this_size;
	}
      }

      if (saved_size <= LZERO)
	error("Something went wrong during part 1 of outlier removal: couldn't find a cluster to merge outlier into");

      // keep save2 - it is likely to be larger than save1
      merge_clusters(&clusters,save2,save1,saved_size);

      // recompute centroid - could avoid this until the end, but only
      // if the current ClusterSizeMethod doesn't use the centroid for
      // its calculations
      if(!set_cluster_centroid(&(*save2),commands[command_index].options.Centroid))
	error("Failed to set_cluster_centroid of newly merged cluster");

      update_merged_sizes_table(&merged_sizes,&clusters,
				commands[command_index].options.ClusterSizeMethod,
				commands[command_index].options.Centroid, save2);


    }
    infoMsg(IM::Tiny,"\n");
    infoMsg(IM::Default,"After low occupancy removal:            ");
    print_cluster_stats(&clusters,IM::Default);
  }

  // -----------------------------------------------------------------------------------
  // Outlier removal part 2: remove small clusters
  // -----------------------------------------------------------------------------------

  if ((unsigned)commands[command_index].options.MinClusterMembers > 1){

    infoMsg(IM::Default,"Starting small cluster removal\n");

    bool no_more_small_clusters=false;
    while (!no_more_small_clusters){

      infoMsg(IM::Tiny,"Small cluster removal: # clusters is %d      ",clusters.size());

      if (clusters.size() < 2){
	infoMsg(IM::Tiny,"\n");
	warning("All items ended up in one cluster; MinClusterMembers may not have been satisfied\n");
	no_more_small_clusters=true;
	continue;
      }

      // find the smallest cluster (in terms of number of members)
      int smallest_cluster_members=-1;
      for(ci=clusters.begin();ci!=clusters.end();ci++)
	if (( (ci->items.size() < (unsigned)smallest_cluster_members) || (smallest_cluster_members<0) ) 
	    && (ci->items.size() < commands[command_index].options.MinClusterMembers) ){
	  smallest_cluster_members = ci->items.size();
	  save1=ci;
	}
    
      if(save1==NULL_ITERATOR){
	infoMsg(IM::Tiny,"\nMinClusterMembers threshold reached:  ");
	print_cluster_stats(&clusters,IM::Tiny);
	no_more_small_clusters=true;
	continue;
      }

      saved_size=LZERO;
      for(ci=clusters.begin();ci!=clusters.end();ci++){
	if(ci==save1) continue;

	this_size=merged_sizes[save1->ident][ci->ident];

	if ( (this_size < saved_size) or (saved_size <= LZERO) ){
	  save2=ci;
	  saved_size=this_size;
	}
      }

      if (saved_size <= LZERO)
	error("Something went wrong during part 2 of outlier removal: couldn't find a cluster to merge outlier into");

      // keep save2 - it is likely to be larger than save1
      merge_clusters(&clusters,save2,save1,saved_size);
      // recompute centroid
      if(!set_cluster_centroid(&(*save2),commands[command_index].options.Centroid))
	error("Failed to set_cluster_centroid of newly merged cluster");

      update_merged_sizes_table(&merged_sizes,&clusters,
				commands[command_index].options.ClusterSizeMethod,
				commands[command_index].options.Centroid, save2);

    }
    infoMsg(IM::Tiny,"\n");
    infoMsg(IM::Default,"After small cluster removal:            ");
    print_cluster_stats(&clusters,IM::Default);
  }
    

  // Set the return value: a vector of vectors of parameter names.
  // The actual tying will be done later
  
  infoMsg(IM::Tiny,"Clusters are:\n");

  clustered_params.resize(clusters.size());
  int c_counter=0, m_counter=0;
  for(ci=clusters.begin();ci!=clusters.end();ci++,c_counter++){
    clustered_params[c_counter].resize(ci->items.size());
    infoMsg(IM::Tiny,"cluster %d  #=%d occ=%1.9e size=%1.9e : ",
	    c_counter , ci->items.size(),ci->occupancy,ci->size);


    m_counter=0;
    
    for(std::list<Clusterable*>::iterator ii=ci->items.begin();ii!=ci->items.end();ii++,m_counter++){
      infoMsg(IM::Tiny," %s",(*ii)->name().c_str());
      clustered_params[c_counter][m_counter]=(*ii)->name();
    }
    infoMsg(IM::Tiny,"\n");

  }
  
  
  // to do : clean up
  
  return true;
}

/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::decision_tree_cluster
 *      top-down clustering 
 * 
 * Preconditions:
 *      GM_Parms must be valid
 *
 * Postconditions:
 *      clustered_params contains the clusters (e.g. for passing to a
 *      tying function)
 *
 * Side Effects:
 *      outliers contains the parameters that were removed before clustering 
 *
 * Results:
 *      true if successful; false otherwise
 *
 *-----------------------------------------------------------------------
 */
bool 
GMTK_Tie::decision_tree_cluster(unsigned command_index, 
				std::vector<std::vector<std::string> > &clustered_params,
				std::vector<std::string> &outliers)
{
  if (commands[command_index].command != CT_DTcluster)
    error("decision_tree_cluster was called, but the command is not DTcluster!");
  
  infoMsg(IM::Huge,"Starting decision_tree_cluster\n");


  if (commands[command_index].options.TreeName == ""){
    //commands[command_index].options.TreeName = new_tree_name("ClusteringDT");
    error("No TreeName supplied for command %d",command_index);
    // but the name may still clash with one in the save file?
  }
  
  if (commands[command_index].options.Filename != "")
    error("Cannot load existing trees for the command DTcluster");
  
  
  DecisionTreeType decision_tree;

  std::map<std::string,DecisionTreeType>::iterator ti=decision_trees.find(commands[command_index].options.TreeName);
  if(ti==decision_trees.end()){
    infoMsg(IM::Huge,"No decision tree called %s has been loaded, so creating it\n",commands[command_index].options.TreeName.c_str());
    
    // create the decision tree - max size is 2^(n+1) for n questions we
    // don't allocate all that here, becuase there may be hundreds of
    // possible questions to ask (we just hope that they aren't all
    // needed...)
    //
    //  (unsigned)pow(2,questions.size()+1)
    decision_tree.tree.assign(2, -2); // default initial size (just a root node)
    decision_tree.tree[1]=-3;
    decision_tree.name=commands[command_index].options.TreeName;

    // set up the pointers to a feature set definition and a set of feature values
    if (feature_definition_sets.find(commands[command_index].options.FeatureSetName) == feature_definition_sets.end())
      error("There is no feature definition set called '%s' loaded into memory",commands[command_index].options.FeatureSetName.c_str());
    decision_tree.feature_definitions = &(feature_definition_sets[commands[command_index].options.FeatureSetName]);

    if (feature_value_sets.find(commands[command_index].options.FeatureValuesName) == feature_value_sets.end())
      error("There is no set of feature values called '%s' loaded into memory",commands[command_index].options.FeatureValuesName.c_str());
    decision_tree.feature_values = &(feature_value_sets[commands[command_index].options.FeatureValuesName]);

    if (question_sets.find(commands[command_index].options.QuestionSetName) == question_sets.end())
      error("There is no question set called '%s' loaded into memory",commands[command_index].options.QuestionSetName.c_str());
    decision_tree.questions = &(question_sets[commands[command_index].options.QuestionSetName]);


  } else {
    infoMsg(IM::Huge,"Using existing decision tree called %s (it may or may not be grown further here)\n",commands[command_index].options.TreeName.c_str());

    // for simpilicity, actually extract a copy of the DT (it will be
    // put back at the end of this function)
    decision_tree = decision_trees[commands[command_index].options.TreeName];


  }

  // clear the temporary maps
  decision_tree.nodes_to_clusters.clear();
  decision_tree.nodes_to_param_names.clear();
  
  
  // sanity checks
  //
  // make sure retrieved decision tree is called the expected name
  assert(decision_tree.name == commands[command_index].options.TreeName);
  // command may not have specified the FeatureSetName or
  // FeatureValuesName if it uses a previously loaded tree
  if ( commands[command_index].options.FeatureSetName != "")
    assert(decision_tree.feature_definitions->name == commands[command_index].options.FeatureSetName);
  if ( commands[command_index].options.FeatureValuesName != "")
    assert(decision_tree.feature_values->name == commands[command_index].options.FeatureValuesName);
  if ( commands[command_index].options.QuestionSetName != "")
    assert(decision_tree.questions->name == commands[command_index].options.QuestionSetName);


  std::list<Cluster> clusters;
  std::list<Cluster>::iterator ci,save1=NULL_ITERATOR,save2=NULL_ITERATOR;
  int p_counter=0;


  // make the initial cluster, which contains all items; put outliers to one side
  
  // this method only clusters single digaonal Gaussian Components, so
  // we must extract them if necessary 
  Cluster new_cluster;
  for(std::vector<std::string>::iterator i=commands[command_index].params.begin();i!=commands[command_index].params.end();i++,p_counter++){
    
    Clusterable *c=NULL;
    DiagGaussian *this_dg=NULL;

    switch(commands[command_index].param_type){
      
    case PT_Unknown:
      error("Parameter type is unknown");
      break;
    case PT_Mixture:
      {
	Mixture *this_mixture = find_Mixture(*i);
	Component *this_component = find_Component_of_Mixture(this_mixture);
	if( !is_DiagGaussian(this_component) )
	  error("In decision_tree_cluster, can only handle DiagGaussian components (offending parameter name is %s)",i->c_str());
	this_dg = (DiagGaussian*)this_component;
      }
      break;
    case PT_Component:
      {
	Component *this_component = find_Component(*i);
	if( !is_DiagGaussian(this_component) )
	  error("In decision_tree_cluster, can only handle DiagGaussian components (offending parameter name is %s)",i->c_str());
	this_dg = (DiagGaussian*)this_component;
      }
      break;
    default:
      error("data_driven_cluster is not available for this parameter type");
      break;
      
    }

    c = new ClusterableDiagGaussian(*i, this_dg, commands[command_index].options.DissimilarityMeasure);
    //infoMsg(IM::Max,"Wrappedparameter called %s (underlying Gaussian component is called %s) in ClusterableDiagGaussian\n",i->c_str(),this_dg->name().c_str());

    // set the feature vector for this item by finding the item in the
    // loaded feature vectors
    if (decision_tree.feature_values->values.find(*i) != decision_tree.feature_values->values.end())
      // set the feature value vector for this clusterable item
      c->features = decision_tree.feature_values->values[*i];
    else
      error("No feature vector could be found for the parameter called '%s'\n",i->c_str());


    if (c->occupancy() > commands[command_index].options.ThresholdOccupancyCount){
      new_cluster.occupancy += c->occupancy();
      infoMsg(IM::Max,"Keeping '%s' for building the tree; occ=%e\n",i->c_str(),c->occupancy());
      new_cluster.items.push_back(c);
    } else {
      infoMsg(IM::Max,"Outlier '%s' removed before building tree; occ=%e\n",i->c_str(),c->occupancy());
      outliers.push_back(*i);
    }
  }

  infoMsg(IM::Mod,"%d outliers (occupancy <= ThresholdOccupancyCount) were removed before building tree; %d parameters are left\n",outliers.size(),new_cluster.items.size());

  
  bool finished=false;

  new_cluster.tree_node=1; // the root of the DT
  new_cluster.size = 0.0; // not needed?
  new_cluster.finished=false;


  new_cluster.questions_used.assign(decision_tree.questions->questions.size(),0);


  // push on to list of clusters (this list will grow, as clusters are
  // split)


  double total_ll;
  if(new_cluster.items.size() == 0){
    warning("All parameters were removed as outliers - all parameters will be tied together.\n");
    new_cluster.finished=true;
    new_cluster.occupancy=0.0;
    finished=true;
    total_ll=LZERO;
    decision_tree.tree[1]=-4; // special flag value
  } else {
    infoMsg(IM::Huge,"Computing initial log likelihood\n");
    total_ll= cluster_scaled_log_likelihood(new_cluster.items,&new_cluster.occupancy);
    infoMsg(IM::Mod,"Initial log likelihood before any splitting is %e\n",total_ll);
    infoMsg(IM::Mod,"Initial total occupancy before any splitting is %e\n", new_cluster.occupancy);
  }

  clusters.push_back(new_cluster);
  infoMsg(IM::Huge,"Initialised tree root node with all items\n");

  // recursively split the cluster using the best available question,
  // until the stopping criterion is met
  //recursive_split(clusters,clusters.begin());

  // print some stats
  /*
  double total_ll=0;
  unsigned cii=0;
    for(ci=clusters.begin();ci!=clusters.end();ci++)
    infoMsg(IM::Max,"Adding log likelihood of cluster %d\n",cii);
    cerr << "cluster size is " << ci->items.size() << endl;
    total_ll += cluster_scaled_log_likelihood(ci->items);
  */



  // a non-recursive implementation
  while(!finished){

    infoMsg(IM::Huge,"Starting an iteration\n");
    //print_tree(decision_tree);
    
    // find a cluster (any cluster) that might need further splitting
    for(ci=clusters.begin();ci!=clusters.end();ci++)
      if(ci->finished == false)
	break;
    // now ci is pointing at a cluster to be split, but let's do a
    // sanity check
    if(ci==clusters.end())
      error("All clusters already finished at start of tree-building iteration\n");
    
    // try to find the best question to split this cluster
    int qindex = -1;

    // scores are of the type specified by the user (this feature not
    // yet implementd: currently hardwired to use HTK-tyle liklihood
    // measure)

    
    double this_score= LZERO, best_score=LZERO;
    double best_ll=LZERO, before_ll = cluster_scaled_log_likelihood(ci->items);

    for (unsigned i=0;i<decision_tree.questions->questions.size();i++){
      
      if(question_previously_used(decision_tree,ci->tree_node,i)){
	//infoMsg(IM::Huge,"Tree says that %d is used\n",i);
	continue; // this question was already used

      } else {
	//infoMsg(IM::Huge,"Tree says that %d is unused\n",i);

	infoMsg(IM::Huge,"Question %d : ",i);
      

	// make a temporary split (just in terms of Clusterable*)
	std::list<Clusterable*> left, right;
	double left_occ=0, right_occ=0;
	for(std::list<Clusterable*>::iterator ci2=ci->items.begin();ci2!=ci->items.end();ci2++){
	  if( decision_tree.questions->questions[i].valueSet.find((*ci2)->features[decision_tree.questions->questions[i].feature_index]) != decision_tree.questions->questions[i].valueSet.end()){
	    left.insert(left.end(),(*ci2));
	    left_occ += (*ci2)->occupancy();
	  } else {
	    right.insert(right.end(),(*ci2));
	    right_occ += (*ci2)->occupancy();
	  }
	}

	infoMsg(IM::Huge,"split %d/%d ",left.size(),right.size());

	if (left.empty() or right.empty()){
	  infoMsg(IM::Huge," results in one empty cluster\n");
	  continue;
	}

	if ( (left.size() < commands[command_index].options.MinClusterMembers) or
	     (right.size() <commands[command_index].options.MinClusterMembers )){
	  infoMsg(IM::Huge," results in a small cluster(s)\n");
	  continue;
	}


	if ( (left_occ < commands[command_index].options.MinOccupancyCount) or
	     (right_occ < commands[command_index].options.MinOccupancyCount)){
	  infoMsg(IM::Huge," results in low occupancy cluster(s)\n");
	  continue;
	}


	// score is increase in log likelihood when going from one
	// cluster to two
	double left_ll   = cluster_scaled_log_likelihood(left);
	double right_ll  = cluster_scaled_log_likelihood(right);


	if( (left_ll > LZERO) and (right_ll > LZERO) ){
	  double after_ll  = left_ll + right_ll; 
	  
	  this_score = after_ll - before_ll;
	  
	  infoMsg(IM::Huge,"improvement for this split %e (%3.1f%)\n",this_score,-100*(this_score/total_ll));
	  
	  if (this_score > best_score) {
	    qindex=i;
	    best_score=this_score;
	    best_ll=after_ll;
	  }
	}

      }
    }

    // TO DO: user should be warned if we ran out of questions before
    // reaching some other stopping criterion

    if ( (best_score <= LZERO) or (qindex == -1)) {
      infoMsg(IM::Mod,"   No remaining question could split this branch\n");
      ci->finished=true;
      decision_tree.tree[ci->tree_node] = -1;
      
    } else if  ( (-100*(best_score/total_ll)) < commands[command_index].options.MinImprovementPercent){
	infoMsg(IM::Mod,"   Percent log likelihood improvement (%3.1f%) too small on this branch\n",-100*(best_score/total_ll));
	ci->finished=true;
	decision_tree.tree[ci->tree_node] = -1;

    } else if ( best_score < commands[command_index].options.MinImprovementAbsolute) {
      infoMsg(IM::Mod,"   Absolute log likelihood improvement (%e) too small on this branch\n",best_score);
      ci->finished=true;
      decision_tree.tree[ci->tree_node] = -1;

    } else {

      infoMsg(IM::Mod,"Splitting using question %d, which improves log likelihood by %e (%3.1f%)\n",qindex,best_score,-100*(best_score/total_ll));


      // we found a question, so go ahead and make the split
      
      // ci points at the cluster being split and qindex indicates the
      // best question that can make the split

      // two new clusters formed
      Cluster left, right;
      left.occupancy=0;
      right.occupancy=0;
      left.size=0.0;
      right.size=0.0;
      for(std::list<Clusterable*>::iterator ci2=ci->items.begin();ci2!=ci->items.end();ci2++){
	if( decision_tree.questions->questions[qindex].valueSet.find((*ci2)->features[decision_tree.questions->questions[qindex].feature_index]) 
	    != decision_tree.questions->questions[qindex].valueSet.end()){
	  left.items.insert(left.items.end(),(*ci2));
	  left.occupancy += (*ci2)->occupancy();
	} else {
	  right.items.insert(right.items.end(),(*ci2));
	  right.occupancy += (*ci2)->occupancy();
	}
      }
      

      // this is inefficient - may need to store qorder in the Cluster
      // object?
      unsigned qorder=1;
      for (unsigned i=0;i<ci->questions_used.size();i++)
	if(ci->questions_used[i] != 0)
	  qorder++;
      ci->questions_used[qindex]=qorder;
      left.questions_used = ci->questions_used;
      right.questions_used = ci->questions_used;

      // insert the question into the tree and send the children down
      // to the next level
      decision_tree.tree[ci->tree_node] = qindex;
      left.tree_node=left_child(ci->tree_node);
      right.tree_node=right_child(ci->tree_node);

      grow_decision_tree_if_needed(decision_tree,right.tree_node);

      // mark these nodes as "under contruction" (for debugging
      // purposes only)
      decision_tree.tree[left.tree_node]=-3;
      decision_tree.tree[right.tree_node]=-3;

      left.finished = true;
      right.finished = true;
      for (unsigned i=0;i<decision_tree.questions->questions.size();i++){
	if(left.questions_used[i] == 0)
	  left.finished = false;
	if(right.questions_used[i] == 0)
	  right.finished = false;
      }

      if (left.finished or right.finished)
	warning("All questions were used up for this branch of the decision tree before any other stopping criterion was met: try adding more questions, since some clusters may benefit from further splitting\n");
      

      // old cluster erased
      clusters.erase(ci);

      // new clusters inserted
      clusters.insert(clusters.end(),left);
      clusters.insert(clusters.end(),right);
    }
    
    // are there any clusters that are not finished?
    finished=true;
    for(ci=clusters.begin();ci!=clusters.end();ci++)
      if(!ci->finished)
	finished=false;
    
    // check if target num of clusters is met
    if( clusters.size() >= commands[command_index].options.MaxNumberParameters){
      infoMsg(IM::Huge,"Maximum number of clusters reached\n",qindex);
      finished=true;
    }
      

    /*
    // debug
    int i=0;
    for(ci=clusters.begin();ci!=clusters.end();ci++,i++){
      cerr << "cluster " << i << "/" << clusters.size() << " size=" << ci->items.size()
	   << " node=" << ci->tree_node
	   << " q=";
      for (unsigned i=0;i<ci->questions_used.size();i++)
	cerr << ci->questions_used[i] << " ";
      if (ci->finished)
	cerr << " finished";
      cerr << endl;
    }
    */

    // print some stats
    //double new_total_ll=0;
    total_ll=0;
    for(ci=clusters.begin();ci!=clusters.end();ci++)
      total_ll += cluster_scaled_log_likelihood(ci->items);
    // it is possible that there is no improvement, in cases where a
    // previously unfinished cluster could not be split further
    //if(new_total_ll > total_ll){
    //  infoMsg(IM::Tiny,"Log likelihood improved from %e to %e (%3.1f%)\n",total_ll,new_total_ll,100*(total_ll-new_total_ll)/total_ll);
    //  total_ll=new_total_ll;
    // }

    infoMsg(IM::Mod,"Total log likelihood is now %e\n",total_ll);

    print_cluster_stats(&clusters,IM::Default);


  }

  // Set the return value: a vector of vectors of parameter names.
  // The actual tying will be done later

  //infoMsg(IM::Tiny,"Final tree is:\n");
  //print_tree(decision_tree);

  infoMsg(IM::Mod,"Final clusters are:\n");

  // clusters are ordered so that the decision tree can note which
  // cluster lives at each leaf node

  clustered_params.resize(clusters.size());
  int c_counter=0, m_counter=0;
  for(ci=clusters.begin();ci!=clusters.end();ci++,c_counter++){

    decision_tree.nodes_to_clusters[ci->tree_node]=c_counter;

    // for Mixtures, no need to map back to the Mixture name (because
    // that is what the user requested to be clustered, even though we
    // worked in terms of DiagGaussian Components in this function)
    // because we named the Clusterable objects after the original
    // Mixture name - cunning, eh?
    
    clustered_params[c_counter].resize(ci->items.size());

    infoMsg(IM::Mod," %d : %d members, occ=%1.9e size=%1.9e",
	    c_counter , ci->items.size(),ci->occupancy,ci->size);


    m_counter=0;
    
    for(std::list<Clusterable*>::iterator ii=ci->items.begin();ii!=ci->items.end();ii++,m_counter++){
      infoMsg(IM::Huge," %s",(*ii)->name().c_str());
      clustered_params[c_counter][m_counter]=(*ii)->name();
    }

    infoMsg(IM::Mod,"\n");

  }


  // copy tree back in
  decision_trees[commands[command_index].options.TreeName] = decision_tree;

  // to do : clean up
  return true;
}



/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::decision_tree_synthesise
 *      synthesise (i.e. tie to existing parameters) using a decision tree
 * 
 * Preconditions:
 *      GM_Parms must be valid
 *      Tree must be loaded into memory
 *
 * Postconditions:
 *      parameters specified are tied to existing parameters
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      true if successful; false otherwise
 *
 *-----------------------------------------------------------------------
 */
bool 
GMTK_Tie::decision_tree_synthesise(unsigned command_index)
{
  if (commands[command_index].command != CT_DTsynthesise)
    error("decision_tree_synthesise was called, but the command is not DTsynthesise!");
  
  infoMsg(IM::Huge,"Starting decision_tree_synthesise\n");

  DecisionTreeType *decision_tree = &(decision_trees[commands[command_index].options.TreeName]);
  assert(commands[command_index].options.TreeName == decision_tree->name);


  if (! decision_tree_synthesise(commands[command_index].params.begin(),commands[command_index].params.end(),
				 decision_tree,commands[command_index].options.CollectionName))
    return false;
  
  return true;

}





// synthesise a vector of parameters using the given tree
bool 
GMTK_Tie::decision_tree_synthesise(std::vector<std::string>::iterator b, std::vector<std::string>::iterator e,
				   DecisionTreeType *decision_tree, std::string collection_name)
{
  // we queue up the parameters being tied so that we can execute a
  // single tie_Mixtures per group (of parameters being tied to the
  // same existing parameter) at the end

  // map is indexed by the existing parameter and contains a vector of
  // parameters to be tied to that existing parameter
  map<string,vector<string> > ties;

  for(std::vector<std::string>::iterator pi=b;pi!=e;pi++){

    string parameter_name=*pi;
    
    // we can only synthesis Mixtures if DiagGaussians at the moment
    if (!is_DiagGaussian(find_Component_of_Mixture(parameter_name)))
      error("Cannot synthesise parameter '%s' because it is not a Mixture containing a single diagonal Gaussian",parameter_name.c_str());
    
    // for ease of code readability
    //FeatureDefinitionSetType *feature_definitions = decision_tree->feature_definitions;
    FeatureValueSetType *feature_values = decision_tree->feature_values;
    QuestionSetType *questions = decision_tree->questions;
    
    // get the feature vector for this parameter
    if(feature_values->values.find(parameter_name) == feature_values->values.end())
      error("No feature vector could be found for the parameter called '%s'\n",parameter_name.c_str());
    
    //features are in decision_tree.feature_values->values[*i];
    
    
    // decend the tree until we find a leaf
    unsigned node=1;
    while(decision_tree->tree[node] != -1){ // -1 means a terminal node
      
      // which question is being asked at this tree node?
      int qindex=decision_tree->tree[node];
      
      // what feature is this question about?
      int findex=questions->questions[qindex].feature_index;
      
      // retrieve the value of that feature for the current parameter
      unsigned fval=feature_values->values[parameter_name][findex];
      
      if( questions->questions[qindex].valueSet.find(fval) != questions->questions[qindex].valueSet.end())
	// go left down the tree
	node=left_child(node);
      else
	// go right down the tree
	node=right_child(node);
      
    }
    
    // must now be at a terminal node
    
    if(decision_tree->nodes_to_param_names.find(node) == decision_tree->nodes_to_param_names.end())
      error("There is no parameter at this leaf: this is a bug");
    
    infoMsg(IM::Huge,"Using '%s' to synthesise '%s'\n",decision_tree->nodes_to_param_names[node].c_str(),parameter_name.c_str());

    ties[decision_tree->nodes_to_param_names[node]].push_back(parameter_name);


  }

  /*
  cerr << "Here are the ties about to be done" << endl;
  for(map<string,vector<string> >::iterator ti=ties.begin();ti!=ties.end();ti++){
    cerr << ti->first << " : ";
    for(vector<string>::iterator vi=ti->second.begin();vi!=ti->second.end();vi++)
      cerr << *vi << " ";
    cerr << endl;
  }
  */

  // now execute all the queued-up ties
  for(map<string,vector<string> >::iterator ti=ties.begin();ti!=ties.end();ti++){

    // put existing parameter at start of vector - this is what
    // tie_Mixtures expects
    ti->second.insert(ti->second.begin(), ti->first);

    if(tie_Mixtures(ti->second,CNT_UseExistingCentroid,collection_name,"",false,false) != ti->first)
      error("tie_Mixtures returned an unexpected centroid: this is a bug");

  } 

  //commit_nc_changes(collection_name);


  return true;
}













/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::question_previously_used
 *      determine if a question was asked higher up in the decision tree
 * 
 * Preconditions:
 *      decision_tree must have been constructed
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      true if question has already been asked
 *
 *-----------------------------------------------------------------------
 */

bool
GMTK_Tie::question_previously_used(DecisionTreeType &decision_tree, unsigned node, unsigned question_index)
{

  assert(node>0);
  assert(node<decision_tree.tree.size());

  if ((int)question_index==decision_tree.tree[node])
    return true;
  else if (node > 1)
    return question_previously_used(decision_tree,parent(node),question_index);
  else
    return false;
}

void
GMTK_Tie::print_tree(DecisionTreeType &decision_tree){
  cerr << "DT " << decision_tree.name << ":";
  for(unsigned i=1;i<decision_tree.tree.size();i++)
    cerr << " " << decision_tree.tree[i];
  cerr << endl;

}

void
GMTK_Tie::grow_decision_tree_if_needed(DecisionTreeType &decision_tree, unsigned node){ 

  if (node > decision_tree.tree.size()-1){
    infoMsg(IM::Huge,"Growing decision tree to new size of %d\n",decision_tree.tree.size());
    assert (decision_tree.tree.size() * 2 < decision_tree.tree.max_size());
    //print_tree(decision_tree);
    decision_tree.tree.resize(decision_tree.tree.size() * 2, -2);
    //print_tree(decision_tree);
  }
}

std::string
GMTK_Tie::new_tree_name(const std::string basename)
{
  std::string name;

  unsigned cloneNo=0; 

  do {
    char buff[256];
    sprintf(buff,"%d",cloneNo);
    name = basename + string("_") + buff;
    cloneNo++;
  } while (decision_trees.find(name) != decision_trees.end());

  return name;
}


void
GMTK_Tie::save_decision_tree(unsigned command_index)
{
  infoMsg(IM::Huge,"Saving tree called '%s' in file %s\n",
	  commands[command_index].options.TreeName.c_str(),
	  commands[command_index].options.Filename.c_str());
			 
  if (decision_trees[commands[command_index].options.TreeName].tree.size()<2)
    error("Empty or unitialised decision tree - cannot save");


  // simple current solution: we will not modify existing files
  //
  // simple current file format: 
  //
  // tree_name % DT name
  // number_of_nodes
  // node_index node_type [question_index|parameter_name] % for terminal nodes, give tied param name
  // node_index node_type [question_index|parameter_name] % for terminal nodes, give tied param name
  // node_type can be "nonterminal", "terminal" or "unused"
  // ... etc



  // for ease of code readability
  DecisionTreeType *decision_tree = &(decision_trees[commands[command_index].options.TreeName]);
  FeatureDefinitionSetType *feature_definitions = decision_tree->feature_definitions;
  FeatureValueSetType *feature_values = decision_tree->feature_values;
  QuestionSetType *questions = decision_tree->questions;

  assert(commands[command_index].options.TreeName == decision_tree->name);


  oDataStreamFile os(commands[command_index].options.Filename.c_str(),false,false); // never in binary format, never append


  os.writeString(decision_tree->name,"tree name");
  os.nl();
  os.writeString(feature_definitions->name,"feature definitions name");
  os.nl();
  os.writeString(feature_values->name,"feature values set name");
  os.nl();
  os.writeString(questions->name,"question set name");
  os.nl();
  os.writeUnsigned(decision_tree->tree.size()-1,"number of nodes");
  os.nl();

  for(unsigned i=1;i<decision_tree->tree.size();i++){
    os.nl();
    os.writeUnsigned(i,"node index");
    int qindex=decision_tree->tree[i];
    if( qindex == -1){
      os.writeString("terminal","node type");
      
      
      // cluster index is stored in nodes_to_clusters[i]
      // each cluster has its centroid in first position
      os.writeString(decision_tree->nodes_to_param_names[i],"parameter name");
    
    } else  if( qindex >= 0){
      os.writeString("nonterminal","node type");
      os.writeUnsigned(qindex,"question index");

      std::string qfeat=questions->questions[qindex].feature_name;

      os.writeString(" % Is "+qfeat+" in {","question's feature",false);
      std::string qvals="";
      for(std::set<unsigned>::iterator qsi=questions->questions[qindex].valueSet.begin();qsi!=questions->questions[qindex].valueSet.end();qsi++)
	if (qsi == questions->questions[qindex].valueSet.begin())
	  os.writeString(feature_definitions->features[questions->questions[qindex].feature_index].string_values[*qsi],"first feature value",false);
	else
	  os.writeString(","+feature_definitions->features[questions->questions[qindex].feature_index].string_values[*qsi],"next feature value",false);
      os.writeString("} ?");

    } else  if( qindex == -2){
      os.writeString("unused","node type");
      continue;

    } else
      error("Ill-formed decision tree - cannot save");
    
  }

  os.nl();

}


void
GMTK_Tie::load_decision_tree(unsigned command_index)
{
  // see if there is a same-named tree already loaded
  if (decision_trees.find(commands[command_index].options.TreeName) != decision_trees.end())
    warning("Attempting to load a tree called '%s' when a tree with this name already exists in memory\n",
	    commands[command_index].options.TreeName.c_str());
  
  // load the tree into temporary space
  DecisionTreeType new_tree;

  if (commands[command_index].options.TreeName == "")
    error("No tree name was given in loadTree command - this is a required parameter");

  iDataStreamFile is(commands[command_index].options.Filename.c_str(),false); // never in binary format
  //  oDataStreamFile os("-",false,false); // never in binary format, never append

  //print_tree(decision_trees[commands[command_index].options.TreeName]);



  std::string s;
  is.readString(s,"tree name");
  if(s != commands[command_index].options.TreeName)
    error("File '%s' contains a tree called '%s' but the loadTree command specified a tree called '%s'",
	  commands[command_index].options.Filename.c_str(),
	  s.c_str(), commands[command_index].options.TreeName.c_str());

  new_tree.name=s;
  infoMsg(IM::High,"Found a decision tree called '%s'\n",s.c_str());

  if(!is.readString(s))
    error("Failed to read name of the feature definitions to use for this tree");
  infoMsg(IM::High," which uses feature definition set '%s'\n",s.c_str());
  if (feature_definition_sets.find(s) == feature_definition_sets.end())
    error("There is no feature definition set called '%s' loaded into memory",s.c_str());
  new_tree.feature_definitions=&(feature_definition_sets[s]);


  if(!is.readString(s))
    error("Failed to read name of the feature value set to use for this decision tree");
  infoMsg(IM::High," and uses feature value set '%s'\n",s.c_str());
  if (feature_value_sets.find(s) == feature_value_sets.end())
    error("There is no feature value set called '%s' loaded into memory",s.c_str());
  new_tree.feature_values=&(feature_value_sets[s]);


  if(!is.readString(s))
    error("Failed to read name of the question set to use for this tree");
  infoMsg(IM::High," and uses question set '%s'\n",s.c_str());
  if (question_sets.find(s) == question_sets.end())
    error("There is no question set called '%s' loaded into memory",s.c_str());
  new_tree.questions=&(question_sets[s]);


  unsigned nnodes;
  is.readUnsigned(nnodes,"number of nodes");
  new_tree.tree.resize(nnodes+1);

  for(unsigned i=1;i<=nnodes;i++){
    unsigned ii;
    is.readUnsigned(ii,"node index");
    if(ii!=i)
      error("Nodes out of order when loading decision tree");

    std::string terminal;
    is.readString(terminal,"node type");

    if (terminal == "terminal"){
      // read the name of the parameter at this leaf
      new_tree.tree[i]=-1;
      std::string pname;
      is.readString(pname,"parameter name");
      new_tree.nodes_to_param_names[i]=pname;

    } else if (terminal == "nonterminal") {
      // read the index of the question being asked at this node
      unsigned qindex;
      is.readUnsigned(qindex,"question index");
      new_tree.tree[i]=qindex;

    } else if (terminal == "unused") {
      new_tree.tree[i]=-2;
      continue;

    } else
      error ("Unknown node type '%s' in decision tree file", terminal.c_str());

  }
  


  std::map<std::string,DecisionTreeType>::iterator old_tree;
  if ( (old_tree = decision_trees.find(commands[command_index].options.TreeName)) != decision_trees.end()){
    // check that the loaded tree is an exact match for the one in
    // memory - if so, continue
    bool are_same=true;

    //print_tree(new_tree);
    //print_tree(old_tree->second);

    if(new_tree.name != old_tree->second.name){
      warning("Different tree names\n");
      are_same=false;

    } else if(new_tree.tree.size() != old_tree->second.tree.size()){
      warning("Different size trees\n");
      are_same=false;

    } else if(new_tree.feature_definitions->name != old_tree->second.feature_definitions->name ){
      warning("Different feature definitions are used\n");
      are_same=false;

    } else if(new_tree.feature_values->name != old_tree->second.feature_values->name ){
      warning("Different sets of feature values are used\n");
      are_same=false;

    } else if(new_tree.questions->name != old_tree->second.questions->name ){
      warning("Different sets of questions are used\n");
      are_same=false;

    } else {
      for(unsigned i=1;i<new_tree.tree.size();i++){
	if(new_tree.tree[i] != old_tree->second.tree[i]){
	  warning("Different node type at node %d : %d %d\n",i,new_tree.tree[i],old_tree->second.tree[i]);
	  are_same=false;
	}
	if( (new_tree.tree[i] == -1) and (new_tree.nodes_to_param_names[i] != old_tree->second.nodes_to_param_names[i]  ) ){
	  warning("Different parameter names at node %d : %s %s\n",i,
		  old_tree->second.nodes_to_param_names[i].c_str(),new_tree.nodes_to_param_names[i].c_str());
	  are_same=false;
	}
      }
    }

    if(!are_same)
      error("The tree called '%s' loaded from file '%s' differs from the tree of the same name already in memory",
	    commands[command_index].options.TreeName.c_str(), commands[command_index].options.Filename.c_str());
    else
      warning("OK: the loaded tree exactly matches the tree with this name already in memory\n");

  } else {
    // move tree from temporary space into main storage
    decision_trees[commands[command_index].options.TreeName]=new_tree;
    
  }
}


/*
void
GMTK_Tie::print_feature_maps()
{
  for(unsigned i=0;i<num_features;i++){
    cerr << "Feature " << i << " is called " << feature_names[i] 
	 << " and has index " << feature_name_to_index[feature_names[i]] << endl;
    

    for(FeatureValueType_to_string_mapType::iterator j=FeatureValueType_to_string_maps[i].begin();
	j!=FeatureValueType_to_string_maps[i].end();j++){
      cerr << j->first;
      cerr << "='" << j->second << "'";
      cerr << " '" << FeatureValueType_to_string_maps[i][j->first] << "'=";
      cerr << string_to_FeatureValueType_maps[i][j->second] << endl;
    }

  }


}
*/
