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

VCID("$Header$")



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
GMTK_Tie::GMTK_Tie(GMParms *GM_Parms_ptr)
{
  GM_Parms=GM_Parms_ptr;

  ////////////////////////////////////////////////////////////////////////
  // make maps between tying command file keywords and the enumerated
  // types

  CommandType_to_string[CT_Unknown]="unknown";
  CommandType_to_string[CT_Tie]="tie";
  CommandType_to_string[CT_Cluster]="cluster";
  CommandType_to_string[CT_DTcluster]="DTcluster";
  CommandType_to_string[CT_Untie]="untie";

  ParameterType_to_string[PT_Unknown]="Unknown";
  ParameterType_to_string[PT_Mixture]="Mixture";
  ParameterType_to_string[PT_Mean]="Mean";
  ParameterType_to_string[PT_Collection]="Collection";

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

    // the parameter type on which it operates
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


    // the list of optional arguments, each of the form key=value
    if (!is.readUnsigned(noptions)) return false;

    fill_in_default_values(commands[i].command, commands[i].options);

    for (unsigned j=0;j<noptions;j++){

      char* opt;
      if (!is.readStr(opt)) return false;
      parse_option(commands[i].command, commands[i].options, std::string(opt));

    }

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
  for (CommandListType::iterator i=commands.begin();i!=commands.end();i++)
    if ( ! (i->options.MinOccupancyCount == 0.0) ) // 0.0 is special flag value
      return true;
  
  return false;
}



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
  default:
    error("Command type %s is not yet supported", CommandType_to_string[ctype].c_str());
  }



}



/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::enter_tie_option
 *      parses an option that applies to the command "tie"
 * 
 * Preconditions:
 *      a 'tie' command must have been encountered in the command
 *      file, and a subsequent option read in
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
GMTK_Tie::enter_tie_option(GMTK_Tie::Options &opt, std::string &key, std::string &value)
{
  if  (key == "Centroid"){
    map<std::string, CentroidType>::iterator i = string_to_CentroidType.find(value);
    if(i != string_to_CentroidType.end())
      opt.Centroid=i->second;
    else
      error("Invalid value of '%s' for option %s",value.c_str(),key.c_str());
  } else
    error("The 'tie' command does not accept the option '%s'",key.c_str());
}



/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::enter_cluster_option
 *      parses an option that applies to the command "cluster"
 *      
 * 
 * Preconditions:
 *      a 'cluster' command must have been encountered in the command
 *      file, and a subsequent option read in
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
GMTK_Tie::enter_cluster_option(GMTK_Tie::Options &opt, std::string &key, std::string &value)
{

  // the strings are hardwired here - this is a bad idea
  // to do: fix it

  if  (key == "MaxClusterSize"){
    opt.MaxClusterSize=atof(value.c_str());

  } else if  (key == "MinClusterMembers"){
    opt.MinClusterMembers=atoi(value.c_str());

  } else if  (key == "MinOccupancyCount"){
    opt.MinOccupancyCount=(double)atof(value.c_str());

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

  } else if  (key == "ClusterSizeMethod"){
    map<std::string, ClusterSizeMethodType>::iterator i = string_to_ClusterSizeMethodType.find(value);
    if(i != string_to_ClusterSizeMethodType.end())
      opt.ClusterSizeMethod=i->second;
    else
      error("Invalid value of '%s' for option %s",value.c_str(),key.c_str());


  } else
    error("The 'cluster' command does not accept the option '%s'",key.c_str());

}


/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::enter_DTcluster_option
 *      parses an option that applies to the command "DTcluster"
 * 
 * Preconditions:
 *      a 'DTcluster' command must have been encountered in the
 *      command file, and a subsequent option read in
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
GMTK_Tie::enter_DTcluster_option(GMTK_Tie::Options &opt, std::string &key, std::string &value)
{
  // not implemented yet
  error("The 'DTcluster' command does not accept the option '%s'",key.c_str());
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
  opt.MaxClusterSize = -1.0;
  opt.MinClusterMembers  = -1;
  opt.MinOccupancyCount = 0.0;
  opt.NewParameterNamePrefix = "";
  opt.Centroid = CNT_Arbitrary;

  switch(ctype){
  case CT_Tie:
  case CT_Untie:
    break;
  case CT_Cluster:
    opt.MaxClusterSize = 20.0; 
    opt.MinClusterMembers = 1;
    opt.DissimilarityMeasure = DMT_EmulateHTK;
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
    case PT_Mean:
      expand_param_names(&(GM_Parms->meansMap), i->param_expression, &(i->params));
      break;
    case PT_Collection:
	error("Collection is not yet a supported parameter type");
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
GMTK_Tie::print_Command(Command &cmd, int index, IM::VerbosityLevels v)
{

  infoMsg(v,"Command %d\n",index);
  infoMsg(v," CommandType: %s\n",CommandType_to_string[cmd.command].c_str());
  infoMsg(v," ParameterType: %s\n",ParameterType_to_string[cmd.param_type].c_str());
  print_Options(cmd.options,v);

  infoMsg(v+IM::Increment," Matching parameters: ");
  for(std::vector<std::string>::iterator i=cmd.params.begin();i!=cmd.params.end();i++)
    infoMsg(v+IM::Increment," %s", i->c_str());
  infoMsg(v,"\n");
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

  infoMsg(v,"  MaxClusterSize: %1.9e\n",opt.MaxClusterSize);
  infoMsg(v,"  MinClusterMembers: %d\n",opt.MinClusterMembers);
  infoMsg(v,"  NewParameterNamePrefix: %s\n",opt.NewParameterNamePrefix.c_str());
  infoMsg(v,"  DissimilarityMeasure: %s\n",DissimilarityMeasureType_to_string[opt.DissimilarityMeasure].c_str());
  infoMsg(v,"  ClusterSizeMethod: %s\n",ClusterSizeMethodType_to_string[opt.ClusterSizeMethod].c_str());
  infoMsg(v,"  Centroid: %s\n",CentroidType_to_string[opt.Centroid].c_str());

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
  print_Command(commands[command_index],command_index,IM::Default);

  if (commands[command_index].params.size() < 1){
    warning("There are no parameters matching this command - doing nothing!");
    return true;
  }

  bool rval=true;

  switch(commands[command_index].command){

  case GMTK_Tie::CT_Unknown:
    error("Command is unknown");
    break;

  case GMTK_Tie::CT_Tie:
    
    switch(commands[command_index].param_type) {
    case GMTK_Tie::PT_Mixture:
      if (!(tie_Mixtures(commands[command_index].params, commands[command_index].options.Centroid,false)))
	error(" Failed to tie Mixtures");
      break;

    case GMTK_Tie::PT_Mean:
      warning("Tying Means is not yet implemented; model not changed");
      break;

    case GMTK_Tie::PT_Collection:
      warning("Tying Collections is not yet implemented; model not changed");
      break;

    default:
      warning("Unsupported parameter type for command 'tie'. Model unchanged");
    }
    break;



  case GMTK_Tie::CT_Cluster:
    {
      std::vector<std::vector<std::string> > clustered_params;
      if(!data_driven_cluster(command_index, clustered_params))
	error(" Failed to do data_driven_cluster");
      
      switch(commands[command_index].param_type) {

      case GMTK_Tie::PT_Mean:
	{
	  for(unsigned j=0;j<clustered_params.size();j++){ // for each cluster
	     infoMsg(IM::Mega,"Tying %d Means within cluster %d\n",clustered_params[j].size(), j);
	    if (!(tie_Means(clustered_params[j], commands[command_index].options.Centroid,false))) 
	      error(" Failed to tie Means after clustering");
	  }
	}
	break;

      case GMTK_Tie::PT_Mixture:
	{
	  for(unsigned j=0;j<clustered_params.size();j++){ // for each cluster
	     infoMsg(IM::Mega,"Tying %d Mixtures within cluster %d\n",clustered_params[j].size(), j);
	    if (!(tie_Mixtures(clustered_params[j], commands[command_index].options.Centroid,false))) 
	      error(" Failed to tie Mixtures after clustering");
	  }
	}
	break;

      default:
	error("Unsupported parameter type for tying after clustering");
      }
    }    
    break;
    
    
  case GMTK_Tie::CT_DTcluster:
    warning("DTcluster not yet implemented");
    break;



  case GMTK_Tie::CT_Untie:
    
    switch(commands[command_index].param_type) {
    case GMTK_Tie::PT_Mixture:
      //the params are expanded by find_matching_command_parameters().
      if (!(untie_Mixtures(commands[command_index].params,false)))
	error(" Failed to untie Mixtures");
      break;

    case GMTK_Tie::PT_Mean:
      warning("Untying Means is not yet implemented; model not changed");
      break;

    case GMTK_Tie::PT_Collection:
      warning("Untying Collections is not yet implemented; model not changed");
      break;

    default:
      warning("Unsupported parameter type for command 'untie'. Model unchanged");
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
bool 
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
    warning("Tying with no parameters makes no sense - doing nothing for this command!");
    return true;
  }

  if ( params.size() == 1){
    warning("Tying a single parameter has no effect; model is unchanged.");
    return true;
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
	    warning("MeanVector '%s' will be tied, although could not find any Components that use it",i->c_str());

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
  GM_Parms->markUsedMixtureComponents(); // does this actually delete anything?
  // and any now-unused  Dense1DPMF from GM_Params
  // TO DO - no function currently exists for this in class GMParms
  // TO DO - other param types too?
    

  return true;
}




/*-
 *-----------------------------------------------------------------------
 * GMTK_Tie::tie_Mixtures
 *      tie the underlying Components and weights of the specified
 *      Mixtures
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
bool 
GMTK_Tie::tie_Mixtures(std::vector<std::string> param_expressions, CentroidType method, bool expand_expressions)
{

  std::vector<std::string>::iterator i;
  std::vector<string> params;

  if(expand_expressions)
    // expand out the list of regex, param_expressions, into actual
    // parameter names. 
    for(i=param_expressions.begin();i!=param_expressions.end();i++)
      expand_param_names(&GM_Parms->mixturesMap,*i,&params);
  else
    params=param_expressions;

  if ( params.size() == 0){
    warning("Tying with no parameters makes no sense - doing nothing for this command!");
    return true;
  }

  if ( params.size() == 1){
    if(expand_expressions)
      warning("Tying a single parameter has no effect; model is unchanged.");
    return true;
  }

  // first, locate the Mixtures concerned in the list of Mixtures and
  // find all their mixture weights and sets of Gaussian components

  
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
      return false;
    }

  }

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


  infoMsg(IM::Mod,"tie_Mixtures to a Centroid selected by method '%s'\n",CentroidType_to_string[method].c_str());

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
      Mixture* tied_mixture = find_Mixture(*i); // GM_Parms.mixtures[GM_Parms.mixturesMap[*i]];
      
      i++; // now pointing to second param
      for(;i!=params.end();i++){
	
	Mixture* this_mixture = find_Mixture(*i); // GM_Parms.mixtures[GM_Parms.mixturesMap[*i]];
	
	infoMsg(IM::Low,"Tying %s to %s\n",i->c_str(),params.begin()->c_str());
	
	// delete existing Components from the Mixture object
	this_mixture->components.clear();
	
	// insert tied components (works for multi-component mixtures too)
	this_mixture->components.assign(tied_mixture->components.begin(),tied_mixture->components.end());
	this_mixture->dense1DPMF = tied_mixture->dense1DPMF ;
      }
    }
    break;



  case CNT_CreateCentroid_averageSingleComponentMixtures:
    {
      // only for single-component Mixtures of Gaussians
      i=params.begin();
      Mixture *this_mixture,*tied_mixture;
      for(;i!=params.end();i++){
	this_mixture = find_Mixture(*i); // GM_Parms.mixtures[GM_Parms.mixturesMap[*i]];
	if (this_mixture->components.size() != 1)
	  error ("CreateCentroid_averageSingleComponentMixtures only available for Mixtures with exactly one Component");
	//if( std::string(typeid(*(this_mixture->components[0])).name()).find("DiagGaussian") == std::string::npos)
	if( this_mixture->components[0]->name() != "Diag Gaussian")
	  error("CreateCentroid_averageSingleComponentMixtures only available for DiagGaussian Components");
	
      }
      
      // keep first Mixture's mean/var/weights, sum all other components
      // into them, then divide by M.
      i=params.begin();
      tied_mixture = find_Mixture(*i); //GM_Parms.mixtures[GM_What must be true before the function is called.Parms.mixturesMap[*i]];
      DiagGaussian *tied_gaussian = (DiagGaussian*)tied_mixture->components[0];
      
      i=params.begin(); i++; // now pointing to second param
      for(;i!=params.end();i++){
	
	Mixture* this_mixture = find_Mixture(*i); //GM_Parms.mixtures[GM_Parms.mixturesMap[*i]];
	
	infoMsg(IM::Low,"Tying %s to %s\n",i->c_str(),params.begin()->c_str());
	
	// sum parameters into first component. Since these are
	// DiagGaussian components, we merge mean and covar
	DiagGaussian *this_gaussian = (DiagGaussian*)this_mixture->components[0];
	
	
	// hmmm - this seems awfully low-level
	for (int k=0; k<tied_gaussian->mean->means.len(); k++){
	  tied_gaussian->mean->means[k] += this_gaussian->mean->means[k];
	  tied_gaussian->covar->covariances[k] += this_gaussian->covar->covariances[k];
	}
	
	
	// next line should be redundant (weight vector should be "1.0"
	// for all Mixtures)
	for (int k=0; k<tied_mixture->dense1DPMF->pmf.len(); k++)
	  tied_mixture->dense1DPMF->pmf[k] += this_mixture->dense1DPMF->pmf[k]; 
	
	// delete existing Components from the Mixture object 
	this_mixture->components.clear();
	
	// second, insert tied component (this is just setting pointers,
	// so can be done even while values of params pointed at are
	// being changed)
	this_mixture->components.assign(tied_mixture->components.begin(),tied_mixture->components.end());
	this_mixture->dense1DPMF = tied_mixture->dense1DPMF;
      }
    
      
      float M=(float)params.size();
      for (int k=0; k<tied_gaussian->mean->means.len(); k++){
	tied_gaussian->mean->means[k] /= M;
	tied_gaussian->covar->covariances[k] /= M;
	
      }
      for (int k=0; k<tied_mixture->dense1DPMF->pmf.len(); k++)
	tied_mixture->dense1DPMF->pmf[k] /= M;
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
    error("EmulateHTKMixturesOfGaussians not yet implemented for tying Mixtures");
    break;

  default:
    error("'%s' is not a valid Centroid selection method for tying Mixtures",CentroidType_to_string[method].c_str());
    break;


  }

  // Clean up
  // --------

  //  delete any now-unused  Components from GM_Params
  GM_Parms->markUsedMixtureComponents(); // does this actually delete anything?
  // and any now-unused  Dense1DPMF from GM_Params
  // TO DO - no function currently exists for this in class GMParms
    

  return true;
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
GMTK_Tie::untie_Mixtures(std::vector<std::string> param_expressions, bool expand_expressions)
{

  // this should work for any type of component - to do: remove any
  // previous checks for Gaussian/DiagGaussian

  // to do - share first part of this function with tie_Mixtures

  std::vector<std::string>::iterator i;
  std::vector<string> params;

  if(expand_expressions)
    // expand out the list of regex, param_expressions, into actual
    // parameter names. 
    for(i=param_expressions.begin();i!=param_expressions.end();i++)
      expand_param_names(&GM_Parms->mixturesMap,*i,&params);
  else
    params=param_expressions;

  if ( params.size() == 0){
    warning("Untying with no parameters makes no sense - doing nothing for this command!");
    return true;
  }

  /*

  this is not true - can untie a single param - it just gets its own copy (in case it was tied to something)

  if ( params.size() == 1){
    if(expand_expressions)
      warning("Untying a single parameter has no effect; model is unchanged.");
    return true;
  }
  */

  // locate the Mixtures concerned in the list of Mixtures and find
  // all their mixture weights and sets of components


  // make an independent copy of Components and weights for each Mixture

  // NOTE: this assumed Mixture tying is ONLY EVER done by tying
  // underlying Components and weights!!!!!

  // we don't care whether this leaves orphaned (unused) parameters -
  // they will get cleaned up somewhere else later; they will remain
  // in GM_Parms for now, so we can find them (no memory leak)
  
  for(i=params.begin();i!=params.end();i++){
    Mixture* this_mixture = find_Mixture(*i);

    for(vector< Component* >::iterator j=this_mixture->components.begin();j!=this_mixture->components.end();j++){

      // clone this component and use the clone to replace the
      // (possibly shared) original; the identicalIndependentClone
      // function takes care of generating a new name and putting the
      // new object into the GM_Parms table

      // the old Component (which we about to lose the pointer to)
      // needs to decrease its usage by 1
      (*j)->adjustNumTimesShared(-1);

      *j = (*j)->identicalIndependentClone();

      // the new Component we just made now gets a usage count of 1
      (*j)->adjustNumTimesShared(1);

      // now adjust the usage counts of any subobjects that have them


   
    }

    // Dense1DPMFs do not have numTimesShared
    this_mixture->dense1DPMF = this_mixture->dense1DPMF->identicalIndependentClone();

  }

  // Clean up
  // --------

  //  delete any now-unused  Components from GM_Params
  GM_Parms->markUsedMixtureComponents(); // does this actually delete anything?
  // and any now-unused  Dense1DPMF from GM_Params
  // TO DO - no function currently exists for this in class GMParms
    

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

  infoMsg(v,"%d clusters with occupancies %1.9e to %1.9e and # members %d to %d\n",
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
 *      none (not actual tying takes place here)
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
    error("data_driven_cluster_Mixtures was called, but the command is not cluster!");

  // simple bottom-up clustering using only the values of the parameters themselves
  
  // TO DO: insist on homogeneous parameter type within the set of
  // parameters being clustered


  std::list<Cluster> clusters;
  std::list<Cluster>::iterator ci;
  std::list<Cluster>::iterator ci2;
  std::list<Cluster>::iterator save1;
  std::list<Cluster>::iterator save2;

  float this_size, saved_size=LZERO;
  int p_counter=0;
  // make the initial set of clusters - one per parameter
  for(std::vector<std::string>::iterator i=commands[command_index].params.begin();i!=commands[command_index].params.end();i++,p_counter++){
    
    Clusterable *c=NULL;
    
    switch(commands[command_index].param_type){
      
    case GMTK_Tie::PT_Unknown:
      error("Parameter type is unknown");
      break;
    case PT_Mixture:
      {
	Mixture *this_mixture = find_Mixture(*i);
	c = new ClusterableMixture(*i, this_mixture, commands[command_index].options.DissimilarityMeasure);
      }
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
    warning("Clustering with fewer than 2 parameters makes no sense - doing nothing for this command!");
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
	warning("All items ended up in one cluster; MinOccupancyCount may not have been satisfied");
	no_more_low_occupancy=true;
	continue;
      }

      // find the lowest occupancy cluster
      save1=clusters.end();
    
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
	save1=clusters.end();
    
      if(save1==clusters.end()){
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
	warning("All items ended up in one cluster; MinClusterMembers may not have been satisfied");
	no_more_small_clusters=true;
	continue;
      }

      // find the smallest cluster (in terms of number of members)
      int smallest_cluster_members=-1;
      for(ci=clusters.begin();ci!=clusters.end();ci++)
	if (( (ci->items.size() < (unsigned)smallest_cluster_members) || (smallest_cluster_members<0) ) 
	    && (ci->items.size() < (unsigned)commands[command_index].options.MinClusterMembers) ){
	  smallest_cluster_members = ci->items.size();
	  save1=ci;
	}
    
      if(save1==clusters.end()){
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

