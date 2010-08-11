/*-
 * GMTK_NameCollection.cc
 *     named colleciton.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */



#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
//#include <iostream.h>
#include <algorithm>

#include "general.h"
#include "error.h"

#include "GMTK_NameCollection.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_Mixture.h"


VCID("$Header$")

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * NameCollection::NameCollection()
 *      Constructor
 *
 * Results:
 *      Constructs the object.
 *
 * Side Effects:
 *      None so far.
 *
 *-----------------------------------------------------------------------
 */
NameCollection::NameCollection()
{
 _is_sorted=false;
}


////////////////////////////////////////////////////////////////////
//        I/O
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * NameCollection::read(is)
 *      read in the table.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the internal table
 *
 *-----------------------------------------------------------------------
 */

void
NameCollection::read(iDataStreamFile& is)
{

  NamedObject::read(is);
  int length;
  is.read(length,"Can't read NameCollection's length");

  if (length <= 0) 
    error("ERROR: reading file '%s' line %d, NameCollection '%s', must have positive number of entries.",
	  is.fileName(),is.lineNo(),name().c_str(),length);
  table.resize(length);
  for (int i=0;i<length;i++) {
    is.read(table[i],"Can't read NameCollection's table entry");
  }
}


/*-
 *-----------------------------------------------------------------------
 * NameCollection::write(os)
 *      write out data to file 'os'. 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effects other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
NameCollection::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.write(table.size(),"NameCollection::write length");
  for (unsigned i=0;i<table.size();i++) {
    os.write(table[i],"NameCollection::write table entry");
    if ((i+1) % 10 == 0)
      os.nl();
  }
  os.nl();
}




/*-
 *-----------------------------------------------------------------------
 * NameCollection::fillMxTable
 *      fill the spmf table with entries from GM_Params
 *      This routine must be called to fill the table
 *      that allow a RV to go from an integer index value
 *      to a pointer to the appropriate object.
 *      Note that if mxTable is already filled (size greater 
 *      or equal to table), then nothing happens.
 *
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      fills the mxTable entry.
 *
 *-----------------------------------------------------------------------
 */
void
NameCollection::fillMxTable()
{
  infoMsg(IM::Huge,"NameCollection::fillMxTable\n");

  if (mxTable.size() >= table.size())
    return;

  mxTable.resize(table.size());
  for (unsigned i=0;i<table.size();i++ ) {
    GMParms::ObjectMapType::iterator it;
    if ((it = GM_Parms.mixturesMap.find(table[i])) == GM_Parms.mixturesMap.end())
      error("Error: collection '%s' has named a mixture '%s' (at table entry %d) that doesn't exist.",
	    name().c_str(),table[i].c_str(),i);
    unsigned indx = (*it).second;
    mxTable[i] = GM_Parms.mixtures[indx];
  }
}



/*-
 *-----------------------------------------------------------------------
 * NameCollection::fillspmfTable
 *      fill the mx table with entries from GM_Params.
 *      This routine must be called to fill the table
 *      that allow a RV to go from an integer index value
 *      to a pointer to the appropriate object.
 *      Note that if mxTable is already filled (size greater 
 *      or equal to table), then nothing happens.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      fills the mxTable entry.
 *
 *-----------------------------------------------------------------------
 */
void
NameCollection::fillSpmfTable()
{
  if (spmfTable.size() >= table.size())
    return;

  spmfTable.resize(table.size());
  for (unsigned i=0;i<table.size();i++ ) {
    GMParms::ObjectMapType::iterator it;
    if ((it = GM_Parms.sPmfsMap.find(table[i])) == GM_Parms.sPmfsMap.end())
      error("Error: collection '%s' has named an SPMF '%s' (at table entry %d) that doesn't exist.",
	    name().c_str(),table[i].c_str(),i);
    unsigned indx = (*it).second;
    spmfTable[i] = GM_Parms.sPmfs[indx];
  }
}



/*-
 *-----------------------------------------------------------------------
 * NameCollection::sort
 *      make sorted_table from table
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      fills the sorted_table vector, makes the unsort_mapping and
 *      sets the _is_sorted flag, mxTable and spmfTable are cleared
 *      since it is invalid now
 *
 *-----------------------------------------------------------------------
 */

bool compfn(const pair<string,unsigned> &a, const pair<string,unsigned> &b)
{
  return (a.first < b.first);
}

void
NameCollection::sort()
{

  sorted_table.resize(table.size());
  unsort_mapping.resize(table.size());
  

  vector<pair<string,unsigned> > tmp;
  tmp.resize(table.size());

  for(unsigned i=0;i<table.size();i++){
    tmp[i].first=table[i];
    tmp[i].second=i;
  }

  std::sort(tmp.begin(),tmp.end(), compfn);
  
    
  int j=0;
  for(vector<pair<string,unsigned> >::iterator i=tmp.begin();i!=tmp.end();i++,j++){
    sorted_table[j]=i->first;
    unsort_mapping[j]=i->second;
  }
    
  mxTable.clear();
  spmfTable.clear();

  _is_sorted=true;
}

/*-
 *-----------------------------------------------------------------------
 * NameCollection::resort
 *      re-sort the sorted_table and the associated unsort_mapping
 * 
 * Results:
 *      No results.
 *
 * Preconditions:
 *      sort() must have been called
 *
 * Side Effects:
 *      resorts the sorted_table vector
 *
 *-----------------------------------------------------------------------
 */
void
NameCollection::resort()
{
  if(!_is_sorted)
    error("Cannot call NameCollection::resort() unless the name collection has already been sorted");


  vector<pair<string,unsigned> > tmp;
  tmp.resize(sorted_table.size());

  for(unsigned i=0;i<sorted_table.size();i++){
    tmp[i].first=sorted_table[i];
    tmp[i].second=unsort_mapping[i];
  }

  std::sort(tmp.begin(),tmp.end(), compfn);
  
    
  int j=0;
  for(vector<pair<string,unsigned> >::iterator i=tmp.begin();i!=tmp.end();i++,j++){
    sorted_table[j]=i->first;
    unsort_mapping[j]=i->second;
  }

}



/*-
 *-----------------------------------------------------------------------
 * NameCollection::unsort
 *      copy sorted_table back into table, unsorting it in the process
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      deletes sorted_table and unsort_mapping and unsets the
 *      _is_sorted flag
 *
 *-----------------------------------------------------------------------
 */

void
NameCollection::unsort()
{
  infoMsg(IM::Huge,"NameCollection::unsort (collection name is %s)\n",name().c_str());

  for(unsigned i=0;i<unsort_mapping.size();i++)
    table[unsort_mapping[i]]=sorted_table[i];

  unsort_mapping.clear();  
  sorted_table.clear();  
  _is_sorted=false;

  // must rebuild the tables
  mxTable.clear();
  spmfTable.clear();

  fillMxTable();
  //fillSpmfTable(); call one or the other, but not both ... but how
  //on earth do we know which?

}



void
NameCollection:: search_and_replace(std::vector<std::string> old_names, std::string new_name){

  //cerr << "search_and_replace_in_name_collection for this many old names:" << old_names.size() << endl;
  infoMsg(IM::Huge,"search_and_replace_in_name_collection\n");
  unsigned nchanged=0;

  //cerr << "nc name is " << name() << endl;

  if(! is_sorted())
    sort();

  //cerr << "sorting" << endl;
  std::sort(old_names.begin(),old_names.end());

  vector<std::string>::iterator nci,i,ncie,ie;

  i=old_names.begin();
  ie=old_names.end();

  nci=sorted_table.begin();
  ncie=sorted_table.end();

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

  resort();

  infoMsg(IM::Huge,"search_and_replace changed %d entries\n",nchanged);

}

void
NameCollection::queue_search_and_replace(std::vector<std::string> old_names, std::string new_name){

  queued_changes.reserve(queued_changes.size() + old_names.size());

  // just note a set of changes, without actually making them just yet
  for (vector<std::string>::iterator i=old_names.begin();i!=old_names.end();i++){
    pair<string,string> c;
    c.first=*i;
    c.second=new_name;
    queued_changes.push_back(c);
  }


}

bool compfn2(const pair<string,string> &a, const pair<string,string> &b)
{
  return (a.first < b.first);
}

void
NameCollection::commit_all_searches_and_replacements(){

  infoMsg(IM::Huge,"Committing all searches and replacements to name collection %s: ",name().c_str());

  // execute searches in sort order
  std::sort(queued_changes.begin(),queued_changes.end(),compfn2);

  // actually execute all the queued-up search/replace operations in
  // such a way that we do not need to resort the sorted_table until
  // the end

  // note: still leaves collection in "sorted" state (i.e. must call
  // unsort() later)


  if(! is_sorted())
    sort();  



  vector<pair<string,string> >::iterator qi=queued_changes.begin(), qie=queued_changes.end();
  vector<std::string>::iterator nci=sorted_table.begin(), ncie=sorted_table.end();
  unsigned nchanged=0;

  //cerr << "running" << endl; 

  
  for(vector<pair<string,string> >::iterator qi=queued_changes.begin(); qi!=queued_changes.end(); qi++){

    //cerr << "processing " << qi->first << endl;

    // jump to first match and change it and all contiguous matching
    // entries in sorted_table
    nci=std::find(nci,ncie,qi->first);
    
    if (nci==ncie)
      error("\n\nThere was a queued change of %s -> %s, but the named collection '%s' does not contain %s. This is probably because an earlier command has already changed this entry in the name collection.\n",
	    qi->first.c_str(),qi->second.c_str(), name().c_str(), qi->first.c_str());

    while ( (nci!=ncie) and (*nci == qi->first) ){
      //cerr << "replacing " << *nci << " with " << qi->second << endl;
      *nci = qi->second;
      nchanged++;
      nci++;
    }
   
    if (nci==ncie)
      break;
  }

  resort();
  queued_changes.clear();

  infoMsg(IM::Huge," changed %d entries\n",nchanged);


  
  
}

