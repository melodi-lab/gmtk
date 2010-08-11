/*-
 * GMTK_NameCollection.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */

/*
 * This class is just a readable/writable table of names of objects
 * (such as Mixtures, SPMFS, and so on).
 * Instances of this class are used to associate integer decision tree leaves
 * to actual object pointers.
 */


#ifndef GMTK_NAMECOLLECTION_H
#define GMTK_NAMECOLLECTION_H

#include <vector>

#include "fileParser.h"
#include "logp.h"

#include "GMTK_NamedObject.h"
#include "GMTK_GMParms.h"

class MixtureCommon;
class Sparse1DPMF;

class NameCollection : public NamedObject  {


  friend class GMParms;
  friend class RV;
  friend class DiscRV;
  friend class FileParser;
  friend class GMTK_Tie;
  friend void search_and_replace_in_name_collection(NameCollection *nc, std::vector<std::string> old_names, std::string new_name);

  // Possible instantiations of this class
  // 0) just table is allocated
  //    (right after being read in from disk)
  // 1) table is allocated
  //    one or both of mxTable and spmfTable are allocated
  //    (after being associated with an object)
  // 2) table is not allocated
  //    one or both of mxTable and spmfTable are allocated
  //    (special global objects)

  // string of names. 
  vector<string> table;

  // a sorted version of table
  vector<string> sorted_table;

  // changes that are waiting to be made to sorted_table
  vector< pair<string,string> > queued_changes;

  // if this is true, then sorted_table has been made from table, and
  // may have subsequently been modified, so we should not use table
  // again until unsort() has been called
  bool _is_sorted;
  // and a place to remember how to unsort it again
  vector<unsigned> unsort_mapping;
  
  // direct pointers to those objects
  // for which this might refer to.
  vector<Mixture*> mxTable;
  vector<Sparse1DPMF*> spmfTable;


public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  NameCollection();
  ~NameCollection() { }

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  //////////////////////////////////////////////
  // routines to fill in tables
  void fillMxTable();
  void fillSpmfTable();

  // access routines to mixtures
  Mixture* mx(int i) { return mxTable[i]; }
  unsigned mxSize() { return mxTable.size(); }
  bool validMxIndex(unsigned u) { return (u < mxSize()); }

  // access routines to sparse 1D PMFs
  Sparse1DPMF* spmf(int i) { return spmfTable[i]; }
  unsigned spmfSize() { return spmfTable.size(); }
  bool validSpmfIndex(unsigned u) { return (u < spmfSize()); }

  //////////////////////////////////////////////
  // sorting and unsorting, replacing, etc
  // to do: clear sorted_table on loading, etc
  inline bool is_sorted() { return _is_sorted; };
  void sort();
  void unsort();
  void resort();
  void search_and_replace(std::vector<std::string> old_names, std::string new_name);
  void queue_search_and_replace(std::vector<std::string> old_names, std::string new_name);
  void commit_all_searches_and_replacements();
};


#endif // defined NAMECOLLECTION
