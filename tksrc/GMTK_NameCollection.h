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
 * (such as Gaussian Mixtures, SPMFS, and so on).
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

class MixGaussiansCommon;
class Sparse1DPMF;

class NameCollection : public NamedObject  {

  friend class GMParms;

  // Possible instantiations of this class
  // 0) just table is allocated
  //    (right after being read in from disk)
  // 1) table is allocated
  //    one or both of mgTable and spmfTable are allocated
  //    (after being associated with an object)
  // 2) table is not allocated
  //    one or both of mgTable and spmfTable are allocated
  //    (special global objects)

  // string of names. 
  vector<string> table;
  // direct pointers to those objects
  // for which this might refer to.
  vector<MixGaussians*> mgTable;
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
  void fillMgTable();
  void fillSpmfTable();

  // access routines to mix Gaussians
  MixGaussians* mg(int i) { return mgTable[i]; }
  unsigned mgSize() { return mgTable.size(); }
  bool validMgIndex(unsigned u) { return (u < mgSize()); }

  // access routines to sparse 1D PMFs
  Sparse1DPMF* spmf(int i) { return spmfTable[i]; }
  unsigned spmfSize() { return spmfTable.size(); }
  bool validSpmfIndex(unsigned u) { return (u < spmfSize()); }

};


#endif // defined NAMECOLLECTION
