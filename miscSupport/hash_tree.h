/*
 * hash_tree.h
 *   Hash Tree, for hashing clique values into a tree.
 *   This structure is a tree, where each node is a linked list hash table .
 *   The structure supports:
 *        - rapid hash'ed search, insertion, etc. for entries
 *        - iterators through the tree, where the path from
 *          root to leaf is stored in external variables
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
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

#ifndef HASHTREE_H
#define HASHTREE_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "hash_map_list.h"


template <class _Key, class _Data> 
class hash_tree : public hash_abstract {

  /////////////////////////////////////////////////////////
  // the depth of the tree, in that there are this
  // many tree levels, where each tree level corresponds
  // to one possible RV value.
  const unsigned treeDepth;

  /////////////////////////////////////////////////////////
  // an array of length treeDepth holding the default hash-table
  // starting sizes.  Note: this is a pointer to an external array
  // managed external to this object, so we don't change anything
  // here. If this is NULL, we ignore it.
  const unsigned* const default_hash_start_sizes;

public:

  ////////////////////
  // constructor
  // TODO: need to add a list of default starting sizes for each hash
  // level.
  hash_tree(unsigned arg_depth,
	    unsigned *arg_dfss) : 
    treeDepth(arg_depth),
    default_hash_start_sizes(arg_dfss) {}

  ////////////////////////////////
  // abstract class for data items.
  class DataItem {
  protected:
    DataItem() {}
  public:
  };

  //////////////////////////////////////////////////////
  // A DataItem for the terminal nodes in the tree. This node contains
  // a true item values.
  class TerminalDataItem : public DataItem {
    friend class hash_tree<_Key,_Data>;
    friend class hash_tree<_Key,_Data>::iterator;
    friend class hash_tree<_Key,_Data>::iterator_vector;
    friend class hash_tree<_Key,_Data>::iterator_vectorp;
    _Data val;
  public:
    TerminalDataItem() {}
  };


  //////////////////////////////////////////////////////
  // A DataItem for the non-terminal nodes in the tree. This node
  // contains a hash-table itself to hash down to the next level
  // DataItems in the tree. At the penultimate level in the tree, this
  // hash table hashes down to TerminalDataItems (which are at the
  // lowest tree level), and at all levels other than the penultimate
  // level, this hash table hashes down to other NonTerminalDataItems.
  class NonTerminalDataItem : public DataItem {
    friend class hash_tree<_Key,_Data>;
    friend class hash_tree<_Key,_Data>::iterator;
    friend class hash_tree<_Key,_Data>::iterator_vector;
    friend class hash_tree<_Key,_Data>::iterator_vectorp;
    hash_map_list< _Key, DataItem* > base_map;
  public:
    NonTerminalDataItem(const unsigned hash_start_size)
      : base_map(hash_start_size) {}
  };

  //////////////////////////////////////////
  // The root level of the tree, the hash table at the root.
  hash_map_list < _Key, DataItem* > base_map;

  ///////////////////////////////////////////////////////
  // insert an item <key vector,val> into the hash table.
  // Return a pointer to the data item in the hash table
  // after it has been inserted. This way, it is possible
  // to insert an empty item (using a _Data() constructor)
  // and then adjust its value later via the returned pointer.
  // The foundp argument is set to true when the key has
  // been found. In this case (when it was found), the existing
  // _Data item is not adjusted to be val.
  _Data* insert( _Key*key_vec, 
		 _Data val,
		 bool&foundp = hash_abstract::global_foundp) {

    unsigned depth = 0;
    NonTerminalDataItem* ntdi;
    // createRest indicates if we need to by default
    // create the rest of the structure during insertion. ON
    // the first not-found, we set this to true. We could
    // also have a variable found = !createRest, but that
    // would be redundant.
    bool createRest = false;

    _Key key = key_vec[depth];
    bool foundp_at_level;
    DataItem** valp = base_map.insert(key,NULL,foundp_at_level);
    if (!foundp_at_level) {
      // key was inserted but not found, meaning it was not there
      // before, so we need to create a new object.
      ntdi = new NonTerminalDataItem(default_hash_start_sizes[depth]);
      *valp = ntdi;
      createRest = true;
    } else 
      ntdi = (NonTerminalDataItem*)(*valp);
  
    for (depth=1;depth<(treeDepth-1);depth++) {
      if (createRest) {
	NonTerminalDataItem* next_ntdi = 
	  new NonTerminalDataItem(default_hash_start_sizes[depth]);
	ntdi->base_map[key_vec[depth]] = next_ntdi;
	ntdi = next_ntdi;
      } else {

	key = key_vec[depth];
	valp = ntdi->base_map.insert(key,NULL,foundp_at_level);
	if (!foundp_at_level) {
	  // not there so we need to create a new object.
	  NonTerminalDataItem* next_ntdi = 
	    new NonTerminalDataItem(default_hash_start_sizes[depth]);
	  *valp = next_ntdi;
	  createRest = true;
	  ntdi = next_ntdi;
	} else {
	  NonTerminalDataItem* next_ntdi = 
	    (NonTerminalDataItem*)(*valp);
	  ntdi = next_ntdi;
	}
      }
    }
    TerminalDataItem* tdi;
    if (createRest) {
      tdi = new TerminalDataItem();
      ntdi->base_map[key_vec[depth]] = tdi;
    } else {
      key = key_vec[depth];
      valp = ntdi->base_map.insert(key,NULL,foundp_at_level);
      if (!foundp_at_level) {
	tdi = new TerminalDataItem();
	*valp = tdi;
	createRest = true;
      } else {
	// this probably shouldn't happen, since it means
	// we are inserting a different val for a key that
	// was already inserted. Might want to add an
	// assert(0) here or some notification.
	tdi = (TerminalDataItem*)(*valp);
      }
    }
    foundp = !createRest;
    if (!foundp)
      numberUniqueEntriesInserted++;
    assert ( ++depth == treeDepth );
    tdi->val = val;
    return &(tdi->val);
  }

  unsigned numEntries() {
    return numberUniqueEntriesInserted;
  }  

  ///////////////////////////////////////////////////
  // insert() but using an sArray()
  _Data* insert(sArray<_Key>& key_vec, 
		_Data val,
		bool&foundp = hash_abstract::global_foundp)
  {
    assert (key_vec.size() == treeDepth);
    return insert(key_vec.ptr,val,foundp);
  }


  ////////////////////////////////////////////////////////
  // search for key returning true if the key is found, otherwise
  // don't change the table. Return pointer to the data item
  // when found, and NULL when not found.
  _Data* find( sArray<_Key>& key_vec) {
    assert ( key_vec.size() == treeDepth );
    unsigned depth = 0;    
    DataItem** valp = base_map.find(key_vec[depth]);
    if (!valp)
      return NULL;
    for (depth=1;depth<treeDepth;depth++) {
      NonTerminalDataItem* ntdi = (NonTerminalDataItem*)(*valp);
      valp = ntdi->base_map.find(key_vec[depth]);
      if (!valp)
	return NULL;
    }
    // get the data item
    TerminalDataItem *tdi = (TerminalDataItem*)(*valp);
    return &(tdi->val);
  }


  //////////////////////////////////
  // a type def for readability below.
  typedef typename hash_map_list <_Key, DataItem*>::iterator 
            base_map_iterator_t;


  ///////////////////////////////////////////
  // a dummy type that is used to simulate STL's
  // way of doing iterator end (i.e., it == obj.end()).
  typedef char iterator_end_object;

  ///////////////////////////////////////////
  // Return an iterator end object, any one will do, just
  // to trigger the end call of the iterator.
  iterator_end_object end() { return (iterator_end_object)0; }

  //////////////////////////////////////////////////////
  // An iterator class for this tree. It iterates over all leaf nodes
  // in the tree by simultaneously iterating through the different
  // hash tables at each tree level.
  class iterator {
    friend class hash_tree<_Key,_Data>;

  protected:
    sArray< base_map_iterator_t > it_array;
    bool endp;
  public:

    iterator& operator ++(int) {
      // test for end.
      unsigned depth = it_array.size()-1;
      const unsigned depthm1 = depth;

      it_array[depth]++;
      while (it_array[depth].end() && depth != 0) {
	depth--;
	it_array[depth]++;
      }
      if (depth == 0 && it_array[depth].end()) {
	endp = true;
	return *this; // we're at the real end.
      }
      // this next bit of code works because the tree is always
      // perfectly balanced (i.e., the length from the root to all
      // leaf nodes is always the same for all leaves).
      while (depth < depthm1) {
	NonTerminalDataItem* ntdi = 
	  (NonTerminalDataItem*)(*it_array[depth]);
	ntdi->base_map.begin(it_array[++depth]);
      }
      return *this;
    }

    _Data& operator*() { 
      TerminalDataItem* tdi = 
	(TerminalDataItem*)
	// TODO: optimize this to use ptr and pre-compute size()-1.
	(*(it_array[it_array.size()-1]));
      return tdi->val;
    }

    bool end() {
      // cheap test for end.
      return endp;
    }

    // use operator == only with a dummy end object.
    bool operator == (iterator_end_object eo) {
      return end();
    }
    bool operator != (iterator_end_object eo) {
      return !end();
    }

  };


  //////////////////////////////////////////////////////
  // Another iterator class for this tree. This one is
  // similar to the above, but also keeps track of 
  // the path from root to leaf node in the array
  // key_vec. Only those values in this vector that
  // change are modified, so that iterating through
  // the tree does not involve multiple repeated writes
  // to the array with the same value.
  class iterator_vector : public iterator {

    friend class hash_tree<_Key,_Data>;
  protected:
    // a length treeDepth vector that stores
    // the current values in the tree during
    // and iteration. Only values that change
    // are re-writen on each iteration.
    _Key* key_vec;
  public:

    iterator_vector& operator ++(int) {

      // test for end.
      unsigned depth = iterator::it_array.size()-1;
      const unsigned depthm1 = depth;

      iterator::it_array[depth]++;
      while (iterator::it_array[depth].end() && depth != 0) {
	depth--;
	iterator::it_array[depth]++;
      }
      if (depth == 0 && iterator::it_array[depth].end()) {
	iterator::endp = true;
	return *this; // we're at the real end.
      }
      while (depth < depthm1) {
	key_vec[depth] = iterator::it_array[depth].key();
	NonTerminalDataItem* ntdi = 
	  (NonTerminalDataItem*)(*iterator::it_array[depth]);
	ntdi->base_map.begin(iterator::it_array[++depth]);
      }
      key_vec[depth] = iterator::it_array[depth].key();
      return *this;
    }
  };


  //////////////////////////////////////////////////////
  // Still another iterator class for this tree. This one is similar
  // to the above, but also keeps track of the path from root to leaf
  // node in the array of pointers to _Keys, key_vecp. Only those
  // values in this vector that change are modified, so that iterating
  // through the tree does not involve multiple repeated writes to the
  // array with the same value.
  class iterator_vectorp : public iterator {

    friend class hash_tree<_Key,_Data>;
  protected:
    // a length treeDepth vector that stores
    // the current values in the tree during
    // and iteration. Only values that change
    // are re-writen on each iteration.
    _Key** key_vecp;
  public:

    iterator_vectorp& operator ++(int) {
      // test for end.
      unsigned depth = iterator::it_array.size()-1;
      const unsigned depthm1 = depth;

      iterator::it_array[depth]++;
      while (iterator::it_array[depth].end() && depth != 0) {
	depth--;
	iterator::it_array[depth]++;
      }
      if (depth == 0 && iterator::it_array[depth].end()) {
	iterator::endp = true;
	return *this; // we're at the real end.
      }
      while (depth < depthm1) {
	(*key_vecp[depth]) = iterator::it_array[depth].key();
	NonTerminalDataItem* ntdi = 
	  (NonTerminalDataItem*)(*iterator::it_array[depth]);
	ntdi->base_map.begin(iterator::it_array[++depth]);
      }
      (*key_vecp[depth]) = iterator::it_array[depth].key();
      return *this;
    }
  };


  ////////////////////////////////////////////////////////////
  // There are multiple ways to begin an iteration.  They all use the
  // begin(iterator) notion rather than the iterator = obj.begin(),
  // thereby avoiding temporaries, copies, and redundant object
  // copying.

  ///////////////////////////////////////////////////////////////
  // create an iterator through the entire tree.
  void begin(iterator& it) {
    it.it_array.resize(treeDepth);
    it.endp = false;
    unsigned depth = 0;
    base_map.begin(it.it_array[depth]);
    for (depth=1;depth<treeDepth;depth++) {
      NonTerminalDataItem* ntdi = 
	(NonTerminalDataItem*)(*it.it_array[depth-1]);
      ntdi->base_map.begin(it.it_array[depth]);
    }
  }

  ///////////////////////////////////////////////////////////////
  // create an iterator but start at the partial path
  // down the tree given by start_key_vec which has length start_depth.
  // If start_depth == 0, then this is a normal iterator. 
  void begin(const _Key*const start_key_vec,
	     const unsigned start_depth,
	     iterator& it) {
    assert(start_depth < treeDepth);
    if (start_depth == 0) {
      // same as a standard iterator
      return begin(it);
    }
    // first find the starting hash table to
    // start iterating at.
    unsigned depth = 0;
    DataItem** valp = base_map.find(start_key_vec[depth]);
    if (!valp) {
      it.endp = true;
      return;
    }
    for (depth = 1;depth <start_depth; depth++) {
      NonTerminalDataItem* ntdi = (NonTerminalDataItem*)(*valp);
      valp = ntdi->base_map.find(start_key_vec[depth]);
      if (!valp) {
	it.endp = true;
	return;
      }
    }
    // now create iterators the rest of the way down.
    it.it_array.resize(treeDepth-start_depth);
    it.endp = false;
    NonTerminalDataItem* ntdi = (NonTerminalDataItem*)(*valp);
    ntdi->base_map.begin(it.it_array[depth-start_depth]);
    for (depth++;depth<treeDepth;depth++) {
      NonTerminalDataItem* ntdi = 
	(NonTerminalDataItem*)(*it.it_array[depth-start_depth-1]);
      ntdi->base_map.begin(it.it_array[depth-start_depth]);
    }
  }


  ///////////////////////////////////////////////////////////////
  // create an iterator through the entire tree.  Keep track of the
  // current path from root to leaf in the array key_vec, the storage
  // of which is managed outside this object (i.e., key_vec is written
  // to and not ever deleted here).
  void begin(iterator_vector& it,_Key*key_vec) {
    it.it_array.resize(treeDepth);
    it.key_vec = key_vec;
    it.endp = false;
    unsigned depth = 0;
    base_map.begin(it.it_array[depth]);
    key_vec[depth] = it.it_array[depth].key();
    for (depth=1;depth<treeDepth;depth++) {
      NonTerminalDataItem* ntdi = 
	(NonTerminalDataItem*)(*it.it_array[depth-1]);
      ntdi->base_map.begin(it.it_array[depth]);
      key_vec[depth] = it.it_array[depth].key();
    }
  }



  ///////////////////////////////////////////////////////////////
  // create an iterator but start at the partial path down the tree
  // given by start_key_vec which has length start_depth.  If
  // start_depth == 0, then this is a normal iterator.  
  // key_vec is the same as and therefore described above.
  void begin(const _Key*const start_key_vec,
	     const unsigned start_depth,
	     iterator_vector& it,
	     _Key*key_vec) {
    assert(start_depth < treeDepth);
    if (start_depth == 0) {
      // same as a standard iterator
      return begin(it,key_vec);
    }
    // first find the starting hash table to
    // start iterating at.
    unsigned depth = 0;
    DataItem** valp = base_map.find(start_key_vec[depth]);
    if (!valp) {
      it.endp = true;
      return;
    }
    for (depth = 1;depth <start_depth; depth++) {
      NonTerminalDataItem* ntdi = (NonTerminalDataItem*)(*valp);
      valp = ntdi->base_map.find(start_key_vec[depth]);
      if (!valp) {
	it.endp = true;
	return;
      }
    }
    // now create iterators the rest of the way down.
    it.it_array.resize(treeDepth-start_depth);
    it.key_vec = key_vec;
    it.endp = false;
    NonTerminalDataItem* ntdi = (NonTerminalDataItem*)(*valp);
    ntdi->base_map.begin(it.it_array[depth-start_depth]);
    key_vec[depth-start_depth] = it.it_array[depth-start_depth].key();
    for (depth++;depth<treeDepth;depth++) {
      NonTerminalDataItem* ntdi = 
	(NonTerminalDataItem*)(*it.it_array[depth-start_depth-1]);
      ntdi->base_map.begin(it.it_array[depth-start_depth]);
      key_vec[depth-start_depth] = it.it_array[depth-start_depth].key();
    }
  }



  ///////////////////////////////////////////////////////////////
  // create an iterator through the entire tree.  Keep track of the
  // current path from root to leaf in the array of pointers key_vecp,
  // the storage of which is managed outside this object (i.e.,
  // key_vec is written to and not ever deleted here).
  void begin(iterator_vectorp& it,_Key**key_vecp) {
    it.it_array.resize(treeDepth);
    it.key_vecp = key_vecp;
    it.endp = false;
    unsigned depth = 0;
    base_map.begin(it.it_array[depth]);
    (*key_vecp[depth]) = it.it_array[depth].key();
    for (depth=1;depth<treeDepth;depth++) {
      NonTerminalDataItem* ntdi = 
	(NonTerminalDataItem*)(*it.it_array[depth-1]);
      ntdi->base_map.begin(it.it_array[depth]);
      (*key_vecp[depth]) = it.it_array[depth].key();
    }
  }


  ///////////////////////////////////////////////////////////////
  // create an iterator but start at the partial path down the tree
  // given by start_key_vec which has length start_depth.  If
  // start_depth == 0, then this is a normal iterator.  
  // key_vecp is the same as and therefore described above.
  void begin(const _Key*const start_key_vec,
	     const unsigned start_depth,
	     iterator_vectorp& it,
	     _Key**key_vecp) {
    assert(start_depth < treeDepth);
    if (start_depth == 0) {
      // same as a standard iterator
      return begin(it,key_vecp);
    }
    // first find the starting hash table to
    // start iterating at.
    unsigned depth = 0;
    DataItem** valp = base_map.find(start_key_vec[depth]);
    if (!valp) {
      it.endp = true;
      return;
    }
    for (depth = 1;depth <start_depth; depth++) {
      NonTerminalDataItem* ntdi = (NonTerminalDataItem*)(*valp);
      valp = ntdi->base_map.find(start_key_vec[depth]);
      if (!valp) {
	it.endp = true;
	return;
      }
    }
    // now create iterators the rest of the way down.
    it.it_array.resize(treeDepth-start_depth);
    it.key_vecp = key_vecp;
    it.endp = false;
    NonTerminalDataItem* ntdi = (NonTerminalDataItem*)(*valp);
    ntdi->base_map.begin(it.it_array[depth-start_depth]);
    (*key_vecp[depth-start_depth]) = it.it_array[depth-start_depth].key();
    for (depth++;depth<treeDepth;depth++) {
      NonTerminalDataItem* ntdi = 
	(NonTerminalDataItem*)(*it.it_array[depth-start_depth-1]);
      ntdi->base_map.begin(it.it_array[depth-start_depth]);
      (*key_vecp[depth-start_depth]) = it.it_array[depth-start_depth].key();
    }
  }



};


#endif

