/*
 * "Copyright 2001, International Business Machines Corporation and University
 * of Washington. All Rights Reserved
 *
 *    Written by Geoffrey Zweig and Jeff Bilmes
 *
 * NO WARRANTY
 * THE PROGRAM IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED INCLUDING, WITHOUT
 * LIMITATION, ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT,
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Each Recipient is
 * solely responsible for determining the appropriateness of using the Program
 * and assumes all risks associated with such use, including but not limited
 * to the risks and costs of program errors, compliance with applicable laws,
 * damage to or loss of data, programs or equipment, and unavailability or
 * interruption of operations.

 * DISCLAIMER OF LIABILITY
 * THE UNIVERSITY OF WASHINGTON, INTERNATIONAL BUSINESS MACHINES CORPORATION,
 * GEOFFREY ZWEIG AND JEFF BILMES SHALL NOT HAVE ANY LIABILITY FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING WITHOUT LIMITATION LOST PROFITS), HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE  OF
 * THE PROGRAM, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES."
*/

#ifndef GMTK_HASH_H
#define GMTK_HASH_H

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>

#include "general.h"
#include "error.h"

extern const unsigned HashTable_SizePrimesArray;
extern const unsigned HashTable_PrimesArray[];

/*
 * Ten random 10 digit primes
 *
 * unsigned long tenDigitPrimes[] = {
 * 5915587277u, 
 * 1500450271u, 
 * 3267000013u, 
 * 5754853343u, 
 * 4093082899u, 
 * 9576890767u, 
 * 3628273133u, 
 * 2860486313u, 
 * 5463458053u, 
 * 3367900313u };
 * 
 */

///////////////////////////////////////////////////////////
// A generic hash table that hashes T objects.  
// Note: T *MUST* be a STL vector< something >, or a C++ 'string'
// or any other object that acts like an array of int types.
//
template <class T>
class HashTable
{
  // uncomment to get test main() stuff working.
public:

  ////////////////////////////////
  // total number of entries in the hash table
  unsigned _totalNumberEntries;

  ///////////////////////////////////////////////////////////////////
  // the possible sizes for the has table. These need to be
  // prime numbers for correctness.

  // the index into the above prime array of the current size.
  unsigned primesArrayIndex; 

  //////////////////////////////////////////////////////////
  // the actual hash table, an array of pointers to T's
  vector < T* > table;

#ifdef COLLECT_COLLISION_STATISTICS
  unsigned maxCollisions;
  unsigned numCollisions;
  unsigned numInserts;
#endif

public:

  // constructors
  HashTable(unsigned approximateStartingSize = 200000) {
    _totalNumberEntries=0;
    for (primesArrayIndex=0;
	 primesArrayIndex<HashTable_SizePrimesArray;
	 primesArrayIndex++) {
      if (HashTable_PrimesArray[primesArrayIndex] >= approximateStartingSize)
	break;
    }
    if (primesArrayIndex == HashTable_SizePrimesArray)
      error("ERROR: Can't create hash table of approximate size %u\n",
	    approximateStartingSize);
#ifdef COLLECT_COLLISION_STATISTICS
    maxCollisions = 0;
    numCollisions = 0;
    numInserts = 0;
#endif
  }
  ~HashTable() {
    clear();
  } 

  //////////////////////////////////////////////////////////////////
  // since this is a double hash table, we define two address functions,
  // addr which gives the starting position, and incr which
  // gives the increment when we have a collision.
  unsigned addr(T & vec, const unsigned size)
  {
    assert ( size > 1 );
    unsigned long a = vec.size();
    int i = (int)(vec.size()-1); do {
      // a =65599*a + vec[i] + 1;
      // a = 402223*a + vec[i] + 1;
      // a = 611953*a + vec[i] + 1;
      a = 3367900313ul*a + vec[i] + 1;
    } while (--i >= 0);
    return a % size;
  }
  // the increment of vec for table of size 'size' when
  // we have a collision.
  // Note: the result of the function must satisfy
  //    1) it must be greater than zero
  //    2) it must be strictly less than size
  //         
  unsigned incr(T &vec, const unsigned size)
  {
    assert ( size > 1 );
    unsigned long a=0;

    int i = (int)(vec.size()-1); do {
      // a =65599*a + vec[i] + 1;
      // a = 402223*a + vec[i] + 1;
      // a = 611953*a + vec[i] + 1;
      // a = 1500450271ul*a + vec[i] + 1;
      a = 3267000013ul*a + vec[i] + 1;
    } while (--i >= 0);

#if 0
    int i=(int)(vec.size()-1); do {
      a = (a << 4) + vec[i] + 1;
      if (a > 0x0fffffff) {
	a ^= (a >> 24) & 0xf0;
	a &= 0x0fffffff;
      }
    } while (--i >= 0);
#endif

    return (
       (a % (size-1)) // this gives [0 : (size-2)]
       +
       1              // this gives [1 : (size-1)] 
       );
  }

  //////////////////////////////////////////////////////////////////
  // return the entry of vec according to table tbl.
  unsigned entryOf(T &vec, vector < T* >& tbl) {

    unsigned a = addr(vec, tbl.size());

#ifdef COLLECT_COLLISION_STATISTICS
    unsigned collisions=0;
#endif

    if ( tbl[a] == NULL || *tbl[a] == vec) {
      return a;
    } 
    // don't compute inc unless we need to.
    const unsigned inc = incr(vec, tbl.size());
    do {
#ifdef COLLECT_COLLISION_STATISTICS
      collisions++;
#endif
      a = (a+inc) % tbl.size(); 
    } while ( tbl[a] != NULL && *tbl[a] != vec);

#ifdef COLLECT_COLLISION_STATISTICS
    if (collisions > maxCollisions)
      maxCollisions = collisions;
    numCollisions += collisions;
#endif

    return a;
  }

  ///////////////////////////////////////////////////////////////
  // insert: insert a new element in the hash table if it is not
  // there, and return a pointer to that element if it is there
  // or not.
  T* insert(T &vec) {
    // can only insert non-zero size entries.
    assert (vec.size() > 0); 

#ifdef COLLECT_COLLISION_STATISTICS
    numInserts++;
#endif

    // make sure the table has entries
    if (table.size() == 0) 
      resize(HashTable_PrimesArray[primesArrayIndex]);
    // compute the address
    const unsigned a = entryOf(vec,table);

    T *nv = table[a];
    if (table[a] == NULL ) {
      // need to insert since not there.
      table[a] = nv = new T(vec);

      // time to resize if getting too big.
      if (++_totalNumberEntries >= table.size()/2) {
	  if (primesArrayIndex == (HashTable_SizePrimesArray-1)) 
	    error("ERROR: Hash table error, table size exceeds max size of %lu",
		  HashTable_PrimesArray[primesArrayIndex]);
	  resize(HashTable_PrimesArray[++primesArrayIndex]);
      }
    }
    return nv;
  }

  //////////////////////////////////////////////////////
  // returns true if vec is contained in table
  bool isContained(T& vec) {
    const unsigned a = entryOf(vec,table);    
    return (table[a] != NULL);
  }

  ////////////////////////////////////////////////////////////
  // resize: resizes the table to be new_size.
  // new_size *MUST* be a prime number.
  void resize(int new_size) {
    // make a new table (nt),  and re-hash everyone in the
    // old table into the new table.

    // the next table, used for table resizing.
    vector < T* > nt;

    nt.resize(new_size);
    for (int i=0; i<new_size; i++) 
      nt[i] = NULL;
    for (unsigned i=0; i<table.size(); i++) {
        if (table[i]) {
	  // find the address of table[i] in the new table
	  unsigned a = entryOf(*table[i],nt);
	  nt[a] = table[i];
        }
    }
    // finally, copy the new table onto the old table
    table=nt;
  }

  /////////////////////////////////////////////////////////
  // clear out the table entirely, including deleting
  // all memory pointed to by the T* pointers.
  void clear() {
    for (unsigned i=0; i<table.size(); i++)
      if (table[i]) 
	delete table[i];
    table.clear();
    _totalNumberEntries=0;
    primesArrayIndex=0;
  }

  ////////////////////////////////////////////////////////////////
  // return the total number of entries that have
  // been inserted into the hash table (this is different
  // than the total allocated size of the table.
  unsigned totalNumberEntries() { return _totalNumberEntries; }

  
};

#endif
