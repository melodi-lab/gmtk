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



#include "general.h"
VCID("$Header$");

#ifdef MAIN
#define COLLECT_COLLISION_STATISTICS
#endif

#include "GMTK_Hash.h"

// the possible sizes for the has table, a list of prime
// numbers each entry roughly doubling in size of the 
// previous entry. Comment the ones at the beginning to get
// the desired starting size of the hash table.
const unsigned HashTable_PrimesArray[] = {
3,
5,
11,
23,
53,
107,
211,
421,
743,
1531,
3709,
6317,
12437,
24547,
47779,
75403,
134053,
262139,
524287,
1048573,
2097143,
4194301,
8388593,
16777213,
33554393,
67108859,
134217689,
268435399,
536870909,
1073741789,
2147483647
};

/////////////////////////////////////////////////////////////
// the size of the above table, determined automatically.
const unsigned HashTable_SizePrimesArray
   = sizeof(HashTable_PrimesArray)/sizeof(unsigned);

#ifdef MAIN

///////////////////////////////////////////
// main driver debugger for hash table.

#include <string>

int main(int argc,char*argv[])
{

  // first insert some numbers.
  int count= 100000;
  int vsize = 3;
  int maxCard = 10;
  
  if (argc > 1)
    count = atoi(argv[1]);
  if (argc > 2)
    vsize = atoi(argv[2]);
  if (argc > 3)
    maxCard = atoi(argv[3]);

  printf("prime table size = %d\n",HashTable_SizePrimesArray);
  printf("Using %d entries, each of length %d\n",count,vsize);

  HashTable< vector<int> > ht;

  printf("inserting %d random entries of size %d\n",count,vsize);
  vector<int> vi;
  vi.resize(vsize);
  for (int i=0;i<count;i++) {
    for (int j=0;j<vsize;j++) {
      vi[j] = rand() % maxCard;
    }
    ht.insert(vi);
  }
  
  printf("max %d cols, %d total collisions with %d inserts, avg = %e, num entries = %d, table size=%d\n",
	 ht.maxCollisions,ht.numCollisions,ht.numInserts,
	 (double)ht.numCollisions/(double)ht.numInserts,
	 ht. totalNumberEntries(),
	 HashTable_PrimesArray[ht.primesArrayIndex]);

  // now go through each entry in the table to test if it is there.
  printf("ensuring inserted entries are all contained\n");
  for (unsigned pos=0;pos<ht.table.size();pos++) {
    if (ht.table[pos] != NULL) {
      // found an entry, make sure it is there.
      if (!ht.isContained(*ht.table[pos]))
	printf("entry should be in hash table but isn't");
    }
  }

  // next go through each entry in table to make sure it
  // returns same pointer to itself.
  printf("ensuring pointers returned are proper\n");
  for (unsigned pos=0;pos<ht.table.size();pos++) {
    if (ht.table[pos] != NULL) {
      // found an entry, make sure it is there.
      if (ht.table[pos] != ht.insert(*ht.table[pos]))
	printf("entry should have same pointers");
    }
  }


  // finally, try some new random entries and count
  // the number of hits and mises.
  printf("checking if any random entries are contained\n");
  int nhits = 0;
  for (int i=0;i<count;i++) {
    for (int j=0;j<vsize;j++) {
      vi[j] = rand() % maxCard;
    }
    if (ht.isContained(vi))
      nhits++;
  }
  printf("Found %d out of %d with random vectors\n",nhits,count);

  // next, try hashing some C++ strings.
  count = 0;
  printf("string hash table: inserting %d random strings of various sizes\n",count);
  HashTable< string > sht;
  string str;
  for (int i=0;i<count;i++) {
    unsigned size = (rand() % 30) + 1;
    str.resize(size);
    for (unsigned j=0;j<size;j++) {
      str[j] = ' ' + (rand() % ('~' - ' '));
    }
    sht.insert(str);
  }



}

#endif
