
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

#ifndef GMTK_CLIQUE_H
#define GMTK_CLIQUE_H 

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include "GMTK_RandomVariable.h"

#include "general.h"


#define MAX_VALUE_HASH_TABLE_SIZE 200000

struct ValueHashTable
{
    vector<vector<RandomVariable::DiscreteVariableType> *> table, nt;
    int table_size, count;
    int addr(vector<RandomVariable::DiscreteVariableType> &vec)
         // return the hash value of the vector
        {unsigned long a = vec.size();
         for (unsigned i=0; i<vec.size(); i++) a = 65599*a + vec[i];
         a = a % table_size;
         return a;}

    vector<RandomVariable::DiscreteVariableType> *insert(
        vector<RandomVariable::DiscreteVariableType> &vec)
        {
	  // if (count >= MAX_VALUE_HASH_TABLE_SIZE/2)
	  //   warning("WARNING: hash table getting large, causing multiple resizes");
	  // if (count >= MAX_VALUE_HASH_TABLE_SIZE)
	  // error("ERROR: hash table too large, needing > %d entries",
	  // MAX_VALUE_HASH_TABLE_SIZE);

	  if (++count>=table_size/2)
	    resize(max(2*table_size,MAX_VALUE_HASH_TABLE_SIZE)); 
	  int a = addr(vec);
	  while (table[a] && *table[a]!=vec) a=(a+1)%table_size;
	  if (!table[a])
	    table[a] = new vector<RandomVariable::DiscreteVariableType>(vec);
	  return table[a];
	}

    void clear() 
        {for (unsigned i=0; i<table.size(); i++) if (table[i]) delete table[i];
         table.clear(); table_size=count=0;}
    ~ValueHashTable() {clear();} 
    ValueHashTable() {count=table_size=0;}
    void resize(int size) 
        {table_size = size;  // used for addressing; not same as table.size()
         nt.resize(table_size);
         for (int i=0; i<table_size; i++) nt[i] = NULL;
         for (unsigned i=0; i<table.size(); i++)
             if (table[i]) {int a=addr(*table[i]); 
              while (nt[a] && *nt[a]!=*table[i]) a=(a+1)%table_size;
              nt[a] = table[i];}
         table=nt;}
};

struct CliqueValue
{
    // A clique has a bunch of possible values, depending on the possible
    // combinations of its constituent variables. With just discrete hidden
    // variables, these can be enumerated out. The probabilities
    // associated with clique values will be stored in these structures.

    logpr lambda, pi;

    // For the clique values associated with non-separators, it is convenient 
    // to keep track of the clique values that are consistent in the surrounding
    // separators. Note that there is just one in each of the surrounding
    // separators. 
    int pred, succ;

    // The underlying variable values corresponding to a particular clique
    // instantiation
    const vector<RandomVariable::DiscreteVariableType> *values;

    // In a dynamic network, the same set of values will occur over and
    // over again in different cliques. To avoid storing them over and over
    // again, keep a global pool. 
    static ValueHashTable global_val_set;
};

struct Clique
{
    vector<RandomVariable *> member;
    // vector of pointers to the constituent variables in the underlying GM.

    vector<RandomVariable *> newMember;
    // vector of pointers to the variables that are present in the clique,
    // but not in its parent.
    // The entries must be organized so that a variable's parents occur
    // before the variable.

    vector<RandomVariable *> discreteMember;
    // A list of the discrete members of the clique.

    vector<RandomVariable *> conditionalProbabilityNode;
    // Each clique has a set of nodes assigned to it, that contribute to
    // its conditional probability. This array stores them.
    // it is set on demand by findConditionalProbabilityNodes()

    vector<RandomVariable::DiscreteVariableType> clampedValues;
    // What are the values of the discrete variables (which should be clamped)
    // set in the course of inference

    void cacheClampedValues();
    // Reads the values of the discrete members and stores them in clampedValues

    void findConditionalProbabilityNodes() 
    {; /* current implementation ignores possible switching efficiencies */}
    // With switching parents, the set of variables assigned to a clique 
    // depends on the value of the clique. This function is called after all
    // the variables have been clamped, and stores the appropriate 
    // RandomVariable pointers in conditionalProbabilityNode.

    bool separator;
    // Is the clique a separator?

    vector<int> instantiation;
    // This stores all the possible instantiations of a clique.
    // Pruning occurs by removing low probability instantiations.
    // The storage is indirect: the integer refers to the index of the
    // instantiation in the global instantiation pool

    map<vector<RandomVariable::DiscreteVariableType>, int>instantiationAddress;
    // A separator clique sums over multiple values from its non-separator
    // parent. instantiationAddress keeps track of the index of the 
    // unique CliqueValue that all the parent values with the same 
    // underlying variable values sum into.

    static vector<CliqueValue> gip;  // the global instantiation pool

    static vector<unsigned> freelist;     
    // keeps track of which CliqueVals are available; the next available may
    // not simply be the next entry in the gip, because pruning frees up
    // CliqueValues everywhere.

    static int nextfree;  // next free thing on the freelist

    static unsigned newCliqueValue();  // return a clique value (index) to use

    static void recycleCliqueValue(unsigned idx); // frees this one for reuse

    void recycleAllCliqueValues() {for (unsigned i=0;i<instantiation.size();i++)
                                       recycleCliqueValue(instantiation[i]);}

    logpr probGivenParents();
    // This computes the conditional probability of the clique.
    // It calls findConditionalProbabilityNodes(), and then for each
    // member of conditionalProbabilityNode it invokes probGivenParents()
    
    void enumerateValues(int new_member_num, int pred_val, bool viterbi=false);
    // This is a recursive function that enumerates all the possible
    // instantiations of a clique.  For consistency with the value of
    // the parent clique, only the new members are instantiated.

    void prune(logpr beam);
    // Removes all instantiations whose forwards probability is less than
    // beam*max.

    void reveal();
    // shows who the members and newMembers are

    ~Clique() {recycleAllCliqueValues();}
};

#endif
