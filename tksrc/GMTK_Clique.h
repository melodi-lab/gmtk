
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
#include "GMTK_Hash.h"

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

    // The instantiation number in the order that the instantiations were
    // enumerated.
    int inum; 

    // In a dynamic network, the same set of values will occur over and
    // over again in different cliques. To avoid storing them over and over
    // again, keep a global pool. 
    /**** 
     * Note: With logspace in place, the overall memory requirements can 
     * be significantly reduced by eliminating this hash table. Instead,
     * the "values" field can be changed from a a pointer into the hash table 
     * to a straight vector. The variable assignments associated with the
     * instantiation can be stored directly in values. Since there are only
     * a small number of instantiations, this is OK. This will help 
     * especially in pruning. The drawback is that it will be horribly 
     * memory-inefficient in the non-logspace case. Doing this will be 
     * ultra-memory efficient, but force the use of the logspace recurions.
    */
    static HashTable< vector<RandomVariable::DiscreteVariableType> >
         global_val_set;
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

    int *instantiation;
    // This stores all the possible instantiations of a clique.
    // Pruning occurs by removing low probability instantiations.
    // The storage is indirect: the integer refers to the index of the
    // instantiation in the global instantiation pool

    int numInstantiations;  // how many are there?
    int instantiation_capacity;  // how big is the instantiation vector?

    void addInstantiation(int inum);  // adds an instantiation

    void clearInstantiationList(); // deletes it

    vector<int> *inum2gnum;
    // A mapping from instantiation number X (before pruning) to 
    // the index in the gip where it's stored
    // Using a self-managed int * would be OK too.
    // It's a vector * here because deleting one gets the memory back,
    // whereas clear()ing one doesn't

    static map<vector<RandomVariable::DiscreteVariableType>, int>
           instantiationAddress;
    // A separator clique sums over multiple values from its non-separator
    // parent. instantiationAddress keeps track of the index of the 
    // unique CliqueValue that all the parent values with the same 
    // underlying variable values sum into.
    // Static because only one copy is needed at a time.

    static vector<CliqueValue> gip;  // the global instantiation pool

    static vector<unsigned> freelist;     
    // keeps track of which CliqueVals are available; the next available may
    // not simply be the next entry in the gip, because pruning frees up
    // CliqueValues everywhere.

    static int nextfree;  // next free thing on the freelist

    static unsigned newCliqueValue();  // return a clique value (index) to use

    static void recycleCliqueValue(unsigned idx); // frees this one for reuse

    void recycleAllCliqueValues() {for (int i=0;i<numInstantiations; i++)
                                       recycleCliqueValue(instantiation[i]);}

    logpr probGivenParents();
    // This computes the conditional probability of the clique.
    // It calls findConditionalProbabilityNodes(), and then for each
    // member of conditionalProbabilityNode it invokes probGivenParents()
    
    void enumerateValues(int new_member_num, int pred_val, vector<int> *in2gn,
    bool viterbi=false);
    // This is a recursive function that enumerates all the possible
    // instantiations of a clique.  For consistency with the value of
    // the parent clique, only the new members are instantiated.

    void prune(logpr beam);
    // Removes all instantiations whose forwards probability is less than
    // beam*max.

    void reveal();
    // shows who the members and newMembers are

    Clique() {instantiation=NULL; numInstantiations=0;
    instantiation_capacity=0; inum2gnum=NULL;}

    ~Clique() {recycleAllCliqueValues(); assert(inum2gnum==NULL);}

     // for incremental use in the logspace algorithm
     void reclaimMemory();
};

#endif
