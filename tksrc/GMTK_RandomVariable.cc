

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


#include <set>
#include "GMTK_RandomVariable.h"

/*-
 *-----------------------------------------------------------------------
 *  setParents
 *     sets up the switching and conditional parents.
 *
 * Results:
 *
 * Side Effects:
 *     sets the switchingParents array
 *     sets the conditionalParentsList array
 *     sets the allPossibleParents Array
 *     adds to the allPossibleChildren array of the variable's parents
 *
 *-----------------------------------------------------------------------
 */

void RandomVariable::setParents(vector<RandomVariable *> &sparents,
         vector<vector<RandomVariable *> > &cpl)
{
    switchingParents = sparents;
    conditionalParentsList = cpl;

    // add this to the child list of each of its parents
    // note that a parent may occur twice, e.g. as a switching and conditional
    // parent, or on twqo different conditional parent lists
    // avoid adding this to the child list twice
    set<RandomVariable *> parents;  // the set of parents this is a child of
    for (unsigned i=0; i<sparents.size(); i++)
    {
         parents.insert(sparents[i]);
         sparents[i]->allPossibleChildren.push_back(this);
    }
    
    for (unsigned i=0; i<cpl.size(); i++)
        for (unsigned j=0; j<cpl[i].size(); j++)
            if (parents.find(cpl[i][j]) == parents.end())
            {
                cpl[i][j]->allPossibleChildren.push_back(this);
                parents.insert(cpl[i][j]);
            }

    // set up the all possible parents array.
    allPossibleParents.resize(parents.size());
    set<RandomVariable *>::iterator si;
    int p=0;
    for (si = parents.begin(); si != parents.end(); si++)
        allPossibleParents[p++] = *si;
}


/*-
 *-----------------------------------------------------------------------
 * clone
 *      copies the data structures necessary for unrolling
 *      all random variables have these structures. Note that
 *      this routine clones *ONLY* those attributes that
 *      are abstract, in that they are associated with 
 *      all derived classes of RandomVariable. 
 *
 * Preconditions:
 *      dtMapper and all parents lists must be set appropriately.
 *
 * Postconditions:
 *      specified data members are copied
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      A cloned RV
 *
 *-----------------------------------------------------------------------
 */

RandomVariable* RandomVariable::clone()
{
  RandomVariable* rv = create();
  rv->label = label;
  rv->hidden = hidden;
  rv->discrete = discrete;
  rv->timeIndex = timeIndex;
  rv->switchingParents = switchingParents;
  rv->allPossibleParents = allPossibleParents;
  rv->allPossibleChildren = allPossibleChildren;
  rv->conditionalParentsList = conditionalParentsList;
  // leave curConditionalParents empty
  rv->dtMapper = dtMapper; 
  // leave cachedIntFromSwitchingState uninitialized since
  // it will be changed individually for each rv.

  return rv;
}



/*-
 *-----------------------------------------------------------------------
 * identicalStructureWith.
 *      Returns true if this rv has identical structure with that of other.
 *      "identical structure" means that the r.v. have the same
 *      number, type, and cardinality parents. If this returns
 *      true, then it will be valid to tie parameters between
 *      these two random variables.
 *
 * 
 * Preconditions:
 *      Both rvs must have parents filled in.
 *
 * Postconditions:
 *      If function returns true, then the variables have
 *      identical structure, otherwise not.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      indicating boolean.
 *
 *-----------------------------------------------------------------------
 */
bool
RandomVariable::identicalStructureWith(const RandomVariable& other)
{

  if (discrete != other.discrete)
    return false;

  // check that switching parents are compatible
  if (switchingParents.size() != other.switchingParents.size())
    return false;
  for (unsigned i=0;i<switchingParents.size();i++) {
    if (switchingParents[i]->discrete !=
	other.switchingParents[i]->discrete)
      return false;
    if (switchingParents[i]->discrete) {
      if (switchingParents[i]->cardinality !=
	  other.switchingParents[i]->cardinality)
	return false;
    }
  }

  // switching parents check out ok, now check all sets
  // of conditional parents
  if (conditionalParentsList.size() !=
      other.conditionalParentsList.size())
    return false;

  for (unsigned i=0;i<conditionalParentsList.size();i++) {
    if (conditionalParentsList[i].size() !=
	other.conditionalParentsList[i].size())
      return false;
    for (unsigned j=0;j<conditionalParentsList[i].size();j++) {
      if (conditionalParentsList[i][j]->discrete !=
	  other.conditionalParentsList[i][j]->discrete)
	return false;
      if (conditionalParentsList[i][j]->discrete) {
	if (conditionalParentsList[i][j]->cardinality !=
	    other.conditionalParentsList[i][j]->cardinality)
	  return false;
      }
    }
  }
  return true;
}




/*-
 *-----------------------------------------------------------------------
 * reveal()
 *      i.e., "print" out the contents of the rv.
 * 
 * Preconditions:
 *      rv must be filled in.
 *
 * Postconditions:
 *      Same as before.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void RandomVariable::reveal(bool show_vals)
{
    cout << label << "-" << timeIndex << " : ";
    if (discrete) cout << "discrete (" << cardinality << ") ";
    if (hidden) cout << "hidden "; else cout << "observed ";
    if (!hidden || show_vals) 
        if (discrete)
            cout << "val (" << val << ")";
        else 
            cout << "continuous";

/*
    cout << " possible parents: ";
    for (unsigned i=0; i<allPossibleParents.size(); i++)
        cout << allPossibleParents[i]->label << "-" 
             << allPossibleParents[i]->timeIndex << " ";
*/
    cout << endl;
}

