
/* 
 * GMTK_CliqueChain.cc
 * Clique chain functions.
 *
 * "Copyright 2001, International Business Machines Corporation and University
 * of Washington. All Rights Reserved

 *    Written by Geoffrey Zweig and Jeff Bilmes

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

#include "GMTK_CliqueChain.h"

/* 
 * Initailize forward pass: free up any memory, and enumerate the values for
 * the root.
*/

map<int,int> copies;

void CliqueChain::initializeForwardPass()
{
// cout << "initializing forward pass\n";

    // memory processing should ensure that nothing is used at this point.
    assert(Clique::nextfree == int(Clique::gip.size())-1);

    // take care of the root
    assert(preorder(0)->inum2gnum == NULL);
    preorder(0)->inum2gnum = new vector<int>;
    preorder(0)->enumerateValues(0, -1, NULL, doingViterbi);
    preorder(0)->prune(beam);
    copies[0]++;
    assert(preorder(0)->instantiationAddress.size() == 0);
}

/* 
 * finishFrowardPass: compute 
 * data/viterbi probabilities. Return true if non-zero probability.
*/
bool CliqueChain::finishForwardPass()
{
// cout << "Finishing forward pass\n";
    // check if nothing had any probability, or pruning was too drastic
    if (postorder(0)->numInstantiations == 0) 
    {
        dataProb = viterbiProb = 0.0;
        return false; 
    }

    // look at the probabilities
    logpr sum=0.0, max=0.0;
    for (int j=0; j<postorder(0)->numInstantiations; j++)
    {
        int inst_idx = postorder(0)->instantiation[j];
        CliqueValue &cv = Clique::gip[inst_idx];
        sum += cv.pi;
        if (cv.pi > max) max = cv.pi;
    }
    dataProb=viterbiProb=0.0;
    if (doingViterbi) 
        viterbiProb=max;
    else 
        dataProb=sum;

    return true;
}

/* 
 * advancePi: advance and store pi values for starting from first_clique
 * and ending at last_clique. first_clique must be a separator and
 * last_clique must be a non-separator. The pis of the predecessor
 * of the first_clique must be known on entry. On exit, the pis of
 * the last clique are known.
*/

void CliqueChain::advancePis(unsigned first_clique, unsigned last_clique, 
bool erase)
{
// cout << "Advancing pis from " << first_clique << " to " << last_clique << "\n";
    assert(preorder(first_clique)->separator);
    assert(!preorder(last_clique)->separator);
    for (int i=int(first_clique)-1; i<int(last_clique); i++)
    {
        // This check against the number of copies of a clique is 
        // necessary because the rightmost clique in a block (which is
        // a non-separator) is saved when the pis are advanced. advancePis
        // may be called more than once with the same rightmost endpoint.
        // (Actually, just twice, once when the block is created, and then
        // once when the rightmost sub-block is processed as a base-case.)
        // It is not possible to delete the existing copy because then the
        // copy that is created will not have a valid succ entry - that is
        // filled in by calling advancePis on the following clique.
        if (copies[i+1]) 
        {
            assert(!preorder(i+1)->separator);
            // reclaim the memory from the clique that just pushed
            if (i-1>=int(first_clique) && erase)
            {
// cout << "reclaiming " << (i-1) << endl;
                preorder(i-1)->reclaimMemory();
                copies[i-1]--;
            }
            continue;
        }
        // propagate the surviving entries to preorder(i+1)
        copies[i+1]++;
        assert(preorder(i+1)->inum2gnum==NULL);
        preorder(i+1)->inum2gnum = new vector<int>;
        assert(preorder(i)->inum2gnum != NULL);
        if (preorder(i+1)->separator)
            assert(Clique::instantiationAddress.size() == 0);
        for (int j=0; j<preorder(i)->numInstantiations; j++)
        {
            if (!preorder(i)->separator)
            {
                // clamp the values of the variables in the clique
                CliqueValue &cv = Clique::gip[preorder(i)->instantiation[j]];
#ifdef NO_HASH_TABLE
                assert(cv.values.size()==preorder(i)->discreteMember.size());
                for (unsigned k=0; k<preorder(i)->discreteMember.size(); k++)
                    preorder(i)->discreteMember[k]->val = cv.values[k];
#else
                assert(cv.values->size()==preorder(i)->discreteMember.size());
                for (unsigned k=0; k<preorder(i)->discreteMember.size(); k++)
                    preorder(i)->discreteMember[k]->val = (*cv.values)[k];
#endif
            }
            else  
            {
                // clamp the values of the previous non-sep clique
                // there are more than in the separator, but this way we don't 
                // have to store values for separators
                vector<int> *in2gn = preorder(i-1)->inum2gnum;
                CliqueValue &cv =
                  Clique::gip[
                    (*in2gn)[Clique::gip[preorder(i)->instantiation[j]].pred]];
                assert(i>0);
#ifdef NO_HASH_TABLE
                assert(cv.values.size()==preorder(i-1)->discreteMember.size());
                for (unsigned k=0; k<preorder(i-1)->discreteMember.size(); k++)
                    preorder(i-1)->discreteMember[k]->val = cv.values[k];
#else
                assert(cv.values->size()==preorder(i-1)->discreteMember.size());
                for (unsigned k=0; k<preorder(i-1)->discreteMember.size(); k++)
                    preorder(i-1)->discreteMember[k]->val = (*cv.values)[k];
#endif
            }

            // compute the clique values of the child clique that are 
            // consistent with this parent instantiation.
            preorder(i+1)->enumerateValues(0, 
               Clique::gip[preorder(i)->instantiation[j]].inum, 
               preorder(i)->inum2gnum, doingViterbi); 
        }

        // we are done enumerating the values for preorder(i+1), and can
        // free the memory used in its instantiationAddress
        if (preorder(i+1)->separator)
            Clique::instantiationAddress.clear();
        else
            assert(Clique::instantiationAddress.size() == 0);

        // prune the low probability entries in preorder(i+1)
        // do not prune separators
        // this ensures that each entry on the instantaition list of a
        // non-separator clique will have a successor, and that li->succ->lambda
        // can  be dereferenced on the backward pass without an extra check.
        if (!preorder(i+1)->separator)
            preorder(i+1)->prune(beam);

        // reclaim the memory from the clique that just pushed
        if (i-1>=int(first_clique) && erase)
        {
// cout << "reclaiming " << (i-1) << endl;
            preorder(i-1)->reclaimMemory();
            copies[i-1]--;
        }
    }
    if (erase)
    {
// cout << "reclaiming " << (last_clique-1) << endl;
        preorder(last_clique-1)->reclaimMemory();
        copies[last_clique-1]--;
    }
}

/* 
 * initializeBackwardPass: compute the lambdas (all 1) for the final
 * clique, and the lambdas for its predecessor non-separator.
*/
void CliqueChain::initializeBackwardPass()
{
// cout << "initializing backward pass\n";
    if (doingViterbi)
    {
        // first find the likeliest instantiation of the last clique
        logpr maxprob = 0.0;
        best = -1;
        Clique *cl = postorder(0);
        for (int j=0; j<cl->numInstantiations; j++)
        {
            CliqueValue &cv = Clique::gip[cl->instantiation[j]];
            if (cv.pi >= maxprob)
            {
                maxprob = cv.pi;
                best = cl->instantiation[j];
            }
        }
        assert(maxprob==viterbiProb); // already computed on forward pass

        // clamp the values of the root variable
        for (unsigned j=0; j<cl->discreteMember.size(); j++)
#ifdef NO_HASH_TABLE
            cl->discreteMember[j]->val = Clique::gip[best].values[j];
#else
            cl->discreteMember[j]->val = (*Clique::gip[best].values)[j];
#endif

 
        if (preorderSize()>1)  // then there must be at least 3
        {
            // set best to the best instantiation of the predecessor of root 
            best = (*(postorder(1)->inum2gnum))[Clique::gip[best].pred];

            // set best to the best instantiation of the second to last 
            // non-separator clique
            best = (*(postorder(2)->inum2gnum))[Clique::gip[best].pred];
        }
    }
    else
    {
        // first do the last clique, where the lambdas are all 1.
        Clique *cl = postorder(0);
        for (int j=0; j<cl->numInstantiations; j++)
        {
            CliqueValue &cv = Clique::gip[cl->instantiation[j]];
            logpr t = cv.lambda;  // retrieve the cached probGivenParents
            cv.lambda = 1.0;
            if (postorderSize() > 1)  // watch out for degenerate case
                Clique::gip[(*(postorder(1)->inum2gnum))[cv.pred]].lambda += t;
        }
        if (doingEM)
            incrementEMStatistics(cl);
    }

    // all done with the last clique
    postorder(0)->reclaimMemory();
    copies[postorderSize()-1]--;
}

/* 
 * recedeLambdas: compute the lambdas for the cliques from right to left.
 * right must be a non-separator and left must be a separator. The
 * lambdas of the successor of the right clique must be known on
 * entry. On exit, the lambdas of the left clique are known.
 * One exit, all the instantiations from left+1 (the first separator)
 * to right+1 (the separator immediately following the chunk) are deleted.
*/
void CliqueChain::recedeLambdas(unsigned right, unsigned left)
{
// cout << "Receding lambdas from " << right << " to " << left << "\n";
    // compute the lambdas:
    // Non-separator clique values simply get the lambda of the separator
    // clique value that is descended from them.
    // Separator clique values get the sum of the lambdas of the non-separator
    // clique values descended from them, multiplied by the probGivenParents
    // of those descended values. Note that this was cleverly cached on the
    // forward pass.
    // separator lambdas are zero already

    assert(!preorder(right)->separator);
    assert(preorder(left)->separator);

    if (doingViterbi)
    {
        // trace backwards and clamp the best values
        for (int i=int(right); i>=int(left); i--)
        {
            Clique *cl = preorder(i);
            // clamp the values -- only care about non-separators
            if (!cl->separator)
                for (unsigned j=0; j<cl->discreteMember.size(); j++)
#ifdef NO_HASH_TABLE
                    cl->discreteMember[j]->val = Clique::gip[best].values[j];
#else
                    cl->discreteMember[j]->val = (*Clique::gip[best].values)[j];
#endif
            // set best to the best instantiation of the predecessor clique
            best = (*(preorder(i-1)->inum2gnum))[Clique::gip[best].pred];
        }
    }
    else
    {
        // now do the middle cliques, whose instantiations must pull in a lambda
        // from the instantiations derived from them, and push a lambda to the 
        // instantiations they derive from.
        for (int i=int(right); i>int(left); i-=2)
        {
            assert(copies[i] == 1);
            Clique *cl = preorder(i);
            assert(!cl->separator);
            for (int j=0; j<cl->numInstantiations; j++)
            {
                CliqueValue &cv = Clique::gip[cl->instantiation[j]];
                // retrieve the cached probGivenParents
                logpr t = cv.lambda;
	        // update non-separator clique value
                cv.lambda = 
                    Clique::gip[(*(preorder(i+1)->inum2gnum))[cv.succ]].lambda;
	        // update separator clique value
	        Clique::gip[(*(preorder(i-1)->inum2gnum))[cv.pred]].lambda += 
                    cv.lambda*t;
            }

            if (doingEM)
                incrementEMStatistics(cl);
        }
    }

    for (unsigned i=left+1; i<=right+1; i++)
    {
        copies[i]--;
        preorder(i)->reclaimMemory();
    }
}

void CliqueChain::finishBackwardPass()
{
// cout << "Finishing backward pass\n";
    // now do the root, which does not do any pushing
    if (preorderSize() > 1) 
    {
        Clique *cl = preorder(0);
        if (doingViterbi)
            // clamp the values of the root variable
            for (unsigned j=0; j<cl->discreteMember.size(); j++)
#ifdef NO_HASH_TABLE
                cl->discreteMember[j]->val = Clique::gip[best].values[j];
#else
                cl->discreteMember[j]->val = (*Clique::gip[best].values)[j];
#endif
        else
        {
            if (doingEM)
                incrementEMStatistics(cl);
        }
        preorder(0)->reclaimMemory();
        copies[0]--;
        preorder(1)->reclaimMemory();
        copies[1]--;
    }
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *     incrementEMStatistics() multiplies the lambdas and the pis for 
 *     each clique instantiation, clamps the nodes in the network, and
 *     increments the EM statistics for each node assigned to the clique.
 *     It assumes that the separators do not have any conditionalProbability
 *     nodes assigned to them. If it 
 *
 * Results:
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */

void CliqueChain::incrementEMStatistics(Clique *cl)
{
    // avoid numerical stability issues by using lambda-pi sum for
    // the current clique. This should be the same as dataProb computed
    // on the forward pass, but with lookup tables, etc...
    logpr cur_dp = 0.0;
    for (int k=0; k<cl->numInstantiations; k++)
    {
        CliqueValue &cv = Clique::gip[cl->instantiation[k]];
        cur_dp += cv.lambda*cv.pi;
    }

    // consider the contribution of each instantiation of the clique
    for (int k=0; k<cl->numInstantiations; k++)
    {
        CliqueValue &cv = Clique::gip[cl->instantiation[k]];
        // clamp the instantiation values
        for (unsigned j=0; j<cl->discreteMember.size(); j++)
#ifdef NO_HASH_TABLE
            cl->discreteMember[j]->val = cv.values[j];
#else
            cl->discreteMember[j]->val = (*cv.values)[j];
#endif

        // find the conditional families, given the switching parent values
        cl->findConditionalProbabilityNodes();

        // do the updates
        logpr posterior = cv.lambda*cv.pi/cur_dp;

        for (unsigned j=0; j<cl->conditionalProbabilityNode.size(); j++)
            cl->conditionalProbabilityNode[j]->emIncrement(posterior);
    }
}

/*
 * left is a separator, right is a non-separator.
*/

bool CliqueChain::recProcess(unsigned left, unsigned right)
{
// cout << "processing from " << left << " to " << right << "\n";
    assert(preorder(left)->separator);
    assert(!preorder(right)->separator);
    if (int(right)-int(left)+1 > baseCaseThreshold)  // must recurse
    {
// cout << "recursively processing from " << left << " to " << right << endl;
        // divide into k parts    
        int block_size = (int(right)-int(left)+1)/numSplits;
        // block size must be even to accomodate sep-nonsep pattern
        if (block_size % 2)  
            block_size = block_size+1;

        vector<int> leftboundary, rightboundary;
        leftboundary.push_back(left);
        int rb;
        while ((rb=leftboundary.back()+block_size-1) <= int(right))
        {
            rightboundary.push_back(rb);
            if (rb < int(right))
                leftboundary.push_back(rb+1);
            else if (rb == int(right))
                break;
        }
        if (!rightboundary.size() || rightboundary.back() < int(right))
            rightboundary.push_back(right);
        assert(leftboundary.size() == rightboundary.size());

        for (int i=0; i<int(leftboundary.size())-1; i++)
            advancePis(leftboundary[i], rightboundary[i], true);

        for (int i=int(leftboundary.size()-1); i>=0; i--)
            if (!recProcess(leftboundary[i], rightboundary[i]))
                return false;
    }
    else
    {
// cout << "entering base case" << endl;
        advancePis(left, right);
        if (right==cliques.size()-1)  // got to the end!!
        {
            if(!finishForwardPass())
                return false;
            initializeBackwardPass();
            recedeLambdas(right-2, left);
        }
        else
            recedeLambdas(right, left);
// cout << "exiting base case\n";
    }
    return true;
}

bool CliqueChain::compute(compute_mode mode, logpr _beam)
{
// cout << "computing\n";
// cout << "cliques range from 0 to " << (cliques.size()-1) << endl;

    map<int,int>::iterator mi;
    for (mi=copies.begin(); mi!=copies.end(); mi++)
        assert((*mi).second==0);  // all copies should have been recycled
    copies.clear();

    if (baseCaseThreshold < 2)
        error("The base case threshold must be at least 2.");
    if (numSplits < 2)
        error("There numst be at least 2 splits.");
    doingEM=doingViterbi=false;
    if (mode==viterbi)
        doingViterbi = true;
    else
        doingEM = true;
    beam = _beam;

    if (cliques.size()==1)  // a special case
    {
        initializeForwardPass();
        if (!finishForwardPass())
        {
            // clean up the memory usage and fail
            for (unsigned i=0; i<preorderSize(); i++)
                if (copies[i])
                {
                    preorder(i)->reclaimMemory();
                    copies[i]--;
                }
            return false;
        }
        initializeBackwardPass();
        finishBackwardPass();
    }
    else
    {    
        initializeForwardPass();
        if (!recProcess(1, cliques.size()-1))
        {
            // clean up the memory usage and fail
            for (unsigned i=0; i<preorderSize(); i++)
                if (copies[i])
                {
                    preorder(i)->reclaimMemory();
                    copies[i]--;
                }
            return false;
        }
        finishBackwardPass();
    }

    // memory processing should ensure that nothing is used at this point.
    assert(Clique::nextfree == int(Clique::gip.size())-1);

    return true;
}
