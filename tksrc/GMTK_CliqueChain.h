
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

#ifndef CLIQUECHAIN_H
#define CLIQUECHAIN_H

#include "GMTK_Clique.h"

enum compute_mode {statistics,viterbi,data_probs};

struct CliqueChain
{
    vector<Clique> cliques;
    // The actual cliques in the chain.

    vector<Clique *> _preorder;
    // Pointers to the cliques in pre and post-order.

    Clique* preorder(const int i) { return _preorder[i]; }
    unsigned preorderSize() { return _preorder.size(); }
    Clique* postorder(const int i) 
         { return _preorder[_preorder.size()-i-1]; }
    unsigned postorderSize() { return _preorder.size(); }

    CliqueChain() {numSplits=3; baseCaseThreshold=10;}

    void setRecursion(int ns, int bct) {numSplits=ns; baseCaseThreshold=bct;}

    void initializeForwardPass();

    bool finishForwardPass();

    void advancePis(unsigned first_clique, unsigned last_clique, 
    bool erase=false);

    void initializeBackwardPass();

    void recedeLambdas(unsigned right, unsigned left);

    void finishBackwardPass();
   
    void incrementEMStatistics(Clique *cl);

    bool recProcess(unsigned left, unsigned right);

    bool compute(compute_mode mode, logpr _beam);

    logpr dataProb, viterbiProb, backwardDataProb;

    int best;  // used in the Viterbi computation to store the best predecessor

    bool doingEM, doingViterbi;

    logpr beam;

    int baseCaseThreshold, numSplits;
};

#endif
