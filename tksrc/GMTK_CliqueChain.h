
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

#include <list>
#include "GMTK_Clique.h"

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

    bool forwardPass(logpr beam=0.0, bool viterbi=false);
    // Computes the alpha probabilities and/or viterbi clique value pointers.
    // Prunes away entries that are less than beam*max.

    void backwardPass();
    // In the backward pass, lambdas are computed for each of the CliqueValues
    // stored in the forward pass. When these are multiplied with the 
    // corresponding pis, and divided by the data prob, the poserior of each
    // clique instantiation results.

    bool doViterbi(logpr beam=0.0);
    // Computes the likeliest value of each clique, and clamps the variables
    // in the underlying network correspondingly.
    // Prunes away entries that are less than beam*max.
    // returns false of the forwards pass fails (0 prob after pruning)

    bool computePosteriors(logpr beam=0.0);
    // Calculates the lambdas and pis for all the cliques.
    // Prunes away entries that are less than beam*max.
    // returns false of the forwards pass fails (0 prob after pruning)

    void incrementEMStatistics();
    // Multiplies the lambdas and pis for each clique instantiation, clamps
    // the network, and increments the EM statistics for each node assigned
    // to a clique.

    logpr dataProb, viterbiProb, backwardDataProb;
};

#endif
