
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

#ifndef GMTK_GM_H
#define GMTK_GM_H

#include <set>

#include "logp.h"
#include "general.h"
#include "bp_range.h"

#include "GMTK_RandomVariable.h"
#include "GMTK_CliqueChain.h"



struct GMTK_GM
{
    vector<RandomVariable *> node;
    // This holds all the variables in the graph.
    // The topology os determined by the Parent and Child arrays associated
    // with each random variable.

    bool emMode;
    // when data probability is computed, this controls whether counts
    // are accumulated. If emMode==true, they are accumulated

    logpr dataProb, viterbiProb;
  
    vector<vector<VariableValue> > *example;
    // In EM there must be a way of iterating over the examples.
    // This store the examples to look at.
    // Note: if files are being used, the globalObservationMatrix stores the
    // data.

    // how far in the example array have we gotten?
    unsigned expos; 

    // Examples can be read from a file, or from an vector of vector of
    // VariableValues. This says which it is.
    bool using_files;

    CliqueChain *chain;
    // A pointer to a clique chain representation of the GM.

    vector<RandomVariable *> gmTemplate;

    int firstChunkFrame, lastChunkFrame;
    // the first and last frames in the template repeating segment

    int obsInTemplate, obsInRepeatSeg, framesInTemplate, framesInRepeatSeg;

    // support for iterating through a subset of the training segments
    BP_Range* trrng;  
    BP_Range::iterator* trrng_it; // iterator

    bool unrollCliqueChain;

    ///////////////////////////////////////////////////////////
    // The dividing line between member variables and member functions
    ///////////////////////////////////////////////////////////

    void verifyTopologicalOrder();
    // verifies that the node array is in topological order
    // copies it to the topological array

    void reveal(vector<RandomVariable *> order, bool show_vals = false);
    // Go through the nodes in the specified order and show them.

    GMTK_GM() {
      example=NULL; 
      chain=NULL; 
      using_files=false;
      trrng = NULL;
      trrng_it = NULL;
      unrollCliqueChain = true;  // appropriate for training and decoding
    }

    ~GMTK_GM() { 
      if (chain) delete chain; 
      deleteObsInVector(node);
      deleteObsInVector(gmTemplate);
      delete trrng;
      delete trrng_it;
    }

    void makeRandom();
    // Goes over each variable in the graph and calls its makeRandom function.

    void makeUniform();
    // Goes over each variable in the graph and calls its makeUniform function.

    void simulate();
    // This could be called 'sample' since it produces a sample from
    // the set of RVs constituting a GM.
    // Goes over the variables in topological order and instantiates them 
    // according to the probability distribution given their parent's values.

    void enumerateProb(int pos=0, logpr p=1.0);
    // A recursive function that enumerates all possible values of the 
    // hidden variables and computes the total data prob by brute force.
    // Good for verifying DP and other inference results.
    // Call with emMode=false to get the data probability.
    // Call with emMode=true and p=1/dataProb to do EM

    void cliqueChainProb(logpr beam=0.0)
    { chain->computePosteriors(beam); dataProb = chain->dataProb; }
    // Computes the likelihood of the observed variable values with
    // dynamic programming on a clique chain

    void cacheValues();
    // Tells all the nodes int the graph to cache the values they are 
    // currently set to.
    // Observation variables can ignore this.

    void restoreCachedValues();
    // Tells all the nodes in the graph to clamp themselves to their cached
    // values.
    // Observation variables can ignore this.

    void storeValues(vector<VariableValue> &vv);
    // stores the values of all the observed variables

    void setValues(vector<VariableValue> &vv);
    // sets the values of all the observed variables

    void enumerateViterbiProb(int pos=0, logpr p=1.0);
    // Computes the probability of the likeliest instantiation of the hidden
    // variables, and stores it in logViterbiProb. 
    // Has the side effect that at termination, the network is clamped to its
    // likeliest value.

    void cliqueChainViterbiProb(logpr beam=0.0) 
    {chain->doViterbi(beam); viterbiProb=chain->viterbiProb;}
    // Computes the probability of the likeliest instantiation of the hidden
    // variables, and stores it in logViterbiProb. 
    // Has the side effect that at termination, the network is clamped to its
    // likeliest value.

    /////////////////////////////////////////////////////////
    //  EM Support
    // Tells each variable to increment its counts by the posterior.
    void emIncrement(logpr posterior);
    ///////////////////////////////////////

    void setExampleStream(vector<vector<VariableValue> > *_example)
    {example=_example;}
    // The examples to be iterated over are held in here

    void setExampleStream(const char *const obs_file_name,
			  const char *const trrng_str = NULL);

    void clampFirstExample(); 
    // Clamps the observation variables according to the first example.

    bool clampNextExample(); 
    // Clamps the observation variables according to the next example.

    // Does EM using brute force inference.
    void enumerativeEM(int iterations);

    // Does EM using dynamic programming on a clique chain.
    void cliqueChainEM(const int iterations=1, 
		       const logpr beam=0.0,
		       const bool writeParametersAfterEachEMIteration=true,
		       const char* const outputParamFile=NULL,
		       const bool binOutFile=false,
		       const char* const loadAccFile=NULL,
		       const char* const loadAccRange=NULL,
		       const char* const storeAccFile=NULL,
		       const bool accFileIsBinary=true,
		       const char* const llStoreFile=NULL,
		       const double lldp=0.001);


    logpr enumerativeExampleProb(vector<vector<VariableValue > > &example);
    // computes the likelihood of the examples

    logpr cliqueChainExampleProb(vector<vector<VariableValue > > &example,
          logpr beam=0.0);
    // computes the likelihood of the examples

    void GM2CliqueChain();
    // this creates and initializes a clique chain corresponding to the GM

    void showCliques();
    // shows the cliques in preorder

    void unroll(int first_frame, int last_frame, int times);
    // duplicates the structure from first_frame to last_frame "times" times
    // e.g. unroll(1,2,5) adds 10 frames -- 5 1,2 chunks

    void cloneVariables(vector<RandomVariable *> &from, 
        vector<RandomVariable *> &to);
    // clones the network structure held in from so that it is reproduced in to.
    // all internal pointers are adjusted correctly

    void setSize(int repeat_segs);
    // sets the network up for infrence with an observation sequence that
    // is time_frames long

    void setupForVariableLengthUnrolling(int first_frame, int last_frame);
};

#endif
