
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

/*-
 * gmtkViterbi.cc
 *     Get the viterbi instantiation for each example in a data set 
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <iostream.h>
#include <fstream.h>
#include <map>
#include <string>

#include "general.h"
#include "error.h"
#include "rand.h"
#include "arguments.h"

#include "ieeeFPsetup.h"

#include "spi.h"

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"

VCID("$Header$");

// the file of observation

float pruneRatio=0.0;
char *obsFileName;
char *strFileName;
char *parmsFileName=NULL;
char *parmsPtrFileName=NULL;
bool showVitVals = false;

char *wordVar=NULL;
char *varMapFile=NULL;
char *transitionLabel=NULL;

ARGS ARGS::Args[] = {

 ARGS("obsFile",ARGS::Req,obsFileName,"File containing observations"),
 ARGS("parmsFile",ARGS::Opt,parmsFileName,"GM Parms File"), 
 ARGS("parmsPtrFile",ARGS::Opt,parmsPtrFileName,"GM Parms File"), 
 ARGS("strFile",ARGS::Req,strFileName,"GM Structure File"),
 ARGS("pruneRatio",ARGS::Opt,pruneRatio,"Pruning Ratio, values less than this*max are pruned"),
 ARGS("showVitVals",ARGS::Opt,showVitVals,"Print the viterbi values??"),

 // These 3 must be used together or not at all
 ARGS("printWordVar",ARGS::Opt,wordVar,"Print the word var - which has this label"),
 ARGS("varMap",ARGS::Opt,varMapFile,"Use this file to map from word-index to string"),
 ARGS("transitionLabel",ARGS::Opt,transitionLabel,"The label of the word transition variable"),

 ARGS()

};


RAND rnd(0);
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;


int
main(int argc,char*argv[])
{
  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();

  ARGS::parse(argc,argv);

  map<int, string> word_map;
  if (wordVar != NULL)
  {
    if (varMapFile==NULL)
      error("File to map from word index to string not specified.");
    if(transitionLabel == NULL)
      error("The label of the transition variable was not specified.");
    ifstream in(varMapFile);
    if (!in) { cout << "Unable to open " << varMapFile << endl; exit(1); }
    string name;
    int val;
    while (!in.eof())
    {
      in >> val >> name >> ws;
      word_map[val] = name;
    }
    in.close();
  }

  /////////
  // read in all parameters
  if (parmsFileName) {
    iDataStreamFile parmsFile(parmsFileName);
    GM_Parms.readAll(parmsFile);
  }
  if (parmsPtrFileName) {
    iDataStreamFile parmsFile(parmsPtrFileName);
    GM_Parms.read(parmsFile);
  }

  /////////////////////////////
  // read in the structure

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName);
  // parse the file
  fp.parseGraphicalModel();
  // create the rv variable objects
  fp.createRandomVariableGraph();
  // make sure that there are no directed loops in the graph
  // by imposing the S,SE,E,NE constrains
  fp.ensureS_SE_E_NE();
  // link the RVs with the parameters that are contained in
  // the bn1_gm.dt file.
  fp.associateWithDataParams();

  // now associate the RVs with a GM
  GMTK_GM gm;
  fp.addVariablesToGM(gm);

  gm.setExampleStream(obsFileName);

  gm.verifyTopologicalOrder();

  gm.GM2CliqueChain();
  // gm.showCliques();

  gm.setupForVariableLengthUnrolling(fp.firstChunkFrame(),fp.lastChunkFrame());

  // and away we go
  gm.clampFirstExample();
  do
  {
    gm.cliqueChainViterbiProb(pruneRatio);
    cout << "Example prob: " << gm.viterbiProb.val() << " : "
         << ((*gm.node.rbegin())->timeIndex+1) << " frames\n";
    if (showVitVals)
      gm.reveal(gm.node, true);
    if (wordVar && gm.viterbiProb!=0.0)
    {
      // print the sequence of values for this variable
      // compress consecutive values into a single instance
      // the times are right if a word transition at time t means there is
      // a new word at t+1
      string pvn = string(wordVar);
      string tl = string(transitionLabel);
      for (int i=0, lv=-1, lf=0; i<int(gm.node.size()); i++)
      {
        if (gm.node[i]->label == pvn)
        {
          if (!gm.node[i]->discrete) 
            error("Can only print Viterbi values for discrete variables");
          if (gm.node[i]->cardinality != int(word_map.size()))
            error("Word-val to string map does not match the number of words.");
          lv = gm.node[i]->val;
        }
        else if (gm.node[i]->label==tl)
        {
          if (gm.node[i]->cardinality != 2) 
            error("Word transition variable should have two values");
          if (gm.node[i]->val==1)  // a word transition
          {
            cout << word_map[lv] << " (" << lf << "-" 
                 << gm.node[i]->timeIndex << ")\n";  
            lf = gm.node[i]->timeIndex+1;
          }
        }
      }
    }
  } while (gm.clampNextExample());

  return 0;  
}
