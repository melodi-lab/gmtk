
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
#include "version.h"

#include "ieeeFPsetup.h"

#include "spi.h"

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"

#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"

VCID("$Header$");

// the file of observation

float beam=-LZERO;
char *strFileName;
char *prmMasterFile=NULL;
char *prmTrainableFile=NULL;
bool binPrmTrainableFile=false;
bool showVitVals = false;

// observation file support

char *obsFileName;

#define MAX_NUM_OBS_FILES (3)
char *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL };
unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
char *frs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
char *irs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
char *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile" };
bool iswps[MAX_NUM_OBS_FILES] = { false, false, false };

int startSkip = 0;
int endSkip = 0;
char *dcdrng_str="all";

char *wordVar=NULL;
char *varMapFile=NULL;
char *transitionLabel=NULL;
double varFloor = GMTK_DEFAULT_VARIANCE_FLOOR;

char *ofilelist = NULL;
char *dumpNames = NULL;

int show_cliques=0;

char *cppCommandOptions = NULL;

int bct=GMTK_DEFAULT_BASECASETHRESHOLD;
int ns=GMTK_DEFAULT_NUM_SPLITS;

bool print_version_and_exit = false;

Arg Arg::Args[] = {

  Arg("of1",Arg::Req,ofs[0],"Observation File 1"),
  Arg("nf1",Arg::Opt,nfs[0],"Number of floats in observation file 1"),
  Arg("ni1",Arg::Opt,nis[0],"Number of ints in observation file 1"),
  Arg("fr1",Arg::Opt,frs[0],"Float range for observation file 1"),
  Arg("ir1",Arg::Opt,irs[0],"Int range for observation file 1"),
  Arg("fmt1",Arg::Opt,fmts[0],"Format (htk,bin,asc,pfile) for observation file 1"),
  Arg("iswp1",Arg::Opt,iswps[0],"Endian swap condition for observation file 1"),


  Arg("of2",Arg::Opt,ofs[1],"Observation File 2"),
  Arg("nf2",Arg::Opt,nfs[1],"Number of floats in observation file 2"),
  Arg("ni2",Arg::Opt,nis[1],"Number of ints in observation file 2"),
  Arg("fr2",Arg::Opt,frs[1],"Float range for observation file 2"),
  Arg("ir2",Arg::Opt,irs[1],"Int range for observation file 2"),
  Arg("fmt2",Arg::Opt,fmts[1],"Format (htk,bin,asc,pfile) for observation file 2"),
  Arg("iswp2",Arg::Opt,iswps[1],"Endian swap condition for observation file 2"),
  Arg("of3",Arg::Opt,ofs[2],"Observation File 3"),
  Arg("nf3",Arg::Opt,nfs[2],"Number of floats in observation file 3"),
  Arg("ni3",Arg::Opt,nis[2],"Number of ints in observation file 3"),
  Arg("fr3",Arg::Opt,frs[2],"Float range for observation file 3"),
  Arg("ir3",Arg::Opt,irs[2],"Int range for observation file 3"),
  Arg("fmt3",Arg::Opt,fmts[2],"Format (htk,bin,asc,pfile) for observation file 3"),
  Arg("iswp3",Arg::Opt,iswps[2],"Endian swap condition for observation file 3"),

  Arg("strFile",Arg::Req,strFileName,"GM Structure File"),
  Arg("inputMasterFile",Arg::Req,prmMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
  Arg("inputTrainableFile",Arg::Opt,prmTrainableFile,"File of only and all trainable parameters"),
  Arg("binInputTrainableFile",Arg::Opt,binPrmTrainableFile,"Binary condition of trainable parameters file"),
  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Command line options to give to cpp"),


  Arg("varFloor",Arg::Opt,varFloor,"Variance Floor"),
  Arg("floorVarOnRead",Arg::Opt,DiagCovarVector::floorVariancesWhenReadIn,
       "Floor the variances to varFloor when they are read in"),
  Arg("dcdrng",Arg::Opt,dcdrng_str,"Range to decode over segment file"),

  Arg("beam",Arg::Opt,beam,"Beam width (less than max*exp(-beam) are pruned away)"),
  Arg("showVitVals",Arg::Opt,showVitVals,"Print the viterbi values??"),

  // These 3 must be used together or not at all
  Arg("printWordVar",Arg::Opt,wordVar,"Print the word var - which has this label"),
  Arg("varMap",Arg::Opt,varMapFile,"Use this file to map from word-index to string"),
  Arg("transitionLabel",Arg::Opt,transitionLabel,"The label of the word transition variable"),
  Arg("startSkip",Arg::Opt,startSkip,"Frames to skip at beginning (i.e., first frame is buff[startSkip])"),
  Arg("endSkip",Arg::Opt,endSkip,"Frames to skip at end (i.e., last frame is buff[len-1-endSkip])"),
  Arg("cptNormThreshold",Arg::Opt,CPT::normalizationThreshold,"Read error if |Sum-1.0|/card > norm_threshold"),

  Arg("showCliques",Arg::Opt,show_cliques,"Show the cliques after the netwok has been unrolled k times."),


  Arg("dumpNames",Arg::Opt,dumpNames,"File containing the names of the variables to save to a file"),
  Arg("ofilelist",Arg::Opt,ofilelist,"List of filenames to dump the hidden variable values to"),

 Arg("numSplits",Arg::Opt,ns,"Number of splits to use in logspace recursion (>=2)."),

  Arg("baseCaseThreshold",Arg::Opt,bct,"Base case threshold to end recursion (>=2)."),

  Arg("version",Arg::Opt,print_version_and_exit,"Print GMTK version number and exit."),

  Arg()

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
  set_new_handler(memory_error);

  Arg::parse(argc,argv);

  if (print_version_and_exit)
    printf("%s\n",gmtk_version_id);

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

  if (dumpNames)
    if (ofilelist==NULL) error("Must specify output files for binary dumping");

  ////////////////////////////////////////////
  // check for valid argument values.
  int nfiles = 0;
  unsigned ifmts[MAX_NUM_OBS_FILES];
  for (int i=0;i<MAX_NUM_OBS_FILES;i++) {
    if (ofs[i] != NULL && nfs[i] == 0 && nis[i] == 0)
      error("ERROR: command line must specify one of nf%d and ni%d not zero",
            i+1,i+1);
    nfiles += (ofs[i] != NULL);
    if (strcmp(fmts[i],"htk") == 0)
      ifmts[i] = HTK;
    else if (strcmp(fmts[i],"binary") == 0)
      ifmts[i] = RAWBIN;
    else if (strcmp(fmts[i],"ascii") == 0)
      ifmts[i] = RAWASC;
    else if (strcmp(fmts[i],"pfile") == 0)
      ifmts[i] = PFILE;
    else
      error("ERROR: Unknown observation file format type: '%s'\n",fmts[i]);
  }

  globalObservationMatrix.openFiles(nfiles,
				    (const char**)&ofs,
				    (const char**)&frs,
				    (const char**)&irs,
				    (unsigned*)&nfs,
				    (unsigned*)&nis,
				    (unsigned*)&ifmts,
				    (bool*)&iswps,
				    startSkip,
				    endSkip);

  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();
  if (beam < 0.0)
    error("beam must be >= 0");
  if (startSkip < 0 || endSkip < 0)
    error("startSkip/endSkip must be >= 0");

  ////////////////////////////////////////////
  // set global variables/change global state from args
  GaussianComponent::setVarianceFloor(varFloor);


  /////////////////////////////////////////////
  // read in all the parameters
  if (prmMasterFile) {
    // flat, where everything is contained in one file, always ASCII
    iDataStreamFile pf(prmMasterFile,false,true,cppCommandOptions);
    GM_Parms.read(pf);
  }
  if (prmTrainableFile) {
    // flat, where everything is contained in one file, always ASCII
    iDataStreamFile pf(prmTrainableFile,binPrmTrainableFile,true,cppCommandOptions);
    GM_Parms.readTrainable(pf);
  }
  GM_Parms.loadGlobal();

  /////////////////////////////
  // read in the structure

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName,cppCommandOptions);

  printf("Finished reading in all parameters and structures\n");

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
  // make sure that all observation variables work
  // with the global observation stream.
  fp.checkConsistentWithGlobalObservationStream();

  // now associate the RVs with a GM
  GMTK_GM gm(&fp);
  fp.addVariablesToGM(gm);

  gm.setExampleStream(obsFileName,dcdrng_str);
  GM_Parms.checkConsistentWithGlobalObservationStream();

  GM_Parms.setFirstUtterance( gm.trrng->min() );
 
  gm.setCliqueChainRecursion(ns, bct);

  gm.verifyTopologicalOrder();

  gm.GM2CliqueChain();
  gm.setupForVariableLengthUnrolling(fp.firstChunkFrame(),fp.lastChunkFrame());

  if (show_cliques)
  {
    cout << "The cliques in the unrolled network are:\n";
    gm.setSize(show_cliques);
    gm.showCliques();
  }

  set<string> dumpVars;
  map<string,int> posFor;
  int ndv=0;
  if (dumpNames)
  {
      ifstream din(dumpNames);
      if (!din) {cout << "Unable to open " << dumpNames << endl; exit(1);}
      string s;
      while (!din.eof())
      {
          din >> s >> ws;
          dumpVars.insert(s);
          posFor[s] = ndv++;  // which position within a frame to put it
      }
      din.close();
  }

  vector<string> ofiles;
  if (ofilelist)
  {
      ifstream oin(ofilelist);
      if (!oin) {cout << "Unable to open " << ofilelist << endl; exit(1);} 
      string s;
      while (!oin.eof())
      {
          oin >> s >> ws;
          ofiles.push_back(s);
      }
      oin.close();
  }

  // and away we go
  int ne=0;
  gm.clampFirstExample();
  do {
    logpr pruneRatio;
    pruneRatio.valref() = -beam;
    gm.cliqueChainViterbiProb(pruneRatio);
    cout << "Example prob: " << gm.viterbiProb.val() << " : "
         << ((*gm.node.rbegin())->timeIndex+1) << " frames\n" << flush;
    if (showVitVals)
      gm.reveal(gm.node, true);
    if (dumpNames)
    {
        if (ne==int(ofiles.size())) error("More utterances than output files");
        FILE *fp = fopen(ofiles[ne++].c_str(), "wb");
        if (!fp) {cout << "Unable to open " << ofiles[ne-1] << endl; exit(1);}
        int *vals = new int[gm.node.size()];
        // just in case a variable isn't present in a slice, write a -1
        for (int i=0; i<int(gm.node.size()); i++) vals[i] = -1;
        int nv=0;
        for (int i=0; i<int(gm.node.size()); i++)
            if (dumpVars.count(gm.node[i]->label)) 
            {
                if (!gm.node[i]->discrete || !gm.node[i]->hidden) 
                    error("variables to dump must be discrete and hidden");
                int p = gm.node[i]->timeIndex*dumpVars.size() 
                        + posFor[gm.node[i]->label];
                assert(p<int(gm.node.size()));
                vals[p] = gm.node[i]->val; 
                nv++;
            }
        fwrite((void *) vals, sizeof(int), nv, fp); 
        delete [] vals;
        fclose(fp);
    }
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
                 << gm.node[i]->timeIndex << ")\n" << flush;  
            lf = gm.node[i]->timeIndex+1;
          }
        }
      }
    }
  } while (gm.clampNextExample());

  exit_program_with_status(0);
}
