/*
 * "Copyright 2001, University of Washington and International Business Machines Corporation. All Rights Reserved
 *
 *    Written by Jeff Bilmes and Geoffrey Zweig
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
 * JEFF BILMES AND GEOFFREY ZWEIG SHALL NOT HAVE ANY LIABILITY FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING WITHOUT LIMITATION LOST PROFITS), HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE  OF
 * THE PROGRAM, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES."
*/

/*
 * This program converts from ascii trainable parameters to binary
 * and vice versa.
 *
 *
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
#include "spi.h"

VCID("$Header$");

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"


/*
 * command line arguments
 */

char *prmTrainableFile=NULL;
bool binPrmTrainableFile=false;

char *prmMasterFile=NULL;
char *prmOutFile="outParms%d.gmp";
bool binPrmOutFile=true;

double varFloor = 1e-10;
char *cppCommandOptions = NULL;


void makeArgs(Argument_List &args)
{
  bool optional=0,required=1;

  args.add("prmMasterFile",optional,&prmMasterFile,
           "Multi-level master CPP processed GM Parms File");

  args.add("prmTrainableFile",required,&prmTrainableFile,
           "File containing Trainable Parameters");
  args.add("binPrmTrainableFile",optional,&binPrmTrainableFile,
           "Is Binary? File containing Trainable Parameters");


  args.add("prmOutFile",required,&prmOutFile,
           "File to place *TRAINABLE* output parametes");
  args.add("binPrmOutFile",optional,&binPrmOutFile,
           "Output parametes binary? (def=false)");
  args.add("cppCommandOptions",optional,&cppCommandOptions,
           "Command line options to give to cpp");

  args.add("varFloor",optional,&varFloor,
           "Variance Floor");
  args.add("floorVarOnRead",optional,&DiagCovarVector::floorVariancesWhenReadIn,
           "Floor the variances to varFloor when they are read in");
  args.add("cptNormThreshold",optional,&CPT::normalizationThreshold,
           "Read error if |Sum-1.0|/card > norm_threshold");
}

/*
 * definition of needed global arguments
 */
RAND rnd(false);
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

  ////////////////////////////////////////////
  // parse arguments
  Argument_List args;
  makeArgs(args);
  args.parse(argc, argv);

  MixGaussiansCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();

  ////////////////////////////////////////////
  // set global variables/change global state from args
  GaussianComponent::setVarianceFloor(varFloor);
  /////////////////////////////////////////////


  ////////////////////////////////////////////
  // finally, pull any trainable parameters out
  // of a master file and send them to a trainable file. 
  if (prmMasterFile != NULL) {
    iDataStreamFile pf(prmMasterFile,false,true,cppCommandOptions);
    GM_Parms.read(pf);
  }
  iDataStreamFile pf(prmTrainableFile,binPrmTrainableFile,true,cppCommandOptions);
  GM_Parms.readTrainable(pf);
  printf("Trainable file '%s' has '%u' total parameters\n",
	 prmTrainableFile,
	 GM_Parms.totalNumberParameters());

  oDataStreamFile of(prmOutFile,binPrmOutFile);
  GM_Parms.writeTrainable(of);

  exit_program_with_status(0);
}

