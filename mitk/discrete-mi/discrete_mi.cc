/** 
 * discrete_mi.cc

 Main functions to calculate the mutual information between 
 two set of discrete-valued vectors.
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <errno.h>
#include <cstring>
#include <values.h>
#include <cmath>
#include <cassert>
#include <ctime>
#include <signal.h>

#include "general.h"
#include "error.h"
#include "range.h"
#include "mixNormal.h"
#include "mixNormalCollection.h"
#include "readRange.h"
#include "GMTK_ObservationMatrix.h"
#include "arguments.h"

#define BUFFER_SIZE 300
#define FREQUENCY 10

#define MAX_LABEL_VAL 16384

#define MIN(a,b) ((a)<(b)?(a):(b))

static const char* program_name;
ObservationMatrix globalObservationMatrix;
ObservationMatrix globalLabelMatrix;
bool gotLabelFile = false;
int min_consec_labels_per_sentence = 0;


#define MAX_NUM_OBS_FILES (3)
const char *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL }; 
unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
const char *frs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
const char *irs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
const char *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile" };
bool iswps[MAX_NUM_OBS_FILES] = { false, false, false };


const char *lofs[] = { NULL}; 
unsigned lnfs[] = { 0};
unsigned lnis[] = { 0};
const char *lfrs[] = { "all"};
const char *lirs[] = { "all"};
const char *lfmts[] = { "pfile"};
bool liswps[] = {false };

//ARGS ARGS::Args[] = {
//  ARGS()
//
//}; 


int readFeatures(Range::iterator krit, size_t &n_frames,
		 size_t &n_samps,
		 Range &lrrng, int labpos,
		 unsigned &frameStart, unsigned &firstFrame);


/**
 * the static usage of this program
 *
 * @param message the message to be outputted
 */
static void usage(const char* message = 0) {
  if (message)
    fprintf(stderr, "%s: %s\n", program_name, message);
  fprintf(stderr, "Usage: %s <options>\n", program_name);
  fprintf(stderr,
	  "Where <options> include:\n"
	  "[-help]           print this message\n"
	  "-i <file-name>  input pfile\n"
	  "-ni1 <number of ints> for file 1\n"
	  "-nf1 <number of floats> for file 1\n"
	  "-ir1 <int range> for file 1\n"
	  "-fr1 <float range> for file 1\n"
	  "-fmt1 <format(htk,bin,asc,pfile)> for file 1\n"
	  "[-i2 <file-name>]  optional second input file to be merged with first.  The total number od discrete features has to be 1.\n"
	  "[-ni2 <number of ints>] for file 2\n"
	  "[-nf2 <number of floats>] for file 2\n"
	  "-ir2 <int range> for file 2\n"
	  "-fr2 <float range> for file 2\n"
	  "-fmt2 <format(htk,bin,asc,pfile)> for file 2\n"
	  "-iswap1         byte-swap input file 1\n"
	  "-iswap2         byte-swap input file 2\n"
	  "-lswap          byte-swap label file\n"
	  "[-o <file-name>]  output pfile (default stdout) \n"
	  "[-rf <file-name>]  range file. Default: rangeFile.dat\n"
	  "[-q]              quiet mode\n"
	  "[-sr <range>]     sentence range. Default: all\n"
	  "[-fr <range>]     feature range. Default: all\n"
	  "[-lr <range>]     label range: Subset of labels to condition on. Default = all.\n\n"
	  "[-lf file]      Pfile containing labels. Use only with none-all -lr option.\n"
	  "-labfmt <format(htk,bin,asc,pfile)> for the label file\n"
	  "[-labpos int] Position of the label.  Default: last discrete feature.\n"
	  "-b              print raw binary data (ints and floats)\n"
	  "-ns             Don't print the frame IDs (i.e., sent and frame #)\n"
	  "-debug <level>  number giving level of debugging output to produce 0=none.\n"
	  "[-mlps d]    Minimum consecutive labels per segment to compute CMI with.\n"
  );
  exit(EXIT_FAILURE);
}


/*
Reads features and labels and performs error checking
 */

#define NO_DATA 0
#define DATA_LEFT 1
#define DONE 2

int readFeatures(Range::iterator krit, size_t &n_frames,
		 size_t &n_samps,
		 Range &lrrng, int labpos,
		 unsigned &frameStart,unsigned &firstFrame) {

  bool segAlreadyLoaded = false;
  if(frameStart != 0) segAlreadyLoaded = true;

  size_t n_labs;
  n_labs = (gotLabelFile?globalLabelMatrix.numDiscrete():globalObservationMatrix.numDiscrete());
  const bool stateCondMI = !lrrng.full();

  if (stateCondMI && n_labs < 1) { 
    error("For conditional MI, number of discrete features per frame in pfile must be at least one.");
  }
  
  ObservationMatrix* obsMat;
  if(gotLabelFile) 
    obsMat = &globalLabelMatrix;
  else
    obsMat = &globalObservationMatrix;
  if(!segAlreadyLoaded)
    obsMat->loadSegment((const unsigned)(*krit));
  n_frames = obsMat->numFrames();
  
  if ( n_frames == SIZET_BAD )
    error("%s couldn't find number of frames at sentence %lu in input pfile.\n", program_name, (unsigned long) *krit);

  if (!stateCondMI) { 
    firstFrame = 0;
    n_samps = n_frames;
    if(n_samps == 0) return NO_DATA;
    return DONE;
  }
  else {
    unsigned label;
    unsigned frameno;
    int pos = (unsigned) labpos;
    if(pos == -1) //use default: last discrete feature
      pos =  obsMat->numFeatures() - 1;
    
    int numFound = 0;
    for(frameno = frameStart; frameno < n_frames; ++frameno) {
      label =  obsMat->unsignedAtFrame(frameno,(const unsigned) pos);
      frameStart = frameno;
      while (lrrng.contains(label)) {
	++numFound;
	++frameno; 
	if(frameno >= n_frames) break;
	label =  obsMat->unsignedAtFrame(frameno,(const unsigned) pos);
      }
      if (numFound == 0 || numFound < min_consec_labels_per_sentence) 
	numFound = 0; //reset
      else break;
    }
    //At this point we have read a continuous chunk
    //Case 1:  We haven't reached the last frame
    //         Everything is good, we return DATA_LEFT
    //Case2:   We've reached the last frame:
    //   Case 2.a: numFound > 0 : return DONE   
    //   Case 2.b: numFound ==0 : return NO_DATA
    //      
    n_samps = numFound;
    firstFrame = frameStart;  
    if(frameno != n_frames) {
      if(!segAlreadyLoaded) {
	globalObservationMatrix.loadSegment((const unsigned)(*krit));
	if(gotLabelFile) { //check that they have the same number of frames
	  unsigned obs_n_frames = globalObservationMatrix.numFrames();
	  if(obs_n_frames > n_frames)
	    warning("There are fewer observation frames than labels\n");
	  else if(obs_n_frames < n_frames)
	    error("There are more observation frames than labels\n");
	}
      }
      frameStart = frameno;  //next time we start from here
      return DATA_LEFT;
    }    
    else if(numFound > 0) {
      if(!segAlreadyLoaded) {
	globalObservationMatrix.loadSegment((const unsigned)(*krit));
	if(gotLabelFile) { //check that they have the same number of frames
	  unsigned obs_n_frames = globalObservationMatrix.numFrames();
	  if(obs_n_frames > n_frames)
	    warning("There are fewer observation frames than labels\n");
	  else if(obs_n_frames < n_frames)
	    error("There are more observation frames than labels\n");
	}
      }
      return DONE;
    }
    else { //numFound == 0
      return NO_DATA;
    }
  } //end else (stateCondMI)
  return DONE;
}


//Iterates over all sentences accumulating the fequency of diff vectors
static void discreteMI(FILE *mi_ofp, // where to store output MI values
		       const char* range_fname,
		       Range &srrng,
		       Range &frrng,
		       Range &lrrng,
		       const bool print_frameid,
		       const bool quiet,
		       const bool binary,
		       int labpos) {
  size_t n_frames, n_samps;
  size_t totalNumSamples = 0;
  size_t n_discrete  = globalObservationMatrix.numDiscrete();
  int readStatus;
  unsigned frameStart,firstFrame;  

  //The label will be included with the discret features.
  //One just needs to be careful to have a correct range specification
  //if(!lrrng.full() && !gotLabelFile) {
  //  n_discrete--;  //because the label is one of the discrete features
                   // and we don't wanna include it
  //}

  cout<<"Starting program...\n";
  RangeSetCollection rngSetCol(range_fname);
  cout<<"Creating collection...\n";
  MixNormalCollection mg(rngSetCol);
  
  //Start
  cout<<"Starting MI calculation...\n";
  mg.startEpoch();   
  for ( Range::iterator srit = srrng.begin(); !srit.at_end(); srit++ ) {
    if ( *srit % FREQUENCY == 0 )
      std::cout << "Processing sentence " << *srit << std::endl;
     frameStart = 0;
     do {
       readStatus = readFeatures(srit, n_frames,n_samps,lrrng,labpos,frameStart,firstFrame);
       if(readStatus == NO_DATA) break;
       totalNumSamples += n_samps;
              mg.addToEpoch(&globalObservationMatrix, n_discrete, n_frames, n_samps, firstFrame, rngSetCol);
     } while(readStatus == DATA_LEFT);
  }   

  cout<<"Finished MI calculation.\n";
  mg.endEpoch( totalNumSamples, rngSetCol, mi_ofp);
  //End

}

/**
 * parse long variable from string
 *
 * @param s the string containing the long value
 * @return the value corresponding to s
 */
static long parse_long(const char*const s) {
  size_t len = strlen(s);
  char *ptr;
  long val;

  val = strtol(s, &ptr, 0);

  if (ptr != (s+len))
    error("Not an integer argument.");

  return val;
}

/**
********************  MAIN  *************
*/
int main(int argc, const char *argv[]) {
  const char *input_fname = 0;  // Input file name.
  const char *input_fname2 = 0;  // Input file name 2.
  const char *output_fname = 0; // Output file name.
  const char *label_fname = 0;   // Label file name (if any).
  char *range_fname = 0; //  range name.
  
  const char *sr_str = 0;   // sentence range string
  Range *sr_rng;
  const char *fr_str = 0;   // feature range string    
  Range *fr_rng;
  const char *lr_str = 0;   // label range string    
  Range *lr_rng;

  int debug_level = 0;
  bool print_frameid = true;
  bool quiet = false;
  bool binary=false;

  bool gotLabs = false;
  bool lswap = false;
  int labpos = -1;  //default: last postion used as the label

  program_name = *argv++;
  argc--;

  while ( argc-- ) {
    char buf[BUFSIZ];
    const char *argp = *argv++;

    if ( strcmp(argp, "-help") == 0 ) {
      usage();
    } else if ( strcmp(argp, "-ns") == 0 ) {
      print_frameid = false;
    } else if ( strcmp(argp, "-q") == 0 ) {
      quiet = true;
    } else if ( strcmp(argp, "-b") == 0 ) {
      binary = true;
    } else if ( strcmp(argp, "-i") == 0 ) {
      // Input file name.
      if ( argc > 0 ) {
	// Next argument is input file name.
	input_fname = *argv++;
	argc--;
      } else
	usage("No input filename given.");
    } else if (strcmp(argp, "-ni1")==0) {
      if (argc>0) {
	nis[0] = (int) parse_long(*argv++);
	argc--;
      } else
	usage("No number of ints specified.");
    } else if (strcmp(argp, "-nf1")==0) {
      if (argc>0) {
	nfs[0] = (int) parse_long(*argv++);
	argc--;
      } else usage("No number of floats specified.");
    } else if (strcmp(argp, "-fr1")==0) {
      if (argc>0) {
	frs[0] = *argv++;
	argc--;
      } else
	usage("No float range given.");
    } else if (strcmp(argp, "-ir1")==0) {
      if (argc>0) {
	irs[0] = *argv++;
	argc--;
      } else
	usage("No int range given for file 1.");
    } else if (strcmp(argp, "-fmt1")==0) {
      if (argc>0) {
	fmts[0] = *argv++;
	argc--;
      } else usage("No format given for file 1.");
    } else if (strcmp(argp, "-ni2")==0) {
      if (argc>0) {
	nis[1] = (int) parse_long(*argv++);
	argc--;
      } else
	usage("No number of ints specified.");
    } else if (strcmp(argp, "-nf2")==0) {
      if (argc>0) {
	nfs[1] = (int) parse_long(*argv++);
	argc--;
      } else
	usage("No number of floats specified.");
    } else if (strcmp(argp, "-ir2")==0) {
      if (argc>0) {
	irs[1] = *argv++;
	argc--;
      } else
	usage("No int range given for file 2.");
    } else if (strcmp(argp, "-fr2")==0) {
      if (argc>0) {
	frs[1] = *argv++;
	argc--;
      } else usage("No float range given.");
    } else if (strcmp(argp, "-fmt2")==0) {
      if (argc>0) {
	fmts[1] = *argv++;
	argc--;
      } else usage("No format given for file 2.");
    } else if (strcmp(argp,"-iswap1") == 0) {
      iswps[0] = true;
    } else if (strcmp(argp,"-iswap2") == 0) {
      iswps[1] = true;
    } else if (strcmp(argp,"-lswap") == 0) {
      lswap = true;
    } else if ( strcmp(argp, "-o") == 0 ) {
      // Output file name.
      if ( argc > 0 ) {
	// Next argument is output file name.
	output_fname = *argv++;
	argc--;
      } else
	usage("No output filename given.");
    }  else if ( strcmp(argp, "-debug") == 0 ) {
      if ( argc > 0 ) {
	// Next argument is debug level.
	debug_level = (int) parse_long(*argv++);
	argc--;
      } else
	usage("No debug level given.");
    } else if ( strcmp(argp, "-sr" )==0) {
      if (argc>0) {
	sr_str = *argv++;
	argc--;
      } else
	usage("No range given.");
    } else if ( strcmp(argp, "-fr")==0 ) {
      if (argc>0) {
	fr_str = *argv++;
	argc--;
      } else
	usage("No range given.");
    } else if ( strcmp(argp, "-lr")==0 ) {
      if (argc>0) {
	gotLabs = true;
	lr_str = *argv++;
	argc--;
      } else
	usage("No range given.");
    } else if (strcmp(argp, "-lf")==0) {
      if (argc>0) {
	label_fname = *argv++;
	argc--;
      }
      else
	usage("No label file name given.");
    }  else if (strcmp(argp, "-labpos")==0) {
      if (argc>0) {
	labpos = parse_long(*argv++);
	argc--;
      } else
	usage("No -labpos *d* value given.");
    } else if (strcmp(argp, "-labfmt")==0) {
      if (argc>0) {
	fmts[2] = *argv++;
	argc--;
      }
      else
	usage("No format given for label file.");
    } else if ( strcmp(argp, "-rf") == 0 ) {
      // range file name.
      if ( argc > 0 ) {
	// Next argument is range file name.
	range_fname = (char*) *argv++;
	argc--;
      } else
	usage("No range filename given.");
    }  else if (strcmp(argp, "-mlps")==0) {
      if (argc>0) {
	min_consec_labels_per_sentence = parse_long(*argv++);
	argc--;
      } else
	usage("No -mlps *d* value given.");
      if (min_consec_labels_per_sentence < 0)
	error("-mlps argument must be positive");
    } else {
      sprintf(buf,"Unrecognized argument (%s).",argp);
      usage(buf);
    }
  }

  // Check all necessary arguments provided before creating objects.
   int fileno = 0;
  if (input_fname == NULL)
    usage("No input file name supplied");
  
  if (input_fname2 == NULL) {
    ofs[fileno++] = input_fname;
  } else {
    ofs[fileno++] = input_fname2;
    printf("NOTE: Merging multiple files frame-by-frame\n");
  }

  char default_range_fname[] = "rangeFile.dat";
  if ( range_fname == NULL ) {
    // usage("No range file name supplied");
    cout<<"No range file name given."<<endl;
    cout<<"Using default range file name "<<"rangeFile.dat"<<endl;
    range_fname = default_range_fname;
  }

  // Create objects
  
    lr_rng = new Range(lr_str,0,MAX_LABEL_VAL);
    if (!lr_rng->full()) { //or better yet if lr_str != NULL
      // only bother to open this if the label range isn't full.
      if (label_fname!=NULL) {
	gotLabelFile = true;
	if(labpos < 0) {
	  printf("Using the last discrete feature of the label file as the label\n");
	}
	lofs[0] = label_fname;
	lnfs[0] = 0; lnis[0] = 1;
	if(lswap == true) liswps[0] = true;
	if (lofs[0] != NULL && lnfs[0] == 0 && lnis[0] == 0)
	  error("ERROR: command line must specify one of lnf and lni not zero");
	unsigned lifmts[1];
	if (strcmp(fmts[0],"htk") == 0)
	  lifmts[0] = HTK;
	else if (strcmp(fmts[0],"binary") == 0)
	  lifmts[0] = RAWBIN;
	else if (strcmp(fmts[0],"ascii") == 0)
	  lifmts[0] = RAWASC;
	else if (strcmp(fmts[0],"pfile") == 0)
	  lifmts[0] = PFILE;
	else
	  error("ERROR: Unknown observation file format type: '%s'\n",lfmts[0]);
	
	globalLabelMatrix.openFiles(1,
				    (const char**)&lofs,
				    (const char**)&lfrs,
				    (const char**)&lirs,
				    (unsigned*)&lnfs,
				    (unsigned*)&lnis,
				    (unsigned*)&lifmts,
				    (bool*)&liswps);

      unsigned numLabelFeatures = globalLabelMatrix.numFeatures();
      unsigned numLabelContinuous = globalLabelMatrix.numContinuous();
      if(labpos != -1 && 
	 (labpos < (int)numLabelContinuous || labpos >= (int)numLabelFeatures) )
      error("labpos (%d) out of range (%d - %d): must be within the range of discrete observations of the label file\n",labpos,numLabelContinuous,numLabelFeatures-1);
      }
      else {
	gotLabelFile = false;
	if(labpos < 0) { //No label position has been given
	  printf("Using the last discrete feature of the last input file as the label\n");
	}
      }
    }

    ////////////////////////////////////////////
    // check for valid argument values.
    int nfiles = 0;
    unsigned ifmts[MAX_NUM_OBS_FILES];
    for (int i=0;i<MAX_NUM_OBS_FILES;i++) {
      if (ofs[i] != NULL && nfs[i] == 0 && nis[i] == 0)
	error("ERROR: command line must specify one of nf%d and ni%d not zero",i+1,i+1);
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
				      (bool*)&iswps);
 
    unsigned numFeatures = globalObservationMatrix.numFeatures();
    unsigned numContinuous = globalObservationMatrix.numContinuous();
    if(!gotLabelFile && labpos != -1 && 
       (labpos < (int)numContinuous || labpos >= (int)numFeatures) )
      error("labpos (%d) out of range (%d - %d): must be within the range of discrete obsevations\n",labpos,numContinuous,numFeatures-1);

    cout<<"nfiles ="<<nfiles<<endl;
    cout<<"numSegs = "<<globalObservationMatrix.numSegments()<<endl;
    cout<<"numFeatures = "<<numFeatures<<endl;
    
  sr_rng = new Range(sr_str,0,globalObservationMatrix.numSegments());  
  fr_rng = new Range(fr_str,0,0);
 
 // If an output pfile name is not supplied, we just
  // compute the statistics.
  FILE *out_fp=NULL;
  if (output_fname==NULL || !strcmp(output_fname,"-"))
    out_fp = stdout; // no output pfile desired.
  else {
    if ((out_fp = fopen(output_fname, "w")) == NULL) {
      error("Couldn't open output file for writing.");
    }
  }     

  cout<<"Parameters setup. Starting...\n"; 
  // Do the work.
  discreteMI(out_fp, 
	     range_fname,
	     *sr_rng, 
	     *fr_rng, 
	     *lr_rng, 
	     print_frameid, 
	     quiet, 
	     binary,  
	     labpos);

  // Clean up and exit.
  if (out_fp && fclose(out_fp))
    error("Couldn't close output file.");

  delete sr_rng;
  delete fr_rng;
  delete lr_rng;

  return EXIT_SUCCESS;
}





