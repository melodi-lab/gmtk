/*
    $Header$

    This program 

    1) reads a set of observation files (input files in PFILE, HTK, ascii or binary format), 

    2) (optionally) performs a number of possible transformations:
    
        - sentence-wide transformations, which are implemented in the
          GMTK Observation library and include selecting a subset of
          the features, frames, or sentences, downsampling,
          upsampling, sentence-wide normalization, filtering,
          combination of the input files
  
        - file-wide transformations, such as normalization (both
          linear and Gaussian), Karhunen-Loeve (PCA), addition of
          silence frames, calculation of various statistics,

    3) writes the transformed or not output in PFILE, HTK, ascii, or
       binary.  In the case of ascii and binary outputs there is the
       additional option of a "flat output" (a single file with all of
       the observation file optionaly index by the sentence and frame
       numbers).
  
    Many features of this program are based on the pfile tool-set by Jeff Bilmes <bilmes@ee.washington.edu>


    
    Added the option to preprocess the input list of file names using CPP -- 29aug2003
    Endian is found dynamically and swapping defaults are set accordingly -- 01dec2003 

    Karim Filali <karim@cs.washington.edu>

*/

#include <stdlib.h>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <values.h>
#include <cmath>
#include <cassert>
#include "pfile.h"
//#include "parse_subset.h"
#include "error.h"
#include "arguments.h"
#include "GMTK_WordOrganization.h"

#include "GMTK_ObsPrint.h"
#include "GMTK_ObsKLT.h"
#include "GMTK_ObsStats.h"
#include "GMTK_ObsNorm.h"
#include "GMTK_ObsGaussianNorm.h"
#include "GMTK_ObsAddSil.h"

#ifdef DEBUG
#define DBGFPRINTF(_x_) fprintf _x_
#else
#define DBGFPRINTF(_x_)
#endif

#define MAXHISTBINS 1000

#define MIN(a,b) ((a)<(b)?(a):(b))

ObservationMatrix globalObservationMatrix;

void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;

char * output_fname = NULL; // Output pfile name.
char * outputList = NULL;
FILE *outputListFp=NULL;
char * outputNameSeparatorStr="_";



int HTK_Param_Kind=USER;
int HTK_Sample_Period=1;


/**

Prints an HTK header given the dimensionality of the feature vector and the number of samples

*/

void printHTKHeader(FILE* ofp, bool oswap, int numInts, int numFloats, int numSamples, short parameterKind=USER, Int32 samplePeriod=1) {

  DBGFPRINTF((stderr,"obsPrint: Printing HTK header.\n"));
  // structure of HTK header
  //  Int32 numSamples;
  //  Int32 samplePeriod=1;
  //  short parameterKind;
  short sampleSize;
  if(numFloats > 0) {
    sampleSize=numFloats*sizeof(float);
    //    parameterKind=USER;  // enum in GMTK_Stream.h
    //    parameterKind=paramKind;  // enum in GMTK_Stream.h
  }
  else {
    sampleSize=numInts*sizeof(short);
    parameterKind=DISCRETE;  // enum in GMTK_Stream.h
  }

  if (oswap) {
    numSamples = swapb_i32_i32(numSamples);
    samplePeriod = swapb_i32_i32(samplePeriod);
    sampleSize = swapb_short_short(sampleSize);
    parameterKind = swapb_short_short(parameterKind);
  }

  if (fwrite(&numSamples,sizeof(Int32),1,ofp) != 1) {
    error("Cannot write HTK number of samples\n");
  }
  if (fwrite((short *)&samplePeriod,sizeof(Int32),1,ofp) != 1) {
    error("Cannot write HTK sample period\n");
  }

  if (fwrite((short *)&sampleSize,sizeof(short),1,ofp) != 1) {
    error("Cannot write HTK sample size\n");
  }

  if (fwrite(&parameterKind,sizeof(short),1,ofp) != 1) {
    error("Cannot write HTK parm kind\n");
  }

  DBGFPRINTF((stderr,"obsPrint: Finished printing HTK header.\n"));

}

void printSegment(unsigned sent_no, FILE* out_fp, float* cont_buf, unsigned num_continuous, UInt32* disc_buf, unsigned num_discrete, unsigned num_frames, const bool dontPrintFrameID,const bool quiet,unsigned ofmt,int debug_level,bool oswap, OutFtrLabStream_PFile* out_stream) {


    if(ofmt==RAWASC || ofmt==RAWBIN || ofmt==HTK) {
      char* current_output_fname = new char[strlen(output_fname)+strlen(outputNameSeparatorStr)+50];
      sprintf(current_output_fname,"%s%s%d",output_fname,outputNameSeparatorStr,sent_no);
      if ((out_fp = fopen(current_output_fname, "w")) == NULL) {
	error("Couldn't open output file (%s) for writing.",current_output_fname);
      }

      if(outputListFp != NULL) {
	fprintf(outputListFp,"%s\n",current_output_fname);
      }

      delete []  current_output_fname;
    }
    
    if(ofmt==HTK) {
      DBGFPRINTF((stderr,"obsPrint: Calling printHTKHeader with numSamples=%d.\n",num_frames));
      printHTKHeader(out_fp,oswap,num_discrete,num_continuous,num_frames);
    }

    float* cont_buf_p = cont_buf;
    UInt32* disc_buf_p = disc_buf;
    //////////  Print the frames ////////////////////////////////////
    for (unsigned frame_no=0; frame_no < num_frames ; ++frame_no) {
      bool ns = false;
      if (!dontPrintFrameID) {
	if (ofmt==FLATBIN || ofmt==RAWBIN) {
	  copy_swap_func_ptr(1,(int*)&sent_no,(int*)&sent_no);
	  copy_swap_func_ptr(1,(int*)&frame_no,(int*)&frame_no);
	  fwrite(&sent_no,sizeof(sent_no),1,out_fp);
	  fwrite(&frame_no,sizeof(frame_no),1,out_fp);
	} else if(ofmt==FLATASC || ofmt==RAWASC ){
	  fprintf(out_fp,"%d %u",sent_no,frame_no);
	  ns = true;
	}
      }
      
      /// Print continuous part of frame /////////////////////////////
      for (unsigned frit=0;frit< num_continuous; ++frit) {
	if (ofmt==FLATBIN || ofmt==RAWBIN || ofmt==HTK) {
	  DBGFPRINTF((stderr,"obsPrint: Printing HTK float %f.\n",cont_buf_p[frit]));
	  copy_swap_func_ptr(1,(int*)&cont_buf_p[frit],(int*)&cont_buf_p[frit]);
	  fwrite(&cont_buf_p[frit], sizeof(cont_buf_p[frit]),  1,out_fp);
	} 
	else if(ofmt==FLATASC || ofmt==RAWASC){
	  if (ns) fprintf(out_fp," ");
	  fprintf(out_fp,"%f",cont_buf_p[frit]);
	  ns = true;
	}
      }
      ///////////////////////////////////////////////////////////////

      /// Print discrete part of the frame ///////////////////////////
      for (unsigned lrit=0;lrit<num_discrete; ++lrit) {
	if (ofmt==FLATBIN || ofmt==RAWBIN || (ofmt==HTK && num_continuous>0) ) {
	  copy_swap_func_ptr(1,(int*)&disc_buf_p[lrit],(int*)&disc_buf_p[lrit]);
	  fwrite(&disc_buf_p[lrit],  sizeof(disc_buf_p[lrit]), 1,out_fp);
	} 
	else if(ofmt==HTK && num_continuous==0) { // in the HTK format we
  // cannot mix floats with discrete data; that's why if there is at
  // least one float component everyting is written out as a float.
	  short short_lab_buf_p=(short)disc_buf_p[lrit];
	  DBGFPRINTF((stderr,"obsPrint: Printing HTK short %d.\n",short_lab_buf_p));
	   if (oswap) {
	     short_lab_buf_p = swapb_short_short(short_lab_buf_p);
	   }
	  fwrite(&short_lab_buf_p,  sizeof(short_lab_buf_p), 1,out_fp);
	}
	else if(ofmt==FLATASC || ofmt==RAWASC) {
	  if (ns) fprintf(out_fp," ");	    
	  fprintf(out_fp,"%d",disc_buf_p[lrit]);
	  ns = true;
	}
      } // end of for (unsigned lrit=0;lrit<num_discrete; ++lrit)
      ////////////////////////////////////////////////////////////////

      if (ofmt==FLATASC || ofmt==RAWASC)
	fprintf(out_fp,"\n");

      cont_buf_p += num_continuous;
      disc_buf_p += num_discrete;
    }  // end of for (unsigned frame_no=0; frame_no < num_frames ; ++frame_no)
    ////////////////////////////////////////////////////////////////////////////

    if(ofmt==PFILE) {
      out_stream->write_ftrslabs(num_frames, cont_buf, disc_buf);
      out_stream->doneseg((SegID) sent_no);
    }
    
    if(ofmt==RAWASC || ofmt==RAWBIN || ofmt==HTK) {
      if (fclose(out_fp)) error("Couldn't close output file.");
    }
}

void obsPrint(FILE* out_fp,Range& srrng,const char * pr_str,const bool dontPrintFrameID,const bool quiet,unsigned ofmt,int debug_level,bool oswap) {

  // Feature and label buffers are dynamically grown as needed.
  size_t buf_size = 300;      // Start with storage for 300 frames.
  const size_t n_labs = globalObservationMatrix.numDiscrete();
  const size_t n_ftrs = globalObservationMatrix.numContinuous();

  float *ftr_buf_p;
  UInt32* lab_buf_p;
  
  float *oftr_buf=NULL,*oftr_buf_p=NULL;
  UInt32* olab_buf=NULL,* olab_buf_p=NULL;
  OutFtrLabStream_PFile* out_stream=NULL;
  
  if(ofmt==PFILE) {
    out_stream = new OutFtrLabStream_PFile(debug_level,"",out_fp,n_ftrs,n_labs,1,oswap);
    oftr_buf = new float[buf_size *  n_ftrs];
    olab_buf = new UInt32[buf_size * n_labs];
  }
  
  // Go through input pfile to get the initial statistics,
  // i.e., max, min, mean, std, etc.
  for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {
    DBGFPRINTF((stderr,"obsPrint: Loading segment no %d\n",*srit));
    globalObservationMatrix.loadSegment((const unsigned)(*srit));
    const size_t n_frames = globalObservationMatrix.numFrames();
    
    // Increase size of buffers if needed.
    if (ofmt==PFILE && n_frames > buf_size) {
      // Free old buffers.
      delete [] oftr_buf;
      delete [] olab_buf;
      // Make twice as big to cut down on future reallocs.
      buf_size = n_frames * 2;
      // Allocate new larger buffers.
      oftr_buf = new float[buf_size *  n_ftrs];
      olab_buf = new UInt32[buf_size * n_labs];
    }
    
    Range prrng(pr_str,0,n_frames);
    
    if (!quiet) {
      if (*srit % 100 == 0)
	printf("Processing sentence %d\n",*srit);
    }
    
    if(ofmt==PFILE) {
      oftr_buf_p = oftr_buf;
      olab_buf_p = olab_buf;
    }
    else if(ofmt==RAWASC || ofmt==RAWBIN || ofmt==HTK) {
      char* current_output_fname = new char[strlen(output_fname)+strlen(outputNameSeparatorStr)+50];
      sprintf(current_output_fname,"%s%s%d",output_fname,outputNameSeparatorStr,*srit);
      if ((out_fp = fopen(current_output_fname, "w")) == NULL) {
	error("Couldn't open output file (%s) for writing.",current_output_fname);
      }

      if(outputListFp != NULL) {
	fprintf(outputListFp,"%s\n",current_output_fname);
      }

      delete []  current_output_fname;
    }
    
    if(ofmt==HTK) {
      DBGFPRINTF((stderr,"obsPrint: Calling printHTKHeader with numSamples=%d.\n",prrng.length()));
      printHTKHeader(out_fp,oswap,n_labs,n_ftrs,prrng.length(),(short)HTK_Param_Kind, HTK_Sample_Period);
    }

    for (Range::iterator prit=prrng.begin(); !prit.at_end() ; ++prit) {
      bool ns = false;
      ftr_buf_p = globalObservationMatrix.floatVecAtFrame((*prit));
      lab_buf_p = globalObservationMatrix.unsignedAtFrame((*prit));
      if (!dontPrintFrameID) {
	//	if (ofmt==FLATBIN || ofmt==RAWBIN) {
	if (ofmt==FLATBIN) {
	  size_t sent_no = *srit;
	  size_t frame_no = *prit;
	  
	  copy_swap_func_ptr(1,(int*)&sent_no,(int*)&sent_no);
	  copy_swap_func_ptr(1,(int*)&frame_no,(int*)&frame_no);

	  fwrite(&sent_no,sizeof(sent_no),1,out_fp);
	  fwrite(&frame_no,sizeof(frame_no),1,out_fp);
	} else if(ofmt==FLATASC || ofmt==RAWASC ){
	  fprintf(out_fp,"%d %u",*srit,*prit);
	  ns = true;
	}
      }
      
      for (unsigned frit=0;frit<globalObservationMatrix.numContinuous(); ++frit) {
	if (ofmt==FLATBIN || ofmt==RAWBIN || ofmt==HTK) {
	  DBGFPRINTF((stderr,"obsPrint: Printing HTK float %f.\n",ftr_buf_p[frit]));
	  copy_swap_func_ptr(1,(int*)&ftr_buf_p[frit],(int*)&ftr_buf_p[frit]);
	  fwrite(&ftr_buf_p[frit], sizeof(ftr_buf_p[frit]),  1,out_fp);
	} 
	else if(ofmt==PFILE) {
	  *oftr_buf_p++ = ftr_buf_p[frit];
	}
	else if(ofmt==FLATASC || ofmt==RAWASC){
	  if (ns) fprintf(out_fp," ");
	  fprintf(out_fp,"%f",ftr_buf_p[frit]);
	  ns = true;
	}
      }
      for (unsigned lrit=0;lrit<globalObservationMatrix.numDiscrete(); ++lrit) {
	if (ofmt==FLATBIN || ofmt==RAWBIN || (ofmt==HTK && n_ftrs>0) ) {
	  copy_swap_func_ptr(1,(int*)&lab_buf_p[lrit],(int*)&lab_buf_p[lrit]);
	  fwrite(&lab_buf_p[lrit],  sizeof(lab_buf_p[lrit]), 1,out_fp);
	} 
	else if(ofmt==HTK && n_ftrs==0) { // in the HTK format we
  // cannot mix floats with discrete data; that's why if there is at
  // least one float component everything is written out as a float.
	  short short_lab_buf=(short)lab_buf_p[lrit];
	  DBGFPRINTF((stderr,"obsPrint: Printing HTK short %d.\n",short_lab_buf));
	   if (oswap) {
	     short_lab_buf = swapb_short_short(short_lab_buf);
	   }
	  fwrite(&short_lab_buf,  sizeof(short_lab_buf), 1,out_fp);
	}
	else if(ofmt==PFILE) {
	  *olab_buf_p++ = lab_buf_p[lrit];
	}
	else if(ofmt==FLATASC || ofmt==RAWASC) {
	  if (ns) fprintf(out_fp," ");	    
	  fprintf(out_fp,"%d",lab_buf_p[lrit]);
	  ns = true;
	}
      }
      if (ofmt==FLATASC || ofmt==RAWASC)
	fprintf(out_fp,"\n");
    }
    if(ofmt==PFILE) {
      out_stream->write_ftrslabs(prrng.length(), oftr_buf, olab_buf);
      out_stream->doneseg((SegID) *srit);
    }
    if(ofmt==RAWASC || ofmt==RAWBIN || ofmt==HTK) {
      if (fclose(out_fp)) error("Couldn't close output file.");
    }
  }
  
  if(ofmt==PFILE) {
    delete [] oftr_buf;
    delete [] olab_buf;
    delete out_stream;
  }

  DBGFPRINTF((stderr,"obsPrint: End.\n"));

}




#define MAX_OBJECTS 10

char *input_fname[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};  // Input file name.
char * ifmtStr[MAX_OBJECTS]={"pfile","pfile","pfile","pfile","pfile","pfile","pfile","pfile","pfile","pfile"};
unsigned ifmt[MAX_OBJECTS];

char * ofmtStr="flatasc";
unsigned ofmt;

unsigned int nis[MAX_OBJECTS];
unsigned int nfs[MAX_OBJECTS];

char  *sr_str               = 0;   // sentence range string
Range *sr_rng;
char  *fr_str[MAX_OBJECTS]  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};   // feature range string    
char  *lr_str[MAX_OBJECTS]  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};   // label range string  
char  *spr_str[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};   // per stream per sentence range string 
char  *pr_str               = 0;   // per-sentence range string

char* actionIfDiffNumFramesStr[MAX_OBJECTS]={"er","er","er","er","er","er","er","er","er","er"};   // 
unsigned actionIfDiffNumFrames[MAX_OBJECTS]={ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR};   // 


char*    actionIfDiffNumSentsStr[MAX_OBJECTS] = {"te","te","te","te","te","te","te","te","te","te"}; 
unsigned actionIfDiffNumSents[MAX_OBJECTS]    = {TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END};   // 

char* perStreamTransforms[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};   // 
char* postTransforms                   = NULL;

int  debug_level = 0;
bool dontPrintFrameID = false;
bool quiet = false;
#ifdef INTV_WORDS_BIGENDIAN
bool iswap[MAX_OBJECTS]={true,true,true,true,true,true,true,true,true,true};
bool oswap = true;
#else
bool iswap[MAX_OBJECTS]= {false,false,false,false,false,false,false,false,false,false};
bool oswap             = false;
#endif 

bool     cppIfAscii        = true;
char*    cppCommandOptions = NULL;

unsigned startSkip = 0;
unsigned endSkip   = 0;


char*    ftrcomboStr = "none";
unsigned ftrcombo    = FTROP_NONE;


bool     Normalize = false;
double   Norm_Mean = 0;
double   Norm_Std  = 1;
UInt32   Norm_Segment_Group_Len      = INT_MAX;  // treat all sentences as one big segment 
char*    Norm_Segment_Group_Len_File = NULL;

bool     Perform_KLT               = false;
char*    KLT_Output_Stat_File_Name = NULL;
char*    KLT_Input_Stat_File_Name  = NULL;
bool     KLT_Ascii_Stat_Files      = false;
char*    KLT_Output_Ftr_Range      = NULL;
bool     KLT_Unity_Variance        = true;

bool     Get_Stats           = false;
unsigned Num_Hist_Bins = 0;

bool     Gaussian_Norm                    = false;
float    Gaussian_Num_Stds                = 5;
bool     Gaussian_Uniform                 = false;
char*    Gauss_Norm_Output_Stat_File_Name = NULL;
char*    Gauss_Norm_Input_Stat_File_Name  = NULL;

bool     Add_Sil                = false;
unsigned Add_Sil_Num_Beg_Frames = 0;
unsigned Add_Sil_Num_End_Frames = 0;
char*    Add_Sil_Beg_Rng_Str    = NULL; 
char*    Add_Sil_End_Rng_Str    = NULL; 
double   Add_Sil_MMF            = 1.0;
double   Add_Sil_MAF            = 0.0;
double   Add_Sil_SMF            = 1.0;
double   Add_Sil_SAF            = 0.0;

char* Usage_Str = NULL;
bool help=false;

Arg Arg::Args[] = {
  Arg("i",      Arg::Req, input_fname,"input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("ifmt",      Arg::Opt,ifmtStr ,"format of input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("o",      Arg::Opt, output_fname,"output file"),
  Arg("ofmt",      Arg::Opt, ofmtStr,"format of output file"),
  Arg("nf",   Arg::Opt, nfs,"number of floats in input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("ni",   Arg::Opt, nis,"number of ints (labels) in input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("q",  Arg::Tog, quiet,"quiet mode"),
  Arg("sr",   Arg::Opt, sr_str,"sentence range"),
  Arg("fr",   Arg::Opt, fr_str,"feature range",Arg::ARRAY,MAX_OBJECTS),
  Arg("spr",  Arg::Opt, spr_str,"per stream per-sentence range",Arg::ARRAY,MAX_OBJECTS),
  //  Arg("pr",   Arg::Opt, pr_str,"per-sentence range"),
  Arg("lr",   Arg::Opt, lr_str,"label range",Arg::ARRAY,MAX_OBJECTS),
  Arg("startskip",   Arg::Opt, startSkip,"start skip"),
  Arg("endskip",   Arg::Opt, endSkip,"end skip"),
  Arg("fdiffact",  Arg::Opt,actionIfDiffNumFramesStr ,"Action if different number of frames in streams: error (er), repeat last frame (rl), first frame (rf), segmentally expand (se), truncate from start (ts), truncate from end (te)",Arg::ARRAY,MAX_OBJECTS),
  Arg("sdiffact",  Arg::Opt,actionIfDiffNumSentsStr ,"Action if different number of sentences in streams: error (er), truncate from end (te), repeat last sent (rl), and wrap around (wa).",Arg::ARRAY,MAX_OBJECTS),
  Arg("trans",  Arg::Opt,perStreamTransforms ,"per stream transformations string",Arg::ARRAY,MAX_OBJECTS),
  Arg("posttrans",  Arg::Opt,postTransforms ,"Final global transformations string"),
  Arg("comb",      Arg::Opt, ftrcomboStr,"Combine float features (none: no combination, add, sub, mul,div"),
  Arg("iswap",Arg::Opt, iswap,"do byte swapping on the input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("oswap",Arg::Opt, oswap,"do byte swapping on the output file"),
  Arg("ns",    Arg::Opt, dontPrintFrameID,"Don't print the frame IDs (i.e., sent and frame #)"),
  Arg("sep",      Arg::Opt, outputNameSeparatorStr,"String to use as separator when outputting raw ascii or binary files (one sentence per file)."),
  Arg("olist",      Arg::Opt, outputList,"output list-of-files name.  Only meaningful if used with the RAW or HTK formats."),
  Arg("cppifascii",Arg::Opt, cppIfAscii,"Pre-process ASCII files using CPP"),
  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),
  Arg("norm",      Arg::Tog, Normalize, "Normalize the observation file"),
  Arg("normMean",  Arg::Opt, Norm_Mean, "NORM: Mean of the resulting output file",Arg::SINGLE,0,true),
  Arg("normStd",   Arg::Opt, Norm_Std,  "NORM: Std of the resulting output file",Arg::SINGLE,0,true),
  Arg("normSl",    Arg::Opt, Norm_Segment_Group_Len,      "NORM: Segment group length",Arg::SINGLE,0,true),
  Arg("normSlf",   Arg::Opt, Norm_Segment_Group_Len_File, "NORM: Ascii file with segment group lengths",Arg::SINGLE,0,true),
  Arg("gauss",           Arg::Tog, Gaussian_Norm, "Normalize the observation file to be Gaussian distributed with zzero mean and unit variance. The features are scaled to be within +/- gaussNumStds standard deviations"),
  Arg("gaussNumStds",    Arg::Opt, Gaussian_Num_Stds, "GAUSS: Number of Gaussian standard deviations",Arg::SINGLE,0,true),
  Arg("gaussUniform",    Arg::Opt, Gaussian_Uniform,  "GAUSS: Output is uniform[0,1] distributed rather than Gaussian",Arg::SINGLE,0,true),
  Arg("gaussOutputStat", Arg::Opt, Gauss_Norm_Output_Stat_File_Name, "GAUSS: Output statistics file",Arg::SINGLE,0,true),
  Arg("gaussInputStat",  Arg::Opt, Gauss_Norm_Input_Stat_File_Name,  "GAUSS: Input statistics file",Arg::SINGLE,0,true),
  Arg("klt",               Arg::Tog, Perform_KLT,"Perform a KLT transform"),
  Arg("kltUnityVar",       Arg::Opt, KLT_Unity_Variance,"KLT: multiply features by inverse eigenvalues (to have unity variance)",Arg::SINGLE,0,true),
  Arg("kltOutputStat",     Arg::Opt,KLT_Output_Stat_File_Name ,"KLT: output stats (covariance, mean, eigenvectors, eigenvalues) matrices binary double precision ('-' for stdout)",Arg::SINGLE,0,true),
  Arg("kltInputStat",      Arg::Opt,KLT_Input_Stat_File_Name ,"KLT: input stats (covariance, mean, eigenvectors, eigenvalues) matrices (i.e., do not compute them)",Arg::SINGLE,0,true),
  Arg("kltAsciiStat",      Arg::Opt,KLT_Ascii_Stat_Files ,"KLT: stat matrices written/read in ascii rather than binary doubles",Arg::SINGLE,0,true),
  Arg("kltOutputFtrRange", Arg::Opt, KLT_Output_Ftr_Range,"KLT: output feature range",Arg::SINGLE,0,true),
  Arg("stats",           Arg::Tog, Get_Stats,"Output statistics of the form:\nfeatnum mean std max @sent# @frame# min @sent# @frame# max/stds min/stds [histogram]"),
  Arg("histBins",   Arg::Opt, Num_Hist_Bins,"STATS/GAUSS: number of histogram bins"),
 Arg("addsil",               Arg::Tog, Add_Sil,"Add silence frames at the begining and end each sentence"), 
 Arg("addsilNumBeg",               Arg::Opt, Add_Sil_Num_Beg_Frames,"Number of new beginning silence frames",Arg::SINGLE,0,true), 
 Arg("addsilPrb",               Arg::Opt, Add_Sil_Beg_Rng_Str,"Per-sentence range to compute beginning silence",Arg::SINGLE,0,true), 
 Arg("addsilNumEnd",               Arg::Opt, Add_Sil_Num_End_Frames,"Number of new ending silence frames",Arg::SINGLE,0,true), 
 Arg("addsilPre",               Arg::Opt, Add_Sil_End_Rng_Str,"Per-sentence range to compute ending silence",Arg::SINGLE,0,true), 
 Arg("addsilMMF",               Arg::Opt, Add_Sil_MMF,"Mean multiplicative factor",Arg::SINGLE,0,true), 
 Arg("addsilMAF",               Arg::Opt, Add_Sil_MAF,"Mean additive factor",Arg::SINGLE,0,true),
 Arg("addsilSMF",               Arg::Opt, Add_Sil_SMF,"Standard deviation multiplicative factor",Arg::SINGLE,0,true),  
 Arg("addsilSAF",               Arg::Opt, Add_Sil_SAF,"Standard deviation additive factor",Arg::SINGLE,0,true),  
  Arg("htkKind",  Arg::Opt, HTK_Param_Kind,"Kind of output HTK parameters"),
  Arg("htkSamplePeriod",  Arg::Opt, HTK_Sample_Period,"Output HTK Sample Period"),
 Arg("debug",  Arg::Opt, debug_level,"Number giving level of debugging output to produce 0=none"),
  Arg("usage",  Arg::Opt, Usage_Str, "Print usage information about one of the follwoing topics: {norm, gauss, klt, addsil}"),
  Arg("help",   Arg::Tog, help,  "Print this message"),
  // The argumentless argument marks the end of the above list.
  Arg()
};


int main(int argc, const char *argv[]) {
  
  int numFiles=0;

  // Figure out the Endian of the machine this is running on and set the swap defaults accordingly
  bool doWeSwap;

  ByteEndian byteEndian = getWordOrganization();
  switch(byteEndian) {
  case BYTE_BIG_ENDIAN:
    doWeSwap=false;
    break;
  case BYTE_LITTLE_ENDIAN:
    doWeSwap=true;
    break;
  default:
    // We weren't able to figure the Endian out.  Leave the swap defaults as they are.
#ifdef INTV_WORDS_BIGENDIAN
    doWeSwap=true;
#else
    doWeSwap=false;
#endif
  }

  oswap=doWeSwap;
  for(int i=0; i<MAX_OBJECTS; ++i) {
    iswap[i]=doWeSwap;
  }
  ///////////////////////////////////////////

  bool parse_was_ok = Arg::parse(argc,(char**)argv);

  if(help) {
    Arg::usage();
    exit(0);
  }
  if(Usage_Str!=NULL) {
    if(strcmp(Usage_Str,"norm")==0)  Arg::usage("norm");
    else if(strcmp(Usage_Str,"gauss")==0)  Arg::usage("gauss");
    else if(strcmp(Usage_Str,"klt")==0)  Arg::usage("klt");
    else if(strcmp(Usage_Str,"stat")==0)  Arg::usage("stat");
    else if(strcmp(Usage_Str,"addsil")==0)  Arg::usage("addsil");
    else {
      error("No usage information on \"%s\"",Usage_Str);
    }
    exit(0);
  }

  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }


  // oswap might been assigned a new value on the command line
  if(oswap) {
    copy_swap_func_ptr=&swapb_vi32_vi32;
  }
  else {
    copy_swap_func_ptr=&copy_vi32_vi32;
  }

  //////////////////////////////////////////////////////////////////////
  // Check all necessary arguments provided before creating objects.
  //////////////////////////////////////////////////////////////////////

 for (int i=0;i<MAX_OBJECTS;i++) {
    numFiles += (input_fname[i] != NULL);
    if (strcmp(ifmtStr[i],"htk") == 0)
      ifmt[i] = HTK;
    else if (strcmp(ifmtStr[i],"binary") == 0)
      ifmt[i] = RAWBIN;
    else if (strcmp(ifmtStr[i],"ascii") == 0)
      ifmt[i] = RAWASC;
    else if (strcmp(ifmtStr[i],"pfile") == 0)
      ifmt[i] = PFILE;
    //    else if (strcmp(ifmtStr[i],"flatbin") == 0)
    //  ifmt[i]=FLATBIN;
    //else if (strcmp(ifmtStr[i],"flatasc") == 0)
    //  ifmt[i]=FLATASC;
    else
      error("ERROR: Unknown observation file format type: '%s'\n",ifmtStr[i]);
  }

 for(int i=0; i < numFiles; ++i)
   if(output_fname!=NULL && strcmp(input_fname[i],output_fname)==0) {
     error("Input and output filenames cannot be the same.");
   }
 
 if (strcmp(ofmtStr,"htk") == 0)
   ofmt = HTK;
 else if (strcmp(ofmtStr,"binary") == 0)
   ofmt = RAWBIN;
 else if (strcmp(ofmtStr,"ascii") == 0)
   ofmt = RAWASC;
 else if (strcmp(ofmtStr,"pfile") == 0)
   ofmt = PFILE;
 else if (strcmp(ofmtStr,"flatbin") == 0)
   ofmt=FLATBIN;
 else if (strcmp(ofmtStr,"flatasc") == 0)
   ofmt=FLATASC;
 else
   error("ERROR: Unknown observation file format type: '%s'\n",ofmtStr);

if (strcmp(ftrcomboStr,"none") == 0)
   ftrcombo = FTROP_NONE;
 else if (strcmp(ftrcomboStr,"add") == 0)
   ftrcombo = FTROP_ADD;
 else if (strcmp(ftrcomboStr,"sub") == 0)
   ftrcombo = FTROP_SUB;
 else if (strcmp(ftrcomboStr,"mul") == 0)
   ftrcombo = FTROP_MUL;
 else if (strcmp(ftrcomboStr,"div") == 0)
   ftrcombo = FTROP_DIV;
 else
   error("ERROR: Unknown feature combination type: '%s'\n",ftrcomboStr);

 for(int i=0; i < MAX_OBJECTS; ++i) {
   if(input_fname[i]!=NULL) {
     if (strcmp(actionIfDiffNumFramesStr[i],"er") == 0)
       actionIfDiffNumFrames[i] = ERROR;
     else if (strcmp(actionIfDiffNumFramesStr[i],"rl") == 0)
       actionIfDiffNumFrames[i] = REPEAT_LAST;
     else if (strcmp(actionIfDiffNumFramesStr[i],"rf") == 0)
       actionIfDiffNumFrames[i] = REPEAT_FIRST;
     else if (strcmp(actionIfDiffNumFramesStr[i],"se") == 0)
       actionIfDiffNumFrames[i] = EXPAND_SEGMENTALLY;
     else if (strcmp(actionIfDiffNumFramesStr[i],"ts") == 0)
       actionIfDiffNumFrames[i] = TRUNCATE_FROM_START;
     else if (strcmp(actionIfDiffNumFramesStr[i],"te") == 0)
       actionIfDiffNumFrames[i] = TRUNCATE_FROM_END;
     else
       error("ERROR: Unknown action when diff num of frames: '%s'\n",actionIfDiffNumFramesStr[i]);
   }
 }

for(int i=0; i < MAX_OBJECTS; ++i) {
   if(input_fname[i]!=NULL) {
     if (strcmp(actionIfDiffNumSentsStr[i],"er") == 0)
       actionIfDiffNumSents[i] = ERROR;
     else if (strcmp(actionIfDiffNumSentsStr[i],"rl") == 0)
       actionIfDiffNumSents[i] = REPEAT_LAST;
     else if (strcmp(actionIfDiffNumSentsStr[i],"wa") == 0)
       actionIfDiffNumSents[i] = WRAP_AROUND;
     else if (strcmp(actionIfDiffNumSentsStr[i],"te") == 0)
       actionIfDiffNumSents[i] = TRUNCATE_FROM_END;
     else
       error("ERROR: Unknown action when diff num of sentences: '%s'\n",actionIfDiffNumSentsStr[i]);
   }
 }


 FILE *out_fp=NULL;
 if (output_fname==0 || !strcmp(output_fname,"-")) {
   if(ofmt==RAWASC || ofmt==RAWBIN || ofmt==HTK) {
     error("Need to specify output filename when output type is RAWBIN, RAWASC, or HTK.\n");
   }
   out_fp = stdout;
 }
 else {
   if(ofmt != RAWASC && ofmt != RAWBIN && ofmt != HTK) {
     if ((out_fp = fopen(output_fname, "w")) == NULL) {
       error("Couldn't open output file for writing.\n");
     }
   }
 }


 if(outputList != NULL) {
   if ((outputListFp = fopen(outputList, "w")) == NULL) {
     error("Couldn't open output list (%s) for writing.\n",outputList);
   }
 }

 //////////////////////////////////////////////////////////////////////
 // Create objects.
 //////////////////////////////////////////////////////////////////////
 // If we have a pfile, we can extract the number if features from the file directly
 for(int i=0; i < MAX_OBJECTS; ++i) {
   if(input_fname[i]!=NULL) {
     if(ifmt[i]==PFILE) {
       FILE *in_fp = fopen(input_fname[i], "r");
       if (in_fp==NULL) error("Couldn't open input pfile for reading.");
       InFtrLabStream_PFile* in_streamp = new InFtrLabStream_PFile(debug_level,"",in_fp,1,iswap[i]);
       nis[i]=in_streamp->num_labs();
       nfs[i]=in_streamp->num_ftrs();
       if (fclose(in_fp)) error("Couldn't close input pfile.");
       delete in_streamp;
     }
     
     if(nis[i]==0 && nfs[i]==0) {
       error("The number of floats and the number of ints cannot be both zero.");
     }
   }
 }
 
 globalObservationMatrix.openFiles(numFiles,  // number of files.   For now we use only one
				   (const char**)&input_fname,
				   (const char**)&fr_str,
				   (const char**)&lr_str,
				   (unsigned*)&nfs,
				   (unsigned*)&nis,
				   (unsigned*)&ifmt,
				   (bool*)&iswap,
				   startSkip,  // startSkip
				   endSkip,  // endSkip  pr_rng takes care of these two  
				   cppIfAscii,
				   cppCommandOptions,
				   (const char**)&spr_str,
				   actionIfDiffNumFrames,
				   actionIfDiffNumSents,
				   perStreamTransforms,
				   postTransforms,
				   ftrcombo);   


     sr_rng = new Range(sr_str,0,globalObservationMatrix.numSegments());

     /////////////////////////////////////////////////////////////////////
     unsigned too_many_switches=0;
     if(Get_Stats)       ++too_many_switches; 
     if(Normalize)       ++too_many_switches; 
     if(Gaussian_Norm)   ++too_many_switches; 
     if(Perform_KLT)     ++too_many_switches; 
     if(Add_Sil)         ++too_many_switches; 

     if(too_many_switches > 1) {
       error("Cannot specify more than one of the following switches: -stats, -norm, -gauss, -klt, addsil");
     }

    //////////////////////////////////////////////////////////////////////
    // Do the work.
    //////////////////////////////////////////////////////////////////////

     if(Get_Stats) {
       Range* fr_rng=new Range(NULL,0,globalObservationMatrix.numContinuous());
       obsStats(out_fp, &globalObservationMatrix, *sr_rng, *fr_rng, NULL, Num_Hist_Bins, quiet);
       delete fr_rng;
     }
     else if(Normalize) {
       obsNorm(out_fp,&globalObservationMatrix,*sr_rng,Norm_Mean, Norm_Std, Norm_Segment_Group_Len_File, Norm_Segment_Group_Len, dontPrintFrameID,quiet,ofmt,debug_level,oswap);
     }
     else if(Gaussian_Norm) {
       Range* fr_rng=new Range(NULL,0,globalObservationMatrix.numContinuous());
       FILE* os_fp=NULL, * is_fp=NULL;
       if(Gauss_Norm_Output_Stat_File_Name != NULL) {
	 if((os_fp=fopen(Gauss_Norm_Output_Stat_File_Name,"w")) == NULL) {
	   error("Could not open output stat file, %s, for writing.",Gauss_Norm_Output_Stat_File_Name);
	 }
       }
       if(Gauss_Norm_Input_Stat_File_Name != NULL) {
	 if((is_fp=fopen(Gauss_Norm_Input_Stat_File_Name,"w")) == NULL) {
	   error("Could not open input stat file, %s, for writing.",Gauss_Norm_Input_Stat_File_Name);
	 }
       }
       gaussianNorm(out_fp,&globalObservationMatrix,is_fp,os_fp,*sr_rng,*fr_rng, pr_str, Num_Hist_Bins, Gaussian_Num_Stds, Gaussian_Uniform, dontPrintFrameID,quiet,ofmt,debug_level,oswap);
       delete fr_rng;
       if(Gauss_Norm_Output_Stat_File_Name != NULL) fclose(os_fp);
       if(Gauss_Norm_Input_Stat_File_Name  != NULL) fclose(is_fp);
     }
     else if(Perform_KLT) {
       Range* ofr_rng = new Range(KLT_Output_Ftr_Range,0,globalObservationMatrix.numContinuous());
       
       FILE* os_fp=NULL, * is_fp=NULL;
       if(KLT_Output_Stat_File_Name != NULL) {
	 if((os_fp=fopen(KLT_Output_Stat_File_Name,"w")) == NULL) {
	   error("Could not open output stat file, %s, for writing.",KLT_Output_Stat_File_Name);
	 }
       }
       if(KLT_Input_Stat_File_Name != NULL) {
	 if((is_fp=fopen(KLT_Input_Stat_File_Name,"w")) == NULL) {
	   error("Could not open input stat file, %s, for writing.",KLT_Input_Stat_File_Name);
	 }
       }
       
       obsKLT(out_fp,&globalObservationMatrix,is_fp,os_fp,*sr_rng,*ofr_rng,KLT_Unity_Variance, KLT_Ascii_Stat_Files, dontPrintFrameID,quiet,ofmt,debug_level,oswap);
       
       delete ofr_rng;
       if(KLT_Output_Stat_File_Name != NULL) fclose(os_fp);
       if(KLT_Input_Stat_File_Name != NULL) fclose(is_fp);
     }
     else if (Add_Sil) {
       addSil(out_fp,&globalObservationMatrix,*sr_rng,Add_Sil_Num_Beg_Frames,Add_Sil_Beg_Rng_Str, Add_Sil_Num_End_Frames,Add_Sil_End_Rng_Str, Add_Sil_MMF, Add_Sil_MAF, Add_Sil_SMF, Add_Sil_SAF, dontPrintFrameID,quiet,ofmt,debug_level,oswap);
     }
     else {
       obsPrint(out_fp,*sr_rng,pr_str,dontPrintFrameID,quiet,ofmt,debug_level,oswap);
     }
    //////////////////////////////////////////////////////////////////////
    // Clean up and exit.
    //////////////////////////////////////////////////////////////////////
    delete sr_rng;
    if(ofmt != RAWASC && ofmt != RAWBIN && ofmt != HTK) {
      if (fclose(out_fp)) error("Couldn't close output file.");
    }

    return EXIT_SUCCESS;
}
