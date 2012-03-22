

/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/*****************************                                     **********************************************/
/*****************************   OBSERVATION INPUT HANDLING   **********************************************/
/*****************************                                     **********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/



#if defined(GMTK_ARG_STREAMING_OUTPUT)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   bool binaryOutputStream = true;
   bool outputNetByteOrder = true;
   bool prettyPrintStream  = false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("\n*** Output observation stream ***\n"),
  Arg("binaryOutputStream", Arg::Opt, binaryOutputStream, "Output observation stream data is binary"),
  Arg("outputNetByteOrder", Arg::Opt, outputNetByteOrder, "For binary output streams, data is big-endian"),
  Arg("prettyPrint", Arg::Opt, prettyPrintStream, "For ASCII output streams, include segment and frame numbers (note that  pretty-printed streams cannot be sent to gmtkOnline)"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (prettyPrintStream && binaryOutputStream) {
    error("%s: cannot pretty-print binary observation streams", argerr);
  }
  if (binaryOutputStream && !outputNetByteOrder) {
    warning("using non-standard byte order for output binary observation stream data");
  }
  // FIXME - error if stream & file(s) both specified
  // warn? if options not relevant to stype specified

#else
#endif
#endif // defined(GMTK_ARG_STREAMING_OUTPUT)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_STREAMING_INPUT)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   bool binaryInputStream = true;
   bool inputNetByteOrder  = true;

   char const *os = "-";

// FIXME - consider adding #float, #int to protocol header
   unsigned snf = 0;
   unsigned sni = 0;
   char const *inputSource = "stream";
   bool streamSource = true;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("\n*** Input observation stream ***\n"),
  Arg("inputSource", Arg::Opt,inputSource,"Observation input source (file, stream)"),
  Arg("os",  Arg::Opt,os,"Input observation stream file name (- is stdin)"),
  Arg("binaryInputStream", Arg::Opt, binaryInputStream, "Input observation stream data is binary"),
  Arg("inputNetByteOrder", Arg::Opt, inputNetByteOrder, "For binary input observation streams, data is big-endian"),
  Arg("snf",  Arg::Opt,snf,"Number of floats in observation stream"),
  Arg("sni",  Arg::Opt,sni,"Number of ints in observation stream"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (binaryInputStream && !inputNetByteOrder) {
    warning("using non-standard byte order for binary input observation stream data");
  }

  // FIXME - should this be a warning?
  if (!strcasecmp(inputSource,"stream")) {
    streamSource = true;
  } else if (!strcasecmp(inputSource, "file")) {
    streamSource = false;
  } else {
    error("%s: unknown inputSource '%s', must be 'stream' or 'file'", argerr, inputSource);
  }

  if (streamSource && snf + sni == 0) {
    error("%s: frames have zero features",argerr);
  }
  // FIXME - error if stream & file(s) both specified
  // warn? if options not relevant to stype specified

#else
#endif
#endif // defined(GMTK_ARG_STREAMING_INPUT)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_OBS_FILES)
#if defined(GMTK_ARGUMENTS_DEFINITION)


// This next code is used by a number of routines to compute and set
// the default endian swapping condition associated with the
// arguments. We figure out the Endian of the machine this is running
// on and set the swap defaults accordingly.

#define DEF_CODE_TO_COMPUTE_ENDIAN(DEFAULT_SWAP_VALUE)   \
  bool doWeSwap; \
  ByteEndian byteEndian = getWordOrganization(); \
  switch(byteEndian) { \
  case BYTE_BIG_ENDIAN: \
    doWeSwap=false; \
    break; \
  case BYTE_LITTLE_ENDIAN: \
     doWeSwap=true; \
     break; \
  default: \
    /* We weren't able to figure the Endian out.  Leave the swap defaults as they are. */ \
    doWeSwap=DEFAULT_SWAP_VALUE; \
  } \
  \
   for(int i=0; i<MAX_NUM_OBS_FILES; ++i) { \
     iswp[i]=doWeSwap; \
  }

#ifdef INTV_WORDS_BIGENDIAN
#define CODE_TO_COMPUTE_ENDIAN DEF_CODE_TO_COMPUTE_ENDIAN(true) 
#else
#define CODE_TO_COMPUTE_ENDIAN DEF_CODE_TO_COMPUTE_ENDIAN(false) 
#endif


   // observation input file handling
#define MAX_NUM_OBS_FILES (5)
   char    *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL, NULL,NULL }; 
   unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
   unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
   const char   *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile","pfile","pfile" };
   const char    *frs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
   const char    *irs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" }; 
   const char    *prefrs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
   const char    *preirs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
   const char    *sr[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
   // per stream frame range string before any tranformations are applied
   char  *prepr[MAX_NUM_OBS_FILES] = {NULL,NULL,NULL,NULL,NULL};   
   // per stream frame range string after per-stream transformations are applied
   char *postpr[MAX_NUM_OBS_FILES] = {NULL,NULL,NULL,NULL,NULL};   
   char *gpr_str                   = NULL;   // global final frame range string

extern bool ObservationsAllowNan;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  // observation input file handling
  Arg("\n*** Observation input file handling ***\n"),
#ifdef OPTIONAL_OBSERVATION_FILES
  Arg("of",  Arg::Opt,ofs,"Observation File.  Replace X with the file number",Arg::ARRAY,MAX_NUM_OBS_FILES),
#else
  Arg("of",  Arg::Req,ofs,"Observation File.  Replace X with the file number",Arg::ARRAY,MAX_NUM_OBS_FILES),
#endif
  Arg("nf",  Arg::Opt,nfs,"Number of floats in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ni",  Arg::Opt,nis,"Number of ints in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fmt", Arg::Opt,fmts,"Format (htk,binary,ascii,pfile) for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("iswp",Arg::Opt,iswp,"Endian swap condition for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("prefr",  Arg::Opt,prefrs,"Float range for observation file X (before transforms)",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("preir",  Arg::Opt,preirs,"Int range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fr",  Arg::Opt,frs,"Float range for observation file X (after transforms)",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ir",  Arg::Opt,irs,"Int range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("sr",  Arg::Opt,sr,"Sentence range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("prepr", Arg::Opt, prepr,"Pre Per-segment frame Range for obs file X before any transforms are applied",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("postpr",Arg::Opt, postpr,"Post Per-segment frame Range for obs file X after per-stream transforms are applied",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("gpr",   Arg::Opt, gpr_str,"Global Per-segment final frame Range"),
  Arg("obsNAN",   Arg::Opt, ObservationsAllowNan," True if observation files allow FP NAN values"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  ////////////////////////////////////////////
  // check for valid argument values.
  int nfiles = 0;
  unsigned ifmts[MAX_NUM_OBS_FILES];
  for (int i=0;i<MAX_NUM_OBS_FILES;i++) {

    if (strcmp(fmts[i],"htk") == 0)
      ifmts[i] = HTK;
    else if (strcmp(fmts[i],"binary") == 0)
      ifmts[i] = RAWBIN;
    else if (strcmp(fmts[i],"ascii") == 0)
      ifmts[i] = RAWASC;
    else if (strcmp(fmts[i],"flatascii") == 0)
      ifmts[i] = FLATASC;
    else if (strcmp(fmts[i], "hdf5") == 0) 
      ifmts[i] = HDF5;
    else if (strcmp(fmts[i],"pfile") == 0)
      ifmts[i] = PFILE;
    else
      error("%s: Unknown observation file format type: '%s'\n",argerr,fmts[i]);

    if (ofs[i] != NULL && ifmts[i]!=PFILE && nfs[i] == 0 && nis[i] == 0)
      error("%s: command line parameters must specify one of nf%d and ni%d as not zero",argerr,
	    i+1,i+1);
    
    if(ofs[i] != NULL && ifmts[i]==PFILE) {
      FILE *in_fp = fopen(ofs[i], "r");
      if (in_fp==NULL) 
	error("Couldn't open input pfile %s for reading.", ofs[i]);
      bool debug_level=0;
      InFtrLabStream_PFile* in_streamp = new InFtrLabStream_PFile(debug_level,"",in_fp,1,iswp[i]);
      unsigned num_labs=in_streamp->num_labs();
      unsigned num_ftrs=in_streamp->num_ftrs();

      ////////////////////////////////////////////////////////////
      // Check consistency between pfile and supplied arguments //
      char search_str[]="nXXXXX";
      sprintf(search_str,"-ni%d",i+1);
      bool found=false;
      for(int j=1; j < argc; ++j) {
	if(strcmp(argv[j],search_str)==0) found=true;
      }
      if(found && nis[i] != num_labs) 
	error("%s: command line parameter ni%d (%d) is different from the one found in the pfile (%d)",argerr,
	      i+1,nis[i],num_labs); 
      sprintf(search_str,"-nf%d",i+1);
      found=false;
      for(int j=1; j < argc; ++j) {
	if(strcmp(argv[j],search_str)==0) found=true;
      }
      if(found && nfs[i] != num_ftrs) 
	error("%s: command line parameter nf%d (%d) is different from the one found in the pfile (%d)",
	      argerr,i+1,nfs[i],num_ftrs); 
      ////////////////////////////////////////////////////////////
      nis[i]=num_labs;
      nfs[i]=num_ftrs;

      if (fclose(in_fp)) 
	error("Couldn't close input pfile %s.", ofs[i]);
      delete in_streamp;
    }
    
    nfiles += (ofs[i] != NULL);
  }

#else
#endif
#endif // defined(GMTK_ARG_OBS_FILES)



/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************      FILE RANGE OPTIONS             ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

#if defined(GMTK_ARG_FILE_RANGE_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** File range options ***\n"),
#endif
#endif

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_DCDRNG)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  const static char *dcdrng_str="all";

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("dcdrng",Arg::Opt,dcdrng_str,"Range to decode over segment file"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_DCDRNG)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_TRRNG)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static const char *trrng_str="all";

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("trrng",Arg::Opt,trrng_str,"Range to decode over segment file"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_TRRNG)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_START_END_SKIP)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static int startSkip = 0;
  static int endSkip = 0;


#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("startSkip",Arg::Opt,startSkip,"Frames to skip at beginning (i.e., first frame is buff[startSkip])"),
  Arg("endSkip",Arg::Opt,endSkip,"Frames to skip at end (i.e., last frame is buff[len-1-endSkip])"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (startSkip < 0 || endSkip < 0)
    error("%s: arguments startSkip=%d/endSkip=%d must both be >= 0",argerr,startSkip,endSkip);

#else
#endif
#endif // defined(GMTK_ARG_START_END_SKIP)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_HELP)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   // 0: no help; HIGHEST_PRIORITY (1) ... LOWEST_PRIORITY (5) : increasing levels of help.  The priority levels are defined in arguments.h 
   static unsigned help = 0;  

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("help",  Arg::Help, help,  "Print this message. Add an argument from 1 to 5 for increasing help info."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if(help) {
    Arg::usage();
    exit(0);
  }

#else
#endif
#endif // defined(GMTK_ARG_HELP)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_VERSION)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool print_version_and_exit = false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("version",Arg::Opt,print_version_and_exit,"Print GMTK version number and exit."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (print_version_and_exit) {
#ifdef HAVE_CONFIG_H
    printf("%s (Mercurial id: %s)\n",gmtk_version_id,HGID);
#else
    printf("%s\n", gmtk_version_id);
#endif
    exit(0);
  }


#else
#endif
#endif // defined(GMTK_ARG_VERSION)




/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/************************                                              ******************************************/
/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
/************************                                              ******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

#if defined(GMTK_ARG_OBS_MATRIX_XFORMATION)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Observation matrix transformation options ***\n"),
#endif
#endif

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_OBS_MATRIX_XFORMATION)
#if defined(GMTK_ARGUMENTS_DEFINITION)

bool     Cpp_If_Ascii        = false;

const char*    Action_If_Diff_Num_Frames_Str[MAX_NUM_OBS_FILES]={"er","er","er","er","er"};   // 
unsigned Action_If_Diff_Num_Frames[MAX_NUM_OBS_FILES]={FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR};   // 
const char*    Action_If_Diff_Num_Sents_Str[MAX_NUM_OBS_FILES]={"te","te","te","te","te"}; 
unsigned Action_If_Diff_Num_Sents[MAX_NUM_OBS_FILES]={SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END};   // 

char    *Per_Stream_Transforms[MAX_NUM_OBS_FILES]={NULL,NULL,NULL,NULL,NULL};   // 
char    *Post_Transforms=NULL;

const char    *Ftr_Combo_Str="none";
unsigned Ftr_Combo=FTROP_NONE;
 
#ifdef INTV_WORDS_BIGENDIAN
bool iswp[MAX_NUM_OBS_FILES] = {true,true,true,true,true};
#else
bool iswp[MAX_NUM_OBS_FILES] = {false,false,false,false,false};
#endif 

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("fdiffact",  Arg::Opt, Action_If_Diff_Num_Frames_Str ,"Action if different number of frames in streams: error (er), repeat last frame (rl), first frame (rf), segmentally expand (se), truncate from start (ts), truncate from end (te)",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("sdiffact",  Arg::Opt, Action_If_Diff_Num_Sents_Str ,"Action if different number of sentences in streams: error (er), truncate from end (te), repeat last sent (rl), and wrap around (wa).",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("cppifascii",Arg::Tog, Cpp_If_Ascii,"Pre-process ASCII files using CPP"),
  Arg("trans",     Arg::Opt,Per_Stream_Transforms ,"per stream transformations string",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("posttrans", Arg::Opt,Post_Transforms ,"Final global transformations string"),
  Arg("comb",      Arg::Opt, Ftr_Combo_Str,"Combine float features (none: no combination, add, sub, mul,div"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (strcmp(Ftr_Combo_Str,"none") == 0)     Ftr_Combo = FTROP_NONE;
  else if (strcmp(Ftr_Combo_Str,"add") == 0) Ftr_Combo = FTROP_ADD;
  else if (strcmp(Ftr_Combo_Str,"sub") == 0) Ftr_Combo = FTROP_SUB;
  else if (strcmp(Ftr_Combo_Str,"mul") == 0) Ftr_Combo = FTROP_MUL;
  else if (strcmp(Ftr_Combo_Str,"div") == 0) Ftr_Combo = FTROP_DIV;
  else error("%s: Unknown feature combination type: '%s'\n",argerr,Ftr_Combo_Str);
  
  for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
    if(ofs[i]!=NULL) {
      if (strcmp(Action_If_Diff_Num_Frames_Str[i],"er") == 0)      Action_If_Diff_Num_Frames[i] = FRAMEMATCH_ERROR;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rl") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_REPEAT_LAST;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rf") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_REPEAT_FIRST;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"se") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_EXPAND_SEGMENTALLY;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"ts") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_TRUNCATE_FROM_START;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"te") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_TRUNCATE_FROM_END;
      else error("%s: Unknown action when diff num of frames: '%s'\n",argerr,Action_If_Diff_Num_Frames_Str[i]);
    }
  }
  
  for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
    if(ofs[i]!=NULL) {
      if (strcmp(Action_If_Diff_Num_Sents_Str[i],"er") == 0)      Action_If_Diff_Num_Sents[i] = SEGMATCH_ERROR;
      else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"rl") == 0) Action_If_Diff_Num_Sents[i] = SEGMATCH_REPEAT_LAST;
      else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"wa") == 0) Action_If_Diff_Num_Sents[i] = SEGMATCH_WRAP_AROUND;
      else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"te") == 0) Action_If_Diff_Num_Sents[i] = SEGMATCH_TRUNCATE_FROM_END;
      else error("%s: Unknown action when diff num of sentences: '%s'\n",argerr,
		 Action_If_Diff_Num_Sents_Str[i]);
    }
  }


#else
#endif
#endif // defined(GMTK_ARG_OBS_MATRIX_XFORMATION)


#if defined(GMTK_ARG_CPP_CMD_OPTS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *cppCommandOptions = NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Command line options to give to 'cpp'"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CPP_CMD_OPTS)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

