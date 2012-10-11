

/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/*****************************                                     **********************************************/
/*****************************   OBSERVATION INPUT HANDLING   **********************************************/
/*****************************                                     **********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

#define MAX_NUM_OBS_FILES (5)

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

#else
#endif
#endif // defined(GMTK_ARG_STREAMING_OUTPUT)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

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


// for gmtkOnline? only accept a single input stream
#if defined(GMTK_ARG_STREAM_INPUT)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#define MAX_NUM_OBS_STREAMS (1)

#define DEFAULT_STREAM_BUFFER_SIZE (sizeof(Data32))

  bool inputNetByteOrder[MAX_NUM_OBS_STREAMS] = {true};
  char *oss[MAX_NUM_OBS_STREAMS] = {NULL};
  unsigned streamBufferSize = DEFAULT_STREAM_BUFFER_SIZE;

   // observation input file handling
   unsigned nfs[MAX_NUM_OBS_STREAMS] = {0};
   unsigned nis[MAX_NUM_OBS_STREAMS] = {0};
   const char   *fmts[MAX_NUM_OBS_STREAMS] = {"binary"};
   unsigned ifmts[MAX_NUM_OBS_STREAMS] = {RAWBIN};

extern bool ObservationsAllowNan;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("\n*** Observation stream input handling ***\n"),
  Arg("streamBufferSize", Arg::Opt,streamBufferSize,"Size in MB of observation frame queue"),
  Arg("os",  Arg::Req,oss,"Input observation stream file name (- is stdin). Replace X with the stream number",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("inputNetByteOrder", Arg::Opt, inputNetByteOrder, "For binary input observation stream X, data is big-endian",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("nf",  Arg::Opt,nfs,"Number of floats in observation stream X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("ni",  Arg::Opt,nis,"Number of ints in observation stream X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("fmt", Arg::Opt,fmts,"Format (for files: htk,binary,ascii,flatascii,hdf5,pfile; for streams: binary,ascii) for observation stream X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("obsNAN",   Arg::Opt, ObservationsAllowNan," True if observation files allow FP NAN values"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if ( streamBufferSize < 1) {
    error("%s: streamBufferSize must be at least 1 (got %u)", streamBufferSize);
  }
#define MEBIBYTE (1048576)
  streamBufferSize *= MEBIBYTE;
  bool gotStdin = false;
  for (int i=0;i<MAX_NUM_OBS_STREAMS;i++) {
    if (oss[i] && strcmp(oss[i],"-") == 0) {
      if (gotStdin) {
	error("%s: can only use stdin as an input stream once", argerr);
      }
      gotStdin = true;
    }
    if (strcmp(fmts[i],"binary") == 0)
      ifmts[i] = RAWBIN;
    else if (strcmp(fmts[i],"ascii") == 0)
      ifmts[i] = RAWASC;
    else
      error("%s: observation stream format -fmt%u must be 'binary' or 'ascii', got '%s'\n",argerr,i, fmts[i]);
    
    if (ifmts[i] == RAWBIN && !inputNetByteOrder[i]) {
      warning("using non-standard byte order for binary input observation stream %u data", i);
    }
  }

#else
#endif
#endif // defined(GMTK_ARG_XX_XX)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_STREAM_AND_FILE_INPUT)
#if defined(GMTK_ARGUMENTS_DEFINITION)

// stream stuff

#define MAX_NUM_OBS_STREAMS (5)

#define DEFAULT_STREAM_BUFFER_SIZE (16)
#define DEFAULT_FILE_BUFFER_SIZE   (16)
#ifndef DEFAULT_FILE_WINDOW_SIZE
#define DEFAULT_FILE_WINDOW_SIZE   (4)
#define DEFAULT_FILE_WINDOW_DELTA  (100)
#endif

  bool inputNetByteOrder[MAX_NUM_OBS_STREAMS] = {true,true,true,true,true};
  char *oss[MAX_NUM_OBS_STREAMS] = {NULL,NULL,NULL,NULL,NULL};
  unsigned streamBufferSize = DEFAULT_STREAM_BUFFER_SIZE;

// file stuff

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
   char    *ofs[MAX_NUM_OBS_STREAMS] = { NULL, NULL, NULL, NULL,NULL }; 
   unsigned nfs[MAX_NUM_OBS_STREAMS] = { 0, 0, 0,0,0 };
   unsigned nis[MAX_NUM_OBS_STREAMS] = { 0, 0, 0,0,0 };
   const char   *fmts[MAX_NUM_OBS_STREAMS] = {"binary","binary","binary","binary","binary"};
   unsigned ifmts[MAX_NUM_OBS_STREAMS];
   const char    *frs[MAX_NUM_OBS_STREAMS] = { "all", "all", "all","all","all" };
   const char    *irs[MAX_NUM_OBS_STREAMS] = { "all", "all", "all","all","all" }; 
   const char    *prefrs[MAX_NUM_OBS_STREAMS] = { NULL, NULL, NULL, NULL,NULL }; 
   const char    *preirs[MAX_NUM_OBS_STREAMS] = { NULL, NULL, NULL, NULL,NULL }; 
   const char    *sr[MAX_NUM_OBS_STREAMS] = { NULL, NULL, NULL, NULL,NULL }; 
   // per stream frame range string before any tranformations are applied
   char  *prepr[MAX_NUM_OBS_STREAMS] = {NULL,NULL,NULL,NULL,NULL};   
   // per stream frame range string after per-stream transformations are applied
   char *postpr[MAX_NUM_OBS_STREAMS] = {NULL,NULL,NULL,NULL,NULL};   
   unsigned fileBufferSize = DEFAULT_FILE_BUFFER_SIZE;
   unsigned fileWindowSize = DEFAULT_FILE_WINDOW_SIZE;
   unsigned fileWindowDelta = DEFAULT_FILE_WINDOW_DELTA;
   bool constantSpace = false;
   char *gpr_str                   = NULL;   // global final frame range string
   const char *justification_str         = "left";
   unsigned justification;

#ifdef INTV_WORDS_BIGENDIAN
   bool iswp[MAX_NUM_OBS_STREAMS] = {true,true,true,true,true};
#else
   bool iswp[MAX_NUM_OBS_STREAMS] = {false,false,false,false,false};
#endif 

extern bool ObservationsAllowNan;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("\n*** Observation stream input handling ***\n"),
  Arg("streamBufferSize", Arg::Opt,streamBufferSize,"Size in MB of the stream observation frame queue"),
  Arg("os",  Arg::Opt,oss,"Input observation stream file name (- is stdin). Replace X with the stream number",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("inputNetByteOrder", Arg::Opt, inputNetByteOrder, "For binary input observation stream X, data is big-endian",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("of",  Arg::Opt,ofs,"Read observation stream X from non-stream file.",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("fmt", Arg::Opt,fmts,"Format (for files: htk,binary,ascii,flatascii,hdf5,pfile; for streams: binary,ascii) for observation stream X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("fileBufferSize", Arg::Opt,fileBufferSize,"Size in MB of the file observation frame buffer"),
  Arg("constantSpace", Arg::Opt,constantSpace,"Use only fileBufferSize memory to hold the observation data"),
  Arg("fileWindowSize", Arg::Opt,fileWindowSize, "Size in MB to load at once if constantSpace is active"),
  Arg("fileWindowDelta", Arg::Opt,fileWindowDelta, "How close (in frames) from the edge of the current window triggers loading more frames"),    
  Arg("nf",  Arg::Opt,nfs,"Number of floats in observation stream X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("ni",  Arg::Opt,nis,"Number of ints in observation stream X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("iswp",Arg::Opt,iswp,"Endian swap condition for observation file X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("prefr",  Arg::Opt,prefrs,"Float range for observation file X (before transforms)",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("preir",  Arg::Opt,preirs,"Int range for observation file X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("fr",  Arg::Opt,frs,"Float range for observation stream X (after transforms)",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("ir",  Arg::Opt,irs,"Int range for observation stream X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("sr",  Arg::Opt,sr,"Sentence range for observation file X",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("prepr", Arg::Opt, prepr,"Pre Per-segment frame Range for obs file X before any transforms are applied",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("postpr",Arg::Opt, postpr,"Post Per-segment frame Range for obs file X after per-stream transforms are applied",Arg::ARRAY,MAX_NUM_OBS_STREAMS),
  Arg("gpr",   Arg::Opt, gpr_str,"Global Per-segment final frame Range"),
  Arg("justification", Arg::Opt, justification_str, "Justification of usable frames (left, center, right)"),
  Arg("obsNAN",   Arg::Opt, ObservationsAllowNan,"True if observation files allow FP NAN values"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#define MEBIBYTE (1048576)
    streamBufferSize *= MEBIBYTE / sizeof(Data32);
    fileBufferSize *= MEBIBYTE / sizeof(Data32);

  /////////////////////////////////////////////////////////
  // check for valid argument values.  File only case 
  int nfiles = 0;
  bool gotStdin = false;
  for (int i=0;i<MAX_NUM_OBS_STREAMS;i++) {

    if (ofs[i] && oss[i]) {
      error("%s: Can't specify both -of%u and -os%u", argerr,i+1,i+1);
    }

    if ( (ofs[i] && strcmp(ofs[i],"-")==0) ||
         (oss[i] && strcmp(oss[i],"-")==0) )
    {
      if (gotStdin) {
	error("%s: can only use stdin as an input once", argerr);
      }
      gotStdin = true;
    }

    if (ofs[i]) {
      int fmtNum = formatStrToNumber(fmts[i]);
      if (fmtNum >= 0) {
	ifmts[i] = fmtNum;
      } else {
	error("%s: Unknown observation file format type '%s' for -of%u\n",argerr,fmts[i],i+1);
      }
    } else {
      if (prefrs[i] || preirs[i]) {
	error("%s: Can't specify -prefrs%u or -preirs%u because -os%u is a stream",argerr,i+1,i+1,i+1);
      }
      if (sr[i]) {
	error("%s: Can't specify -sr%u because -os%u is a stream", argerr,i+1,i+1);
      }
      if (prepr[i] || postpr[i]) {
	error("%s: Can't specify -prepr%u or -postpr%u because -os%u is a stream", argerr,i+1,i+1,i+1);
      }
      if (strcmp(fmts[i],"binary") == 0)
	ifmts[i] = RAWBIN;
      else if (strcmp(fmts[i],"ascii") == 0)
	ifmts[i] = RAWASC;
      else
	error("%s: observation stream format -fmt%u must be 'binary' or 'ascii', got '%s'\n",argerr,i+1, fmts[i]);
      
      if (ifmts[i] == RAWBIN && !inputNetByteOrder[i]) {
	warning("using non-standard byte order for binary input observation stream %u data", i+1);
      }
    }

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
    
    nfiles += (ofs[i] != NULL || oss[i] != NULL);
  }

  if (strcmp(justification_str, "left") == 0) {
    justification = FRAMEJUSTIFICATION_LEFT;
  } else if (strcmp(justification_str, "center") == 0) {
    justification = FRAMEJUSTIFICATION_CENTER;
  } else if (strcmp(justification_str, "right") == 0) {
    justification = FRAMEJUSTIFICATION_RIGHT;
  } else {
    error("%s: must specify 'left', 'center', or 'right' for -justification",argerr);
  }


#else
#endif
#endif // defined(GMTK_ARG_STREAM_AND_FILE_INPUT)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_OBS_FILES)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#define DEFAULT_FILE_BUFFER_SIZE   (16)
#ifndef DEFAULT_FILE_WINDOW_SIZE
#define DEFAULT_FILE_WINDOW_SIZE   (4)
#define DEFAULT_FILE_WINDOW_DELTA  (100)
#endif

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
#ifndef MAX_NUM_OBS_FILES
#define MAX_NUM_OBS_FILES (5)
#endif
   char    *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL, NULL,NULL }; 
   unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
   unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
   const char   *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile","pfile","pfile" };
   unsigned ifmts[MAX_NUM_OBS_FILES];
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
   const char *justification_str         = "left";
   unsigned justification;
   unsigned fileBufferSize = DEFAULT_FILE_BUFFER_SIZE;
   unsigned fileWindowSize = DEFAULT_FILE_WINDOW_SIZE;
   unsigned fileWindowDelta = DEFAULT_FILE_WINDOW_DELTA;
   bool constantSpace = false;

#ifdef INTV_WORDS_BIGENDIAN
   bool iswp[MAX_NUM_OBS_FILES] = {true,true,true,true,true};
#else
   bool iswp[MAX_NUM_OBS_FILES] = {false,false,false,false,false};
#endif 


extern bool ObservationsAllowNan;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  // observation input file handling
  Arg("\n*** Observation input file handling ***\n"),
  Arg("of",  Arg::Req,ofs,"Observation File.  Replace X with the file number",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("nf",  Arg::Opt,nfs,"Number of floats in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ni",  Arg::Opt,nis,"Number of ints in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fmt", Arg::Opt,fmts,"Format (htk,binary,ascii,flatascii,hdf5,pfile) for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fileBufferSize", Arg::Opt,fileBufferSize,"Size in MB of the file observation frame buffer"),
  Arg("constantSpace", Arg::Opt,constantSpace,"Use only fileBufferSize memory to hold the observation data"),
  Arg("fileWindowSize", Arg::Opt,fileWindowSize, "Size in MB to load at once if constantSpace is active"),
  Arg("fileWindowDelta", Arg::Opt,fileWindowDelta, "How close (in frames) from the edge of the current window triggers loading more frames"),    
  Arg("iswp",Arg::Opt,iswp,"Endian swap condition for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("prefr",  Arg::Opt,prefrs,"Float range for observation file X (before transforms)",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("preir",  Arg::Opt,preirs,"Int range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fr",  Arg::Opt,frs,"Float range for observation file X (after transforms)",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ir",  Arg::Opt,irs,"Int range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("sr",  Arg::Opt,sr,"Sentence range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("prepr", Arg::Opt, prepr,"Pre Per-segment frame Range for obs file X before any transforms are applied",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("postpr",Arg::Opt, postpr,"Post Per-segment frame Range for obs file X after per-stream transforms are applied",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("gpr",   Arg::Opt, gpr_str,"Global Per-segment final frame Range"),
  Arg("justification", Arg::Opt, justification_str, "Justification of usable frames (left, center, right)"),
  Arg("obsNAN",   Arg::Opt, ObservationsAllowNan,"True if observation files allow FP NAN values"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#define MEBIBYTE (1048576)

    fileBufferSize *= MEBIBYTE / sizeof(Data32);

  /////////////////////////////////////////////////////////
  // check for valid argument values.  File only case 
  int nfiles = 0;
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

  if (strcmp(justification_str, "left") == 0) {
    justification = FRAMEJUSTIFICATION_LEFT;
  } else if (strcmp(justification_str, "center") == 0) {
    justification = FRAMEJUSTIFICATION_CENTER;
  } else if (strcmp(justification_str, "right") == 0) {
    justification = FRAMEJUSTIFICATION_RIGHT;
  } else {
    error("%s: must specify 'left', 'center', or 'right' for -justification",argerr);
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

  int startSkip = 0;
  int endSkip = 0;


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

  Arg("\n*** Help options ***\n"),
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

#if !defined(MAX_NUM_OBS_FILES)
#  if defined(MAX_NUM_OBS_STREAMS)
#    define MAX_NUM_OBS_FILES MAX_NUM_OBS_STREAMS
#  endif
#endif

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

  char *cppCommandOptions = NULL;

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

