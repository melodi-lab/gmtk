/*-
 *-----------------------------------------------------------------------
 * A simple library to parse command lines in C++ using
 * an easy interface to quickly define arguments.
 * 
 *   Jeff Bilmes <bilmes@ee.washington.edu>
 *
 *  $Header$
 *
 *-----------------------------------------------------------------------
 */

// #include <strstream.h>
#include <iostream.h>
#include <fstream.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "error.h"
#include "arguments.h"


/*-
 *-----------------------------------------------------------------------
 * static entities
 *-----------------------------------------------------------------------
 */

bool Arg::EMPTY_ARGS_FLAG;
char* const Arg::NOFLAG = "noflg";
char* const Arg::NOFL_FOUND = "nofl_fnd";
const char Arg::COMMENTCHAR = '#';
int Arg::Num_Arguments = 0;
bool* Arg::Argument_Specified = NULL;
char* Arg::Program_Name = "";
const char* const ARGS_FILE_NAME = "argsFile";
const char* const ArgsErrStr = "Argument Error:";

/*-
 *-----------------------------------------------------------------------
 * CONSTRUCTORS
 *-----------------------------------------------------------------------
 */



/*-
 *-----------------------------------------------------------------------
 * Arg::initialize
 *      initialize the structures.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void Arg::initialize(char*m,ArgDisposition r,char *d) 
{
  flag = m;
  arg_kind = r;
  if (arg_kind == Tog) { 
    // type must be boolean. Toggle is never a required arg.
    if (mt.type != MultiType::bool_type) {
      error("%s A toggle argument must boolean",ArgsErrStr);
    }
  }
  description = d;
}


/*-
 *-----------------------------------------------------------------------
 * Arg::Arg()
 *      creates a special empty argument structure
 *      for using at the end of a list.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
Arg::Arg() : mt(EMPTY_ARGS_FLAG) {
  initialize(NULL,Arg::Opt,"");
}



/*-
 *-----------------------------------------------------------------------
 * Arg::Arg()
 *      creates a complete args object, with all bells and whistles.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
Arg::Arg(char*m ,ArgDisposition r,MultiType ucl,char* d)
: mt(ucl) 
{
  initialize(m,r,d);
}


/*-
 *-----------------------------------------------------------------------
 * Arg::Arg()
 *      creates a flagless argument.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
Arg::Arg(ArgDisposition r,MultiType ucl,char* d)
: mt(ucl) 
{
  initialize(NOFLAG,r,d);
}


/*-
 *-----------------------------------------------------------------------
 * Arg::Arg()
 *      copy constructor.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
Arg::Arg(const Arg& a)
: mt(a.mt) {
  flag = a.flag;
  arg_kind = a.arg_kind;
  description = a.description;
}

Arg::~Arg() {}



/*-
 *-----------------------------------------------------------------------
 * countAndClearArgBits()
 *      count the number of arguments that are used, and clear any of the
 *      used bits (if this hasn't been called before)
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void Arg::countAndClearArgBits() {
  if (Argument_Specified != NULL)
    return;
  Num_Arguments = 0;
  Arg* arg_ptr = Args;
  while (arg_ptr->flag != NULL) {
    Num_Arguments ++;
    arg_ptr++;
  }
  Argument_Specified = new bool[Num_Arguments];
  for (int i=0;i<Num_Arguments;i++) {
    Argument_Specified[i] = false;
  }
}


/*-
 *-----------------------------------------------------------------------
 * checkMissing()
 *      check if any required arguments are missing.
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      return true if any are missing.
 *
 *-----------------------------------------------------------------------
 */
bool Arg::checkMissing(bool printMessage) {
  Arg* arg_ptr;
  arg_ptr = Args;
  bool missing = false;
  while (arg_ptr->flag != NULL) {
    if (!Argument_Specified[arg_ptr-Args] && 
	arg_ptr->arg_kind == Req) {
      if (printMessage) {
	char brackets[2];
	fprintf(stderr,"%s Missing REQUIRED argument:",ArgsErrStr);
	if (!noFlagP(arg_ptr->flag))
	  fprintf(stderr," -%s ",arg_ptr->flag);
	// bool_type args are always optional. (i.e. -b T, -b F, or -b)
	if (arg_ptr->mt.type == MultiType::bool_type) {  
	  brackets[0] = '[';  brackets[1] = ']';
	} else {
	  brackets[0] = '<'; brackets[1] = '>';
	}
	fprintf(stderr," %c%s%c\n",
		brackets[0],
		arg_ptr->mt.printable(arg_ptr->mt.type),
		brackets[1]);
      }
      missing = true;
    }
    arg_ptr++;
  }
  return missing;
}


/*-
 *-----------------------------------------------------------------------
 * searchArgs()
 *      search the arg list one that matches a given flag.
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      Return NULL for unknown switch, return (Arg*)(-1) for ambiguous switch,
 *      otherwise returns the Arg* that matches.
 *
 *-----------------------------------------------------------------------
 */
Arg* Arg::searchArgs(Arg* ag,char *flag) {
  int flaglen = ::strlen(flag);
  Arg* arg_ptr = ag;
  int numTaged = 0;  
  int lastTaged = -1;

  // find the one that best matches.
  while (arg_ptr->flag != NULL) {
    if (!noFlagP(arg_ptr->flag)) {
      if (!::strncmp(arg_ptr->flag,flag,flaglen)) {
	numTaged++;
	lastTaged = (arg_ptr - ag);
      }
    }
    arg_ptr++;
  }
  // include check for args file since we shouldn't
  // be ambiguous with respect to this as well.
  if (!::strncmp(ARGS_FILE_NAME,flag,flaglen))
    numTaged++;

  if (numTaged == 0)
    return NULL;
  else if (numTaged > 1)
    return (Arg*)(-1);
  else 
    return &ag[lastTaged];
}


/*-
 *-----------------------------------------------------------------------
 * argsSwitch()
 *      checks the argument to see if it has the appropriate value
 *      for a given flag (which already was found out).
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      only returns Arg::ARG_OK or ARG_ERROR 
 *
 *-----------------------------------------------------------------------
 */
Arg::ArgsRetCode 
Arg::argsSwitch(Arg* arg_ptr,char *arg,int& index,bool& found,char*flag)
{
  switch(arg_ptr->mt.type) {
  // =============================================================
  case MultiType::int_type: {
    if (arg == NULL) { // end of arguments.
      warning("%s Integer argument needed: %s",
	      ArgsErrStr,
	      flag);
      return ARG_ERROR;
    }
    int n;
    char *endp;
    n = (int)strtol(arg,&endp,0); 
    if ( endp == arg ) {
      warning("%s Integer argument needed: %s %s",
	      ArgsErrStr,
	      flag,
	      arg);
      return ARG_ERROR;
    }
    arg_ptr->mt.ptr->integer = n;
    found = true;
  }
  break;
  // =============================================================
  case MultiType::uint_type: {
    if (arg == NULL) { // end of arguments.
      warning("%s Unsigned integer argument needed: %s",
	      ArgsErrStr,
	      flag);
      return ARG_ERROR;
    }
    unsigned int n;
    char *endp;
    n = (unsigned)strtoul(arg,&endp,0);
    if ( endp == arg) {
      warning("%s Unsigned integer argument needed: %s %s",
	      ArgsErrStr,
	      flag,
	      arg);
      return ARG_ERROR;
    }
    arg_ptr->mt.ptr->uinteger = n;
    found = true;
  }
  break;
  // =============================================================
  case MultiType::float_type: {
    if (arg == NULL) { // end of arguments.
      warning("%s Real number argument needed: %s",
	      ArgsErrStr,
	      flag);
      return ARG_ERROR;
    }
    float f;
    char *endp;
    f = strtof(arg,&endp);
    if ( arg == endp ) {
      warning("%s Floating point number argument needed: %s %s",
	      ArgsErrStr,
	      flag,
	      arg);
      return ARG_ERROR;
    }
    arg_ptr->mt.ptr->single_prec = f;
    found = true;
  }
  break;
  // =============================================================
  case MultiType::double_type: {
    if (arg == NULL) { // end of arguments.
      warning("%s Real number argument needed: %s",
	      ArgsErrStr,
	      flag);
      return ARG_ERROR;
    }
    double d;
    char *endp;
    d = strtod(arg,&endp);
    if (arg == endp) {
      warning("%s Integer argument needed: %s %s",
	      ArgsErrStr,
	      flag,
	      arg); 
      return ARG_ERROR;
    }
    arg_ptr->mt.ptr->double_prec = d;
    found = true;
  }
  break;
  // =============================================================
  case MultiType::str_type: {
    if (arg == NULL) { // end of arguments.
      warning("%s String argument needed: %s",
	      ArgsErrStr,
	      flag);
      return ARG_ERROR;
    }
    arg_ptr->mt.ptr->string =
      ::strcpy(new char[strlen(arg)+1],arg);
    found = true;
  }
  break;
  // =============================================================
  case MultiType::char_type: {
    if (arg == NULL || // end of arguments.
	::strlen(arg) != 1) {
      warning("%s Character argument needed: %s %s",
	      ArgsErrStr,
	      flag,arg);
      return ARG_ERROR;
    }
    arg_ptr->mt.ptr->ch = arg[0];
    found = true;
  }
  break;
  // =============================================================
  case MultiType::bool_type: {
    bool b; 
    if (arg == NULL || // end of arguments.
	arg[0] == '-') { // for optionless flag, just turn on.
      // for bool case, just turn it on.
      if (arg_ptr->arg_kind != Tog)
	arg_ptr->mt.ptr->boolean = true;
      else
	arg_ptr->mt.ptr->boolean = 
	  ((arg_ptr->mt.ptr->boolean == true) ? false : true);
      if (arg != NULL)
	index--;
      found = true;
      break;	      
    }
    // If kind  is Tog and we have a valid boolean 
    // argument, then treat argument as normal explicit boolean argument.
    if (!validBoolean(arg,b)) {
      warning("%s Boolean argument needed: %s %s",
	      ArgsErrStr,
	      flag,arg);
      return ARG_ERROR;
    }
    arg_ptr->mt.ptr->boolean = b;
    found = true;
  }
  break;
  // =============================================================
  default:
    error("%s Unknown internal argument",ArgsErrStr);
    break;
  }
  return ARG_OK;
}


/*-
 *-----------------------------------------------------------------------
 * validBoolean()
 *      returns true (and its value) if string is a valid boolean value
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      boolean value
 *
 *-----------------------------------------------------------------------
 */
bool Arg::validBoolean(char *string,bool& value)
{
  bool rc;
  int arglen = strlen(string);
  char *upcasearg = new char[arglen+1];
  ::strcpy(upcasearg,string);
  for (int i=0;i<arglen;i++)
    upcasearg[i] = toupper(upcasearg[i]);
  if (!::strncmp("TRUE",upcasearg,arglen) 
      || !::strncmp("YES",upcasearg,arglen)  
      || !::strncmp("ON",upcasearg,arglen)  
      || (arglen == 1 && string[0] == '1'))  {
    value = true;
    rc = true;
  } 
  else if (!::strncmp("FALSE",upcasearg,arglen) 
	   || !::strncmp("NO",upcasearg,arglen) 
	   || !::strncmp("OFF",upcasearg,arglen) 
	   || (arglen == 1 && string[0] == '0')) {
    value = false;
    rc = true;
  } else {
    rc = false;
  }
  delete [] upcasearg;
  return rc;
}


// return true if there is no flag.
bool Arg::noFlagP(char *flg) {
  // cant test here equality with NOFLAG, so do the following instread.
  if (flg == NOFLAG || flg == NOFL_FOUND)
    return true;
  else
    return false;
}



/*-
 *-----------------------------------------------------------------------
 * Printing operations
 *-----------------------------------------------------------------------
 */


/*-
 *-----------------------------------------------------------------------
 *  MultiType::print(FILE* f)
 *      prints the argument value for the current arg.
 * 
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void MultiType::print(FILE* f) {
  switch (type) {
  case MultiType::bool_type:
    fprintf(f,"%c",(ptr->boolean?'T':'F'));
    break;
  case MultiType::char_type:
    fprintf(f,"%c",(ptr->ch));
    break;
  case MultiType::str_type:
    fprintf(f,"%s",(ptr->string == NULL ? "" : ptr->string));
    break;
  case MultiType::int_type:
    fprintf(f,"%d",(ptr->integer));
    break;
  case MultiType::uint_type:
    fprintf(f,"%u",(ptr->uinteger));
    break;
  case MultiType::float_type:
    fprintf(f,"%e",(double)(ptr->single_prec));
    break;
  case MultiType::double_type:
    fprintf(f,"%e",(ptr->double_prec));
    break;
  default:
    error("%s Internal argument error",ArgsErrStr);
    break;
  }
}

/*-
 *-----------------------------------------------------------------------
 * Arg::printArgs()
 *      prints all the args in the list.

 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void Arg::printArgs(Arg*args,FILE* f) {
  Arg *arg_ptr = args;
  while (arg_ptr->flag != NULL) {
    arg_ptr->print(f);
    arg_ptr++;
  }
}

/*-
 *-----------------------------------------------------------------------
 * Arg::print()
 *      prints the current argument entry.

 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void Arg::print(FILE* f) {
  if (!noFlagP(flag))
    fprintf(f,"%s",flag);
  fprintf(f,":%s = ",MultiType::printable(mt.type));
  mt.print(f);
  fprintf(f,"; # %s\n",
	  ((description == NULL) ? "" : description));
}


/*-
 *-----------------------------------------------------------------------
 * printable()
 *      returns a printable version of an argument type
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
char *MultiType::printable(MultiType::ArgumentType at) {
  switch (at) {
  case MultiType::bool_type:
    return "bool";
  case MultiType::char_type:
    return "char";
  case MultiType::str_type:
    return "str";
  case MultiType::int_type:
    return "integer";
  case MultiType::uint_type:
    return "unsigned";
  case MultiType::float_type:
    return "float";
  case MultiType::double_type:
    return "double";
  default:
    return "error: unknown type";
  }
}


/*-
 *-----------------------------------------------------------------------
 * usage()
 *      prints a usage message of the program. This will print out
 *      the arguments, their types, documentation, and default values.
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      stuff is printed out.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void Arg::usage() {

  fprintf(stderr,"Usage: %s  [[[-flag] [option]] ...]\n",Program_Name);
  fprintf(stderr,"Required: <>; Optional: []; Flagless arguments must be in order.\n");

  Arg* arg_ptr = Args;
  int longest_variation = 0;

  while (arg_ptr->flag != NULL) {
    int len = 0;
    if (!noFlagP(arg_ptr->flag)) {
      // add one for the '-', as in "-flag"
      len += ::strlen(arg_ptr->flag)+1;
      len ++; //  add one for the ' ' in "-flag "
    }
    len += ::strlen(MultiType::printable(arg_ptr->mt.type));
    len += 2; // add two for brackets. '[',']', or '<','>' around type.
    if (len  > longest_variation)
      longest_variation = len;
    arg_ptr++;
  }

  for (int printOptional=0;printOptional<2;printOptional++) {
    arg_ptr = Args;
    while (arg_ptr->flag != NULL) {
      int this_variation = 0;
      char brackets[2];
      if (arg_ptr->arg_kind == Req) {
	if (printOptional)
	  goto skip;
	brackets[0] = '<'; brackets[1] = '>';
      } else {
	if (!printOptional)
	  goto skip;
	brackets[0] = '['; brackets[1] = ']';
      }
      fprintf(stderr," %c",brackets[0]);

      if (!noFlagP(arg_ptr->flag)) {
	fprintf(stderr,"-%s",arg_ptr->flag);
	// add one for the '-', as in "-flag"
	this_variation = ::strlen(arg_ptr->flag) + 1;
	fprintf(stderr," ");
	this_variation ++; //  add one for the ' ' in "-flag "
      }
      fprintf(stderr,"%s",arg_ptr->mt.printable(arg_ptr->mt.type));

      this_variation += ::strlen(MultiType::printable(arg_ptr->mt.type));
      // add two for brackets. '[',']', or '<','>' around type.
      this_variation += 2;

      fprintf(stderr,"%c",brackets[1]);

      while (this_variation++ < longest_variation)
	fprintf(stderr," ");
      fprintf(stderr,"   %s {",
	      ((arg_ptr->description == NULL) ? "" : arg_ptr->description));
      arg_ptr->mt.print(stderr);      
      fprintf(stderr,"}\n");
    skip:
      arg_ptr++;
    }
  }

  fprintf(stderr," [-%s <str>]",ARGS_FILE_NAME);
  int this_variation = 9 + strlen(ARGS_FILE_NAME);
  while (this_variation++ < longest_variation)
    fprintf(stderr," ");
  fprintf(stderr,"   File to obtain additional arguments from {}\n");
}


/*-
 *-----------------------------------------------------------------------
 * parseArgsFromFile()
 *      more direct parses from directly from a file
 * 
 * Preconditions:
 *      file must be valid format.
 *
 * Postconditions:
 *      file is parsed, any previously defined argument values are lost.
 *
 * Side Effects:
 *      modifies internal object static variables.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
Arg::ArgsRetCode Arg::parseArgsFromFile(char *fileName)
{
  countAndClearArgBits();
  ifstream ifile(fileName);
  if (!ifile) {
    warning("%s Can't file argument file: %s",
	    ArgsErrStr,
	    fileName);
    return ARG_ERROR;
  } else {
    // get the remaining from file
    const unsigned max_line_length = 32*1024;
    char buffer[max_line_length];

    while (ifile.getline(buffer,max_line_length,'\n')) {
      // printf("Got line (%s)\n",buffer);
      const unsigned line_length = strlen(buffer);
      if (line_length+1 >= max_line_length) {
	// input string is too long.
 	warning("%s Line length too long in command line parameter arguments file: %s",
 		ArgsErrStr,
 		fileName);
 	return ARG_ERROR; // give up
      }
      char* buffp = buffer;
      while (*buffp) {
	if (*buffp == COMMENTCHAR) {
	  *buffp = '\0';
	  break;
	}
	buffp++;
      }

      buffp = buffer;
      while (*buffp && isspace(*buffp))
        buffp++;       // skip space
      if (!*buffp) {
	continue; // empty line
      }
      char *flag = buffp; // get flag
      char *arg;
      // get command up to space or ':'
      while (*buffp && *buffp != ' ' && *buffp != '\t' &&  *buffp != ':')
	buffp++;
      if (buffp == flag)
	continue; // empty line or empty flag
      if (*buffp) { 
	// get ':' and position buffp to start of arg.
	if (*buffp == ':')
	  // we have the flag
	  *buffp++ = '\0'; 
	else {
	  // we have the flag, but need to get rid of ':' if there.
	  *buffp++ = '\0'; 
	  while (*buffp && *buffp != ':')
	    buffp++;
	  if (*buffp == ':')
	    buffp++;
	}
	while (*buffp == ' ' || *buffp == '\t')
	  buffp++; 	// skip space
      }

      if (!*buffp)
	arg = NULL;
      else {
	arg = buffp;
	// get command up to space
	while (*buffp && *buffp != ' ' && *buffp != '\t')
	  buffp++;	
	*buffp = '\0';
      }      

      // check to see if it is the special parse from file
      // argument name, which is not valid here.
      if (!::strncmp(ARGS_FILE_NAME,flag,::strlen(flag)))
	error("%s Can not recursively parse arguments from files\n",ArgsErrStr);

      Arg* arg_ptr = searchArgs(Args,flag);
      if (arg_ptr == NULL) {
	warning("%s Unknown switch in %s : %s",
		ArgsErrStr,
		fileName,flag);
	return (ARG_ERROR);
      } else if (arg_ptr == (Arg*)(-1)) {
	warning("%s Ambiguous switch in %s : %s",
		ArgsErrStr,
		fileName,flag);
	return (ARG_ERROR);
      } else {
	int i;
	if (argsSwitch(arg_ptr,arg,i,Argument_Specified[arg_ptr - Args],flag) 
	    != ARG_OK) {
	  warning("%s Error in %s",
		  ArgsErrStr,
		  fileName);
	  return (ARG_ERROR);
	}
      }
    }
  }
  if (checkMissing())
    return ARG_MISSING;
  return ARG_OK;
}



/*-
 *-----------------------------------------------------------------------
 * parseArgsFromCommandLine()
 *      more direct parses from the command line
 * 
 * Preconditions:
 *      command line must be specified
 *
 * Postconditions:
 *      command line is parsed.
 *
 * Side Effects:
 *      modifies internal object static variables.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
Arg::ArgsRetCode 
Arg::parseArgsFromCommandLine(int argc,char**argv)
{
  if (argv[0] != NULL)
    Program_Name = argv[0];
  countAndClearArgBits();
  Arg* arg_ptr;
  for (int i=1;i<argc;i++) {
    if (argv[i][0] == '-') {
      char *flag = argv[i];
      int flaglen = ::strlen(flag);
      arg_ptr = searchArgs(Args,&flag[1]);
      if (arg_ptr == NULL) {
	warning("%s Unknown switch: %s",
		ArgsErrStr,
		argv[i]);
	return (ARG_ERROR);
      } else if (arg_ptr == (Arg*)(-1)) {
	warning("%s Ambiguous switch: %s",
		ArgsErrStr,
		argv[i]);
	return (ARG_ERROR);
      } else {
	// first check to see if it is the special parse from file
	// argument name, which is always valid.
	if (!::strncmp(ARGS_FILE_NAME,&flag[1],flaglen-1)) {
	  // so argument is presumably a file name.
	  char *arg=NULL;
	  if ((i+1) >= argc)
	    error("%s Expecting file name after %s argument flag\n",
		  ArgsErrStr,
		  ARGS_FILE_NAME);
	  else
	    arg = argv[++i];
	  parseArgsFromFile(arg);
	} else {
	  char *arg;
	  if ((i+1) >= argc)
	    arg = NULL;
	  else
	    arg = argv[++i];
	  argsSwitch(arg_ptr,arg,i,Argument_Specified[arg_ptr - Args],flag);
	}
      }
    } else { // Go through and look for no flag case, in order.
      arg_ptr = Args;
      while (arg_ptr->flag != NULL) {
	if (arg_ptr->flag == NOFLAG) {
	  // assume string type
	  char *arg = argv[i];
	  int dummy;
	  argsSwitch(arg_ptr,arg,dummy,Argument_Specified[arg_ptr-Args],"");
	  arg_ptr->flag = NOFL_FOUND;
	  break;
	}
	arg_ptr++;
      }
    }
  }
  if (checkMissing())
    return ARG_MISSING;
  return ARG_OK;
}


/*-
 *-----------------------------------------------------------------------
 * parse()
 *      parses the arguments given to the program
 * 
 * Preconditions:
 *      arguments must be specified
 *
 * Postconditions:
 *      arguments are parsed
 *
 * Side Effects:
 *      modifies internal object static variables.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void Arg::parse(int argc,char** argv)
{
  ArgsRetCode rc;
  rc = parseArgsFromCommandLine(argc,argv);
  if (rc == ARG_MISSING) { // still missing ??
    Arg::checkMissing(true);
  }
  if (rc != ARG_OK) {
    Arg::usage();
    ::exit(1);
  }
}


/*-
 *-----------------------------------------------------------------------
 * Test driver routine, -DMAIN to compile into program.
 *-----------------------------------------------------------------------
 */


#ifdef MAIN

/*
 * arguments without flags
 */

char *string_fl="This is a string";
char char_fl = 'C';
float float_fl = 3.4;
double double_fl = 4.5;
bool bool_fl = false;
int int_fl = 343;

/*
 * arguments with flags
 */
char *myString = "BARSTR";
float aSingle = 2.3;
double aDouble = .4;
int int1 = 3;
int int2 = 3;
bool bvalue = true;
bool aToggle = true;
bool anotherToggle = false;
char aChar = 'c';

/*
 * the argument list
 */
Arg Arg::Args[] = {

 // Arguments with flags. 
 Arg("myString",  Arg::Opt, myString,"A string"),
 Arg("aSingle",   Arg::Opt, aSingle,"A single precision floating point num"),
 Arg("aDouble",   Arg::Req, aDouble,"A double precision floating point num"),
 Arg("int1",      Arg::Req, int1,"An integer"),
 Arg("int2",      Arg::Req, int2,"A differnet integer"),
 Arg("rbvalue",   Arg::Req, bvalue,"A required boolean"),
 Arg("bvalue",    Arg::Opt, bvalue,"A boolean"),
 Arg("aToggle",   Arg::Tog, aToggle,"a toggle"),
 Arg("anotherToggle",Arg::Tog, anotherToggle,"another toggle"),
 Arg("aChar",     Arg::Opt, aChar,"A character"),

 // Arguments without flags. The order on the command line must be in
 // the the same as the order given here. i.e. string first, then int,
 // then bool, etc.
 Arg(Arg::Req, string_fl,  "string"),
 Arg(Arg::Opt, char_fl,    "char"),
 Arg(Arg::Opt, float_fl,   "float"),
 Arg(Arg::Opt, double_fl,  "double"),
 Arg(Arg::Opt, bool_fl,    "bool"),
 Arg(Arg::Opt, int_fl,     "int"),

 // The argumentless argument marks the end
 // of the above list.
 Arg()

};


int main(int argc,char*argv[])
{
  Arg::parse(argc,argv);
  Arg::printArgs(Arg::Args,stdout);
}

#endif // #ifdef MAIN
