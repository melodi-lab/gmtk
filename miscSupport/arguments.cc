//
// *************************************************************************
// * A general mechanism to parse command line and file arguments for C++.
// * Jeff Bilmes
// * <bilmes@media.mit.edu> Aug 1992
// * See main() in arguments.cc for example.
// * Compile with -DMAIN to demonstrate example.
// *************************************************************************
// 

// TODO and Known Bugs: 
//   o A printDefault message should wait to print until
//    all args are set (including those in files).
//    Add another argument to parseFromCommandLine so it can tell what
//    is the file to get input from. Do we want multiple such files?
//   o Add ability to get parameters from stdin.
//   o Fix bug with flagless and boolean arguments where
//     we can do: program -boolarg flaglessargument
//   o add the ability to have multiple arguments for one flag.
//     i.e -flag <arg1> <arg2> ... <argn>
//       ARGS("foo",ARGS::Req,a,b,c,d,e,"Takes (requires) 5 arguments");
//     (or even an unknown number of arguments of one type,
//         i.e. an array of ints, or of char*s).
//       char **ptrp;
//       int size_ptrp;
//       ARGS("foo",ARGS::Req || ARGS:Unlimited,ptrp,size_ptrp,"Takes many char*s");
//     isnumArg(char*) and isoptionArg(char*) routines will be useful.
//   o Add ability to unsatisfy arguments (i.e. if a -q flag implies 
//     needing a -p flag).
//   o Add actions that get called after all arguments are set from
//     command line, file, etc. so that user can check for special 
//     requirements.
//   o add types short,long,unsigned short, unsigned long, etc.
//   o add ability for usage message to also print out default values.
//   o Order of parsing arguments should be correct. I.e. in order
//     of appearance on command line, and if command line is a file
//     argument, that file should be parsed before subsequent command line
//     arguments.
//   o add bit arguments, where boolean flags can set or reset just one
//     specified bit of an integer argument.
//   o add ability to do "one of" arguments, e.g. a direction
//        argument, "progr -dir up" where dir is checked to
//        be up,down,left, or right and an error results if it isn't.
//     possible oneof syntax:
//         ARGS(
//              ARGS("arg1",...),
//              ARGS("arg2",...),
//               ..
//                 );
//   o Arguments file should be piped through cpp. Comments should start with ";"
//   o Make the entire package non-static so multiple args structs
//     can be created.
//   o Allow for different sets of valid arguments to be selected
//     based on previous command line parameters.

#include <strstream.h>
#include <iostream.h>
#include <fstream.h>
#include <stddef.h>
#include <stdlib.h>

/* files included by std.h
#include <_G_config.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
*/

#include <string.h>
#include <ctype.h>


#include "arguments.h"

bool ARGS::__args__dummy__;
ARGS ARGS::END(NULL,ARGS::Opt,ARGS::__args__dummy__, "");

char* const ARGS::NOFLAG = "noflg";
char* const ARGS::NOFL_FOUND = "nofl_fnd";
bool ARGS::usage_called = false;
const char ARGS::CommentChar = '#';
int ARGS::numArgs = 0;
bool* ARGS::found = NULL;
char* ARGS::progName = "";
argsAction ARGS::no_action;
bool ARGS::ignoreUnknownSwitch = false;

void ARGS::init(char*m,argKind r,char *d,argsAction& a) 
{
  flag = m;
  arg_kind = r;
  if (arg_kind == Tog) { // type must be boolean. Toggle is never a required arg.
    if (uc.type != unionClass::bool_t) {
      Error("Internal: ARGS::ARGS() toggle argument must be bool");
      ::exit(1);
    }
  }
  description = d;
  if (&a == &no_action)
    action = NULL;
  else
    action = &a;
}

ARGS::ARGS() : uc(__args__dummy__) {
  init(NULL,ARGS::Opt,"",no_action);
}


ARGS::ARGS(char*m ,argKind r,unionClass ucl,char* d,argsAction& a)
: uc(ucl) 
{
  init(m,r,d,a);
}


ARGS::ARGS(argKind r,unionClass ucl,char* d,argsAction& a)
: uc(ucl) 
{
  init(NOFLAG,r,d,a);
}


ARGS::ARGS(unionClass ucl,char* d,argsAction& a)
: uc(ucl) 
{
  init(NOFLAG,Opt,d,a);
}


ARGS::ARGS(const ARGS& a)
: uc(a.uc) {
  flag = a.flag;
  arg_kind = a.arg_kind;
  description = a.description;
  action = a.action;
}

void ARGS::makeFound() {
  if (found != NULL)
    return;
  numArgs = 0;
  ARGS* args_iter = Args;
  while (args_iter->flag != NULL) {
    numArgs ++;
    args_iter++;
  }
  found = new bool[numArgs];
  for (int i=0;i<numArgs;i++) {
    found[i] = false;
  }
}


ARGS::ArgsRetCode 
ARGS::parseFromCommandLine(int argc,char**argv)
{
  if (argv[0] != NULL)
    progName = argv[0];
  makeFound();
  ARGS* args_iter;
  for (int i=1;i<argc;i++) {
    if (argv[i][0] == '-') {
      char *flag = argv[i];
      args_iter = findMatch(Args,&flag[1]);
      if (args_iter == NULL) {
	if (ignoreUnknownSwitch == false) {
	  Error("Unknown switch: ",argv[i]);
	  return (ARG_ERROR);
	}
      } else if (args_iter == (ARGS*)(-1)) {
	Error("Ambiguous switch: ",argv[i]);
	return (ARG_ERROR);
      } else {
	char *arg;
	if ((i+1) >= argc)
	  arg = NULL;
	else
	  arg = argv[++i];
	argsSwitch(args_iter,arg,i,found[args_iter - Args],flag);
      }
    } else { // Go through and look for no flag case, in order.
      args_iter = Args;
      while (args_iter->flag != NULL) {
	if (args_iter->flag == NOFLAG) {
	  // assume string type
	  char *arg = argv[i];
	  int dummy;
	  argsSwitch(args_iter,arg,dummy,found[args_iter-Args],"");
	  args_iter->flag = NOFL_FOUND;
	  break;
	}
	args_iter++;
      }
    }
  }
  if (checkMissing())
    return ARG_MISSING;
  return ARG_OK;
}

bool ARGS::checkMissing(bool printMessage) {
  ARGS* args_iter;
  args_iter = Args;
  bool missing = false;
  while (args_iter->flag != NULL) {
    if (!found[args_iter-Args] && args_iter->arg_kind == Req) {
      if (printMessage) {
	ostrstream oss;
	if (noFlagP(args_iter->flag))
	  oss << "";
	else
	  oss << " -" << args_iter->flag;
	oss << " " << args_iter->uc.type;
	Error("Required argument missing:",oss.str());
      }
      missing = true;
    }
    args_iter++;
  }
  return missing;
}

// LIBG++ 2.2 dependant.
// THIS FUNCTION WILL NEED RE-WRITTING AS libg++ 2.2 fixes bugs.
// Still currently can't have lines with no characters in them. It messes
// up libg++ 2.2's internal state.
ARGS::ArgsRetCode ARGS::parseFromFile(char *fileName)
{
  makeFound();
  ifstream ifile(fileName);
  if (!ifile) {
    Error("Can't file argument file: ",fileName);
    return ARG_ERROR;
  } else {
    // get the arguments from the file.
    char buffer[512];
    char buffer_orig[512];
    char c;
    while (ifile.good()) { 
      ifile.get(buffer,512,'\n');
      // ****************************
      // ** LIBG++ VERSION 2.2 BUG WORK-ARROUND 
      // ****************************
      if (ifile.flags() | ios::failbit && ::strlen(buffer) == 0)
	ifile.unsetf(ios::failbit);
      // ****************************
      // ** END OF LIBG++ BUG WORK-ARROUND  
      // ****************************
      if (ifile.get(c) && c != '\n') {
	// input string is longer than expected.
	Error("Line length too long in file: ",fileName);
	return ARG_ERROR; // give up
      }
      ::strcpy(buffer_orig,buffer);
      char* buffp = buffer;
      while (*buffp) {
	if (*buffp == CommentChar)
	  *buffp = '\0';
	buffp++;
      }

      buffp = buffer;
      // skip space
      while (*buffp == ' ' || *buffp == '\t')
        buffp++;
      // get flag
      char *flag = buffp;
      char *arg;
      // get command up to space or ':'
      while (*buffp && *buffp != ' ' && *buffp != '\t' &&  *buffp != ':')
	buffp++;
      if (buffp == flag)
	continue; // empty line.
      
      if (*buffp) { 
	// get ':' and position buffp to start of arg.
	if (*buffp == ':')
	  *buffp++ = '\0'; // we have the flag
	else {
	  *buffp++ = '\0'; // we have the flag, but need to get rid of ':' if there.
	  while (*buffp && *buffp != ':')
	    buffp++;
	  if (*buffp == ':')
	    buffp++;
	}
	// skip space
	while (*buffp == ' ' || *buffp == '\t')
	  buffp++;
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
      ARGS* args_iter = findMatch(Args,flag);
      if (args_iter == NULL) {
	Error("Unknown switch in ",fileName,": ",flag);
	return (ARG_ERROR);
      } else if (args_iter == (ARGS*)(-1)) {
	Error("Ambiguous switch in ",fileName,": ",flag);
	return (ARG_ERROR);
      } else {
	// In a file, only set arg if it hasn't been set before. This
	// allows command line arguments to take precedence, but arguments
	// in files are only set once.
	if (found[args_iter-Args] != true) {
	  int i;
	  if (argsSwitch(args_iter,arg,i,found[args_iter - Args],flag) 
	      != ARG_OK) {
	    Error("Error in ",fileName);
	    return (ARG_ERROR);
	  }
	}
      }
    }
  }
  if (checkMissing())
    return ARG_MISSING;
  return ARG_OK;
}



// Return NULL for unknown switch, return (ARGS*)(-1) for ambiguous switch,
// otherwise returns the ARGS* that matches.
ARGS* ARGS::findMatch(ARGS* ag,char *flag) {
  int flaglen = ::strlen(flag);
  ARGS* args_iter = ag;
  int numTaged = 0;  
  int lastTaged = -1;

  // find the one that best matches.
  while (args_iter->flag != NULL) {
    if (!noFlagP(args_iter->flag)) {
      if (!::strncmp(args_iter->flag,flag,flaglen)) {
	numTaged++;
	lastTaged = (args_iter - ag);
      }
    }
    args_iter++;
  }
  if (numTaged == 0)
    return NULL;
  else if (numTaged > 1)
    return (ARGS*)(-1);
  else 
    return &ag[lastTaged];
}


// only returns ARGS::ARG_OK or ARG_ERROR
ARGS::ArgsRetCode 
ARGS::argsSwitch(ARGS* args_iter,char *arg,int& index,bool& found,char*flag)
{
  switch(args_iter->uc.type) {
  // ---------------------------------------------------------------------
  case unionClass::bool_t: {
    bool b; 
    if (arg == NULL || // end of arguments.
	arg[0] == '-') { // for optionless flag, just turn on.
      // for bool case, just turn it on.
      if (args_iter->arg_kind != Tog)
	args_iter->uc.ptr->boolean = true;
      else
	args_iter->uc.ptr->boolean = 
	  ((args_iter->uc.ptr->boolean == true) ? false : true);
      if (arg != NULL)
	index--;
      found = true;
      if (args_iter->action != NULL)
	args_iter->action->act(args_iter->uc.ptr->boolean);
      break;	      
    }

    // If kind  is Tog and we have a valid boolean 
    // argument, then treat argument as normal explicit boolean argument.
    if (!boolable(arg,b)) {
      Error("Boolean argument needed: ",flag," ",arg);
      return ARG_ERROR;
    }
    args_iter->uc.ptr->boolean = b;
    found = true;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->boolean);
  }
    break;
  // ---------------------------------------------------------------------
  case unionClass::char_t: {
    if (arg == NULL || // end of arguments.
	::strlen(arg) != 1) {
      Error("Character argument needed: ",flag," ",arg);
      return ARG_ERROR;
    }
    args_iter->uc.ptr->ch = arg[0];
    found = true;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->ch);
  }
    break;
  // ---------------------------------------------------------------------
  case unionClass::str_t: {
    if (arg == NULL) { // end of arguments.
      Error("String argument needed: ",flag);
      return ARG_ERROR;
    }
    args_iter->uc.ptr->string =
      ::strcpy(new char[strlen(arg)+1],arg);
    found = true;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->string);
  }
    break;
  // ---------------------------------------------------------------------
  case unionClass::int_t: {
    if (arg == NULL) { // end of arguments.
      Error("Integer argument needed: ",flag);
      return ARG_ERROR;
    }
    int n;
    istrstream ist(arg);
    if (!(ist >> n)) {
      Error("Integer argument needed: ",flag," ",arg);
      return ARG_ERROR;
    }
    args_iter->uc.ptr->integer = n;
    found = true;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->integer);
  }
    break;
  // ---------------------------------------------------------------------
  case unionClass::uint_t: {
    if (arg == NULL) { // end of arguments.
      Error("Unsigned integer argument needed: ",flag);
      return ARG_ERROR;
    }
    unsigned int n;
    istrstream ist(arg);
    if (!(ist >> n)) {
      Error("Unsigned integer argument needed: ",flag," ",arg);
      return ARG_ERROR;
    }
    args_iter->uc.ptr->uinteger = n;
    found = true;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->uinteger);
  }
    break;
  // ---------------------------------------------------------------------
  case unionClass::float_t: {
    if (arg == NULL) { // end of arguments.
      Error("Real number argument needed: ",flag);
      return ARG_ERROR;
    }
    float f;
    istrstream ist(arg);
    if (!(ist >> f)) {
      Error("Floating point number argument needed: ",flag," ",arg);
      return ARG_ERROR;
    }
    args_iter->uc.ptr->single_prec = f;
    found = true;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->single_prec);
  }
    break;
  // ---------------------------------------------------------------------
  case unionClass::double_t: {
    if (arg == NULL) { // end of arguments.
      Error("Real number argument needed: ",flag);
      return ARG_ERROR;
    }
    double d;
    istrstream ist(arg);
    if (!(ist >> d)) {
      Error("Integer argument needed: ",flag," ",arg); 
      return ARG_ERROR;
    }
    args_iter->uc.ptr->double_prec = d;
    found = true;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->double_prec);
  }
    break;
  // ---------------------------------------------------------------------
  case unionClass::action_t:
    // the argsAction* is stored in ptr itself.
    ((argsAction*)args_iter->uc.ptr)->act();
    if (arg != NULL)
      index--;
    found = true;
    if (args_iter->action != NULL)
      args_iter->action->act();
    break;
  // ---------------------------------------------------------------------
  default:
    Error("Unknown internal argument: Internal error.");
    ::exit(1);
    break;
  }
  return ARG_OK;
}

bool ARGS::boolable(char *string,bool& value)
{
  bool rc;
  int arglen = strlen(string);
  char *upcasearg = new char[arglen+1];
  ::strcpy(upcasearg,string);
  for (int i=0;i<arglen;i++)
    upcasearg[i] = toupper(upcasearg[i]);
  if (!::strncmp("TRUE",upcasearg,arglen) ||
      !::strncmp("YES",upcasearg,arglen)  ||
      !::strncmp("ON",upcasearg,arglen)  ||
      (arglen == 1 && string[0] == '1')) {
    value = true;
    rc = true;
  } 
  else if (!::strncmp("FALSE",upcasearg,arglen) || 
	   !::strncmp("NO",upcasearg,arglen) || 
	   !::strncmp("OFF",upcasearg,arglen) || 
	   (arglen == 1 && string[0] == '0')) {
    value = false;
    rc = true;
  } else {
    rc = false;
  }
  delete [] upcasearg;
  return rc;
}

// return true if there is no flag.
bool ARGS::noFlagP(char *flg) {
  // cant test here equality with NOFLAG, so do the following instread.
  if (flg == NOFLAG || flg == NOFL_FOUND)
    return true;
  else
    return false;
}




void ARGS::Error(char* s1,char*s2,char*s3,char*s4,char*s5)
{
  cerr << "Argument error: " << s1;
  if (s2 != NULL)
    cerr << s2;
  if (s3 != NULL)
    cerr << s3;
  if (s4 != NULL)
    cerr << s4;
  if (s5 != NULL)
    cerr << s5;
  cerr << "\n";
}

void ARGS::usage(bool force) {

  if (usage_called && !force)
    return;
  usage_called = true;

  cerr << "Usage: " << progName << " [[[-flag] [option]] ...]\n";
  cerr << "Required: <>; Optional: []; Flagless arguments must be in order.\n";

  ARGS* args_iter = Args;
  int longest_variation = 0;

  while (args_iter->flag != NULL) {
    int len = 0;
    if (!noFlagP(args_iter->flag)) {
      // add one for the '-', as in "-flag"
      len += ::strlen(args_iter->flag)+1;
      if (args_iter->uc.type != unionClass::action_t)
	len ++; //  add one for the ' ' in "-flag "
    }
    len += ::strlen(unionClass::typeStr(args_iter->uc.type));
    if (args_iter->uc.type != unionClass::action_t)
      len += 2; // add two for brackets. '[',']', or '<','>' around type.
    if (len  > longest_variation)
      longest_variation = len;
    args_iter++;
  }
  
  args_iter = Args;
  while (args_iter->flag != NULL) {
    int this_variation = 0;
    char brackets[2];
    if (args_iter->arg_kind == Req) {
      brackets[0] = '<'; brackets[1] = '>';
    } else {
      brackets[0] = '['; brackets[1] = ']';
    }
    cerr << "  " << brackets[0]; 

    if (!noFlagP(args_iter->flag)) {
      cerr << "-" << args_iter->flag;
      // add one for the '-', as in "-flag"
      this_variation = ::strlen(args_iter->flag) + 1;
      if (args_iter->uc.type != unionClass::action_t) {
	cerr << " ";
	this_variation ++; //  add one for the ' ' in "-flag "
      }
    }
    cerr << args_iter->uc.type;
    this_variation += ::strlen(unionClass::typeStr(args_iter->uc.type));
    if (args_iter->uc.type != unionClass::action_t) {
      // add two for brackets. '[',']', or '<','>' around type.
      this_variation += 2;
    }

    cerr << brackets[1];
    while (this_variation++ < longest_variation)
      cerr << " ";
    cerr << "     ";
    cerr << ((args_iter->description == NULL) ? "" : args_iter->description) 
	     << '\n';
    args_iter++;
  }  
}

ostream& operator << (ostream& os,unionClass::argType& t) {
  if (t == unionClass::action_t)
    return os;
  char brackets[2];
  if (t == unionClass::bool_t) {  // bool_t arg is always optional. (i.e. -b T, -b F, or -b)
    brackets[0] = '['; 
    brackets[1] = ']';
  } else {
    brackets[0] = '<'; 
    brackets[1] = '>';
  }
  os << brackets[0] << unionClass::typeStr(t) << brackets[1];
  return os;
}


char *unionClass::typeStr(unionClass::argType at) {
  switch (at) {
  case unionClass::bool_t:
    return "bool";
  case unionClass::char_t:
    return "char";
  case unionClass::str_t:
    return "string";
  case unionClass::int_t:
    return "integer";
  case unionClass::float_t:
    return "float";
  case unionClass::double_t:
    return "double";
  case unionClass::action_t:
    return ""; // "ACTION";
  default:
    return "error";
  }
}

ostream& operator << (ostream&os,ARGS*args) { 
  ARGS *args_iter = args;
  while (args_iter->flag != NULL) {
    os << *args_iter;
    args_iter ++;
  }
  return os; 
}

ostream& operator << (ostream&os,ARGS& arg) { 
  if (!ARGS::noFlagP(arg.flag))
    os << arg.flag;
  os << ':';
  os << unionClass::typeStr(arg.uc.type) << " = ";
  switch (arg.uc.type) {
  case unionClass::bool_t:
    os << arg.uc.ptr->boolean; break;
  case unionClass::char_t:
    os << arg.uc.ptr->ch; break;
  case unionClass::str_t: 
    os << (arg.uc.ptr->string == NULL ? "" : arg.uc.ptr->string); break;
  case unionClass::int_t:
    os << arg.uc.ptr->integer; break;
  case unionClass::float_t:
    os << arg.uc.ptr->single_prec; break;
  case unionClass::double_t:
    os << arg.uc.ptr->double_prec; break;
  case unionClass::action_t:
    break;
  default:
    break;
  }
  os << ";  // " << ((arg.description == NULL) ? "" : arg.description) 
     << '\n';
  return os;
}


// A typical way one might want to parse their arguments.
// exit the program if any Required arguments are missing or if there is
// some other kind of error.
void ARGS::parse(int argc,char** argv,char* &argsFilep)
{
  ArgsRetCode rc;
  rc = parseFromCommandLine(argc,argv);
  if (argsFilep != NULL)
    rc = parseFromFile(argsFilep);
  if (rc == ARG_MISSING) { // still missing ??
    ARGS::printMissing();
  }
  if (rc != ARG_OK) {
    ARGS::usage();
    ::exit(1);
  }
}


argsAction::ActionRetCode
argsAction::satisfy(char*flag) {
  ARGS* args_iter;
  args_iter = ARGS::findMatch(ARGS::Args,flag);
  if (args_iter == NULL || args_iter == (ARGS*)(-1))
    return ERROR;
  ARGS::found[args_iter - ARGS::Args] = true;
  return OK;
}

argsAction::ActionRetCode
argsAction::setArg(char* flag,unionClass::argType t,unionClass::argU au) 
{
  ARGS* args_iter;
  args_iter = ARGS::findMatch(ARGS::Args,flag);
  if (args_iter == NULL || args_iter == (ARGS*)(-1))
    return ERROR;
  if (args_iter->uc.type != t)
    return ERROR;
  ARGS::found[args_iter - ARGS::Args] = true;
  switch(t) {
  case unionClass::bool_t:
    args_iter->uc.ptr->boolean = au.boolean;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->boolean);
    break;
  case unionClass::char_t:
    args_iter->uc.ptr->ch = au.ch;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->ch);
    break;
  case unionClass::str_t:
    args_iter->uc.ptr->string = au.string;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->string);
    break;
  case unionClass::int_t:
    args_iter->uc.ptr->integer = au.integer;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->integer);
    break;
  case unionClass::float_t:
    args_iter->uc.ptr->single_prec = au.single_prec;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->single_prec);
    break;
  case unionClass::double_t:
    args_iter->uc.ptr->double_prec = au.double_prec;
    if (args_iter->action != NULL)
      args_iter->action->act(args_iter->uc.ptr->double_prec);
    break;
  case unionClass::action_t:
    break;
  default:
    break;
  }
  return OK;
}

argsAction::ActionRetCode
argsAction::setArg(char* flag,bool b) {
  unionClass::argU au;
  au.boolean = b;
  return setArg(flag,unionClass::bool_t,au);
}


argsAction::ActionRetCode
argsAction::setArg(char* flag,char c) {
  unionClass::argU au;
  au.ch = c;
  return setArg(flag,unionClass::char_t,au);
}


argsAction::ActionRetCode
argsAction::setArg(char* flag,char* str) {
  unionClass::argU au;
  au.string = str;
  return setArg(flag,unionClass::str_t,au);
}


argsAction::ActionRetCode
argsAction::setArg(char* flag,int i) {
  unionClass::argU au;
  au.integer = i;
  return setArg(flag,unionClass::int_t,au);
}


argsAction::ActionRetCode
argsAction::setArg(char* flag,float f) {
  unionClass::argU au;
  au.single_prec = f;
  return setArg(flag,unionClass::float_t,au);
}


argsAction::ActionRetCode
argsAction::setArg(char* flag,double d) {
  unionClass::argU au;
  au.double_prec = d;
  return setArg(flag,unionClass::double_t,au);
}



#ifdef MAIN

class myAction : public argsAction {
  char *str;
public:
  myAction(char *a) { str = a; }
  void act() { cerr << "In void myAction::act() :" << str << '\n'; }
  void act(bool& b) { cerr << "In void myAction::act(bool&), b = " << b
		      << ": " << str << '\n'; }
  void act(float& f) { cerr << "In void myAction::act(float&), f = " << f
			    << ": " << str << '\n';
		       setArg("doub",(double)f);
		     }

};


int bufSize = 3;
int buffSize = 3;
char *bla = "BARSTR";
bool truth = false;
char chr = 'c';
float flt = 2.3;
double dob = .4;
bool toggle1 = false;
bool toggle2 = true;

char *nfstr="OSTRTesting";
int nfint = 5;
bool nfbool = false;
char nfchar = 'C';
float nffloat = 3.4;
double nfdouble = 4.5;

char *argsFile = NULL;

myAction ma1("toggle 2 action");
myAction ma2("action 1 action");
myAction ma3("set doub action");

ARGS ARGS::Args[] = {

 ARGS("buffSize",ARGS::Req,buffSize,"The buffer size"),
 ARGS("bufSize",ARGS::Req,bufSize,"The bufer size"),
 ARGS("bla",ARGS::Opt,bla,"The bla"),
 ARGS("truth",ARGS::Opt,truth,"A truth value"),
 ARGS("chars",ARGS::Opt,chr,"A char"),
 ARGS("sing",ARGS::Opt,flt,"A single prec",ma3),
 ARGS("doub",ARGS::Req,dob,"A double prec"),
 ARGS("toggle1",ARGS::Tog,toggle1,"toggle1"),
 ARGS("toggle2",ARGS::Tog,toggle2,"toggle2",ma1),
 ARGS("action",ARGS::Opt,ma2,"action 1"),
 ARGS("argsFile",ARGS::Opt,argsFile,"Arguments file"),


 // Flagless (No-flag) arguments. The order on the command line must be in the
 //   the same as the order given here. i.e. string first, then int, then
 //   bool, etc.
 ARGS(ARGS::Req,nfstr,"No flag string"),
 ARGS(ARGS::Opt,nfint,"No flag int"),
 ARGS(ARGS::Opt,nfbool,"No flag bool"),
 ARGS(ARGS::Opt,nfchar,"No flag char"),
 ARGS(ARGS::Opt,nffloat,"No flag float"),
 ARGS(ARGS::Opt,nfdouble,"No flag double"),


 // the following END marker must be present.
 // ARGS::END
 ARGS()

};


int
main(int argc,char*argv[])
{
  ARGS::parse(argc,argv,argsFile);
  cout << ARGS::Args;
}

/* 
** Easy and small template to copy from:
**

int buffSize = 10;
ARGS ARGS::Args[] = {
 ARGS("buffSize",ARGS::Req,buffSize,"The buffer size"),
 ARGS::END
};

*/

#endif // def MAIN
