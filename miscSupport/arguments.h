// *************************************************************************
// * A general mechanism to parse command line and file arguments for 
// * C++ programs.
// * <bilmes@media.mit.edu> Aug 1992
// *
// *************************************************************************
//

// TODO:
//    1. Usage message should print required arguments first.  
//    2. Parse from file should have capability to read multi word string args.
//    3. Make error reports better.
//    4. Fix bug with empty line in parse args from file (libg++ 2.2).
//    5. negative action, action called if argument is not present.
//    6. Add ability to specify controls in documentation string.
//       i.e. to we can do:  ... "Parm X, default = %d",default);

#ifndef ARGS_h
#define ARGS_h

// ... OR ...
// The following bool definition will work fine.
#include <iostream.h>
// enum bool {true = 1, false = 0};
// inline ostream& operator << (ostream&os,bool b) {
//   return (os << ((b == true) ? "T" : "F"));
// }

class ARGS;
class argsAction;
extern ostream& operator << (ostream&os,ARGS*);
extern ostream& operator << (ostream&os,ARGS&);


class unionClass {
  friend class ARGS;
  friend class argsAction;
  friend ostream& operator << (ostream&,ARGS*);
  friend ostream& operator << (ostream&,ARGS&);
  enum argType { bool_t, char_t, str_t, int_t, uint_t, float_t, 
		   double_t, action_t };
  union argU {
    bool boolean;
    char ch;
    char * string;
    int integer;
    unsigned int uinteger;
    float single_prec;
    double double_prec;
    argsAction* action;  // dummy field.
  };
  argU* ptr;
  argType type;
  static char* typeStr(unionClass::argType);
  friend ostream& operator << (ostream&,unionClass::argType&);
 public:
  unionClass(bool& b) { ptr = (argU*)&b; type = bool_t; }
  unionClass(char& c) { ptr = (argU*)&c; type = char_t; }
  unionClass(char*& s) { ptr = (argU*)&s; type = str_t; }
  unionClass(int& i)  { ptr = (argU*)&i; type = int_t; }
  unionClass(unsigned int& i)  { ptr = (argU*)&i; type = uint_t; }
  unionClass(float& f) { ptr = (argU*)&f; type = float_t; }
  unionClass(double& d) { ptr = (argU*)&d; type = double_t; }
  unionClass(argsAction& a) { ptr = (argU*)&a; type = action_t; }
};

extern ostream& operator << (ostream&,unionClass::argType&);

class argsAction {
  enum ActionRetCode { OK, ERROR };
 private:
  ActionRetCode setArg(char*,unionClass::argType,unionClass::argU);
 protected:
  ActionRetCode setArg(char*,bool);
  ActionRetCode setArg(char*,char);
  ActionRetCode setArg(char*,char*);
  ActionRetCode setArg(char*,int);
  ActionRetCode setArg(char*,float);
  ActionRetCode setArg(char*,double);
  ActionRetCode satisfy(char*);
 public:
  virtual void act() {}
  virtual void act(bool& b) { b=b;}
  virtual void act(char& c) { c=c; }
  virtual void act(char*& str) { str=str; }
  virtual void act(int& i) { i=i; }
  virtual void act(unsigned int& ui) { ui=ui; }
  virtual void act(float& f) { f=f; }
  virtual void act(double& d) { d=d; }
};



class ARGS {
  friend class unionClass;
  friend class argsAction;
 public:

  enum argKind { Opt, Req, Tog }; // optional, required, or toggle
  enum ArgsRetCode { ARG_MISSING, ARG_OK, ARG_ERROR };
  static ARGS Args[];
  static ARGS END;

 private:

  static char* const NOFLAG;
  static char* const NOFL_FOUND;
  static const char CommentChar;  // for argument files.
  static bool usage_called;
  static int numArgs;
  static bool* found;
  static char* progName;
  static bool ignoreUnknownSwitch;

  // instance data members.
  char *flag; // name to match on command line, NULL when end.
  argKind arg_kind;  // optional, required, toggle
  unionClass uc;
  char *description;
  argsAction* action;


  static argsAction no_action;
  static bool __args__dummy__;

 public:
  ARGS(char*,argKind,unionClass,char*d=NULL,argsAction&a=no_action);
  ARGS(argKind,unionClass,char*d=NULL,argsAction&a=no_action);
  ARGS(unionClass,char*d=NULL,argsAction&a=no_action);
  ARGS(const ARGS&);
  ARGS();

  static ArgsRetCode parseFromCommandLine(int,char**);
  static ArgsRetCode parseFromFile(char*f="argsFile");
  static void parse(int,char**,char*&);
  static void parse(int i,char**c) { char*d=NULL; parse(i,c,d); }
  static void printMissing() { checkMissing(true); }
  static void usage(bool force=false);
  static void ignoreUnknownFlag() { ignoreUnknownSwitch = true; }
  static void regardUnknownFlag() { ignoreUnknownSwitch = false; }

 private:
  void init(char*,argKind,char*,argsAction&);
  static void Error(char*,char*s2=NULL,char*s3=NULL,char*s4=NULL,
		    char*s5=NULL);

  static bool noFlagP(char *);
  static ArgsRetCode argsSwitch(ARGS*,char *,int&,bool&,char*);
  static ARGS* findMatch(ARGS*,char*);
  static void makeFound();
  static bool checkMissing(bool printMessage=false);
  static bool boolable(char* string,bool&value);
  friend ostream& operator << (ostream&,ARGS*);
  friend ostream& operator << (ostream&,ARGS&);

};

#endif
