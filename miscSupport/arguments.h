/*-
 *-----------------------------------------------------------------------
 * A simple library to parse command lines in C++ using
 * an easy interface to quickly define arguments.
 * 
 *       Jeff Bilmes <bilmes@ee.washington.edu>
 *
 *  $Header$
 *
 *-----------------------------------------------------------------------
 */


#ifndef ARGS_h
#define ARGS_h

#include <iostream>
#include <stdio.h>

class MultiType {
  friend class Arg;
  enum ArgumentType { 
    int_type,    // integer
    uint_type,   // unsigned integer
    float_type,  // single precision float
    double_type, // double precision float
    str_type,    // string precision float
    char_type,   // char precision float
    bool_type   // boolean precision float
  };
  union MultiArg {
    int integer;
    unsigned int uinteger;
    float single_prec;
    double double_prec;
    char * string;
    char ch;
    bool boolean;
  };
  MultiArg* ptr;
  ArgumentType type;
  static char* printable(MultiType::ArgumentType);
 public:
  MultiType(bool& b) { ptr = (MultiArg*)&b; type = bool_type; }
  MultiType(char& c) { ptr = (MultiArg*)&c; type = char_type; }
  MultiType(char*& s) { ptr = (MultiArg*)&s; type = str_type; }
  MultiType(int& i)  { ptr = (MultiArg*)&i; type = int_type; }
  MultiType(unsigned int& i)  { ptr = (MultiArg*)&i; type = uint_type; }
  MultiType(float& f) { ptr = (MultiArg*)&f; type = float_type; }
  MultiType(double& d) { ptr = (MultiArg*)&d; type = double_type; }
  void print(FILE*);
};


class Arg {
  friend class MultiType;
 public:

  // the argument disposition, optional, required, or toggle
  enum ArgDisposition { Opt, Req, Tog };
  // the return codes, missing, ok, or in error.
  enum ArgsRetCode { ARG_MISSING, ARG_OK, ARG_ERROR };
  // the actual argument array itself.
  static Arg Args[];

 private:

  static char* const NOFLAG;
  static char* const NOFL_FOUND;
  static const char COMMENTCHAR;  // for argument files.
  // the total number of arguments in a given program.
  static int Num_Arguments;
  // an array of bits set if a particular argument is used
  static bool* Argument_Specified;
  // the program name, saved for usage messages.
  static char* Program_Name;

  // instance data members.
  char *flag; // name to match on command line, NULL when end.
  ArgDisposition arg_kind;  // optional, required, toggle
  MultiType mt;
  char *description;
  
  // special lvalue boolean to be able to construct empty argument
  // entry for end of array.
  static bool EMPTY_ARGS_FLAG;

 public:
  Arg(char*,ArgDisposition,MultiType,char*d=NULL);
  Arg(ArgDisposition,MultiType,char*d=NULL);
  Arg(const Arg&);
  Arg();
  ~Arg();

  static ArgsRetCode parseArgsFromCommandLine(int,char**);
  static ArgsRetCode parseArgsFromFile(char*f="argsFile");
  static void parse(int i,char**c);
  static void usage();
  static void printArgs(Arg*args,FILE*f);

 private:
  void initialize(char*,ArgDisposition,char*);

  static bool noFlagP(char *);
  static ArgsRetCode argsSwitch(Arg*,char *,int&,bool&,char*);
  static Arg* searchArgs(Arg*,char*);
  static void countAndClearArgBits();
  static bool checkMissing(bool printMessage=false);
  static bool validBoolean(char* string,bool&value);
  void print(FILE*);

};

#endif
