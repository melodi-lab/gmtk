/*
    $Header$
  
    Simple verbosity,informational,debugging error function.
    There are two ways to use this facility:

    1) as a class-specific debug message thing, where you use inheritance and do:
          class FOO : public IM { ... }
       In this case, the debug levels are specific to the class, and info messages
       can be called as infoMesg(Low,"message");

    2) A global info message handler, just call infoMsg(IM::Low,"message");
       (see below).

    TODO: change the name of this from debug.h to infoMesg.h

    Jeff Bilmes <bilmes@ee.washington.edu>
*/



#ifndef DEBUG_H
#define DEBUG_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

// set this to zero to make all the routines nop stubs.
#define INFO_MESSAGES_ON 1

class IM {
  friend class InferenceMaxClique;
  static unsigned globalMessageLevel;
  unsigned messageLevel;
  bool messagesOn;
  bool flush;

public:

  IM() {
    messageLevel = globalMessageLevel;
    messagesOn = true;
    flush = true;
  }

  // levels for general informational and/or debugging messages.
  // In general, there are levels between 0 and 100, where
  // larger means be more verbose, and lower means be less
  // verbose. A info message gives a tag that says
  // how verbose it is (i.e., a tag of 0 means that this is a
  // warning that could indicate a problem, a tag of 100 means
  // that this is something that should not exist.
  enum VerbosityLevels { 
    Nothing = 0, // At this level nothing is intended to be printed (note
                 // that user should never do infoMsg(Nothing,"") for this to work).
    Min = 1,     // Min level at which anything is printed.
    Warning = 1, // Anything lower than this should be an error
                 // called with the error() function which will kill
                 // the program. Also, this level is the
                 // threshold at which things get sent to stderr.
    SoftWarning = 6, // warnings that aren't too bad.
    Nano = 10,   // This should be the default.
    Default = 10, // Another name for default. Prints using default command line option.
    Info = 11, // Informative but not crucially important commands.
    Tiny = 20,
    Low  = 30,
    Moderate = 40,
    Mod  = 40, // alias for Moderate
    Med  = 50,
    High = 60,
    Huge = 70,
    Mega = 80,
    Giga = 90,
    Max  = 100,    //  Maximum setable level.

    // We also define an "increment" variable here, so that a user
    // can do something like 
    //   infoMsg( (Nano + 0.3*Increment), "bar" );
    Increment = 10
    // Note also that the user could do something
    // like 
    //   infoMsg( (Tiny + Low)/2, "baz" );
    // to get something between tiny and low.
  };


  inline bool message(unsigned v) {
#if INFO_MESSAGES_ON
    return (messagesOn && ((v <= messageLevel) || (v <= globalMessageLevel)));
#else 
    return false;
#endif
  }

  static inline bool messageGlb(unsigned v) {
#if INFO_MESSAGES_ON
    return (v <= globalMessageLevel);
#else
    return false;
#endif    
  }

  inline void infoMsg(unsigned v,char*format, ...) 
  {
#if INFO_MESSAGES_ON
    if (message(v)) {
      va_list ap;
      va_start(ap,format);
      if (v == Warning)
	(void) vfprintf(stderr, format, ap);
      else {
	(void) vfprintf(stdout, format, ap);
	if (flush) fflush(stdout);
      }
      va_end(ap);
    }
#endif
  }

  inline void infoMsg(char*format, ...) 
  {
#if INFO_MESSAGES_ON
    if (message(Default)) {
      va_list ap;
      va_start(ap,format);
      if (Default == Warning)
	(void) vfprintf(stderr, format, ap);
      else {
	(void) vfprintf(stdout, format, ap);
	if (flush) fflush(stdout);
      }
      va_end(ap);
    }
#endif
  }

  // already checked/force message case.
  // I.e., this routine can be used with an
  // external check, so we don't pay to do a
  // routine call case this doesn't get inlined.
  inline void infoMsgForce(char*format, ...) 
  {
#if INFO_MESSAGES_ON
    va_list ap;
    va_start(ap,format);
    if (Default == Warning)
      (void) vfprintf(stderr, format, ap);
    else {
      (void) vfprintf(stdout, format, ap);
      if (flush) fflush(stdout);
    }
    va_end(ap);
#endif
  }


  void msgsOn() { messagesOn = true; }
  void msgsOff() { messagesOn = false; }
  unsigned msgLevel() { return messageLevel; }
  unsigned setMsgLevel(const unsigned ml) { messageLevel = ml; return ml; }
  static unsigned glbMsgLevel() { return globalMessageLevel; }
  static unsigned setGlbMsgLevel(const unsigned ml) { 
    globalMessageLevel = ml; 
    return ml; 
  }

};

// Global version
inline void infoMsg(unsigned v,char*format, ...) 
{
#if INFO_MESSAGES_ON
  if (IM::messageGlb(v)) {
    va_list ap;
    va_start(ap,format);
    if (v == IM::Warning)
      (void) vfprintf(stderr, format, ap);
    else 
      (void) vfprintf(stdout, format, ap);
    va_end(ap);
  }
#endif
}

inline void infoMsg(char*format, ...) 
{
#if INFO_MESSAGES_ON
  if (IM::messageGlb(IM::Default)) {
    va_list ap;
    va_start(ap,format);
    if (IM::Default == IM::Warning)
      (void) vfprintf(stderr, format, ap);
    else 
      (void) vfprintf(stdout, format, ap);
    va_end(ap);
  }
#endif
}


inline void infoMsgForce(char*format, ...) 
{
#if INFO_MESSAGES_ON
    va_list ap;
    va_start(ap,format);
    if (IM::Default == IM::Warning)
      (void) vfprintf(stderr, format, ap);
    else 
      (void) vfprintf(stdout, format, ap);
    va_end(ap);
#endif
}



#endif
