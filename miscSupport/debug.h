/*
    $Header$
  
    Simple verbosity,informational,debugging error function.

    Jeff Bilmes <bilmes@cs.berkeley.edu>
*/



#ifndef DEBUG_H
#define DEBUG_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

// set this to zero to make all the routines nop stubs.
#define INFO_MESSAGES_ON 1

class IM {

  static unsigned globalMessageLevel;
  unsigned messageLevel;
  bool messagesOn;

public:

  IM() {
    messageLevel = globalMessageLevel;
    messagesOn = true;
  }

  // levels for general informational and/or debugging messages.

  enum VerbosityLevels { 
    Min = 0,
    Warning = 0, // Anything lower than this should be an error
                 // called with the error() function which will kill
                 // the program. Also, this level is the
                 // threshold at which things get sent to stderr.
                 //
    Nano = 10,   // This should be the default.
    Default = 10,
    Tiny = 20,
    Low  = 30,
    Med  = 40,
    High = 50,
    Huge = 60,
    Mega = 70,
    Rediculous = 80,
    Max  = 80    //  keep this set at the maximum
  };


  inline bool message(unsigned v) {
    return (messagesOn && ((v <= messageLevel) || (v <= globalMessageLevel)));
  }

  static inline bool messageGlb(unsigned v) {
    return (v <= globalMessageLevel);
  }

  inline void infoMsg(unsigned v,char*format, ...) 
  {
#if INFO_MESSAGES_ON
    if (message(v)) {
      va_list ap;
      va_start(ap,format);
      if (v == Warning)
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
    if (message(Default)) {
      va_list ap;
      va_start(ap,format);
      if (Default == Warning)
	(void) vfprintf(stderr, format, ap);
      else 
	(void) vfprintf(stdout, format, ap);
      va_end(ap);
    }
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



#endif
