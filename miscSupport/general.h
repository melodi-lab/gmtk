//
// General miscellaneous stuff that belongs nowhere else.
//
// $Header$
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu


#ifndef GENERAL_H
#define GENERAL_H

using namespace std;

// added this to keep gcc -pedantic flag from complaining
// about incompatibilities with not defining certain
// extern "C" functions with the 'throw()' directive
// in the function prototypes. Such functions include
// popen, pclose, drand48, seed48, etc.
#ifndef __THROW
#define __THROW
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>

#include <vector>

#include "machine-dependent.h"

#define VCID(x) static char * version_control_id = x; static char *___tmp___ = version_control_id;

char *copyToNewStr(const char *const str);

bool strIsInt(const char *const str, int* i=NULL,int* len=NULL);
bool strIsInt(const char *const str, unsigned* i=NULL,int* len=NULL);


// a general swapping routine.
template <class T>
void genSwap(T& v1, T& v2) 
{
   T tmp = v1;
   v1 = v2;
   v2 = tmp;
}

// general max/min class
template <class T>
T& max(T& v1, T& v2)
{
  if (v1 > v2)
    return v1;
  else
    return v2;
}

template <class T>
T& min(T& v1, T& v2)
{
  if (v1 < v2)
    return v1;
  else
    return v2;
}


template <class T>
void deleteObsInVector(vector < T* >& v) {
  for (unsigned i=0;i<v.size();i++) 
    delete v[i];
}


#define CSWT_EMPTY_TAG (~0)
// Copies input over to result and if 
// input has any occurence of '%d' in it, replace it with
// the string version of the integer tag. Assume
// result is plenty big.
void copyStringWithTag(char *result,const char *const input,
		       const int tag, const int maxLen);


unsigned long fsize(FILE*stream);
unsigned long fsize(const char* const filename);


void print_date_string(FILE* f);

void exit_program_with_status(const int stat);

void memory_error();

// return the log10 add of v1 and v2,i.e., res = log10(10^v1 + 10^v2)
double log10add(double v1,double v2);

// report timing for getrusage
void reportTiming(// input 
		  const struct rusage& rus,
		  const struct rusage& rue,
		  // output
		  double& userTime, 
		  double& sysTime,
		  // input
		  FILE* outputf=NULL);


int stringprintf(string& str,char *format, ...);


/*
 * returns number of bits required to represent
 * a value v between [0 <= v < val], where val > 0..
 */
unsigned bitsRequiredUptoNotIncluding(unsigned val); 

#endif
