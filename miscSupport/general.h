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
#include <string>

#include "machine-dependent.h"

#define VCID(x) static const char * const version_control_id = x;

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


int stringprintf(string& str,const char *format, ...);


/*
 * returns number of bits required to represent
 * a value v between [0 <= v < val], where val > 0..
 */
unsigned bitsRequiredUptoNotIncluding(unsigned val); 


/*
 * returns integer 2^ceil(log2(val)), least
 * power of two that is >= val. I.e., returns
 * 2^k* where k* = argmin_k { k : 2^k >= val }
 */
unsigned nextPower2(unsigned val); 

/*
 * returns ceil(log2(val))
 */
unsigned ceilLog2(unsigned val);

/*
 * returns the number of bits set in the unsigned
 *
 */
unsigned int numBitsSet(unsigned u);


/*
 * returns a char* to the cpp command to use.
 * It returns a string that should not be freed.
 */
const char* CPP_Command();

bool freadUntilEOF(FILE *f);

/*
 * MACHINE DEPENDENT typedefs.
 */

// 8-bit (1 byte) unsigned integer type
typedef unsigned char UInt8;

// 8-bit (1 byte) signed integer type
typedef char Int8;

// 16-bit (2 byte) unsigned integer type
typedef unsigned short UInt16;

// 16-bit (2 byte) signed integer type
typedef int short Int16;

// 32-bit (4 byte) unsigned integer type.
typedef unsigned int UInt32;

// 32-bit (4 byte) signed integer type.
typedef int Int32;

// 64-bit (8 byte) unsigned integer type.
// typedef long long unsigned UInt64;

// 64-bit (8 byte) signed integer type.
// typedef long long int Int64;


#endif
