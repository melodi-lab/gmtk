//
// General miscellaneous stuff that belongs nowhere else.
//
// $Header$
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu


#ifndef GENERAL_H
#define GENERAL_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "machine-dependent.h"

#define VCID(x) static char * version_control_id = x; static char *___tmp___ = version_control_id;

char *copyToNewStr(const char *const str);

bool strIsInt(const char *const str, int* i=NULL,int* len=NULL);


// a general swapping routine.
template <class T>
void genSwap(T& v1, T& v2) 
{
   T tmp = v1;
   v1 = v2;
   v2 = tmp;
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
void copyStringWithTag(char *result,char *input,
		       const int tag, const int maxLen);


unsigned long fsize(FILE*stream);
unsigned long fsize(const char* const filename);




#endif
