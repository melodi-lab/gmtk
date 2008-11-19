//
// General miscellaneous stuff that belongs nowhere else.
// 
// Written by: Jeff Bilmes
//             bilmes@ee.wasington.edu

#include <ctype.h>
#include <sys/time.h>
#include <time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <string.h>

#include "general.h"
#include "rand.h"
#include "error.h"


VCID("$Header$")

char *copyToNewStr(const char *const str)
{
  char *rc = new char[::strlen(str)+1];
  ::strcpy(rc,str);
  return rc;
}

bool strIsInt(const char*const str, int* i,int* len) 
{
  char *p;
  long l = strtol(str,&p,0);
  if (p == str)
    return false;
  else {
    if (i != NULL)
      *i = l;
    if (len != NULL)
      *len = (p - str);
    return true;
  }
}
bool strIsInt(const char*const str, unsigned* i,int* len) 
{
  char *p;
  long l = strtol(str,&p,0);
  if (p == str || l < 0)
    return false;
  else {
    if (i != NULL)
      *i = (unsigned)l;
    if (len != NULL)
      *len = (p - str);
    return true;
  }
}



// Copies input over to result and if 
// input has any occurence of '%d' in it, replace it with
// the string version of the integer tag. Assume
// result is plenty big.
void copyStringWithTag(char *result,const char *const input, 
		       int tag, const int maxLen)
{

  char buff[1024];
  if (tag == CSWT_EMPTY_TAG)
    sprintf(buff,"%s","");
  else 
    sprintf(buff,"%d",tag);
  int buflen = strlen(buff);
  
  char *result_p = result;
  char *result_endp = result + maxLen;
  const char *input_p = input;
  const char *input_endp = input+strlen(input)+1;
  while (input_p < input_endp) {
    if ( (*input_p != '@') || !(input_p+1 != input_endp && input_p[1] == 'D')) {
      if (result_p < result_endp)
	*result_p++ = *input_p++;
      else {
	error("copyStringWithTag: input to long for length %d",maxLen);
	return;
      }
    } else if (input_p+1 != input_endp && input_p[1] == 'D') {
      input_p += 2;
      for (int i=0;i<buflen;i++) {
	if (result_p < result_endp)
	  *result_p++ = buff[i];
	else {
	  error("copyStringWithTag: input to long for length %d",maxLen);
	  return;
	}
      }
    }
  }
}


//
// Return the file size without moving the file
// pointer.
unsigned long
fsize(FILE*stream) 
{
  long curpos = ftell(stream);
  // rewind
  (void) fseek (stream, 0L, SEEK_END);
  long filesize = ftell(stream);
  (void) fseek (stream, curpos, SEEK_SET);
  return filesize;
}
 

//
// Return the file size of file 'filename' or 
// return 0 if it doesn't exist (i.e., this will
// return 0 if either the file doesn't exist or if
// the file exists and it has zero size).
unsigned long
fsize(const char *const filename) 
{
  if (filename == NULL)
    return 0l;
  FILE*stream;
  if ((stream = fopen(filename,"r")) == NULL)
    return 0l;
  (void) fseek (stream, 0L, SEEK_END);
  long filesize = ftell(stream);
  fclose(stream);
  return filesize;
}

void print_date_string(FILE*f)
{
  time_t tloc;
  struct tm*tms;
  char buf[BUFSIZ];
  time(&tloc);
  tms = localtime(&tloc);
  strftime(buf,BUFSIZ,"%A %B %d %Y, %H:%M:%S %Z",tms);
  fprintf(f,"%s",buf);
}
 
void exit_program_with_status(const int stat)
{
  printf("____ PROGRAM ENDED %sWITH STATUS %d AT ", 
	 (stat == 0 ? "SUCCESSFULLY " : ""),
	 stat);
  print_date_string(stdout);
  printf(" ____\n");
  exit (stat);
}
 
void memory_error()
{
  fprintf(stderr,"ERROR: can't allocate any more memory, malloc() failed. Either decrease your model size, find a better triangulation, use more and/or better pruning, use low-memory options, and/or use a machine with more memory.\n");
  exit_program_with_status(-1);
}

// returns log10(10^v1 + 10^v2)
double log10add(double v1, double v2) 
{
  double small,big;
  if (v1 < v2) {
    small = v1; big = v2;
  } else {
    small = v2; big = v1;
  }
  return (big + ::log10(1+::pow(10.0,small-big)));
}

// report timing in seconds for two calls of
// getrusage()
void reportTiming(// input 
		  const struct rusage& rus,
		  const struct rusage& rue,
		  // output
		  double& userTime, 
		  double& sysTime,
		  // input
		  FILE* outputf)
{

  struct timeval utime;
  double utimef;
  struct timeval stime;
  double stimef;

  /* user time */
  utime.tv_sec = rue.ru_utime.tv_sec - rus.ru_utime.tv_sec ;
  if ( rue.ru_utime.tv_usec < rus.ru_utime.tv_usec ) {
    utime.tv_sec--;
    utime.tv_usec = 1000000l - rus.ru_utime.tv_usec +
      rue.ru_utime.tv_usec;
  } else
    utime.tv_usec = rue.ru_utime.tv_usec -
      rus.ru_utime.tv_usec ;
  utimef = (double)utime.tv_sec + (double)utime.tv_usec/1e6;

  /* system time */
  stime.tv_sec = rue.ru_stime.tv_sec - rus.ru_stime.tv_sec ;
  if ( rue.ru_stime.tv_usec < rus.ru_stime.tv_usec ) {
    stime.tv_sec--;
    stime.tv_usec = 1000000l - rus.ru_stime.tv_usec +
      rue.ru_stime.tv_usec;
  } else
    stime.tv_usec = rue.ru_stime.tv_usec -
      rus.ru_stime.tv_usec ;
  
  stimef = (double)stime.tv_sec + (double)stime.tv_usec/1e6;
  if (outputf != NULL)
    fprintf(outputf,"User: %f, System: %f, CPU %f\n", utimef, stimef, utimef+stimef);
  
  // if (userTime)
  userTime = utimef;
  // if (sysTime)
  sysTime = stimef;
}


/*
 * a version of printf that goes to a C++ string 
 */
int stringprintf(string& str,const char *format, ...)
{
  char buff[512];
  va_list ap;
  va_start(ap,format);
#ifdef __CYGWIN__
  // TODO: as of Fri Sep 07 08:04:46 2007, cygwin doesn't have
  // a vsnprintf, so we use vsprintf for now, but we know
  // that this might cause an error.
  int rc = vsprintf(buff, format, ap);
#else
  int rc = vsnprintf(buff,sizeof(buff), format, ap);
#endif
  va_end(ap);
  str = buff;
  return rc;
}



/*
 * returns number of bits required to represent
 * a value v between [0 <= v < val], where val > 0..
 */
unsigned
bitsRequiredUptoNotIncluding(unsigned val) 
{
  assert ( val > 0 );
  val --;
  unsigned res = 0;
  while (val) {
    res++;
    val >>= 1;
  }
  return res;
}

/*
 * returns ceil(log2(val))
 */
unsigned ceilLog2(unsigned val)
{
  unsigned numBitsSet=0;
  unsigned highBit=0;
  if (val == 0) return 0;
  while (val) { 
    numBitsSet += (val&0x1); 
    highBit++;
    val >>= 1; 
  }
  if (numBitsSet == 1)
    return (highBit-1);
  else
    return (highBit);
}

/*
 * returns integer 2^ceil(log2(val)), least
 * power of two that is >= val. I.e., returns
 * 2^k* where k* = argmin_k { k : 2^k >= val }
 */
unsigned nextPower2(unsigned val)
{
  unsigned numBitsSet=0;
  unsigned highBit=0;
  if (val == 0) return 0;
  while (val) { 
    numBitsSet += (val&0x1); 
    highBit++;
    val >>= 1; 
  }
  if (numBitsSet == 1)
    return (1<<(highBit-1));
  else
    return (1<<(highBit));
}


/*
 * returns the number of bits set in the unsigned
 *
 */
unsigned numBitsSet(unsigned u) {
  unsigned count=0;
  while (u) { count += (u&0x1); u >>= 1; }
  return count;
}


char *CPP_Command() 
{
  char *rc = getenv("GMTK_CPP_CMD");
  if (rc == NULL)
    rc = "cpp";
  return rc;
}
