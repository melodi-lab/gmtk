//
// General miscellaneous stuff that belongs nowhere else.
// 
// Written by: Jeff Bilmes
//             bilmes@ee.wasington.edu

#include <ctype.h>
#include <sys/time.h>
#include <time.h>


#include "general.h"
#include "rand.h"
#include "error.h"


VCID("$Header$");

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
    if (*input_p != '%') {
      if (result_p < result_endp)
	*result_p++ = *input_p++;
      else {
	error("copyStringWithTag: input to long for length %d",maxLen);
	return;
      }
    }
    else if (input_p+1 != input_endp && input_p[1] == 'd') {
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
  strftime(buf,BUFSIZ,"%A %B %d %Y, %T %Z",tms);
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
