//
// Simple code to read/write data files in either ASCII or binary
// consisting of either chars, ints, floats, or doubles. If it's an
// ASCII file, supports reading/writting comments starting with COMMENTCHAR, 
// and simple formating with an explicit new-line routine and an 
// indent feature. In binary mode, the formating routines have no effect.
// 
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <ctype.h>

#include "general.h"
VCID("$Header$");
#include "error.h"
#include "fileParser.h"
#include "sArray.h"


#define MAXLINSIZE 131072
#define COMMENTCHAR '%'

#define FLOATWRITESTR   "%0.10e "
#define DOUBLEWRITESTR   "%0.17e "


/////////////////////////////////////////
// Define if you want to pipe all ASCII
// files through the C pre-processor to
// get to use it's macro facilities. Might
// need to change the comment character above.
#define PIPE_ASCII_FILES_THROUGH_CPP

bool 
ioDataStreamFile::errorReturn(char *from,char *msg)
{
  if (msg != NULL) {
    error("%s occurred in %s, file '%s': %s\n",
	  (feof(fh) ? "EOF" : (ferror(fh) ? "Error" : "Strange Error")),
	  from,fileName(),
	  msg);
  }
  return false;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////  iDataStreamFile //////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


#ifdef PIPE_ASCII_FILES_THROUGH_CPP
extern "C" {
  FILE     *popen(const char *, const char *);
  int pclose(FILE *stream);
};
#endif


iDataStreamFile::iDataStreamFile(const char *const _name, bool _Binary)
  : ioDataStreamFile(_name,_Binary)
{
  if (_name == NULL)
    error("Error: Can't open null file for reading.");
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  if (!Binary) {
    if (!strcmp("-",_name)) {
      fh = ::popen("cpp","r");
      if (fh == NULL) {
	error("ERROR: unable to open standard input via cpp");
      }
    }  else {
      if ((fh = ::fopen(_name,"r")) == NULL) {
	error("ERROR: unable to open file (%s) for reading",_name);
      }
      fclose(fh);
      string str = (string)"cpp " + (string)_name;
      fh = ::popen(str.c_str(),"r");    
      if (fh == NULL)
	error("ERROR, can't open file stream from (%s)",_name);
    }
  } else {
    if (!strcmp("-",_name)) {
      fh = stdin;
    } else if ((fh=fopen(_name,"r")) == NULL) {
      error("Error: Can't open file (%s) for reading.",_name);
    }
  }
#else
  if (!strcmp("-",_name)) {
    fh = stdin;
  } else if ((fh=fopen(_name,"r")) == NULL) {
    error("Error: Can't open file (%s) for reading.",_name);
  }
#endif
  if (!Binary) {
    buff = new char[MAXLINSIZE];
    buffp = buff;
    state = GetNextLine;
  }
}

iDataStreamFile::~iDataStreamFile()
{
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  if (pclose(fh) != 0) {
    error("Error: Can't close file.");
  }
#else
  if (fclose(fh) != 0) {
    error("Error: Can't close file.");
  }
#endif
  if (!Binary)
    delete [] buff;
}

bool
iDataStreamFile::prepareNext()
{

  if (state==GetNextLine) {
    bool haveData = false;

    do {
      char *s = fgets(buff,MAXLINSIZE,fh);
      if (s == NULL)
	return false;
      char *cstart = ::index(s,COMMENTCHAR);
      if (cstart != NULL) {
	if (cstart == buff)
	  continue;
	*cstart = '\0';
      }
      buffp = buff;
      while (*buffp && isspace(*buffp)) {
	buffp++;
      }
      if (!*buffp)
	continue;
      haveData = true;
    } while (!haveData);
    state = UseCurLine;
  } else {
    while (*buffp && isspace(*buffp)) {
      buffp++;
    }
    if (!*buffp) {
      state = GetNextLine;
      return prepareNext();
    }
  }
  return true;
}



bool 
iDataStreamFile::readChar(char& c, char *msg) 
{
  if (Binary) {
    size_t rc = fread(&c, sizeof(char), 1,fh);
    if (rc != 1)
      return errorReturn("readChar",msg);
    return true;
  } else {
    if (!prepareNext())
      return errorReturn("readChar",msg);
    c = *buffp++;
    return true;
  }
}


bool 
iDataStreamFile::readStr(char*& str, char *msg) 
{
  sArray<char> tmp(20);
  int len=0;
  if (Binary) {
    char c;
    // read a string up to the next NULL character.
    do {
      size_t rc = fread(&c, sizeof(char), 1,fh);
      if (rc != 1)
	return errorReturn("readChar",msg);
      tmp.growByNIfNeededAndCopy(2,len+1);
      tmp.ptr[len++] = c;
    } while (c != '\0');
  } else {
    if (!prepareNext())
      return errorReturn("readChar",msg);
    // read until a space. Add a null character
    // onto the end of the string.
    char c;
    do {
      tmp.growByNIfNeededAndCopy(2,len+1);
      tmp.ptr[len++] = c = *buffp++;
    } while (!isspace(c) &&  c != '\n');
    tmp.ptr[len-1] = '\0';
  }
  tmp.resizeAndCopy(len);
  str = copyToNewStr(tmp.ptr);
  return true;
}


bool 
iDataStreamFile::readString(string& str, char *msg) 
{
  str.erase();
  if (Binary) {
    char c;
    // read a string up to the next NULL character.
    do {
      size_t rc = fread(&c, sizeof(char), 1,fh);
      if (rc != 1)
	return errorReturn("readChar",msg);
      str += c;
    } while (c != '\0');
  } else {
    if (!prepareNext())
      return errorReturn("readChar",msg);
    // read until a space. Add a null character
    // onto the end of the string.
    char c = *buffp++;
    while (!isspace(c) &&  c != '\n') {
      str += c;
      c = *buffp++;
    }
  }
  return true;
}


bool 
iDataStreamFile::readToken(string& str, const string& tokenChars, char *msg) 
{
  str.erase();
  if (Binary) {
    char c;
    // read a string up to the next NULL character.
    do {
      size_t rc = fread(&c, sizeof(char), 1,fh);
      if (rc != 1)
	return errorReturn("readChar",msg);
      if (c == '\0' || tokenChars.find(c,0) == string::npos)
	break;
      str += c;
    } while (1);
  } else {
    if (!prepareNext())
      return errorReturn("readChar",msg);
    // read until a space. Add a null character
    // onto the end of the string.
    char c = *buffp;
    while (!isspace(c) &&  c != '\n' && tokenChars.find(c,0) != string::npos) {
      str += c;
      c = *++buffp;
    }
  }
  return true;
}



bool 
iDataStreamFile::readInt(int& i, char *msg) 
{
  if (Binary) {
    size_t rc = fread(&i, sizeof(int), 1,fh);
    if (rc != 1)
      return errorReturn("readInt",msg);
    return true;
  } else {
    if (!prepareNext())
      return errorReturn("readInt",msg);
    char *ptr;
    long l = strtol(buffp,&ptr,0);
    if (ptr == buffp) {
      error("readInt: Can't form int at (%s) in file (%s). %s",buffp,
	    fileName(),
	    (msg?msg:""));
    } else
      buffp = ptr;
    i = (int)l;
    return true;
  }
}


bool 
iDataStreamFile::readFloat(float& f, char *msg) 
{
  if (Binary) {
    size_t rc = fread(&f, sizeof(float), 1,fh);
    if (rc != 1)
      return errorReturn("readFloat",msg);
    return true;
  } else {
    if (!prepareNext())
      return errorReturn("readFloat",msg);
    char *ptr;
    f = strtod(buffp,&ptr);
    if (ptr == buffp) {
      error("readFloat: Can't form float at (%s) in file (%s). %s",buffp,
	    fileName(),
	    (msg?msg:""));
    } else
      buffp = ptr;
    return true;
  }
}


bool 
iDataStreamFile::readDouble(double& d, char *msg) 
{
  if (Binary) {
    size_t rc = fread(&d, sizeof(double), 1,fh);
    if (rc != 1)
      return errorReturn("readDouble",msg);
    return true;
  } else {
    if (!prepareNext())
      return errorReturn("readDouble",msg);
    char *ptr;
    d = strtod(buffp,&ptr);
    if (ptr == buffp) {
      error("readDouble: Can't form double at (%s) in file (%s). %s",buffp,
	    fileName(),
	    (msg?msg:""));
    } else
      buffp = ptr;
    return true;
  }
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////  oDataStreamFile //////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


oDataStreamFile::oDataStreamFile(const char *const _name,bool _Binary)
  : ioDataStreamFile(_name,_Binary)
{
  if (_name == NULL)
    error("Error: Can't open null file for reading.");
  if (!strcmp("-",_name)) {
    fh = stdout;
  } else if (!strcmp("--",_name)) {
    fh = stderr;
  } else if ((fh=fopen(_name,"w")) == NULL) {
    error("Error: Can't open file (%s) for writing.",_name);
  }
  char buff[2048];
  sprintf(buff,FLOATWRITESTR,0.0);
  float_space = strlen(buff);
  sprintf(buff,DOUBLEWRITESTR,0.0);
  double_space = strlen(buff);
}

oDataStreamFile::~oDataStreamFile()
{
  if (fclose(fh) != 0) {
    error("Error: Can't close file.");
  }
}

bool oDataStreamFile::writeStr(const char *const str,char *msg)
{
  if (Binary) {
    const int len = strlen(str)+1;
    size_t rc = fwrite(str, sizeof(char), len,fh);
    if ((int)rc != len)
      return errorReturn("writeStr",msg);
    return true;
  } else {
    if (fprintf(fh,"%s ",str) == 0) 
      return errorReturn("writeStr",msg);
    return true;
  }
}


bool oDataStreamFile::writeString(const string& str,char *msg)
{
  if (Binary) {
    const int len = str.length()+1;
    size_t rc = fwrite(str.c_str(), sizeof(char), len,fh);
    if ((int)rc != len)
      return errorReturn("writeStr",msg);
    return true;
  } else {
    if (fprintf(fh,"%s ",str.c_str()) == 0) 
      return errorReturn("writeStr",msg);
    return true;
  }
}


bool oDataStreamFile::writeChar(const char c,char *msg)
{
  if (Binary) {
    size_t rc = fwrite(&c, sizeof(char), 1,fh);
    if (rc != 1)
      return errorReturn("writeChar",msg);
    return true;
  } else {
    if (fprintf(fh,"%c ",c) == 0) 
      return errorReturn("writeChar",msg);
    return true;
  }
}


bool oDataStreamFile::writeInt(const int i,char *msg)
{
  if (Binary) {
    size_t rc = fwrite(&i, sizeof(int), 1,fh);
    if (rc != 1)
      return errorReturn("writeInt",msg);
    return true;
  } else {
    if (fprintf(fh,"%d ",i) == 0) 
      return errorReturn("writeInt",msg);
    return true;
  }
}


bool oDataStreamFile::writeFloat(const float f,char *msg)
{
  if (Binary) {
    size_t rc = fwrite(&f, sizeof(float), 1,fh);
    if (rc != 1)
      return errorReturn("writeFloat",msg);
    return true;
  } else {
    if (fprintf(fh,FLOATWRITESTR,f) == 0) 
      return errorReturn("writeFloat",msg);
    return true;
  }
}


bool oDataStreamFile::writeDouble(const double d,char *msg)
{
  if (Binary) {
    size_t rc = fwrite(&d, sizeof(double), 1,fh);
    if (rc != 1)
      return errorReturn("writeDouble",msg);
    return true;
  } else {
    if (fprintf(fh,DOUBLEWRITESTR,d) == 0) 
      return errorReturn("writeDouble",msg);
    return true;
  }
}


bool oDataStreamFile::writeComment(char *comment, ...)
{
  if (Binary) {
    // do nothing.
    return true;
  } else {
    if (fprintf(fh,"%c ",COMMENTCHAR) == 0)
      return errorReturn("writeComment",comment);
    va_list ap;
    va_start(ap,comment);
    (void) vfprintf(fh, comment, ap);
    va_end(ap);
    return true;
  }
}



bool oDataStreamFile::indent(const int i,const bool doubSpace,char *msg)
{
  if (Binary) {
  } else {
    int tmp = i;
    while (tmp--) {
      int tmp2=(doubSpace?double_space:float_space);
      while (tmp2--) {
	if (fprintf(fh," ") == 0) 
	  return errorReturn("indent",msg);
      }
    }
  }
  return true;
}


bool oDataStreamFile::space(const int numSpaceChars,char *msg)
{
  if (Binary) {
  } else {
    int tmp = numSpaceChars;
    while (tmp--) {
      if (fprintf(fh," ") == 0) 
	return errorReturn("space",msg);
    }
  }
  return true;
}


bool oDataStreamFile::nl(char *msg)
{
  if (Binary) {
  } else {
    if (fprintf(fh,"\n") == 0) 
      return errorReturn("indent",msg);
  }
  return true;
}


bool oDataStreamFile::flush(char *msg)
{
  if (fflush(fh) != 0)
    return errorReturn("flush",msg);
  return true;
}

void oDataStreamFile::rewind()
{
  (void) fseek (fh, 0L, SEEK_SET);
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////  Test Code  ///////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

#ifdef MAIN


int main()
{
  char c;
  int i;
  float f;
  double d;
  char *str;
  string cppstr = "This_is_a_stupid_slow_C++_string.";

  float far1[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };

  float far2[sizeof(far1)/sizeof(float)];

  bool bin=false;
  {
    oDataStreamFile of("/tmp/foo_ascii",bin);
    of.writeChar('c'); of.nl();
    of.writeInt(10003); of.nl();
    of.writeFloat(10003.043); of.nl();
    of.indent(1,false); of.writeDouble(10003.043343434); of.nl();
    of.writeStr("This_is_a_String"); of.nl();
    of.write(far1,sizeof(far1)/sizeof(float),"writing array");
    of.write(cppstr);
  }

  bin=false;
  {
    iDataStreamFile inf("/tmp/foo_ascii",bin);
    inf.readChar(c);
    inf.readInt(i);
    inf.readFloat(f);
    inf.readDouble(d);
    inf.readStr(str);
    inf.read(far2,sizeof(far1)/sizeof(float),"reading array");
    inf.read(cppstr);
    printf("%c %d %f %f %s\n",c,i,f,d,str);
    for (int i=0;i<(int)(sizeof(far1)/sizeof(float));i++) {
      printf("%f ",far2[i]);
    }
    printf("\n");
    printf("cppstr=(%s)\n",cppstr.c_str());
    printf("\n");
  }

  bin = true;
  {
    oDataStreamFile of("/tmp/foo_bin",bin);
    of.writeChar(c); of.nl();
    of.writeInt(i); of.nl();
    of.writeFloat(f); of.nl();
    of.indent(1,false); of.writeDouble(d); of.nl();
    of.writeStr(str); of.nl();
    of.write(far1,sizeof(far1)/sizeof(float),"writing array");
    of.write(cppstr);
  }

  bin=true;
  {
    iDataStreamFile inf("/tmp/foo_bin",bin);
    inf.readChar(c);;
    inf.readInt(i);
    inf.readFloat(f);
    inf.readDouble(d);
    inf.readStr(str);
    inf.read(far2,sizeof(far2)/sizeof(float),"reading array");
    inf.read(cppstr);
    printf("%c %d %f %f %s\n",c,i,f,d,str);
    for (int i=0;i<(int)(sizeof(far1)/sizeof(float));i++) {
      printf("%f ",far2[i]);
    }
    printf("\n");
    printf("cppstr=(%s)\n",cppstr.c_str());
    printf("\n");

  }

}



#endif
