//
// Simple code to read/write data files in either ASCII or binary
// consisting of either chars, ints, floats, or doubles. If it's an
// ASCII file, supports reading/writting comments starting with COMMENTCHAR, 
// and simple formating with an explicit new-line routine and an 
// indent feature. In binary mode, the formating routines have no effect.
// 
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <ctype.h>
#include <errno.h>
#include <string>

#ifdef __CYGWIN__
// added for cygwin Tue May 14 12:30:46 2002
extern "C" { char *index(const char* str, int c); }
#endif

#include "general.h"
VCID("$Header$");
#include "error.h"
#include "fileParser.h"
#include "sArray.h"


#define MAXLINSIZEPLUS1 (262144)
#define COMMENTCHAR '%'

#define FLOATWRITESTR   "%0.10e "
#define DOUBLEWRITESTR   "%0.17e "


#ifdef PIPE_ASCII_FILES_THROUGH_CPP
#define CPP_DIRECTIVE_CHAR '#'
#endif


bool 
ioDataStreamFile::errorReturn(char *from,char *msg)
{
  if (msg != NULL) {
    error("%s occurred in %s, file '%s' line %d: %s\n",
	  (feof(fh) ? "EOF" : (ferror(fh) ? "Error" : "Strange Error")),
	  from,fileName(),lineNo(),
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
#ifndef DECLARE_POPEN_FUNCTIONS_EXTERN_C
extern "C" {
#ifdef __CYGWIN__
 FILE     *popen(const char *, const char *) __THROW;
 int pclose(FILE *stream) __THROW;
#endif
};
#endif
#endif

#ifdef PIPE_ASCII_FILES_THROUGH_CPP
iDataStreamFile::iDataStreamFile(const char *const _name, bool _Binary, bool _cppIfAscii, const char *const _cppCommandOptions)
  : ioDataStreamFile(_name,_Binary), cppIfAscii(!_Binary && _cppIfAscii)
#else
iDataStreamFile::iDataStreamFile(const char *const _name, bool _Binary)
  : ioDataStreamFile(_name,_Binary)
#endif
{
  if (_name == NULL)
    error("Error: Can't open null file for reading.");
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  if (!Binary) {
    if (cppIfAscii) {
      string cppCommand = string("cpp");
      if (_cppCommandOptions != NULL) {
	cppCommand = cppCommand + string(" ") + string(_cppCommandOptions);
      }
      if (!strcmp("-",_name)) {
	fh = ::popen(cppCommand.c_str(),"r");
	if (fh == NULL) {
	  error("ERROR: unable to open standard input via cpp");
	}
      }  else {
	// make sure the file  exists first.
	if ((fh = ::fopen(_name,"r")) == NULL) {
	  error("ERROR: unable to open file (%s) for reading",_name);
	}
	fclose(fh);

	// add path of file to include directory paths.
	string path = _name;
	unsigned slashPos = path.rfind("/");
	if (slashPos != string::npos) {
	  // then '/' is found
	  cppCommand = cppCommand + " -I" + path.substr(0,slashPos);
	}
	// Lastly, add CWD to default CPP command options for include files
	// (i.e., we look for include files in CWD only if all previous
	// ones fail, cpp has this behavior.
	cppCommand = cppCommand + " -I.";

	cppCommand = cppCommand + string(" ") + string(_name);
	fh = ::popen(cppCommand.c_str(),"r");    
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
    buff = new char[MAXLINSIZEPLUS1];
    buffp = buff;
    state = GetNextLine;
  }
}

iDataStreamFile::~iDataStreamFile()
{
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  if (cppIfAscii) {
    if (pclose(fh) != 0) {
      warning("WARNING: Can't close pipe 'cpp %s'.",fileName());
    }
  } else {
    if (fclose(fh) != 0) {
      warning("WARNING: Can't close file '%s'.",fileName());
    }
  }
#else
  if (fclose(fh) != 0) {
    warning("WARNING: Can't close file '%s'.",fileName());
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
      char *s = fgets(buff,MAXLINSIZEPLUS1,fh);
      if (s == NULL)
	return false;
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
      if (*s == CPP_DIRECTIVE_CHAR) {
	// check if there is a filename/lineNo change.
	{
	  if (strcmp("#line",s) == 0) {
	    // then this is probabilty of the form: ^"#line"{ws}{int}{ws}{string}.* 
	    // skip over '#line' string
	    s += 5;
	  } else { 
	    // then might be of the form: ^"#"{ws}{int}{ws}{string}.*
	    // skip over '#' character
	    s += 1; 
	  }
	  // skip over ws
	  while (*s && isspace(*s))
	    s++;
	  // try to get an int.
	  char *ss;
	  long newLineNo = strtol(s,&ss,0);
	  if (ss == s) {
	    // then not an int, so we just skip this line
	    continue;
	  }
	  // ok, got int, presumably line number is in newLineNo
	  s = ss;
	  // skip over ws
	  while (*s && isspace(*s))
	    s++;
	  if (!(*s)) {
	    // then at end of line, so no file name, so skip line
	    continue;
	  }
	  // scan filename until next ws or eol.
	  ss = s;
	  while (*ss && !isspace(*ss) && (*ss != '\n'))
	    ss++;
	  // end the line if not ended already.
	  *ss = '\0';
	  // we're done. We've got new line number and get remainder
	  // of string the presumably new file name.
	  _curLineNo = newLineNo;
	  _fileName = s;
	  // now remove any " characters from filename.
	  if (_fileName[0] == '"') {
	    _fileName.erase(0,1);
	  }
	  if (_fileName[_fileName.size()-1] == '"') {	  
	    _fileName.erase(_fileName.size()-1,1);
	  }
	}
	continue;
      }
#endif
      // otherwise, we've successfully read the next line.
      _curLineNo ++;

      if (::strlen(s) == (MAXLINSIZEPLUS1-1))
	error("ERROR: maximum line length of %d reached in ASCII file '%s', line %d\n",
	      fileName(),lineNo());

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


void iDataStreamFile::rewind()
{
  assert ( ! cppIfAscii );
  if (::fseek (fh, 0L, SEEK_SET) != 0)
    error("ERROR: trouble seeking to beginning of file '%s', %s\n",
	  fileName(),strerror(errno));
  _curLineNo = 0;
  state = GetNextLine;
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
	return errorReturn("readStr",msg);
      tmp.growByNIfNeededAndCopy(2,len+1);
      tmp.ptr[len++] = c;
    } while (c != '\0');
  } else {
    if (!prepareNext())
      return errorReturn("readStr",msg);
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
    size_t rc = fread(&c, sizeof(char), 1,fh);
    if (rc != 1)
      return errorReturn("readString",msg);
    while (c != '\0') {
      str += c;
      size_t rc = fread(&c, sizeof(char), 1,fh);
      if (rc != 1)
	return errorReturn("readString",msg);
    }
    if (str.size() == 0)
	return errorReturn("readString, zero length string",msg);
  } else {
    if (!prepareNext())
      return errorReturn("readString",msg);
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
iDataStreamFile::readStringUntil(
  string& str, 
  const char delimiter, 
  bool spaceIsDelimiter, 
  char *msg ) 
{
  str.erase();
  if (Binary) {
    char c;
    // read a string up to the next delimiter character.
    size_t rc = fread(&c, sizeof(char), 1,fh);
    if (rc != 1)
      return errorReturn("readStringUntil",msg);
    while ((c != delimiter) && ((!spaceIsDelimiter) || (!isspace(c)))) {
      str += c;
      size_t rc = fread(&c, sizeof(char), 1,fh);
      if (rc != 1)
	return errorReturn("readStringUntil",msg);
    }
    if (str.size() == 0)
	return errorReturn("readStringUntil, zero length string",msg);
  } else {
    if (!prepareNext()) {
      return errorReturn("readStringUntil_aaa",msg);
    }
    char c = *buffp++;

    while ((c != delimiter) && ((!spaceIsDelimiter) || (!isspace(c)))) {
      str += c;

      if (c == '\n') {
        if (!prepareNext()) {
          return errorReturn("readStringUntil_whee",msg);
        }
      }
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
	return errorReturn("readToken",msg);
      if (c == '\0' || tokenChars.find(c,0) == string::npos)
	break;
      str += c;
    } while (1);
  } else {
    if (!prepareNext())
      return errorReturn("readToken",msg);
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
      error("readInt: Can't form int at (%s) in file (%s) line %d. %s",buffp,
	    fileName(),lineNo(),
	    (msg?msg:""));
    } else
      buffp = ptr;
    i = (int)l;
    return true;
  }
}


bool 
iDataStreamFile::readUnsigned(unsigned& i, char *msg) 
{
  if (Binary) {
    size_t rc = fread(&i, sizeof(unsigned), 1,fh);
    if (rc != 1)
      return errorReturn("readUnsigned",msg);
    return true;
  } else {
    if (!prepareNext())
      return errorReturn("readUnsigned",msg);
    char *ptr;
    unsigned long l = strtoul(buffp,&ptr,0);
    if (ptr == buffp) {
      error("readUnsigned: Can't form unsigned at (%s) in file (%s) line %d. %s",buffp,
	    fileName(),lineNo(),
	    (msg?msg:""));
    } else
      buffp = ptr;
    i = (unsigned)l;
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
      error("readFloat: Can't form float at (%s) in file (%s) line %d. %s",buffp,
	    fileName(),lineNo(),
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
      error("readDouble: Can't form double at (%s) in file (%s) line %d. %s",buffp,
	    fileName(),lineNo(),
	    (msg?msg:""));
    } else
      buffp = ptr;
    return true;
  }
}


bool
iDataStreamFile::readLine(char *&line, size_t n, char *msg) {
  if ( Binary )
    error("getline: Can't getline in a binary file.");

  /*
  ssize_t rc = ::getline(&line, &n, fh);
  if ( rc < 0 )
	error("getline", msg);
  
  buffp += rc;
  */

  int len=0;

  if (!prepareNext())
    return errorReturn("getline",msg);

  // read until '\n'. Add a null character
  // onto the end of the string.
  char c;
  do {
    line[len++] = c = *buffp++;
  } while ( c != '\n' );
  line[len-1] = '\n';
  line[len] = '\0';

  return true;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////  oDataStreamFile //////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


oDataStreamFile::oDataStreamFile(const char *const _name,bool _Binary, bool _Append)
  : ioDataStreamFile(_name,_Binary)
{
  if (_name == NULL)
    error("Error: Can't open null file for reading.");
  if (!strcmp("-",_name)) {
    fh = stdout;
  } else if (!strcmp("--",_name)) {
    fh = stderr;
  } else if (_Append) {
    if ((fh=fopen(_name,"a")) == NULL)
      error("Error: Can't open file (%s) for appending.",_name);
  } else { 
    if ((fh=fopen(_name,"w")) == NULL)
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
    warning("WARNING: Can't close file '%s'.",fileName());
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

bool oDataStreamFile::writeUnsigned(const unsigned int u,char *msg)
{
  if (Binary) {
    size_t rc = fwrite(&u, sizeof(int), 1,fh);
    if (rc != 1)
      return errorReturn("writeUnsigned",msg);
    return true;
  } else {
    if (fprintf(fh,"%u ",u) == 0) 
      return errorReturn("writeUnsigned",msg);
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
  _curLineNo++;
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
  _curLineNo = 0;
  if (::fseek (fh, 0L, SEEK_SET) != 0)
    error("ERROR: trouble seeking to beginning of file '%s', %s\n",
	  fileName(),strerror(errno));
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
  string cppstr = "This_is_a_C++_string.";

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
