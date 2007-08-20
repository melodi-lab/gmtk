// 
// Simple code to parse a file
//
// $Header$
//
//
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu


#ifndef FILEPARSER_H
#define FILEPARSER_H


/////////////////////////////////////////
// Define if you want to pipe all ASCII
// files through the C pre-processor to
// get to use it's macro facilities. Might
// need to change the comment character above.
#define PIPE_ASCII_FILES_THROUGH_CPP


#include <string>
#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "general.h"

class ioDataStreamFile {

 protected:
  FILE *fh;
  bool Binary;
  bool errorReturn(const char *from,const char *msg);
  // the current file name
  string _fileName;
  // the current line number in the current file name (always 0 for binary file).
  int _curLineNo;

 public:

  ioDataStreamFile(const char *name,bool _Binary = false) : 
    Binary(_Binary), _fileName(name), _curLineNo(0) {}
  ~ioDataStreamFile() { }

  bool binary() { return Binary; }
  const char *const fileName() { return _fileName.c_str(); }
  int lineNo() { return _curLineNo; }

  long ftell() const { return(::ftell(fh)); }
  int  fseek ( long offset , int origin ) { return(::fseek(fh,offset,origin)); }
};


class iDataStreamFile : public ioDataStreamFile {

  char *buff;
  char *buffp;
  enum State { GetNextLine, UseCurLine } state;
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  const bool cppIfAscii;
#endif
  const char extraCommentChar;

 public:

#define IDATASTREAMFILE_DEFAULT_EXTRA_COMMENT_CHAR '\1'

#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  iDataStreamFile(const char *_name, bool _Binary = false, bool _cppIfAscii = true, const char * const _cppCommandOptions = NULL,const char extraCommentChar = IDATASTREAMFILE_DEFAULT_EXTRA_COMMENT_CHAR);
#else
  iDataStreamFile(const char *_name, bool _Binary = false,const char extraCommentChar = IDATASTREAMFILE_DEFAULT_EXTRA_COMMENT_CHAR);
#endif

  ~iDataStreamFile();

  bool prepareNext();
  
  void rewind();

  bool isEOF() { 
    if(feof(fh)) return true;
    else         return false;
  }

  // These return true if a successful value is returned,
  // false otherwise. If 'msg' is provided and an error
  // occurs, msg will be printed and the program will die.

  // type explicit
  bool readChar(char& c,const char *msg=NULL);
  bool readStr(char*& str,const char *msg=NULL);
  bool readInt(int& i,const char *msg=NULL);
  bool readUnsigned(unsigned& i,const char *msg=NULL);

  bool readFloat(float& f,const char *msg=NULL);
  bool readDouble(double& d,const char *msg=NULL);

  bool readFloatVec(float* fp,const unsigned len,const char *msg=NULL);
  bool readDoubleVec(double* dp,const unsigned len,const char *msg=NULL);


  bool readString(string& str,const char *msg=NULL);
  bool readToken(string& str,const string& tokenChars,const char *msg=NULL);

  bool readIfMatch(const string& matchTokenStr,const char *msg=NULL);

  // type implicit
  bool read(char*& str,const char *msg=NULL) { return readStr(str,msg); }
  bool read(char& c,const char *msg=NULL) { return readChar(c,msg); }
  bool read(int& i,const char *msg=NULL) { return readInt(i,msg); }
  bool read(unsigned& i,const char *msg=NULL) { return readUnsigned(i,msg); }
  bool read(float& f,const char *msg=NULL) { return readFloat(f,msg); }
  bool read(double& d,const char *msg=NULL) { return readDouble(d,msg); }
  bool read(float* fp, const unsigned len, const char *msg=NULL) { return readFloatVec(fp,len,msg); }
  bool read(double* dp, const unsigned len, const char *msg=NULL) { return readDoubleVec(dp,len,msg); }
  bool read(string& str,const char *msg=NULL) { return readString(str,msg); }
  bool read(string& str,const string& tokenChars,const char *msg=NULL) {
    return readToken(str,tokenChars,msg);
  }

  bool readStringUntil(
    string& str, 
    const char delimiter, 
    bool spaceIsDelimiter, 
    const char *msg=NULL );

  // read a line consisting of a total of n characters including the final '\0'.
  // Memory allocation is done externally.
  bool readLine(char * lineptr, size_t n, const char *msg = NULL);

  char peekChar(const char *msg = NULL);

  template <class T>
  bool readArray(T* location, const int length, const char *msg = NULL) 
  {
    assert ( length >= 0 && location != NULL);
    
    if (length == 0)
      return true;

    bool rc;
    int i=0; do {
      rc = read(location[i],msg);
      i++;
    } while (rc == true && i < length );
    return rc;
  }

  template <class T>
  bool read(T* location, const int length, const char *msg = NULL)
  { return readArray(location,length,msg); }


  template <class T>
  bool readVector(vector<T>& location, const int length, const char *msg = NULL) 
  {
    assert ( length >= 0 );
    
    if (length == 0)
      return true;

    bool rc;
    location.reserve(length);
    int i=0; do {
      rc = read(location[i],msg);
      i++;
    } while (rc == true && i < length );
    return rc;
  }

  template <class T>
  bool read(vector<T>& location, const int length, const char *msg = NULL)
  { return readVector(location,length,msg); }



};



class oDataStreamFile : public ioDataStreamFile {
  int float_space;
  int double_space;

 public:
  oDataStreamFile(const char *_name, bool _Binary = false, bool _Append = false);
  ~oDataStreamFile();

  // type explicit
  bool writeStr(const char * const str, const char *msg=NULL,const bool writeSpaceSuffixAscii=true);
  bool writeString(const string& str,const char *msg=NULL,const bool writeSpaceSuffixAscii=true);
  bool writeChar(const char c, const char *msg=NULL,const bool writeSpaceSuffixAscii=true);
  bool writeInt(const int i,const char *msg=NULL);
  bool writeUnsigned(const unsigned int u,const char *msg=NULL);
  bool writeFloat(const float f,const char *msg=NULL);
  bool writeDouble(const double d,const char *msg=NULL);
  bool writeFloatVec(const float* fp,unsigned len, const char *msg=NULL);
  bool writeDoubleVec(const double* dp,unsigned len,const char *msg=NULL);

  bool writeComment(const char *comment, ...);
  bool indent(const int i,const bool doubSpace, const char *msg=NULL);
  bool space(const int numSpaceChars, const char *msg=NULL);
  bool nl(const char *msg=NULL);
  bool flush(const char *msg=NULL);
  void rewind();

  // type implicit
  bool write(const char *const str,const char *msg=NULL,const bool writeSpaceSuffixAscii=true) { return writeStr(str,msg,writeSpaceSuffixAscii); }
  bool write(const string& str,const char *msg=NULL,const bool writeSpaceSuffixAscii=true) { return writeString(str,msg,writeSpaceSuffixAscii); }
  bool write(const char c, const char *msg=NULL,const bool writeSpaceSuffixAscii=true) { return writeChar(c,msg,writeSpaceSuffixAscii); }

  bool write(const int i,const char *msg=NULL) { return writeInt(i,msg); }
  bool write(const unsigned int u,const char *msg=NULL) { return writeUnsigned(u,msg); }
#ifdef _AIX
  bool write(const size_t i,const char *msg=NULL) { return writeUnsigned(i,msg); }
#endif
  bool write(const float f,const char *msg=NULL) { return writeFloat(f,msg); }
  bool write(const double d,const char *msg=NULL) { return writeDouble(d,msg); }
  bool write(const float* fp,unsigned len,const char *msg=NULL) { return writeFloatVec(fp,len,msg); }
  bool write(const double* dp,unsigned len,const char *msg=NULL) { return writeDoubleVec(dp,len,msg); }


  template <class T>
  bool writeArray(T* location, const int length, const char *msg = NULL) 
  {
    assert ( length >= 0 && location != NULL);
    
    if (length == 0)
      return true;

    bool rc;
    int i=0; do {
      rc = write(location[i],msg);
      i++;
    } while (rc == true && i < length );
    return rc;
  }

  // template <class T>
  // bool write(T* location, const int length, const char *msg = NULL)
  // { return writeArray(location,length,msg); }


  template <class T>
  bool writeVector(vector<T> location, const char *msg = NULL) 
  {
    if (location.length() == 0)
      return true;

    bool rc;
    int i=0; do {
      rc = write(location[i],msg);
      i++;
    } while (rc == true && i < location.length() );
    return rc;
  }

  // template <class T>
  // bool write(vector<T> location, const char *msg = NULL)
  // { return writeVector(location,msg); }



};



#endif
