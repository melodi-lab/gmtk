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
  bool errorReturn(char *from,char *msg);
  const char *const _fileName;

 public:

  ioDataStreamFile(const char *name,bool _Binary = false) : 
    Binary(_Binary), _fileName(copyToNewStr(name)) {}
  ~ioDataStreamFile() { delete [] _fileName; }

  bool binary() { return Binary; }
  const char *const fileName() { return _fileName; }

};


class iDataStreamFile : public ioDataStreamFile {

  char *buff;
  char *buffp;
  enum State { GetNextLine, UseCurLine } state;
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  const bool cppIfAscii;
#endif

 public:
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  iDataStreamFile(const char *_name, bool _Binary = false, bool _cppIfAscii = true, const char * const _cppCommandOptions = NULL);
#else
  iDataStreamFile(const char *_name, bool _Binary = false);
#endif

  ~iDataStreamFile();

  bool prepareNext();
  
  void rewind();


  // These return true if a successful value is returned,
  // false otherwise. If 'msg' is provided and an error
  // occurs, msg will be printed and the program will die.

  // type explicit
  bool readChar(char& c,char *msg=NULL);
  bool readStr(char*& str,char *msg=NULL);
  bool readInt(int& i,char *msg=NULL);
  bool readUnsigned(unsigned& i,char *msg=NULL);
  bool readFloat(float& f,char *msg=NULL);
  bool readDouble(double& d,char *msg=NULL);
  bool readString(string& str,char *msg=NULL);
  bool readToken(string& str,const string& tokenChars,char *msg=NULL);


  // type implicit
  bool read(char*& str,char *msg=NULL) { return readStr(str,msg); }
  bool read(char& c,char *msg=NULL) { return readChar(c,msg); }
  bool read(int& i,char *msg=NULL) { return readInt(i,msg); }
  bool read(unsigned& i,char *msg=NULL) { return readUnsigned(i,msg); }
  bool read(float& f,char *msg=NULL) { return readFloat(f,msg); }
  bool read(double& d,char *msg=NULL) { return readDouble(d,msg); }
  bool read(string& str,char *msg=NULL) { return readString(str,msg); }
  bool read(string& str,const string& tokenChars,char *msg=NULL) {
    return readToken(str,tokenChars,msg);
  }


  template <class T>
  bool readArray(T* location, const int length, char *msg = NULL) 
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
  bool read(T* location, const int length, char *msg = NULL)
  { return readArray(location,length,msg); }


  template <class T>
  bool readVector(vector<T>& location, const int length, char *msg = NULL) 
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
  bool read(vector<T>& location, const int length, char *msg = NULL)
  { return readVector(location,length,msg); }



};



class oDataStreamFile : public ioDataStreamFile {
  int float_space;
  int double_space;

 public:
  oDataStreamFile(const char *_name, bool _Binary = false);
  ~oDataStreamFile();

  // type explicit
  bool writeStr(const char * const str, char *msg=NULL);
  bool writeString(const string& str,char *msg=NULL);
  bool writeChar(const char c, char *msg=NULL);
  bool writeInt(const int i,char *msg=NULL);
  bool writeUnsigned(const unsigned int u,char *msg=NULL);
  bool writeFloat(const float f,char *msg=NULL);
  bool writeDouble(const double d,char *msg=NULL);
  bool writeComment(char *comment, ...);
  bool indent(const int i,const bool doubSpace, char *msg=NULL);
  bool space(const int numSpaceChars, char *msg=NULL);
  bool nl(char *msg=NULL);
  bool flush(char *msg=NULL);
  void rewind();

  // type implicit
  bool write(const char *const str,char *msg=NULL) { return writeStr(str,msg); }
  bool write(const string& str,char *msg=NULL) { return writeString(str,msg); }
  bool write(const char c, char *msg=NULL) { return writeChar(c,msg); }
  bool write(const int i,char *msg=NULL) { return writeInt(i,msg); }
  bool write(const unsigned int u,char *msg=NULL) { return writeUnsigned(u,msg); }
#ifdef _AIX
  bool write(const size_t i,char *msg=NULL) { return writeUnsigned(i,msg); }
#endif
  bool write(const float f,char *msg=NULL) { return writeFloat(f,msg); }
  bool write(const double d,char *msg=NULL) { return writeDouble(d,msg); }


  template <class T>
  bool writeArray(T* location, const int length, char *msg = NULL) 
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

  template <class T>
  bool write(T* location, const int length, char *msg = NULL)
  { return writeArray(location,length,msg); }


  template <class T>
  bool writeVector(vector<T> location, char *msg = NULL) 
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

  template <class T>
  bool write(vector<T> location, char *msg = NULL)
  { return writeVector(location,msg); }



};



#endif
