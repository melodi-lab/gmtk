// 
// Simple code to parse a file
//
// $Header$
//
//
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu


#ifndef FILEPARSER_H
#define FILEPARSER_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

class ioDataStreamFile {

 protected:
  FILE *fh;
  bool Binary;
  bool errorReturn(char *from,char *msg);

 public:

  void rewind() { ::rewind(fh); }
  ioDataStreamFile(bool _Binary = false) : Binary(_Binary) {}

};


class iDataStreamFile : public ioDataStreamFile {

  char *buff;
  char *buffp;
  enum State { GetNextLine, UseCurLine } state;

 public:
  iDataStreamFile(char *_name, bool _Binary = false);
  ~iDataStreamFile();

  bool prepareNext();

  // These return true if a successful value is returned,
  // false otherwise. If 'msg' is provided and an error
  // occurs, msg will be printed and the program will die.

  // type explicit
  bool readChar(char& c,char *msg=NULL);
  bool readStr(char*& str,char *msg=NULL);
  bool readInt(int& i,char *msg=NULL);
  bool readFloat(float& f,char *msg=NULL);
  bool readDouble(double& d,char *msg=NULL);


  // type implicit
  bool read(char*& str,char *msg=NULL) { return readStr(str,msg); }
  bool read(char& c,char *msg=NULL) { return readChar(c,msg); }
  bool read(int& i,char *msg=NULL) { return readInt(i,msg); }
  bool read(unsigned& i,char *msg=NULL) { return readInt((int&)i,msg); }
  bool read(float& f,char *msg=NULL) { return readFloat(f,msg); }
  bool read(double& d,char *msg=NULL) { return readDouble(d,msg); }


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

};



class oDataStreamFile : public ioDataStreamFile {
  int float_space;
  int double_space;

 public:
  oDataStreamFile(char *_name, bool _Binary = false);
  ~oDataStreamFile();

  // type explicit
  bool writeStr(const char * const str, char *msg=NULL);
  bool writeChar(const char c, char *msg=NULL);
  bool writeInt(const int i,char *msg=NULL);
  bool writeFloat(const float f,char *msg=NULL);
  bool writeDouble(const double d,char *msg=NULL);
  bool writeComment(char *comment, ...);
  bool indent(const int i,const bool doubSpace, char *msg=NULL);
  bool nl(char *msg=NULL);
  bool flush(char *msg=NULL);
  void rewind();

  // type implicit
  bool write(const char *const str,char *msg=NULL) { return writeStr(str,msg); }
  bool write(const char c, char *msg=NULL) { return writeChar(c,msg); }
  bool write(const int i,char *msg=NULL) { return writeInt(i,msg); }
  bool write(const unsigned i,char *msg=NULL) { return writeInt((int)i,msg); }
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



};



#endif
