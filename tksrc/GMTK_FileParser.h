/*
 * GMTK_GM.h
 * Parses a text file giving the basic GM structure
 * over hidden varialbes.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 */

#ifndef GMTK_FILEPARSER_H
#define GMTK_FILEPARSER_H

#include <stdio.h>
#include <stdlib.h>

//  #include "sArray.h"
//  #include "logp.h"

/*
#include "GMTK_GM.h" 
*/

class FileParser
{


  //  GMTK_GM& gm; 

  //////////////////////////////////
  // The file we are currently parsing.
  FILE *ifile;



  //////////////////////////////////////////////
  // prepares the next input token and fills the
  // structure
  void prepareNextToken();

  //////////////////////////////////////////////
  // ensures that we are not at EOF, and dies with
  // a parse error message if we are.
  void ensureNotEOF(const char *const str = NULL);
  // returns true if we're at EOF
  bool endOfFile() { return (tokenInfo.rc == 0); }

  // parse error thing.
  void parseError(const char *const str = NULL);

public:

  enum TokenType {
    Token_EOF=0,
    Token_Integer=1,
    Token_Real=2,
    Token_Colon=3,
    Token_SemiColon=4,
    Token_LeftBrace=5,
    Token_RightBrace=6,
    Token_LeftParen=7,
    Token_RightParen=8,
    Token_VirtBar=9,
    Token_Keyword=10,
    Token_Identifier=11,
    Token_Undefined=12
  };

  struct TokenInfo {
    // the line in the file of the most recently processed token
    int srcLine; 
    // the char on the current line of most rec. proc. token
    int srcChar; 

    // the previously returned yylex return code
    int rc;

    // the token in string form
    char *tokenStr; 

    // the type of the token
    TokenType tokenType;

    ///////////////////////////////////////////
    // if of the appropriate type, the value
    // of the read token.
    union {
      int int_val;
      float float_val;
      double doub_val;
    };

    ///////////////////////////////////////////////////////////////
    // a simple string equality check
    bool operator == (const char *const s) { return !::strcmp(s,tokenStr); }
    bool operator != (const char *const s) { return ::strcmp(s,tokenStr); }

    ///////////////////////////////////////////////////////////////
    // a simple token type equality check
    bool operator == (const enum TokenType t) { return t == tokenType; }
    bool operator != (const enum TokenType t) { return t != tokenType; }

  };

  /////////////////////////////////////
  // filled in automatically by the scanner.
  static TokenInfo tokenInfo;

  void setTokenPosition();


  //////////////////////////////////////////////
  //  FileParser(const char *const fileName,
  //	     GMTK_GM& gm);

  //////////////////////////////////////////////
  // constructor opens and parses file, or dies if an
  // error occurs.
  FileParser(const char *const fileName);

  void parseGraphicalModel();
  void parseFrameList();
  void parseFrame();
  void parseRandomVariableList();
  void parseRandomVariable();
  void parseRandomVariableAttributeList();
  void parseRandomVariableAttribute();
  void parseRandomVariableType();
  void parseRandomVariableDisposition();
  void parseSwitchingParentList();
  void parseParentList();
  void parseMappingSpec();
  void parseConditionalParentListList();
  void parseConditionalParentList();
  void parseCPT_SPEC();
  void parseCPT_TYPE();
  void parseParent();
  void parseContinousImplementation();
  void parseChunkSpecifier();

};

#endif
