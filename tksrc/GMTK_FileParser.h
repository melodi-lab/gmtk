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

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_GM.h" 
#include "GMTK_CPT.h"
#include "GMTK_MixGaussiansCommon.h"

class RandomVariable;

class FileParser
{


private:

  ///////////////////////////////////////////////////
  typedef pair<string,int> rvParent;

  ///////////////////////////////////////////////////
  // A class where information about a RV can
  // be placed as it is being parsed. 
  class RVInfo {
    friend class FileParser;


    ////////////////////////////////////////////////////////////
    // define a bunch of types that are used in RVs
    enum Type { t_discrete, t_continuous, t_unknown };
    enum Disposition { d_hidden, d_observed, d_unknown };
    struct FeatureRange {
      bool filled;
      unsigned firstFeatureElement;
      unsigned lastFeatureElement;
      FeatureRange() { filled = false; }
      void clear() { filled = false; }
    };
    struct ListIndex {
      enum ListIndexType { li_String, li_Index, li_Unknown } liType;
      // if this is an integer index, this is used.
      unsigned intIndex;
      // Otherwise, if this is an string index, this is used.
      string nameIndex;
      void clear() { liType = li_Unknown; }
    };

    ///////////////////////////////////////////////////////////
    // data associated with a RV

    // the frame where it was defined (so is part of name really)
    unsigned frame;
    // line number of file where this RV was first declared
    unsigned fileLineNumber;
    // rv's name
    string name;    

    // fields from type
    Type rvType;
    Disposition rvDisp;
    unsigned rvCard;
    // if it is an observed variable, it must have a feature range.
    FeatureRange rvFeatureRange;

    // switching parents stuff
    vector< rvParent > switchingParents;
    ListIndex switchMapping;

    // conditional parents stuff
    vector<vector< rvParent > > conditionalParents;
    // if discrete, then the list if discrete implementations
    vector< CPT::DiscreteImplementaton > discImplementations;
    // if continuous, the list of continuous implementations
    vector< MixGaussiansCommon::ContinuousImplementation > contImplementations;
    // in either case, a low-level parameter index
    vector< ListIndex > listIndices;

    /////////////////////////////////////////////////////////
    // An actual pointer to the RV once we instantiate it. 
    RandomVariable* rv;

  public:
    /////////////////////////////////////////////////////////
    // constructor
    RVInfo() { rvType = t_unknown; rv = NULL; };
    // copy constructor
    RVInfo(const RVInfo&);

    // clear out the current RV structure when we 
    // are parsing and encounter a new RV.
    void clear() { 
      name.erase();

      rvType = t_unknown;
      rvDisp = d_unknown;
      rvFeatureRange.clear();

      switchingParents.clear();
      switchMapping.clear();

      conditionalParents.clear();
      discImplementations.clear();
      contImplementations.clear();
      listIndices.clear();

    }

    // checks the contents of the object to ensure
    // that this is a validly instantiated RV (i.e.,
    // all the information is there, it all makes sense, etc.)
    // Die if there is a problem.
    void checkConsistency();

  };

  ////////////////////////////////////////////////////////////////
  // The current pre-allocated random variable that is being
  // parsed and filled in as we go. 
  RVInfo curRV;


  ///////////////////////////////////
  // the current frame we are parsing.
  int curFrame;

  ///////////////////////////////////////////////////
  // Mapping from the name of the random variable
  // to its pointer.
  // map<pair<string, unsigned>, RandomVariable *> variableNamed;
  map < rvParent , unsigned > nameRVmap;

  //////////////////////////////////////////////
  // This is where the parser puts partially 
  // completed RVs as it is parsing them.
  vector < RVInfo > rvInfoVector;

  //////////////////////////////////////////////
  // the result of the chunk parse
  unsigned _firstChunkframe;
  unsigned _lastChunkframe;
  unsigned _maxFrame;  

public:

  //////////////////////////////////////////////
  //  FileParser(const char *const fileName,
  //	     GMTK_GM& gm);

  // the different type of tokens
  enum TokenType {
    TT_EOF=0,
    TT_Integer=1,
    TT_Real=2,
    TT_Colon=3,
    TT_SemiColon=4,
    TT_LeftBrace=5,
    TT_RightBrace=6,
    TT_LeftParen=7,
    TT_RightParen=8,
    TT_VirtBar=9,
    TT_Keyword=10,
    TT_Identifier=11,
    TT_Comma=12,
    TT_String=13,
    TT_Undefined=15
  };

  // when token type is a keyword, the different types of keywords
  enum TokenKeyword {
    KW_Frame = 0,
    KW_Variable = 1,
    KW_Type=2,
    KW_Cardinality=3,
    KW_Switchingparents=4,
    KW_Conditionalparents=5,
    KW_Discrete=6,
    KW_Continuous=7,
    KW_Hidden=8,
    KW_Observed=9,
    KW_Nil=10,
    KW_Using=11,
    KW_Mapping=12,
    KW_MDCPT=13,
    KW_MSCPT=14,
    KW_MTCPT=15,
    KW_MixGaussian=16,
    KW_GausSwitchMixGaussian=17,
    KW_LogitSwitchMixGaussian=18,
    KW_MlpSwitchMixGaussin=19,
    KW_Chunk=20,
    KW_GRAPHICAL_MODEL=21
  };

  // list of token keyword strings.
  static const vector < string >  KeywordTable;
  // the function that fills this in.
  static const vector < string >  fillKeywordTable();


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
    bool operator == (const string& s) { return (s == tokenStr); }
    bool operator != (const string& s) { return (s != tokenStr); }

    ///////////////////////////////////////////////////////////////
    // a simple token type equality check
    bool operator == (const enum TokenType t) { return t == tokenType; }
    bool operator != (const enum TokenType t) { return t != tokenType; }

    ///////////////////////////////////////////////////////////////
    // a simple token keyword equality check, these are
    // true only if both 1) the token is a keyword and the token
    // keyword equality holds true.

    bool operator == (const enum TokenKeyword kw) { 
      return ((*this == TT_Keyword) && (*this == KeywordTable[kw]));
    }
    bool operator != (const enum TokenKeyword kw) { 
      return !(*this == kw);
    }

  };


private:

  //  GMTK_GM& gm; 

  //////////////////////////////////
  // The file we are currently parsing.
  FILE *ifile;

  //////////////////////////////////////////////
  // prepares the next input token and fills the
  // structure
  void prepareNextToken();

  //////////////////////////////////////////////
  // prepares the next input token and fills the
  // structure
  void consumeToken();

  // returns true if we're at EOF
  bool endOfFile() { return (tokenInfo.rc == 0); }
  //////////////////////////////////////////////
  // ensures that the next position in the file
  // is not the EOF. It dies with
  // a parse error message if we are.
  void ensureNotAtEOF(const char *const str = NULL);
  // special parse error for key words
  void ensureNotAtEOF(const TokenKeyword kw);

  // parse error thing.
  void parseError(const char *const str = NULL);
  // special parse error for key words
  void parseError(const TokenKeyword kw);

  ////////////////////////////////////////////////////////////
  // the actual routines for the recursive descent parser 
  void parseFrameList();
  void parseFrame();
  void parseRandomVariableList();
  void parseRandomVariable();
  void parseRandomVariableAttributeList();
  void parseRandomVariableAttribute();
  void parseRandomVariableTypeAttribute();
  void parseRandomVariableType();
  void parseRandomVariableDiscreteType();
  void parseRandomVariableContinuousType();
  void parseRandomVariableParentAttribute();
  void parseSwitchingParentAttribute();
  void parseConditionalParentSpecList();
  void parseConditionalParentSpec();
  void parseConditionalParentList();

  void parseImplementation();
  void parseDiscreteImplementation();
  void parseContinuousImplementation();
  void parseContObsDistType();

  void parseChunkSpecifier();

  ///////////////////////////////////////////////////////////
  // parsing routines that are used by multiple
  // other parsing routines and therefore need their
  // own location to place data.

  void parseParentList();
  void parseParent();
  vector < rvParent > parentList;

  void parseMappingSpec();
  void parseListIndex();
  RVInfo::ListIndex listIndex;

public:

  /////////////////////////////////////////////
  // filled in automatically by the scanner.
  // This member acts as a lookahead token, to
  // the next "unconsumed" token.
  static TokenInfo tokenInfo;

  //////////////////////////////////////////////
  // constructor opens and parses file, or dies if an
  // error occurs.
  FileParser(const char *const fileName);
  ~FileParser();
  void parseGraphicalModel();
  void createRandomVariableGraph();
  void associateWithDataParams(bool allocateIfNotThere=false);
  // ensure links are "south", "south east",
  // "east", or "north east", meaning that there
  // is a numeric ordering on the nodes such that 
  // any parents of a node at a particular
  // numeric position have their position
  // earlier in the ordering.
  void ensureS_SE_E_NE();

  // add all the variables inself to the GM. Presumably
  // gm has not been added to yet.
  void addVariablesToGM(GMTK_GM& gm);

  //////////////////////////////////////////////////////////////
  // access to the chunk information.
  unsigned firstChunkFrame() { return _firstChunkframe; }
  unsigned lastChunkFrame() { return _firstChunkframe; }
  unsigned maxFrame() { return _maxFrame; }
  unsigned numFrames() { return _maxFrame+1; }
};

#endif
