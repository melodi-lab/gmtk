/*
 * GMTK_FileParser.h
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
 *
 * $Header$
 *
 */

#ifndef GMTK_FILEPARSER_H
#define GMTK_FILEPARSER_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_CPT.h"
// #include "GMTK_GMTemplate.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_GraphicalModel.h"
#include "GMTK_RVInfo.h"

#include "fileParser.h"

class RV;

class FileParser
{
 private:
  friend class RV;
  friend class StructPage;

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
  map < RVInfo::rvParent , unsigned > nameRVmap;

  //////////////////////////////////////////////
  // This is where the parser puts partially 
  // completed RVs as it is parsing them. This array
  // will be filled in in order of variables encountered
  // in the file. This means that all frame 0 variables
  // will be seen (and positioned in the array) first, 
  // then frame 1, and so on.
  vector < RVInfo > rvInfoVector;

  //////////////////////////////////////////////
  // the result of the chunk parse
  unsigned _firstChunkframe;
  unsigned _lastChunkframe;
  unsigned _maxFrame;  

  ////////////////////////////////////////////////
  // more information filled in after the parse
  unsigned numVarsInPrologue;
  unsigned numVarsInChunk;
  unsigned numVarsInEpilogue;

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
    TT_multiLineString=14,
    TT_Undefined=15
  };


  // ******************************************************************
  // This enum must always be order consistent with the table 
  // in fillKeywordTable in the .cc file, and consistent with the
  // list of keywords in the .lex file.
  // ******************************************************************
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
    KW_Collection=13,
    KW_MDCPT=14,
    KW_MSCPT=15,
    KW_MTCPT=16,
    KW_NGRAMCPT=17,
    KW_FNGRAMCPT=18,
    KW_Mixture=19,
    KW_GausSwitchMixture=20,
    KW_LogitSwitchMixture=21,
    KW_MlpSwitchMixture=22,
    KW_Chunk=23,
    KW_GRAPHICAL_MODEL=24,
    KW_Value=25,
    KW_Weight=26,
    KW_Scale=27,
    KW_Penalty=28,
    KW_Shift=29,
    KW_EliminationHint=30,
    KW_FrameNum=31,
    KW_NumFrames=32,
    KW_SegmentNum=33,
    KW_NumSegments=34,
    KW_VECPT=35
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
      // float float_val;
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
  void parseErrorExpecting(const char *const str = NULL);
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

  void parseRandomVariableWeightAttribute();
  void parseRandomVariableWeightAttributeSpecList();
  void parseRandomVariableWeightAttributeSpec();
  void parseRandomVariableWeightOptionList();

  void parseRandomVariableEliminationHintAttribute();
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
  vector < RVInfo::rvParent > parentList;

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
  static string fileNameParsing;
  FileParser(const char *const fileName, 
	     const char *const cppCommandOptions = NULL);
  ~FileParser();
  void parseGraphicalModel();
  void createRandomVariableGraph();
  // ensure no loops in graph for all possible unrollings.
  void ensureValidTemplate();

  enum MdcptAllocStatus { noAllocate, allocateRandom, allocateUniform };
  void associateWithDataParams(MdcptAllocStatus allocate = noAllocate);
  // ensure links are "south", "south east",
  // "east", or "north east", meaning that there
  // is a numeric ordering on the nodes such that 
  // any parents of a node at a particular
  // numeric position have their position
  // earlier in the ordering.
  void ensureS_SE_E_NE();

  // this function checks to make sure that variable parents do not
  // span multiple regions.  I.e., variables in epilogue can't have
  // parents in the prologue, parents in prolog can't have (future)
  // parents in the epilogue. Note that in some cases, this would be
  // valid, but it would break some constrained triangulation schemes.
  // For example, if the network is fully unrolled for each T and then
  // triangulated, this routine is not required.
  void ensureVariablesDoNotReachAcrossRegion() {}

  // add all the variables to a template, essentially
  // keeping all the rv information but removing the
  // file specific information.
  // void addVariablesToTemplate(GMTemplate&);

  void checkConsistentWithGlobalObservationStream();

  // unroll the template chunk k times, and place
  // the result in the existing vector of random
  // variables node.
  void unroll(unsigned k,vector<RV*> &unrolledVarSet);
  void unroll(unsigned k,
	      vector<RV*> &unrolledVarSet,
	      map < RVInfo::rvParent, unsigned >& ppf);

  // A routine to write out the graph template (P,C,E) in condensed
  // form but in sufficient detail so that it can be used to quickly
  // ID the current template (e.g., so that a given elimination
  // order can be checked with a current graph).
  void writeGMId(oDataStreamFile& os);

  // A routine that reads in a graph id that was written
  // in condensed by writeGMId(), and it verifies that
  // the GM ID written matches the current template.
  bool readAndVerifyGMId(iDataStreamFile& is);


  //////////////////////////////////////////////////////////////
  // access to the chunk information.
  unsigned firstChunkFrame() { return _firstChunkframe; }
  unsigned lastChunkFrame() { return _lastChunkframe; }
  unsigned maxFrame() { return _maxFrame; }
  unsigned numFrames() { return _maxFrame+1; }
  // number of frames in prologue
  unsigned numFramesInP() { return _firstChunkframe; }
  // number of frames in chunk
  unsigned numFramesInC() { return _lastChunkframe - _firstChunkframe +1; }
  // number of frames in epilogue  
  unsigned numFramesInE() { return _maxFrame - _lastChunkframe; }
};

#endif
