/*-
 * GMTK_FileParser.cc
 *     structure file parsing
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_USCPT.h"
#include "GMTK_MixGaussians.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_GraphicalModel.h"
#include "GMTK_RVInfo.h"

VCID("$Header$");

#ifndef DECLARE_POPEN_FUNCTIONS_EXTERN_C
extern "C" {
  FILE     *popen(const char *, const char *) __THROW;
  int pclose(FILE *stream) __THROW;
};
#endif

/*
***********************************************************************
***********************************************************************

The GM Grammar: 
        This is a valid grammer for the parser below.
	Please try to keep this grammer up to date if any
	changes are made to the parser.

GM = "GRAPHICAL_MODEL" identifier FrameList ChunkSpecifier

FrameList = Frame FrameList | NULL

Frame = "frame" ":" integer "{" RandomVariableList "}"

RandomVariableList = RandomVariable RandomVariableList | NULL

RV = "variable" ":" name "{" RandomVariableAttribute "}"

RandomVariableAttributeList =
        RandomVariableAttribute RandomVariableAttributeList | NULL

RandomVariableAttribute = TypeAttribute | ParentsAttribute

TypeAttribute = "type" ":" RandomVariableType ";"

ParentsAttribute =
      ( "switchingparents" ":" SwitchingParentAttribute ";" )
   |  ( "conditionalparents" ":" ConditionalParentSpecList  ";" )
    # Semanatics requires that we have both
    # switchingparents & conditionalparents in a RV. Note that
    # the parse grammer allows this not to be the case.

RandomVariableType = RandomVariableDiscreteType | 
                     RandomVariableContinuousType

RandomVariableDiscreteType = 
      "discrete" 
      ( "hidden" | "observed" (integer ":" integer | "value" integer) ) 
     "cardinality" integer

RandomVariableContinuousType = 
       "continuous" 
       ("hidden" | "observed" integer ":" integer)

SwitchingParentAttribute = "nil" | ParentList "using" MappingSpec

ConditionalParentSpecList =
    ConditionalParentSpec "|" ConditionalParentSpecList
  | ConditionalParentSpec

ConditionalParentSpec = ConditionalParentList using Implementation

ConditionalParentList = "nil" | ParentList

ParentList = 
       Parent "," ParentList
     | Parent

Parent = identifier "(" integer ")"

Implementation = DiscreteImplementation | ContinuousImplementation

DiscreteImplementation = ( "DenseCPT" | "SparseCPT" | "DeterministicCPT" )  "(" ListIndex ")"

ContinuousImplementation = ContObsDistType
        (
             "(" ListIndex ")"
           |
              "collection" "(" ListIndex ")" MappingSpec
        )
       #
       # A ContinousImplementation has two cases:
       # 1) the first case uses the syntax
       #
       #          "(" ListIndex ")"
       #
       # and is used if conditional parents 
       # are nil in which case we select only one dist. In 
       # this case, 'ListIndex' refers directly to 
       # a particular single gaussian mixture, using
       # the normal 'ListIndex' format (either an int
       # or a string name). Note that in this case, 
       # the number of conditional parents (as given
       # in the 'conditionalParents' array, will be zero.
       # This is analogous to the discrete case where you directly
       # select a CPT with the appropriate parents
       #
       # 2) in the second case, which uses the syntax
       #
       #           "collection" "(" ListIndex ")" MappingSpec
       #
       # this is when we have multiple
       # conditional parents, and we need another decision tree to map
       # from the conditional parents values to the appropriate
       # distribution. Therefore we use the mapping syntax.
       # In this case, the 'MappingSpec' must refer to one of
       # the decision trees.  The collection referes to the 
       # collection of objects in global arrays, and the mappingspec
       # gives the mapping from random variable parent indices to
       # the relative offset in the collection array.
       # 


ContObsDistType = "mixGaussian" | "gausSwitchMixGaussian" 
  | "logitSwitchMixGaussian" | "mlpSwitchMixGaussian"

MappingSpec = "mapping" "(" ListIndex ")"
     # A MappingSpec always indexes into one of the decision trees.
     # The integer (or string) is used to index into a table
     # of decision trees to choose the decision tree
     # that will map from the switching parents to one of the
     # conditional parent lists.

ChunkSpecifier = "chunk"  integer ":" integer

ListIndex = integer | string

======================================================================


Example of a grammatical GM file
-----------------------------------

#
# An actual model definition (that parses).
# 
GRAPHICAL_MODEL FHMM

frame:0 {

     variable : word {
          type: discrete hidden cardinality 3 ;
          switchingparents: nil ;
          conditionalparents: nil using DenseCPT("initcpt") ;
        }
     variable : phone {
          type: discrete hidden cardinality 4 ;
          switchingparents : nil ;
          conditionalparents :
                     word(0) using DenseCPT("dep_on_word") ;
        }
     variable : state1 {
          type: discrete hidden cardinality 4 ;
          switchingparents: phone(0) using mapping("phone2state1") ;
          conditionalparents :
                  nil using SparseCPT("f1")
                | nil using DenseCPT("f2") ;
       }
       variable : state2 {
          type: discrete hidden cardinality 4 ;
          switchingparents: nil ;
          conditionalparents:
                     nil using DenseCPT("f5") ;
       }

       variable : obs1 {
          type: continuous observed 0:5 ;
          switchingparents: state1(0), state2(0)
                   using mapping("state2obs6") ;
          conditionalparents: 
                 nil using mixGaussian("the_forth_gaussian")
               | state1(0) using mixGaussian collection("gaussianCollection") mapping("gausmapping");
       }
       variable : obs2 {
          type: continuous observed 6:25  ;
          switchingparents: nil ;
          conditionalparents: 
                state1(0) using mlpSwitchMixGaussian collection("mlpgaussianCollection")
                          mapping("gaussmapping3") ;
       }
}

frame:1 {

     variable : word {
          type: discrete hidden cardinality 3 ;
          switchingparents: nil ;
          conditionalparents: word(-1) using DenseCPT("wordbigram") ;
        }

     variable : phone {
          type: discrete hidden cardinality 4 ;
          switchingparents : nil ;
          conditionalparents :
                phone(-1),word(0) using DenseCPT("dep_on_word_phone") ;
        }

     variable : state1 {
          type: discrete hidden cardinality 4 ;
          switchingparents: phone(0) using mapping("phone2state1") ;
          conditionalparents :
	          # in this first case, state1 is dep on prev time
                  state1(-1) using SparseCPT("f3")
	          # in this second case, it is cond. indep. of prev. time.
                | nil using DenseCPT("f2") ;
       }

       variable : state2 {
          type: discrete hidden cardinality 4 ;
          switchingparents: nil ;
          conditionalparents:
	             # this is a sep. hidden markov chain
                     state2(-1) using DenseCPT("f9") ;
       }

       # like in the first frame, the obs. are only dep.
       # on RVs from the current frame.
       variable : obs1 {
          type: continuous observed 0:5 ;
          switchingparents: state1(0), state2(0) 
                   using mapping("state2obs6") ;
          conditionalparents: 
                 nil using mixGaussian("the_forth_gaussian")
               | state1(0) using mixGaussian collection("gaussianCollection") mapping("gausmapping");
       }

       variable : obs2 {
          type: continuous observed 6:25  ;
          switchingparents: nil ;
          conditionalparents: 
                state1(0) using mlpSwitchMixGaussian collection("gaussianCollection")
                          mapping("gaussmapping3") ;
       }
}

chunk: 1:1

------------------------------------------------------------

*********************************************************************** 
*********************************************************************** 
*/


////////////////////////////////////////////////////////////////////
//        lex declarations
////////////////////////////////////////////////////////////////////

extern FILE *yyin, *yyout;
extern int yylex();


////////////////////////////////////////////////////////////////////
//        static data declarations
////////////////////////////////////////////////////////////////////

FileParser::TokenInfo FileParser::tokenInfo;
string FileParser::fileNameParsing;

////////////////////////////////////////////////////////////////////
//        Internal Keyword Symbols
////////////////////////////////////////////////////////////////////


const vector<string> 
FileParser::KeywordTable = FileParser::fillKeywordTable();

const vector<string> 
FileParser::fillKeywordTable()
{
  // ******************************************************************
  // This table must always be consistent with the enum TokenKeyword in 
  // the .h file and with 'keyword' in the .lex file.
  // ******************************************************************
  const char*const kw_table[] = {
    "frame",
    "variable",
    "type",
    "cardinality",
    "switchingparents",
    "conditionalparents",
    "discrete",
    "continuous",
    "hidden",
    "observed",
    "nil",
    "using",
    "mapping",
    "collection",
    "DenseCPT",
    "SparseCPT",
    "DeterministicCPT",
    "mixGaussian",
    "gausSwitchMixGaussian",
    "logitSwitchMixGaussian",
    "mlpSwitchMixGaussian",
    "chunk",
    "GRAPHICAL_MODEL",
    "value",
    "weight",
    "observation" 
  };
  vector<string> v;
  const unsigned len = sizeof(kw_table)/sizeof(char*);
  v.resize(len);
  for (unsigned i=0;i<len;i++) {
    v[i] = kw_table[i];
  }
  return v;
}


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * FileParser
 *      constructor of an object, takes file argument, program
 *      dies if there is an error or file does not parse.
 *
 *      *******************************************************
 *      NOTE: CAN ONLY PARSE ONE FILE AT A TIME since static
 *      members are used to hold current state of object.
 *      *******************************************************
 *
 * Preconditions:
 *      file should contain a valid GM structure file.
 *
 * Postconditions:
 *
 * Side Effects:
 *      program might die, totally changes internal structure of object.
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
FileParser::FileParser(const char*const file,
		       const char*const cppCommandOptions)
{
  FILE* f;
  if (file == NULL)
    error("FileParser::FileParser, can't open NULL file");
  string cppCommand = string("cpp");
  if (cppCommandOptions != NULL)
    cppCommand = cppCommand + string(" ") + cppCommandOptions;
  if (!strcmp("-",file)) {
    f = ::popen(cppCommand.c_str(),"r");
    if (f == NULL) {
      error("ERROR: unable to open with standard input structure file");
    }
  }  else {
    if ((f = ::fopen(file,"r")) == NULL) {
      error("ERROR: unable to open file (%s) for reading",file);
    }
    fclose(f);
    cppCommand = cppCommand + string(" ") + (string)file;
    f = ::popen(cppCommand.c_str(),"r");
    if (f == NULL)
      error("FileParser::FileParser, can't open file stream from (%s)",file);
  }
  yyin = f;
  fileNameParsing = file;
}



/*-
 *-----------------------------------------------------------------------
 * ~FileParser
 *      destructor
 *
 * Preconditions:
 *      should be constructed.
 *
 * Postconditions:
 *
 * Side Effects:
 *      Just closes the input file.
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
FileParser::~FileParser()
{
  pclose(yyin);
}




/*-
 *-----------------------------------------------------------------------
 * prepareNextToken
 *   prepares the next input token and fills the tokeninfo structure.
 *   Note that this depends on lex. 
 * 
 * Preconditions:
 *      Object must be in the midst of parsing a file.
 *
 * Postconditions:
 *      Same as before, new token parsed, possible EOF condition.
 *
 * Side Effects:
 *      Changes the interal token info function.
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
void
FileParser::prepareNextToken()
{
  tokenInfo.rc = yylex();
}



/*-
 *-----------------------------------------------------------------------
 * ConsumeToken
 *   consumes the current token 
 * 
 * Preconditions:
 *      Object must be in the midst of parsing a file.
 *
 * Postconditions:
 *      Same as before, new token parsed, possible new 
 *      EOF condition might exist.
 *
 * Side Effects:
 *      Changes the interal token info function to point
 *      to the next unconsumed token (or EOF if that occurs).
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
void
FileParser::consumeToken()
{
#if 0  
  printf("Consuming token (%s) of type %d from src line (%d)\n",
	 tokenInfo.tokenStr,tokenInfo.tokenType,tokenInfo.srcLine);
#endif
  prepareNextToken();
}



/*-
 *-----------------------------------------------------------------------
 * ensureNotAtEOF:
 *   makes sure that we are not at the EOF.
 * 
 * Preconditions:
 *      Object must be in the midst of parsing a file.
 *
 * Postconditions:
 *      same as before.
 *
 * Side Effects:
 *      program might die as result of error.
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
void
FileParser::ensureNotAtEOF(const char *const msg)
{
  if (tokenInfo.rc == TT_EOF) {
    fprintf(stderr,"Unexpected EOF Error: expecting %s at line %d\n",
	    msg,
	    tokenInfo.srcLine);
    error("Exiting Program");
  }
}
void
FileParser::ensureNotAtEOF(const TokenKeyword kw)
{
  if (tokenInfo.rc == TT_EOF) {
    fprintf(stderr,"Unexpected EOF Error: expecting keyword (%s) at line %d\n",
	    KeywordTable[kw].c_str(),
	    tokenInfo.srcLine);
    error("Exiting Program");
  }
}




/*-
 *-----------------------------------------------------------------------
 * parseError
 *   exists with a parse error condition.
 * 
 * Preconditions:
 *      some error
 *
 * Postconditions:
 *      program dies.
 *
 * Side Effects:
 *      program dies.
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
void
FileParser::parseError(const char* const str)
{
  fprintf(stderr,"Parse Error in file '%s': %s at or before line %d, near (%s)\n",
	  fileNameParsing.c_str(),
	  str,
	  tokenInfo.srcLine,
	  tokenInfo.tokenStr);
  error("Exiting Program");
}

void
FileParser::parseError(const TokenKeyword kw)
{
  fprintf(stderr,"Parse Error in file '%s': expecting keyword (%s) at or before line %d, near (%s)\n",
	  fileNameParsing.c_str(),
	  KeywordTable[kw].c_str(),
	  tokenInfo.srcLine,
	  tokenInfo.tokenStr);
  error("Exiting Program");
}




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//                  PARSING FUNCTIONS                             //
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


void
FileParser::parseGraphicalModel()
{
  tokenInfo.srcLine = 1;
  tokenInfo.srcChar = 1;
  curFrame = -1;
  rvInfoVector.clear();

  // prepare the first lookahead token
  prepareNextToken();

  // now start parsing.
  ensureNotAtEOF("GM magic keyword");
  if (tokenInfo != KW_GRAPHICAL_MODEL)
    parseError(KW_GRAPHICAL_MODEL);
  consumeToken();

  ensureNotAtEOF("GM name");
  if (tokenInfo != TT_Identifier)
    parseError("expecting GM name identifier");
  consumeToken();

  // look ahead to make sure it is a frame
  ensureNotAtEOF("frame keyword");
  if (tokenInfo != KW_Frame)
    parseError(KW_Frame);
  // now this will consume the token
  parseFrameList();
  _maxFrame = curFrame;

  parseChunkSpecifier();

  if (_lastChunkframe > _maxFrame)
    parseError("last chunk integer must be no greater than number of frames");

#if 0
  printf("Found %d random variables.\n",
	 rvInfoVector.size());
  for (unsigned i=0;i<rvInfoVector.size();i++ ) {
    printf("RV(%d) is named (%s)\n",
	   i,rvInfoVector[i].name.c_str());
  }
#endif


}


void
FileParser::parseFrameList()
{
  if (tokenInfo != KW_Frame)
    return;
  parseFrame();
  parseFrameList();
}


void
FileParser::parseFrame()
{
  ensureNotAtEOF(KW_Frame);
  if (tokenInfo != KW_Frame)
    parseError(KW_Frame);
  consumeToken();

  ensureNotAtEOF("frame colon");
  if (tokenInfo != TT_Colon)
    parseError("frame colon");
  consumeToken();

  ensureNotAtEOF("frame number");
  if (tokenInfo != TT_Integer)
    parseError("frame number");
  if (tokenInfo.int_val != (curFrame+1)) {
    parseError("frame number out of order");
  }
  curFrame++;
  consumeToken();


  ensureNotAtEOF("open frame {");
  if (tokenInfo != TT_LeftBrace)
    parseError("open frame {");
  consumeToken();

  parseRandomVariableList();

  ensureNotAtEOF("close frame }");
  if (tokenInfo != TT_RightBrace)
    parseError("close frame }");
  consumeToken();
}

void
FileParser::parseRandomVariableList()
{
  if (tokenInfo != KW_Variable)
    return;

  curRV.clear();
  parseRandomVariable();

  //////////////////////////////
  // make sure that there was a type and switching/cont parents
  if (curRV.rvType == RVInfo::t_unknown)
    error("Random variable type unknown for variable %s at frame %d, line %d\n",
	  curRV.name.c_str(),curRV.frame,curRV.fileLineNumber);

  // we need to at least have specified one conditional parents set,
  // which might be 'nil' 
  if (curRV.conditionalParents.size() == 0)  
    error("Conditional parents unknown for random variable %s at frame %d, line %d\n",
	  curRV.name.c_str(),curRV.frame,curRV.fileLineNumber);


  // check if we've already seen this RV
  map < RVInfo::rvParent , unsigned >::iterator it;
  it = nameRVmap.find(RVInfo::rvParent(curRV.name,curRV.frame));
  if (it != nameRVmap.end()) {
    error("Error: random variable (%s) at frame (%d) defined twice, "
	  "on both line %d and %d\n",curRV.name.c_str(),
	  curRV.frame,
	  curRV.fileLineNumber,
	  rvInfoVector[(*it).second].fileLineNumber);
  }
  
  curRV.variablePositionInStrFile = rvInfoVector.size();

  // everything looks ok, insert it in our tables.
  rvInfoVector.push_back(curRV);
  // add to map
  nameRVmap[
	    RVInfo::rvParent(curRV.name,curRV.frame)
            ] 
    = rvInfoVector.size()-1;

  parseRandomVariableList();
}


void
FileParser::parseRandomVariable()
{
  ensureNotAtEOF(KW_Variable);
  if (tokenInfo != KW_Variable)
    parseError(KW_Variable);
  curRV.frame = curFrame;
  curRV.fileLineNumber = tokenInfo.srcLine;
  curRV.rvFileName  = fileNameParsing;
  consumeToken();


  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseError(":");
  consumeToken();

  ensureNotAtEOF("variable name");
  if (tokenInfo != TT_Identifier)
    parseError("variable name");
  curRV.name = tokenInfo.tokenStr;
  consumeToken();

  ensureNotAtEOF("open RV {");
  if (tokenInfo != TT_LeftBrace)
    parseError("open RV {");
  consumeToken();

  parseRandomVariableAttributeList();

  ensureNotAtEOF("close frame }");
  if (tokenInfo != TT_RightBrace)
    parseError("close RV }");
  consumeToken();

}

void
FileParser::parseRandomVariableAttributeList()
{
  if (tokenInfo == KW_Type || 
      tokenInfo == KW_Weight ||
      tokenInfo == KW_Switchingparents ||
      tokenInfo == KW_Conditionalparents) 
    {
      parseRandomVariableAttribute();
      parseRandomVariableAttributeList();
    }
    return;
}

void
FileParser::parseRandomVariableAttribute()
{
  ensureNotAtEOF("variable attribute");
  if (tokenInfo == KW_Type)
    return parseRandomVariableTypeAttribute();
  else if (tokenInfo == KW_Weight)
    return parseRandomVariableWeightAttribute();
  else if (tokenInfo == KW_Switchingparents ||
	   tokenInfo == KW_Conditionalparents) {
    if (curRV.rvType == RVInfo::t_unknown)
      parseError("type must be first attribute");
    return parseRandomVariableParentAttribute();
  }  else
    parseError("variable attribute");
}

void
FileParser::parseRandomVariableTypeAttribute()
{

  ensureNotAtEOF(KW_Type);
  if (tokenInfo != KW_Type)
    parseError(KW_Type);
  consumeToken();

  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseError(":");
  consumeToken();

  parseRandomVariableType();

  ensureNotAtEOF(";");
  if (tokenInfo != TT_SemiColon)
    parseError(";");
  consumeToken();

}

void
FileParser::parseRandomVariableType()
{
  ensureNotAtEOF("variable type");
  if (tokenInfo == KW_Discrete)
    parseRandomVariableDiscreteType(); 
  else if (tokenInfo == KW_Continuous)
    parseRandomVariableContinuousType();
  else 
    parseError("variable type");
}



void
FileParser::parseRandomVariableDiscreteType()
{
  ensureNotAtEOF(KW_Discrete);
  if (tokenInfo != KW_Discrete)
    parseError(KW_Discrete);
  if (curRV.rvType != RVInfo::t_unknown)
    error("Error parsing file: RV (%s) already has a type, frame %d, line %d",
	  curRV.name.c_str(),
	  curRV.frame,
	  curRV.fileLineNumber);
  curRV.rvType = RVInfo::t_discrete;
  consumeToken();

  ensureNotAtEOF("variable disposition (hidden|discrete)");
  if (tokenInfo == KW_Hidden) {
    // hidden discrete RV
    curRV.rvDisp = RVInfo::d_hidden;
    consumeToken();
  } else if (tokenInfo == KW_Observed) {
    // observed discrete RV    
    curRV.rvDisp = RVInfo::d_observed;
    consumeToken();
    
    ensureNotAtEOF("first feature range|value keyword");
    if (tokenInfo == TT_Integer) {
      // should be n:m syntax, so values come from file

      if (tokenInfo != TT_Integer)
	parseError("first feature range");
      curRV.rvFeatureRange.firstFeatureElement = tokenInfo.int_val;
      consumeToken();

      ensureNotAtEOF("feature range separator");
      if (tokenInfo != TT_Colon)
	parseError("feature range separator");
      consumeToken();
    
      ensureNotAtEOF("second feature range");
      if (tokenInfo != TT_Integer)
	parseError("second feature range");
      curRV.rvFeatureRange.lastFeatureElement = tokenInfo.int_val;
      if (curRV.rvFeatureRange.lastFeatureElement < curRV.rvFeatureRange.firstFeatureElement)
	parseError("first range num must be < second range num");
      curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_Range;
      // A discrete random variable is, at this time, only a scalar.
      // Therefore, the feature range spec should be of the
      // form n:n (i.e., it may only specify a single number).
      // Ultimately, a discrete RV will be a vector so we 
      // keep the n:m notation, but for we check this
      if (curRV.rvFeatureRange.lastFeatureElement != curRV.rvFeatureRange.firstFeatureElement)
	parseError("for scalar RV, first feature range num be same as second range num");
      consumeToken();
    }  else if (tokenInfo == KW_Value) {
      // should be "value n" syntax
      consumeToken(); // consume the 'value' token

      if (tokenInfo != TT_Integer)
	parseError("value integer");
      if (tokenInfo.int_val < 0)
	parseError("non-negative value integer");
      
      curRV.rvFeatureRange.firstFeatureElement = tokenInfo.int_val;
      curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_FirstIsValue;
      consumeToken();

    } else {
      // parse error
      parseError("first feature range n:m | value keyword");
    }
  } else 
    parseError("variable disposition (hidden|discrete)");

  ensureNotAtEOF(KW_Cardinality);
  if (tokenInfo != KW_Cardinality)
    parseError(KW_Cardinality);
  consumeToken();

  ensureNotAtEOF("cardinality value");
  if (tokenInfo != TT_Integer)
    parseError("cardinality value");
  curRV.rvCard = tokenInfo.int_val;
  if (curRV.rvCard == 0)
    parseError("positive cardinality value");
  consumeToken();
}


void
FileParser::parseRandomVariableContinuousType()
{
  ensureNotAtEOF(KW_Continuous);
  if (tokenInfo != KW_Continuous)
    parseError(KW_Continuous);
  if (curRV.rvType != RVInfo::t_unknown)
    error("Error parsing file: RV (%s) already has a type, frame %d, line %d",
	  curRV.name.c_str(),
	  curRV.frame,
	  curRV.fileLineNumber);
  curRV.rvType = RVInfo::t_continuous;
  consumeToken();

  ensureNotAtEOF("variable disposition (hidden|discrete)");
  if (tokenInfo == KW_Hidden) {
    // hidden continuous RV
    curRV.rvDisp = RVInfo::d_hidden;
    // ... fill in later
    consumeToken();
    
    ///////////////////////////////////////////////
    // TODO: removing the following is safe when
    // hidden continuous variables are implemented.
    parseError("sorry, hidden continuous RVs not yet implemented");

  } else if (tokenInfo == KW_Observed) {
    // observed continuous RV    
    curRV.rvDisp = RVInfo::d_observed;
    consumeToken();
    
    ensureNotAtEOF("first feature range");
    if (tokenInfo != TT_Integer)
      parseError("first feature range");
    curRV.rvFeatureRange.firstFeatureElement = tokenInfo.int_val;
    consumeToken();

    ensureNotAtEOF("feature range separator");
    if (tokenInfo != TT_Colon)
      parseError("feature range separator");
    consumeToken();
    
    ensureNotAtEOF("second feature range");
    if (tokenInfo != TT_Integer)
      parseError("second feature range");
    curRV.rvFeatureRange.lastFeatureElement = tokenInfo.int_val;
    if (curRV.rvFeatureRange.lastFeatureElement < curRV.rvFeatureRange.firstFeatureElement)
      parseError("first range num must be <= second range num");
    curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_Range;
    consumeToken();

  } else
    parseError("variable disposition (hidden|discrete)");

}

void
FileParser::parseRandomVariableWeightAttribute()
{


  ensureNotAtEOF(KW_Weight);
  if (tokenInfo != KW_Weight)
    parseError(KW_Weight);
  if (curRV.rvWeightInfo.wt_Status != RVInfo::WeightInfo::wt_NoWeight) {
    parseError("RV already has weight attribute");
  }
  consumeToken();

  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseError(":");
  consumeToken();

  ensureNotAtEOF("weight value or observation");
  if (tokenInfo == KW_Value) {
    consumeToken();

    ensureNotAtEOF("weight floating-point value");
    // allow an int to be treated as a float value.
    if (tokenInfo != TT_Real && tokenInfo != TT_Integer)
      parseError("weight floating-point value");
    if (tokenInfo == TT_Real)
      curRV.rvWeightInfo.weight_value = tokenInfo.doub_val;
    else 
      curRV.rvWeightInfo.weight_value = (double)tokenInfo.int_val;
    consumeToken();

    curRV.rvWeightInfo.wt_Status = RVInfo::WeightInfo::wt_Constant;    
    
  } else if (tokenInfo == KW_Observation) {
    consumeToken();
    
    ensureNotAtEOF("first feature range");
    if (tokenInfo != TT_Integer)
      parseError("first feature range");
    curRV.rvWeightInfo.firstFeatureElement = tokenInfo.int_val;
    consumeToken();

    ensureNotAtEOF("feature range separator");
    if (tokenInfo != TT_Colon)
      parseError("feature range separator");
    consumeToken();
    
    ensureNotAtEOF("second feature range");
    if (tokenInfo != TT_Integer)
      parseError("second feature range");
    curRV.rvWeightInfo.lastFeatureElement = tokenInfo.int_val;
    if (curRV.rvWeightInfo.lastFeatureElement != curRV.rvWeightInfo.firstFeatureElement
	|| (curRV.rvWeightInfo.lastFeatureElement < 0))
      parseError("first range num must be == second range num for observation scalar weight");
    consumeToken();

    curRV.rvWeightInfo.wt_Status = RVInfo::WeightInfo::wt_Observation;

  } else 
    parseError("weight value or observation");


  ensureNotAtEOF(";");
  if (tokenInfo != TT_SemiColon)
    parseError(";");
  consumeToken();

}



void
FileParser::parseRandomVariableParentAttribute()
{

  ensureNotAtEOF("parent attribute");
  if (tokenInfo == KW_Switchingparents) {
    consumeToken();

    ensureNotAtEOF("attribute seperator");
    if (tokenInfo != TT_Colon)
      parseError("attribute separator");
    consumeToken();

    parseSwitchingParentAttribute();

    ensureNotAtEOF(";");
    if (tokenInfo != TT_SemiColon)
      parseError(";");
    consumeToken();

  } else if (tokenInfo == KW_Conditionalparents) {
    consumeToken();

    ensureNotAtEOF("attribute seperator");
    if (tokenInfo != TT_Colon)
      parseError("attribute separator");
    consumeToken();

    parseConditionalParentSpecList(); 

    ensureNotAtEOF(";");
    if (tokenInfo != TT_SemiColon)
      parseError(";");
    consumeToken();

  } else 
    parseError("parent attribute");

}

void
FileParser::parseSwitchingParentAttribute()
{

  ensureNotAtEOF("list of switching parents");
  parentList.clear();
  if (tokenInfo == KW_Nil) {
    // set parent list to the empty array.
    curRV.switchingParents = parentList;
    consumeToken();
  } else {

    parseParentList();
    curRV.switchingParents = parentList;
    ensureNotAtEOF(KW_Using);
    if (tokenInfo != KW_Using)
      parseError(KW_Using);
    consumeToken();

    parseMappingSpec();
    curRV.switchMapping = listIndex;
  }
}

void
FileParser::parseConditionalParentSpecList()
{
  parseConditionalParentSpec();
  if (tokenInfo == TT_VirtBar) {
    consumeToken();
    parseConditionalParentSpecList();
  }
}

void
FileParser::parseConditionalParentSpec()
{
  parseConditionalParentList();
  curRV.conditionalParents.push_back(parentList);

  ensureNotAtEOF(KW_Using);
  if (tokenInfo != KW_Using)
    parseError(KW_Using);
  consumeToken();

  parseImplementation();



}

void
FileParser::parseConditionalParentList()
{
  parentList.clear();
  if (tokenInfo == KW_Nil) {
    consumeToken();
  } else {
    parseParentList();
  }
}

void
FileParser::parseParentList()
{
  parseParent();
  if (tokenInfo == TT_Comma) {
    consumeToken();
    parseParentList();
  }
}

void
FileParser::parseParent()
{
  RVInfo::rvParent p;

  ensureNotAtEOF("parent RV name");
  if (tokenInfo != TT_Identifier) 
    parseError("parent RV name");
  p.first = tokenInfo.tokenStr;
  consumeToken();

  ensureNotAtEOF("(");
  if (tokenInfo != TT_LeftParen) 
    parseError("(");
  consumeToken();

  ensureNotAtEOF("parent RV offset");
  if (tokenInfo != TT_Integer) 
    parseError("parent RV offset");
  p.second = tokenInfo.int_val;;
  // Simple topology check right here, make
  // sure we are not pointing to ourselves.
  if ((p.first == curRV.name) &&
      (p.second == 0)) {
    parseError("parent variable must not refer to self");
  }

  parentList.push_back(p);
  consumeToken();

  ensureNotAtEOF(")");
  if (tokenInfo != TT_RightParen) 
    parseError(")");
  consumeToken();

}

void
FileParser::parseImplementation()
{
  ensureNotAtEOF("implementation");  
  if (curRV.rvType == RVInfo::t_discrete)
    parseDiscreteImplementation();
  else
    parseContinuousImplementation();

#if 0

  if (tokenInfo == KW_MDCPT || tokenInfo == KW_MSCPT 
      || tokenInfo == KW_MTCPT) {
    if (curRV.rvType != RVInfo::t_discrete) 
      parseError("need discrete implementations in discrete RV");
    parseDiscreteImplementation();
  } else {
    if (curRV.rvType != RVInfo::t_continuous) 
      parseError("need continuous implementations in continuous RV");
    parseContinuousImplementation();
  }
#endif

}


void
FileParser::parseDiscreteImplementation()
{
  ensureNotAtEOF("discrete implementation");  
  if (tokenInfo == KW_MDCPT || tokenInfo == KW_MSCPT
      || tokenInfo == KW_MTCPT) {

    if (tokenInfo == KW_MDCPT)
      curRV.discImplementations.push_back(CPT::di_MDCPT);
    else if (tokenInfo == KW_MSCPT)
      curRV.discImplementations.push_back(CPT::di_MSCPT);
    else // (tokenInfo == KW_MTCPT)
      curRV.discImplementations.push_back(CPT::di_MTCPT);
    consumeToken();


    ensureNotAtEOF("(");
    if (tokenInfo != TT_LeftParen) {
      parseError("(");
    }
    consumeToken();

    parseListIndex();
    curRV.listIndices.push_back(listIndex);

    ensureNotAtEOF(")");
    if (tokenInfo != TT_RightParen) {
      parseError(")");
    }
    consumeToken();


  } else {
    parseError("need discrete implementations in discrete RV");
  }

}


void
FileParser::parseContinuousImplementation()
{
  
  parseContObsDistType();

  ensureNotAtEOF("remainder of continuous random variable spec");  
  if (tokenInfo == TT_LeftParen) {

    consumeToken();

    parseListIndex();
    curRV.listIndices.push_back(listIndex);

    ensureNotAtEOF(")");
    if (tokenInfo != TT_RightParen) {
      parseError(")");
    }
    consumeToken();

    //////////////////////////////////////
    // AAA: some semantic checking here.
    // The current implementation presumably comes from a RV 
    // with nil conditional parents, as in:
    //
    // ... | nil using mixGaussian("the_forth_gaussian") | ...
    //
    // this means that we are specifying one and only one particular
    // gaussian "the_forth_gaussian" in the gaussian file,
    // since there are no conditional parents. Note that
    // this situation doesn't arrise in the discrete case here
    // because it is assumed that the pointer to the CPT will
    // be such that, if nil exists, then the CPT will not
    // have any "parents".
    // 
    // In any case, we need to make sure here that there is indeed
    // 'nil' as the current parents. Do this by checking the length
    // of the most recently pushed parrent list, and making
    // sure that it is zero.
    if (curRV.conditionalParents[curRV.conditionalParents.size()-1].size() > 0)
      parseError("decision tree 'mapping' needed when conditional parents are given");


  } else {

    // parse a string of the form 'collection("foo")'

    ensureNotAtEOF(KW_Collection);
    if (tokenInfo != KW_Collection)
      parseError(KW_Collection);
    consumeToken();

    ensureNotAtEOF("(");
    if (tokenInfo != TT_LeftParen) {
      parseError("(");
    }
    consumeToken();

    ensureNotAtEOF("collection name");
    string collectionName;
    if (tokenInfo == TT_String) {
      collectionName = tokenInfo.tokenStr;
      // remove the double quotes
      collectionName.replace(0,1,"");
      collectionName.replace(collectionName.length()-1,1,"");
      consumeToken();
    } else
      parseError("expecting collection name");

    ensureNotAtEOF(")");
    if (tokenInfo != TT_RightParen) {
      parseError(")");
    }
    consumeToken();

    // parse a string of the form 'mapping("bar")'
    parseMappingSpec();

    listIndex.collectionName = collectionName;

    // save the data.
    curRV.listIndices.push_back(listIndex);

    //////////////////////////////////////////////////////////////////
    // semantic check, this is the dual check of the
    // check at AAA above.
    if (curRV.conditionalParents[curRV.conditionalParents.size()-1].size() == 0)
      parseError("decision tree 'mapping' needs > 0 conditional parents");
  }
}


void
FileParser::parseContObsDistType()
{
  ensureNotAtEOF("continuous observation distribution type");
  if (tokenInfo == KW_MixGaussian) {
    curRV.contImplementations.push_back(MixGaussiansCommon::ci_mixGaussian);
    consumeToken();
  } else if (tokenInfo == KW_GausSwitchMixGaussian) {
    curRV.contImplementations.push_back(MixGaussiansCommon::ci_gausSwitchMixGaussian);
    consumeToken();
  }  else if (tokenInfo == KW_LogitSwitchMixGaussian) {
    curRV.contImplementations.push_back(MixGaussiansCommon::ci_logitSwitchMixGaussian);
    consumeToken();
  }  else if (tokenInfo == KW_MlpSwitchMixGaussin) {
    curRV.contImplementations.push_back(MixGaussiansCommon::ci_mlpSwitchMixGaussian);
    consumeToken();
  } else 
    parseError("continuous observation distribution type");
}

void
FileParser::parseChunkSpecifier()
{
  ensureNotAtEOF(KW_Chunk);
  if (tokenInfo != KW_Chunk) 
    parseError(KW_Chunk);
  consumeToken();

  ensureNotAtEOF("first chunk integer");
  if (tokenInfo != TT_Integer)
    parseError("first chunk integer");
  if (tokenInfo.int_val < 0) 
    parseError("non-negative chunk range");
  _firstChunkframe = (unsigned)tokenInfo.int_val;
  consumeToken();

  ensureNotAtEOF("chunk colon");
  if (tokenInfo != TT_Colon) 
    parseError("expecting chunk colon");
  consumeToken();

  ensureNotAtEOF("second chunk integer");
  if (tokenInfo != TT_Integer) 
    parseError("expecting second chunk integer");
  if (tokenInfo.int_val < 0) 
    parseError("non-negative chunk range");
  _lastChunkframe = (unsigned)tokenInfo.int_val;

  if (_firstChunkframe > _lastChunkframe)
    parseError("last chunk integer must be no less than first chunk integer");


  consumeToken();

}


void
FileParser::parseMappingSpec()
{

  ensureNotAtEOF(KW_Mapping);
  if (tokenInfo != KW_Mapping)
    parseError(KW_Mapping);
  consumeToken();

  ensureNotAtEOF("(");
  if (tokenInfo != TT_LeftParen) {
    parseError("(");
  }
  consumeToken();

  parseListIndex();

  ensureNotAtEOF(")");
  if (tokenInfo != TT_RightParen) {
    parseError(")");
  }
  consumeToken();

}

void
FileParser::parseListIndex()
{
  ensureNotAtEOF("list index");
  if (tokenInfo == TT_Integer) {
    listIndex.liType = RVInfo::ListIndex::li_Index;
    listIndex.intIndex = tokenInfo.int_val;
    consumeToken();
  } else if (tokenInfo == TT_String) {
    listIndex.liType = RVInfo::ListIndex::li_String;
    listIndex.nameIndex = tokenInfo.tokenStr;
    // remove the double quotes
    // listIndex.nameIndex.replace(listIndex.nameIndex.find("\""),1,"");
    // listIndex.nameIndex.replace(listIndex.nameIndex.find("\""),1,"");
    listIndex.nameIndex.replace(0,1,"");
    listIndex.nameIndex.replace(listIndex.nameIndex.length()-1,1,"");
    consumeToken();
  } else
    parseError("expecting list index");
}





////////////////////////////////////////////////////////////////////
//        Create Random Variables
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * createRandomVariableGraph()
 *      Takes the list of random variables that have been created
 *      by the parser and produces a valid fully linked 
 *      RV graph.
 *
 * Preconditions:
 *      successful parse from before.
 *      parseGraphicalModel() must have been called.
 *
 * Postconditions:
 *      New RV graph created.
 *
 * Side Effects:
 *      Changes the internal structure of the parser
 *      objects to create random variables.
 *
 * Results:
 *      nothing
 *
 *-----------------------------------------------------------------------
 */

void
FileParser::createRandomVariableGraph()
{

  // first create the RV objects and fill in var
  // counts as we go.
  numVarsInPrologue = numVarsInChunk = numVarsInEpilogue = 0;
  for (unsigned i=0;i<rvInfoVector.size();i++) {
    const unsigned frame = rvInfoVector[i].frame;
    if (frame < _firstChunkframe)
      numVarsInPrologue++;
    else if (frame >= _firstChunkframe &&
	     frame <= _lastChunkframe)
      numVarsInChunk++;
    else if (frame <= _maxFrame)
      numVarsInEpilogue++;
    else // shouldn't happen
      assert(0);

    if (rvInfoVector[i].rvType == RVInfo::t_discrete) {
      DiscreteRandomVariable*rv = 
	new DiscreteRandomVariable(rvInfoVector[i],
				   rvInfoVector[i].name,
				   rvInfoVector[i].rvCard);
      rv->hidden = (rvInfoVector[i].rvDisp == RVInfo::d_hidden);
      if (!rv->hidden) {
	if (rvInfoVector[i].rvFeatureRange.filled == RVInfo::FeatureRange::fr_Range) {
	  rv->featureElement = 
	    rvInfoVector[i].rvFeatureRange.firstFeatureElement;
	} else if (rvInfoVector[i].rvFeatureRange.filled == RVInfo::FeatureRange::fr_FirstIsValue) {
	  if (rvInfoVector[i].rvFeatureRange.firstFeatureElement >=
	      rvInfoVector[i].rvCard)
	    error("Error: RV \"%s\" at frame %d (line %d) specifies a value observation (%d) incompatible with its cardinality (%d)\n",
		  rvInfoVector[i].name.c_str(),
		  rvInfoVector[i].frame,
		  rvInfoVector[i].fileLineNumber,
		  rvInfoVector[i].rvFeatureRange.firstFeatureElement,
		  rvInfoVector[i].rvCard);
	  rv->val = rvInfoVector[i].rvFeatureRange.firstFeatureElement;
	  rv->featureElement = DRV_USE_FIXED_VALUE_FEATURE_ELEMENT;
	} else {
	  // this shouldn't happen (the parser should not let this case occur) 
	  // but we keep the check just in case.
	  error("ERROR: internal parser error, feature range is not filled in");
	}
      }
      rvInfoVector[i].rv = rv;
    } else {
      ContinuousRandomVariable*rv = 
	new ContinuousRandomVariable(rvInfoVector[i],
				     rvInfoVector[i].name);
      rv->hidden = (rvInfoVector[i].rvDisp == RVInfo::d_hidden);
      if (!rv->hidden) {
	rv->firstFeatureElement = 
	  rvInfoVector[i].rvFeatureRange.firstFeatureElement;
	rv->lastFeatureElement = 
	  rvInfoVector[i].rvFeatureRange.lastFeatureElement;
      }
      rvInfoVector[i].rv = rv;
    }
    rvInfoVector[i].rv->timeIndex = rvInfoVector[i].frame;

    // add weight value stuff
    rvInfoVector[i].rv->wtStatus = RandomVariable::wt_NoWeight;
    if (rvInfoVector[i].rvWeightInfo.wt_Status == RVInfo::WeightInfo::wt_Constant) {
      rvInfoVector[i].rv->wtStatus = RandomVariable::wt_Constant;
      rvInfoVector[i].rv->wtWeight = rvInfoVector[i].rvWeightInfo.weight_value;
    } else if (rvInfoVector[i].rvWeightInfo.wt_Status == RVInfo::WeightInfo::wt_Observation) {
      rvInfoVector[i].rv->wtStatus = RandomVariable::wt_Observation;
      rvInfoVector[i].rv->wtFeatureElement = rvInfoVector[i].rvWeightInfo.firstFeatureElement;
    }

  }

  // now set up all the parents of each random variable.
  for (unsigned i=0;i<rvInfoVector.size();i++) {

    const unsigned frame = rvInfoVector[i].frame;

    // build the switching parents list.
    vector<RandomVariable *> sparents;
    for (unsigned j=0;j<rvInfoVector[i].switchingParents.size();j++) {

      RVInfo::rvParent pp(rvInfoVector[i].switchingParents[j].first,
		  frame+rvInfoVector[i].switchingParents[j].second);

      ////////////////////////////////////////////////////
      // Make sure the rv at the time delta from the current
      // frame exists.
      if (nameRVmap.find(pp) == nameRVmap.end())
	error("Error: RV \"%s\" at frame %d (line %d) specifies a parent \"%s\" at frame %d that does not exist\n",
	      rvInfoVector[i].name.c_str(),
	      rvInfoVector[i].frame,
	      rvInfoVector[i].fileLineNumber,
	      pp.first.c_str(),pp.second);
      // error("Error: parent random variable \"%s\" at frame %d does not exist\n",
      // pp.first.c_str(),pp.second);

      const RVInfo& par = rvInfoVector[ nameRVmap[ pp ] ];

      ////////////////////////////////////////////////////
      // make sure we have no continuous parents.
      if (par.rvType == RVInfo::t_continuous)
	error("Error: RV \"%s\" at frame %d (line %d) specifies a continuous parent \"%s\" at frame %d (line %d)\n",
	      rvInfoVector[i].name.c_str(),
	      rvInfoVector[i].frame,
	      rvInfoVector[i].fileLineNumber,
	      par.name.c_str(),
	      par.frame,
	      par.fileLineNumber);

      // add
      sparents.push_back(par.rv);
    }
    
    // now build conditional parent list
    vector<vector<RandomVariable * > > cpl(rvInfoVector[i].conditionalParents.size());
    for (unsigned j=0;j<rvInfoVector[i].conditionalParents.size();j++) {
      for (unsigned k=0;k<rvInfoVector[i].conditionalParents[j].size();k++) {

	RVInfo::rvParent pp(rvInfoVector[i].conditionalParents[j][k].first,
		    frame+rvInfoVector[i].conditionalParents[j][k].second);

	////////////////////////////////////////////////////
	// Make sure the rv at the time delta from the current
	// frame exists.
	if (nameRVmap.find(pp) == nameRVmap.end())
	  error("Error: RV \"%s\" at frame %d (line %d) specifies a parent \"%s\" at frame %d that does not exist\n",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		pp.first.c_str(),pp.second);
	// error("Error: parent random variable \"%s\" at frame %d does not exist\n",
	// pp.first.c_str(),pp.second);

	const RVInfo& par = rvInfoVector[ nameRVmap[ pp ] ];

	////////////////////////////////////////////////////
	// make sure we have no continuous parents.
	if (par.rvType == RVInfo::t_continuous)
	  error("Error: RV \"%s\" at frame %d (line %d) specifies a continuous parent \"%s\" at frame %d (line %d)\n",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		par.name.c_str(),
		par.frame,
		par.fileLineNumber);
	// add
	cpl[j].push_back(par.rv);
      }

    }

    // finally, add all the parents.
    rvInfoVector[i].rv->setParents(sparents,cpl);
  }
}


/*-
 *-----------------------------------------------------------------------
 * ensureS_SE_E_NE,
 *    ensure links are "south", "south east",
 *    "east", or "north east", meaning that there
 *    is a numeric ordering on the nodes such that 
 *    any parents of a node at a particular
 *    numeric position have their position
 *    earlier in the ordering.
 *
 * Preconditions:
 *      parseGraphicalModel() must have been called.
 *
 * Postconditions:
 *      ordering exists if program is still running.
 *
 * Side Effects:
 *      none other than possibly killing the program.
 *
 * Results:
 *      nothing.
 *
 *-----------------------------------------------------------------------
 */

void
FileParser::ensureS_SE_E_NE()
{

  // now set up all the parents of each random variable.
  for (unsigned i=0;i<rvInfoVector.size();i++) {
    for (unsigned j=0;j<rvInfoVector[i].switchingParents.size();j++) {

      RVInfo::rvParent pp(rvInfoVector[i].switchingParents[j].first,
		  rvInfoVector[i].frame
		  +rvInfoVector[i].switchingParents[j].second);

      ////////////////////////////////////////////////////
      // Make sure the rv at the time delta from the current
      // frame exists.
      if (nameRVmap.find(pp) == nameRVmap.end())
	error("Error: RV \"%s\" at frame %d (line %d) specifies a parent \"%s\" at frame %d that does not exist\n",
	      rvInfoVector[i].name.c_str(),
	      rvInfoVector[i].frame,
	      rvInfoVector[i].fileLineNumber,
	      pp.first.c_str(),pp.second);
      // error("Error: parent random variable \"%s\" at frame %d does not exist\n",
      // pp.first.c_str(),pp.second);

      unsigned parent_position = nameRVmap[ pp ];

      if (parent_position > i) {
	const RVInfo& par = rvInfoVector[ nameRVmap[ pp ] ];
	error("Error: parent variable \"%s\", frame %d (line %d) is later than child \"%s\" frame %d (line %d)",
	      par.name.c_str(),
	      par.frame,
	      par.fileLineNumber,
	      rvInfoVector[i].name.c_str(),
	      rvInfoVector[i].frame,
	      rvInfoVector[i].fileLineNumber
	      );
      }
    }

    for (unsigned j=0;j<rvInfoVector[i].conditionalParents.size();j++) {
      for (unsigned k=0;k<rvInfoVector[i].conditionalParents[j].size();k++) {

	RVInfo::rvParent pp(rvInfoVector[i].conditionalParents[j][k].first,
		    rvInfoVector[i].frame		    
		    +rvInfoVector[i].conditionalParents[j][k].second);

	////////////////////////////////////////////////////
	// Make sure the rv at the time delta from the current
	// frame exists.
	if (nameRVmap.find(pp) == nameRVmap.end())
	  error("Error: RV \"%s\" at frame %d (line %d) specifies a parent \"%s\" at frame %d that does not exist\n",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		pp.first.c_str(),pp.second);
	// error("Error: parent random variable \"%s\" at frame %d does not exist\n",
	// pp.first.c_str(),pp.second);
	unsigned parent_position = nameRVmap[ pp ];

	if (parent_position > i) {
	  const RVInfo& par = rvInfoVector[ nameRVmap[ pp ] ];
	  error("Error: parent variable \"%s\", frame %d (line %d) is later than child \"%s\" frame %d (line %d)",
		par.name.c_str(),
		par.frame,
		par.fileLineNumber,
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber
		);
	}
      }
    }
  }

}


/*-
 *-----------------------------------------------------------------------
 * associateWithDataParams 
 *      This sets up the set of random variables 
 *      were created in createRandomVariableGraph() and associates their
 *      paramers with that which is contained in GM_Params.
 *
 * Preconditions:
 *      createRandomVariableGraph() and parseGraphicalModel() must have been called.
 *
 * Postconditions:
 *      All RVs have parameters now.
 *
 * Side Effects:
 *      changes all RV objects.
 *
 * Results:
 *      nothing.
 *
 *-----------------------------------------------------------------------
 */
void
FileParser::associateWithDataParams(MdcptAllocStatus allocate)
{
  // now set up rest of info about RV such as
  // - feature range
  // - implementations
  // - pointers to global parameter arrays & checks that they are consistent
  //   with each other.
  // Note that all the other data in GMParms must be allocated.
  for (unsigned i=0;i<rvInfoVector.size();i++) {

    // first set up the switching parent's DT mapper.
    if (rvInfoVector[i].switchingParents.size() == 0) {
      rvInfoVector[i].rv->dtMapper = NULL;
    } else {
      if (rvInfoVector[i].switchMapping.liType == 
	  RVInfo::ListIndex::li_String) {
	if (GM_Parms.dtsMap.find(rvInfoVector[i].switchMapping.nameIndex) == GM_Parms.dtsMap.end())
	  error("Error: RV \"%s\" at frame %d (line %d), switching parent DT \"%s\" doesn't exist\n",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		rvInfoVector[i].switchMapping.nameIndex.c_str());

	rvInfoVector[i].rv->dtMapper = 
	  GM_Parms.dts[
		       GM_Parms.dtsMap[rvInfoVector[i].switchMapping.nameIndex
		       ]];
      } else if  (rvInfoVector[i].switchMapping.liType == 
		  RVInfo::ListIndex::li_Index) {
	if ((rvInfoVector[i].switchMapping.intIndex < 0) ||
	    (rvInfoVector[i].switchMapping.intIndex > GM_Parms.dts.size()))
	  error("Error: RV \"%s\" at frame %d (line %d), switching parent num %d out of range\n",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		rvInfoVector[i].switchMapping.intIndex);
	rvInfoVector[i].rv->dtMapper = GM_Parms.dts[rvInfoVector[i].switchMapping.intIndex];
      } else {
	// this shouldn't happen, unless the parser has a bug.
	assert ( 0 );
      }

      // now check to make sure that the decision tree matches
      // the set of parents that were set up as the switching parents.
      if (rvInfoVector[i].rv->dtMapper->numFeatures() != 
	  rvInfoVector[i].switchingParents.size()) {
	error("Error: RV \"%s\" at frame %d (line %d), num switching parents different than required by decision tree named \"%s\".\n",
	      rvInfoVector[i].name.c_str(),
	      rvInfoVector[i].frame,
	      rvInfoVector[i].fileLineNumber,
	      rvInfoVector[i].rv->dtMapper->name().c_str());
      }
    }

    //////////////////////////////////////////////////////////
    // NOW set up the conditional parent's parameters.
    //
    if (rvInfoVector[i].rvType == RVInfo::t_discrete) {
      ///////////////////////////////////////////////////////
      //
      // DISCRETE form of the current rv
      //
      ///////////////////////////////////////////////////////
      DiscreteRandomVariable* rv = 
	(DiscreteRandomVariable*) rvInfoVector[i].rv;

      vector<CPT*> cpts(rvInfoVector[i].conditionalParents.size());
      for (unsigned j=0;j<rvInfoVector[i].conditionalParents.size();j++) {

	////////////
	// check each implementation of CPT, and add the
	// appropriate one of the appropriate type to the the 
	// virtual CPTS array "cpts". The following code
	// is essentially the same code but duplicated 
	// one time for each CPT implementation.

	if (rvInfoVector[i].discImplementations[j] == CPT::di_MDCPT) {

	  //////////////////////////////////////////////////////
	  // set the CPT to a MDCPT, depending on if a string
	  // or integer index was used in the file.
	  if (rvInfoVector[i].listIndices[j].liType 
	      == RVInfo::ListIndex::li_String) {
	    if (GM_Parms.mdCptsMap.find(
		      rvInfoVector[i].listIndices[j].nameIndex) ==
		GM_Parms.mdCptsMap.end()) {
	      if (allocate == noAllocate) {
		error("Error: RV \"%s\" at frame %d (line %d), conditional parent MDCPT \"%s\" doesn't exist\n",
		      rvInfoVector[i].name.c_str(),
		      rvInfoVector[i].frame,
		      rvInfoVector[i].fileLineNumber,
		      rvInfoVector[i].listIndices[j].nameIndex.c_str());
	      }
	      else {
		// allocate the MDCPT with name and install it.
		MDCPT* mdcpt = new MDCPT();
		mdcpt->setName(rvInfoVector[i].listIndices[j].nameIndex);
		mdcpt->
		  setNumParents
		      (rvInfoVector[i].conditionalParents[j].size());

		for (unsigned k=0;k<rvInfoVector[i].conditionalParents[j].size();k++) {


		  RVInfo::rvParent pp(rvInfoVector[i].conditionalParents[j][k].first,
			      rvInfoVector[i].frame
			      +rvInfoVector[i].conditionalParents[j][k].second);
		  map < RVInfo::rvParent , unsigned >::iterator it;
		  it = nameRVmap.find(pp);
		  if (it == nameRVmap.end()) {
		    // this really shouldn't happen at this point since
		    // it should have been checked somewhere else,
		    // but we include the check nonetheless
		    error("Error: RV \"%s\" at frame %d (line %d) specifies a parent \"%s\" at frame %d that does not exist\n",
			  rvInfoVector[i].name.c_str(),
			  rvInfoVector[i].frame,
			  rvInfoVector[i].fileLineNumber,
			  pp.first.c_str(),pp.second,
			  rvInfoVector[i].conditionalParents[j][k].first.c_str(),
			  rvInfoVector[i].conditionalParents[j][k].second);
		    // error("Error: parent random variable \"%s\" at" 
		    // "frame %d does not exist\n",
		    // rvInfoVector[i].conditionalParents[j][k].first.c_str(),
		    // rvInfoVector[i].conditionalParents[j][k].second);
		  }
		  mdcpt->setNumCardinality(k,
					   rvInfoVector[(*it).second].rvCard);
		}
		mdcpt->setNumCardinality(rvInfoVector[i].conditionalParents[j].size(),
					 rvInfoVector[i].rvCard);
		mdcpt->allocateBasicInternalStructures();
		if (allocate == allocateRandom)
		  mdcpt->makeRandom();
		else if (allocate == allocateUniform)
		  mdcpt->makeUniform();
		else
		  assert(0);

		GM_Parms.mdCpts.push_back(mdcpt);
		GM_Parms.mdCptsMap[rvInfoVector[i].listIndices[j].nameIndex]
		  = GM_Parms.mdCpts.size()-1;
		cpts[j] = mdcpt;
	      }
	    } else {
	      // otherwise add it
	      cpts[j] = 
		GM_Parms.mdCpts[
		 GM_Parms.mdCptsMap[
			  rvInfoVector[i].listIndices[j].nameIndex
		 ]
		];

	    }
	  } else {
	    // need to remove the integer index code.
	    assert (0);
#if 0
	    if (rvInfoVector[i].listIndices[j].intIndex >= 
		GM_Parms.mdCpts.size()) {
	      if (!allocateIfNotThere) {
		error("Error: RV \"%s\" at frame %d (line %d), conditional parent index (%d) too large\n",
		      rvInfoVector[i].name.c_str(),
		      rvInfoVector[i].frame,
		      rvInfoVector[i].fileLineNumber,
		      rvInfoVector[i].listIndices[j].intIndex);
	      } else {
		error("Can't allocate with integer cpt index");
	      }
	    } else {
	      // otherwise add it
	      cpts[j] =
		GM_Parms.mdCpts[rvInfoVector[i].listIndices[j].intIndex];
	    }
#endif
	  }

	} else 
	  if (rvInfoVector[i].discImplementations[j] == CPT::di_MSCPT) {

	    /////////////////////////////////////////////////////////
	    // same code as above, but using MSCPTs rather
	    // then MDCPTs. There should be a better way to do this.

	    //////////////////////////////////////////////////////
	    // set the CPT to a MSCPT, depending on if a string
	    // or integer index was used in the file.
	    if (rvInfoVector[i].listIndices[j].liType 
		== RVInfo::ListIndex::li_String) {
	      if (GM_Parms.msCptsMap.find(
					  rvInfoVector[i].listIndices[j].nameIndex) ==
		  GM_Parms.msCptsMap.end()) {
		error("Error: RV \"%s\" at frame %d (line %d), conditional parent MSCPT \"%s\" doesn't exist\n",
		      rvInfoVector[i].name.c_str(),
		      rvInfoVector[i].frame,
		      rvInfoVector[i].fileLineNumber,
		      rvInfoVector[i].listIndices[j].nameIndex.c_str());
	      } else {
		// otherwise add it
		cpts[j] = 
		  GM_Parms.msCpts[
				  GM_Parms.msCptsMap[
						     rvInfoVector[i].listIndices[j].nameIndex
				  ]
		  ];
	      }
	    } else {
	      // need to remove the integer index code.
	      assert(0);
#if 0
	      if (rvInfoVector[i].listIndices[j].intIndex >= 
		  GM_Parms.msCpts.size()) {
		if (!allocateIfNotThere) {
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent index (%d) too large\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].intIndex);
		} else {
		  error("Can't allocate with integer cpt index");
		}
	      } else {
		// otherwise add it
		cpts[j] =
		  GM_Parms.msCpts[rvInfoVector[i].listIndices[j].intIndex];
	      }
#endif
	    }

	} else 
	  if (rvInfoVector[i].discImplementations[j] == CPT::di_MTCPT) {

	    /////////////////////////////////////////////////////////
	    // Once again, same code as above, but using MTCPTs rather
	    // then MDCPTs or MSCPTs. 

	    //////////////////////////////////////////////////////
	    // set the CPT to a MTCPT, depending on if a string
	    // or integer index was used in the file.
	    if (rvInfoVector[i].listIndices[j].liType 
		== RVInfo::ListIndex::li_String) {
	      if (GM_Parms.mtCptsMap.find(
					  rvInfoVector[i].listIndices[j].nameIndex) ==
		  GM_Parms.mtCptsMap.end()) {
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent MTCPT \"%s\" doesn't exist\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].nameIndex.c_str());
	      } else {
		// otherwise add it
		cpts[j] = 
		  GM_Parms.mtCpts[
				  GM_Parms.mtCptsMap[
						     rvInfoVector[i].listIndices[j].nameIndex
				  ]
		  ];
	      }
	    } else {
	      // need to remove the integer index code.
	      assert(0);
#if 0
	      if (rvInfoVector[i].listIndices[j].intIndex >= 
		  GM_Parms.mtCpts.size()) {
		if (!allocateIfNotThere) {
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent index (%d) too large\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].intIndex);
		} else {
		  error("Can't allocate with integer cpt index");
		}
	      } else {
		// otherwise add it
		cpts[j] =
		  GM_Parms.mtCpts[rvInfoVector[i].listIndices[j].intIndex];
	      }
#endif
	    }

	} else {
	  // Again, this shouldn't happen. If it does, something is wrong
	  // with the parser code or some earlier code, and it didn't correctly
	  // set the CPT type. 
	  assert ( 0 );
	}


	// do checking to ensure that CPTs specified are
	// compatible with parent cardinalities.

	// first get the name for better error reporting.
	string cptType;
	if (rvInfoVector[i].discImplementations[j] == CPT::di_MDCPT) {
	  cptType = "MDCPT";
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_MSCPT) {
	  cptType = "MSCPT";
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_MTCPT) {
	  cptType = "MTCPT";
	}

	// check to make sure this cpt matches this
	// number of parents.
	if (cpts[j]->numParents() != 
	    rvInfoVector[i].conditionalParents[j].size()) {
	  error("Error: RV \"%s\" at frame %d (line %d), num parents cond. %d different than required by %s \"%s\".\n",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		j,
		cptType.c_str(),
		cpts[j]->name().c_str());
	}

	if (cpts[j]->cptType == CPT::di_USCPT) {
	  // then we need to do special checking for an USCPT, meaning
	  // we need to check only if that this discrete random
	  // variable is observed. Note we've already checked for
	  // the number of parents above, so since a USCPT always
	  // has zero parents, we don't need to do further
	  // checking of parent cardinality below.
	  if (rvInfoVector[i].rvDisp != RVInfo::d_observed) {
	    error("Error: RV '%s' at frame %d (line %d), wants to use a unity score DenseCPT which requires an observed variable with no parents, but variable '%s' is not specified as observed.\n",
		  rvInfoVector[i].name.c_str(),
		  rvInfoVector[i].frame,
		  rvInfoVector[i].fileLineNumber,
		  rvInfoVector[i].name.c_str());
	  }
	} else {
	  // regular non USCPT checking below.
	  if ((unsigned)cpts[j]->card() != rvInfoVector[i].rvCard) {
	    error("Error: RV \"%s\" at frame %d (line %d), cardinality of RV is %d, but %s \"%s\" requires cardinality of %d.\n",
		  rvInfoVector[i].name.c_str(),
		  rvInfoVector[i].frame,
		  rvInfoVector[i].fileLineNumber,
		  rvInfoVector[i].rvCard,
		  cptType.c_str(),		
		  cpts[j]->name().c_str(),
		  cpts[j]->card());
	  }
	  for (unsigned par=0;par<cpts[j]->numParents();par++) {
	    if (rv->conditionalParentsList[j][par]->cardinality !=
		cpts[j]->parentCardinality(par))
	      error("Error: RV \"%s\" at frame %d (line %d), cardinality of parent '%s' is %d, but %d'th parent of %s \"%s\" requires cardinality of %d.\n",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rv->conditionalParentsList[j][par]->name().c_str(),
		    rv->conditionalParentsList[j][par]->cardinality,
		    par,
		    cptType.c_str(),
		    cpts[j]->name().c_str(),
		    cpts[j]->parentCardinality(par));
	  }
	}
      }

      // finally, give the the cpts to the rv
      rv->setCpts(cpts);

    } else { 
      ///////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////
      // This is a CONTINUOUS RV, so 
      // get a cont. form of the current rv
      ///////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////
      ContinuousRandomVariable* rv = 
	(ContinuousRandomVariable*) rvInfoVector[i].rv;

      rv->conditionalGaussians.resize
	(
	 rvInfoVector[i].conditionalParents.size()
	 );
      for (unsigned j=0;j<rvInfoVector[i].conditionalParents.size();j++) {

	// TODO:: implement the other continuous implementations
	// and allow the other cases here.
	if (rvInfoVector[i].contImplementations[j] !=
	    MixGaussiansCommon::ci_mixGaussian) {
	  error("ERROR: Only supports mixGaussian implementations for now");
	  // ultimately, we can remove this and we'll have to duplicate
	  // the below code for each continuous implementation, just
	  // like what was done above for the discrete implementations.
	}

	if (rvInfoVector[i].conditionalParents[j].size() == 0) {
	  // then under this value for switching parents, there
	  // are no conditional parents and we have a direct
	  // pointer to some Gaussian mixture.
	  rv->conditionalGaussians[j].direct = true;
	  if (rvInfoVector[i].listIndices[j].liType 
		== RVInfo::ListIndex::li_String) {
	    if (GM_Parms.mixGaussiansMap.find(
		    rvInfoVector[i].listIndices[j].nameIndex) ==
		GM_Parms.mixGaussiansMap.end()) {
	      error("Error: RV \"%s\" at frame %d (line %d), Gaussian mixture \"%s\" doesn't exist\n",
		      rvInfoVector[i].name.c_str(),
		      rvInfoVector[i].frame,
		      rvInfoVector[i].fileLineNumber,
		      rvInfoVector[i].listIndices[j].nameIndex.c_str());
	    } else {
	      // Mix Gaussian is there, so add it.
	      rv->conditionalGaussians[j].gaussian =
		GM_Parms.mixGaussians[
		      GM_Parms.mixGaussiansMap[
			  rvInfoVector[i].listIndices[j].nameIndex
		      ]];
	    }
	  } else {
	    // need to remove the integer index code.
	    assert(0);
#if 0
	    // the list index is an integer
	    if (rvInfoVector[i].listIndices[j].intIndex >= 
		GM_Parms.mixGaussiansMap.size()) {
		if (!allocateIfNotThere) {
		  error("Error: RV \"%s\" at frame %d (line %d), Gaussian mixture index (%d) too large\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].intIndex);
		} else {
		  error("Can't allocate with integer Gaussian mixture index");
		}
	      } else {
		// otherwise add it
		rv->conditionalGaussians[j].gaussian =
		  GM_Parms.mixGaussians[
                        rvInfoVector[i].listIndices[j].intIndex
		      ];
	      }
#endif
	  }
	} else {
	  // there are > 0 conditional parents for this
	  // set of switching values. The index should
	  // specify a DT which maps from the set of
	  // conditional parents to a GM collection.
	  rv->conditionalGaussians[j].direct = false;
	  if (rvInfoVector[i].listIndices[j].liType 
		== RVInfo::ListIndex::li_String) {
	    // string name index to DT
	    if (GM_Parms.dtsMap.find(
		    rvInfoVector[i].listIndices[j].nameIndex) ==
		GM_Parms.dtsMap.end()) {
		error("Error: RV \"%s\" at frame %d (line %d), conditional parent %d specifies a DT \"%s\" that doesn't exist\n",
		      rvInfoVector[i].name.c_str(),
		      rvInfoVector[i].frame,
		      rvInfoVector[i].fileLineNumber,
		      j,
		      rvInfoVector[i].listIndices[j].nameIndex.c_str());
	    } else {
	      rv->conditionalGaussians[j].mapping.dtMapper =
		GM_Parms.dts[
		   GM_Parms.dtsMap[
				rvInfoVector[i].listIndices[j].nameIndex
		  ]];
	    }
	  } else {
	    // int index to DT
	    if (rvInfoVector[i].listIndices[j].intIndex >= 
		GM_Parms.dts.size()) {
	      error("Error: RV \"%s\" at frame %d (line %d), conditional parent %d specifies a DT index (%d) that is too large\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			j,
			rvInfoVector[i].listIndices[j].intIndex);
	      } else {
		// otherwise add it
		rv->conditionalGaussians[j].mapping.dtMapper =
		  GM_Parms.dts[rvInfoVector[i].listIndices[j].intIndex];
	      }
	  }

	  // get the collection as well.
	  if (GM_Parms.nclsMap.find(
		    rvInfoVector[i].listIndices[j].collectionName) ==
		GM_Parms.nclsMap.end()) {
		error("Error: RV \"%s\" at frame %d (line %d), conditional parent %d specifies a collection name \"%s\" that doesn't exist\n",
		      rvInfoVector[i].name.c_str(),
		      rvInfoVector[i].frame,
		      rvInfoVector[i].fileLineNumber,
		      j,
		      rvInfoVector[i].listIndices[j].collectionName.c_str());
	    } else {
	      rv->conditionalGaussians[j].mapping.collection =
		GM_Parms.ncls[
		   GM_Parms.nclsMap[
				rvInfoVector[i].listIndices[j].collectionName
		  ]];
	      // make sure that the pointer table is filled in for this name.
	      rv->conditionalGaussians[j].mapping.collection->fillMgTable();
	      // Might as well do this here. Check that all Gaussians
	      // in the collection are of the right dimensionality.
	      // This check is important for ContinuousRandomVariable::probGivenParents().
	      const unsigned rvDim = 
		(rvInfoVector[i].rvFeatureRange.lastFeatureElement - 
		 rvInfoVector[i].rvFeatureRange.firstFeatureElement)+1;
	      // make sure all Gaussians in the collection have the
	      // same dimensionality as the RV.
	      for (unsigned u=0;u<
		     rv->conditionalGaussians[j].mapping.collection->mgSize();
		   u++) {
		if ( // we only check the dimension for those Gaussians that have it.
		    (rv->conditionalGaussians[j].mapping.collection->mg(u)->mixType
		     != MixGaussiansCommon::ci_zeroScoreMixGaussian)
		    &&
		    (rv->conditionalGaussians[j].mapping.collection->mg(u)->mixType
		     != MixGaussiansCommon::ci_unityScoreMixGaussian)
		    && 
		    (rv->conditionalGaussians[j].mapping.collection->mg(u)->dim() 
		     != rvDim)) {
		  error("Error: RV \"%s\" at frame %d (line %d), dimensionality %d, conditional parent %d, collection \"%s\" specifies Gaussian \"%s\" at pos %d with wrong dimension %d\n",
		      rvInfoVector[i].name.c_str(),
		      rvInfoVector[i].frame,
		      rvInfoVector[i].fileLineNumber,
		      rvDim,
		      j,
		      rvInfoVector[i].listIndices[j].collectionName.c_str(),
		      rv->conditionalGaussians[j].mapping.collection->mg(u)->name().c_str(),
		      u,
		      rv->conditionalGaussians[j].mapping.collection->mg(u)->dim());
		}
	      }
	    }
	}
      }
    }
  }
}


/*-
 *-----------------------------------------------------------------------
 * checkConsistentWithGlobalObservationStream
 *      Make sure that the obs variables defined here are 
 *      ok with the current global observation stream object.
 *
 * Preconditions:
 *      Global observation stream should have been set up (i.e.,
 *      number features per frame, etc. all set up).
 *      createRandomVariableGraph() MUST HAVE BEEN CALLED FIRST.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      might kill the program.
 *
 * Results:
 *      nothing.
 *
 *-----------------------------------------------------------------------
 */
void
FileParser::checkConsistentWithGlobalObservationStream()
{
  for (unsigned i=0;i<rvInfoVector.size();i++) { 
    if (rvInfoVector[i].rvType == RVInfo::t_discrete) {
      DiscreteRandomVariable* rv = 
	(DiscreteRandomVariable*) rvInfoVector[i].rv;
      if (!rv->hidden) {
	if (rv->featureElement != DRV_USE_FIXED_VALUE_FEATURE_ELEMENT) {
	  if (!globalObservationMatrix.elementIsDiscrete(rv->featureElement)) {
	    if (globalObservationMatrix.numDiscrete() > 0) 
	      error("ERROR: discrete observed random variable '%s', frame %d, line %d, specifies a feature element %d:%d that is out of discrete range ([%d:%d] inclusive) of observation matrix",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rv->featureElement,
		    rv->featureElement,
		    globalObservationMatrix.numContinuous(),
		    globalObservationMatrix.numFeatures()-1);
	    else
	      error("ERROR: discrete observed random variable '%s', frame %d, line %d, specifies a feature element %d:%d for an observation matrix with zero discrete features.",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rv->featureElement,
		    rv->featureElement,
		    globalObservationMatrix.numContinuous(),
		    globalObservationMatrix.numFeatures()-1);
	  }
	}
      }
    } else { // (rvInfoVector[i].rvType == RVInfo::t_continuous) {
      ContinuousRandomVariable* rv = 
	(ContinuousRandomVariable*) rvInfoVector[i].rv;
      if (!rv->hidden) {
	if (rv->lastFeatureElement >=  globalObservationMatrix.numContinuous())
	      error("ERROR: continuous observed random variable '%s', frame %d, line %d, specifies feature elements %d:%d that are out of continuous range ([%d:%d] inclusive) of observation matrix",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rv->firstFeatureElement,
		    rv->lastFeatureElement,
		    0,
		    globalObservationMatrix.numContinuous()-1);
      }
    }

    // checks common for both discrete and continuous RVs.

    // if weight comes from observation matrix, make sure it indexes into
    // a valid index and a float value.
    if (rvInfoVector[i].rvWeightInfo.wt_Status == RVInfo::WeightInfo::wt_Observation) {
      if (rvInfoVector[i].rvWeightInfo.lastFeatureElement >= globalObservationMatrix.numContinuous()) {
	error("ERROR: random variable '%s', frame %d, line %d, specifies weight's observation feature element %d:%d that are out of continuous range ([%d:%d] inclusive) of observation matrix",
	      rvInfoVector[i].name.c_str(),
	      rvInfoVector[i].frame,
	      rvInfoVector[i].fileLineNumber,
	      rvInfoVector[i].rvWeightInfo.firstFeatureElement,
	      rvInfoVector[i].rvWeightInfo.lastFeatureElement,
	      0,
	      globalObservationMatrix.numContinuous()-1);
      }
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * addVariablesToGM()
 *      Adds all the variables that have been created by the
 *      parser to the GM.
 *
 * Preconditions:
 *      createRandomVariableGraph() and parseGraphicalModel(),
 *      and associateWithDataParams() must have been called.
 *
 * Postconditions:
 *      GM object contains all of the RVs.
 *
 * Side Effects:
 *      none internally, but changes gm object.
 *
 * Results:
 *      nothing.
 *
 *-----------------------------------------------------------------------
 */
void
FileParser::addVariablesToGM(GMTK_GM& gm)
{
  gm.node.resize(rvInfoVector.size());
  for (unsigned i=0;i<rvInfoVector.size();i++) {
    gm.node[i] = rvInfoVector[i].rv;
  }
  // set some of the GMs variables.
  gm.firstChunkFrame = firstChunkFrame();
  gm.lastChunkFrame = lastChunkFrame();
  gm.framesInTemplate = numFrames();
  gm.framesInRepeatSeg = lastChunkFrame()- firstChunkFrame() + 1;
}



#if 0
/*-
 *-----------------------------------------------------------------------
 * addVariablesToGMTemplate()
 *      Adds all the variables that have been created by the
 *      parser to the GMTemplate argument. This routine is like
 *      an "export" command, in that it exports only the
 *      graph information contained in the file w/o any of the
 *      file specific information. After calling this routine,
 *      the template is entirely dissociated with the fileparser,
 *      and the fileparser information may be deleted.
 *
 * Preconditions:
 *      createRandomVariableGraph() and parseGraphicalModel(),
 *      and associateWithDataParams() must have been called.
 *
 * Postconditions:
 *      GMTemplate contains a stripped down version
 *      of the template contained in the fileparser, but without
 *      all of the file/line number information.
 *
 * Side Effects:
 *      none internally, but changes GMTemplate object.
 *
 * Results:
 *      nothing.
 *
 *-----------------------------------------------------------------------
 */
void
FileParser::addVariablesToTemplate(GMTemplate& gm_template)
{
  gm_template.numFrames = numFrames();

  // TODO: rhs, computed already?
  gm_template.prologueNumFrames = firstChunkFrame();
  gm_template.chunkNumFrames = lastChunkFrame() - firstChunkFrame() + 1;
  gm_template.epilogueNumFrames = numFrames() - firstChunkFrame() - 1;

  gm_template.firstChunkFrame = firstChunkFrame();
  gm_template.lastChunkFrame = lastChunkFrame();

  // unroll the file parser's template 1 time and place
  // it in the GM template.
  unroll(2,gm_template.rvs);

//   gm_template.frames.resize(numFrames());
//   ///////////////////////////////////////////////////////////
//   // Assume that the rvs inserted into rvInfoVector are 
//   // insereted in frame order.
//   unsigned rvindex=0;
//   for (unsigned frameIndex=0;frameIndex<numFrames();frameIndex++) {
//     while (rvInfoVector[rvindex].frame == frameIndex) {
//       gm_template.frames[frameIndex].rvs.push_back(rvInfoVector[rvindex].rv);
//       rvindex++;
//     }
//  }

}

#endif



/*-
 *-----------------------------------------------------------------------
 * unroll()
 *      Unrolls the internal stored template some number of times
 *      and returns the result in the vector of random variables.
 *      Note, this is a general unrolling routine, it allows
 *      for variables to have parents both in the past and in
 *      the future. It also allows for networks that have
 *      cycles. It does not allow for 
 *
 *      In general, unrolling is valid only if all parent variables in
 *      the unrolled network are {\em compatible} with those parents
 *      that exist in the template. A variable is said to be
 *      compatible if has 1) the same name, 2) the same type (i.e.,
 *      discrete or continuous), and 3) the same cardinality (number
 *      of possible values).  For example, suppose that $A$ is a
 *      random variable in the template having $B$ as a parent that is
 *      $k$ frames to the left of $A$ in the template. After
 *      unrolling, each variable $A'$ in the unrolled network that is
 *      derived from $A$ in the template must have a parent having the
 *      same name as $B$, it must be $k$ frames to the left of each
 *      $A'$, and must have the same type and (if discrete)
 *      cardinality as $B$. If these conditions are not met, the
 *      template is invalid.
 *
 *      Note: unrolling 0 times means return a structure like the
 *      template. Unrolling 1 times means duplicate the chunk once.
 *      
 *      Note: This routine is part of FileParser becuase it uses
 *      the parent offset information that is contained in 
 *      the template.
 *
 *      There are TWO interfaces to this routine. One of them
 *      also returns a map from rvParent structures to the position
 *      in the returned array of the corresponding variable (which
 *      might be useful to the caller).
 *
 * Preconditions:
 *      createRandomVariableGraph() and parseGraphicalModel(),
 *      and associateWithDataParams() must have been called.
 *      Assumes that rvInfoVector exists and template r.v.s
 *      have been stored in frame order (frame 0 first, 1 second, etc.)
 *
 * Postconditions:
 *      argument 'unrolledVarSet' contains unrolled network.
 *      Order of vars in 'unrolledVarSet' is GUARANTEED to be
 *      in same order as before (i.e., in frame order
 *      which will be the order the variables are presented in 
 *      the file).
 *
 * Side Effects:
 *      Nothing internal is changed, but changes argument.
 *
 * Results:
 *      nothing.
 *
 *----------------------------------------------------------------------- 
 */
void
FileParser::unroll(unsigned timesToUnroll,
		   vector<RandomVariable*> &unrolledVarSet)
{
  // a map from the r.v. name and new frame number to
  // position in the new unrolled array of r.v.'s
  map < RVInfo::rvParent, unsigned > posOfParentAtFrame;
  return unroll(timesToUnroll,unrolledVarSet,posOfParentAtFrame);
}
void
FileParser::unroll(unsigned timesToUnroll,
		   vector<RandomVariable*> &unrolledVarSet,
		   map < RVInfo::rvParent, unsigned >& posOfParentAtFrame)
{

  // a map from the new r.v. to the corresponding rv info
  // vector position in the template.
  map < RandomVariable *, unsigned > infoOf;

  // first go though and create all of the new r.v., without parents,
  // but with parameters shared

  // clear out the old var set if any.
  for (unsigned i=0;i<unrolledVarSet.size();i++)
    delete unrolledVarSet[i];
  unrolledVarSet.clear();
  // set up the size for the new one in advance.
  unrolledVarSet.resize(numVarsInPrologue+
			(timesToUnroll+1)*numVarsInChunk+numVarsInEpilogue);
  // uvsi: unrolledVarSet index
  unsigned uvsi=0;
  // tvi: template variable indexa
  unsigned tvi=0;
  for (;tvi<numVarsInPrologue;tvi++) {

    const unsigned frame = rvInfoVector[tvi].frame;

    unrolledVarSet[uvsi] = rvInfoVector[tvi].rv->cloneWithoutParents();

    // not necessary for prologue variables, since
    // time index retained in clone.
    // unrolledVarSet[uvsi]->timeIndex = frame;

    infoOf[unrolledVarSet[uvsi]] = tvi;
    RVInfo::rvParent p(unrolledVarSet[uvsi]->name(),frame);
    posOfParentAtFrame[p] = uvsi;

    uvsi++;
  }

  const unsigned numFramesInChunk = (_lastChunkframe-_firstChunkframe+1);
  for (unsigned i=0;i<(timesToUnroll+1);i++) {
    tvi = numVarsInPrologue;
    for (unsigned j=0;j<numVarsInChunk;j++) {
      const unsigned frame = 
	i*numFramesInChunk+rvInfoVector[tvi].frame;

      unrolledVarSet[uvsi] = rvInfoVector[tvi].rv->cloneWithoutParents();
      unrolledVarSet[uvsi]->timeIndex = frame;

      infoOf[unrolledVarSet[uvsi]] = tvi;
      RVInfo::rvParent p(unrolledVarSet[uvsi]->name(),frame);
      posOfParentAtFrame[p] = uvsi;

      uvsi++;
      tvi++;
    }
  }

  for (unsigned i=0;i<numVarsInEpilogue;i++) {
    const unsigned frame = 
      timesToUnroll*numFramesInChunk+rvInfoVector[tvi].frame;

    unrolledVarSet[uvsi] = rvInfoVector[tvi].rv->cloneWithoutParents();
    unrolledVarSet[uvsi]->timeIndex = frame;

    infoOf[unrolledVarSet[uvsi]] = tvi;
    RVInfo::rvParent p(unrolledVarSet[uvsi]->name(),frame);
    posOfParentAtFrame[p] = uvsi;

    uvsi++;
    tvi++;
  }
  
  // now we have all the rv,s but with the wrong parents, fix that.

  for (uvsi=0;uvsi<unrolledVarSet.size();uvsi++) {

    const unsigned frame = unrolledVarSet[uvsi]->timeIndex;

    // get the info objectd for this one.
    RVInfo* info = &rvInfoVector[infoOf[unrolledVarSet[uvsi]]];
    
    // build the switching parents list.
    vector<RandomVariable *> sparents;
    for (unsigned j=0;j<info->switchingParents.size();j++) {

      // grab a pointer to the parent in the template
      RandomVariable* template_parent_rv =
	info->rv->switchingParents[j];

      // pointer to parent in unrolled network
      RVInfo::rvParent pp(info->switchingParents[j].first,
		  frame+info->switchingParents[j].second);


      ///////////////////////////////////////////
      // next set of checks ensure compatibility.
      // 1. first the name check at that frame.

      ////////////////////////////////////////////////////////
      // Make sure the rv at the time delta from the current
      // frame exists.
      map < RVInfo::rvParent , unsigned >::iterator it;      
      if ((it = posOfParentAtFrame.find(pp)) == posOfParentAtFrame.end()) {
	error("Error: random variable '%s' (template frame %d, line %d) specifies a parent '%s(%d)' that does not exist when unrolling template %d times (in unrolled network, variable at frame %d asked for parent '%s' at frame '%d' which does not exist)\n",
	      info->name.c_str(),
	      info->frame,
	      info->fileLineNumber,
	      info->switchingParents[j].first.c_str(),
	      info->switchingParents[j].second,
	      timesToUnroll,
	      frame,
	      info->switchingParents[j].first.c_str(),
	      frame+info->switchingParents[j].second);
      }

      ///////////////////////////////////////////
      // 2. next the type compatibility check
      RandomVariable *unrolled_parent_rv = unrolledVarSet[(*it).second];
      if (!unrolled_parent_rv->discrete) {
	error("Error: random variable '%s' (template frame %d, line %d) specifies a parent '%s(%d)' that is continous when unrolling template %d times (in unrolled network, variable at frame %d asked for parent '%s' at frame '%d' which is continuous)\n",
	      info->name.c_str(),
	      info->frame,
	      info->fileLineNumber,
	      info->switchingParents[j].first.c_str(),
	      info->switchingParents[j].second,
	      timesToUnroll,
	      frame,
	      info->switchingParents[j].first.c_str(),
	      frame+info->switchingParents[j].second);
      }

      //////////////////////////////////////////////////////////////
      // 3. lastly, check that the cardinality of the parent matches 
      // since it is a discrete r.v.
      DiscreteRandomVariable* dunrolled_parent_rv = 
	(DiscreteRandomVariable*) unrolled_parent_rv;
      
      if (dunrolled_parent_rv->cardinality != 
	  template_parent_rv->cardinality) {
	error("Error: random variable '%s' (template frame %d, line %d) specifies a parent '%s(%d)' that has the wrong cardinality when unrolling template %d times (in unrolled network, variable at frame %d asked for parent '%s' at frame '%d' which is has cardinality %d, but cardinality in template parent is %d)\n",
	      info->name.c_str(),
	      info->frame,
	      info->fileLineNumber,
	      info->switchingParents[j].first.c_str(),
	      info->switchingParents[j].second,
	      timesToUnroll,
	      frame,
	      info->switchingParents[j].first.c_str(),
	      frame+info->switchingParents[j].second,
	      unrolled_parent_rv->cardinality,
	      template_parent_rv->cardinality);
      }

      // add
      sparents.push_back(unrolled_parent_rv);

    }

    // now build the conditional parents list.
    vector<vector<RandomVariable * > > cpl(info->conditionalParents.size());
    for (unsigned j=0;j<info->conditionalParents.size();j++) {
      for (unsigned k=0;k<info->conditionalParents[j].size();k++) {


	// grab a pointer to the parent in the template
	RandomVariable* template_parent_rv =
	  info->rv->conditionalParentsList[j][k];

	RVInfo::rvParent pp(info->conditionalParents[j][k].first,
		    frame+info->conditionalParents[j][k].second);


	///////////////////////////////////////////
	// next set of checks ensure compatibility.
	// 1. first the name check at that frame.

	////////////////////////////////////////////////////////
	// Make sure the rv at the time delta from the current
	// frame exists.
	map < RVInfo::rvParent , unsigned >::iterator it;      
	if ((it = posOfParentAtFrame.find(pp)) == posOfParentAtFrame.end()) {
	  error("Error: random variable '%s' (template frame %d, line %d) specifies a parent '%s(%d)' that does not exist when unrolling template %d times (in unrolled network, variable at frame %d asked for parent '%s' at frame '%d' which does not exist)\n",
		info->name.c_str(),
		info->frame,
		info->fileLineNumber,
		info->conditionalParents[j][k].first.c_str(),
		info->conditionalParents[j][k].second,
		timesToUnroll,
		frame,
		info->conditionalParents[j][k].first.c_str(),
		frame+info->conditionalParents[j][k].second);
	}

	///////////////////////////////////////////
	// 2. next the type compatibility check
	RandomVariable *unrolled_parent_rv = unrolledVarSet[(*it).second];
	if (!unrolled_parent_rv->discrete) {
	  error("Error: random variable '%s' (template frame %d, line %d) specifies a parent '%s(%d)' that is continous when unrolling template %d times (in unrolled network, variable at frame %d asked for parent '%s' at frame '%d' which is continuous)\n",
		info->name.c_str(),
		info->frame,
		info->fileLineNumber,
		info->conditionalParents[j][k].first.c_str(),
		info->conditionalParents[j][k].second,
		timesToUnroll,
		frame,
		info->conditionalParents[j][k].first.c_str(),
		frame+info->conditionalParents[j][k].second);
	}

	//////////////////////////////////////////////////////////////
	// 3. lastly, check that the cardinality of the parent matches 
	// since it is a discrete r.v.
	DiscreteRandomVariable* dunrolled_parent_rv = 
	  (DiscreteRandomVariable*) unrolled_parent_rv;
      
	if (dunrolled_parent_rv->cardinality != 
	    template_parent_rv->cardinality) {
	  error("Error: random variable '%s' (template frame %d, line %d) specifies a parent '%s(%d)' that has the wrong cardinality when unrolling template %d times (in unrolled network, variable at frame %d asked for parent '%s' at frame '%d' which is has cardinality %d, but cardinality in template parent is %d)\n",
		info->name.c_str(),
		info->frame,
		info->fileLineNumber,
		info->conditionalParents[j][k].first.c_str(),
		info->conditionalParents[j][k].second,
		timesToUnroll,
		frame,
		info->conditionalParents[j][k].first.c_str(),
		frame+info->conditionalParents[j][k].second,
		unrolled_parent_rv->cardinality,
		template_parent_rv->cardinality);
	}

	// add
	cpl[j].push_back(unrolled_parent_rv);
    
      }
    }
    unrolledVarSet[uvsi]->setParents(sparents,cpl);
  }

  // can't call topological sort for now since clique
  // sizes get too big for frontier algorithm.

  // vector<RandomVariable*> res;
  // GraphicalModel::topologicalSort(unrolledVarSet,res);
  // unrolledVarSet = res;

}






/*-
 *-----------------------------------------------------------------------
 * writeGMId()
 *  A routine to write out the graph template (P,C,E) in condensed
 *  form but in sufficient detail so that it can be used to quickly
 *  ID the current template (e.g., so that a given elimination
 *  order can be checked with a current graph).
 *
 * Preconditions:
 *     The routines parseGraphicalModel() and createRandomVariableGraph()
 *     must have been called at some point, thereby creating
 *     a valid random variable graph template in rvInfoVector.
 *     It is also assumed that the order that the rvs in rvInfoVector
 *     occur in is the same order as that in the file. If a file
 *     has re-ordered variables in the file, then that will
 *     create a differnet ID (i.e., there is no attempt by this
 *     routine to normalize the order of the random variable presentation
 *     and everything is done in strict structure file (.str) order.
 *     'os' must point to a valid open output file.
 *
 * Postconditions:
 *     An "ID" is written to the file that gives sufficient
 *     detail about the graph (but not the parameters, etc.)
 *
 * Side Effects:
 *     None, other than change the file pointer in os.
 *
 * Results:
 *     None
 *
 *----------------------------------------------------------------------- 
 */
void
FileParser::writeGMId(oDataStreamFile& os)
{


  os.writeComment("Structure File Identification Information\n");
  for (unsigned i=0;i<rvInfoVector.size();i++) {

    // For each variable in template, we write out one per line the following:
    //   0. Position of in file of variable
    //   1. rv name
    //   2. rv frame number
    //   3. rv cardinality
    //   4. rv type
    //   5. number of switching parents
    //   6. list of switching parents, each a <name,offset> pair
    //   7. Number of sets of conditional parents
    //   8. For each set of conditional parents
    //        A. number of conditional parents in this set
    //        B. list of conditional parents, each a <name,offset> pair

    os.write(i,"position");
    os.write(rvInfoVector[i].name,"name");
    os.write(rvInfoVector[i].frame,"frame");
    os.write(rvInfoVector[i].rvCard,"card");
    os.write(((rvInfoVector[i].rvType == RVInfo::t_discrete)?"D":"C"),"type");
    os.write(rvInfoVector[i].switchingParents.size(),"num sp");
    for (unsigned j=0;j<rvInfoVector[i].switchingParents.size();j++) {
      os.write(rvInfoVector[i].switchingParents[j].first,"s par n");
      os.write(rvInfoVector[i].switchingParents[j].second,"s par f");
    }
    os.write(rvInfoVector[i].conditionalParents.size(),"num cps");
    for (unsigned j=0;j<rvInfoVector[i].conditionalParents.size();j++) {
      os.write(rvInfoVector[i].conditionalParents[j].size(),"num cp");
      for (unsigned k=0;k<rvInfoVector[i].conditionalParents[j].size();k++) {
	os.write(rvInfoVector[i].conditionalParents[j][k].first,"c par n");
	os.write(rvInfoVector[i].conditionalParents[j][k].second,"c par f");	
      }
    }
    os.nl();
  }
}




/*-
 *-----------------------------------------------------------------------
 * readAndVerifyGMId()
 *  A routine that reads in a graph id that was written
 *  in condensed by writeGMId(), and it verifies that
 *  the GM ID written matches the current template.
 *  See the routine 'writeGMId()' for further information.
 *
 * Preconditions:
 *     The routines parseGraphicalModel() and createRandomVariableGraph()
 *     must have been called at some point, thereby creating
 *     a valid random variable graph template in rvInfoVector.
 *     It is also assumed that the order that the rvs in rvInfoVector
 *     occur in is the same order as that in the file. If a file
 *     has re-ordered variables in the file, then that will
 *     create a differnet ID (i.e., there is no attempt by this
 *     routine to normalize the order of the random variable presentation
 *     and everything is done in strict structure file (.str) order.
 *     'is' must point to a valid open input file.
 *
 * Postconditions:
 *
 * Side Effects:
 *     none, other than moving file position.
 *
 * Results:
 *     If ID in file matches, true is returned, otherwise false.a
 *
 *----------------------------------------------------------------------- 
 */
bool
FileParser::readAndVerifyGMId(iDataStreamFile& is)
{
  // just go through and make sure everything is the same.

  // variables for checking
  int ival;
  unsigned uval;
  string nm;
  
  for (unsigned i=0;i<rvInfoVector.size();i++) {

    if (!is.read(uval)) return false;
    if (uval != i) return false;

    if (!is.read(nm)) return false;
    if (nm != rvInfoVector[i].name) return false;

    if (!is.read(uval)) return false;
    if (uval != rvInfoVector[i].frame) return false;

    if (!is.read(uval)) return false;
    if (uval != rvInfoVector[i].rvCard) return false;

    if (!is.read(nm)) return false;
    if (rvInfoVector[i].rvType == RVInfo::t_discrete) {
      if (nm != "D") return false;
    } else {
      if (nm != "C") return false;
    }

    if (!is.read(uval)) return false;    
    if (uval != rvInfoVector[i].switchingParents.size()) return false;

    for (unsigned j=0;j<rvInfoVector[i].switchingParents.size();j++) {
      if (!is.read(nm)) return false;
      if (nm != rvInfoVector[i].switchingParents[j].first) return false;


      if (!is.read(ival)) return false;
      if (ival != rvInfoVector[i].switchingParents[j].second) return false;

    }

    if (!is.read(uval)) return false;
    if (uval != rvInfoVector[i].conditionalParents.size()) return false;

    for (unsigned j=0;j<rvInfoVector[i].conditionalParents.size();j++) {

      if (!is.read(uval)) return false;      
      if (uval != rvInfoVector[i].conditionalParents[j].size()) return false;

      for (unsigned k=0;k<rvInfoVector[i].conditionalParents[j].size();k++) {

	if (!is.read(nm)) return false;
	if (nm != rvInfoVector[i].conditionalParents[j][k].first) return false;

	if (!is.read(ival)) return false;	
	if (ival != rvInfoVector[i].conditionalParents[j][k].second) return false;

      }
    }
  }
  // all checked out ok
  return true;
}







////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

GMParms GM_Parms;

int
main(int argc,char*argv[])
{
  if (argc > 1) {
    FileParser fp(argv[1]);    
    fp.parseGraphicalModel();
    fp.createRandomVariableGraph();
  } else {
    FileParser fp("-");
    fp.parseGraphicalModel();
    fp.createRandomVariableGraph();
  }


}


#endif
