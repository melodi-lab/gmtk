/*-
 * GMTK_MDCPT.cc
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
#include "GMTK_GMParms.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_MixGaussians.h"


VCID("$Header$");

extern "C" {
  FILE     *popen(const char *, const char *);
  int pclose(FILE *stream);
};

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

DiscreteImplementation = ( "MDCPT" | "MSCPT" | "MTCPT" )  "(" ListIndex ")"

ContinuousImplementation = ContObsDistType
        (
             "(" ListIndex ")"
           |
              MappingSpec
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
       #           MappingSpec
       #
       # this is when we have multiple
       # conditional parents, and we need another decision tree to map
       # from the conditional parents values to the appropriate
       # distribution. Therefore we use the mapping syntax.
       # In this case, the 'MappingSpec' must refer to one of
       # the decision trees. 
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
          conditionalparents: nil using MDCPT("initcpt") ;
        }
     variable : phone {
          type: discrete hidden cardinality 4 ;
          switchingparents : nil ;
          conditionalparents :
                     word(0) using MDCPT("dep_on_word") ;
        }
     variable : state1 {
          type: discrete hidden cardinality 4 ;
          switchingparents: phone(0) using mapping("phone2state1") ;
          conditionalparents :
                  nil using MSCPT("f1")
                | nil using MDCPT("f2") ;
       }
       variable : state2 {
          type: discrete hidden cardinality 4 ;
          switchingparents: nil ;
          conditionalparents:
                     nil using MDCPT("f5") ;
       }

       variable : obs1 {
          type: continuous observed 0:5 ;
          switchingparents: state1(0), state2(0)
                   using mapping("state2obs6") ;
          conditionalparents: 
                 nil using mixGaussian("the_forth_gaussian")
               | state1(0) using mixGaussian mapping("gausmapping");
       }
       variable : obs2 {
          type: continuous observed 6:25  ;
          switchingparents: nil ;
          conditionalparents: 
                state1(0) using mlpSwitchMixGaussian
                          mapping("gaussmapping3") ;
       }
}

frame:1 {

     variable : word {
          type: discrete hidden cardinality 3 ;
          switchingparents: nil ;
          conditionalparents: word(-1) using MDCPT("wordbigram") ;
        }

     variable : phone {
          type: discrete hidden cardinality 4 ;
          switchingparents : nil ;
          conditionalparents :
                phone(-1),word(0) using MDCPT("dep_on_word_phone") ;
        }

     variable : state1 {
          type: discrete hidden cardinality 4 ;
          switchingparents: phone(0) using mapping("phone2state1") ;
          conditionalparents :
	          # in this first case, state1 is dep on prev time
                  state1(-1) using MSCPT("f3")
	          # in this second case, it is cond. indep. of prev. time.
                | nil using MDCPT("f2") ;
       }

       variable : state2 {
          type: discrete hidden cardinality 4 ;
          switchingparents: nil ;
          conditionalparents:
	             # this is a sep. hidden markov chain
                     state2(-1) using MDCPT("f9") ;
       }

       # like in the first frame, the obs. are only dep.
       # on RVs from the current frame.
       variable : obs1 {
          type: continuous observed 0:5 ;
          switchingparents: state1(0), state2(0) 
                   using mapping("state2obs6") ;
          conditionalparents: 
                 nil using mixGaussian("the_forth_gaussian")
               | state1(0) using mixGaussian mapping("gausmapping");
       }

       variable : obs2 {
          type: continuous observed 6:25  ;
          switchingparents: nil ;
          conditionalparents: 
                state1(0) using mlpSwitchMixGaussian
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


////////////////////////////////////////////////////////////////////
//        Internal Keyword Symbols
////////////////////////////////////////////////////////////////////


const vector<string> 
FileParser::KeywordTable = FileParser::fillKeywordTable();

const vector<string> 
FileParser::fillKeywordTable()
{
  // ******************************************************************
  // This table must always be consistent with the enum TokenKeyword.
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
    "MDCPT",
    "MSCPT",
    "MTCPT",
    "mixGaussian",
    "gausSwitchMixGaussian",
    "logitSwitchMixGaussian",
    "mlpSwitchMixGaussian",
    "chunk",
    "GRAPHICAL_MODEL",
    "value"
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
//        RVInfo functions
////////////////////////////////////////////////////////////////////


FileParser::RVInfo::RVInfo(const RVInfo& v)
{
  frame = v.frame;
  fileLineNumber = v.fileLineNumber;
  name = v.name;
  rvType = v.rvType;
  rvDisp = v.rvDisp;
  rvCard = v.rvCard;
  rvFeatureRange = v.rvFeatureRange;
  switchingParents = v.switchingParents;
  switchMapping = v.switchMapping;
  conditionalParents = v.conditionalParents;
  discImplementations = v.discImplementations;
  contImplementations = v.contImplementations;
  listIndices = v.listIndices;
  rv = v.rv;
}


//
// 
// Check the consistency of the RV information
// that couldn't be checked while it was parsing (if anything)
void
FileParser::RVInfo::checkConsistency()
{
  error("not implemented");
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
FileParser::FileParser(const char*const file)
{
  FILE* f;
  if (file == NULL)
    error("FileParser::FileParser, can't open NULL file");
  if (!strcmp("-",file)) {
    f = ::popen("cpp","r");
    if (f == NULL) {
      error("ERROR: unable to open with standard input structure file");
    }
  }  else {
    if ((f = ::fopen(file,"r")) == NULL) {
      error("ERROR: unable to open file (%s) for reading",file);
    }
    fclose(f);

    string str = (string)"cpp " + (string)file;
    f = ::popen(str.c_str(),"r");    
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
  fprintf(stderr,"Parse Error in file '%s': %s at line %d, near (%s)\n",
	  fileNameParsing.c_str(),
	  str,
	  tokenInfo.srcLine,
	  tokenInfo.tokenStr);
  error("Exiting Program");
}

void
FileParser::parseError(const TokenKeyword kw)
{
  fprintf(stderr,"Parse Error in file '%s': expecting keyword (%s) at line %d, near (%s)\n",
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

  printf("Found %d random variables.\n",
	 rvInfoVector.size());
  for (unsigned i=0;i<rvInfoVector.size();i++ ) {
    printf("RV(%d) is named (%s)\n",
	   i,rvInfoVector[i].name.c_str());
  }


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
  map < rvParent , unsigned >::iterator it;
  it = nameRVmap.find(rvParent(curRV.name,curRV.frame));
  if (it != nameRVmap.end()) {
    error("Error: random variable (%s) at frame (%d) defined twice, "
	  "on both line %d and %d\n",curRV.name.c_str(),
	  curRV.frame,
	  curRV.fileLineNumber,
	  rvInfoVector[(*it).second].fileLineNumber);
  }
  

  // everything looks ok, insert it in our tables.
  rvInfoVector.push_back(curRV);
  // add to map
  nameRVmap[
	    rvParent(curRV.name,curRV.frame)
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
  rvParent p;

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
  // make sure we are not pointing to ourselves.
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
    parseError("discrete implementation");
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

    parseMappingSpec();
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
  // first create the RV objects
  for (unsigned i=0;i<rvInfoVector.size();i++) {
    if (rvInfoVector[i].rvType == RVInfo::t_discrete) {
      DiscreteRandomVariable*rv = 
	new DiscreteRandomVariable(rvInfoVector[i].name,
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
	new ContinuousRandomVariable(rvInfoVector[i].name);
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
  }

  // now set up all the parents of each random variable.
  for (unsigned i=0;i<rvInfoVector.size();i++) {

    unsigned frame = rvInfoVector[i].frame;

    // build the switching parents list.
    vector<RandomVariable *> sparents;
    for (unsigned j=0;j<rvInfoVector[i].switchingParents.size();j++) {

      rvParent pp(rvInfoVector[i].switchingParents[j].first,
		  frame+rvInfoVector[i].switchingParents[j].second);

      ////////////////////////////////////////////////////
      // Make sure the rv at the time delta from the current
      // frame exists.
      if (nameRVmap.find(pp) == nameRVmap.end())
	error("Error: parent random variable \"%s\" at frame %d does not exist\n",
	      pp.first.c_str(),pp.second);

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
    
    // now build continous parent list
    vector<vector<RandomVariable * > > cpl(rvInfoVector[i].conditionalParents.size());
    for (unsigned j=0;j<rvInfoVector[i].conditionalParents.size();j++) {
      for (unsigned k=0;k<rvInfoVector[i].conditionalParents[j].size();k++) {

	rvParent pp(rvInfoVector[i].conditionalParents[j][k].first,
		    frame+rvInfoVector[i].conditionalParents[j][k].second);

	////////////////////////////////////////////////////
	// Make sure the rv at the time delta from the current
	// frame exists.
	if (nameRVmap.find(pp) == nameRVmap.end())
	  error("Error: parent random variable \"%s\" at frame %d does not exist\n",
		pp.first.c_str(),pp.second);

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

      rvParent pp(rvInfoVector[i].switchingParents[j].first,
		  rvInfoVector[i].frame
		  +rvInfoVector[i].switchingParents[j].second);

      ////////////////////////////////////////////////////
      // Make sure the rv at the time delta from the current
      // frame exists.
      if (nameRVmap.find(pp) == nameRVmap.end())
	error("Error: parent random variable \"%s\" at frame %d does not exist\n",
	      pp.first.c_str(),pp.second);

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

	rvParent pp(rvInfoVector[i].conditionalParents[j][k].first,
		    rvInfoVector[i].frame		    
		    +rvInfoVector[i].conditionalParents[j][k].second);

	////////////////////////////////////////////////////
	// Make sure the rv at the time delta from the current
	// frame exists.
	if (nameRVmap.find(pp) == nameRVmap.end())
	  error("Error: parent random variable \"%s\" at frame %d does not exist\n",
		pp.first.c_str(),pp.second);
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
FileParser::associateWithDataParams(bool allocateIfNotThere)
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
	      if (!allocateIfNotThere) {
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


		  rvParent pp(rvInfoVector[i].conditionalParents[j][k].first,
			      rvInfoVector[i].frame
			      +rvInfoVector[i].conditionalParents[j][k].second);
		  map < rvParent , unsigned >::iterator it;
		  it = nameRVmap.find(pp);
		  if (it == nameRVmap.end()) {
		    // this really shouldn't happen at this point since
		    // it should have been checked somewhere else,
		    // but we include the check nonetheless
		    error("Error: parent random variable \"%s\" at" 
			  "frame %d does not exist\n",
			  rvInfoVector[i].conditionalParents[j][k].first.c_str(),
			  rvInfoVector[i].conditionalParents[j][k].second);
		  }
		  mdcpt->setNumCardinality(k,
					   rvInfoVector[(*it).second].rvCard);
		}
		mdcpt->setNumCardinality(rvInfoVector[i].conditionalParents[j].size(),
					 rvInfoVector[i].rvCard);
		mdcpt->allocateBasicInternalStructures();

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
	  }

	  // now check to make sure this cpt matches this
	  // number of parents.
	  if (cpts[j]->numParents() != 
	      rvInfoVector[i].conditionalParents[j].size()) {
	    error("Error: RV \"%s\" at frame %d (line %d), num parents cond. %d different than required by MDCPT \"%s\".\n",
		  rvInfoVector[i].name.c_str(),
		  rvInfoVector[i].frame,
		  rvInfoVector[i].fileLineNumber,
		  j,
		  cpts[j]->name().c_str());
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
		if (!allocateIfNotThere) {
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent MSCPT \"%s\" doesn't exist\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].nameIndex.c_str());
		}
		else {
		  // allocate the MSCPT with name and install it.
		  MSCPT* mscpt = new MSCPT();
		  mscpt->setName(rvInfoVector[i].listIndices[j].nameIndex);
		  mscpt->
		    setNumParents
		    (rvInfoVector[i].conditionalParents[j].size());

		  for (unsigned k=0;k<rvInfoVector[i].conditionalParents[j].size();k++) {


		    rvParent pp(rvInfoVector[i].conditionalParents[j][k].first,
				rvInfoVector[i].frame
				+rvInfoVector[i].conditionalParents[j][k].second);
		    map < rvParent , unsigned >::iterator it;
		    it = nameRVmap.find(pp);
		    if (it == nameRVmap.end()) {
		      // this really shouldn't happen at this point since
		      // it should have been checked somewhere else,
		      // but we include the check nonetheless
		      error("Error: parent random variable \"%s\" at" 
			    "frame %d does not exist\n",
			    rvInfoVector[i].conditionalParents[j][k].first.c_str(),
			    rvInfoVector[i].conditionalParents[j][k].second);
		    }
		    mscpt->setNumCardinality(k,
					     rvInfoVector[(*it).second].rvCard);
		  }
		  mscpt->setNumCardinality(rvInfoVector[i].conditionalParents[j].size(),
					   rvInfoVector[i].rvCard);
		  mscpt->allocateBasicInternalStructures();

		  GM_Parms.msCpts.push_back(mscpt);
		  GM_Parms.msCptsMap[rvInfoVector[i].listIndices[j].nameIndex]
		    = GM_Parms.msCpts.size()-1;
		  cpts[j] = mscpt;
		}
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
	    }

	    // now check to make sure this cpt matches this
	    // number of parents.
	    if (cpts[j]->numParents() != 
		rvInfoVector[i].conditionalParents[j].size()) {
	      error("Error: RV \"%s\" at frame %d (line %d), num parents cond. %d different than required by MSCPT \"%s\".\n",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    j,
		    cpts[j]->name().c_str());
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
		if (!allocateIfNotThere) {
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent MTCPT \"%s\" doesn't exist\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].nameIndex.c_str());
		}
		else {
		  // allocate the MTCPT with name and install it.
		  MTCPT* mtcpt = new MTCPT();
		  mtcpt->setName(rvInfoVector[i].listIndices[j].nameIndex);
		  mtcpt->
		    setNumParents
		    (rvInfoVector[i].conditionalParents[j].size());

		  for (unsigned k=0;k<rvInfoVector[i].conditionalParents[j].size();k++) {


		    rvParent pp(rvInfoVector[i].conditionalParents[j][k].first,
				rvInfoVector[i].frame
				+rvInfoVector[i].conditionalParents[j][k].second);
		    map < rvParent , unsigned >::iterator it;
		    it = nameRVmap.find(pp);
		    if (it == nameRVmap.end()) {
		      // this really shouldn't happen at this point since
		      // it should have been checked somewhere else,
		      // but we include the check nonetheless
		      error("Error: parent random variable \"%s\" at" 
			    "frame %d does not exist\n",
			    rvInfoVector[i].conditionalParents[j][k].first.c_str(),
			    rvInfoVector[i].conditionalParents[j][k].second);
		    }
		    mtcpt->setNumCardinality(k,
					     rvInfoVector[(*it).second].rvCard);
		  }
		  mtcpt->setNumCardinality(rvInfoVector[i].conditionalParents[j].size(),
					   rvInfoVector[i].rvCard);
		  mtcpt->allocateBasicInternalStructures();

		  GM_Parms.mtCpts.push_back(mtcpt);
		  GM_Parms.mtCptsMap[rvInfoVector[i].listIndices[j].nameIndex]
		    = GM_Parms.mtCpts.size()-1;
		  cpts[j] = mtcpt;
		}
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
	    }

	    // now check to make sure this cpt matches this
	    // number of parents.
	    if (cpts[j]->numParents() != 
		rvInfoVector[i].conditionalParents[j].size()) {
	      error("Error: RV \"%s\" at frame %d (line %d), num parents cond. %d different than required by MTCPT \"%s\".\n",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    j,
		    cpts[j]->name().c_str());
	    }

	} else {
	  // Again, this shouldn't happen. If it does, something is wrong
	  // with the parser code or some earlier code, and it didn't correctly
	  // set the CPT type. 
	  assert ( 0 );
	}
      }
      // now add the cpts to the rv
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
	  }
	} else {
	  // there are > 0 conditional parents for this
	  // set of switching values. The index should
	  // specify a DT.
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
	      rv->conditionalGaussians[j].dtMapper =
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
		rv->conditionalGaussians[j].dtMapper =
		  GM_Parms.dts[rvInfoVector[i].listIndices[j].intIndex];
	      }
	  }
	}
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
