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
#include "GMTK_DiscreteVariable.h"
#include "GMTK_ContinuousVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMParms.h"


VCID("$Header$");


/*
***********************************************************************
***********************************************************************

# The GM Grammar:

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
      ( "hidden" | "observed" integer ":" integer ) 
     "cardinality" integer

RandomVariableContinuousType = 
       "continuous" 
       ("hidden" | "observed" integer:integer)

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

Implementation = DiscreteImplementation | ContinousImplementation

DiscreteImplementation = ( "MDCPT" | "MSCPT" | "MTCPT" )  "(" ListIndex ")"

ContinousImplementation = ContObsDistType
        (
             "(" ListIndex ")"
           |
              MappingSpec
        )
       # A ContinousImplementation has two cases.
       # 1) the first case is used if conditional parents 
       # are nil in which case we select only one dist. This
       # is analogous to the discrete case where you directly
       # select a CPT with the appropriate parents
       # 2) in the second case, this is when we have multiple
       #  conditional parents, and we need another decision tree to map
       # from the conditional parents values to the appropriate
       # distribution.


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
# Actual model definition (that parsed).
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

Notes on the semantic stuff:

o   Check for unconnected RVs?
o   Allow non DAGs (untimately for loopy inference) but include a routine 
    that provides a DAG check.
o   


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
    "GRAPHICAL_MODEL"
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

  if (file == NULL)
    error("FileParser::FileParser, can't open NULL file");
  if (!strcmp("-",file))
    yyin = stdin;
  else {
    if ((yyin = fopen (file,"r")) == NULL)
      error("FileParser::FileParser, can't open file (%s)",file);
  }

}


/*-
 *-----------------------------------------------------------------------
 * prepareNextToken
 *   prepares the next input token and fills the tokeninfo structure.
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
  printf("Consuming token (%s) of type %d from src line (%d)\n",
	 tokenInfo.tokenStr,tokenInfo.tokenType,tokenInfo.srcLine);
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
  fprintf(stderr,"Parse Error: %s at line %d, near (%s)\n",
	  str,
	  tokenInfo.srcLine,
	  tokenInfo.tokenStr);
  error("Exiting Program");
}

void
FileParser::parseError(const TokenKeyword kw)
{
  fprintf(stderr,"Parse Error: expecting keyword (%s) at line %d, near (%s)\n",
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

  parseChunkSpecifier();

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
  map < pair<string,unsigned> , RVInfo* >::iterator it;
  it = nameRVmap.find(pair<string,unsigned>(curRV.name,curRV.frame));
  if (it != nameRVmap.end()) {
    error("Error: random variable (%s) at frame (%d) defined twice, "
	  "on both line %d and %d\n",curRV.name.c_str(),
	  curRV.frame,
	  curRV.fileLineNumber,
	  (*it).second->fileLineNumber);
  }
  

  // everything looks ok, insert it in our tables.
  rvInfoVector.push_back(curRV);
  // add to map with invalid pointers for now.
  nameRVmap[
	    pair<string,unsigned>(curRV.name,curRV.frame)
  ] = &rvInfoVector[rvInfoVector.size()-1];

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
  else if (tokenInfo == KW_Continous)
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
    
    ensureNotAtEOF("first feature range");
    if (tokenInfo != TT_Integer)
      parseError("first feature range");
    curRV.rvFeatureRange.start = tokenInfo.int_val;
    consumeToken();

    ensureNotAtEOF("feature range separator");
    if (tokenInfo != TT_Colon)
      parseError("feature range separator");
    consumeToken();
    
    ensureNotAtEOF("second feature range");
    if (tokenInfo != TT_Integer)
      parseError("second feature range");
    curRV.rvFeatureRange.stop = tokenInfo.int_val;
    if (curRV.rvFeatureRange.stop < curRV.rvFeatureRange.start)
      parseError("first range num must be < second range num");
    curRV.rvFeatureRange.filled = true;
    consumeToken();

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
  ensureNotAtEOF(KW_Continous);
  if (tokenInfo != KW_Continous)
    parseError(KW_Continous);
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
    curRV.rvFeatureRange.start = tokenInfo.int_val;
    consumeToken();

    ensureNotAtEOF("feature range separator");
    if (tokenInfo != TT_Colon)
      parseError("feature range separator");
    consumeToken();
    
    ensureNotAtEOF("second feature range");
    if (tokenInfo != TT_Integer)
      parseError("second feature range");
    curRV.rvFeatureRange.stop = tokenInfo.int_val;
    curRV.rvFeatureRange.filled = true;
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
      curRV.discImplementations.push_back(RVInfo::di_MDCPT);
    else if (tokenInfo == KW_MSCPT)
      curRV.discImplementations.push_back(RVInfo::di_MSCPT);
    else // (tokenInfo == KW_MTCPT)
      curRV.discImplementations.push_back(RVInfo::di_MTCPT);
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
    curRV.contImplementations.push_back(RVInfo::ci_mixGaussian);
    consumeToken();
  } else if (tokenInfo == KW_GausSwitchMixGaussian) {
    curRV.contImplementations.push_back(RVInfo::ci_gausSwitchMixGaussian);
    consumeToken();
  }  else if (tokenInfo == KW_LogitSwitchMixGaussian) {
    curRV.contImplementations.push_back(RVInfo::ci_logitSwitchMixGaussian);
    consumeToken();
  }  else if (tokenInfo == KW_MlpSwitchMixGaussin) {
    curRV.contImplementations.push_back(RVInfo::ci_mlpSwitchMixGaussian);
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
  firstChunkframe = (unsigned)tokenInfo.int_val;
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
  lastChunkframe = (unsigned)tokenInfo.int_val;
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
}


//
// 
// Check the consistency of the RV information
// that couldn't be checked while it was parsing (if anything)
void
FileParser::RVInfo::checkConsistency()
{
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
  } else {
    FileParser fp("-");
    fp.parseGraphicalModel();
  }

}


#endif
