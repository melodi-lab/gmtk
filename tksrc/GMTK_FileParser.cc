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
#include "GMTK_GM.h"
#include "GMTK_GMParms.h"

VCID("$Header$");


/*
***********************************************************************
***********************************************************************

# The GM Grammar:

GM = "GRAPHICAL_MODEL" identifier Frame_List Chunk_Specifier

Frame_List = Frame Frame_List | NULL

Frame = "frame" ":" integer "{" RV_List "}"

RV_List = RV RV_List | NULL

RV = "variable" ":" name "{" RV_Attribute "}"

RV_Attribute_List =
        RV_Attribute RV_Attribute_List | NULL

RV_Attribute = Type_Attribute | Parents_Attribute

Type_Attribute = "type" ":" RV_Type ";"

Parents_Attribute =
      ( "switchingparents" ":" Switching_Parent_LIST ";" )
   |  ( "conditionalparents" ":" Conditional_Parent_Spec_List  ";" )
    # Semanatics requires that we have both
    # switchingparents & conditionalparents in a RV. Note that
    # the parse grammer allows this not to be the case.

RV_Type = Discrete_RV_TYPE | Continuous_RV_TYPE

Discrete_RV_TYPE = 
      "discrete" 
      ( "hidden" | "observed" integer ":" integer ) 
     "cardinality" integer

Continuous_RV_TYPE = 
       "continuous" 
       ("hidden" | "observed" integer:integer)

Switching_Parent_LIST = "nil" | Parent_List "using" Mapping_Spec

Conditional_Parent_Spec_List =
    Conditional_Parent_Spec "|" Conditional_Parent_Spec_List
  | Conditional_Parent_Spec

Conditional_Parent_Spec = Conditional_Parent_List using Implementation

Cond_Parent_List = "nil" | Parent_List

Parent_List = 
       Parent "," Parent_List
     | Parent

Parent = identifier "(" integer ")"

Implementation = DiscreteImplementation | ContinousImplementation

DiscreteImplementation = "MDCPT" | "MSCPT"

ContinousImplementation = DirectContObsDist | MappingToAContObsDist

DirectContObsDist = ContObsDistType "(" List_Index ")"
       # direct cont. observation dist is used if conditional parents 
       # are nil
       # in which case we select only one dist. This
       # is analogous to the discrete case where you directly
       # select a CPT with the appropriate parents

MappingToAContObsDist = ContObsDistType Mapping_Spec
       # this is when we have multiple conditional parents,
       # and we need another decision tree to map
       # from the conditional parents values to the appropriate
       # distribution.

ContObsDistType = "mixGaussian" | "gausSwitchMixGaussian" 
  | "logitSwitchMixGaussian" | "mlpSwitchMixGaussian"

Mapping_Spec = "mapping" "(" List_Index ")"
     # A Mapping_Spec always indexes into one of the decision trees.
     # The integer (or string) is used to index into a table
     # of decision trees to choose the decision tree
     # that will map from the switching parents to one of the
     # conditional parent lists.

Chunk_Specifier = "chunk"  integer ":" integer

List_Index = integer | string

======================================================================


Example of a grammatical GM file
-----------------------------------

# Actual model definition
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
          type: continous observed 0:5 ;
          switchingparents: state1(0), state2(0) 
                   using mapping("state2obs6") ;
          conditionalparents: 
                 nil using mixGaussian("the_forth_gaussian");
               | state1(0) using mixGaussian mapping("gausmapping");
       }
       variable : obs2 {
          type: continous observed 6:25  ;
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
          type: continous observed 0:5 ;
          switchingparents: state1(0), state2(0) 
                   using mapping("state2obs6") ;
          conditionalparents: 
                 nil using mixGaussian("the_forth_gaussian");
               | state1(0) using mixGaussian mapping("gausmapping");
       }

       variable : obs2 {
          type: continous observed 6:25  ;
          switchingparents: nil ;
          conditionalparents: 
                state1(0) using mlpSwitchMixGaussian
                          mapping("gaussmapping3") ;
       }
}

chunk: 1:1



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
  tokenInfo.srcLine = 1;
  tokenInfo.srcChar = 1;
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
  parseRandomVariable();
  parseRandomVariableList();
}


void
FileParser::parseRandomVariable()
{
  ensureNotAtEOF(KW_Variable);
  if (tokenInfo != KW_Variable)
    parseError(KW_Variable);
  consumeToken();

  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseError(":");
  consumeToken();

  ensureNotAtEOF("variable name");
  if (tokenInfo != TT_Identifier)
    parseError("variable name");
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
	   tokenInfo == KW_Conditionalparents)
    return parseRandomVariableParentAttribute();
  else
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
  consumeToken();

  ensureNotAtEOF("variable disposition (hidden|discrete)");
  if (tokenInfo == KW_Hidden) {
    // hidden discrete RV
    consumeToken();
  } else if (tokenInfo == KW_Observed) {
    // observed discrete RV    
    consumeToken();
    
    ensureNotAtEOF("first feature range");
    if (tokenInfo != TT_Integer)
      parseError("first feature range");
    consumeToken();

    ensureNotAtEOF("feature range separator");
    if (tokenInfo != TT_Colon)
      parseError("feature range separator");
    consumeToken();
    
    ensureNotAtEOF("second feature range");
    if (tokenInfo != TT_Integer)
      parseError("second feature range");
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
  consumeToken();
}


void
FileParser::parseRandomVariableContinuousType()
{
  ensureNotAtEOF(KW_Continous);
  if (tokenInfo != KW_Continous)
    parseError(KW_Continous);
  consumeToken();

  ensureNotAtEOF("variable disposition (hidden|discrete)");
  if (tokenInfo == KW_Hidden) {
    // hidden discrete RV
    consumeToken();
  } else if (tokenInfo == KW_Observed) {
    // observed discrete RV    
    consumeToken();
    
    ensureNotAtEOF("first feature range");
    if (tokenInfo != TT_Integer)
      parseError("first feature range");
    consumeToken();

    ensureNotAtEOF("feature range separator");
    if (tokenInfo != TT_Colon)
      parseError("feature range separator");
    consumeToken();
    
    ensureNotAtEOF("second feature range");
    if (tokenInfo != TT_Integer)
      parseError("second feature range");
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

    parseSwitchingParentList();

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
FileParser::parseSwitchingParentList()
{
}

void
FileParser::parseConditionalParentSpecList()
{
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
  consumeToken();

  ensureNotAtEOF("chunk colon");
  if (tokenInfo != TT_Colon) 
    parseError("expecting chunk colon");
  consumeToken();

  ensureNotAtEOF("second chunk integer");
  if (tokenInfo != TT_Integer) 
    parseError("expecting second chunk integer");
  consumeToken();

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
