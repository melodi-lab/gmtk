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

Frame_List = Frame | Frame Frame_List

Frame = "frame" ":" integer "{" RV_List "}"

RV_List = RV | RV RV_List

RV = "variable" ":" name "{" RV_Attribute_List "}"

RV_Attribute_List = Type_Attribute | Parents_Attribute

Type_Attribute = "type" ":" RV_Type ";"

Parents_Attribute =
    "switchingparents" ":" Switching_Parent_LIST ";"
    "conditionalparents" ":" Conditional_Parent_List_List  ";"

RV_Type = Discrete_RV | Continuous_RV

Discrete_RV = "discrete" ( "hidden" | "observed" integer ":" integer ) "cardinality" integer

Continuous_RV = "continuous" ("hidden" | "observed" "features" integer:integer)

Switching_Parent_LIST = "nil" | Parent_List "using" Mapping_Spec

Conditional_Parent_List_List =
    Conditional_Parent_List using Implementation |
    Conditional_Parent_List using Implementation "|" Conditional_Parent_List_List

Cond_Parent_List = "nil" | Parent_List

Parent_List = Parent | Parent "," Parent_List

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
          type: continous observed features 0:5 ;
          switchingparents: state1(0), state2(0) 
                   using mapping("state2obs6") ;
          conditionalparents: 
                 nil using mixGaussian("the_forth_gaussian");
               | state1(0) using mixGaussian mapping("gausmapping");
       }
       variable : obs2 {
          type: continous observed features 6:25  ;
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
          type: continous observed features 0:5 ;
          switchingparents: state1(0), state2(0) 
                   using mapping("state2obs6") ;
          conditionalparents: 
                 nil using mixGaussian("the_forth_gaussian");
               | state1(0) using mixGaussian mapping("gausmapping");
       }

       variable : obs2 {
          type: continous observed features 6:25  ;
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
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////



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
 * ensureNotEOF:
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
FileParser::ensureNotEOF(const char *const msg)
{
  if (tokenInfo.rc == Token_EOF) {
    fprintf(stderr,"Unexpected EOF Error: expecting %s at line %d\n",
	    msg,
	    tokenInfo.srcLine);
    error("Exiting Program");
  }
}


void
FileParser::parseError(const char* const str)
{
  fprintf(stderr,"Parse Error: %s at line %d, near (%s)\n",
	  str,
	  tokenInfo.srcLine,
	  tokenInfo.tokenStr);
  error("Exiting Program");
}



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
  parseGraphicalModel();
}


void
FileParser::parseGraphicalModel()
{

  prepareNextToken();
  ensureNotEOF("GM magic keyword");
  if (tokenInfo != "GRAPHICAL_MODEL")
    parseError("expecting key word GRAPHICAL_MODEL");

  prepareNextToken();
  ensureNotEOF("GM name");
  if (tokenInfo != Token_Identifier)
    parseError("expecting GM name ID");

  prepareNextToken();
  ensureNotEOF("frame keyword");


  if (tokenInfo != "frame")
    parseError("expecting frame keyword");
  parseFrameList();
  parseChunkSpecifier();

}


void
FileParser::parseFrameList()
{
  if (tokenInfo != "frame")
    return;
  parseFrame();
  parseFrameList();
}



void
FileParser::parseFrame()
{

  prepareNextToken();
  ensureNotEOF("frame colon");
  if (tokenInfo != Token_Colon)
    parseError("frame colon");


  prepareNextToken();
  ensureNotEOF("frame number");
  if (tokenInfo != Token_Integer)
    parseError("frame number");

  prepareNextToken();

}



void
FileParser::parseChunkSpecifier()
{
  if (tokenInfo != "chunk") 
    parseError("expecting chunk keyword");

  prepareNextToken();
  ensureNotEOF("first chunk integer");

  if (tokenInfo != Token_Integer)
    parseError("expecting first chunk integer");

  prepareNextToken();
  ensureNotEOF("chunk colon");

  if (tokenInfo != Token_Colon) 
    parseError("expecting chunk colon");

  prepareNextToken();
  ensureNotEOF("second chunk integer");

  if (tokenInfo != Token_Integer) 
    parseError("expecting second chunk integer");

}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

GMParms GM_Parms;

int
main(int argc,char*argv[])
{
  if (argc > 1) 
    FileParser fp(argv[1]);    
  else 
    FileParser fp("-");
}


#endif
