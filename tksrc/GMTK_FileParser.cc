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

// include all random variable classes.
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_Sw_HidDiscRV.h"
#include "GMTK_ScPnSh_HidDiscRV.h"
#include "GMTK_ScPnSh_Sw_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_Sw_ObsDiscRV.h"
#include "GMTK_ScPnSh_ObsDiscRV.h"
#include "GMTK_ScPnSh_Sw_ObsDiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_ObsContRV.h"
#include "GMTK_Sw_ObsContRV.h"
#include "GMTK_ScPnSh_ObsContRV.h"
#include "GMTK_ScPnSh_Sw_ObsContRV.h"

#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_NGramCPT.h"
#include "GMTK_FNGramCPT.h"
#include "GMTK_USCPT.h"
#include "GMTK_VECPT.h"
#include "GMTK_LatticeNodeCPT.h"
#include "GMTK_LatticeEdgeCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_GraphicalModel.h"
#include "GMTK_RVInfo.h"

VCID("$Header$")

#ifndef DECLARE_POPEN_FUNCTIONS_EXTERN_C
extern "C" {
  //   FILE     *popen(const char *, const char *) __THROW;
  //   int pclose(FILE *stream) __THROW;
}
#endif

#define TRIFILE_END_OF_ID_STRING "@@@!!!TRIFILE_END_OF_ID_STRING!!!@@@"


/*
***********************************************************************
***********************************************************************

The GM Grammar: 

     This is a valid grammer for the parser below. The grammar
     specifies the language that is excepted by the parser, but other
     "semantic" errors are checked either after parsing, or sometimes
     while parsing when possible.

     NOTE: Please keep this grammer up to date if any changes are made
     to the parser.

GM = "GRAPHICAL_MODEL" identifier FrameList ChunkSpecifier

FrameList = Frame FrameList | NULL

Frame = "frame" ":" integer "{" RandomVariableList "}"

RandomVariableList = RandomVariable RandomVariableList | NULL

RandomVariable = "variable" ":" name "{" RandomVariableAttributeList "}"
               |  "factor" ":" name "{" FactorAttributeList "}"


RandomVariableAttributeList =
        RandomVariableAttribute RandomVariableAttributeList | NULL

RandomVariableAttribute = TypeAttribute | 
                          EliminationHintAttribute | 
                          WeightAttribute | 
                          ParentsAttribute

WeightAttribute = "weight" : WeightAttributeSpecList ";"

WeightAttributeSpecList =
    WeightAttributeSpec "|" WeightAttributeSpecList
  | WeightAttributeSpec

WeightAttributeSpec =  WeightOptionList
                     | "nil"   

WeightOptionList = 
     ( WeightType WeightOption ) WeightOptionList
   | ( WeightType WeightOption )

WeightType =  "scale" | "penalty" | "shift"

WeightOption = ( "value" number )
             | ( integer ":" integer )

EliminationHintAttribute = "elimination_hint" : number ";"


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

DiscreteImplementation = ( "DenseCPT" | "SparseCPT" | "DeterministicCPT" | 
                           "NGramCPT" | "FNGramCPT" )  "(" ListIndex ")"

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


ContObsDistType = "mixture" | "gausSwitchMixture" 
  | "logitSwitchMixture" | "mlpSwitchMixture"

MappingSpec = "mapping" "(" ListIndex ")"
     # A MappingSpec always indexes into one of the decision trees.
     # The integer (or string) is used to index into a table
     # of decision trees to choose the decision tree
     # that will map from the switching parents to one of the
     # conditional parent lists.

ChunkSpecifier = "chunk"  integer ":" integer

ListIndex = integer | string


number = integer | floating_point_value

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
                 nil using mixture("the_forth_gaussian_mixture")
               | state1(0) using mixture collection("gaussianCollection") mapping("gausmapping");
       }
       variable : obs2 {
          type: continuous observed 6:25  ;
          switchingparents: nil ;
          conditionalparents:
                state1(0) using mlpSwitchMixture collection("mlpgaussianCollection")
                          mapping("gaussmapping3") ;
       }
}

frame:1 {

     variable : word {
          type: discrete hidden cardinality 3 ;
          switchingparents: nil ;
          conditionalparents: word(-1) using DenseCPT("wordbigram") ;
          weight: penalty -10.0 ; 
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
                 nil using mixture("the_forth_gaussian_mixture")
               | state1(0) using mixture collection("gaussianCollection") mapping("gausmapping");
       }

       variable : obs2 {
          type: continuous observed 6:25  ;
          switchingparents: nil ;
          conditionalparents:
                state1(0) using mlpSwitchMixture collection("gaussianCollection")
                          mapping("gaussmapping3") ;
       }
}

chunk: 1:1

------------------------------------------------------------

Example factor syntax for reference.

   factor: firstFactor {
      variables: fooA(0),fooA(-1);

      // only one type of constraint can be defined at a time??
      symmetricConstraint:
	allVarsEqual;
        allVarsUnequal;
        varsNotEqual;
        varsSumTo(n);  // only values that sum to n
        varsMultiplyTo(n); // only values that multiply to n
        varsSumMod(m,n);  // only values such that (sum(vars) % m) = n
                     // e.g., force even parity done by sumMod(2,0)
        varsSatisfy using mapping("bla");
           // where mapping is a DT

     directionalConstraint: fooA(0) = functionOf(fooA(-1)) using mapping("bla");
        // where 'bla' is the name of a decision tree

     softConstraint: using table("foo"); // table is a new object, EM learnable, will also
                                            be used for TableCPT.
     softConstraint: using logLinear("bar"); // logLinear is a new loglinear object, EM learnable
  }


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
  // This table must always be order consistent with the enum TokenKeyword in 
  // the .h file and consistent with 'keyword' in the .lex file.
  // ******************************************************************
  const char*const kw_table[] = {
    /* 0  */ "frame",
    /* 1  */ "variable",
    /* 2  */ "type",
    /* 3  */ "cardinality",
    /* 4  */ "switchingparents",
    /* 5  */ "conditionalparents",
    /* 6  */ "discrete",
    /* 7  */ "continuous",
    /* 8  */ "hidden",
    /* 9  */ "observed",
    /* 10 */ "nil",
    /* 11 */ "using",
    /* 12 */ "mapping",
    /* 13 */ "collection",
    /* 14 */ "DenseCPT",
    /* 15 */ "SparseCPT",
    /* 16 */ "DeterministicCPT",
    /* 17 */ "NGramCPT",
    /* 18 */ "FNGramCPT",
    /* 19 */ "mixture",
    /* 20 */ "gausSwitchMixture",
    /* 21 */ "logitSwitchMixture",
    /* 22 */ "mlpSwitchMixture",
    /* 23 */ "chunk",
    /* 24 */ "GRAPHICAL_MODEL",
    /* 25 */ "value",
    /* 26 */ "weight",
    /* 27 */ "scale",
    /* 28 */ "penalty",
    /* 29 */ "shift",
    /* 30 */ "elimination_hint",
    /* 31 */ "frameNum",
    /* 32 */ "numFrames",
    /* 33 */ "segmentNum",
    /* 34 */ "numSegments",
    /* 35 */ "VirtualEvidenceCPT",
    /* 36 */ "LatticeNodeCPT",
    /* 37 */ "LatticeEdgeCPT",
    /* 38 */ "factor",
    /* 39 */ "variables",
    /* 40 */ "symmetricConstraint",
    /* 41 */ "allVarsEqual",
    /* 42 */ "allVarsUnequal",
    /* 43 */ "varsNotEqual",
    /* 44 */ "varsSumTo",
    /* 45 */ "varsMultiplyTo",
    /* 46 */ "varsSumMod",
    /* 47 */ "varsSatisfy",
    /* 48 */ "directionalConstraint",
    /* 49 */ "functionOf",
    /* 50 */ "softConstraint",
    /* 51 */ "table",
    /* 52 */ "logLinear",
    /* 53 */ "emarfNum"
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
  string cppCommand = CPP_Command();
  if (cppCommandOptions != NULL)
    cppCommand = cppCommand + string(" ") + cppCommandOptions;

  if (!strcmp("-",file)) {
    f = ::popen(cppCommand.c_str(),"r");
    if (f == NULL) {
      error("ERROR: unable to open with standard input structure file");
    }
  }  else {

    // check that file exists.
    if ((f = ::fopen(file,"r")) == NULL) {
      error("ERROR: unable to open file (%s) for reading",file);
    }
    fclose(f);

    // add path of file to include directory paths.
    string path = file;
    unsigned long slashPos = path.rfind("/");
    if (slashPos != string::npos) {
      // then '/' is found
      cppCommand = cppCommand + " -I" + path.substr(0,slashPos);
    }
    // Lastly, add CWD to default CPP command options for include files
    // (i.e., we look for include files in CWD only if all previous
    // ones fail, cpp has this behavior.
    cppCommand = cppCommand + " -I.";

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
FileParser::parseErrorExpecting(const char* const str)
{
  fprintf(stderr,"Parse Error in file '%s': expecting %s at or before line %d, near (%s)\n",
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
  factorList.clear();

  // prepare the first lookahead token
  prepareNextToken();

  // now start parsing.
  ensureNotAtEOF("GM magic keyword");
  if (tokenInfo != KW_GRAPHICAL_MODEL)
    parseError(KW_GRAPHICAL_MODEL);
  consumeToken();

  ensureNotAtEOF("GM name");
  if (tokenInfo != TT_Identifier)
    parseErrorExpecting("GM name identifier");
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

  pclose(yyin);

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
    parseErrorExpecting("frame colon");
  consumeToken();

  ensureNotAtEOF("frame number");
  if (tokenInfo != TT_Integer)
    parseErrorExpecting("frame number");
  if (tokenInfo.int_val != (curFrame+1)) {
    parseError("frame number out of order");
  }
  curFrame++;
  consumeToken();


  ensureNotAtEOF("open frame {");
  if (tokenInfo != TT_LeftBrace)
    parseErrorExpecting("open frame {");
  consumeToken();

  parseFrameEntryList();

  ensureNotAtEOF("close frame }");
  if (tokenInfo != TT_RightBrace)
    parseErrorExpecting("close frame }");
  consumeToken();
}


void
FileParser::parseFrameEntryList()
{
  if (tokenInfo != KW_Variable && tokenInfo != KW_Factor)
    return;

  if (tokenInfo == KW_Variable) {
    // the current frame entry must be a random variable definition.

    curRV.clear();
    parseRandomVariable();

    ////////////////////////////////////////////////////////////
    // A number of "semantic" errors can be checked right
    // here before going on, so we do that.
    ///////////////////////////////////////////////////////////////

    //////////////////////////////
    // make sure that there was a type and switching/cont parents
    if (curRV.rvType == RVInfo::t_unknown)
      error("Random variable type unknown for variable %s at frame %d, line %d\n",
	    curRV.name.c_str(),curRV.frame,curRV.fileLineNumber);

    // we need to at least have specified one conditional parents set,
    // which might be 'nil' 
    if (curRV.conditionalParents.size() == 0)
      error("Conditional parents unknown/unspecified for random variable '%s' at frame %d, line %d\n",
	    curRV.name.c_str(),curRV.frame,curRV.fileLineNumber);

    // Make sure that if we have switching weights, then we also have a
    // sets of conditional parents, and that we either
    //    1) have only one weight (scale,penalty,shift)
    // or 2) have as many weights as we have sets of conditional parents.
    if (curRV.rvWeightInfo.size() > 1) {
      // then we have switching weights, make sure that we have
      // switching parents and list of conditional parent sets, where
      // the list is the same length.
      if (curRV.conditionalParents.size() != curRV.rvWeightInfo.size()) {
	error("Random variable '%s', frame %d, line %d of file %s has %d switching weights but "
	      "only %d set(s) of conditional parents. Must either have the same number, or a single weight item.\n",
	      curRV.name.c_str(),
	      curRV.frame,
	      curRV.fileLineNumber,
	      curRV.rvFileName.c_str(),
	      curRV.rvWeightInfo.size(),
	      curRV.conditionalParents.size());
      }
    }

    // check that if we have conditionalparents > 1 we also have switching parents.
    if (curRV.conditionalParents.size() > 1 && curRV.switchingParents.size() == 0) {
      error("Random variable '%s', frame %d, line %d of file %s has %d sets of conditional parents but does not list any switching parents.\n",
	    curRV.name.c_str(),
	    curRV.frame,
	    curRV.fileLineNumber,
	    curRV.rvFileName.c_str(),
	    curRV.conditionalParents.size());
    }

    // check if we've already seen this RV
    map < RVInfo::rvParent , unsigned >::iterator it;
    it = nameRVmap.find(RVInfo::rvParent(curRV.name,curRV.frame));
    if (it != nameRVmap.end()) {
      error("Error: random variable '%s' at frame (%d) defined twice, "
	    "on both line %d and line %d\n",curRV.name.c_str(),
	    curRV.frame,
	    curRV.fileLineNumber,
	    rvInfoVector[(*it).second].fileLineNumber);
    }

    curRV.variablePositionInStrFile = rvInfoVector.size();

    // compute the RV's internal status variables
    curRV.computeAndReturnDeterministicStatus();
    curRV.computeAndReturnSparseStatus();

    // everything looks ok, insert it in our tables.
    rvInfoVector.push_back(curRV);
    // add to map
    nameRVmap[
	      RVInfo::rvParent(curRV.name,curRV.frame)
    ]
      = rvInfoVector.size()-1;

  } else {
    // this must be a factor definition.
    curFactor.clear();
    parseFactor();

    // Do dumb linear search since number of factors will probably be
    // very small.
    for (unsigned i=0;i<factorList.size();i++) {
      if (factorList[i].frame == curFactor.frame &&
	  factorList[i].name == curFactor.name) {
	error("Error: factor '%s' at frame (%d) defined twice in same frame, "
	      "on both line %d and line %d\n",
	      curFactor.name.c_str(),
	      curFactor.frame,
	      factorList[i].fileLineNumber,
	      curFactor.fileLineNumber);
      }
    }

    if (curFactor.variables.size() == 0) {
      error("Error: factor '%s' frame (%d), line %d of file %s, must define more than zero variables\n",
	    curFactor.name.c_str(),
	    curFactor.frame,
	    curFactor.fileLineNumber,
	    curFactor.fileName.c_str());
    }

    factorList.push_back(curFactor);
  }

  parseFrameEntryList();
}



void
FileParser::parseFactor()
{

  ensureNotAtEOF(KW_Factor);
  if (tokenInfo != KW_Factor)
    parseError(KW_Factor);
  
  curFactor.frame = curFrame;
  curFactor.fileLineNumber = tokenInfo.srcLine;
  curFactor.fileName  = fileNameParsing;
  consumeToken();


  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseErrorExpecting("':'");
  consumeToken();

  ensureNotAtEOF("factor name");
  if (tokenInfo != TT_Identifier)
    parseErrorExpecting("factor name");
  curFactor.name = tokenInfo.tokenStr;
  consumeToken();

  ensureNotAtEOF("open factor {");
  if (tokenInfo != TT_LeftBrace)
    parseErrorExpecting("open factor {");
  consumeToken();

  parseFactorAttributeList();

  ensureNotAtEOF("close factor }");
  if (tokenInfo != TT_RightBrace)
    parseErrorExpecting("close factor }");
  consumeToken();

}


void
FileParser::parseFactorAttributeList()
{
  if (tokenInfo == KW_Variables  
      || tokenInfo == KW_SymmetricConstraint 
      || tokenInfo == KW_DirectionalConstraint
      || tokenInfo == KW_SoftConstraint) 
    {
      parseFactorAttribute();
      parseFactorAttributeList();
    }
  return;
}


void
FileParser::parseFactorAttribute()
{
  ensureNotAtEOF("factor attribute");
  if (tokenInfo == KW_Variables) {
    return parseFactorVariablesAttribute();
  } else if (tokenInfo == KW_SymmetricConstraint) {
    return parseFactorSymmetricConstraintAttribute();
  } else if (tokenInfo == KW_DirectionalConstraint) {
    return parseFactorDirectionalConstraintAttribute();
  } else if (tokenInfo == KW_SoftConstraint) {
    return parseFactorSoftConstraintAttribute();
  } else
    parseErrorExpecting("factor attribute");

}

void
FileParser::parseFactorVariablesAttribute()
{
  ensureNotAtEOF(KW_Variables);
  if (tokenInfo != KW_Variables)
    parseError(KW_Variables);
  consumeToken();

  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseErrorExpecting("':'");
  consumeToken();

  // use same format as a list of parents to define the list of
  // parents (and use the relative offset notation, where the offset
  // is relative to the current frame number for convenience).
  
  ensureNotAtEOF("list of factor variables");
  rvDeclarationList.clear();
  parseRVDeclarationList();
  curFactor.variables = rvDeclarationList;

  ensureNotAtEOF(";");
  if (tokenInfo != TT_SemiColon)
    parseErrorExpecting("';'");
  consumeToken();

}


void
FileParser::parseFactorSymmetricConstraintAttribute()
{

  /*
      symmetricConstraint: // followed by one of:
	allVarsEqual;
        allVarsUnequal;
        varsNotEqual;
        varsSumTo(n);  // only values that sum to n
        varsMultiplyTo(n); // only values that multiply to n
        varsSumMod(m,n);  // only values such that (sum(vars) % m) = n
	// e.g., force even parity done by sumMod(2,0)
	varsSatisfy using mapping("bla");
	// where mapping is a DT
  */

  ensureNotAtEOF(KW_SymmetricConstraint);
  if (tokenInfo != KW_SymmetricConstraint)
    parseError(KW_SymmetricConstraint);
  consumeToken();

  if (curFactor.fType != FactorInfo::ft_unknown)
    parseError("already defined a type for this factor, can't have more than one,");
  curFactor.fType = FactorInfo::ft_symmetricConstraint;

  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseErrorExpecting("':'");
  consumeToken();

  ensureNotAtEOF("type of symmetric constraint");
  if (tokenInfo == KW_AllVarsEqual) {
    curFactor.symmetricConstraintInfo.symmetricConstraintType = FactorInfo::sct_allVarsEqual;
    consumeToken();

  } else if (tokenInfo == KW_AllVarsUnequal) {
    curFactor.symmetricConstraintInfo.symmetricConstraintType = FactorInfo::sct_allVarsUnequal;
    consumeToken();

  } else if (tokenInfo == KW_VarsNotEqual) {
    curFactor.symmetricConstraintInfo.symmetricConstraintType = FactorInfo::sct_varsNotEqual;
    consumeToken();

  } else if (tokenInfo == KW_VarsSumTo) {
    curFactor.symmetricConstraintInfo.symmetricConstraintType = FactorInfo::sct_varsSumTo;
    consumeToken();

    ensureNotAtEOF("(");
    if (tokenInfo != TT_LeftParen)
      parseErrorExpecting("'('");
    consumeToken();

    ensureNotAtEOF("int to sum");
    if (tokenInfo != TT_Integer)
      parseErrorExpecting("sum-to integer");
    if (tokenInfo.int_val < 0) 
      parseError("sum-to integer must be non-negative");
    curFactor.symmetricConstraintInfo.n = (unsigned) tokenInfo.int_val;
    consumeToken();

    ensureNotAtEOF(")");
    if (tokenInfo != TT_RightParen)
      parseErrorExpecting("')'");
    consumeToken();

  } else if (tokenInfo == KW_VarsMultiplyTo) {
    curFactor.symmetricConstraintInfo.symmetricConstraintType = FactorInfo::sct_varsMultiplyTo;
    consumeToken();

    ensureNotAtEOF("(");
    if (tokenInfo != TT_LeftParen)
      parseErrorExpecting("'('");
    consumeToken();

    ensureNotAtEOF("int to multiply");
    if (tokenInfo != TT_Integer)
      parseErrorExpecting("multiply-to integer");
    if (tokenInfo.int_val < 0) 
      parseError("multiply-to integer must be non-negative");
    curFactor.symmetricConstraintInfo.n = (unsigned) tokenInfo.int_val;
    consumeToken();

    ensureNotAtEOF(")");
    if (tokenInfo != TT_RightParen)
      parseErrorExpecting("')'");
    consumeToken();


  } else if (tokenInfo == KW_VarsSumMod) {
    curFactor.symmetricConstraintInfo.symmetricConstraintType = FactorInfo::sct_varsSumMod;
    consumeToken();


    ensureNotAtEOF("(");
    if (tokenInfo != TT_LeftParen)
      parseErrorExpecting("'('");
    consumeToken();

    ensureNotAtEOF("int to mod");
    if (tokenInfo != TT_Integer)
      parseErrorExpecting("mod-by integer");
    if (tokenInfo.int_val < 0) 
      parseError("mod-by integer must be non-negative");
    curFactor.symmetricConstraintInfo.m = (unsigned) tokenInfo.int_val;
    consumeToken();

    ensureNotAtEOF(",");
    if (tokenInfo != TT_Comma)
      parseErrorExpecting("','");
    consumeToken();

    ensureNotAtEOF("int result of mod");
    if (tokenInfo != TT_Integer)
      parseErrorExpecting("mod-result integer");
    if (tokenInfo.int_val < 0) 
      parseError("mod-result integer must be non-negative");
    curFactor.symmetricConstraintInfo.n = (unsigned) tokenInfo.int_val;
    consumeToken();

    ensureNotAtEOF(")");
    if (tokenInfo != TT_RightParen)
      parseErrorExpecting("')'");
    consumeToken();


  } else if (tokenInfo == KW_VarsSatisfy) {
    curFactor.symmetricConstraintInfo.symmetricConstraintType = FactorInfo::sct_varsSatisfy;
    consumeToken();

    ensureNotAtEOF(KW_Using);
    if (tokenInfo != KW_Using)
      parseError(KW_Using);
    consumeToken();

    parseMappingSpec();
    // name went into listIndex variable, grab it here.
    curFactor.symmetricConstraintInfo.mappingDT = listIndex.nameIndex;
    

  } else 
    parseErrorExpecting("symmetric constraint type");


  ensureNotAtEOF(";");
  if (tokenInfo != TT_SemiColon)
    parseErrorExpecting("';'");
  consumeToken();
}


void
FileParser::parseFactorDirectionalConstraintAttribute()
{
  // directionalConstraint: fooA(0) = functionOf(fooA(-1),fooB(-1)) using mapping("bla");

  if (tokenInfo != KW_DirectionalConstraint)
    parseError(KW_DirectionalConstraint);
  consumeToken();

  if (curFactor.fType != FactorInfo::ft_unknown)
    parseError("already defined a type for this factor");
  curFactor.fType = FactorInfo::ft_directionalConstraint;

  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseErrorExpecting("':'");
  consumeToken();

  ensureNotAtEOF("child variable");
  rvDeclarationList.clear();
  parseRVDeclaration();
  curFactor.directionalConstraintInfo.child = rvDeclarationList[0];

  ensureNotAtEOF("=");
  if (tokenInfo != TT_Equals)
    parseErrorExpecting("'='");
  consumeToken();

  ensureNotAtEOF(KW_FunctionOf);
  if (tokenInfo != KW_FunctionOf)
    parseError(KW_FunctionOf);
  consumeToken();

  ensureNotAtEOF("(");
  if (tokenInfo != TT_LeftParen)
    parseErrorExpecting("'('");
  consumeToken();

  ensureNotAtEOF("list of variables");
  rvDeclarationList.clear();
  parseRVDeclarationList();
  curFactor.directionalConstraintInfo.parents = rvDeclarationList;

  ensureNotAtEOF(")");
  if (tokenInfo != TT_RightParen)
    parseErrorExpecting("')'");
  consumeToken();


  ensureNotAtEOF(KW_Using);
  if (tokenInfo != KW_Using)
    parseError(KW_Using);
  consumeToken();

  parseMappingSpec();
  // name went into listIndex variable, grab it here.
  curFactor.directionalConstraintInfo.mappingDT = listIndex.nameIndex;

  ensureNotAtEOF(";");
  if (tokenInfo != TT_SemiColon)
    parseErrorExpecting("';'");
  consumeToken();
}


void
FileParser::parseFactorSoftConstraintAttribute()
{
  // softConstraint: using table("foo"); // table is a new object, EM learnable, will also
  //                             be used for TableCPT.
  // softConstraint: using logLinear("bar"); // logLinear is a new loglinear object, EM learnable

  if (tokenInfo != KW_SoftConstraint)
    parseError(KW_SoftConstraint);
  consumeToken();

  if (curFactor.fType != FactorInfo::ft_unknown)
    parseError("already defined a type for this factor");
  curFactor.fType = FactorInfo::ft_softConstraint;

  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseErrorExpecting("':'");
  consumeToken();

  ensureNotAtEOF(KW_Using);
  if (tokenInfo != KW_Using)
    parseError(KW_Using);
  consumeToken();

  if (tokenInfo == KW_Table) {
    curFactor.softConstraintInfo.softConstraintType = FactorInfo::fct_table;
  } else if (tokenInfo == KW_LogLinear) {
    curFactor.softConstraintInfo.softConstraintType = FactorInfo::fct_logLinear;
  } else
    parseErrorExpecting("table or logLinear");
  

  ensureNotAtEOF("(");
  if (tokenInfo != TT_LeftParen)
    parseErrorExpecting("'('");
  consumeToken();

  parseListIndex();
  curFactor.softConstraintInfo.name = listIndex.nameIndex;

  ensureNotAtEOF(")");
  if (tokenInfo != TT_RightParen)
    parseErrorExpecting("')'");
  consumeToken();


  ensureNotAtEOF(";");
  if (tokenInfo != TT_SemiColon)
    parseErrorExpecting("';'");
  consumeToken();
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
    parseErrorExpecting("':'");
  consumeToken();

  ensureNotAtEOF("variable name");
  if (tokenInfo != TT_Identifier)
    parseErrorExpecting("variable name");
  curRV.name = tokenInfo.tokenStr;
  consumeToken();

  ensureNotAtEOF("open RV {");
  if (tokenInfo != TT_LeftBrace)
    parseErrorExpecting("open RV {");
  consumeToken();

  parseRandomVariableAttributeList();

  ensureNotAtEOF("close RV }");
  if (tokenInfo != TT_RightBrace)
    parseErrorExpecting("close RV }");
  consumeToken();

}

void
FileParser::parseRandomVariableAttributeList()
{
  if (tokenInfo == KW_Type || 
      tokenInfo == KW_Weight ||
      tokenInfo == KW_EliminationHint ||
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
  else if (tokenInfo == KW_EliminationHint)
    return parseRandomVariableEliminationHintAttribute();
  else if (tokenInfo == KW_Switchingparents ||
	   tokenInfo == KW_Conditionalparents) {
    if (curRV.rvType == RVInfo::t_unknown)
      parseError("type must be first attribute");
    return parseRandomVariableParentAttribute();
  }  else
    parseErrorExpecting("variable attribute");
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
    parseErrorExpecting("':'");
  consumeToken();

  parseRandomVariableType();

  ensureNotAtEOF(";");
  if (tokenInfo != TT_SemiColon)
    parseErrorExpecting("';'");
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
    parseErrorExpecting("variable type");
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
	parseErrorExpecting("first feature range");
      curRV.rvFeatureRange.firstFeatureElement = tokenInfo.int_val;
      consumeToken();

      ensureNotAtEOF("feature range separator");
      if (tokenInfo != TT_Colon)
	parseErrorExpecting("feature range separator");
      consumeToken();
    
      ensureNotAtEOF("second feature range");
      if (tokenInfo != TT_Integer)
	parseErrorExpecting("second feature range");
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

      if (tokenInfo == TT_Integer) {
	if (tokenInfo.int_val < 0)
	  parseError("non-negative value integer");
	curRV.rvFeatureRange.firstFeatureElement = tokenInfo.int_val;
	curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_FirstIsValue;
      } else if (tokenInfo == KW_FrameNum) {
	curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_FrameNumIsValue;
      } else if (tokenInfo == KW_EmarfNum) {
	curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_EmarfNumIsValue;
      } else if (tokenInfo == KW_NumFrames) {
	curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_NumFramesIsValue;
      } else if (tokenInfo == KW_SegmentNum) {
	curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_SegmentNumIsValue;
      } else if (tokenInfo == KW_NumSegments) {
	curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_NumSegmentsIsValue;
      } else {
	parseErrorExpecting("value {integer|frameNum|emarfNum|numFrames|segmentNum|numSegments}");
      }
      // consume whatever the value was.
      consumeToken();

    } else {
      // parse error
      parseErrorExpecting("first feature range n:m | value keyword");
    }
  } else 
    parseErrorExpecting("variable disposition (hidden|discrete)");

  ensureNotAtEOF(KW_Cardinality);
  if (tokenInfo != KW_Cardinality)
    parseError(KW_Cardinality);
  consumeToken();

  ensureNotAtEOF("cardinality value");
  if (tokenInfo != TT_Integer)
    parseErrorExpecting("cardinality value");
  curRV.rvCard = tokenInfo.int_val;
  if (curRV.rvCard <= 1)
    parseError("cardinality must be greater than one (1), for cardinality 1 use an observed variable,");

  if (curRV.rvDisp == RVInfo::d_observed &&
      curRV.rvFeatureRange.filled == RVInfo::FeatureRange::fr_FirstIsValue) {
    // check here that observed value is within cardinality.
    if (curRV.rvFeatureRange.firstFeatureElement >= curRV.rvCard) {
      parseError("given cardinality not large enough for immediate observed value");
    }
  }

  // TODO: allow cardinality 1 variables, and check here for zero
  // cardinality.  cardinality 1 might be reasonable to have
  // "observations" which use a DT that alwyas returns a fixed value
  // (such as frame number). Also, card 1 vars should probably never
  // have any parents.

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
      parseErrorExpecting("first feature range");
    curRV.rvFeatureRange.firstFeatureElement = tokenInfo.int_val;
    consumeToken();

    ensureNotAtEOF("feature range separator");
    if (tokenInfo != TT_Colon)
      parseErrorExpecting("feature range separator");
    consumeToken();
    
    ensureNotAtEOF("second feature range");
    if (tokenInfo != TT_Integer)
      parseErrorExpecting("second feature range");
    curRV.rvFeatureRange.lastFeatureElement = tokenInfo.int_val;
    if (curRV.rvFeatureRange.lastFeatureElement < curRV.rvFeatureRange.firstFeatureElement)
      parseError("first range num must be <= second range num");
    curRV.rvFeatureRange.filled = RVInfo::FeatureRange::fr_Range;
    consumeToken();

  } else
    parseErrorExpecting("variable disposition (hidden|discrete)");

}

void
FileParser::parseRandomVariableEliminationHintAttribute()
{
  ensureNotAtEOF(KW_EliminationHint);
  if (tokenInfo != KW_EliminationHint)
    parseError(KW_EliminationHint);
  consumeToken();  

  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseErrorExpecting("':'");
  consumeToken();

  static const char *const tmp_str = "elimination hint numeric value";
  ensureNotAtEOF(tmp_str);
  // allow an int to be treated as a float value.
  if (tokenInfo != TT_Real && tokenInfo != TT_Integer)
    parseErrorExpecting(tmp_str);
  if (tokenInfo == TT_Real)
    curRV.eliminationOrderHint = tokenInfo.doub_val;
  else 
    curRV.eliminationOrderHint = (double)tokenInfo.int_val;
  consumeToken();

  ensureNotAtEOF(";");
  if (tokenInfo != TT_SemiColon)
    parseErrorExpecting("';'");
  consumeToken();


}


void
FileParser::parseRandomVariableWeightAttribute()
{
  ensureNotAtEOF(KW_Weight);
  if (tokenInfo != KW_Weight)
    parseError(KW_Weight);

  if (curRV.rvWeightInfo.size() > 0 ) {
    parseError("RV already has previously specified weight attribute");
  }
  consumeToken();

  ensureNotAtEOF(":");
  if (tokenInfo != TT_Colon)
    parseErrorExpecting("':'");
  consumeToken();
  parseRandomVariableWeightAttributeSpecList();

  ensureNotAtEOF(";");
  if (tokenInfo != TT_SemiColon) {
    parseErrorExpecting("'}'");
  } else {
    // we've got a semi, so time to end the option list.
    consumeToken();
  }

}


void
FileParser::parseRandomVariableWeightAttributeSpecList()
{
  parseRandomVariableWeightAttributeSpec();
  if (tokenInfo == TT_VirtBar) {
    consumeToken();
    parseRandomVariableWeightAttributeSpecList();
  }
}

void
FileParser::parseRandomVariableWeightAttributeSpec()
{
  // create a new empty weight spec and push
  // it on the list.
  curRV.rvWeightInfo.push_back(RVInfo::WeightInfo());
  ensureNotAtEOF("nil | scale, penalty, or shift");
  if (tokenInfo == KW_Nil) {
    // We've just inserted a weight type that always does nothing.
    consumeToken();
  } else {
    // there is actual weight application here.
    parseRandomVariableWeightOptionList();
  }

}


void
FileParser::parseRandomVariableWeightOptionList()
{
  
  assert ( curRV.rvWeightInfo.size() > 0 );
  const unsigned curWI = curRV.rvWeightInfo.size()-1;

  ensureNotAtEOF("scale, penalty, or shift");

  RVInfo::WeightInfo::WeightItem* curWtItem = NULL;
  if (tokenInfo == KW_Scale) {
    curWtItem = &curRV.rvWeightInfo[curWI].scale;
  } else if (tokenInfo == KW_Penalty) {
    curWtItem = &curRV.rvWeightInfo[curWI].penalty;
  } else if (tokenInfo == KW_Shift) {
    curWtItem = &curRV.rvWeightInfo[curWI].shift;
  } else
    parseErrorExpecting("scale, penalty, or shift");
  consumeToken();
  
  ensureNotAtEOF("integer or floating-point value");
  // we now either have an int which is part of an int:int
  // observation file specification, or we have an int or a float
  // which is an immediate value specification. 
  if (tokenInfo == TT_Real) {
    // we definitely have an immediate value since a real
    // value can never specify a feature position in an observation file.
    curWtItem->wt_Status =RVInfo::WeightInfo::WeightItem::wt_Constant;
    curWtItem->weight_value = tokenInfo.doub_val;;
    consumeToken();
  } else if (tokenInfo == TT_Integer) {
    // we now either have an int which is part of an int:int
    // observation file specification, or we have an int 
    // which is an immediate value specification. 
      
    // store int value for now.
    int firstElement = tokenInfo.int_val;
    consumeToken();

    if (tokenInfo == TT_Colon) {
      // then assume we have a feature range separator
      consumeToken();

      if (tokenInfo != TT_Integer)
	parseErrorExpecting("integer value int:int");	

      // ok, we've got the int:int form.
      curWtItem->wt_Status =RVInfo::WeightInfo::WeightItem::wt_Observation;	

      curWtItem->firstFeatureElement = firstElement;
      curWtItem->lastFeatureElement = tokenInfo.int_val;	

      if ((curWtItem->firstFeatureElement !=
	   curWtItem->lastFeatureElement)
	  || (curWtItem->lastFeatureElement < 0))
	parseError("first range num must be == second range num for observation file weight");
      
      consumeToken();
    } else {
      // then the integer was really the immedate value
      curWtItem->wt_Status =RVInfo::WeightInfo::WeightItem::wt_Constant;
      curWtItem->weight_value  = (double)firstElement;
      // we do not consume current token since it is part of the next lexeme.
    }
  } else
    parseErrorExpecting("integer or floating-point value");
  
  ensureNotAtEOF("; or another weight (scale,penalty,shift) option");
  if (tokenInfo == KW_Scale || tokenInfo == KW_Penalty || tokenInfo == KW_Shift) {
    // it is another weight option.
    parseRandomVariableWeightOptionList();
  } else {
    // we've got something else, time to end.
  }
}


void
FileParser::parseRandomVariableParentAttribute()
{

  ensureNotAtEOF("parent attribute");
  if (tokenInfo == KW_Switchingparents) {
    consumeToken();

    ensureNotAtEOF("attribute seperator");
    if (tokenInfo != TT_Colon)
      parseErrorExpecting("attribute separator");
    consumeToken();

    parseSwitchingParentAttribute();

    ensureNotAtEOF(";");
    if (tokenInfo != TT_SemiColon)
      parseErrorExpecting("';'");
    consumeToken();

  } else if (tokenInfo == KW_Conditionalparents) {
    consumeToken();

    ensureNotAtEOF("attribute seperator");
    if (tokenInfo != TT_Colon)
      parseErrorExpecting("attribute separator");
    consumeToken();

    parseConditionalParentSpecList();

    ensureNotAtEOF(";");
    if (tokenInfo != TT_SemiColon)
      parseErrorExpecting("';'");
    consumeToken();

  } else
    parseErrorExpecting("parent attribute");

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
    parseErrorExpecting("parent RV name");
  p.first = tokenInfo.tokenStr;
  consumeToken();

  ensureNotAtEOF("(");
  if (tokenInfo != TT_LeftParen) 
    parseErrorExpecting("'('");
  consumeToken();

  ensureNotAtEOF("parent RV offset");
  if (tokenInfo != TT_Integer) 
    parseErrorExpecting("parent RV offset");
  p.second = tokenInfo.int_val;;
  // Simple topology check right here, make
  // sure we are not pointing to ourselves.
  if ((p.first == curRV.name) &&
      (p.second == 0)) {
    parseError("parent variable must not refer to self");
  }

  // make sure that the parent wasn't given before.
  for (unsigned pp=0;pp<parentList.size();pp++) {
    if ((p.first == parentList[pp].first)
	&&
	(p.second == parentList[pp].second)) {
      parseError("parent name and frame specified in parent list more than one time");
    }
  }

  parentList.push_back(p);
  consumeToken();

  ensureNotAtEOF(")");
  if (tokenInfo != TT_RightParen) 
    parseErrorExpecting("')'");
  consumeToken();

}



void
FileParser::parseRVDeclarationList()
{
  parseRVDeclaration();
  if (tokenInfo == TT_Comma) {
    consumeToken();
    parseRVDeclarationList();
  }
}

void
FileParser::parseRVDeclaration()
{
  RVInfo::rvParent p;

  ensureNotAtEOF("RV name");
  if (tokenInfo != TT_Identifier) 
    parseErrorExpecting("RV name");
  p.first = tokenInfo.tokenStr;
  consumeToken();

  ensureNotAtEOF("(");
  if (tokenInfo != TT_LeftParen) 
    parseErrorExpecting("'('");
  consumeToken();

  ensureNotAtEOF("RV offset");
  if (tokenInfo != TT_Integer) 
    parseErrorExpecting("RV offset");
  p.second = tokenInfo.int_val;;


  // unlike parents, it is not illegal if a variable in this list
  // occurs twice (although it might be a mistake. Perhaps
  // issue a warning).

  rvDeclarationList.push_back(p);
  consumeToken();

  ensureNotAtEOF(")");
  if (tokenInfo != TT_RightParen) 
    parseErrorExpecting("')'");
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

}


void
FileParser::parseDiscreteImplementation()
{
  ensureNotAtEOF("discrete implementation");
  if (tokenInfo == KW_MDCPT || tokenInfo == KW_MSCPT
      || tokenInfo == KW_MTCPT || tokenInfo == KW_NGRAMCPT || tokenInfo == KW_FNGRAMCPT
      || tokenInfo == KW_VECPT || tokenInfo == KW_LATTICENODECPT || tokenInfo == KW_LATTICEEDGECPT ) {

    if (tokenInfo == KW_MDCPT)
      curRV.discImplementations.push_back(CPT::di_MDCPT);
    else if (tokenInfo == KW_MSCPT)
      curRV.discImplementations.push_back(CPT::di_MSCPT);
    else if (tokenInfo == KW_MTCPT)
      curRV.discImplementations.push_back(CPT::di_MTCPT);
    else if (tokenInfo == KW_NGRAMCPT)
      curRV.discImplementations.push_back(CPT::di_NGramCPT);
    else if (tokenInfo == KW_FNGRAMCPT)
      curRV.discImplementations.push_back(CPT::di_FNGramCPT);
    else if (tokenInfo == KW_VECPT)
      curRV.discImplementations.push_back(CPT::di_VECPT);
    else if (tokenInfo == KW_LATTICENODECPT)
	    curRV.discImplementations.push_back(CPT::di_LatticeNodeCPT);
    else if (tokenInfo == KW_LATTICEEDGECPT)
	    curRV.discImplementations.push_back(CPT::di_LatticeEdgeCPT);

    consumeToken();


    ensureNotAtEOF("(");
    if (tokenInfo != TT_LeftParen) {
      parseErrorExpecting("'('");
    }
    consumeToken();

    parseListIndex();

    ensureNotAtEOF(") or ,");
    if (tokenInfo != TT_RightParen) {
      if ( tokenInfo != TT_Comma )
	parseErrorExpecting("')' or ','");
      // parse CPT with fewer parents like
      // using FNGramCPT("fngram", 0, 1);
      
      while ( tokenInfo != TT_RightParen ) {
	if ( tokenInfo != TT_Comma )
	  parseErrorExpecting("','");
	consumeToken();
	
	if ( tokenInfo != TT_Integer )
	  parseErrorExpecting("integer");
	listIndex.fnparents.push_back(tokenInfo.int_val);
	consumeToken();
      }
      
      // before push_back curRV.listIndices.size() is the last index which we are looking for
      if ( listIndex.fnparents.size() != curRV.conditionalParents[curRV.listIndices.size()].size() ) {
	error("Error: RV \"%s\" at frame %d (line %d), conditional parent set %d has %d parents but FNGramCPT sub-parent list specifices %d parents indices. They must match\n",
	      curRV.name.c_str(),
	      curRV.frame,
	      curRV.listIndices.size(),
	      curRV.conditionalParents[curRV.listIndices.size()].size(),
	      listIndex.fnparents.size());
      }
    }

    // consume what now must be the right paren.
    consumeToken();

    curRV.listIndices.push_back(listIndex);

    // we need to clear this up.
    listIndex.fnparents.resize(0);
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
      parseErrorExpecting("')'");
    }
    consumeToken();

    //////////////////////////////////////
    // AAA: some semantic checking here (see AAA: tag below)
    // The current implementation presumably comes from a RV
    // with nil conditional parents, as in:
    //
    // ... | nil using mixture("the_forth_gaussian_mixture") | ...
    //
    // this means that we are specifying one and only one particular
    // gaussian mixture "the_forth_gaussian_mixture" in the gaussian file,
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
      parseErrorExpecting("'('");
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
      parseErrorExpecting("collection name");

    ensureNotAtEOF(")");
    if (tokenInfo != TT_RightParen) {
      parseErrorExpecting("')'");
    }
    consumeToken();

    // parse a string of the form 'mapping("bar")'
    parseMappingSpec();

    listIndex.collectionName = collectionName;

    // save the data.
    curRV.listIndices.push_back(listIndex);

    //////////////////////////////////////////////////////////////////
    // semantic check, this is the dual check of the
    // check at "AAA:" above.
    if (curRV.conditionalParents[curRV.conditionalParents.size()-1].size() == 0)
      parseError("decision tree 'mapping' needs > 0 conditional parents");
  }
}


void
FileParser::parseContObsDistType()
{
  ensureNotAtEOF("continuous observation distribution type");
  if (tokenInfo == KW_Mixture) {
    curRV.contImplementations.push_back(MixtureCommon::ci_mixture);
    consumeToken();
  } else if (tokenInfo == KW_GausSwitchMixture) {
    curRV.contImplementations.push_back(MixtureCommon::ci_gausSwitchMixture);
    consumeToken();
  }  else if (tokenInfo == KW_LogitSwitchMixture) {
    curRV.contImplementations.push_back(MixtureCommon::ci_logitSwitchMixture);
    consumeToken();
  }  else if (tokenInfo == KW_MlpSwitchMixture) {
    curRV.contImplementations.push_back(MixtureCommon::ci_mlpSwitchMixture);
    consumeToken();
  } else
    parseErrorExpecting("continuous observation distribution type");
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
    parseErrorExpecting("first chunk integer");
  if (tokenInfo.int_val < 0) 
    parseError("non-negative chunk range");
  _firstChunkframe = (unsigned)tokenInfo.int_val;
  consumeToken();

  ensureNotAtEOF("chunk colon");
  if (tokenInfo != TT_Colon) 
    parseErrorExpecting("chunk colon");
  consumeToken();

  ensureNotAtEOF("second chunk integer");
  if (tokenInfo != TT_Integer) 
    parseErrorExpecting("second chunk integer");
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
    parseErrorExpecting("'('");
  }
  consumeToken();

  parseListIndex();

  ensureNotAtEOF(")");
  if (tokenInfo != TT_RightParen) {
    parseErrorExpecting("')'");
  }
  consumeToken();

}

void
FileParser::parseListIndex()
{
  ensureNotAtEOF("name of object");
  if (tokenInfo == TT_Integer) {
    // TODO: need to remove the integer index code.
    error("Integers no longer supported as object specifiers\n");
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
    parseErrorExpecting("name of object");
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

    RV*rv;
    // go through and consider all possible RV types instantiating the
    // currect one here.
    if (rvInfoVector[i].rvType == RVInfo::t_discrete) {
      // discrete
      if (rvInfoVector[i].rvDisp == RVInfo::d_hidden) {
	// discrete hidden
	if (rvInfoVector[i].switchingParents.size() > 0) {
	  // discrete hidden switching
	  if (rvInfoVector[i].rvWeightInfo.size() > 0) {
	    // discrete hidden switching weighted
	    rv = new ScPnSh_Sw_HidDiscRV(rvInfoVector[i],rvInfoVector[i].frame,rvInfoVector[i].rvCard);
	  } else {
	    // discrete hidden switching not-weighted
	    rv = new Sw_HidDiscRV(rvInfoVector[i],rvInfoVector[i].frame,rvInfoVector[i].rvCard);
	  }
	} else {
	  // discrete hidden no-switching
	  if (rvInfoVector[i].rvWeightInfo.size() > 0) {
	    // discrete hidden no-switching weighted
	    rv = new ScPnSh_HidDiscRV(rvInfoVector[i],rvInfoVector[i].frame,rvInfoVector[i].rvCard);
	  } else {
	    // discrete hidden no-switching not-weighted
	    rv = new HidDiscRV(rvInfoVector[i],rvInfoVector[i].frame,rvInfoVector[i].rvCard);
	  }
	}
      } else {
	// discrete observed
	if (rvInfoVector[i].switchingParents.size() > 0) {
	  // discrete observed switching
	  if (rvInfoVector[i].rvWeightInfo.size() > 0) {
	    // discrete observed switching weighted
	    rv = new ScPnSh_Sw_ObsDiscRV(rvInfoVector[i],rvInfoVector[i].frame,rvInfoVector[i].rvCard);
	  } else {
	    // discrete observed switching not-weighted
	    rv = new Sw_ObsDiscRV(rvInfoVector[i],rvInfoVector[i].frame,rvInfoVector[i].rvCard);
	  }
	} else {
	  // discrete observed no-switching
	  if (rvInfoVector[i].rvWeightInfo.size() > 0) {
	    // discrete observed no-switching weighted
	    rv = new ScPnSh_ObsDiscRV(rvInfoVector[i],rvInfoVector[i].frame,rvInfoVector[i].rvCard);
	  } else {
	    // discrete observed no-switching not-weighted
	    rv = new ObsDiscRV(rvInfoVector[i],rvInfoVector[i].frame,rvInfoVector[i].rvCard);
	  }
	}
      }
    } else {
      // continuous
      if (rvInfoVector[i].rvDisp == RVInfo::d_hidden) {
	// continuous hidden
	error("ERROR: GMTK does not yet support hidden continuous RVs.");
	rv = NULL; // suppress compiler warning.
      } else {
	// continuous observed
	if (rvInfoVector[i].switchingParents.size() > 0) {
	  // continuous observed switching
	  if (rvInfoVector[i].rvWeightInfo.size() > 0) {
	    // continuous observed switching weighted
	    rv = new ScPnSh_Sw_ObsContRV(rvInfoVector[i],rvInfoVector[i].frame);
	  } else {
	    // continuous observed switching not-weighted
	    rv = new Sw_ObsContRV(rvInfoVector[i],rvInfoVector[i].frame);
	  }
	} else {
	  // continuous observed no-switching
	  if (rvInfoVector[i].rvWeightInfo.size() > 0) {
	    // continuous observed no-switching weighted
	    rv = new ScPnSh_ObsContRV(rvInfoVector[i],rvInfoVector[i].frame);
	  } else {
	    // continuous observed no-switching not-weighted
	    rv = new ObsContRV(rvInfoVector[i],rvInfoVector[i].frame);
	  }
	}
      }
    }
    rvInfoVector[i].rv = rv;    
    // rv->printSelf(stdout,true);
  }
  // now set up all the parents of each random variable.
  for (unsigned i=0;i<rvInfoVector.size();i++) {

    const unsigned frame = rvInfoVector[i].frame;

    // build the (possibly empty) switching parents list.
    vector<RV *> sparents;
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
    vector<vector<RV * > > cpl(rvInfoVector[i].conditionalParents.size());
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


  // now go through factor list and added edges to neighbors of each
  // random variable (and make sure that all variables in each
  // factor exists according to the template).

  for (unsigned factorNo=0;factorNo<factorList.size();factorNo++) {
    FactorInfo& factor = factorList[factorNo];
    factor.rvs.resize(factor.variables.size());

    set < RV* > factorVarsSet; // possibly include this in the FactorInfo itself

    for (unsigned varNo=0;varNo<factor.variables.size();varNo++) {
      RVInfo::rvParent pp(factor.variables[varNo].first,
			  factor.variables[varNo].second 
			  + factor.frame);

      if (nameRVmap.find(pp) == nameRVmap.end()) {
	error("Error: RV \"%s(%d)\" doesn't exist, declared as \"%s(%d)\" in factor '%s' at frame %d (line %d of file '%s')\n",
	      pp.first.c_str(),pp.second,
	      pp.first.c_str(),factor.variables[varNo].second,
	      factor.name.c_str(),
	      factor.frame,
	      factor.fileLineNumber,
	      factor.fileName.c_str());
      } else {
	factor.rvs[varNo] = rvInfoVector[ nameRVmap[ pp ] ].rv;
	// factors only defined over discrete variables for now.
	if (factor.rvs[varNo]->continuous()) {
	  error("Error: RV \"%s(%d)\" declared as \"%s(%d)\" in factor '%s' at frame %d (line %d of file '%s'), must be discrete (for now)\n",
		pp.first.c_str(),pp.second,
		pp.first.c_str(),factor.variables[varNo].second,
		factor.name.c_str(),
		factor.frame,
		factor.fileLineNumber,
		factor.fileName.c_str());
	}
	factorVarsSet.insert(factor.rvs[varNo]);
      }
    }

    // do a bit more sanity checking of factors.
    if (factor.fType == FactorInfo::ft_directionalConstraint) {
      // check that each of child and parent are defined in the factor.
      
      // check child
      RVInfo::rvParent pp(factor.directionalConstraintInfo.child.first,
			  factor.directionalConstraintInfo.child.second 
			  + factor.frame);
      if (nameRVmap.find(pp) == nameRVmap.end())
	error("Error: RV \"%s(%d)\" doesn't exist, declared as dependent variable \"%s(%d)\" in factor '%s' with directional constraint at frame %d (line %d of file '%s')\n",
	      pp.first.c_str(),pp.second,
	      pp.first.c_str(),factor.directionalConstraintInfo.child.second,
	      factor.name.c_str(),
	      factor.frame,
	      factor.fileLineNumber,
	      factor.fileName.c_str());

      RV* child = rvInfoVector[ nameRVmap[ pp ] ].rv;
      if (factorVarsSet.find(child) == factorVarsSet.end()) {
	error("Error: RV \"%s(%d)\" must exist in current factor, declared as dependent variable \"%s(%d)\" in factor '%s' with directional constraint at frame %d (line %d of file '%s')\n",
	      pp.first.c_str(),pp.second,
	      pp.first.c_str(),factor.directionalConstraintInfo.child.second,
	      factor.name.c_str(),
	      factor.frame,
	      factor.fileLineNumber,
	      factor.fileName.c_str());
      }

      // now go through and do the same deal for the parents.
      for (unsigned varNo=0;varNo<factor.directionalConstraintInfo.parents.size();varNo++) {      
	RVInfo::rvParent pp(factor.directionalConstraintInfo.parents[varNo].first,
			    factor.directionalConstraintInfo.parents[varNo].second 
			    + factor.frame);

	if (nameRVmap.find(pp) == nameRVmap.end())
	  error("Error: RV \"%s(%d)\" doesn't exist, declared as independent variable \"%s(%d)\" in factor '%s' with directional constraint at frame %d (line %d of file '%s')\n",
		pp.first.c_str(),pp.second,
		pp.first.c_str(),factor.directionalConstraintInfo.child.second,
		factor.name.c_str(),
		factor.frame,
		factor.fileLineNumber,
		factor.fileName.c_str());

	RV* child = rvInfoVector[ nameRVmap[ pp ] ].rv;
	if (factorVarsSet.find(child) == factorVarsSet.end()) {
	  error("Error: RV \"%s(%d)\" must exist in current factor, declared as independent variable \"%s(%d)\" in factor '%s' with directional constraint at frame %d (line %d of file '%s')\n",
		pp.first.c_str(),pp.second,
		pp.first.c_str(),factor.directionalConstraintInfo.child.second,
		factor.name.c_str(),
		factor.frame,
		factor.fileLineNumber,
		factor.fileName.c_str());
	}

      }
      
    }

    // possible other checks of factors here.

    // We add neighbor members for each template random variable
    // according to factors at the point just before we need to
    // triangulate.

  }

}



/*-
 *-----------------------------------------------------------------------
 * ensureValidTemplate()
 *      Takes the current graph instantiation, and makes sure
 *      that it is valid (i.e., directed, a valid template, etc.)
 *      die with an error if not.
 *
 * Preconditions:
 *      parseGraphicalModel() must have been called.
 *      createRandomVariableGraph() must have been called.
 *
 * Postconditions:
 *      It is true that graph is valid.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      nothing
 *
 *-----------------------------------------------------------------------
 */
void
FileParser::ensureValidTemplate(bool longCheck)
{
  vector <RV*> vars;
  vector <RV*> vars2;

  const unsigned longCheckStart = 3;

  // always do short check
  for (unsigned unrollAmount=0;unrollAmount<longCheckStart;unrollAmount++) {
    infoMsg(Max,"Ensuring Valid Template when unrolling %d out of %d times\n",unrollAmount,longCheckStart-1);
    unroll(unrollAmount,vars);
    // TODO: fix error messages to give indication as to where loop is.
    if (!GraphicalModel::topologicalSort(vars,vars2))
      error("ERROR. Graph is not acyclic, contains a directed loop when unrolled %d times.\n",unrollAmount);
    // need to delete the variables here.
    for (unsigned i=0;i<vars.size();i++)
      delete vars[i];
    vars.clear(); vars2.clear();
  }

  if (longCheck) {
    // note that it is possible for the directed loop to show up only when the graph is unrolled numVarsInChunk
    // times. Unrolling this amount is sufficient for any further unrolling though, but this
    // can take a long time for big graphs.
    for (unsigned unrollAmount=longCheckStart;unrollAmount<=numVarsInChunk;unrollAmount++) 
      {
	infoMsg(Max,"Longcheck: Ensuring Valid Template when unrolling %d out of %d times\n",unrollAmount,numVarsInChunk);
	unroll(unrollAmount,vars);
	// TODO: fix error messages to give indication as to where loop is.
	if (!GraphicalModel::topologicalSort(vars,vars2))
	  error("ERROR. Graph is not directed, contains a directed loop when unrolled %d times.\n",unrollAmount);
	// need to delete the variables here.
	for (unsigned i=0;i<vars.size();i++)
	  delete vars[i];
	vars.clear(); vars2.clear();
      }
  }
  infoMsg(Max,"Done ensuring valid template, longCheck = %d\n",longCheck);

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
      // do nothing, since there is no switching going on.
    } else {

      RngDecisionTree *dtMapper;
      if (rvInfoVector[i].switchMapping.liType == 
	  RVInfo::ListIndex::li_String) {
	if (GM_Parms.dtsMap.find(rvInfoVector[i].switchMapping.nameIndex) == GM_Parms.dtsMap.end())
	  error("Error: RV \"%s\" at frame %d (line %d), switching parent DT \"%s\" doesn't exist\n",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		rvInfoVector[i].switchMapping.nameIndex.c_str());

	dtMapper = 
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
	dtMapper = GM_Parms.dts[rvInfoVector[i].switchMapping.intIndex];
      } else {
	// this shouldn't happen, unless the parser has a bug.
	assert ( 0 );
	dtMapper = NULL; // to avoid compiler warning.
      }

      // now check to make sure that the decision tree matches
      // the set of parents that were set up as the switching parents.
      if (dtMapper->numFeatures() != rvInfoVector[i].switchingParents.size()) {
	error("Error: RV \"%s\" at frame %d (line %d), num switching parents different than required by decision tree named \"%s\".\n",
	      rvInfoVector[i].name.c_str(),
	      rvInfoVector[i].frame,
	      rvInfoVector[i].fileLineNumber,
	      dtMapper->name().c_str());
      }
      // everything ok, so we set it. If this is not a RV that accepts
      // a DT, this will cause a dynamic error.
      rvInfoVector[i].rv->setDTMapper(dtMapper);
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
      DiscRV* rv = 
	(DiscRV*) rvInfoVector[i].rv;

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
		error("Error: RV \"%s\" at frame %d (line %d), conditional parent DenseCPT \"%s\" doesn't exist\n",
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
	    // TODO: need to remove the integer index code.
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
		error("Error: RV \"%s\" at frame %d (line %d), conditional parent SparseCPT \"%s\" doesn't exist\n",
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
	      // TODO: need to remove the integer index code.
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
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent DeterministicCPT \"%s\" doesn't exist\n",
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
	      // TODO: need to remove the integer index code.
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

	} else
	  if (rvInfoVector[i].discImplementations[j] == CPT::di_NGramCPT) {

	    /////////////////////////////////////////////////////////
	    // Once again, same code as above, but using NGramCPTs rather
	    // then MDCPTs, MSCPTs or MTCPTs.

	    //////////////////////////////////////////////////////
	    // set the CPT to a NGramCPT, depending on if a string
	    // or integer index was used in the file.
	    if (rvInfoVector[i].listIndices[j].liType
		== RVInfo::ListIndex::li_String) {
	      if (GM_Parms.ngramCptsMap.find(
					  rvInfoVector[i].listIndices[j].nameIndex) ==
		  GM_Parms.ngramCptsMap.end()) {
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent NGramCPT \"%s\" doesn't exist\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].nameIndex.c_str());
	      } else {
		// otherwise add it
		cpts[j] = (CPT*)
		  GM_Parms.ngramCpts[
				  GM_Parms.ngramCptsMap[
						     rvInfoVector[i].listIndices[j].nameIndex
				  ]
		  ];
	      }
	    } else {
	      // TODO: need to remove the integer index code.
	      assert(0);
#if 0
	      if (rvInfoVector[i].listIndices[j].intIndex >=
		  GM_Parms.ngramCpts.size()) {
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
		  GM_Parms.ngramCpts[rvInfoVector[i].listIndices[j].intIndex];
	      }
#endif
	    }

	} else
	  if (rvInfoVector[i].discImplementations[j] == CPT::di_FNGramCPT) {
	    /////////////////////////////////////////////////////////
	    // Once again, same code as above, but using FNGramCPTs rather
	    // then MDCPTs, MSCPTs or MTCPTs.

	    //////////////////////////////////////////////////////
	    // set the CPT to a FNGramCPT, depending on if a string
	    // or integer index was used in the file.
		if (rvInfoVector[i].listIndices[j].liType == RVInfo::ListIndex::li_String) {
			if ( GM_Parms.fngramImpsMap.find(rvInfoVector[i].listIndices[j].nameIndex) == GM_Parms.fngramImpsMap.end() )
				error("Error: RV \"%s\" at frame %d (line %d), conditional parent FNGramCPT \"%s\" doesn't exist\n",
						rvInfoVector[i].name.c_str(), rvInfoVector[i].frame, rvInfoVector[i].fileLineNumber, rvInfoVector[i].listIndices[j].nameIndex.c_str());

			// the naming convention is if there is full parents, then use fngramCPT name
			// otherwise use name like "fngram:0,1,3"
			string fngramCptName = rvInfoVector[i].listIndices[j].nameIndex;
			if ( rvInfoVector[i].listIndices[j].fnparents.size() != 0 ) {
				fngramCptName += ":";
				for ( unsigned k = 0; k < rvInfoVector[i].listIndices[j].fnparents.size() - 1; k++ ) {
					char tmp[20];
					sprintf(tmp, "%d,", rvInfoVector[i].listIndices[j].fnparents[k]);
					fngramCptName += string(tmp);
				}
				char tmp[20];
				sprintf(tmp, "%d", rvInfoVector[i].listIndices[j].fnparents[rvInfoVector[i].listIndices[j].fnparents.size() - 1]);
				fngramCptName += string(tmp);
			}

			if ( GM_Parms.fngramCptsMap.find(fngramCptName) == GM_Parms.fngramCptsMap.end() ) {
				// Here we will contruct the object for FNGramCPT based on FNGramImp
				// we will check whether it is "ftrigram" or "ftrigram:0,1"
				if ( rvInfoVector[i].listIndices[j].fnparents.size() == 0 ) {
					// create new FNGramCPT
					FNGramCPT *ob = new FNGramCPT();
					ob->setNumParents(rvInfoVector[i].conditionalParents[j].size());
					ob->setFNGramImp(GM_Parms.fngramImps[GM_Parms.fngramImpsMap[rvInfoVector[i].listIndices[j].nameIndex]]);
					ob->setName(rvInfoVector[i].listIndices[j].nameIndex);

					// add it to the map list
					GM_Parms.fngramCptsMap[fngramCptName] = GM_Parms.fngramCpts.size();
					GM_Parms.fngramCpts.push_back(ob);
					cpts[j] = (CPT*)ob;
				} else {
					// create new FNGramCPT
					FNGramCPT *ob = new FNGramCPT();
					ob->setNumParents(rvInfoVector[i].conditionalParents[j].size());
					ob->setFNGramImp(GM_Parms.fngramImps[GM_Parms.fngramImpsMap[rvInfoVector[i].listIndices[j].nameIndex]]);

					// set the parent index
					ob->setName(fngramCptName);
					ob->setParentsPositions(rvInfoVector[i].listIndices[j].fnparents);

					// add it into the map list
					GM_Parms.fngramCptsMap[fngramCptName] = GM_Parms.fngramCpts.size();
					GM_Parms.fngramCpts.push_back(ob);
					cpts[j] = (CPT*)ob;
				}
			} else {
				// otherwise add it
				cpts[j] = (CPT*)GM_Parms.fngramCpts[GM_Parms.fngramCptsMap[rvInfoVector[i].listIndices[j].nameIndex]];
			}

	    } else {
	      // TODO: need to remove the integer index code.
	      assert(0);
#if 0
	      if (rvInfoVector[i].listIndices[j].intIndex >=
		  GM_Parms.fngramCpts.size()) {
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
		  GM_Parms.fngramCpts[rvInfoVector[i].listIndices[j].intIndex];
	      }
#endif
	    }

	} else
	  if (rvInfoVector[i].discImplementations[j] == CPT::di_VECPT) {

	    /////////////////////////////////////////////////////////
	    // Once again, same code as above, but using VECPTs.

	    //////////////////////////////////////////////////////
	    // set the CPT to a VECPT, depending on if a string
	    // or integer index was used in the file.
	    if (rvInfoVector[i].listIndices[j].liType
		== RVInfo::ListIndex::li_String) {
	      if (GM_Parms.veCptsMap.find(
					  rvInfoVector[i].listIndices[j].nameIndex) ==
		  GM_Parms.veCptsMap.end()) {
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent VirtualEvidenceCPT \"%s\" doesn't exist\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].nameIndex.c_str());
	      } else {
		// otherwise add it
		cpts[j] = 
		  GM_Parms.veCpts[
				  GM_Parms.veCptsMap[
						     rvInfoVector[i].listIndices[j].nameIndex
				  ]
		  ];
	      }
	    } else {
	      // TODO: need to remove the integer index code.
	      assert(0);
	    }
	  } else if (rvInfoVector[i].discImplementations[j] == CPT::di_LatticeNodeCPT) {

	    /////////////////////////////////////////////////////////
	    // Once again, same code as above, but using LatticeNodeCPTs rather
	    // then MDCPTs or MSCPTs. 

	    //////////////////////////////////////////////////////
	    // set the CPT to a LatticeNodeCPT, depending on if a string
	    // or integer index was used in the file.
	    if (rvInfoVector[i].listIndices[j].liType == RVInfo::ListIndex::li_String) {
	      if (GM_Parms.latticeNodeCptsMap.find(
					  rvInfoVector[i].listIndices[j].nameIndex) ==
		  GM_Parms.latticeNodeCptsMap.end()) {
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent LatticeNodeCPT \"%s\" doesn't exist\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].nameIndex.c_str());
	      } else {
		// otherwise add it
		cpts[j] = 
		  GM_Parms.latticeNodeCpts[
				  GM_Parms.latticeNodeCptsMap[
						     rvInfoVector[i].listIndices[j].nameIndex
				  ]
		  ];
	      }
	    } else {
	      // TODO: need to remove the integer index code.
	      assert(0);
	    }
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_LatticeEdgeCPT) {

	    /////////////////////////////////////////////////////////
	    // Once again, same code as above, but using LatticeEdgeCPTs rather
	    // then MDCPTs or MSCPTs. 

	    //////////////////////////////////////////////////////
	    // set the CPT to a LatticeEdgeCPT, depending on if a string
	    // or integer index was used in the file.
	    if (rvInfoVector[i].listIndices[j].liType == RVInfo::ListIndex::li_String) {
	      if (GM_Parms.latticeEdgeCptsMap.find(
					  rvInfoVector[i].listIndices[j].nameIndex) ==
		  GM_Parms.latticeEdgeCptsMap.end()) {
		  error("Error: RV \"%s\" at frame %d (line %d), conditional parent LatticeNodeCPT \"%s\" doesn't exist\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].nameIndex.c_str());
	      } else {
		// otherwise add it
		cpts[j] = 
		  GM_Parms.latticeEdgeCpts[
				  GM_Parms.latticeEdgeCptsMap[
						     rvInfoVector[i].listIndices[j].nameIndex
				  ]
		  ];
	      }
	    } else {
	      // TODO: need to remove the integer index code.
	      assert(0);
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
	  cptType = "DenseCPT";
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_MSCPT) {
	  cptType = "SparseCPT";
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_MTCPT) {
	  cptType = "DeterministicCPT";
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_NGramCPT) {
		cptType = "NGramCPT";
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_FNGramCPT) {
		cptType = "FNGramCPT";
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_VECPT) {
		cptType = "VirtualEvidenceCPT";
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_LatticeNodeCPT) {
		cptType = "LatticeNodeCPT";
	} else if (rvInfoVector[i].discImplementations[j] == CPT::di_LatticeEdgeCPT) {
		cptType = "LatticeEdgeCPT";
	}

	// check to make sure this cpt matches this
	// number of parents.
	// Special case for NGramCPT because we only restrict number of parents is less than ngram order.
	if ( rvInfoVector[i].discImplementations[j] == CPT::di_NGramCPT ) {
		if ( cpts[j]->numParents() < rvInfoVector[i].conditionalParents[j].size()) {
			error("Error: RV \"%s\" at frame %d (line %d), number of parents (at switching condition %d) is %d, but that is different than what is required by %s \"%s\" which is %d.\n",
			      rvInfoVector[i].name.c_str(), 
			      rvInfoVector[i].frame, 
			      rvInfoVector[i].fileLineNumber, 
			      j, 
			      rvInfoVector[i].conditionalParents[j].size(),
			      cptType.c_str(), 
			      cpts[j]->name().c_str(),
			      cpts[j]->numParents()
			      );
		}
	} else if (cpts[j]->numParents() !=
	    rvInfoVector[i].conditionalParents[j].size()) {
	  error("Error: RV \"%s\" at frame %d (line %d), number of parents (at switching condition %d) is %d, but that is different than required by %s \"%s\" which is %d.\n",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		j,
		rvInfoVector[i].conditionalParents[j].size(),
		cptType.c_str(),
		cpts[j]->name().c_str(),
		cpts[j]->numParents()
		);
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
	} else if ( cpts[j]->cptType == CPT::di_NGramCPT ) {
	  // Because we allow fewer parents in using ngram cpt, we use cpts[j]->numParents() instead.
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
	  for ( unsigned par = 0; par < rv->condParentsVec(j).size(); par++ ) {
	    if ( RV2DRV(rv->condParentsVec(j)[par])->cardinality != cpts[j]->parentCardinality(par) )
	      error("Error: RV \"%s\" at frame %d (line %d), cardinality of parent '%s' is %d, but %d'th parent of %s \"%s\" requires cardinality of %d.\n",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rv->condParentsVec(j)[par]->name().c_str(),
		    RV2DRV(rv->condParentsVec(j)[par])->cardinality,
		    par,
		    cptType.c_str(),
		    cpts[j]->name().c_str(),
		    cpts[j]->parentCardinality(par));
	  }
	} else if ( cpts[j]->cptType == CPT::di_LatticeNodeCPT ) {
	  // lattice node cpt only has three or two parent
	  if ( ((LatticeNodeCPT*)cpts[j])->useTimeParent() ) {
	    if ( cpts[j]->numParents() != 3 )
	    error("Error: RV \"%s\" at frame %d (line %d), should have only three parents",
		  rvInfoVector[i].name.c_str(),
		  rvInfoVector[i].frame,
		  rvInfoVector[i].fileLineNumber);
	  } else if ( cpts[j]->numParents() != 2 )
	    error("Error: RV \"%s\" at frame %d (line %d), should have only two parents",
		  rvInfoVector[i].name.c_str(),
		  rvInfoVector[i].frame,
		  rvInfoVector[i].fileLineNumber);
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
	    if (RV2DRV(rv->condParentsVec(j)[par])->cardinality !=
		cpts[j]->parentCardinality(par))
	      error("Error: RV \"%s\" at frame %d (line %d), cardinality of parent '%s' is %d, but %d'th parent of %s \"%s\" requires cardinality of %d.\n",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rv->condParentsVec(j)[par]->name().c_str(),
		    RV2DRV(rv->condParentsVec(j)[par])->cardinality,
		    par,
		    cptType.c_str(),
		    cpts[j]->name().c_str(),
		    cpts[j]->parentCardinality(par));
	  }
	} else if ( cpts[j]->cptType == CPT::di_LatticeEdgeCPT ) {
		// lattice edge cpt only has two parents
		if ( cpts[j]->numParents() != 2 )
			error("Error: RV \"%s\" at frame %d (line %d), should have only one parent",
				rvInfoVector[i].name.c_str(),
				rvInfoVector[i].frame,
				rvInfoVector[i].fileLineNumber);
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
	    if (RV2DRV(rv->condParentsVec(j)[par])->cardinality !=
		cpts[j]->parentCardinality(par))
	      error("Error: RV \"%s\" at frame %d (line %d), cardinality of parent '%s' is %d, but %d'th parent of %s \"%s\" requires cardinality of %d.\n",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rv->condParentsVec(j)[par]->name().c_str(),
		    RV2DRV(rv->condParentsVec(j)[par])->cardinality,
		    par,
		    cptType.c_str(),
		    cpts[j]->name().c_str(),
		    cpts[j]->parentCardinality(par));
	  }
	} else {
	  // regular general case checking below.
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
	    if (RV2DRV(rv->condParentsVec(j)[par])->cardinality !=
		cpts[j]->parentCardinality(par))
	      error("Error: RV \"%s\" at frame %d (line %d), cardinality of parent '%s' is %d, but %d'th parent of %s \"%s\" requires cardinality of %d.\n",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rv->condParentsVec(j)[par]->name().c_str(),
		    RV2DRV(rv->condParentsVec(j)[par])->cardinality,
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
      ContRV* rv =
	(ContRV*) rvInfoVector[i].rv;

      rv->conditionalMixtures.resize
	(
	 rvInfoVector[i].conditionalParents.size()
	 );
      // set the current one to position 0. This will for for
      // non-switching variables, and switching variables will modify
      // it as appropriate.
      rv->curMappingOrDirect = &rv->conditionalMixtures[0];
      for (unsigned j=0;j<rvInfoVector[i].conditionalParents.size();j++) {

	// TODO:: implement the other continuous implementations and
	// allow the other cases here.
	if (rvInfoVector[i].contImplementations[j] !=
	    MixtureCommon::ci_mixture) {
	  error("ERROR: Only supports mixture implementations for now");
	  // ultimately, we can remove this and we'll have to duplicate
	  // the below code for each continuous implementation, just
	  // like what was done above for the discrete implementations.
	}

	if (rvInfoVector[i].conditionalParents[j].size() == 0) {
	  // then under this value for switching parents, there
	  // are no conditional parents and we have a direct
	  // pointer to some mixture.
	  rv->conditionalMixtures[j].direct = true;
	  if (rvInfoVector[i].listIndices[j].liType 
		== RVInfo::ListIndex::li_String) {
	    if (GM_Parms.mixturesMap.find(
		    rvInfoVector[i].listIndices[j].nameIndex) ==
		GM_Parms.mixturesMap.end()) {
	      error("Error: RV \"%s\" at frame %d (line %d), mixture \"%s\" doesn't exist\n",
		      rvInfoVector[i].name.c_str(),
		      rvInfoVector[i].frame,
		      rvInfoVector[i].fileLineNumber,
		      rvInfoVector[i].listIndices[j].nameIndex.c_str());
	    } else {
	      // Mixture is there, so add it.
	      rv->conditionalMixtures[j].mixture =
		GM_Parms.mixtures[
		      GM_Parms.mixturesMap[
			  rvInfoVector[i].listIndices[j].nameIndex
		      ]];
	    }
	  } else {
	    // TODO: need to remove the integer index code.
	    assert(0);
#if 0
	    // the list index is an integer
	    if (rvInfoVector[i].listIndices[j].intIndex >= 
		GM_Parms.mixturesMap.size()) {
		if (!allocateIfNotThere) {
		  error("Error: RV \"%s\" at frame %d (line %d), mixture index (%d) too large\n",
			rvInfoVector[i].name.c_str(),
			rvInfoVector[i].frame,
			rvInfoVector[i].fileLineNumber,
			rvInfoVector[i].listIndices[j].intIndex);
		} else {
		  error("Can't allocate with integer mixture index");
		}
	      } else {
		// otherwise add it
		rv->conditionalMixtures[j].mixture =
		  GM_Parms.mixtures[
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
	  rv->conditionalMixtures[j].direct = false;
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
	      rv->conditionalMixtures[j].mapping.dtMapper =
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
		rv->conditionalMixtures[j].mapping.dtMapper =
		  GM_Parms.dts[rvInfoVector[i].listIndices[j].intIndex];
	      }
	  }
	  // 
	  // check that the DT matches the number of current parents and their cardinalities.

	  // Check to make sure that the decision tree (which is used to map sets of parent values down
	  // to the integer that indexes into the collection) has a number of features/parents that
	  // matches the current set of conditional parents.
	  if (rv->conditionalMixtures[j].mapping.dtMapper->numFeatures() != 
	      rvInfoVector[i].conditionalParents[j].size()) {
	    error("Error: RV \"%s\" at frame %d (line %d), Decision Tree '%s' wants %d parents, but number of current parents is %d.\n",
		  rvInfoVector[i].name.c_str(),
		  rvInfoVector[i].frame,
		  rvInfoVector[i].fileLineNumber,
		  rv->conditionalMixtures[j].mapping.dtMapper->name().c_str(),
		  rv->conditionalMixtures[j].mapping.dtMapper->numFeatures(),
		  rvInfoVector[i].conditionalParents[j].size());
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
	      rv->conditionalMixtures[j].mapping.collection =
		GM_Parms.ncls[
		   GM_Parms.nclsMap[
				rvInfoVector[i].listIndices[j].collectionName
		  ]];
	      // make sure that the pointer table is filled in for this name.
	      rv->conditionalMixtures[j].mapping.collection->fillMxTable();
	      // Might as well do this here. Check that all component densities
	      // in the collection are of the right dimensionality.
	      // This check is important for ContRV::probGivenParents().
	      const unsigned rvDim = 
		(rvInfoVector[i].rvFeatureRange.lastFeatureElement - 
		 rvInfoVector[i].rvFeatureRange.firstFeatureElement)+1;
	      // make sure all component distributions in the collection have the
	      // same dimensionality as the RV.
	      for (unsigned u=0;u<
		     rv->conditionalMixtures[j].mapping.collection->mxSize();
		   u++) {
		if ( // we only check the dimension for those component distributions that have it.
		    (rv->conditionalMixtures[j].mapping.collection->mx(u)->mixType
		     != MixtureCommon::ci_zeroScoreMixture)
		    &&
		    (rv->conditionalMixtures[j].mapping.collection->mx(u)->mixType
		     != MixtureCommon::ci_unityScoreMixture)
		    && 
		    (rv->conditionalMixtures[j].mapping.collection->mx(u)->dim() 
		     != rvDim)) {
		  error("Error: RV \"%s\" at frame %d (line %d), dimensionality %d, conditional parent %d, collection \"%s\" specifies Mixture \"%s\" at pos %d with wrong dimension %d\n",
		      rvInfoVector[i].name.c_str(),
		      rvInfoVector[i].frame,
		      rvInfoVector[i].fileLineNumber,
		      rvDim,
		      j,
		      rvInfoVector[i].listIndices[j].collectionName.c_str(),
		      rv->conditionalMixtures[j].mapping.collection->mx(u)->name().c_str(),
		      u,
		      rv->conditionalMixtures[j].mapping.collection->mx(u)->dim());
		}
	      }
	    }
	}
      }
    }
  }

  // now go through factor list and associate the 
  // factor with appropriate GMTK objects (DTs, tables, etc.) when
  // needed.
  for (unsigned factorNo=0;factorNo<factorList.size();factorNo++) {
    FactorInfo& factor = factorList[factorNo];

    // TODO: finish this function.


    if (factor.fType == FactorInfo::ft_directionalConstraint) {
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
      if (rvInfoVector[i].rvDisp != RVInfo::d_hidden) {
	// then observed
	if (rvInfoVector[i].rvFeatureRange.filled == RVInfo::FeatureRange::fr_Range) {
	  // then observed, value from a feature range. Need to check to make sure
	  // it corresponds to a true discrete value.
	  if (!globalObservationMatrix.elementIsDiscrete(rvInfoVector[i].rvFeatureRange.firstFeatureElement)) {
	    if (globalObservationMatrix.numDiscrete() > 0) 
	      error("ERROR: discrete observed random variable '%s', frame %d, line %d, specifies a feature element %d:%d that is out of discrete range ([%d:%d] inclusive) of observation matrix",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rvInfoVector[i].rvFeatureRange.firstFeatureElement,
		    rvInfoVector[i].rvFeatureRange.firstFeatureElement,
		    globalObservationMatrix.numContinuous(),
		    globalObservationMatrix.numFeatures()-1);
	    else
	      error("ERROR: discrete observed random variable '%s', frame %d, line %d, specifies a feature element %d:%d for an observation matrix with zero discrete features.",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rvInfoVector[i].rvFeatureRange.firstFeatureElement,
		    rvInfoVector[i].rvFeatureRange.firstFeatureElement,
		    globalObservationMatrix.numContinuous(),
		    globalObservationMatrix.numFeatures()-1);
	  }
	}
      }
    } else { // (rvInfoVector[i].rvType == RVInfo::t_continuous) {
      if (rvInfoVector[i].rvDisp != RVInfo::d_hidden) {
	if (rvInfoVector[i].rvFeatureRange.lastFeatureElement >=  globalObservationMatrix.numContinuous())
	      error("ERROR: continuous observed random variable '%s', frame %d, line %d, specifies feature elements %d:%d that are out of continuous range ([%d:%d] inclusive) of observation matrix",
		    rvInfoVector[i].name.c_str(),
		    rvInfoVector[i].frame,
		    rvInfoVector[i].fileLineNumber,
		    rvInfoVector[i].rvFeatureRange.firstFeatureElement,
		    rvInfoVector[i].rvFeatureRange.lastFeatureElement,
		    0,
		    globalObservationMatrix.numContinuous()-1);
      }
    }

    // checks common for both discrete and continuous RVs.

    // if weight comes from observation matrix, make sure it indexes into
    // a valid index and a float value.
    for (unsigned wt=0;wt<rvInfoVector[i].rvWeightInfo.size();wt++) {

      if (rvInfoVector[i].rvWeightInfo[wt].penalty.wt_Status 
	  == RVInfo::WeightInfo::WeightItem::wt_Observation) {
	if (rvInfoVector[i].rvWeightInfo[wt].penalty.lastFeatureElement >= globalObservationMatrix.numContinuous()) {
	  error("ERROR: random variable '%s', frame %d, line %d, weight attribute at position %d has penalty observation feature element %d:%d that is out of continuous range ([%d:%d] inclusive) of observation matrix",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		wt,
		rvInfoVector[i].rvWeightInfo[wt].penalty.firstFeatureElement,
		rvInfoVector[i].rvWeightInfo[wt].penalty.lastFeatureElement,
		0,
		globalObservationMatrix.numContinuous()-1);
	}
      }
      if (rvInfoVector[i].rvWeightInfo[wt].scale.wt_Status 
	  == RVInfo::WeightInfo::WeightItem::wt_Observation) {
	if (rvInfoVector[i].rvWeightInfo[wt].scale.lastFeatureElement >= globalObservationMatrix.numContinuous()) {
	  error("ERROR: random variable '%s', frame %d, line %d, weight attribute at position %d has scale observation feature element %d:%d that is out of continuous range ([%d:%d] inclusive) of observation matrix",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		wt,
		rvInfoVector[i].rvWeightInfo[wt].scale.firstFeatureElement,
		rvInfoVector[i].rvWeightInfo[wt].scale.lastFeatureElement,
		0,
		globalObservationMatrix.numContinuous()-1);
	}
      }
      if (rvInfoVector[i].rvWeightInfo[wt].shift.wt_Status 
	  == RVInfo::WeightInfo::WeightItem::wt_Observation) {
	if (rvInfoVector[i].rvWeightInfo[wt].shift.lastFeatureElement >= globalObservationMatrix.numContinuous()) {
	  error("ERROR: random variable '%s', frame %d, line %d, weight attribute at position %d has shift observation feature element %d:%d that is out of continuous range ([%d:%d] inclusive) of observation matrix",
		rvInfoVector[i].name.c_str(),
		rvInfoVector[i].frame,
		rvInfoVector[i].fileLineNumber,
		wt,
		rvInfoVector[i].rvWeightInfo[wt].shift.firstFeatureElement,
		rvInfoVector[i].rvWeightInfo[wt].shift.lastFeatureElement,
		0,
		globalObservationMatrix.numContinuous()-1);
	}
      }

    }
  }
}





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
		   vector<RV*> &unrolledVarSet)
{
  // a map from the r.v. name and new frame number to
  // position in the new unrolled array of r.v.'s
  map < RVInfo::rvParent, unsigned > posOfParentAtFrame;
  return unroll(timesToUnroll,unrolledVarSet,posOfParentAtFrame);
}
void
FileParser::unroll(unsigned timesToUnroll,
		   vector<RV*> &unrolledVarSet,
		   map < RVInfo::rvParent, unsigned >& posOfParentAtFrame)
{

  // a map from the new r.v. to the corresponding rv info
  // vector position in the template.
  map < RV *, unsigned > infoOf;

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

    unrolledVarSet[uvsi] = rvInfoVector[tvi].rv->cloneRVShell();

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

      unrolledVarSet[uvsi] = rvInfoVector[tvi].rv->cloneRVShell();
      unrolledVarSet[uvsi]->timeFrame = frame;

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

    unrolledVarSet[uvsi] = rvInfoVector[tvi].rv->cloneRVShell();
    unrolledVarSet[uvsi]->timeFrame = frame;

    infoOf[unrolledVarSet[uvsi]] = tvi;
    RVInfo::rvParent p(unrolledVarSet[uvsi]->name(),frame);
    posOfParentAtFrame[p] = uvsi;

    uvsi++;
    tvi++;
  }
  
  // now we have all the rv,s but with the wrong parents, fix that.

  for (uvsi=0;uvsi<unrolledVarSet.size();uvsi++) {

    const unsigned frame = unrolledVarSet[uvsi]->timeFrame;

    // get the info objectd for this one.
    RVInfo* info = &rvInfoVector[infoOf[unrolledVarSet[uvsi]]];
    
    // build the switching parents list.
    vector<RV *> sparents;
    for (unsigned j=0;j<info->switchingParents.size();j++) {

      // grab a pointer to the parent in the template
      RV* template_parent_rv =
	info->rv->switchingParentsVec()[j];

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
      RV *unrolled_parent_rv = unrolledVarSet[(*it).second];
      if (!unrolled_parent_rv->discrete()) {
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
      DiscRV* dunrolled_parent_rv = 
	(DiscRV*) unrolled_parent_rv;
      
      if (dunrolled_parent_rv->cardinality != 
	  RV2DRV(template_parent_rv)->cardinality) {
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
	      dunrolled_parent_rv->cardinality,
	      RV2DRV(template_parent_rv)->cardinality);
      }

      // add
      sparents.push_back(unrolled_parent_rv);

    }

    // now build the conditional parents list.
    vector<vector<RV * > > cpl(info->conditionalParents.size());
    for (unsigned j=0;j<info->conditionalParents.size();j++) {
      for (unsigned k=0;k<info->conditionalParents[j].size();k++) {


	// grab a pointer to the parent in the template
	RV* template_parent_rv =
	  info->rv->condParentsVec(j)[k];

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
	RV *unrolled_parent_rv = unrolledVarSet[(*it).second];
	if (!unrolled_parent_rv->discrete()) {
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
	DiscRV* dunrolled_parent_rv = 
	  (DiscRV*) unrolled_parent_rv;
      
	if (dunrolled_parent_rv->cardinality != 
	    RV2DRV(template_parent_rv)->cardinality) {
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
		dunrolled_parent_rv->cardinality,
		RV2DRV(template_parent_rv)->cardinality);
	}

	// add
	cpl[j].push_back(unrolled_parent_rv);
    
      }
    }
    unrolledVarSet[uvsi]->setParents(sparents,cpl);
  }

  // We *may* at some point in the future want to add argument to
  // optionally call the following if we want the neighbors members of
  // RVS to respect the factor's that were defined in the .str
  // file. Note that for long unrollings, this can add considerable
  // expense (sihce this will compute the list of factors), so we
  // leave this commented out and assume the caller of unroll() will
  // do the right thing if need be.
  // 
  // addUndirectedFactorEdges(unrolledVarSet,posOfParentAtFrame);

}


/*-
 *-----------------------------------------------------------------------
 * addUndirectedFactorEdges()
 *
 *  Given an unrolled graph (and the corresponding data structures
 *  returned by the unroll routine), add any undirected edges
 *  according to any undirected factors specified in the structure file.
 *
 *  Note that the semantics of factors in a graph is that the
 *  factor is added relative to the frame in the *unrolled*
 *  graph. This means that if a factor is defined with a random
 *  variable foo(-3) in the template, and that factor when
 *  unrolled exists at frame n, then it expects (and affects) random
 *  variable foo(n-3) in the unrolled graph.  Note that this is
 *  similar behavior to that of random variables that are defined in a
 *  frame, i.e., it is the parents of a random variable that are
 *  relative to the child in the unrolled graph (so we first create
 *  the graph of children, and then create the parents relative to
 *  those children). Here, we first create the graph of nodes in
 *  frames, and then create the factors relative to those frames.
 *
 * Preconditions:
 *     Must be called from unroll() after all the random variales have been set up.
 *     All RV's neighbors members are presumably empty.
 *     'timesToUnroll' argument must be same as corresponding call with unroll().
 *
 * Postconditions:
 *     The neighbors members of the random variables in the unrolled graph
 *     have been adjusted to reflect the factors given in the structure file.
 *
 * Side Effects:
 *     Changes neighbors members of RVs.
 *
 * Results:
 *     None
 *
 *----------------------------------------------------------------------- 
 */
void
FileParser::addUndirectedFactorEdges(unsigned timesToUnroll,
				     vector<RV*> &rvs,
				     map < RVInfo::rvParent, unsigned >& pos,
				     // return value
				     vector < set < RV* > >& factorArray)
{
  // the set of rvs corresponding to a given factor.
  set<RV*> rv_factor;
  // 
  // Add neighbors structures in unrolled graph based on factors.
  for (unsigned factorNo=0;factorNo<factorList.size();factorNo++) {
    FactorInfo& factor = factorList[factorNo];
    if (frameInTemplateP(factor.frame)) {
      const unsigned offset = 0;
      completeRVsInFactor(factor,offset,rvs,pos,rv_factor);
      factorArray.push_back(rv_factor);
    } else if (frameInTemplateC(factor.frame)) {
      // do this for all original C partitions.
      for (unsigned i=0;i<(timesToUnroll+1);i++) {
	const unsigned offset = 
	  i*numFramesInC();
	completeRVsInFactor(factor,offset,rvs,pos,rv_factor);	
	factorArray.push_back(rv_factor);
      }
    } else {
      assert (frameInTemplateE(factor.frame));
      const unsigned offset = 
	timesToUnroll*numFramesInC();
      completeRVsInFactor(factor,offset,rvs,pos,rv_factor);
      factorArray.push_back(rv_factor);
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * completeRVsInFactor()
 *  Given a file parser factor (which, as you should know, is not nec. a max clique),
 *  and an offset, complete (connect all neighbors) of all the random variables
 *  in the unrolled graph according to that factor.
 *
 *
 *
 * Preconditions:
 *     Must be called from unroll() after all the random variales have been set up.
 *     All RV's neighbors members are presumably empty.
 *
 * Postconditions:
 *     The neighbors members of the random variables in the unrolled graph
 *     have been adjusted to reflect the factors given in the structure file.
 *
 * Side Effects:
 *     Changes neighbors members of RVs.
 *
 * Results:
 *     None
 *
 *----------------------------------------------------------------------- 
 */
void
FileParser::completeRVsInFactor(FactorInfo& factor,
				const unsigned offset,
				vector<RV*> &rvs,
				map < RVInfo::rvParent, unsigned >& pos,
				// return value
				set<RV*>& rv_factor)
{

  // form a set of the random variables themselves.
  rv_factor.clear();

  map < RVInfo::rvParent , unsigned >::iterator it;
  for (unsigned varNo=0;varNo<factor.variables.size();varNo++) {
    RVInfo::rvParent pp(factor.variables[varNo].first,
			factor.variables[varNo].second 
			+ factor.frame
			+ offset);

    if ((it = pos.find(pp)) == pos.end()) {
      coredump("INTERNAL ERROR: Can't find random variable %s(%d) in unrolled collection.\n",
	       pp.first.c_str(),pp.second);
    }
    rv_factor.insert(rvs[(*it).second]);
  }
  MaxClique::makeComplete(rv_factor);

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
  
  for (unsigned i=0;i<factorList.size();i++) {
    // for each factor, write out
    //    1. Name of factor
    //    2. frame where defined
    //    3. Number of random variables
    //    4. the list of random variables (name,frame) pairs
    os.write(factorList[i].name,"name");
    os.write(factorList[i].frame,"frame");
    os.write(factorList[i].variables.size(),"num vars");
    for (unsigned j=0;j<factorList[i].variables.size();j++) {
      os.write(factorList[i].variables[j].first,"nm");
      os.write(factorList[i].variables[j].second,"frm");
    }
    
    // could write out more factor information, but we only write
    // out what could effect the triangulation (and the only
    // constraint is that a factor is complete, so the rest is
    // ok to change under a current triangulation).

  }
  

  os.write(_firstChunkframe);
  os.write(_lastChunkframe);
  os.nl();
  os.write(TRIFILE_END_OF_ID_STRING);
  os.nl();

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
FileParser::readAndVerifyGMId(iDataStreamFile& is,const bool checkCardinality)
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
    if (checkCardinality) {
      if (uval != rvInfoVector[i].rvCard) return false;
    }

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


  for (unsigned i=0;i<factorList.size();i++) {
    // for each factor, check
    //    1. Name of factor
    //    2. frame where defined
    //    3. Number of random variables
    //    4. the list of random variables (name,frame) pairs


    if (!is.read(nm)) return false;
    // don't bother checking that name matches, since it curently doesn't matter.
    // if (nm != factorList[i].name) return false;

    if (!is.read(uval)) return false;
    if (uval != factorList[i].frame) return false;

    if (!is.read(uval)) return false;
    if (uval != factorList[i].variables.size()) return false;

    for (unsigned j=0;j<factorList[i].variables.size();j++) {
      if (!is.read(nm)) return false;
      if (nm != factorList[i].variables[j].first) return false;


      if (!is.read(ival)) return false;
      if (ival != factorList[i].variables[j].second) return false;

    }

    // we could check the rest of the factor, but since it doesn't
    // effect the triangulation, we don't (meaning the user
    // can change some of the information in the factor w/o needing to 
    // retriangulate).

  }


  if (!is.read(uval)) return false;
  if (uval != _firstChunkframe) return false;

  if (!is.read(uval)) return false;
  if (uval != _lastChunkframe) return false;

  if (!is.read(nm)) return false;
  if (nm != TRIFILE_END_OF_ID_STRING) return false;

  // all checked out ok
  return true;
}


/*-
 *-----------------------------------------------------------------------
 * ensureS_SE_E_NE,
 *
 *    ensure links are "south", "south east", "east", or "north east",
 *    meaning that there is a numeric ordering on the nodes such that
 *    any parents of a node at a particular numeric position have
 *    their position earlier in the ordering. Note that this routine
 *    is not necessary, as the variables in a .str file may be in any
 *    order (as long as there are no directed cycles).
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
