/*-
 * GMTK_RVInfo.cc
 *     RV generic information
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2007, < fill in later >
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

#include "GMTK_RVInfo.h"


#include "GMTK_RV.h"

VCID("$Header$")


/*-
 *-----------------------------------------------------------------------
 * RVInfo::clear()
 *   clear out the current RV structure to a generic known state, and 
 *   also clear up any memory used by this object. Useful for construction/destruction,
 *   or when we are parsing and encounter a new RV.
 * 
 * Preconditions:
 *      None
 *
 * Postconditions:
 *      None
 *
 * Side Effects:
 *      changes internal member variables.
 *
 * Results:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void 
RVInfo::clear() {
    name.erase();

    rvType = t_unknown;
    rvDisp = d_unknown;
    rvFeatureRange.clear();
    rvWeightInfo.clear();
    eliminationOrderHint = 0.0;
    variablePositionInStrFile = -1;

    switchingParents.clear();
    switchMapping.clear();

    conditionalParents.clear();
    discImplementations.clear();
    contImplementations.clear();
    listIndices.clear();

    isDeterministic = isSparse = false;

    symbolTableCollectionName.clear();
    symbolTable = NULL;

    if (rv) delete rv;
    rv = NULL;
}



/*-
 *-----------------------------------------------------------------------
 * RVInfo::rvParent
 *   parse the string given as an argument into an RV parent object.
 *   
 *   
 * 
 * Preconditions:
 *      The string should be of the form "foo(3)" or "Variableame(-7)",
 *      i.e., a string followed by an integer (zero, positive, or negative)
 *      enclosed in parentheses. Note that this parsing code exists in
 *      the lexer and parser for structure files, but is given here
 *      to make it easy for command line user input of random variable specs.
 *
 * Postconditions:
 *      None
 *
 * Side Effects:
 *      None
 *
 * Results:
 *      creates a new object and returns it.
 *
 *-----------------------------------------------------------------------
 */
RVInfo::rvParent
RVInfo::parseRVParent(const char* const str)
{
  // A RV spec consists if a case-sensitive identifier (starting with
  // an apha character), followed by optional space, then an signed integer
  // surrounded by parents. Examples include
  //  "wordTransition(3)"
  //  "foo(+3)"
  //  "bar(0)"
  //  "baz(-3)"
  // In this routine, moreover, we allow an optional missing integer in parents,
  // so "foo" is the same as "foo(0)"

  const char *buffp = str;

  // skip any initial space
  while (!isspace(*buffp) &&  (*buffp) != '\n') {
    buffp++;
  }   

  // read the identifier, accept any alphanum char.
  const char * const startingPosition = buffp;
  if (isalpha(*buffp)) {
    buffp++;
  } else {
    error("ERROR: parsing random variable specification '%s', must begin with an alphabetical character.\n");
  }
  while (isalnum(*buffp)) {
    buffp++;
  }
  // we should be at end of identifier
  if (buffp == str) {
    error("ERROR: parsing random variable specification '%s', must have non-zero length.\n");
  }

  string res_name;
  res_name.append(startingPosition,buffp-startingPosition);

  // skip space
  while (!isspace(*buffp) &&  (*buffp) != '\n') {
    buffp++;
  }   

  // we should be at end of string or have a '(' characer.

  int offset = 0;
  if (*buffp) {

    // we're not at the end of the string

    if (*buffp == '(') 
      buffp++;
    else 
      error("ERROR: parsing random variable spec '%s', must be form 'id(int)'.\n",str);

    // skip space
    while (!isspace(*buffp) &&  (*buffp) != '\n') {
      buffp++;
    }   
    // read the integer part
    char *end_ptr;
    long l = strtol(buffp,&end_ptr,0);
    if (end_ptr == buffp) {
      error("ERROR: parsing random variable spec '%s' and can't form integer.\n",str);
    } else
      buffp = end_ptr;
    offset = (int)l;

    // check for final paren.
    while (!isspace(*buffp) &&  (*buffp) != '\n') {
      buffp++;
    }   
    if (*buffp == ')') 
      buffp++;
    else 
      error("ERROR: parsing random variable spec '%s', must be form 'id(int)'.\n",str);
  }

  return rvParent (res_name,offset);

}
