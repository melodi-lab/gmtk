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

/*
#include "sArray.h"
#include "logp.h"
#include "GMTK_GM.h" 
*/

class FileParser
{


  //  GMTK_GM& gm; 

public:

  enum {
    TT_INT=0,
    TT_FLT=1,
    TT_KWD=2,
    TT_ID=3,
    TT_OP=4,
    TT_UNDEF=5
  };


  //  FileParser(const char *const fileName,
  //	     GMTK_GM& gm);


  FileParser(const char *const fileName);

  parseGraphicalModel();
  parseFrameList();
  parseFrame();
  parseRandomVariableList();
  parseRandomVariable();
  parseRandomVariableBody();


	     


};

#endif
