/*
 * GMTK_GMTemplate.h
 * Provides a class that supports GMTK templates.
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

#ifndef GMTK_GMTEMPLATE_H
#define GMTK_GMTEMPLATE_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_RandomVariable.h"

class FileParser;
class GraphicalModel;

class GMTemplate
{
  friend class FileParser;
  friend class GraphicalModel;

  // The number of frames in this template
  unsigned numFrames;
  // First chunk frame
  unsigned firstChunkFrame;
  // Last Chunk Frame
  unsigned lastChunkFrame;

  struct frame {
    // list of random variables for this frame
    vector <RandomVariable*> rvs;

  };


  // the set of frames;
  vector < frame > frames;
public:

  GMTemplate() {numFrames=0;}
  ~GMTemplate() {}


  ///////////////////////////////////////////////////////////////////
  // prints out the contents of this template (i.e., variables and
  // their parents.
  void print() {
    for (unsigned i=0;i<numFrames;i++) {
      for (unsigned j=0;j<frames[i].rvs.size();j++) {
	printf("frame %d, variable %d: %s\n",
	       i,j,frames[i].rvs[j]->name().c_str());
      }
    }
  }
};

#endif

