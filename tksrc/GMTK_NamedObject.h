/*-
 * GMTK_NamedObject.h
 *      parent class for objects with string names.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#ifndef GMTK_NAMEDOBJECT_H
#define GMTK_NAMEDOBJECT_H

#include <string>

#include "fileParser.h"

#define NAMED_OBJECT_DEFAULT_NAME "-"

class NamedObject {

protected:
  ///////////////////////////////////////////////////////////  
  // the name
  string _name;


public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  NamedObject(const string& nm) : _name(nm) {} 
  NamedObject() { _name = NAMED_OBJECT_DEFAULT_NAME; } 


  ///////////////////////////////////////////////////////////    
  void read(iDataStreamFile& is) {
    is.read(_name,"NamedObject::read");
  }

  ///////////////////////////////////////////////////////////    
  void write(oDataStreamFile& os) {
    if (_name.size() == 0)
      os.write("-","NamedObject::write");
    else
      os.write(_name,"NamedObject::wriite");
  }

  ////////////////////////////////////////
  void setName(const string& nm) { _name = nm; }
  ////////////////////////////////////////
  const string& name() { return _name; }

};



#endif 
