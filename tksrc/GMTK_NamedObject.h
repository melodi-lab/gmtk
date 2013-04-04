/*-
 * GMTK_NamedObject.h
 *      parent class for objects with string names.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
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
    is.read(_name,"Can't read name");
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
