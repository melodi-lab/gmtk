/*
    Simple informational message system
    Jeff Bilmes <bilmes@ee.washington.edu>
    $Header$
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "debug.h"

unsigned IM::globalMessageLevel = IM::Nano;
bool IM::globalFlush = true;
#ifdef GMTK_ARG_VERB_DEF_VAL
unsigned IM::globalModuleLevel[(unsigned)ModuleCount] = {GMTK_ARG_VERB_DEF_VAL,...};
#else
unsigned IM::globalModuleLevel[(unsigned)ModuleCount] = {IM::Default,};
#endif
const char*IM::moduleString[(unsigned)ModuleCount] = 
  {"default","inference","inference-memory", "training", "triangulation","boundary","unrolling","printing","modelinfo","obsfile","obsstream"};

#ifdef MAIN

class FOO : public IM {
  
public:
  int i;
  
  FOO() { i = 3; }
  void bar() { infoMsg(Low,"This is a test, ml = %d\n",msgLevel()); }

};

int main()
{
  

  FOO f;
  
  f.setMsgLevel(IM::Nano);
  f.setMsgLevel(IM::Huge);
  f.msgsOff();  
  f.msgsOn();
  f.bar();

  return 0;
}

#endif


