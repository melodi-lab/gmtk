//static char *rcsid = "$Id$";

#include "sArray.h"

// basic structure for data record: contiguous float and int arrays

class GM_rec {

  sArray<float> fval; // float part of input record
  sArray<int> ival;  // int part of input record

};
