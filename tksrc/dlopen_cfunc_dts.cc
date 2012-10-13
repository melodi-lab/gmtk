
// to get helper macro definitions
#include "GMTK_dlopenDeterministicMappings.h"


/*
 *
 * 
 * 
 * Define the registerMappers function, containing a REGISTER_MAPPER
 * for each mapping function.
 *
 * Compile with:
 *   g++ -O3 -I. -I../miscSupport -I../featureFileIO -I../IEEEFloatingpoint -DHAVE_CONFIG_H -fPIC -c dlopen_cfunc_dts.cc
 *   g++ -O3 -shared -o dlopen_cfunc_dts.so dlopen_cfunc_dts.o 
 *
 * This builds a shared library from your source code.
 * Then run gmtkCMD ... -map1 dlopen_cfunc_dts.so
 * This will cause the GMTK program to load the shared library and
 * make your mapping functions available for use in models. You can
 * load up to 5 libraries (-map1 through -map5), and each library
 * can contain multiple mapping functions.
 * 
 */

// Use DEFINE_DETERMINISTIC_MAPPER_C_CODE to define mapping functions.
// Pass in the name of the mapping function and the # of input features.
// You can pass CDT_VARIABLE_NUMBER_FEATURES to indicate the function
// takes a variable # of input features. numParents will indicate the
// number of input features. The mapping functions will be available as 
// "user_internal:name" (eg, user_internal:ajitMapping).

// You can access the value of the features as pX, where 0 <= X <= 31.
// You can access an arbitrary feature using par(i), but you must 
// ensure i < numParents. Feature cardinality is available as cpX or
// cardPar(i). cc gives the cardinality of the mapping function's output.
// These values are all unsigned ints.

#define AJIT_MAPPING_NUM_FEATURES 10
DEFINE_DETERMINISTIC_MAPPER_C_CODE(ajitMapping,AJIT_MAPPING_NUM_FEATURES)
{ 
  const DiscRVType SHIFT_ZERO=0;

  return (p0 > 0 ? ((  ( (p1+1 - 17*p7 - 18*p8 + p9) + (p0>>1))/p0) + 1*(p6-SHIFT_ZERO) >= 3 ? ( ( ((p1 + 1 - 17*p7 - 18*p8 + p9) + (p0>>1))/p0) + 1*(p6-SHIFT_ZERO) < 4 ? 1 : 0) : 0) : 0) 
    ||
    (p3 > 0 ? (( ((p4+19 - 17*p7 - 18*p8 + p9) + (p3>>1))/p3) + 1*(p6-SHIFT_ZERO) >= 3 ? ( ( ((p4+19 - 17*p7 - 18*p8 + p9) + (p3>>1))/p3 + 1*(p6-SHIFT_ZERO)) < 4 ? 1 : 0) : 0) : 0);
}


// Here's another simple exmaple mapping function

#define RICHARD_MAPPING_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(richardMapping,RICHARD_MAPPING_NUM_FEATURES)
{
  return (p0 + cp0) % cc;
}


// After loading the shared library, the GMTK program will call the 
// registerMappers function to make the defined mapping functions 
// available for use by models. Thus each library must define registerMappers()
// with a REGISTER_MAPPER for each mapping function.

extern "C" void
registerMappers() {
  REGISTER_MAPPER(ajitMapping,AJIT_MAPPING_NUM_FEATURES)
  REGISTER_MAPPER(richardMapping,RICHARD_MAPPING_NUM_FEATURES)
}

