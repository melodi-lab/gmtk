/*-
 * GMTK_dlopenDeterministicMappings
 *      Include this to write mapping functions to be dynamically linked
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu> & Richard Rogers <rprogers@uw.edu>
 * 
 * 
 * Copyright (c) 2012, < fill in later >
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


#ifndef GMTK_DLOPEN_DETERMINISTIC_MAPPINGS_H
#define GMTK_DLOPEN_DETERMINISTIC_MAPPINGS_H

#include "GMTK_DiscRV.h"
#include "GMTK_GMParms.h"

// this defines the name of the routine in the C++ namespace
#define DETERMINISTIC_MAPPER_C_CODE_NAME(name) cFunctionDeterministicMapping_##name

// use this value if a C decision tree (CDT) has a variale number of
// features, in which case the C routine needs to check the size of
// the variable array (macros for this are given below).
#define CDT_VARIABLE_NUMBER_FEATURES (~0x0)



vector<CFunctionMapperType>   mapperFunctions;
vector<unsigned> mapperNumFeatures;
vector<char const*>   mapperNames;

#define REGISTER_MAPPER(name,numFeatures)       \
  mapperNames.push_back("user_internal:"#name); \
  mapperNumFeatures.push_back(numFeatures);     \
  mapperFunctions.push_back(&cFunctionDeterministicMapping_##name);

#define DEFINE_DETERMINISTIC_MAPPER_C_CODE(name,numFeatures) \
DiscRVType                                                   \
cFunctionDeterministicMapping_##name(const vector< RV* >& parent_variables,const RV* const child_rv)


///// Macros for easy access to variable names.

// parent values
#define p0 (RV2DRV(parent_variables[0])->val)
#define p1 (RV2DRV(parent_variables[1])->val)
#define p2 (RV2DRV(parent_variables[2])->val)
#define p3 (RV2DRV(parent_variables[3])->val)
#define p4 (RV2DRV(parent_variables[4])->val)
#define p5 (RV2DRV(parent_variables[5])->val)
#define p6 (RV2DRV(parent_variables[6])->val)
#define p7 (RV2DRV(parent_variables[7])->val)
#define p8 (RV2DRV(parent_variables[8])->val)
#define p9 (RV2DRV(parent_variables[9])->val)
#define p10 (RV2DRV(parent_variables[10])->val)
#define p11 (RV2DRV(parent_variables[11])->val)
#define p12 (RV2DRV(parent_variables[12])->val)
#define p13 (RV2DRV(parent_variables[13])->val)
#define p14 (RV2DRV(parent_variables[14])->val)
#define p15 (RV2DRV(parent_variables[15])->val)
#define p16 (RV2DRV(parent_variables[16])->val)
#define p17 (RV2DRV(parent_variables[17])->val)
#define p18 (RV2DRV(parent_variables[18])->val)
#define p19 (RV2DRV(parent_variables[19])->val)
#define p20 (RV2DRV(parent_variables[20])->val)
#define p21 (RV2DRV(parent_variables[21])->val)
#define p22 (RV2DRV(parent_variables[22])->val)
#define p23 (RV2DRV(parent_variables[23])->val)
#define p24 (RV2DRV(parent_variables[24])->val)
#define p25 (RV2DRV(parent_variables[25])->val)
#define p26 (RV2DRV(parent_variables[26])->val)
#define p27 (RV2DRV(parent_variables[27])->val)
#define p28 (RV2DRV(parent_variables[28])->val)
#define p29 (RV2DRV(parent_variables[29])->val)
#define p30 (RV2DRV(parent_variables[30])->val)
#define p31 (RV2DRV(parent_variables[31])->val)

// Grab a parent value using a variable index, must make sure that
// index is in range, or a run-time error will occur.
#define numParents ((parent_variables.size()))

// Grab a parent value using a variable index, must make sure that
// index is in range, or a run-time error will occur.
#define par(i) (RV2DRV(parent_variables[(i)])->val)

// cardinality of parent (cp).
#define cp0 (RV2DRV(parent_variables[0])->cardinality)
#define cp1 (RV2DRV(parent_variables[1])->cardinality)
#define cp2 (RV2DRV(parent_variables[2])->cardinality)
#define cp3 (RV2DRV(parent_variables[3])->cardinality)
#define cp4 (RV2DRV(parent_variables[4])->cardinality)
#define cp5 (RV2DRV(parent_variables[5])->cardinality)
#define cp6 (RV2DRV(parent_variables[6])->cardinality)
#define cp7 (RV2DRV(parent_variables[7])->cardinality)
#define cp8 (RV2DRV(parent_variables[8])->cardinality)
#define cp9 (RV2DRV(parent_variables[9])->cardinality)
#define cp10 (RV2DRV(parent_variables[10])->cardinality)
#define cp11 (RV2DRV(parent_variables[11])->cardinality)
#define cp12 (RV2DRV(parent_variables[12])->cardinality)
#define cp13 (RV2DRV(parent_variables[13])->cardinality)
#define cp14 (RV2DRV(parent_variables[14])->cardinality)
#define cp15 (RV2DRV(parent_variables[15])->cardinality)
#define cp16 (RV2DRV(parent_variables[16])->cardinality)
#define cp17 (RV2DRV(parent_variables[17])->cardinality)
#define cp18 (RV2DRV(parent_variables[18])->cardinality)
#define cp19 (RV2DRV(parent_variables[19])->cardinality)
#define cp20 (RV2DRV(parent_variables[20])->cardinality)
#define cp21 (RV2DRV(parent_variables[21])->cardinality)
#define cp22 (RV2DRV(parent_variables[22])->cardinality)
#define cp23 (RV2DRV(parent_variables[23])->cardinality)
#define cp24 (RV2DRV(parent_variables[24])->cardinality)
#define cp25 (RV2DRV(parent_variables[25])->cardinality)
#define cp26 (RV2DRV(parent_variables[26])->cardinality)
#define cp27 (RV2DRV(parent_variables[27])->cardinality)
#define cp28 (RV2DRV(parent_variables[28])->cardinality)
#define cp29 (RV2DRV(parent_variables[29])->cardinality)
#define cp30 (RV2DRV(parent_variables[30])->cardinality)
#define cp31 (RV2DRV(parent_variables[31])->cardinality)

// Grab a parent cardinality value using a variable index, must make
// sure that index is in range, or a run-time error will occur.
#define cardPar(i) (RV2DRV(parent_variables[(i)])->cardinality)

// child cardinality
#define cc (RV2DRV(child_rv)->cardinality)

#endif 
