/*-
 * GMTK_RngDecisionTree
 *      General class to map from vectors of integers to some 
 *      basic type (such as int, float, etc.). Uses bp_ranges
 *      to form the queries.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu> & Chris Bartels <bartels@ee.washington.edu>
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


#ifndef GMTK_C_FUNCTION_DETERMINISTIC_MAPPINGS_H
#define GMTK_C_FUNCTION_DETERMINISTIC_MAPPINGS_H

wid#define DETERMINISTIC_MAPPER_C_CODE_NAME(name) cFunctionDeterministicMapping_##name

#define DEFINE_DETERMINISTIC_MAPPER_C_CODE(name,numFeatures) \
DiscRVType cFunctionDeterministicMapping_##name(const vector< RV* >& parent_variables,const RV* const child_rv)


///// Macros for easy access to variable names.

// parent values
#define p0 (RV2DRV(variables[0])->val)
#define p1 (RV2DRV(variables[1])->val)
#define p2 (RV2DRV(variables[2])->val)
#define p3 (RV2DRV(variables[3])->val)
#define p4 (RV2DRV(variables[4])->val)
#define p5 (RV2DRV(variables[5])->val)
#define p6 (RV2DRV(variables[6])->val)
#define p7 (RV2DRV(variables[7])->val)
#define p8 (RV2DRV(variables[8])->val)
#define p9 (RV2DRV(variables[9])->val)
#define p10 (RV2DRV(variables[10])->val)
#define p11 (RV2DRV(variables[11])->val)
#define p12 (RV2DRV(variables[12])->val)
#define p13 (RV2DRV(variables[13])->val)
#define p14 (RV2DRV(variables[14])->val)
#define p15 (RV2DRV(variables[15])->val)
#define p16 (RV2DRV(variables[16])->val)
#define p17 (RV2DRV(variables[17])->val)
#define p18 (RV2DRV(variables[18])->val)
#define p19 (RV2DRV(variables[19])->val)
#define p20 (RV2DRV(variables[20])->val)
#define p21 (RV2DRV(variables[21])->val)
#define p22 (RV2DRV(variables[22])->val)
#define p23 (RV2DRV(variables[23])->val)
#define p24 (RV2DRV(variables[24])->val)
#define p25 (RV2DRV(variables[25])->val)
#define p26 (RV2DRV(variables[26])->val)
#define p27 (RV2DRV(variables[27])->val)
#define p28 (RV2DRV(variables[28])->val)
#define p29 (RV2DRV(variables[29])->val)
#define p30 (RV2DRV(variables[30])->val)
#define p31 (RV2DRV(variables[31])->val)

// parent cardinality
#define cp0 (RV2DRV(variables[0])->cardinality)
#define cp1 (RV2DRV(variables[1])->cardinality)
#define cp2 (RV2DRV(variables[2])->cardinality)
#define cp3 (RV2DRV(variables[3])->cardinality)
#define cp4 (RV2DRV(variables[4])->cardinality)
#define cp5 (RV2DRV(variables[5])->cardinality)
#define cp6 (RV2DRV(variables[6])->cardinality)
#define cp7 (RV2DRV(variables[7])->cardinality)
#define cp8 (RV2DRV(variables[8])->cardinality)
#define cp9 (RV2DRV(variables[9])->cardinality)
#define cp10 (RV2DRV(variables[10])->cardinality)
#define cp11 (RV2DRV(variables[11])->cardinality)
#define cp12 (RV2DRV(variables[12])->cardinality)
#define cp13 (RV2DRV(variables[13])->cardinality)
#define cp14 (RV2DRV(variables[14])->cardinality)
#define cp15 (RV2DRV(variables[15])->cardinality)
#define cp16 (RV2DRV(variables[16])->cardinality)
#define cp17 (RV2DRV(variables[17])->cardinality)
#define cp18 (RV2DRV(variables[18])->cardinality)
#define cp19 (RV2DRV(variables[19])->cardinality)
#define cp20 (RV2DRV(variables[20])->cardinality)
#define cp21 (RV2DRV(variables[21])->cardinality)
#define cp22 (RV2DRV(variables[22])->cardinality)
#define cp23 (RV2DRV(variables[23])->cardinality)
#define cp24 (RV2DRV(variables[24])->cardinality)
#define cp25 (RV2DRV(variables[25])->cardinality)
#define cp26 (RV2DRV(variables[26])->cardinality)
#define cp27 (RV2DRV(variables[27])->cardinality)
#define cp28 (RV2DRV(variables[28])->cardinality)
#define cp29 (RV2DRV(variables[29])->cardinality)
#define cp30 (RV2DRV(variables[30])->cardinality)
#define cp31 (RV2DRV(variables[31])->cardinality)

// child cardinality

#define cc (RV2DRV(rv)->cardinality)

#endif 




