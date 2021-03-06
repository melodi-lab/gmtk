
/*
 *
 * Copyright (C) 2007 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

/* The ISO C99 standard specifies that in C++ implementations some
   macros should only be defined if explicitly requested.  */
#define __STDC_LIMIT_MACROS 1
#define __STDC_CONSTANT_MACROS 1
   // The ISO C99 standard specifies that the macros in inttypes.h must
   //  only be defined if explicitly requested. 
#define __STDC_FORMAT_MACROS 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>
#endif
#if HAVE_STDINT_H
#  include <stdint.h>
#endif

#include <assert.h>
#include <cctype>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>


#include "vbyteswapping.h"

// Byte swapping/copying routines. Also included are a few that
// perform the four operations, addition, substraction, multiplication
// and division on 32-bit floating-point types.  In the case of
// division, no warning is given if the denominator is zero; instead
// the division is not performed.  The routines are used for feature
// stream combination in the GMTK Observation Matrix library. 
// Much of the code was taken originally from Quicknet's fltvec and intvec
// library, but only those routines that were needed here were used. The
// swap and copy routines were written by Karim Filali (karim@cs.washington.edu).
// This was then changed to a separate file, and then updated to work with
// C++ type-punning and aliasing rules by Bilmes.


// 32-bit unions for byte swapping 32-bit ints and 32-bit floats.
union union_float_int32 {
  intv_int32_t i;
  float f;
};

union union_double_int64 {
  intv_int64_t i;
  double       f;
};


// 32-bit unsigned Endian byte swap
intv_int32_t swapb_i32_i32(intv_int32_t val)
{
  REGISTER intv_uint32_t uval;
  REGISTER intv_uint32_t res;
  REGISTER intv_int32_t b0, b1, b2, b3;

  uval = (intv_uint32_t) val;
  b0 = uval >> 24;
  b1 = (uval >> 8) & 0x0000ff00;
  b2 = (uval << 8) & 0x00ff0000;
  b3 = uval << 24;

  res = b0 | b1 | b2 | b3;
  return (intv_int32_t) res;
}


// 64-bit unsigned Endian byte swap
intv_int64_t swapb_i64_i64(intv_int64_t val)
{
  REGISTER intv_uint64_t uval;
  REGISTER intv_uint64_t res;
  REGISTER intv_int64_t b0, b1, b2, b3, b4, b5, b6, b7;

  uval = (intv_uint64_t) val;
  b0 =  uval >> 56;
  b1 = (uval >> 40) & INT64_C(0x000000000000ff00);
  b2 = (uval >> 24) & INT64_C(0x0000000000ff0000);
  b3 = (uval >>  8) & INT64_C(0x00000000ff000000);
  b4 = (uval <<  8) & INT64_C(0x000000ff00000000);
  b5 = (uval << 24) & INT64_C(0x0000ff0000000000);
  b6 = (uval << 40) & INT64_C(0x00ff000000000000);
  b7 =  uval << 56;

  res = b0 | b1 | b2 | b3 | b4 | b5 | b6 | b7;
  return (intv_int32_t) res;
}


double swapb_f64_f64(double fval) {
  union_double_int64 diu;

  diu.f = fval;
  diu.i = swapb_i64_i64(diu.i);
  return diu.f;
}


// TODO: this should be changed to i16 at some point.
short swapb_short_short(short sval) {
 
  short res;
  short s0,s1;

  unsigned short usval = (unsigned short) sval;  

  s0 = usval >> 8;
  s1 = (usval << 8) & 0x0000ff00;
  res = s0 | s1 ;
  
  return res;
}

// 32-bit float byte swap
float
swapb_f32_f32(const float fval)
{
  union_float_int32 fiu;

  fiu.f = fval;
  fiu.i = swapb_i32_i32(fiu.i);
  return fiu.f;
}

/////////////////////////////////////////////////////////////////////////////////
	
void
copy_i32_vi32(size_t len, intv_int32_t from, intv_int32_t* to)
{
  size_t i;
  for (i=0; i<len; i++)
    *to++ = from;
}

void
swapb_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
  for (i=0; i<len; i++)
    *to++ = swapb_i32_i32(*from++);
}

void
copy_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
  for (i=0; i<len; i++)
    *to++ = *from++;
}

void
copy_i_vi(const size_t len, intv_int32_t from, intv_int32_t* to)
{
  size_t i;
  for (i=0; i<len; i++)
    *to++ = from;
}

/////////////////////////////////////////////////////////////////////////////////
// 32-bit float versions of the above
	
void
copy_f32_vf32(const size_t len, float from, float* to)
{
  size_t i;
  for (i=0; i<len; i++)
    *to++ = from;
}

void
swapb_vf32_vf32(const size_t len, const float* from, float* to)
{
  size_t i;
  for (i=0; i<len; i++)
    *to++ = swapb_f32_f32(*from++);
}

void
copy_vf32_vf32(const size_t len, const float* from, float* to)
{
  size_t i;
  for (i=0; i<len; i++)
    *to++ = *from++;
}

void
copy_f_vf(const size_t len, float from, float* to)
{
  size_t i;
  for (i=0; i<len; i++)
    *to++ = from;
}


/////////////////////////////////////////////////////////////////////////////////

#if 0

//old routines

void copy_add_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
  float * flto=(float*) to;
  float * flfrom=(float*) from;
  for (i=0; i<len; i++) {
    *flto += (*flfrom++);
    ++flto;
  }
  
}


void copy_mul_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
  float * flto=(float*) to;
  float * flfrom=(float*) from;
  for (i=0; i<len; i++) {
    *flto *= (*flfrom++);
    ++flto;
  }

}


void copy_sub_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;

  float * flto=(float*) to;
  float * flfrom=(float*) from;
  for (i=0; i<len; i++) {
    *flto -= (*flfrom++);
    ++flto;
  }

}


void copy_div_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
  float * flto=(float*) to;
  float * flfrom=(float*) from;
  for (i=0; i<len; i++) {
    if(*flfrom !=0) {
      *flto /= (*flfrom);
    }
    ++flto;++flfrom;
  }
}


void swapb_add_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
  float * float_to=(float*) to;
  intv_int32_t swapped_from;
  float* float_from;
  for (i=0; i<len; i++) {
    swapped_from=swapb_i32_i32(*from++);
    float_from =  (float*)&swapped_from;
    *float_to += (float)*float_from;
    ++float_to;
  }
}


void swapb_mul_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
  float * float_to=(float*) to;
  intv_int32_t swapped_from;
  float* float_from;
  for (i=0; i<len; i++) {
    swapped_from=swapb_i32_i32(*from++);
    float_from =  (float*)&swapped_from;
    *float_to *= (float)*float_from;
    ++float_to;
  }

}


void swapb_sub_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
  float * float_to=(float*) to;
  intv_int32_t swapped_from;
  float* float_from;
  for (i=0; i<len; i++) {
    swapped_from=swapb_i32_i32(*from++);
    float_from =  (float*)&swapped_from;
    *float_to -= (float)*float_from;
    ++float_to;
  }

}


void swapb_div_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
  float * float_to=(float*) to;
  intv_int32_t swapped_from;
  float* float_from;
  for (i=0; i<len; i++) {
    swapped_from=swapb_i32_i32(*from++);
    float_from =  (float*)&swapped_from;
    if(float_from!=0)
      *float_to /= (float)*float_from;
    ++float_to;
  }

}


#endif

/////////////////////////////////////////////////////////////////////



void copy_add_vf32_vf32(size_t len, const float* from, float* to)
{
  if (len == 0) return;
  const float* to_endp = to + len;
  while (to != to_endp) {
    (*to++) += (*from++);
  }
}


void copy_mul_vf32_vf32(size_t len, const float* from, float* to)
{
  if (len == 0) return;
  const float* to_endp = to + len;
  while (to != to_endp) {
    (*to++) *= (*from++);
  }
}


void copy_sub_vf32_vf32(size_t len, const float* from, float* to)
{
  if (len == 0) return;
  const float* to_endp = to + len;
  while (to != to_endp) {
    (*to++) -= (*from++);
  }
}


void copy_div_vf32_vf32(size_t len, const float* from, float* to)
{
  if (len == 0) return;
  const float* to_endp = to + len;
  while (to != to_endp) {
    if(*from != 0.0f) {
      *to /= (*from);
    }
    to++;
    from++;
  }
}


void swapb_add_vf32_vf32(size_t len, const float* from, float* to)
{
  if (len == 0) return;
  const float* to_endp = to + len;
  while (to != to_endp) {
    (*to++) += swapb_f32_f32(*from++);
  }
}


void swapb_mul_vf32_vf32(size_t len, const float* from, float* to)
{
  if (len == 0) return;
  const float* to_endp = to + len;
  while (to != to_endp) {
    (*to++) *= swapb_f32_f32(*from++);
  }
}


void swapb_sub_vf32_vf32(size_t len, const float* from, float* to)
{
  if (len == 0) return;
  const float* to_endp = to + len;
  while (to != to_endp) {
    (*to++) -= swapb_f32_f32(*from++);
  }
}


void swapb_div_vf32_vf32(size_t len, const float* from, float* to)
{
  if (len == 0) return;
  const float* to_endp = to + len;
  while (to != to_endp) {
    if(*from != 0.0f) {
      *to /= swapb_f32_f32(*from);
    }
    to++;
    from++;
  }
}

