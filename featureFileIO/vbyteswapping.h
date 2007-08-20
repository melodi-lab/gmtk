/*
 * Vector copying, conversion, and byte swapping code for GMTK
 *
 *
 */

#ifndef VBYTESWAPPING_H
#define VBYTESWAPPING_H

#include "pfile.h"

// basic byte swapping (swapb) routines
extern intv_int32_t swapb_i32_i32(intv_int32_t val);
extern short swapb_short_short(short sval);
extern float swapb_f32_f32(float fval);

// static intv_int32_t copy_i32_i32(intv_int32_t from) { return from; }
extern void copy_i32_vi32(const size_t len, intv_int32_t from, intv_int32_t* to);
extern void swapb_vi32_vi32(const size_t len, const intv_int32_t* from, intv_int32_t* to);
extern void copy_vi32_vi32(const size_t len, const intv_int32_t* from, intv_int32_t* to);
extern void copy_i_vi(const size_t len, intv_int32_t from, intv_int32_t* to);


// static float copy_f32_f32(float from) { return from; }
extern void copy_f32_vf32(const size_t len, float from, float* to);
extern void swapb_vf32_vf32(const size_t len, const float* from, float* to);
extern void copy_vf32_vf32(const size_t len, const float* from, float* to);
extern void copy_f_vf(const size_t len, float from, float* to);


#if 0

//old routines

extern void copy_add_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to);
extern void copy_mul_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to);
extern void copy_sub_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to);
extern void copy_div_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to);

extern void swapb_add_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to);
extern void swapb_mul_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to);
extern void swapb_sub_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to);
extern void swapb_div_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to);

#endif

// 32-bit float versions of the above

extern void copy_add_vf32_vf32(size_t len, const float * from, float* to);
extern void copy_mul_vf32_vf32(size_t len, const float* from, float* to);
extern void copy_sub_vf32_vf32(size_t len, const float* from, float* to);
extern void copy_div_vf32_vf32(size_t len, const float* from, float* to);

extern void swapb_add_vf32_vf32(size_t len, const float* from, float* to);
extern void swapb_mul_vf32_vf32(size_t len, const float* from, float* to);
extern void swapb_sub_vf32_vf32(size_t len, const float* from, float* to);
extern void swapb_div_vf32_vf32(size_t len, const float* from, float* to);



#endif
