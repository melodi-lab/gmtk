#ifndef ENDIAN_H_INCLUDED
#define ENDIAN_H_INCLUDED

typedef int intv_int32_t;
typedef unsigned int intv_uint32_t;

intv_int32_t
swapb_i32_i32(intv_int32_t val)
{
  intv_uint32_t uval;
  intv_uint32_t res;
  intv_int32_t b0, b1, b2, b3;

  uval = (intv_uint32_t) val;
  b0 = uval >> 24;
  b1 = (uval >> 8) & 0x0000ff00;
  b2 = (uval << 8) & 0x00ff0000;
  b3 = uval << 24;

  res = b0 | b1 | b2 | b3;

  return (intv_int32_t) res;
}

intv_int32_t
copy_i32_i32(intv_int32_t from) {
  return from;
}

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
copy_f_vf(size_t len, float from, float* to)
{
  size_t i;

  for (i=0; i<len; i++)
    *to++ = from;
}




#endif




