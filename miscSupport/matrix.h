//
// Various matrix related routines.
//
// $Header$
//
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu


#ifndef MATRIX_H
#define MATRIX_H

#include "general.h"


////////////////////////////////////////////////////
// MATRIX TRANSPOSE: mTranspose
//   Transposes the m X n matrix 'in' and
//   places the result into 'out'.
// ASSUMPTIONS:
//     m > 0 
//     n > 0
// If these are not met, unpredictable results will occur.
////////////////////////////////////////////////////
template <class inType, class outType>
void
mTranspose(const inType *const in,
	   const int m, // rows of 'in'
	   const int n, // cols of 'in'
	   outType *const out)
{
  assert( m > 0 );
  assert( n > 0 );
  const inType *inp = in;
  outType *outp = out;
  int i;
  i = 0;
  do {
    outType *outpp = outp;
    outType *const outp_endp = outp + m;
    const inType *inpp = inp;
    do {
      *outpp = *inpp;
      inpp += n;
      outpp++;
    } while (outpp != outp_endp);
    inp++;
    outp = outp_endp;
    i++;
  } while (i<n);
}

////////////////////////////////////////////////////
// matrixSelfOuterProduct: compute Z = VV'
// 
// Note that result is a symetric part, so we only compute
// the upper triangular portion.
// Assume that:
//     n > 0
//     v is nX1
//     z is upper triangular part of an nXn matrix
//       (i.e., z has only n(n+1)/2 entries.
////////////////////////////////////////////////////
template <class mType>
void matrixSelfOuterProduct(const mType *const v,
			    const int n,
			    mType *const z)
{
  switch (n) {
  case 1:
    {
      mType t0 = v[0];
      z[0] = t0*t0;
    }
    break;
  case 2:
    {
      mType t0,t1;
      t0 = v[0];
      t1 = v[1];
      z[0] = t0*t0; 
      z[1] = t0*t1;
      z[2] = t1*t1;
    }
    break;
  case 3:
    {
      mType t0,t1,t2;
      t0 = v[0];
      t1 = v[1];
      t2 = v[2];
      z[0] = t0*t0;
      z[1] = t0*t1;
      z[2] = t0*t2;
      z[3] = t1*t1;
      z[4] = t1*t2;
      z[5] = t2*t2;
    }
    break;
  case 4:
    {
      mType t0,t1,t2,t3;
      t0 = v[0];
      t1 = v[1];
      t2 = v[2];
      t3 = v[3];
      z[0] = t0*t0;
      z[1] = t0*t1;
      z[2] = t0*t2;
      z[3] = t0*t3;
      z[4] = t1*t1;
      z[5] = t1*t2;
      z[6] = t1*t3;
      z[7] = t2*t2;
      z[8] = t2*t3;
      z[9] = t3*t3;
    }
    break;
  case 5:
    {
      mType t0,t1,t2,t3,t4;
      t0 = v[0];
      t1 = v[1];
      t2 = v[2];
      t3 = v[3];
      t4 = v[4];
      z[0] = t0*t0;
      z[1] = t0*t1;
      z[2] = t0*t2;
      z[3] = t0*t3;
      z[4] = t0*t4;
      z[5] = t1*t1;
      z[6] = t1*t2;
      z[7] = t1*t3;
      z[8] = t1*t4;
      z[9] = t2*t2;
      z[10] = t2*t3;
      z[11] = t2*t4;
      z[12] = t3*t3;
      z[13] = t3*t4;
      z[14] = t4*t4;
    }
    break;
  case 6:
    {
      mType t0,t1,t2,t3,t4,t5;
      t0 = v[0];
      t1 = v[1];
      t2 = v[2];
      t3 = v[3];
      t4 = v[4];
      t5 = v[5];
      z[0] = t0*t0;
      z[1] = t0*t1;
      z[2] = t0*t2;
      z[3] = t0*t3;
      z[4] = t0*t4;
      z[5] = t0*t5;
      z[6] = t1*t1;
      z[7] = t1*t2;
      z[8] = t1*t3;
      z[9] = t1*t4;
      z[10] = t1*t5;
      z[11] = t2*t2;
      z[12] = t2*t3;
      z[13] = t2*t4;
      z[14] = t2*t5;
      z[15] = t3*t3;
      z[16] = t3*t4;
      z[17] = t3*t5;
      z[18] = t4*t4;
      z[19] = t4*t5;
      z[20] = t5*t5;
    }
    break;
  case 7:
    {
      mType t0,t1,t2,t3,t4,t5,t6;
      t0 = v[0];
      t1 = v[1];
      t2 = v[2];
      t3 = v[3];
      t4 = v[4];
      t5 = v[5];
      t6 = v[6];
      z[0] = t0*t0;
      z[1] = t0*t1;
      z[2] = t0*t2;
      z[3] = t0*t3;
      z[4] = t0*t4;
      z[5] = t0*t5;
      z[6] = t0*t6;
      z[7] = t1*t1;
      z[8] = t1*t2;
      z[9] = t1*t3;
      z[10] = t1*t4;
      z[11] = t1*t5;
      z[12] = t1*t6;
      z[13] = t2*t2;
      z[14] = t2*t3;
      z[15] = t2*t4;
      z[16] = t2*t5;
      z[17] = t2*t6;
      z[18] = t3*t3;
      z[19] = t3*t4;
      z[20] = t3*t5;
      z[21] = t3*t6;
      z[22] = t4*t4;
      z[23] = t4*t5;
      z[24] = t4*t6;
      z[25] = t5*t5;
      z[26] = t5*t6;
      z[27] = t6*t6;
    }
    break;
  case 8:
    {
      mType t0,t1,t2,t3,t4,t5,t6,t7;
      t0 = v[0];
      t1 = v[1];
      t2 = v[2];
      t3 = v[3];
      t4 = v[4];
      t5 = v[5];
      t6 = v[6];
      t7 = v[7];
      z[0] = t0*t0;
      z[1] = t0*t1;
      z[2] = t0*t2;
      z[3] = t0*t3;
      z[4] = t0*t4;
      z[5] = t0*t5;
      z[6] = t0*t6;
      z[7] = t0*t7;
      z[8] = t1*t1;
      z[9] = t1*t2;
      z[10] = t1*t3;
      z[11] = t1*t4;
      z[12] = t1*t5;
      z[13] = t1*t6;
      z[14] = t1*t7;
      z[15] = t2*t2;
      z[16] = t2*t3;
      z[17] = t2*t4;
      z[18] = t2*t5;
      z[19] = t2*t6;
      z[20] = t2*t7;
      z[21] = t3*t3;
      z[22] = t3*t4;
      z[23] = t3*t5;
      z[24] = t3*t6;
      z[25] = t3*t7;
      z[26] = t4*t4;
      z[27] = t4*t5;
      z[28] = t4*t6;
      z[29] = t4*t7;
      z[30] = t5*t5;
      z[31] = t5*t6;
      z[32] = t5*t7;
      z[33] = t6*t6;
      z[34] = t6*t7;
      z[35] = t7*t7;
    }
    break;
  default:
    {
      const mType *v_rp=v; // v's row pointer
      const mType *const v_endp = &v[n]; // v end pointer
      mType *zp = z;
      do {
	const mType rv = *v_rp; // row value
	const mType *v_cp=v_rp; // v's column pointer
	do {
	  *zp++ = rv*(*v_cp);
	  v_cp++;
	} while (v_cp != v_endp);
	v_rp ++;
      } while (v_rp != v_endp);
    }
    break;
  }
}

#endif

