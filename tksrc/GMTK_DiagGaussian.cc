/*-
 * GMTK_DiagGaussian.cc
 *        Code for plain vanilla diagonal Gaussians.
 *
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ieeefp.h>
#include <float.h>
#include <assert.h>

#include <string>

#include "general.h"
VCID("$Header$");
#include "error.h"
#include "rand.h"


#include "GMTK_DiagGaussian.h"
#include "GMTK_GMParms.h"

#ifndef M_PI
#define M_PI               3.14159265358979323846  /* pi */
#endif


DiagGaussian::~DiagGaussian()
{
}


void
DiagGaussian::read(iDataStreamFile& is)
{
  // read name
  NamedObject::read(is);

  // read mean vector
  string str;
  is.read(str);
  



  // read covariance vector



}


void
DiagGaussian::write(oDataStreamFile& os)
{
  os.writeComment("burying values");
  os.nl();
  float *burValsp = burVals;
  for (int i=0;i<parent->nFeats;i++) {
    const int nComs_p1 = parent->numComs[i]+1;
    for (int j=0;j<nComs_p1;j++) {
      os.writeFloat(*burValsp++,"DiagGaussian::write bvs");
    }
    os.nl();
  }
  os.writeInt(0,"DiagGaussian::write cvtyp");
  os.writeComment("diagonal covariance type");
  os.nl();
  for (int i=0;i<parent->nFeats;i++) 
    os.writeFloat(variances[i],"DiagGaussian::write vars");
  os.nl();
}


void
DiagGaussian::preCompute()
{
  double det = 1.0;
  for (int i=0;i<parent->nFeats;i++) {
    if (variances[i] <= varianceFloor) {
      // Theoretically, this shouldn't happen unless you are reading
      // in a file that was computed from a previous run with a different threshold.
      coredump("DiagGaussian::preCompute, variance %d hit floor, cmp %d, parent id %d. Dumping.",i,parent->whichComponentAmI(this),parent->id());
    }
    variances_inv[i] = 1.0/variances[i];
    det *= variances[i];
  }
  if (det <= DBL_MIN) {
    coredump("DiagGaussian::preCompute, determinant hit minimum, cmp %d, parent id %d. Dumping.",parent->whichComponentAmI(this),parent->id());
  }
  const double tmp = (pow(2*M_PI,parent->nFeats/2.0)*sqrt(det));
  if (tmp <= DBL_MIN) {
    coredump("DiagGaussian::preCompute, norm const hit maximum, cmp %d, parent id %d. Dumping.",parent->whichComponentAmI(this),parent->id());
  }
  // log_inv_normConst = -log(tmp);
  log_inv_normConst = -0.5*(parent->nFeats*log(2*M_PI) + log(det));
}

void
DiagGaussian::randomize()
{
  assert ( bitmask & bm_basicAllocated );
  means -> makeRandom();
  variance -> makeRandom();
}



// This can be changed from 'float' to 'double' to
// provide extra range for temporary accumulators. Alternatively,
// decreasing the program's mixCoeffVanishRatio at the beginning
// of training should eliminate any component that produces
// such low scores.
#define TMP_ACCUMULATOR_TYPE double

//
// compute the log probability of x with stride 'stride'
// 
logpr
DiagGaussian::log_p(const float *const x,
		    const ptr32* const base,
		    const int stride)
{

  TMP_ACCUMULATOR_TYPE d=0.0;
  const float *xp = x;
  float *var_invp = variance->variances_inv;
  // we assume that parent->nFeats > 0
  const int nFeats = parent->nFeats;

  const float * burValsp = burVals;
  const int *numComsp = parent->numComs;
  const int *numComs_endp = parent->numComs + nFeats;
  const int *lagStrideOffsetsp = parent->lagStrideOffsets;

  do {
    float u=0.0;
    const int nComs = *numComsp++;

    if (nComs > 0) {
      const int *lagStrideOffsets_endp = lagStrideOffsetsp + nComs;
      do {
	u += (*burValsp) *
	  *(x + *lagStrideOffsetsp);
	lagStrideOffsetsp++;
	burValsp++;
      } while (lagStrideOffsetsp != lagStrideOffsets_endp);
    }

    // get the last Z=1 burying value.
    u += (*burValsp++);
    
    TMP_ACCUMULATOR_TYPE tmp = (u - *xp);
    d += tmp*tmp*(*var_invp);
    xp++;
    var_invp++;
  } while (numComsp != numComs_endp);
  d *= -0.5;
  return logpr(0,(log_inv_normConst + d));
}



//
// Compute the log probability of 'x' with stride 'stride'
// using comentaries 'z' of length parent->z_cache.len()
// We assume that parent->nFeats > 0
// 
logpr
DiagGaussian::log_p(const float *const x,const int stride,
		    const float *const z)
{

  assert ( parent->bitmask & GaussianMixture::bm_strideAllocated );

  const int nFeats = parent->nFeats;
  TMP_ACCUMULATOR_TYPE d=0.0;
  const float *xp = x;
  const float *x_endp = x + nFeats;
  const float *zp = z;
  float *var_invp = variances_inv;
  const float * burValsp = burVals;

  if (parent->z_cache.len() > 0) {
    // then this is a real BMM model with commentary variables.
    const int *numComsp = parent->numComs;
    do {
      float u=0.0;
      const int nComs = *numComsp++;
      if (nComs > 0) {
	const float * const z_endp = zp + nComs;
	do {
	  u += (*burValsp++)  * (*zp++);
	} while (zp != z_endp);
      }
      // get the last Z=1 burying value.
      u += (*burValsp++);
      const TMP_ACCUMULATOR_TYPE tmp = (u - *xp);
      d += tmp*tmp*(*var_invp);
      xp++;
      var_invp++;
    } while (xp != x_endp);
  } else {
    // then this is a degenerate HMM
    do {
      // the last Z=1 burying value only
      const float u = (*burValsp++);
      const TMP_ACCUMULATOR_TYPE tmp = (u - *xp);
      d += tmp*tmp*(*var_invp);
      xp++;
      var_invp++;
    } while (xp != x_endp);
  }

#if 0
#if (TMP_ACCUMULATOR_TYPE == double)
  // no need to do this if the type was "float"
  // as we should have died with an overflow exception by now.
  if (d > MAXFLOAT) {
    warning("Warning: very large d=%e value. Gaus. Mix num %d, cmp %d",d,
	    parent->id(),parent->whichComponentAmI(this));
  }
#endif
#endif

  d *= -0.5;
  return logpr(0,(log_inv_normConst + d));
}


/////////////////
// EM routines //
/////////////////

void
DiagGaussian::emInit()
{
  if (bitmask & bm_emAllocated)
    return;

  n_variances = new float[parent->nFeats];
  n_burVals = new float[parent->sum_nComsp1];

  gamma_oo = new float[parent->nFeats];
  gamma_oz = new float[parent->sum_nComsp1];
  gamma_zz = new float[parent->sum_nComsp1SqH];
  gamma_zz_expanded = new 
    float[2*parent->sum_nComsp1SqH - parent->sum_nComsp1] ;
  bitmask |= bm_emAllocated;
}


void
DiagGaussian::startEmIteration()
{
  // the accumulators    
  ::memset(gamma_oz,0,sizeof(gamma_oz[0])*parent->sum_nComsp1);
  ::memset(gamma_zz,0,sizeof(gamma_zz[0])*parent->sum_nComsp1SqH);
  ::memset(gamma_oo,0,sizeof(gamma_oo[0])*parent->nFeats);
  preCompute();
}


void
DiagGaussian::emSwapCurAndNew()
{
  assert ( bitmask & 
	   (bm_basicAllocated | bm_emAllocated) );

  if (bitmask & bm_swapped)
    return;
  else 
    bitmask |= bm_swapped;

  genSwap(burVals,n_burVals);
  genSwap(variances,n_variances);

}


void
DiagGaussian::emIncrement(const float prob,
			   const float *const oz_array,
			   const float *const zz_array,
			   const float *const oo_array)
{
  assert (bitmask & bm_emAllocated);

  const float *srcp;
  const float *src_endp; 
  float *dstp;

  // do three saxpy's 

  dstp = gamma_oo;
  srcp = oo_array;
  src_endp = oo_array + parent->nFeats;
  do {
    *dstp++ += prob*(*srcp++);
  } while (srcp != src_endp);

  dstp = gamma_oz;
  srcp = oz_array;
  src_endp = oz_array + parent->sum_nComsp1;
  do {
    *dstp++ += prob*(*srcp++);
  } while (srcp != src_endp);


  dstp = gamma_zz;
  srcp = zz_array;
  src_endp = zz_array + parent->sum_nComsp1SqH;
  do {
    *dstp++ += prob*(*srcp++);
  } while (srcp != src_endp);

}


void
DiagGaussian::emEndIteration(logpr cmpSop_acc)
{
  assert (bitmask & bm_emAllocated);
  if (cmpSop_acc.zero())
  {
    printf("Error: zero cmpSop_acc in DiagGaussian::endEmIteration\n");
    return;
  }

  const double cmpSop_acc_inv = cmpSop_acc.inverse().unlog();

  int flooredVariances=0;

  expandGammaZZ();
  // the accumulators
  float *a_oz = gamma_oz;
  float *a_zz = gamma_zz_expanded;

  float *n_burValsp = n_burVals;
  for (int feat=0;feat<parent->nFeats;feat++) {
    const int nComs_p1 = parent->numComs[feat]+1;

    // Solve for the burying coefficients.
    // Solves Ax = b for x where
    // A is nXn, x is nX1, and b is nX1
    // here,
    //     A = a_zz,
    //     x = the output
    //     b = a_oz (which gets destroyed)
    // First copy over a_oz since we need it later

    ::memcpy(n_burValsp,a_oz,sizeof(float)*nComs_p1);
    // do the multiply, putting the real burying results in bv
    ::lineqsolve(nComs_p1,1,
		 a_zz,n_burValsp);
    
    // now solve for the variances
    float tmp=0;
    for (int i=0;i<nComs_p1;i++) {
      tmp += ( n_burValsp[i] *a_oz[i]);
    }
    
    n_variances[feat] = (gamma_oo[feat] - tmp)*cmpSop_acc_inv;

    /////////////////////////////////////////////////////
    // When variances hit zero:
    // There could be several reasons for the variances hitting the floor:
    //   1) BMM prediction of means is very good which leads to low
    //      variances.  In this case, we shouldn't drop the component,
    //      instead we should just hard-limit the variance (i.e., here
    //      mixCoeffs[this] is not too small).  
    //   2) Very small quantity of training data (i.e., mixCoeffs[this] is
    //      very small). In this case we should remove the component
    //      completely.
    //   3) If only one mixture is left and this happens, then it could be
    //      that the sop is small. In this case, there's probably a problem
    //      with the HMM topology. I.e., we should remove the state. For now,
    //      however, if this happens, the variance will be floored like
    //      in case 1.
    if (n_variances[feat] < varianceFloor) {
      // Don't let variances go less than variance floor. We
      // assume this is happening becuase the burying links
      // are doing such a good job predicting the mean that
      // it is making the variances very small.
      flooredVariances++;

      // Could either keep old variance or set to varianceFloor.

      // keep old variance
      n_variances[feat] = variances[feat];

      // set to variance floor. 
      // n_variances[feat] = varianceFloor;

    }

    a_oz += nComs_p1;
    n_burValsp += nComs_p1;
    a_zz += nComs_p1*nComs_p1;
  }
  if (flooredVariances > 0)
    fprintf(stderr,"WARNING: n_variance(s) hit floor for %d/%d feats, Gaus. Mix num %d, cmp %d, adjusting...\n",flooredVariances,parent->nFeats,
	    parent->id(),parent->whichComponentAmI(this));

}


void
DiagGaussian::emStoreAccumulators(oDataStreamFile& ofile)
{
  int i;

  ofile.writeComment("DiagGaussian oz values");
  ofile.nl();
  for (i=0;i<parent->sum_nComsp1;i++) {
    ofile.write(gamma_oz[i]);
  }
  ofile.nl();

  ofile.writeComment("DiagGaussian oo values");
  ofile.nl();
  for (i=0;i<parent->nFeats;i++) {
    ofile.write(gamma_oo[i]);
  }
  ofile.nl();

  ofile.writeComment("DiagGaussian zz values");
  ofile.nl();
  for (i=0;i<parent->sum_nComsp1SqH;i++) {
    ofile.write(gamma_zz[i]);
  }
  ofile.nl();
}

void
DiagGaussian::emLoadAccumulators(iDataStreamFile& ifile)
{
  int i;

  for (i=0;i<parent->sum_nComsp1;i++) {
    ifile.read(gamma_oz[i]);
  }

  for (i=0;i<parent->nFeats;i++) {
    ifile.read(gamma_oo[i]);
  }

  for (i=0;i<parent->sum_nComsp1SqH;i++) {
    ifile.read(gamma_zz[i]);
  }
}


void
DiagGaussian::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  int i;
  float f;

  for (i=0;i<parent->sum_nComsp1;i++) {
    ifile.read(f);
    gamma_oz[i] += f;
  }

  for (i=0;i<parent->nFeats;i++) {
    ifile.read(f);
    gamma_oo[i] += f;
  }

  for (i=0;i<parent->sum_nComsp1SqH;i++) {
    ifile.read(f);
    gamma_zz[i] += f;
  }
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////


void DiagGaussian::sampleGenerate(float *const sample,
				  const int numBefore,
				  const sArray<float>& means)
{
  
  // first set variables to N(0,sigma)
  for (int feat=0;feat<parent->nFeats;feat++) {
    sample[feat] = sqrt(variances[feat])*
      ::inverse_normal_func(rnd.drand48pe());
  }

  // compute and add the means.
  float *burValsp = burVals;
  for (int feat=0;feat<parent->nFeats;feat++) {
    float u=0.0;

    GaussianMixture::bindex* bIndicesp = 
      (parent->bIndices[feat]);
    for (int i=0;i<parent->numComs[feat];i++) {
      float burCoefficient = *burValsp++;
      float burFeature;
      const int stride = means.len();
      
      assert (bIndicesp->lag < 0);
      assert (bIndicesp->offset < means.len());

      if (-bIndicesp->lag > numBefore)
	burFeature = means[feat];
      else 
        burFeature =
	  *(sample + bIndicesp->lag*stride + bIndicesp->offset);

      u += burFeature*burCoefficient;
      bIndicesp++;
    }
    // get the last Z=1 burying value.
    u += *burValsp++;

    sample[feat] += u;

  }
}









