#include <iostream>
#include "mixNormal.h"

#include "matrix-ops.h"
#include "rand.h"

#ifndef M_2PI
#define M_2PI 6.28318530717958647692
#endif

//////////////////// sampleUniform ////////////////////

/** 
 *  sample from a uniform distribution between [a,b] 
 */
unsigned MixNormal::sampleUniform(unsigned a, unsigned b) const {
  return (unsigned)(a + rnd.drand48()*(b-a+1));
}

//////////////////// invertCov ////////////////////

/** 
 * returns the inverse of positive definite matrix and invDet is the
 * determinant of the inverse.  Check that this is indeed faster than
 * computing the inverse directly.  cov can be the same as invCov
 */
PARAM_DATA_TYPE MixNormal::invertCov(PARAM_DATA_TYPE* cov, PARAM_DATA_TYPE* invCov, unsigned d) const {

  PARAM_DATA_TYPE vars[MAX_POINTER_SET_DIM];

  diagCholeskyDecomp(cov,_tmpM,vars,d); // cov = _tmpM'* vars * _tmpM, 
  
  inverse(_tmpM,_tmpM,d); //_tmpM = _tmpM.inverse();       
  
  PARAM_DATA_TYPE invDet = 1.0 / prod(vars,d); //invDet = 1.0 / Vars.prod();
  
  for( unsigned i = 0; i < d*d; i++ ) 
    *(invCov + i) = 0.0;

  // The following loop is equivalent to the operations:
  // for ( unsigned i = 0; i < d; i++) S[i][i] = 1.0 / vars[i];
  // outM =  b * S * b.transpose();     // invCov = b^-1 * S^-1 * b'^-1
  for ( unsigned i = 0; i < d; i++) 
      for ( unsigned j = 0; j < d; j++) 
	for ( unsigned k = 0; k< d; k++) 
	  *(invCov + i*d + j) += *(_tmpM + i*d + k) / vars[k] *  *(_tmpM + j*d + k);

  return invDet;
}


//////////////////// getNumberofMixtures ////////////////////

unsigned MixNormal::getNumberofMixtures() const {
  return _numMixtures;
}

//////////////////// marginalize ////////////////////

/** 
 * Find the marginal distribution of the first dX components of 
 * vector [X' Y']' modeled by a mixture of Gaussians
 */
void MixNormal::marginalize(unsigned dX) {

 unsigned dY = _numVariables - dX;
 unsigned dim = _numVariables;

 if( ( dX < 1 ) || ( dX > _numVariables-1 ) )
   error("ERROR: the reduced dimension is less or more");

 for(unsigned l = 0; l < _numMixtures; l++){

   if (_fullCoVar) {
     // In this instance, _invCovsXY will contain the covariance matrix not its inverse despite its name.  It will get inverted later.

     LDU2moment(_invCovsXY+l*dim*dim, (_b + l*dim*dim), (_invVars + l*dim), dim);     //_invCovsXY = LDU2moment(_b[l], _invVars[l]);

   }
   else {
     for(unsigned i = 0; i < dim; i++)
       //*(_tmpM+i*dim+i) = 1.0 /(PARAM_DATA_TYPE) *(_invVars+l*dim + i);  // fill the diagonal of _tmpM
       //*(_invCovsXY+l*dim*dim+i) = 1.0 /(PARAM_DATA_TYPE) *(_invVars+l*dim + i);
       // we are not inverting _invCovsXY so store _invVars directly
       *(_invCovsXY+l*dim*dim+i) = *(_invVars+l*dim + i);
   }

   // Might also wanna remove this function altogether.  It is easy to access blocks of the full covariance
   //partitionCov(dX, _tmpM, _invCovsX+l*dX*dX, _invCovsY+l*dY*dY,dim);
   partitionCov(dX, _invCovsXY+l*dim*dim, _invCovsX+l*dX*dX, _invCovsY+l*dY*dY,dim);

   // now this is done only when we have full covariance
   if (_fullCoVar) _invDets[l] = invertCov(_invCovsXY+l*dim*dim, _invCovsXY+l*dim*dim, dim);

   for(unsigned i = 0; i < dX; i++)
     *(_meansX +l*dX +i) = *(_means+l*dim+i); 


   DBGFPRINTF((stderr,"In MixNormal::marginalize, calling moment2LDU.\n"));
   moment2LDU(_invCovsX+l*dX*dX,_bX+l*dX*dX,_invVarsX+l*dX,dX);


   // do we invert here even if we are using a diagonal covariance matrix?
   _invDetsX[l] = invertCov(_invCovsX+l*dX*dX, _invCovsX+l*dX*dX, dX);    //_invCovsX[l] = invertCov(_invCovsX[l], _invDetsX[l]);

   for(unsigned i = 0; i < dY; i++)
     *(_meansY +l*dY +i) = *(_means+l*dim+i+dX); 

   moment2LDU(_invCovsY+l*dY*dY,_bY+l*dY*dY,_invVarsY+l*dY,dY);

   _invDetsY[l] = invertCov(_invCovsY+l*dY*dY, _invCovsY+l*dY*dY,dY);    //_invCovsY[l] = invertCov(_invCovsY[l], _invDetsY[l]);
   


 }

 _inv_pow_sqrt_2piX = 1.0 / (double)pow(sqrt(M_2PI), (int)dX);
 _inv_pow_sqrt_2piY = 1.0 / (double)pow(sqrt(M_2PI), (int)dY);
}



//////////////////// marginalizeWithFirstParentOut ////////////////////

/** 
 * Find the marginal distribution of the first dX components of 
 * vector [X' Y']' modeled by a mixture of Gaussians
 */
void MixNormal::marginalizeWithFirstParentOut(unsigned dX) {

 unsigned dY = _numVariables - dX-1;
 unsigned dim = _numVariables;

 if( ( dX < 1 ) || ( dX > _numVariables-1 ) )
   error("ERROR: the reduced dimension is less or more");

 for(unsigned l = 0; l < _numMixtures; l++){

   if (_fullCoVar) {
     // In this instance, _invCovsXY will contain the covariance matrix not its inverse despite its name.  It will get inverted later.

     LDU2moment(_invCovsXY+l*dim*dim, (_b + l*dim*dim), (_invVars + l*dim), dim);     //_invCovsXY = LDU2moment(_b[l], _invVars[l]);

   }
   else {
     for(unsigned i = 0; i < dim; i++)
       //*(_tmpM+i*dim+i) = 1.0 /(PARAM_DATA_TYPE) *(_invVars+l*dim + i);  // fill the diagonal of _tmpM
       //*(_invCovsXY+l*dim*dim+i) = 1.0 /(PARAM_DATA_TYPE) *(_invVars+l*dim + i);
       // we are not inverting _invCovsXY so store _invVars directly
       *(_invCovsXY+l*dim*dim+i) = *(_invVars+l*dim + i);
   }

   // Might also wanna remove this function altogether.  It is easy to access blocks of the full covariance
   //partitionCov(dX, _tmpM, _invCovsX+l*dX*dX, _invCovsY+l*dY*dY,dim);
   partitionCov(dX, _invCovsXY+l*dim*dim, _invCovsX+l*dX*dX, _invCovsY+l*dY*dY,dim);

   // now this is done only when we have full covariance
   if (_fullCoVar) _invDets[l] = invertCov(_invCovsXY+l*dim*dim, _invCovsXY+l*dim*dim, dim);

   for(unsigned i = 0; i < dX; i++)
     *(_meansX +l*dX +i) = *(_means+l*dim+i); 


   DBGFPRINTF((stderr,"In MixNormal::marginalize, calling moment2LDU.\n"));
   moment2LDU(_invCovsX+l*dX*dX,_bX+l*dX*dX,_invVarsX+l*dX,dX);


   // do we invert here even if we are using a diagonal covariance matrix?
   _invDetsX[l] = invertCov(_invCovsX+l*dX*dX, _invCovsX+l*dX*dX, dX);    //_invCovsX[l] = invertCov(_invCovsX[l], _invDetsX[l]);

   for(unsigned i = 0; i < dY; i++)
     *(_meansY +l*dY +i) = *(_means+l*dim+i+dX+1); 

   moment2LDU(_invCovsY+l*dY*dY,_bY+l*dY*dY,_invVarsY+l*dY,dY);

   _invDetsY[l] = invertCov(_invCovsY+l*dY*dY, _invCovsY+l*dY*dY,dY);    //_invCovsY[l] = invertCov(_invCovsY[l], _invDetsY[l]);
   


 }

 _inv_pow_sqrt_2piX = 1.0 / (double)pow(sqrt(M_2PI), (int)dX);
 _inv_pow_sqrt_2piY = 1.0 / (double)pow(sqrt(M_2PI), (int)dY);
}



//////////////////// prob_x_Gaussian ////////////////////

/**
 * Compute the pdf a Gaussian at x by natural parameters
 * @return invSqrtTwoPiD * sqrt(invDet) * exp(-.5 * (invCov*(x-mean)).dot(x-mean)); 
*/
PARAM_DATA_TYPE MixNormal::prob_x_Gaussian(const PointerSetToDataPoints& ps,
					   const unsigned sampleNum,
					   const PARAM_DATA_TYPE* mean, 
					   const PARAM_DATA_TYPE* invCov, 
					   const PARAM_DATA_TYPE invDet,
					   const double invSqrtTwoPiD) const {
  unsigned dim = _numVariables;
  PARAM_DATA_TYPE sum=0.0;

  //  cout<<"vec = "<<*(ps.start[0]+sampleNum*ps.skip)<<" mean = "<<*mean<<" invCov = "<<*invCov<<" invDet = "<<invDet<<" invSqrtTwoPiD = "<<invSqrtTwoPiD<<endl;

  for(unsigned i=0; i<dim;++i) {
    for(unsigned j=0; j<dim;++j) {
      sum += *(invCov + j*dim +i) *  
	( *(ps.start[j]+sampleNum*ps.skip) - *(mean + j) ) * 
	( *(ps.start[i]+sampleNum*ps.skip)  - *(mean + i) );
    }
  }
  return (PARAM_DATA_TYPE) invSqrtTwoPiD * sqrt(invDet) * exp(-.5 * sum);
}

//////////////////// prob_x_Gaussian ////////////////////

/**
 * Compute the pdf a Gaussian at x by natural parameters
 */
PARAM_DATA_TYPE MixNormal::prob_x_Gaussian(const PointerSetToDataPoints& ps,
				const unsigned sampleNum,
				const unsigned start, // start of subVector.  Used only for ps.  Might change that and get rid of cov partioning
				const unsigned end,   // end of subVector 
				const PARAM_DATA_TYPE* mean, 
				const PARAM_DATA_TYPE* invCov, 
				const PARAM_DATA_TYPE invDet,
				const double invSqrtTwoPiD) const {
  unsigned vecDim = end - start + 1;
  PARAM_DATA_TYPE sum=0.0;

  for(unsigned i=0; i< vecDim;++i) {
      for(unsigned j=0; j< vecDim;++j) {
	sum += *(invCov + j*vecDim +i) *  
	  ( *(ps.start[j+start]+sampleNum*ps.skip) - *(mean + j) ) * 
	  ( *(ps.start[i+start]+sampleNum*ps.skip)  - *(mean + i) );
      }
  }
  return (PARAM_DATA_TYPE) invSqrtTwoPiD * sqrt(invDet) * exp(-.5 * sum);
}


//////////////////// prob_x_Gaussian ////////////////////

/**
 * Compute the pdf a Gaussian at vec by natural parameters
 */
PARAM_DATA_TYPE MixNormal::prob_x_Gaussian(const PARAM_DATA_TYPE *vec,
					   const unsigned vecSize,
					   const PARAM_DATA_TYPE* mean, 
					   const PARAM_DATA_TYPE* invCov, 
					   const PARAM_DATA_TYPE invDet,
					   const double invSqrtTwoPiD) const {
  PARAM_DATA_TYPE sum=0.0;

  for(unsigned i=0; i<vecSize;++i) {
      for(unsigned j=0; j<vecSize;++j) {
	sum += *(invCov + j*vecSize +i) *  ( *(vec+j) - *(mean + j) ) *  ( *(vec+i)  - *(mean + i) );
      }
  }
  return (PARAM_DATA_TYPE) invSqrtTwoPiD * sqrt(invDet) * exp(-.5 * sum);
}


//////////////////// prob_x_GM ////////////////////

/**
 * Compute the pdf of Gaussian Mixture at x by natural parameters 
 */
PARAM_DATA_TYPE MixNormal::prob_x_GM(const PointerSetToDataPoints& ps,
			      const unsigned sampleNum,
			      const PARAM_DATA_TYPE *alphas,
			      const PARAM_DATA_TYPE *means,  
			      const PARAM_DATA_TYPE *invCovs, 
			      const PARAM_DATA_TYPE *invDets,
			      const double invSqrtTwoPiD) const {
  PARAM_DATA_TYPE sum = 0.0;
  unsigned d = _numVariables;
  for ( unsigned l = 0; l < _numMixtures; l++ )
    sum += *(alphas+l) * prob_x_Gaussian(ps,sampleNum, means+l*d, invCovs+l*d*d, *(invDets+l), invSqrtTwoPiD);
  return sum;
} 

//////////////////// prob_x_GM ////////////////////

/**
 * Compute the pdf of Gaussian Mixture at ps by natural parameters 
 */
PARAM_DATA_TYPE MixNormal::prob_x_GM(const PointerSetToDataPoints& ps,
				     const unsigned sampleNum,
				     const unsigned start, // start of subVector
				     const unsigned end,   // end of subVector 
				     const PARAM_DATA_TYPE *alphas,
				     const PARAM_DATA_TYPE *means,  
				     const PARAM_DATA_TYPE *invCovs, 
				     const PARAM_DATA_TYPE *invDets,
				     const double invSqrtTwoPiD) const {
  PARAM_DATA_TYPE sum = 0;
  unsigned vecDim = end - start + 1;
  assert (vecDim > 0);
  for ( unsigned l = 0; l < _numMixtures; l++ )
    sum += *(alphas+l) * prob_x_Gaussian(ps,sampleNum, start, end, means+l*vecDim, invCovs+l*vecDim*vecDim, *(invDets+l), invSqrtTwoPiD);
  return sum;
} 


//////////////////// prob_x_GM ////////////////////

/**
 * Compute the pdf of Gaussian Mixture at vec by natural parameters
 */
PARAM_DATA_TYPE MixNormal::prob_x_GM(const PARAM_DATA_TYPE* vec,
				     const unsigned vecSize,
				     const PARAM_DATA_TYPE *alphas,
				     const PARAM_DATA_TYPE *means,  
				     const PARAM_DATA_TYPE *invCovs, 
				     const PARAM_DATA_TYPE *invDets,
				     const double invSqrtTwoPiD) const {
  PARAM_DATA_TYPE sum = 0;
  for ( unsigned l = 0; l < _numMixtures; l++ )
    sum += *(alphas+l) * prob_x_Gaussian(vec,vecSize, means+l*vecSize, invCovs+l*vecSize*vecSize, *(invDets+l), invSqrtTwoPiD);
  return sum;
} 

//////////////////// partitionCov ////////////////////

/**
 * Partitions a [dxd] matrice into
 *            
 *             [CovX]_dXxdX |  [CovXY]_dXxdY
 * [Cov]_dxd =  -----------------------------
 *	             A      |  [CovY]_dYxdY
 */
void MixNormal::partitionCov(unsigned dX,
			     PARAM_DATA_TYPE* cov,
			     PARAM_DATA_TYPE* covX,
			     PARAM_DATA_TYPE* covY,
			     unsigned d) const{
  unsigned dY = d - dX;

  for( unsigned i = 0; i < dX; i++)
    for( unsigned j = 0; j < dX; j++)
      *(covX + i*dX + j) = *(cov+i*d+j);
  for( unsigned i = dX; i < d; i++)
    for( unsigned j = dX; j < d; j++)
      *(covY+(i-dX)*dY+(j-dX)) = *(cov+i*d+j);

#if DEBUG
  DBGFPRINTF((stderr,"Full covariance matrix covXY: \n"));
  for( unsigned i = 0; i < d; i++) {
    for( unsigned j = 0; j < d; j++) {
      DBGFPRINTF((stderr,"%f ",*(cov+i*d+j)));
    }
    DBGFPRINTF((stderr,"\n"));
  }
  DBGFPRINTF((stderr,"\n"));
  DBGFPRINTF((stderr,"Partitioned covariance matrix: \n"));
  DBGFPRINTF((stderr,"covX:\n"));
  for( unsigned i = 0; i < dX; i++) {
    for( unsigned j = 0; j < dX; j++) {
      DBGFPRINTF((stderr,"%f ",*(covX+i*dX+j)));
    }
    DBGFPRINTF((stderr,"\n"));
  }
  DBGFPRINTF((stderr,"\n"));
  DBGFPRINTF((stderr,"covY:\n"));
  for( unsigned i = 0; i < dY; i++) {
    for( unsigned j = 0; j < dY; j++) {
      DBGFPRINTF((stderr,"%f ",*(covY+i*dY+j)));
    }
    DBGFPRINTF((stderr,"\n"));
  }
  DBGFPRINTF((stderr,"\n"));
#endif
}

//////////////////// moment2LDU ////////////////////

/**

 * Performs LDU decomposition of a covariance matrix, cov, into b'*D*b
 * where b is a cholesky matrix with 1s in the diagonal and D is a
 * diagonal matrix.  b and D are inverted before being returned. b is
 * also transposed.  Therefore the final b and D are such that b'*D*b = cov-1
 *
 * @param cov -- a covariance matrix (non inverted)
 * @param b -- output parameter that holds the inverse of the B matrix
 * @param invVars -- output parameter taht holds the inverse of the variances
 */
bool MixNormal::moment2LDU(PARAM_DATA_TYPE*  cov, PARAM_DATA_TYPE* b,PARAM_DATA_TYPE* invVars, unsigned d) const {
  bool malFormed = false;
  malFormed = diagCholeskyDecomp(cov,b,invVars,d);   //invVars = Cov.diagCholeskyDecomp(b).inverse();
  if(malFormed) return true;
  invertVec(invVars,invVars,d);
  //  b = (b^-1)'
  inverse(b,_tmpM,d);     // _tmpM = b^-1;
  transpose(b,_tmpM,d);   // b = _tmpM';
  return false;
}

//////////////////// LDU2moment ////////////////////

/**
 * Converts from b and S to outM (covariance matrix). outM^-1 = b'*S*b 
 * @return outM = b^-1 * invVars^-1 * (b^-1)'
 */
void MixNormal::LDU2moment(PARAM_DATA_TYPE* outM, PARAM_DATA_TYPE* b, PARAM_DATA_TYPE* invVars, unsigned d) const {

  PARAM_DATA_TYPE tmp_buffer[MAX_POINTER_SET_DIM*MAX_POINTER_SET_DIM];

  for( unsigned i = 0; i < d*d; i++ ) 
    *(outM + i) = 0.0;

  inverse(b,tmp_buffer,d);
  for ( unsigned i = 0; i < d; i++) 
    for ( unsigned j = 0; j < d; j++) 
      for ( unsigned k = 0; k< d; k++) 
	*(outM + i*d + j) += *(tmp_buffer + i*d + k) / *(invVars+k) *  *(tmp_buffer + j*d + k);
}

