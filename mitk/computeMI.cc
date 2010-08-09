#include <string.h>

#include "mixNormalCollection.h"
#include "matrix-ops.h"


#ifndef M_2PI
#define M_2PI 6.28318530717958647692
#endif


//////////////////// MixNormalCollection::computeMI ////////////////////

/**
 *  computes MI between the first set and second set of each Gaussian
 *  Mixture using the law of large numbers (LLN) by sampling data 
 */
void MixNormalCollection::computeMI(FILE *out_fp, unsigned nSamples,
				    const RangeSetCollection &rangeCol,
				    FILE* rangeFileFP,
				    bool marginalizeFirstParentOut) {
  RangeSet range;
  PARAM_DATA_TYPE I, Hx, Hy, Hxy;
  unsigned firstSetSize, scndSetSize, k = 0;
  char line[MAX_LINE_LEN];
  char *pos = line;
  for ( MixNormal *p = _ftrMI; p != _ftrMI_endp; p++ ) {
    range = rangeCol.rs[k];
    firstSetSize = range.getSize(0);
    scndSetSize = range.getSize(1);
    if( scndSetSize == 0 )  // only one set?
      I = p->computeEntropy(nSamples,Hx);  // compute entropy
    else {
      if(marginalizeFirstParentOut) {
	//DBGFPRINTF((stderr,"Marginalizing the first parent.\n"));
	I = p->computeMIWithFirstParentOut(nSamples, firstSetSize, Hx, Hy, Hxy); // compute MI
      }
      else {
	I = p->computeMI(nSamples, firstSetSize, Hx, Hy, Hxy); // compute MI
      }
    }
    // remove comment lines
    do {
      fgets(line,MAX_LINE_LEN,rangeFileFP); 
      if( ( pos = strchr(line,'#') ) != NULL ) *pos = '\0';  //remove comments
      else pos = line;
      while(isspace(*pos)) ++pos;  //eat spaces
    } while(*pos=='\0' || *pos =='\n');
    // --------------------
    if(line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
    if( scndSetSize != 0 ) {
      fprintf(out_fp,"%s %e %d %d %e %e %e\n", 
	      line, I,p->getNumSamples(),(int)p->getNumEMIterToFinish(),Hx,Hy,Hxy);
    }
    else {
      fprintf(out_fp,"%s %e %d %d\n", 
	      line, I,p->getNumSamples(),(int)p->getNumEMIterToFinish());
    }
    k++;
  }
  fflush(out_fp);
}


//////////////////// MixNormal::computeMI ////////////////////

/** 
 * Compute the mutual information between X and Y by the law of large
 * numbers. First generate nSamples samples from the joint
 * distribution of X and Y and then calculate the mutual information as:
 *  I(X;Y) = H(X) + H(Y) - H(X,Y) where H() denotes the entropy
 *  H(X) = sum log( 1/p(x) ) / nSamples
 * Equivalently, 
 *  I = [ sum log( p(x,y)/p(x)p(y) ) ] / nSamples
 */
PARAM_DATA_TYPE MixNormal::computeMI(unsigned  nSamples, unsigned  dX, PARAM_DATA_TYPE& Hx, PARAM_DATA_TYPE& Hy, PARAM_DATA_TYPE& Hxy) {
  unsigned dY = _numVariables - dX;
  PARAM_DATA_TYPE probX, probY, probXY; 
  PARAM_DATA_TYPE sampleXY[MAX_POINTER_SET_DIM];
  //, sampleX[MAX_POINTER_SET_DIM], sampleY[MAX_POINTER_SET_DIM];

  if( ( dX < 1 ) || ( dX >= _numVariables ) ){
    error("ERROR:Inappropriate dX (%d).  Should be between 1 and %d.\n",dX,_numVariables-1);
  }


#ifdef DEBUG
    DBGFPRINTF((stderr,"------------- B matrix is  -----------\n"));
    for (unsigned ii=0; ii<_numVariables; ++ii) {
      for (unsigned jj=0; jj<_numVariables; ++jj) {
	DBGFPRINTF((stderr,"%f ",*(_b+ii*_numVariables+jj)));
      }
      DBGFPRINTF((stderr,"\n"));
    }
    DBGFPRINTF((stderr,"\n"));
#endif

  DBGFPRINTF((stderr,"In MixNormal::computeMI, marginalizing.\n"));
  marginalize(dX);


  PARAM_DATA_TYPE* sampleX=sampleXY;
  PARAM_DATA_TYPE* sampleY=sampleXY+dX;

  Hx = 0; 
  Hy = 0; 
  Hxy = 0;
  DBGFPRINTF((stderr,"Starting sampling to compute MI\n"));
  for( unsigned i = 0; i < nSamples; i++){

    sample(sampleXY);   //    sampleXY = sample();   


#ifdef DEBUG
#define PRINT_SAMPLES
#define MAX_SAMPLES 10
#ifdef PRINT_SAMPLES
    if(i< MAX_SAMPLES) {
      for (unsigned i=0; i<_numVariables; ++i) {
	DBGFPRINTF((stderr,"%f ",*(sampleXY+i)));
      }
      DBGFPRINTF((stderr,"\n"));
    }
#endif
#endif
    //probX  = prob_x_GM(sampleX, dX, _alphas, _meansX, _invCovsX, _invDetsX,_inv_pow_sqrt_2piX);
    probX  = prob_x(sampleX, dX, _meansX, _bX, _invVarsX, _invDetsX,_inv_pow_sqrt_2piX);
    //probY  = prob_x_GM(sampleY, dY,_alphas, _meansY, _invCovsY, _invDetsY,_inv_pow_sqrt_2piY);
    probY  = prob_x(sampleY, dY, _meansY, _bY, _invVarsY, _invDetsY,_inv_pow_sqrt_2piY);
    probXY = prob_x(sampleXY,_numVariables);
    assert( (probX >= 0) && (probY >= 0) && (probXY >= 0) );
    if( (probX > 0) && (probY > 0) && (probXY > 0) ){
      Hx  = Hx  - log(probX);
      Hy  = Hy  - log(probY);
      Hxy = Hxy - log(probXY);
    }
  }

  PARAM_DATA_TYPE invNlog2 = 1.0 / (PARAM_DATA_TYPE)(nSamples*log(2.0));
  Hx  = Hx  * invNlog2;
  Hy  = Hy  * invNlog2;
  Hxy = Hxy * invNlog2;

  // I = H(X) + H(Y) - H(X;Y) 
  return Hx+Hy-Hxy;
}



//////////////////// MixNormal::computeMIWithFirstParentOut ////////////////////

/** 
 * Compute the mutual information between X and Y by the law of large
 * numbers. First generate nSamples samples from the joint
 * distribution of X and Y and then calculate the mutual information as:
 *  I(X;Y) = H(X) + H(Y) - H(X,Y) where H() denotes the entropy
 *  H(X) = sum log( 1/p(x) ) / nSamples
 * Equivalently, 
 *  I = [ sum log( p(x,y)/p(x)p(y) ) ] / nSamples
 */
PARAM_DATA_TYPE MixNormal::computeMIWithFirstParentOut(unsigned  nSamples, unsigned  dX, PARAM_DATA_TYPE& Hx, PARAM_DATA_TYPE& Hy, PARAM_DATA_TYPE& Hxy) {
  if(_numVariables > 2) {
    DBGFPRINTF((stderr,"In MixNormal::computeMIWithFirstParent, _numVariables=%d\n",_numVariables));
    // Remove first parent
    unsigned origDim=_numVariables;
    _numVariables-=1;
    unsigned dim=_numVariables;
    for(unsigned l=0; l<_numMixtures; ++l) {
      if ( _fullCoVar ) {
	LDU2moment(_invCovsXY+l*origDim*origDim, (_b + l*origDim*origDim), (_invVars + l*origDim), origDim);     //_invCovsXY = LDU2moment(_b[l], _invVars[l]);
	
	// remove first parent here
	unsigned newIndex=0;
	for(unsigned j=0; j<origDim; ++j) {
	  if(j==dX) { // row == first parent position -> delete whole row
	    continue;
	  }
	  for(unsigned k=0; k<origDim; ++k) {
	  if(k==dX) { // column = parent pos -> delete element
	    continue;
	  }
	  *(_invCovsXY+l*dim*dim+newIndex) = *(_invCovsXY+l*origDim*origDim+j*origDim+k);
	  newIndex++;
	  }
	}
	
	// Reconvert to B matrix
	moment2LDU(_invCovsXY+l*dim*dim,_b+l*dim*dim,_invVars+l*dim,dim);
	
	*(_invDets+l) = prod(_invVars+l*dim,dim); //*(_invDets+l) = _invVars[l].prod();
	
	_inv_pow_sqrt_2pi = 1.0 / pow(sqrt(M_2PI), (int)_numVariables);
	
      }
      unsigned newIndex=0;
      for(unsigned j=0; j<origDim; ++j) {
	if(j==dX) continue;
	*(_means+l*dim +newIndex)=*(_means+l*origDim +j);
	newIndex++;
      }
      // the below is already done when calling moment2LDU
      //    for(unsigned j=dX+1; j<origDim; ++j) {
      //  *(_invVars+l*dim +j-1)=*(_invVars+l*origDim +j);
      //}
    }
  }

  unsigned dY = _numVariables - dX;
  PARAM_DATA_TYPE probX, probY, probXY; 
  PARAM_DATA_TYPE sampleXY[MAX_POINTER_SET_DIM];

  if( ( dX < 1 ) || ( dX >= _numVariables )){
    error("ERROR:Inappropriate dX (%d).  Should be between 1 and %d.\n",dX,_numVariables-1);
  }


  //  marginalizeWithFirstParentOut(dX);
  marginalize(dX);

  PARAM_DATA_TYPE* sampleX=sampleXY;
  PARAM_DATA_TYPE* sampleY=sampleXY+dX;

  Hx = 0; 
  Hy = 0; 
  Hxy = 0;
  for( unsigned i = 0; i < nSamples; i++){
    sample(sampleXY);   //    sampleXY = sample();   

    probX  = prob_x(sampleX, dX, _meansX, _bX, _invVarsX, _invDetsX,_inv_pow_sqrt_2piX);
    probY  = prob_x(sampleY, dY, _meansY, _bY, _invVarsY, _invDetsY,_inv_pow_sqrt_2piY);
    probXY = prob_x(sampleXY,_numVariables);
    assert( (probX >= 0) && (probY >= 0) && (probXY >= 0) );
    if( (probX > 0) && (probY > 0) && (probXY > 0) ){
      Hx  = Hx  - log(probX);
      Hy  = Hy  - log(probY);
      Hxy = Hxy - log(probXY);
    }
  }

  PARAM_DATA_TYPE invNlog2 = 1.0 / (PARAM_DATA_TYPE)(nSamples*log(2.0));
  Hx  = Hx  * invNlog2;
  Hy  = Hy  * invNlog2;
  Hxy = Hxy * invNlog2;

  // I = H(X) + H(Y) - H(X;Y) 
  return Hx+Hy-Hxy;
}
