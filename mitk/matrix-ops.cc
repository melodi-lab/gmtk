#include "matrix-ops.h"

/// transpose ///
/**
 * transposes a matrix. out CANNOT be the same as in.
 * 
 */
void transpose(PARAM_DATA_TYPE *  out, PARAM_DATA_TYPE * in, unsigned n) {
  PARAM_DATA_TYPE* pOut= out, *pIn=in;
  for(unsigned i=0; i<n; ++i) 
      for(unsigned j=0; j<n; ++j)
	*pOut++ = *(pIn + i + j*n); 
}

//////////////////// det ////////////////////

/**
 * find the determinent of the matrix
 * in this code, the Gaussian elimination with full pivoting
 * note:
 *     if you are dealing with positive definite matrices,
 * you might want to use Cholesky decomposition to find
 * out the determinant as well.
 *
 * @return the determinent of the matrix
 * @throws OutOfMemoryError when out of memory
 * Taken from Gang Gi's code (gmatrix.h)
*/
double determinant(double* M, unsigned n) {

  unsigned i;
  unsigned rowSizeMinusOne = n - 1;

  double det = 1.0;
  double pivot, tmp;

  double * tmpBuf = new double[n];		// this is for quick swapping
  if ( tmpBuf == NULL ) {
    cout<<"Not enough memory to allocate buffer";
    exit(-1);
  }
  double * tmpM = new double[n*n];	// we don't want to change the original one
  if ( tmpM == NULL ) {
    cout<<"Not enough memory to allocate buffer";
    exit(-1);
  }
  
  
  for(unsigned j = 0; j<n*n ; ++j) *(tmpM + j) = *(M + j); 

  double * p_b = tmpM;
  double * p_end;
  double * const p_g_end = tmpM + n*n;
  double * p;
  double * p_max;

  double * p_s;
  double * p_s_b;

  for ( i = 0; i < rowSizeMinusOne; i++ ) {
    // first we find the pivot
    pivot = fabs(*(p_b + i));
    p_max = p_b;
    p = p_b + n;
    while ( p != p_g_end ) {
      tmp = fabs(*(p + i));
      if ( tmp > pivot ) {
	p_max = p;
	pivot = tmp;
      }
      p += n;
    }

    if ( pivot == 0 )	// singular matrix
      return 0;

    // now p_max is the pointer to the pivoting row
    // we exchange the rows if necesary
    // note we need to change the sign of determinent
    // when we make the swap
    if ( p_max != p_b ) {
      memcpy(tmpBuf, p_b, n * sizeof(double));
      memcpy(p_b, p_max, n * sizeof(double));
      memcpy(p_max, tmpBuf, n * sizeof(double));
      det = -det;
    }

    // now we do the substraction to make the front be zeros
    p_end = p_b + n;
    p_s_b = p_end;
    while ( p_s_b != p_g_end ) {
      p = p_b + i + 1;	// we don't need to do the front
      p_s = p_s_b + i;
      tmp = *p_s / (*(p_b+i));
      *p_s++ = 0;
      while ( p != p_end ) {
	*p_s++ -= *p++ * tmp;
      }
      p_s_b += n;
    }

    p_b += n;
  }

  delete [] tmpBuf;

  // now we multiply the diagonal terms and get the determinant
  p = tmpM;
  for ( i = 0; i < n; i++ ) {
    det *= *(p + i);
    p += n;
  }

  delete [] tmpM;

  return det;
} // end det

//////////////////// prod ////////////////////

/**
 * calculate the product of all the elements of a vector
 *
 * @return the product of all elements
 */
PARAM_DATA_TYPE prod(PARAM_DATA_TYPE* v, unsigned n) {
  PARAM_DATA_TYPE *p = v;
  PARAM_DATA_TYPE prod = 1;
  for ( unsigned i = 0; i < n; i++)
    prod *= *p++;

  return prod;
} // end prod


/**
 * applies a transformation to a vector: multiplies a square matrix M
 * (nxn) by a vector V (n). The result is MxV and is stored in V itself.
 *  */
void matVecProduct(PARAM_DATA_TYPE * M, PARAM_DATA_TYPE* V, unsigned n) {
  PARAM_DATA_TYPE tmpV[MAX_POINTER_SET_DIM];

  for(unsigned i=0; i<n;++i)       { tmpV[i] = *(V+i);  *(V+i) = 0.0; }

  for(unsigned i=0; i<n;++i) 
    for(unsigned j=0; j<n;++j)
        *(V+i) +=  *(M+i*n+j) * tmpV[j];
} 


/**
 * multiplies a vector V (n) by a square matrix M
 * (nxn). The result is MxV and is stored in V itself.
 *  */
void vecMatProduct(PARAM_DATA_TYPE * V, PARAM_DATA_TYPE* M, unsigned n) {
  PARAM_DATA_TYPE tmpV[MAX_POINTER_SET_DIM];

  for(unsigned i=0; i<n;++i)       { tmpV[i] = *(V+i);  *(V+i) = 0.0; }

  for(unsigned i=0; i<n;++i) 
    for(unsigned j=0; j<n;++j)
        *(V+i) +=  *(M+j*n+i) * tmpV[j];
} 


//////////////////// dot ////////////////////

/**
 * dot two vectors.
 *
 * @param v1 and v2 the vectors to dot
 * @return the dot product of v1 and v2
 */
PARAM_DATA_TYPE dot(const PARAM_DATA_TYPE* v1,const PARAM_DATA_TYPE* v2,unsigned n) {
  PARAM_DATA_TYPE sum = 0.0;
  const PARAM_DATA_TYPE *p1 = v1;
  const PARAM_DATA_TYPE *p2 = v2;

  for(unsigned i=0; i<n;++i) sum += (*p1++) * (*p2++);

  return sum;

} // end dot two

/**
 * dot three vectors
 *
 * @param v1 the vector to dot
 * @param v2 the vector to dot
 * @param v3 the vector to dot
 * @return the sum of each coresponding production
 * @throws NoSuchMethodException when the three vectors have different sizes
 */
PARAM_DATA_TYPE dot(const PARAM_DATA_TYPE * v1, const PARAM_DATA_TYPE * v2,const PARAM_DATA_TYPE * v3, unsigned n) {

  PARAM_DATA_TYPE sum = 0.0;
  const PARAM_DATA_TYPE *p1 = v1;
  const PARAM_DATA_TYPE *p2 = v2;
  const PARAM_DATA_TYPE *p3 = v3;

  for(unsigned i=0; i<n;++i) sum += (*p1++) * (*p2++) * (*p3++);

  return sum;
} // end dot three


////////////////// maxVal ////////////////////

/**
 * the maximum of the elements in the vector v, of size n
 *
 * @return the maximum of the elements
 */
PARAM_DATA_TYPE maxVal(PARAM_DATA_TYPE* v, unsigned n)  {

  PARAM_DATA_TYPE max = *v;
  PARAM_DATA_TYPE *p = v;

  for(unsigned i=0; i<n;++i,p++) if( *p > max) max = *p;

  return max;
} // end maxVal


/**
 * unify the matrix
 *
 */
void unify(PARAM_DATA_TYPE* M, unsigned n) {

  PARAM_DATA_TYPE * p = M;

  unsigned shift=0;
  for (unsigned i=0; i<n; ++i) {
    for (unsigned j=0; j<n; ++j)  *(p+j)=0.0;
    *(p + shift) = 1.0;
    shift++;
    p += n;
  }

} // end unify



//////////////////// inverse ////////////////////

/**
 * find the inverse of the matrix using LU decomposition
 * by Crout's method with partial pivoting
 *
 * @param M the matrix to invert
 * @param INV output parameter result of the inversion.  INV can be the same as M
 * @return whether the matrix is singular or not (true -> matrix is singular)
 */

bool inverse(PARAM_DATA_TYPE * M, PARAM_DATA_TYPE * INV, unsigned n)  {

  PARAM_DATA_TYPE det = 1.0;

  int iMax;
  unsigned *indx = new unsigned [n];
  if ( indx == NULL )
    error("out of memory");
  PARAM_DATA_TYPE big, dum, sum, temp;
  PARAM_DATA_TYPE *vv = new PARAM_DATA_TYPE[n];
  if ( vv == NULL )
    error("out of memory");

  PARAM_DATA_TYPE tiny = (PARAM_DATA_TYPE)1.0e-100;
  PARAM_DATA_TYPE **lu = new PARAM_DATA_TYPE* [n];
  if ( lu == NULL )
    error("out of memory");
  for ( unsigned i = 0; i < n; i++ ) {
    lu[i] = new PARAM_DATA_TYPE [n];
    if ( lu[i] == NULL )
      error("out of memory");
  }

  /* first, we copy the data */
  for ( unsigned i = 0; i < n ; i++ )
    memcpy(lu[i], M + i * n, n * sizeof(PARAM_DATA_TYPE));

  for ( unsigned i = 0; i < n; i++ ) {
    big = 0;
    for ( unsigned j = 0; j < n; j++ )
      if ( (temp = fabs(lu[i][j])) > big )
	big = temp;
    if ( big == 0 ) {
      det = 0;		// singular matrix
      //return det;
      return true;  //matrix is singular
    }

    vv[i] = 1.0 / big;
  }

  for ( unsigned j = 0; j < n; j++ ) {
    for ( unsigned i = 0; i < j; i++ ) {
      sum = lu[i][j];
      for ( unsigned k = 0; k < i; k++ )
	sum -= lu[i][k] * lu[k][j];
      lu[i][j] = sum;
    }
    big = 0.0;
    iMax = (int) j;
    for ( unsigned i = j; i < n; i++ ) {
      sum = lu[i][j];
      for ( unsigned k = 0; k < j; k++ )
	sum -= lu[i][k] * lu[k][j];

      lu[i][j] = sum;
      if ( (dum = vv[i] * fabs(sum)) >= big ) {
	big = dum;
	iMax = i;
      }
    }
    if ( (int)j != iMax ) {
      for ( unsigned k = 0; k < n; k++ ) {
	dum = lu[iMax][k];
	lu[iMax][k] = lu[j][k];
	lu[j][k] = dum;
      }
      det = -det;
      vv[iMax] = vv[j];
    }
    indx[j] = iMax;
    if ( lu[j][j] == 0.0 )
      lu[j][j] = tiny;
    if ( j != n - 1 ) {
      dum = 1.0 / (lu[j][j]);
      for (unsigned i = j + 1; i < n; i++ )
	lu[i][j] *= dum;
    }
  }

  for ( unsigned i = 0; i < n; i++ )
    det *= lu[i][i];


  // now use the LU matrices to solve linear equations

  for ( unsigned j = 0; j < n; j++ ) {
    //memset(vv, 0, _rowSize * sizeof(T));
     for(unsigned int ii=0;ii<n;++ii) (vv[ii])=0;
    vv[j] = (PARAM_DATA_TYPE)1;
    iMax = -1;

    for ( unsigned i = 0; i < n; i++ ) {
      unsigned k = indx[i];
      sum = vv[k];
      vv[k] = vv[i];
      if ( iMax >= 0 )
	for ( unsigned k = iMax; k <= i - 1; k++)
	  sum -= lu[i][k] * vv[k];
      else if ( sum )
	iMax = i;
      vv[i] = sum;
    }
    for ( int i = n - 1; i >= 0; i-- ) {
      sum = vv[i];
      for ( unsigned k = (unsigned) i + 1; k < n; k++)
	sum -= lu[i][k] * vv[k];
      vv[i] = sum / lu[i][i];
    }

    for ( unsigned i = 0; i < n; i++ ) {
      *(INV + i * n + j) = vv[i];
    }
  }

  for ( unsigned i = 0; i < n; i++ )
    delete [] lu[i];
  delete [] lu;
  delete [] indx;
  delete [] vv;

  return false; //not singular
} // end inverse


//////////////////// invertVec ////////////////////

/**
 * inverse the vector as v^(-1)
 *
 * @param v vector to invert
 * @param output parameter inv vector with each element is reverse of v. inv can be the same as v
 * @param n the size of the vector
 * @exception OverFlowException when divided by zero
 */
void invertVec(PARAM_DATA_TYPE* v, PARAM_DATA_TYPE * inv, unsigned n)  {

  PARAM_DATA_TYPE *p = v;
  PARAM_DATA_TYPE *inv_p = inv;

  for(unsigned i=0; i<n;++i) {
    if ( *p == 0 )
      error("invertVec(): cannot inverse vector because of a zero in the vector to invert.");
    *inv_p++ = 1.0 / (*p++);
  }

} // end invertVec



//////////////////// pointProdAdd ////////////////////

/**
 * dot product a vector with itself without summing, multiply by a constant and add to a 2nd vector
 *
 * @param out the vector to add the pairwise product to
 * @param ps the set of pointers to the vector we want to dot product with itself without summation
 * @param sampleNum the sample number 
 * @param multiplier the constant to multiply the dot product by
 * @parm n the dimension of the pointer set is equal to ps.dim
 * @throws NoSuchMethodException when the two vectors have different sizes
 */
void pointProdAdd(PARAM_DATA_TYPE * out, const PointerSetToDataPoints& ps, unsigned sampleNum, PARAM_DATA_TYPE multiplier,unsigned n)  {

  PARAM_DATA_TYPE *pOut = out;

  for ( unsigned i = 0; i < n; i++ ) 
    *pOut++ += *(ps.start[i]+sampleNum*ps.skip) * *(ps.start[i]+sampleNum*ps.skip) * multiplier;

} // end pointProd


//////////////////// pointProdSub ////////////////////

/**
 * dot product of two vectors without summing and substract from the third vector
 *
 * @param v1 the vector to pairwise product
 * @param v2 the vector to pairwise product
 * @param v3 the vector to substract the pairwise product from
 * @throws NoSuchMethodException when the two vectors have different sizes
 */
void pointProdSub(PARAM_DATA_TYPE *v1, PARAM_DATA_TYPE * v2, PARAM_DATA_TYPE * v3, unsigned n)  {

  PARAM_DATA_TYPE *p1 = v1;
  PARAM_DATA_TYPE *p2 = v2;
  PARAM_DATA_TYPE *p3 = v3;

  for ( unsigned i = 0; i < n; i++ ) *p3++ -= (*p1++) * (*p2++);

} // end pointProd


/**
 * dot producting two vectors without summing, mul my constant and add to ouput matrix
 *
 * @param out the output matrix
 * @param ps the pointer set to the vector wich is mutiplied with itself (nx1 * 1xn -> nxn)
 * @param sampleNum sample number in the observation matrci to access the vector at
 * @param multiplier
 * @param n dimension of pointer set == ps->dim
 */
void vecProdAdd( PARAM_DATA_TYPE * out, const PointerSetToDataPoints& ps, unsigned sampleNum, PARAM_DATA_TYPE multiplier, unsigned n)  {

  PARAM_DATA_TYPE *pOut = out;

   for(unsigned j = 0; j < n; ++j)
    for(unsigned k = 0; k < n; ++k) {
      *(pOut + j*n + k)  +=   *(ps.start[j]+sampleNum*ps.skip) * *(ps.start[k]+sampleNum*ps.skip) * multiplier;
    }

} // vecProdSub


/**
 * dot producting two vectors without summing and substract from the third vector
 *
 * @param v the vector to pairwise product
 */
void vecProdSub( PARAM_DATA_TYPE * out, PARAM_DATA_TYPE *v1, PARAM_DATA_TYPE * v2, unsigned n)  {

  PARAM_DATA_TYPE *p1 = v1;
  PARAM_DATA_TYPE *p2 = v2;
  PARAM_DATA_TYPE *p3 = out;

   for(unsigned j = 0; j < n; ++j)
    for(unsigned k = 0; k < n; ++k) {
      *(p3 + j*n + k)  -= ( *(p1 + j) * *(p2 + k) );
    }

} // vecPrdSub



//////////////////// choleskyDecomp ////////////////////

/**
 * do Cholesky decomposition on a symetrix matrix
 *    M = L' L where R is a lower triangular matrix
 *    M = U U' where U is an upper triangular matrix
 * the details of the program goes like this
 *    U_{i,i} = (a_{i,i} - \sum_{k=0}^{i-1}U_{k,i}^2)^{1/2}
 *					for i = 0, ..., n-1
 *    U_{i,j} = (a_{i,j} - \sum-{k=0}^{i-1}U_{k,i}U_{k,j})
 *              / U_{i,i}
 *					for j = i+1, ..., n-1
 *
 * @param M the matrix to decompose
 * @param chol the output upper triangular Cholesky matrix
 * @return true when the matrix is NOT positive definite, false otherwise
 */
bool choleskyDecomp(PARAM_DATA_TYPE * M, PARAM_DATA_TYPE * chol, unsigned n)  {

  unsigned i, j, k;
  double tmp, tmp_diag;
  PARAM_DATA_TYPE* p_b = M;
  PARAM_DATA_TYPE* p;

  PARAM_DATA_TYPE* ch_p_b = chol;
  PARAM_DATA_TYPE* ch_p;
  PARAM_DATA_TYPE* ch_i_p;
  PARAM_DATA_TYPE* ch_j_p;

  // Initializes the output cholesky matrix to zeros
  for(unsigned i=0; i < n; ++i) 
    for(unsigned j=0; j < n; ++j) 
      *(chol + i*n +j) = 0.0;

  for ( i = 0; i < n; i++ ) {
    p = p_b + i;
    ch_p = ch_p_b + i;

    // first treat the U_{i,i} term
    tmp = (double) *p++;
    ch_i_p = chol + i;

    for ( k = 0; k < i; k++ ) {
      tmp -= (*ch_i_p) * (*ch_i_p);
      ch_i_p += n;
    }
    if ( tmp <= 0 ) {
       return true;
    }
    tmp_diag = (PARAM_DATA_TYPE)sqrt(tmp);
    *ch_p++ = tmp_diag;

    // now treat the U_{i,j} term
    for ( j = i + 1; j < n; j++ ) {
      tmp = (double) *p++;
      ch_i_p = chol + i;
      ch_j_p = chol + j;

      for ( k = 0; k < i; k++ ) {
	tmp -= (*ch_i_p) * (*ch_j_p);
	ch_i_p += n;
	ch_j_p += n;
      }
      tmp /= tmp_diag;

      *ch_p++ = (PARAM_DATA_TYPE)tmp;
    }

    p_b += n;
    ch_p_b += n;
  }
  return false;  
} // end choleskyDecomp



//////////////////// diagCholeskyDecomp ////////////////////

/**
 * do diagonal Cholesky decomposition of a symmetric matrix
 *     S = C'DC
 *   where C is Cholesky matrix with diagonal 1 and
 *   D is a diagonal matrix
 *
 * @param S is the input matrix
 * @param chol is the output Cholesky matrix with diagonal 1
 * @param diagonal is the output diagonal matrix 
 * @return true if the input matrix is NOT positive definite, false otherwise
 */
bool diagCholeskyDecomp(PARAM_DATA_TYPE * S, PARAM_DATA_TYPE * chol, PARAM_DATA_TYPE * diagonal,unsigned n)  {
  //cout<<"Cov BEFORE:\n";
  //for(unsigned i=0;i<n*n;++i) cout<<*(S+i)<<" ";
  //cout<<endl; 
  bool malFormed = choleskyDecomp(S,chol,n);
  if(malFormed == true) {
    //cout<<"Cov AFTER:\n";
    //for(unsigned i=0;i<n*n;++i) cout<<*(S+i)<<" ";
    //cout<<endl;
    //    error("ERRROR");
    return true;
  }
  double tmp;
  PARAM_DATA_TYPE* p_b = chol;
  PARAM_DATA_TYPE* p_end;
  PARAM_DATA_TYPE* p;
  PARAM_DATA_TYPE* v_p = diagonal;

  for ( unsigned i = 0; i < n; i++ ) {
    p = p_b + i;
    p_end = p_b + n;
    tmp = 1.0 / *p;
    *v_p++ = *p * *p;
    *p++ = (PARAM_DATA_TYPE)1.0;

    while ( p != p_end ) {
      *p++ *= tmp;
    }

    p_b += n;
  }

  return false; // No problem encountered. The input matrix is
  // positive definite.
} // end diagCholeskyDecomp
