#include <iostream>
#include <string.h>

#include "data-points.h"


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
 * @throws NoSuchMethodException when the matrix is not square
 * @throws OutOfMemoryError when out of memory
 */
double determinant(double* M, unsigned n);

void transpose(PARAM_DATA_TYPE *  out, PARAM_DATA_TYPE * in, unsigned n);
void matVecProduct(PARAM_DATA_TYPE * M, PARAM_DATA_TYPE* V, unsigned n);
void vecMatProduct(PARAM_DATA_TYPE * V, PARAM_DATA_TYPE* M, unsigned n);
void unify(PARAM_DATA_TYPE* M, unsigned n);
PARAM_DATA_TYPE prod(PARAM_DATA_TYPE* v, unsigned n);
void pointProdAdd(PARAM_DATA_TYPE * out, const PointerSetToDataPoints& ps, unsigned sampleNum, PARAM_DATA_TYPE multiplier,unsigned n);
void pointProdSub(PARAM_DATA_TYPE *v1, PARAM_DATA_TYPE * v2, PARAM_DATA_TYPE * v3, unsigned n);
void vecProdSub(PARAM_DATA_TYPE *v1, PARAM_DATA_TYPE * v2, PARAM_DATA_TYPE * v3, unsigned n);
void vecProdAdd(PARAM_DATA_TYPE * out, const PointerSetToDataPoints& ps, unsigned sampleNum, PARAM_DATA_TYPE multiplier, unsigned n);
PARAM_DATA_TYPE maxVal(PARAM_DATA_TYPE* v, unsigned n);
PARAM_DATA_TYPE dot(const PARAM_DATA_TYPE* v1,const PARAM_DATA_TYPE* v2,unsigned n);
PARAM_DATA_TYPE dot(const PARAM_DATA_TYPE * v1, const PARAM_DATA_TYPE * v2,const PARAM_DATA_TYPE * v3, unsigned n);
bool inverse(PARAM_DATA_TYPE * M, PARAM_DATA_TYPE * INV, unsigned n);
void invertVec(PARAM_DATA_TYPE* v, PARAM_DATA_TYPE * inv, unsigned n);
bool choleskyDecomp(PARAM_DATA_TYPE * M, PARAM_DATA_TYPE * chol, unsigned n);
bool diagCholeskyDecomp(PARAM_DATA_TYPE * S, PARAM_DATA_TYPE * chol, PARAM_DATA_TYPE * diagonal,unsigned n);

