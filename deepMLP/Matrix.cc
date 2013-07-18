
#if defined(HAVE_CONFIG_H)
#  include <config.h>
#endif

#include <assert.h>

#if defined(HAVE_MKL)
#  include "mkl.h"
#  include "mkl_lapacke.h"
#  include "mkl_spblas.h"
#  include "mkl_trans.h"
#elif defined(HAVE_BLAS)
extern "C" {            /* Assume C declarations for C++ */
#  include <cblas.h>
}
#endif

#if defined(USE_PHIPAC)
extern "C" {            /* Assume C declarations for C++ */
  void
    phipac_dgemm(char* transA, char* transB,
    int* M, int* N, int* K,
    double* alpha,
    double* A, int* Astride,
    double* B, int* Bstride,
    double* beta,
    double* C, int* Cstride);
}
#endif

#include "Matrix.h"

#if !defined(HAVE_MKL)
void 
  my_Domatcopy(char trans, int rows, int cols, 
  const double alpha, double const *A, int lda, 
  double * B, int ldb)
{
  assert(trans == 'n' || trans == 't');
  unsigned Ainc1 = (trans == 'n' || trans == 'N') ? 1   : lda;
  unsigned Ainc2 = (trans == 'n' || trans == 'N') ? lda : 1;

  double       *Bcol = B;
  double       *pastB = B + cols * ldb;
  double const *Apos = A;
  do {
    double       *Bp = Bcol;
    double       *colEnd = Bp + rows;
    double const *Ap = Apos;
    do {
      *(Bp++) = alpha * *Ap;
      Ap += Ainc1;
    } while (Bp != colEnd);
    Bcol += lda;
    Apos += Ainc2;
  } while (Bcol != pastB);
}

void
  my_Domatadd(char aTrans, char bTrans,
  int M, int N, 
  const double alpha, double const *A, int lda,
  const double beta,  double const *B, int ldb,
  double *C, int ldc)
{
  assert(aTrans == 'n' || aTrans == 't');
  assert(bTrans == 'n' || bTrans == 't');
  // MKL docs say A, B, C must not overlap, but it appears to work
  // even when one of the sourc matrices is also the destination.
  // The training code will call this with A=C or B=C.

  unsigned Ainc1 = (aTrans == 'n' || aTrans == 'N') ? 1   : lda;
  unsigned Ainc2 = (aTrans == 'n' || aTrans == 'N') ? lda : 1;
  unsigned Binc1 = (bTrans == 'n' || bTrans == 'N') ? 1   : ldb;
  unsigned Binc2 = (bTrans == 'n' || bTrans == 'N') ? ldb : 1;

  double       *Ccol = C;
  double       *pastC = C + N * ldc;
  double const *Bpos = B;
  double const *Apos = A;
  do {
    double       *Cp = Ccol;
    double       *colEnd = Cp + M;
    double const *Bp = Bpos;
    double const *Ap = Apos;
    do {
      *(Cp++) = alpha * *Ap + beta * *Bp;
      Bp += Binc1;
      Ap += Ainc1;
    } while (Cp != colEnd);
    Ccol += ldc;
    Bpos += Binc2;
    Apos += Ainc2;
  } while (Ccol != pastC);
}
#endif

Matrix Vector::AsMatrix(int numR, int numC) const {
  assert(numR * numC == _len && _inc == 1);
  return Matrix(_start, numR, numC, numR, false);
}

void MutableVector::CopyFrom(const Vector & vec) const {
  assert (_len == vec.Len());
  cblas_dcopy(_len, vec.Start(), vec.Inc(), Start(), Inc());
}

MutableMatrix MutableVector::AsMatrix(int numR, int numC) const {
  return Vector::AsMatrix(numR, numC);
}

MutableMatrix MutableVector::AsMatrix(int numR) const {
  return AsMatrix(numR, Len() / numR);
}

Matrix Vector::AsMatrix(int numR) const {
  return AsMatrix(numR, Len() / numR);
}

void MutableVector::Axpby(const VecScal & expr, double a, double b) const {
  assert (_len == expr.Len());
#if HAVE_CBLAS_DAXPBY
  cblas_daxpby(_len, a * expr.a, expr.v.Start(), expr.v.Inc(), b, Start(), Inc());
#else
  cblas_dscal(_len, b, Start(), Inc());
  cblas_daxpy(_len, a * expr.a, expr.v.Start(), expr.v.Inc(), Start(), Inc());
#endif
}

const MutableMatrix & MutableMatrix::operator=(const MatMatMult & expr) const {
  return Dgemm(expr.a, expr.A, expr.B, 0.0);
}

const MutableMatrix & MutableMatrix::operator+=(const MatMatMult & expr) const {
  return Dgemm(expr.a, expr.A, expr.B, 1.0);
}

const MutableMatrix & MutableMatrix::Dgemm(double a, const Matrix & A, const Matrix & B, double b) const {
  if (_trans) {
    Trans().Dgemm(a, B.Trans(), A.Trans(), b);
    return *this;
  }

  assert(_numR == A.NumR() && _numC == B.NumC() && A.NumC() == B.NumR());
#if USE_PHIPAC

  // I _think_ PHiPAC assumes row-major (C) order, but Galen's code
  // seems to be column-major.

  int numR   = _numR;
  int numC   = _numC;
  int ANumC  = A.NumC();
  int ALd    = A.Ld();
  int BLd    = B.Ld();
  int thisLd = Ld();
  char Aop   = A.IsTrans() ? 'T' : 'N';
  char Bop   = B.IsTrans() ? 'T' : 'N';
  phipac_dgemm(&Aop, &Bop, 
    &numR, &numC, &ANumC, 
    &a, 
    const_cast<double *>(A.Start()), &ALd, 
    const_cast<double *>(B.Start()), &BLd, 
    &b, 
    Start(), &thisLd);
#else
  cblas_dgemm(CblasColMajor, A.IsTrans() ? CblasTrans : CblasNoTrans, B.IsTrans() ? CblasTrans : CblasNoTrans,
    _numR, _numC, A.NumC(), a, A.Start(), A.Ld(), B.Start(), B.Ld(), b, Start(), Ld());
#endif
  return *this;
}

const MutableMatrix & MutableMatrix::operator=(const MatScaledSum & expr) const {
  assert (NumR() == expr.NumR());
  assert (NumC() == expr.NumC());

  bool aTrans = expr.A.IsTrans() ^ IsTrans(), bTrans = expr.B.IsTrans() ^ IsTrans();
#if HAVE_MKL
  MKL_Domatadd ('c', aTrans ? 't' : 'n', bTrans ? 't' : 'n',
    _numR, _numC, expr.a, expr.A.Start(), expr.A.Ld(),
    expr.b, expr.B.Start(), expr.B.Ld(), Start(), _ld);
#else
  my_Domatadd (aTrans ? 't' : 'n', bTrans ? 't' : 'n',
    _numR, _numC, expr.a, expr.A.Start(), expr.A.Ld(),
    expr.b, expr.B.Start(), expr.B.Ld(), Start(), _ld);  
#endif
  return *this;
}

const MutableMatrix & MutableMatrix::operator+=(const MatScal & expr) const {
  return operator=(MatScaledSum(*this, expr));
}

const MutableMatrix & MutableMatrix::operator=(const MatScal & expr) const {
  assert (NumR() == expr.NumR());
  assert (NumC() == expr.NumC());

  bool trans = expr.A.IsTrans() ^ _trans;
#if HAVE_MKL
  MKL_Domatcopy('c', trans ? 't' : 'n', expr.A.DeepNumR(), expr.A.DeepNumC(), expr.a, expr.A.Start(), expr.A.Ld(), Start(), Ld());
#else
  my_Domatcopy(trans ? 't' : 'n', expr.A.DeepNumR(), expr.A.DeepNumC(), expr.a, expr.A.Start(), expr.A.Ld(), Start(), Ld());
#endif
  return *this;
}

const MutableVector & MutableVector::operator*=(double a) const {
  return operator=(*this * a);
}

const MutableMatrix & MutableMatrix::operator*=(double a) const {
  return operator=(*this * a);
}

const MutableMatrix & MutableMatrix::operator-=(const MatScal & expr) const {
  return operator+= (-1 * expr);
}

const MutableMatrix & MutableMatrix::operator-=(const MatMatMult & expr) const {
  return operator+= (-1 * expr);
}

const MutableVector & MutableVector::operator-=(const VecScal & expr) const {
  return operator+= (-1.0 * expr);
}

const MutableVector & MutableVector::operator=(const VecScaledSum & expr) const {
  assert (Len() == expr.Len());
#if HAVE_MKL
  MKL_Domatadd('c', 'n', 'n', 1, Len(), expr.a, expr.x.Start(), expr.x.Inc(), expr.b, expr.y.Start(), expr.y.Inc(), Start(), Inc());
#else
  my_Domatadd('n', 'n', 1, Len(), expr.a, expr.x.Start(), expr.x.Inc(), expr.b, expr.y.Start(), expr.y.Inc(), Start(), Inc());
#endif
  return *this;
}

const MutableVector & MutableVector::operator=(const VecScal & expr) const {
  assert (Len() == expr.Len());
#if HAVE_MKL
  MKL_Domatcopy('c', 'n', 1, Len(), expr.a, expr.v.Start(), expr.v.Inc(), Start(), Inc());
#else
  my_Domatcopy('n', 1, Len(), expr.a, expr.v.Start(), expr.v.Inc(), Start(), Inc());
#endif
  return *this;
}

const MutableVector & MutableVector::operator+=(const VecScal & expr) const {
  assert (Len() == expr.Len());
  cblas_daxpy(Len(), expr.a, expr.v.Start(), expr.v.Inc(), Start(), Inc());
  return *this;
}

VecScal Vector::operator-() const {
  return *this * -1;
}

MatScal Matrix::operator-() const {
  return *this * -1;
}

MatScal operator*(double b, const MatScal & ms) {
  return MatScal(ms.A, ms.a * b);
}

MatScal operator*(const MatScal & ms, double b) {
  return b * ms;
}

MatMatMult operator*(const MatScal & A, const MatScal & B) {
  return MatMatMult(A.A, B.A, A.a * B.a);
}

MatScaledSum operator*(double c, const MatScaledSum & mss) {
  return mss.A * mss.a * c + mss.B * mss.b * c;
}

MatScaledSum operator*(const MatScaledSum & mss, double c) {
  return c * mss;
}

MatMatMult operator*(double b, const MatMatMult & mmm) {
  return MatMatMult(mmm.A, mmm.B, mmm.a * b);
}

MatMatMult operator*(const MatMatMult & mmm, double b) {
  return b * mmm;
}

MatScaledSum operator+(const MatScal & ax, const MatScal & by) {
  return MatScaledSum(ax, by);
}

MatScaledSum operator-(const MatScal & ax, const MatScal & by) {
  return MatScaledSum(ax, -by);
}

VecScal operator*(double b, const VecScal & vs) {
  return VecScal(vs.v, b * vs.a);
}

VecScal operator*(const VecScal & vs, double b) {
  return b * vs;
}

double operator*(const VecScal & ax, const VecScal & by) {
  assert(ax.v.Len() == by.v.Len());

  return ax.a * by.a * cblas_ddot(ax.v.Len(), ax.v.Start(), ax.v.Inc(), by.v.Start(), by.v.Inc());
}

VecScaledSum operator+(const VecScal & ax, const VecScal & by) {
  return VecScaledSum(ax, by);
}

VecScaledSum operator-(const VecScal & ax, const VecScal & by) {
  return VecScaledSum(ax, -by);
}

VecScaledSum operator*(double c, const VecScaledSum & vss) {
  return vss.x * vss.a * c + vss.y * vss.b * c;
}

VecScaledSum operator*(const VecScaledSum & vss, double c) {
  return c * vss;
}

void MutableMatrix::CopyFrom(const Matrix & mat, double scale) const {
  assert (NumR() == mat.NumR() && NumC() == mat.NumC());

  char trans = (IsTrans() ^ mat.IsTrans()) ? 't' : 'n';
#if HAVE_MKL
  mkl_domatcopy('c', trans, mat.DeepNumR(), mat.DeepNumC(), scale, mat.Start(), mat.Ld(), Start(), _ld);
#else
  my_Domatcopy(trans, mat.DeepNumR(), mat.DeepNumC(), scale, mat.Start(), mat.Ld(), Start(), _ld);
#endif
}
