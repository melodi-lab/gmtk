#pragma once

/*
 * Written by Galen Andrew gmandrew@uw.edu
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#if defined(HAVE_CONFIG_H)
#  include <config.h>
#endif
#if HAVE_INTTYPES_H
  // The ISO C99 standard specifies that the macros in inttypes.h must
  //  only be defined if explicitly requested. 
#  ifndef __STDC_FORMAT_MACROS
#    define __STDC_FORMAT_MACROS 1
#  endif
#  include <inttypes.h>
#endif
#if HAVE_STDINT_H
#  include <stdint.h>
#endif

#include <vector>
#include <assert.h>
#include <string>
#include <stdio.h>

#if defined(HAVE_MKL)
#  include "mkl.h"
#endif

#include "Globals.h"

class Matrix;
class MutableMatrix;
struct MatScal;
struct MatScaledSum;
struct MatMatMult;
struct VecScal;
struct VecScaledSum;

/**********************
The Matrix and Vector classes were created to allow matrix/vector
operations to be performed simply and intuitively with arithmetic
operators, while being backed with a fast library like Intel MKL.

The major classes fall in two hierarchies:
  Vector < MutableVector < AllocatingVector
and
  Matrix < MutableMatrix < AllocatingMatrix

The "allocating" versions create space for the data on
construction and destroy it on destruction. The other versions
are "views", containing only pointers to the data. The difference
between them is that the "mutable" classes allow write operations
on the data, while the basic versions treat the data as readonly.

Matrices are stored in column order, and maintain their transpose
state, as well as the "ld" distance between starts of adjacent
columns. Vectors also maintain the vector increment.

Common math operators (multiplication by a scalar, dot-product,
matrix-vector product, matrix-matrix product, addition, etc.) are
implemented using a handful of utility classes and operator
overloading.

NOTE ABOUT INDEXING.
Several operators like Vector::operator[](int i) and Matrix::
At(r,c) index into the structures. If any of the arguments i are
negative, the result is the same as if the value (~i + n - 1)
were used, where n is the relevant size (vector length, matrix
number of rows or columns). This way for example operator[](~k)
indexes the kth-from-last element of a Vector.

Negative arguments are interpreted slightly differently for
operations that specify a range (SubVector, SubMatrix, GetCols).
Here, if any argument i is negative, it is interpreted
as (n - ~a). So for example to get the last two elements of a
Vector, use Vector::SubVector(~2, ~0).

 **********************/

// Base class for Vectors. Contains read-only pointer to data.
// Maintains increment between vector elements and increment
// between adjacent vectors in a larger matrix, if applicable
class Vector {
protected:
  const double *_start;
  int _len, _inc, _ld;

public:
	// construct a null vector
  Vector() : _start(NULL), _len(0), _inc(1), _ld(0) { }

	// construct a vector with the given data, size, and increment
	// ld is optional, but if provided allows the vector to skip
	// through rows or columns of a matrix it is part of
  Vector(const double * start, int len, int inc = 1, int ld = 0)
    :
  _start(start), _len(len), _inc(inc), _ld(ld) { }

	// copy constructor
  Vector(const std::vector<double> & vec) : _start(&vec[0]), _len(vec.size()), _inc(1), _ld(0) { }

	// get start pointer
  const double *Start() const { return _start; }

	// get vector length
  int Len() const { return _len; }

	// get end pointer (one past last entry of vector)
  const double *End() const { return _start + _inc * _len; }

	// get distance between elements
  int Inc() const { return _inc; }

	// advance to next row or column, if vector is part of a
	// matrix (if ld was provided to constructor)
  void operator++() {
    _start += _ld;
  }

	// go back to last row or column, if vector is part of a
	// matrix (if ld was provided to constructor)
  void operator--() {
    _start -= _ld;
  }

	// negate vector. This is an implicit operation, no data
	// is copied until the resulting VecScal is assigned to
	// another vector
  VecScal operator-() const;

	// get the ith element. See NOTE ABOUT INDEXING above
  const double & operator[](int i) const {
    if (i < 0) i += _len;
    assert (0 <= i && i < _len);
    return *(_start + i * _inc);
  }

	// Returns a new vector view that begins with the given
	// element position and is of length (end-begin)
	// See NOTE ABOUT INDEXING above
  Vector SubVector(int begin, int end) const {		
    if (begin < 0) begin = _len - ~begin;
    if (end < 0) end = _len - ~end;
    assert (0 <= begin && begin <= end && end <= _len);
    return Vector(_start + begin, (end - begin), _inc, _ld);
  }

	// Interpret the vector data as a matrix with the given
	// number of rows and columns. It is required that
	// numR * numC == _len
  Matrix AsMatrix(int numR, int numC) const;

	// Interpret the vector data as a matrix with the given
	// number of rows. It is required that _len % numR == 0
  Matrix AsMatrix(int numR) const;

	// Apply a visitor object to each element
	// A Visitor is assumed to implement operator()(double x)
	// Very handy with lambda expressions!
  template <class Visitor>
  void Visit(Visitor visitor) const {
    const double *p = _start;
    while (p != End()) {
      visitor(*p);
      p += _inc;
    }
  }

	// Apply a visitor object to each pair of elements, one
	// from *this* vector and one from the argument vector.
	// It is assumed that the two vectors have the same length.
	// A Visitor is assumed to implement operator()(double x, double y)
	// where the first argument (x) is an element of *this* vector
	// and the second is an element of the argument vector.
	// Very handy with lambda expressions!
  template <class Visitor>
  void Visit(Visitor visitor, const Vector & arg) const {
    assert (arg.Len() == Len());
    const double *p = _start;
    const double *pA = arg.Start();
    while (p != End()) {
      visitor(*p, *pA);
      p += _inc;
      pA += arg.Inc();
    }
  }

	// Print a representation of the vector to the console.
  void Print(int prec = 3) const {
    if (prec < 0) prec = 0;
    else if (prec > 5) prec = 5;
    int len = 4 + prec;
    std::string s = "%-l.pf";
    s[2] = '0' + len;
    s[4] = '0' + prec;
    for (int i = 0; i < Len(); ++i) {
      printf(s.c_str(), operator[](i));
    }
    printf("\n");
  }
};

// Base class for Matrices. Contains read-only pointer to data.
// Matrices are stored in column order. Maintains whether the
// matrix is transposed, and the distance between first elements
// of each column
class Matrix {
protected:
  const double *_start;
  int _numR, _numC, _ld;
  bool _trans;

public:
	// create an empty matrix
  Matrix() : _start(NULL), _numR(0), _numC(0), _ld(0), _trans(false) { }

	// creates a matrix with the given pointer to data (start),
	// number of rows and columns, distance between starts of adjacent
	// columns ld, and possibly transposed
  Matrix(const double *start, int numR, int numC, int ld, bool trans)
    :
  _start(start), _numR(numR), _numC(numC), _ld(ld), _trans(trans) { }

	// pointer to start of matrix
  const double *Start() const { return _start; }

	// pointer to end of matrix (location which would be start of
	// next column after last)
  const double *End() const { return Start() + (_numC * _ld); }

	// Number of underlying rows, not considering that matrix may
	// be carrying a transpose tag
  int DeepNumR() const { return _numR; }

	// Number of underlying columns, not considering that matrix may
	// be carrying a transpose tag
  int DeepNumC() const { return _numC; }

	// Number of rows
  int NumR() const { return _trans ? _numC : _numR; }

	// Number of columns
  int NumC() const { return _trans ? _numR : _numC; }

	// distance between initial elements of each column
  int Ld() const { return _ld; }

	// returns true iff the matrix is contiguous in memory
	// that is, whether ld==numR, or each column starts immediately
	// following the last
  bool IsVec() const { return (_ld == _numR); }

	// Total number of elements IF IsVec() is true
	// otherwise causes program abort
  int VecLen() const { assert(IsVec()); return _numR * _numC; }

	// Returns vector representation of matrix IF IsVec() is true
	// otherwise causes program abort
  Vector Vec() const { return Vector(_start, VecLen(), 1, 0); }

	// returns true iff the matrix is transposed
  bool IsTrans() const { return _trans; }

	// returns a matrix pointing to the same data, but transposed
  Matrix Trans() const { return Matrix(_start, _numR, _numC, _ld, !_trans); }

	// Return the element at the given row and column
	// See NOTE ABOUT INDEXING above
  const double & At(int r, int c) const {
    if (_trans) std::swap(r, c);
    if (r < 0) r += _numR;
    if (c < 0) c += _numC;
    assert (r >= 0 && r < _numR);
    assert (c >= 0 && c < _numC);
    return *(_start + c * _ld + r);
  }

	// Prints a representation of the Matrix to the console
  void Print(int prec = 3) const {
    if (prec < 0) prec = 0;
    else if (prec > 5) prec = 5;
    int len = 4 + prec;
    std::string s = "%-l.pf";
    s[2] = '0' + len;
    s[4] = '0' + prec;
    for (int r = 0; r < NumR(); ++r) {
      for (int c = 0; c < NumC(); ++c) {
        printf(s.c_str(), At(r,c));
      }
      printf("\n");
    }
  }

	// Return a Vector representation of the given column
	// See NOTE ABOUT INDEXING above
  Vector GetCol(int c) const {
    if (_trans) return Trans().GetRow(c);

    if (c < 0) c += _numC;
    assert (0 <= c && c < _numC);

    return Vector(Start() + c * _ld, _numR, 1, _ld);
  }

	// Return a Vector representation of the given row
	// See NOTE ABOUT INDEXING above
  Vector GetRow(int r) const {
    if (_trans) return Trans().GetCol(r);

    if (r < 0) r += _numR;
    assert (0 <= r && r < _numR);

    return Vector(Start() + r, _numC, _ld, 1);
  }

	// Return a Matrix consisting of the columns from beginCol
	// to endCol (not inclusive: the result will have
	// endCol-beginCol total columns)
	// See NOTE ABOUT INDEXING above
  Matrix GetCols(int beginCol, int endCol) const {
    return SubMatrix(0, -1, beginCol, endCol);
  }

	// Return a Matrix consisting of the rows from beginRow
	// to endRow (not inclusive: the result will have
	// endRow-beginRow total rows)
	// See NOTE ABOUT INDEXING above
  Matrix GetRows(int beginRow, int endRow) const {
    return SubMatrix(beginRow, endRow, 0, -1);
  }

	// negate the matrix. This is an implicit operation, no data
	// is copied until the resulting MatScal is assigned to
	// another matrix
  MatScal operator-() const;

	// Return a Matrix consisting of the rows from beginRow
	// to endRow, and the columns from beginCol to endCol
	// (not inclusive: the result will have (endRow-beginRow)
	// total rows and (endCol-beginCol) total columns)
	// See NOTE ABOUT INDEXING above
  Matrix SubMatrix(int beginRow, int endRow, int beginCol, int endCol) const {
    if (_trans) {
      std::swap(beginRow, beginCol);
      std::swap(endRow, endCol);
    }

    if (beginRow < 0) beginRow = _numR - ~beginRow;
    if (endRow < 0) endRow = _numR - ~endRow;
    assert (0 <= beginRow && beginRow <= endRow && endRow <= _numR);

    if (beginCol < 0) beginCol = _numC - ~beginCol;
    if (endCol < 0) endCol = _numC - ~endCol;
    assert (0 <= beginCol && beginCol <= endCol && endCol <= _numC);

    return Matrix(_start + (int64_t)beginCol * _ld + beginRow, endRow - beginRow, endCol - beginCol, _ld, _trans);
  }
};

// Subclass of Vector that adds methods that can alter the
// the underlying data
class MutableVector : public Vector {
  friend class MutableMatrix;
  MutableVector(const Vector & vec) : Vector(vec) { }

public:
	// Make an empty vector
  MutableVector() : Vector() { }

	// construct a vector with the given data, size, and increment
	// ld is optional, but if provided allows the vector to skip
	// through rows or columns of a matrix it is part of
  MutableVector(double * start, int len, int inc = 1, int ld = 0) : Vector(start, len, inc, ld) { }

	// create a vector pointing to the data in this std::vector
  MutableVector(std::vector<double> & vec) : Vector(vec) { }

	// This appears to have no implementation and is not being used
	//  MutableVector(const MutableMatrix & m);

	// Get the pointer to the beginning of the data
  double *Start() const { return const_cast<double*>(Vector::Start()); }

	// get the ith element. See NOTE ABOUT INDEXING above
  double & operator[](int i) const {
    return const_cast<double &>(Vector::operator[](i));
  }

	// Returns a new vector view that begins with the given
	// element position and is of length (end-begin)
	// See NOTE ABOUT INDEXING above
  MutableVector SubVector(int begin, int end) const {
    return Vector::SubVector(begin, end);
  }

	// Interpret the vector data as a matrix with the given
	// number of rows and columns. It is required that
	// numR * numC == _len
  MutableMatrix AsMatrix(int numR, int numC) const;

	// Interpret the vector data as a matrix with the given
	// number of rows. It is required that _len % numR == 0
  MutableMatrix AsMatrix(int numR) const;

	// Overwrite this vector's data with the results of evaluating
	// the linear combination of two other vectors
  const MutableVector & operator=(const VecScaledSum & expr) const;

	// Overwrite this vector's data with the results of multiplying
	// another vector with a scalar
  const MutableVector & operator=(const VecScal & expr) const;

	// Overwrite this vector's data by adding the results of
	// multiplying another vector with a scalar
  const MutableVector & operator+=(const VecScal & expr) const;

	// Overwrite this vector's data by subtracting the results of
	// multiplying another vector with a scalar
  const MutableVector & operator-=(const VecScal & expr) const;

	// Overwrite this vector's data by with the results of
	// multiplying it with a scalar
  const MutableVector & operator*=(double a) const;

	// Overwrite this vector's data by with the results of
	// dividing it by a scalar
  const MutableVector & operator/=(double a) const { return operator*=(1.0 / a); }

	// Overwrite this vector's data with a copy of vec's data
	// It is a precondition that the vector lengths are the same
  void CopyFrom(const Vector & vec) const;

	// Overwrite this vector's data by applying a mutator object to each element
	// A Mutator is assumed to implement operator()(double x)->double
	// Very handy with lambda expressions!
  template <class Mutator>
  void Apply(Mutator mut) const {
    double *p = Start();
    while (p != End()) {
      *p = mut(*p);
      p += _inc;
    }
  }

	// Overwrite this vector's data by applying a mutator object
	// to each pair of elements, one from *this* vector and one from
	// the argument vector.
	// It is assumed that the two vectors have the same length.
	// A Mutator is assumed to implement
	// operator()(double x, double y)->double
	// where the first argument (x) is an element of *this* vector
	// and the second is an element of the argument vector.
	// Very handy with lambda expressions!
  template <class Mutator>
  void Apply(Mutator mut, const Vector & arg) const {
    assert (arg.Len() == Len());
    double *p = Start();
    const double *pA = arg.Start();
    while (p != End()) {
      *p = mut(*p, *pA);
      p += _inc;
      pA += arg.Inc();
    }
  }

	// Overwrite this vector's data by applying a Func object to each element
	// A Func is assumed to implement operator()()->double
	// Note the difference with Apply is that the current data is not used
	// Very handy with lambda expressions!
  template <class Func>
  void Replace(Func trans) const {
    double *p = Start();
    while (p != End()) {
      *p = trans();
      p += _inc;
    }
  }

	// Overwrite this vector's data by applying a Func object to each element
	// of another vector arg. A Func is assumed to implement
	// operator()(double)->double
	// Note the difference with Apply is that the current data is not used
	// Very handy with lambda expressions!
  template <class Func>
  void Replace(Func trans, const Vector & arg) const {
    assert (arg.Len() == Len());
    double *p = Start();
    const double *pA = arg.Start();
    while (p != End()) {
      *p = trans(*pA);
      p += _inc;
      pA += arg.Inc();
    }
  }

	// Overwrite this vector's data by applying a Func object to each pair of
	// elements one from each of the two argument vectors. A Func is assumed
	// to implement operator()(double x, double y)->double
	// where x is the value from the first arg, and y from the second.
	// Note the difference with Apply is that the current data is not used
	// Very handy with lambda expressions!
  template <class Func>
  void Replace(Func trans, const Vector & arg, const Vector & arg2) const {
    assert (arg.Len() == Len());
    double *p = Start();
    const double *pA = arg.Start(), *pA2 = arg2.Start();
    while (p != End()) {
      *p = trans(*pA, *pA2);
      p += _inc;
      pA += arg.Inc();
      pA2 += arg2.Inc();
    }
  }

	// Overwrite every element of this vector's data with the given value
  void Assign(double val) const {
    Replace([val] () { return val; });
  }

	// Overwrite this vector's data by applying
	// an Intel Vector Matrix Library function taking a single vector
	// argument (such as vdexp) to the argument vector
  template <class VMLFunc>
  void ApplyVML(VMLFunc Func, const Vector & arg) const {
    assert (_inc == 1 && arg.Inc() == 1);
    Func(_len, arg.Start(), Start());
  }

	// Overwrite this vector's data by applying
	// an Intel Vector Matrix Library function taking two vector
	// arguments (such as vdadd) to the argument vectors
  template <class VMLFunc>
  void ApplyVML(VMLFunc Func, const Vector & arg1, const Vector & arg2) const {
    assert (_inc == 1 && arg1.Inc() == 1 && arg2.Inc() == 1);
    Func(_len, arg1.Start(), arg2.Start(), Start());
  }

	// Overwrite this vector's data by applying
	// an Intel Vector Matrix Library function taking a single vector
	// argument (such as vdexp) to the current data
  template <class VMLFunc>
  void ApplyVML(VMLFunc Func) const {
    assert (_inc == 1);
    Func(_len, Start(), Start());
  }

	// Overwrite this vector's data as a linear combination of
	// itself with another vector
	// this = a * expr + b * this
  void Axpby(const VecScal & expr, double a, double b) const;
};

// Subclass of Matrix that adds methods that can alter the
// the underlying data
class MutableMatrix : public Matrix {
  friend class MutableVector;
  MutableMatrix(const Matrix & mat) : Matrix(mat) { }

public:
	// create an empty matrix
  MutableMatrix() : Matrix() { }

	// creates a matrix with the given pointer to data (start),
	// number of rows and columns, distance between starts of adjacent
	// columns ld, and possibly transposed
  MutableMatrix(double *start, int numR, int numC, int ld, bool trans) : Matrix(start, numR, numC, ld, trans) { }

	// returns a matrix pointing to the same data, but transposed
  MutableMatrix Trans() const { return MutableMatrix(Start(), _numR, _numC, _ld, !_trans); }

	// pointer to start of matrix
  double *Start() const { return const_cast<double*>(_start); }

	// pointer to end of matrix (location which would be start of
	// next column after last)
  double *End() const { return Start() + (_numC * _ld); }

	// Return the element at the given row and column
	// See NOTE ABOUT INDEXING above
  double & At(int r, int c) const {
    return const_cast<double &>(Matrix::At(r,c));
  }

	// Returns vector representation of matrix IF IsVec() is true
	// otherwise causes program abort
  MutableVector Vec() const {
    return MutableVector(Start(), VecLen(), 1, 0);
  }

	// Overwrite this matrix's data with a copy of mat's data (possibly scaled)
  void CopyFrom(const Matrix & mat, double scale = 1.0) const;

	// Return a MutableVector representation of the given row
	// See NOTE ABOUT INDEXING above
  MutableVector GetRow(int r) const {
    return Matrix::GetRow(r);
  }

	// Return a MutableVector representation of the given column
	// See NOTE ABOUT INDEXING above
  MutableVector GetCol(int c) const {
    return Matrix::GetCol(c);
  }

	// Return a Matrix consisting of the columns from beginCol
	// to endCol (not inclusive: the result will have
	// endCol-beginCol total columns)
	// See NOTE ABOUT INDEXING above
  MutableMatrix GetCols(int beginCol, int endCol) const {
    return SubMatrix(0, -1, beginCol, endCol);
  }

	// Return a Matrix consisting of the rows from beginRow
	// to endRow (not inclusive: the result will have
	// endRow-beginRow total rows)
	// See NOTE ABOUT INDEXING above
  MutableMatrix GetRows(int beginRow, int endRow) const {
    return SubMatrix(beginRow, endRow, 0, -1);
  }

	// Return a Matrix consisting of the rows from beginRow
	// to endRow, and the columns from beginCol to endCol
	// (not inclusive: the result will have (endRow-beginRow)
	// total rows and (endCol-beginCol) total columns)
	// See NOTE ABOUT INDEXING above
  MutableMatrix SubMatrix(int beginRow, int endRow, int beginCol, int endCol) const {
    return Matrix::SubMatrix(beginRow, endRow, beginCol, endCol);
  }

	// Overwrite this matrix's data with the result of the
	// matrix/matrix multiplication op specified in expr
  const MutableMatrix & operator=(const MatMatMult & expr) const;

	// Overwrite this matrix's data with the results of evaluating
	// the linear combination of two other matrices
  const MutableMatrix & operator=(const MatScaledSum & expr) const;

	// Overwrite this matrix's data with the result of the given
	// scaling operation on another matrix
  const MutableMatrix & operator=(const MatScal & expr) const;

	// Overwrite this matrix's data by adding a multiple of another
	// matrix (supplied in expr) to the current data
  const MutableMatrix & operator+=(const MatScal & expr) const;

	// Overwrite this matrix's data by adding a product of two other
	// matrices (supplied in expr) to the current data
  const MutableMatrix & operator+=(const MatMatMult & expr) const;

	// Overwrite this matrix's data by subtracting a multiple of another
	// matrix (supplied in expr) to the current data
  const MutableMatrix & operator-=(const MatScal & expr) const;

	// Overwrite this matrix's data by subtracting a product of two other
	// matrices (supplied in expr) to the current data
  const MutableMatrix & operator-=(const MatMatMult & expr) const;

	// Overwrite this matrix's data with the results of
	// multiplying it by a scalar
  const MutableMatrix & operator*=(double a) const;

	// Overwrite this matrix's data with the results of
	// dividing it by a scalar
  const MutableMatrix & operator/=(double a) const { return operator*=(1.0 / a); }

	// The complete Dgemm specification
	// If this matrix represents the matrix X, overwrite this matrix's
	// data with the result of evaluating aAB+bX
	// Note transpose is handled by calling Trans() on A or B
  const MutableMatrix & Dgemm(double a, const Matrix & A, const Matrix & B, double b) const;

	// Overwrite every element of the matrix with the given value
  void Assign(double val) const {
    if (IsVec()) Vec().Assign(val);
    else if (_trans) Trans().Assign(val);
    else for (int c = 0; c < _numC; ++c) {
      GetCol(c).Assign(val);
    }
  }
};

// Subclass of MutableVector that allocates the underlying array at
// construction and destroyes it at destruction
class AllocatingVector : public MutableVector {
  std::vector<double> _arr;
  void ResetStart() {
    _start = (_arr.size() > 0) ? &_arr[0] : NULL;
  }

public:
	// Make an empty (size-0) vector
  AllocatingVector() : MutableVector() { }

	// construct a vector of length len
	// all values are initialized to val
  AllocatingVector(int len, double val = 0) : MutableVector(NULL, len), _arr(len, val) {
    ResetStart();
  }

	// construct a new vector with length and initial values taken from vec
  AllocatingVector(const Vector & vec) : MutableVector(NULL, vec.Len()), _arr(vec.Len()) {
    ResetStart();
    MutableVector::CopyFrom(vec);
  }

	// construct a new vector with length and initial values taken from vec
  AllocatingVector(const MutableVector & vec) : MutableVector(NULL, vec.Len()), _arr(vec.Len()) {
    ResetStart();
    MutableVector::CopyFrom(vec);
  }

	// construct a new vector with length and initial values taken from vec
  AllocatingVector(const AllocatingVector & vec) : MutableVector(vec), _arr(vec._arr) {
    ResetStart();
  }

	// Grow or shrink the vector
	// Because it is backed with std::vector, if the new size is smaller
	// the old values at those positions are guaranteed to remain
	// if the new size is larger, there are no guarantees about the data
  void Resize(int len) {
    _arr.resize(len);
    ResetStart();
    _len = len;
  }

	// Change the size if necessary and copy the data from vec
  void CopyFrom(const Vector & vec) {
    Resize(vec.Len());
    MutableVector::CopyFrom(vec);
  }

	// Resize to length len and assign all elements the value val
  void Assign(int len, double val) {
    Resize(len);
    MutableVector::Assign(val);
  }

	// assign all elements the value val
  void Assign(double val) {
    MutableVector::Assign(val);
  }

	// Exchange underlying data with another AllocatingVector
	// without any allocation or copying
  void Swap(AllocatingVector & other) {		
    std::swap(_len, other._len);
    std::swap(_start, other._start);
    _arr.swap(other._arr);
  }

	// Change the size if necessary and copy the data from vec
  const Vector & operator=(const Vector & other) {
    CopyFrom(other);
    return *this;
  }

	// Change the size if necessary and copy the data from vec
  const Vector & operator=(const MutableVector & other) {
    CopyFrom(other);
    return *this;
  }

	// Change the size if necessary and copy the data from vec
  const Vector & operator=(const AllocatingVector & other) {
    CopyFrom(other);
    return *this;
  }

	// Resize the vector if necessary to match the size of the
	// result of evaluating the expression expr, then overwrite
	// the vector's data with the results of the evaluation
  template <class Expr>
  const Vector & operator=(const Expr & expr) {
    Resize(expr.Len());
    return MutableVector::operator=(expr);
  }
};

// Subclass of MutableMatrix that allocates the underlying array at
// construction and destroyes it at destruction
class AllocatingMatrix : public MutableMatrix {
  std::vector<double> _arr;
  void ResetStart() {
    _start = (_arr.size() > 0) ? &_arr[0] : NULL;
  }

public:
	// creates an empty (0x0) matrix
  AllocatingMatrix() : MutableMatrix() { }

	// creates a matrix with the supplied number of rows and columns
	// and initializes the value of each element to val
  AllocatingMatrix(int numR, int numC, double val = 0) : MutableMatrix(NULL, numR, numC, numR, false), _arr(_numR * _numC, val) {
    ResetStart();
  }

	// creates a matrix whose dimensions are the same as mat
	// and initializes the data by copying the values of mat
  AllocatingMatrix(const Matrix & mat) : MutableMatrix(NULL, mat.NumR(), mat.NumC(), mat.NumR(), false), _arr(_numR * _numC) {
    ResetStart();
    MutableMatrix::CopyFrom(mat);
  }

	// creates a matrix whose dimensions are the same as mat
	// and initializes the data by copying the values of mat
  AllocatingMatrix(const MutableMatrix & mat) : MutableMatrix(NULL, mat.NumR(), mat.NumC(), mat.NumR(), false), _arr(_numR * _numC) {
    ResetStart();
    MutableMatrix::CopyFrom(mat);
  }

	// creates a matrix whose dimensions are the same as mat
	// and initializes the data by copying the values of mat
  AllocatingMatrix(const AllocatingMatrix & mat) : MutableMatrix(mat), _arr(mat._arr) {
    ResetStart();
  }

	// Change the dimensions of the matrix
	// there are no guarantees about the values of the matrix entries
	// after calling Resize
  void Resize(int numR, int numC) {
    _arr.resize(numR * numC);
    _trans = false;
    ResetStart();
    _numR = numR;
    _numC = numC;
    _ld = _numR;
  }

	// Change the dimensions of the matrix to match the dimensions of A
  void Resize(const Matrix & A) { Resize(A.NumR(), A.NumC()); }

	// Overwrite all entries of the matrix with val
  void Assign(double val) {
    MutableMatrix::Assign(val);
  }

	// Change the dimensions of the matrix and set all values to val
  void Assign(int numR, int numC, double val) {
    Resize(numR, numC);
    MutableMatrix::Assign(val);
  }

	// Change the dimensions of the matrix (if necessary) to match
	// those of mat, then copy the data from mat, including an optional
	// scaling factor
  void CopyFrom(const Matrix & mat, double scale = 1.0) {
    Resize(mat);
    MutableMatrix::CopyFrom(mat, scale);
  }

	// Exchange underlying data with another AllocatingMatrix
	// without any allocation or copying
  void Swap(AllocatingMatrix & other) {		
    std::swap(_numR, other._numR);
    std::swap(_numC, other._numC);
    std::swap(_start, other._start);
    std::swap(_trans, other._trans);
    std::swap(_ld, other._ld);
    _arr.swap(other._arr);
  }

	// Change the dimensions of the matrix (if necessary) to match
	// those of mat, then copy the data from mat
  const Matrix & operator=(const Matrix & other) {
    CopyFrom(other);
    return *this;
  }

	// Change the dimensions of the matrix (if necessary) to match
	// those of mat, then copy the data from mat
  const Matrix & operator=(const MutableMatrix & other) {
    CopyFrom(other);
    return *this;
  }

	// Change the dimensions of the matrix (if necessary) to match
	// those of mat, then copy the data from mat
  const Matrix & operator=(const AllocatingMatrix & other) {
    CopyFrom(other);
    return *this;
  }

	// Change the dimensions of the matrix (if necessary) to match
	// the dimensions of the result of evaluating expr, then overwrite
	// the data of the matrix with the results of the evaluation
  template <class Expr>
  const Matrix & operator=(const Expr & expr) {
    Resize(expr.NumR(), expr.NumC());
    return MutableMatrix::operator=(expr);
  }
};

// structure to store the parts of the expression a * A
// where A is a matrix and a is a scalar
// The expression is evaluated if the resulting MatScal is
// assigned to a matrix object
struct MatScal {
  const Matrix A;
  const double a;

  MatScal(const Matrix & A, const double a = 1.0) : A(A), a(a) { }
  MatScal operator-() const { return MatScal(A, -a); }

	// number of rows in the result
  int NumR() const { return A.NumR(); }

	// number of columns in the result
  int NumC() const { return A.NumC(); }
};

MatScal operator*(double b, const MatScal & ms);
MatScal operator*(const MatScal & ms, double b);


// structure to store the parts of the expression a * A + b * B
// where A and B are matrices and a and b are scalars
// The expression is evaluated if the resulting MatScaledSum is
// assigned to a matrix object
struct MatScaledSum {
  const Matrix A, B;
  const double a, b;

  MatScaledSum(const MatScal & ax, const MatScal & by) : A(ax.A), B(by.A), a(ax.a), b(by.a) {
    assert (A.NumC() == B.NumC() && A.NumR() == B.NumR());
  }

	// number of rows in the result
  int NumR() const { return A.NumR(); }

	// number of columns in the result
  int NumC() const { return A.NumC(); }
};

MatScaledSum operator+(const MatScal & ax, const MatScal & by);
MatScaledSum operator-(const MatScal & ax, const MatScal & by);
MatScaledSum operator*(double c, const MatScaledSum & mss);
MatScaledSum operator*(const MatScaledSum & mss, double c);


// structure to store the parts of the expression a * A * B
// where A and B are matrices and a is a scalar scalars
// The expression is evaluated if the resulting MatMatMult
// assigned to a matrix object
struct MatMatMult {
  const Matrix A;
  const Matrix B;
  const double a;

  MatMatMult(const Matrix & A, const Matrix & B, double a) : A(A), B(B), a(a) { 
    assert (A.NumC() == B.NumR());
  }

	// number of rows in the result
  int NumR() const { return A.NumR(); }

	// number of columns in the result
  int NumC() const { return B.NumC(); }
};

MatMatMult operator*(const MatScal & A, const MatScal & B);
MatMatMult operator*(double b, const MatMatMult & mmm);
MatMatMult operator*(const MatMatMult & mmm, double b);


// structure to store the parts of the expression a * v
// where v is a vector and a is a scalar
// The expression is evaluated if the resulting VecScal is
// assigned to a vector object
struct VecScal {
  const Vector v;
  const double a;

  VecScal(const Vector & v, const double a = 1.0) : v(v), a(a) { }
  VecScal operator-() const { return VecScal(v, -a); }

	// the length of the result 
  int Len() const { return v.Len(); }
};

VecScal operator*(double b, const VecScal & ms);
VecScal operator*(const VecScal & ms, double b);
double operator*(const VecScal & ax, const VecScal & by);


// structure to store the parts of the expression a * x + b * y
// where x and y are vectors and a and b are scalars
// The expression is evaluated if the resulting VecScal is
// assigned to a vector object
struct VecScaledSum {
  const Vector x, y;
  const double a, b;

  VecScaledSum(const VecScal & ax, const VecScal & by) : x(ax.v), y(by.v), a(ax.a), b(by.a) {
    assert(x.Len() == y.Len());
  }

	// the length of the result 
  int Len() const { return x.Len(); }
};

VecScaledSum operator+(const VecScal & ax, const VecScal & by);
VecScaledSum operator-(const VecScal & ax, const VecScal & by);
VecScaledSum operator*(double c, const VecScaledSum & mss);
VecScaledSum operator*(const VecScaledSum & mss, double c);

