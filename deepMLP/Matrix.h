#pragma once

#if defined(HAVE_CONFIG_H)
#  include <config.h>
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

class Vector {
protected:
	const double *_start;
	int _len, _inc, _ld;

public:
	Vector() : _start(NULL), _len(0), _inc(1), _ld(0) { }

	Vector(const double * start, int len, int inc = 1, int ld = 0)
		:
	_start(start), _len(len), _inc(inc), _ld(ld) { }

	Vector(const std::vector<double> & vec) : _start(&vec[0]), _len(vec.size()), _inc(1), _ld(0) { }

	const double *Start() const { return _start; }
	int Len() const { return _len; }
	const double *End() const { return _start + _inc * _len; }
	int Inc() const { return _inc; }

	void operator++() {
		_start += _ld;
	}

	void operator--() {
		_start -= _ld;
	}

	VecScal operator-() const;

	const double & operator[](int i) const {
		if (i < 0) i += _len;
		assert (0 <= i && i < _len);
		return *(_start + i * _inc);
	}

	Vector SubVector(int begin, int end) const {		
		if (begin < 0) begin = _len - ~begin;
		if (end < 0) end = _len - ~end;
		assert (0 <= begin && begin <= end && end <= _len);
		return Vector(_start + begin, (end - begin), _inc, _ld);
	}

	Matrix AsMatrix(int numR, int numC) const;
	
	Matrix AsMatrix(int numR) const;
	
	template <class Visitor>
	void Visit(Visitor visitor) const {
		const double *p = _start;
		while (p != End()) {
			visitor(*p);
			p += _inc;
		}
	}
	
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

class Matrix {
protected:
	const double *_start;
	int _numR, _numC, _ld;
	bool _trans;

public:
	Matrix() : _start(NULL), _numR(0), _numC(0), _ld(0), _trans(false) { }

	Matrix(const double *start, int numR, int numC, int ld, bool trans)
		:
	_start(start), _numR(numR), _numC(numC), _ld(ld), _trans(trans) { }
	
	const double *Start() const { return _start; }
	const double *End() const { return Start() + (_numC * _ld); }
	int DeepNumR() const { return _numR; }
	int DeepNumC() const { return _numC; }
	int NumR() const { return _trans ? _numC : _numR; }
	int NumC() const { return _trans ? _numR : _numC; }
	int Ld() const { return _ld; }

	bool IsVec() const { return (_ld == _numR); }
	int VecLen() const { ASSERT (IsVec()); return _numR * _numC; }
	Vector Vec() const { return Vector(_start, VecLen(), 1, 0); }

	bool IsTrans() const { return _trans; }
	Matrix Trans() const { return Matrix(_start, _numR, _numC, _ld, !_trans); }

	const double & At(int r, int c) const {
		if (_trans) std::swap(r, c);
		if (r < 0) r += _numR;
		if (c < 0) c += _numC;
		assert (r >= 0 && r < _numR);
		assert (c >= 0 && c < _numC);
		return *(_start + c * _ld + r);
	}

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
	
	Vector GetCol(int c) const {
		if (_trans) return Trans().GetRow(c);

		if (c < 0) c += _numC;
		assert (0 <= c && c < _numC);

		return Vector(Start() + c * _ld, _numR, 1, _ld);
	}

	Vector GetRow(int r) const {
		if (_trans) return Trans().GetCol(r);

		if (r < 0) r += _numR;
		assert (0 <= r && r < _numR);

		return Vector(Start() + r, _numC, _ld, 1);
	}
	
	Matrix GetCols(int beginCol, int endCol) const {
		return SubMatrix(0, -1, beginCol, endCol);
	}

	Matrix GetRows(int beginRow, int endRow) const {
		return SubMatrix(beginRow, endRow, 0, -1);
	}

	MatScal operator-() const;

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

		return Matrix(_start + beginCol * _ld + beginRow, endRow - beginRow, endCol - beginCol, _ld, _trans);
	}
};

class MutableVector : public Vector {
	friend class MutableMatrix;
	MutableVector(const Vector & vec) : Vector(vec) { }

public:
	MutableVector() : Vector() { }

	MutableVector(double * start, int len, int inc = 1, int ld = 0) : Vector(start, len, inc, ld) { }

	MutableVector(std::vector<double> & vec) : Vector(vec) { }

	MutableVector(const MutableMatrix & m);

	double *Start() const { return const_cast<double*>(Vector::Start()); }
	
	double & operator[](int i) const {
		return const_cast<double &>(Vector::operator[](i));
	}

	MutableVector SubVector(int begin, int end) const {
		return Vector::SubVector(begin, end);
	}

	MutableMatrix AsMatrix(int numR, int numC) const;

	MutableMatrix AsMatrix(int numR) const;
	
	const MutableVector & operator=(const VecScaledSum & expr) const;
	
	const MutableVector & operator=(const VecScal & expr) const;
	
	const MutableVector & operator+=(const VecScal & expr) const;
	
	const MutableVector & operator-=(const VecScal & expr) const;
	
	const MutableVector & operator*=(double a) const;
	
	const MutableVector & operator/=(double a) const { return operator*=(1.0 / a); }

	void CopyFrom(const Vector & vec) const;
	
	template <class Mutator>
	void Apply(Mutator mut) const {
		double *p = Start();
		while (p != End()) {
			*p = mut(*p);
			p += _inc;
		}
	}
	
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
	
	template <class Func>
	void Replace(Func trans) const {
		double *p = Start();
		while (p != End()) {
			*p = trans();
			p += _inc;
		}
	}
	
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

	void Assign(double val) const {
		Replace([val] () { return val; });
	}
	
	template <class VMLFunc>
	void ApplyVML(VMLFunc Func, const Vector & arg) const {
		assert (_inc == 1 && arg.Inc() == 1);
		Func(_len, arg.Start(), Start());
	}

	template <class VMLFunc>
	void ApplyVML(VMLFunc Func, const Vector & arg1, const Vector & arg2) const {
		assert (_inc == 1 && arg1.Inc() == 1 && arg2.Inc() == 1);
		Func(_len, arg1.Start(), arg2.Start(), Start());
	}

	template <class VMLFunc>
	void ApplyVML(VMLFunc Func) const {
		assert (_inc == 1);
		Func(_len, Start(), Start());
	}

	void Axpby(const VecScal & expr, double a, double b) const;
};

class MutableMatrix : public Matrix {
	friend class MutableVector;
	MutableMatrix(const Matrix & mat) : Matrix(mat) { }

public:
	MutableMatrix() : Matrix() { }

	MutableMatrix(double *start, int numR, int numC, int ld, bool trans) : Matrix(start, numR, numC, ld, trans) { }

	MutableMatrix Trans() const { return MutableMatrix(Start(), _numR, _numC, _ld, !_trans); }

	double *Start() const { return const_cast<double*>(_start); }
	double *End() const { return Start() + (_numC * _ld); }
	
	double & At(int r, int c) const {
		return const_cast<double &>(Matrix::At(r,c));
	}
	
	MutableVector Vec() const {
		return MutableVector(Start(), VecLen(), 1, 0);
	}

	void CopyFrom(const Matrix & mat, double scale = 1.0) const;
	
	MutableVector GetRow(int r) const {
		return Matrix::GetRow(r);
	}

	MutableVector GetCol(int c) const {
		return Matrix::GetCol(c);
	}
	
	MutableMatrix GetCols(int beginCol, int endCol) const {
		return SubMatrix(0, -1, beginCol, endCol);
	}

	MutableMatrix GetRows(int beginRow, int endRow) const {
		return SubMatrix(beginRow, endRow, 0, -1);
	}

	MutableMatrix SubMatrix(int beginRow, int endRow, int beginCol, int endCol) const {
		return Matrix::SubMatrix(beginRow, endRow, beginCol, endCol);
	}
	
	const MutableMatrix & operator=(const MatMatMult & expr) const;

	const MutableMatrix & operator=(const MatScaledSum & expr) const;
	
	const MutableMatrix & operator=(const MatScal & expr) const;
	
	const MutableMatrix & operator+=(const MatScal & expr) const;

	const MutableMatrix & operator+=(const MatMatMult & expr) const;
	
	const MutableMatrix & operator-=(const MatScal & expr) const;

	const MutableMatrix & operator-=(const MatMatMult & expr) const;
	
	const MutableMatrix & operator*=(double a) const;

	const MutableMatrix & operator/=(double a) const { return operator*=(1.0 / a); }

	const MutableMatrix & Dgemm(double a, const Matrix & A, const Matrix & B, double b) const;

	void Assign(double val) const {
		if (IsVec()) Vec().Assign(val);
		else if (_trans) Trans().Assign(val);
		else for (int c = 0; c < _numC; ++c) {
			GetCol(c).Assign(val);
		}
	}
};

class AllocatingVector : public MutableVector {
	std::vector<double> _arr;
	void ResetStart() {
		_start = (_arr.size() > 0) ? &_arr[0] : NULL;
	}

public:
	AllocatingVector() : MutableVector() { }

	AllocatingVector(int len, double val = 0) : MutableVector(NULL, len), _arr(len, val) {
		ResetStart();
	}
	
	AllocatingVector(const Vector & vec) : MutableVector(NULL, vec.Len()), _arr(vec.Len()) {
		ResetStart();
		MutableVector::CopyFrom(vec);
	}

	AllocatingVector(const MutableVector & vec) : MutableVector(NULL, vec.Len()), _arr(vec.Len()) {
		ResetStart();
		MutableVector::CopyFrom(vec);
	}

	AllocatingVector(const AllocatingVector & vec) : MutableVector(vec), _arr(vec._arr) {
		ResetStart();
	}

	void Resize(int len) {
		_arr.resize(len);
		ResetStart();
		_len = len;
	}

	void CopyFrom(const Vector & vec) {
		Resize(vec.Len());
		MutableVector::CopyFrom(vec);
	}
	
	void Assign(int len, double val) {
		Resize(len);
		MutableVector::Assign(val);
	}

	void Assign(double val) {
		MutableVector::Assign(val);
	}

	void Swap(AllocatingVector & other) {		
		std::swap(_len, other._len);
		std::swap(_start, other._start);
		_arr.swap(other._arr);
	}
	
	const Vector & operator=(const Vector & other) {
		CopyFrom(other);
		return *this;
	}

	const Vector & operator=(const MutableVector & other) {
		CopyFrom(other);
		return *this;
	}

	const Vector & operator=(const AllocatingVector & other) {
		CopyFrom(other);
		return *this;
	}

	template <class Expr>
	const Vector & operator=(const Expr & expr) {
		Resize(expr.Len());
		return MutableVector::operator=(expr);
	}
};

class AllocatingMatrix : public MutableMatrix {
	std::vector<double> _arr;
	void ResetStart() {
		_start = (_arr.size() > 0) ? &_arr[0] : NULL;
	}

public:
	AllocatingMatrix() : MutableMatrix() { }
	
	AllocatingMatrix(int numR, int numC, double val = 0) : MutableMatrix(NULL, numR, numC, numR, false), _arr(_numR * _numC, val) {
		ResetStart();
	}
	
	AllocatingMatrix(const Matrix & mat) : MutableMatrix(NULL, mat.NumR(), mat.NumC(), mat.NumR(), false), _arr(_numR * _numC) {
		ResetStart();
		MutableMatrix::CopyFrom(mat);
	}

	AllocatingMatrix(const MutableMatrix & mat) : MutableMatrix(NULL, mat.NumR(), mat.NumC(), mat.NumR(), false), _arr(_numR * _numC) {
		ResetStart();
		MutableMatrix::CopyFrom(mat);
	}

	AllocatingMatrix(const AllocatingMatrix & mat) : MutableMatrix(mat), _arr(mat._arr) {
		ResetStart();
	}

	void Resize(int numR, int numC) {
		_arr.resize(numR * numC);
		_trans = false;
		ResetStart();
		_numR = numR;
		_numC = numC;
		_ld = _numR;
	}

	void Resize(const Matrix & A) { Resize(A.NumR(), A.NumC()); }

	void Assign(double val) {
		MutableMatrix::Assign(val);
	}

	void Assign(int numR, int numC, double val) {
		Resize(numR, numC);
		MutableMatrix::Assign(val);
	}

	void CopyFrom(const Matrix & mat) {
		Resize(mat);
		MutableMatrix::CopyFrom(mat);
	}
	
	void Swap(AllocatingMatrix & other) {		
		std::swap(_numR, other._numR);
		std::swap(_numC, other._numC);
		std::swap(_start, other._start);
		std::swap(_trans, other._trans);
		std::swap(_ld, other._ld);
		_arr.swap(other._arr);
	}
	
	const Matrix & operator=(const Matrix & other) {
		CopyFrom(other);
		return *this;
	}

	const Matrix & operator=(const MutableMatrix & other) {
		CopyFrom(other);
		return *this;
	}

	const Matrix & operator=(const AllocatingMatrix & other) {
		CopyFrom(other);
		return *this;
	}

	template <class Expr>
	const Matrix & operator=(const Expr & expr) {
		Resize(expr.NumR(), expr.NumC());
		return MutableMatrix::operator=(expr);
	}
};

struct MatScal {
	const Matrix A;
	const double a;

	MatScal(const Matrix & A, const double a = 1.0) : A(A), a(a) { }
	MatScal operator-() const { return MatScal(A, -a); }

	int NumR() const { return A.NumR(); }
	int NumC() const { return A.NumC(); }
};

MatScal operator*(double b, const MatScal & ms);
MatScal operator*(const MatScal & ms, double b);

struct MatScaledSum {
	const Matrix A, B;
	const double a, b;

	MatScaledSum(const MatScal & ax, const MatScal & by) : A(ax.A), a(ax.a), B(by.A), b(by.a) {
		assert (A.NumC() == B.NumC() && A.NumR() == B.NumR());
	}

	int NumR() const { return A.NumR(); }
	int NumC() const { return A.NumC(); }
};

MatScaledSum operator+(const MatScal & ax, const MatScal & by);
MatScaledSum operator-(const MatScal & ax, const MatScal & by);
MatScaledSum operator*(double c, const MatScaledSum & mss);
MatScaledSum operator*(const MatScaledSum & mss, double c);

struct MatMatMult {
	const Matrix A;
	const Matrix B;
	const double a;

	MatMatMult(const Matrix & A, const Matrix & B, double a) : A(A), B(B), a(a) { 
		assert (A.NumC() == B.NumR());
	}

	int NumR() const { return A.NumR(); }
	int NumC() const { return B.NumC(); }
};

MatMatMult operator*(const MatScal & A, const MatScal & B);
MatMatMult operator*(double b, const MatMatMult & mmm);
MatMatMult operator*(const MatMatMult & mmm, double b);

struct VecScal {
	const Vector v;
	const double a;

	VecScal(const Vector & v, const double a = 1.0) : v(v), a(a) { }
	VecScal operator-() const { return VecScal(v, -a); }

	int Len() const { return v.Len(); }
};

VecScal operator*(double b, const VecScal & ms);
VecScal operator*(const VecScal & ms, double b);
double operator*(const VecScal & ax, const VecScal & by);

struct VecScaledSum {
	const Vector x, y;
	const double a, b;

	VecScaledSum(const VecScal & ax, const VecScal & by) : x(ax.v), a(ax.a), y(by.v), b(by.a) {
		assert(x.Len() == y.Len());
	}

	int Len() const { return x.Len(); }
};

VecScaledSum operator+(const VecScal & ax, const VecScal & by);
VecScaledSum operator-(const VecScal & ax, const VecScal & by);
VecScaledSum operator*(double c, const VecScaledSum & mss);
VecScaledSum operator*(const VecScaledSum & mss, double c);

