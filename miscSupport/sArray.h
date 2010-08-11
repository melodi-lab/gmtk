//
//
// sArray: a Simple Array class that is basically a glorified type*.
// 
// This class is NOT meant for protection. It is mean to have type*
// (such as an int*, float*, double*, etc.)  and where the underlying
// pointer is easy to access for low level loops. But this class also
// provides for convenient managing of lengths, allocation,
// deallocation, and resizing. It also includes a simple sort routine
// that will work as long as the object T has an operator<(), copy
// constructor and operator=().
//
// ******************* NOTE *****************************************
// There is no COPY CONSTRUCTOR for this object, so all copies are
// bit-wise copies. This means that sArrays work *COMPLETELY*
// differently then say STL vector types, etc.. sArrays are meant to
// be used when you wish not to create lots of temporaries during
// construction/destruction in routine return and call. If you include
// sArrays in a class, you will need to explicitly manage any
// copying/constructing.
// ******************* NOTE *****************************************
// 
//
// $Header$
//
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu
//
//  Added a method that initializes an array to zero in addition to
//  resizing it.  -- Karim Filali (karim@cs.washington.edu)
//
// added various convenience operators and functions: +=, /=,
// dot_product, elementwise_product. -- Simon King


#ifndef SARRAY_H
#define SARRAY_H

#include <string.h>

#include "error.h"
#include "assert.h"

template <class T>
class sArray {

  int _size;

  // simple dirty ascending quicksort with dumb pivoting. Current type
  // must have assignment, comparison (operator <()), and be swapable
  // (i.e., copy constructor and/or writable)
  void internalSort(int begin,int end) {
    if (end > begin) {
      // last element is pivot. Could also choose a random median-like
      // element and swap it to the end.
      T& pivot = ptr[end];
      int l = begin;
      int r = end - 1;
      while (l < r) {
	if (ptr[l] < pivot) {
	  l++;
	} else {
	  T tmp(ptr[l]);
	  ptr[l] = ptr[r];
	  ptr[r] = tmp;
	  r--;
	}
      }
      if (ptr[l] < pivot) {
	r++;
      } else {
	T tmp(ptr[end]);
	ptr[end] = ptr[l];
	ptr[l] = tmp;
	l--;
	r++;
      }
      internalSort(begin, l);
      internalSort(r, end);
    }
  }

 public:

  T *ptr;

  sArray(int arg_size=0) {
    _size = arg_size;
    ptr = NULL;
    if (_size < 0)
      coredump("Error: sArray::sArray arg_size < 0");
    if (_size > 0)
      ptr = new T[_size];
  }

  ~sArray() {
    delete [] ptr;
  }

  // We don't create an operator= since we don't
  // use this array like a regular container class.
  sArray<T>* makeCopyOfSelf()  {
    // makes a copy of self with duplicated memory
    sArray<T>* cpy = new sArray<T>(_size);
    ::memcpy((void*)cpy->ptr,(void*)ptr,sizeof(T)*_size);
    return cpy;
  }

  void copyOtherIntoSelf(sArray<T>& other)  {
    delete [] ptr;
    _size  = other._size;
    ptr = new T[_size];
    ::memcpy((void*)ptr,(void*)other.ptr,sizeof(T)*_size);
  }


  void resize(int arg_size) {
    if (arg_size < 0)
      coredump("Error: Sarray:resize arg_size < 0");
    delete [] ptr;
    _size = arg_size;
    ptr = new T[_size];
  }

  void resizeAndZero(int arg_size) {
    if (arg_size < 0)
      error("Error: Sarray:resize arg_size < 0");
    delete [] ptr;
    _size = arg_size;
    ptr = new T[_size];
    for(int i=0; i < _size; ++i)
      *(ptr+i)=0;
  }


  void resizeIfDifferent(int arg_size) {
    // note that this does not call destructors if the sizes are the
    // same.
    if (arg_size != _size)
      resize(arg_size);
  }
  void growIfNeeded(const int arg_size) {
    if (arg_size > _size)
      resize(arg_size);
  }

  void growByNIfNeeded(const int n,const int arg_size) {
    assert ( n >= 1 );
    if (arg_size > _size)
      resize(n*arg_size);
  }

  void growByFIfNeeded(const float f,const int arg_size) {
    assert ( f >= 1.0 );
    if (arg_size > _size)
      resize((int)(f*arg_size+1.0));
  }

  // see comment within this routine before using it.
  void resizeAndCopy(const int arg_size) {
    if (arg_size < 0)
      coredump("Error: Sarray:resize arg_size < 0");
    T* tmp = new T[arg_size];
    const int nsize = (_size<arg_size?_size:arg_size);
    ::memcpy((void*)tmp,(void*)ptr,sizeof(T)*nsize);
    // note that this could be a problem for objects that have
    // destructors, as this will call the destructor for the object
    // that still has live pointers. Care should be used when using an
    // sArray for non pointer types.
    delete [] ptr;
    _size = arg_size;
    ptr = tmp;
  }

	

  void growIfNeededAndCopy(int arg_size) {
    if (arg_size > _size)
      resizeAndCopy(arg_size);
  }

  void growByNIfNeededAndCopy(const int n,int arg_size) {
    if (arg_size > _size)
      resizeAndCopy(n*arg_size);
  }

  void growByFIfNeededAndCopy(const float f,int arg_size) {
    assert ( f >= 1.0 );
    if (arg_size > _size)
      resizeAndCopy((int)(f*arg_size+1.0));
  }

  void swapPtrs(sArray<T>& sa) {
    if (_size != sa._size)
      coredump("Error: Sarray:swapPtrs, different sizes");
    T *tmp = ptr;
    ptr = sa.ptr;
    sa.ptr = tmp;
  }


  void swap(sArray<T>& sa) {
    T *tmp = ptr;
    ptr = sa.ptr;
    sa.ptr = tmp;
    int tmp_size = _size;
    _size = sa._size;
    sa._size = tmp_size;
  }

  // Append array y to the end of this one
  void concatenate(sArray<T>& y)
  {
    // Resize ourselves
    if (y._size <= 0)
      return;

    int oldsize = _size;
    resizeAndCopy(_size + y._size);
    ::memcpy((void*)(ptr+oldsize),(void*)y.ptr,sizeof(T)*y._size);

  }
    
  inline void clear() {
    delete [] ptr;
    ptr = NULL;
    _size = 0;
  }

  inline bool contains(T& x)
  {
    for (int i=0; i<_size; i++)
      if (ptr[i] == x)
        return true;
    return false;
  }
    
  inline int len() const { return _size; }
  inline unsigned size() const { return (unsigned)_size; }
  inline T& operator[](const int i) { 
    assert ( i >= 0 && i < _size );
    return ptr[i]; 
  }
  inline T operator[] (const int i) const { 
    assert ( i >= 0 && i < _size );
    return ptr[i]; 
  }
  
  void assignAllToValue(const T x) {
    for (int i=0; i<_size; i++)
      ptr[i] = x;
  }

  // simple naive implementation of ascending quicksort.  Current
  // type must have assignment, comparison (operator <()), and be
  // swapable (i.e., copy constructor or writable)
  void sort() { internalSort(0,_size-1); }
  // sort an subset range, inclusive [start,end]
  void sort(unsigned start,unsigned end) { 
    assert (start >= 0 && start < _size);
    assert (end >= 0 && end < _size);
    internalSort(start,end); 
  }


  inline sArray<T>& operator += (sArray<T> &s){
    assert(_size == s.len());
    for (int i=0; i<_size; i++)
      ptr[i] += s[i];
    return *this;
  }

  inline sArray<T>& operator -= (sArray<T> &s){
    assert(_size == s.len());
    for (int i=0; i<_size; i++)
      ptr[i] -= s[i];
    return *this;
  }

  inline sArray<T>& operator /= (const T n){
    for (int i=0; i<_size; i++)
      ptr[i] /= n;
    return *this;
  }


};


template<class T>
inline sArray<T> elementwise_product(const sArray<T> &s1, const sArray<T> &s2){
  assert(s1.len() == s2.len());
  sArray<T> rval;
  rval.resize(s1.len());
  for (int i=0; i<rval.len(); i++)
    rval[i] = s1[i] * s2[i];
  return rval;
}

template<class T>
inline sArray<T> operator + (const sArray<T> &s1, const sArray<T> &s2){
  assert(s1.len() == s2.len());
  sArray<T> rval;
  rval.resize(s1.len());
  for (int i=0; i<rval.len(); i++)
    rval[i] = s1[i] + s1[i];
  return rval;
}


template<class T>
inline sArray<T> operator * (const T a, const sArray<T> &s){
  sArray<T> rval;
  rval.resize(s.len());
  for (int i=0; i<rval.len(); i++)
    rval[i] = a * s[i];
  return rval;
}

template<class T>
inline sArray<T> operator * (const sArray<T> &s, const T a){ return a * s; }


inline sArray<double> operator * (const double a, const sArray<float> &s){
  sArray<double> rval;
  rval.resize(s.len());
  for (int i=0; i<rval.len(); i++)
    rval[i] = a * s[i];
  return rval;
}

#endif
