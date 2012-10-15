//
//
// mArray: a Simple Array class that is basically a glorified type*.
//
// This version uses mmap to allocate memory.
// 
// This class is NOT meant for protection or abstraction. It is mean
// to have type* (such as an int*, float*, double*, etc.)  and where
// the underlying pointer is easy to access for low level loops (and
// where a few simple ones a provided for you). But this class also
// provides for convenient managing of lengths, allocation,
// deallocation, and resizing. It also includes a simple sort routine
// that will work as long as the object T has an operator<(), copy
// constructor and operator=().
//
// ******************* NOTE *****************************************
// There is no COPY CONSTRUCTOR for this object, so all copies are
// bit-wise copies. This means that mArrays work *COMPLETELY*
// differently then say STL vector types, etc.. mArrays are meant to
// be used when you wish not to create lots of temporaries during
// construction/destruction in routine return and call. If you include
// mArrays in a class, you will need to explicitly manage any
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
//
//  mmap() version -- Richard Rogers


#ifndef MARRAY_H
#define MARRAY_H

#include <new>
#include <string.h>
#include <sys/mman.h>

#include "error.h"
#include "assert.h"

template <class T>
class mArray_nd {

 protected:
  int _size;
  bool _usemmap;


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

// see https://j.ee.washington.edu/trac/gmtk/ticket/371
#ifndef   MAP_ANONYMOUS
#  define MAP_ANONYMOUS MAP_ANON
#endif

  mArray_nd(int arg_size=0, bool usemmap=true) {
    _size = arg_size;
    _usemmap = usemmap;
    ptr = NULL;
    if (_size < 0)
      coredump("Error: mArray::mArray arg_size < 0");
    if (_size > 0) {
#ifdef HAVE_MMAP
      if (_usemmap) {
	void *block = mmap(NULL, _size * sizeof(T), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
	assert(block);
	ptr = new(block) T[_size];
      } else {
	ptr = new T[_size];
      }
#else
      ptr = new T[_size];
#endif
      assert(ptr);
    }
  }

  ~mArray_nd() {
    // deallocate();
  }

  void allocate(int arg_size=0) {
#ifdef HAVE_MMAP
    if (_usemmap) {
      void *block = mmap(NULL, arg_size * sizeof(T), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
      assert(block);
      ptr = new(block) T[arg_size];
    } else {
      ptr = new T[arg_size];
    }
#else
    ptr = new T[arg_size];
#endif
    assert(ptr);
    _size = arg_size;
  }

  void deallocate() {
#ifdef HAVE_MMAP
    if (_usemmap) {
      while (_size) 
	ptr[--_size].~T();
      if (ptr) munmap(ptr, _size * sizeof(T));
    } else {
      delete [] ptr;
    }
#else
    delete [] ptr;
#endif
    _size = 0;
    ptr = NULL;
  }

  // We don't create an operator= since we don't
  // use this array like a regular container class.
  mArray_nd<T>* makeCopyOfSelf()  {
    // makes a copy of self with duplicated memory
    mArray_nd<T>* cpy = new mArray_nd<T>(_size, _usemmap);
    ::memcpy((void*)cpy->ptr,(void*)ptr,sizeof(T)*_size);
    return cpy;
  }

  void copyOtherIntoSelf(mArray_nd<T>& other)  {
    deallocate();
    allocate(other._size);
    ::memcpy((void*)ptr,(void*)other.ptr,sizeof(T)*_size);
  }


  void resize(int arg_size) {
    if (arg_size < 0)
      coredump("Error: mArray:resize arg_size < 0");
    deallocate();
    allocate(arg_size);
  }

  void resizeAndZero(int arg_size) {
    if (arg_size < 0)
      error("Error: mArray:resize arg_size < 0");
    deallocate();
    allocate(arg_size);
    if (!_usemmap) {
      // MAP_ANONYMOUS zeros the memory
      for(int i=0; i < _size; ++i)
	*(ptr+i)=0;
    }
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
      coredump("Error: mArray:resize arg_size < 0");
    T* tmp;
#ifdef HAVE_MMAP
    if (_usemmap) {
      void *block = mmap(NULL, arg_size * sizeof(T), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
      assert(block);
      tmp = new(block) T[arg_size];
    } else {
      tmp = new T[arg_size];
    }
#else
    tmp = new T[arg_size];
#endif
    assert(tmp);
    const int nsize = (_size<arg_size?_size:arg_size);
    ::memcpy((void*)tmp,(void*)ptr,sizeof(T)*nsize);
    // note that this could be a problem for objects that have
    // destructors, as this will call the destructor for the object
    // that still has live pointers. Care should be used when using an
    // mArray_nd for non pointer types.
    deallocate();
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

  void swapPtrs(mArray_nd<T>& sa) {
    if (_size != sa._size)
      coredump("Error: mArray:swapPtrs, different sizes");
    T *tmp = ptr;
    ptr = sa.ptr;
    sa.ptr = tmp;

    bool tmpmmap = _usemmap;
    _usemmap = sa._usemmap;
    sa._usemmap = tmpmmap;
  }


  void swap(mArray_nd<T>& sa) {
    T *tmp = ptr;
    ptr = sa.ptr;
    sa.ptr = tmp;
    int tmp_size = _size;
    _size = sa._size;
    sa._size = tmp_size;
    bool tmpmmap = _usemmap;
    _usemmap = sa._usemmap;
    sa._usemmap = tmpmmap;
  }

  // Append array y to the end of this one
  void concatenate(mArray_nd<T>& y)
  {
    // Resize ourselves
    if (y._size <= 0)
      return;

    int oldsize = _size;
    resizeAndCopy(_size + y._size);
    ::memcpy((void*)(ptr+oldsize),(void*)y.ptr,sizeof(T)*y._size);

  }
    
  inline void clear() {
    deallocate();
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


  inline mArray_nd<T>& operator += (mArray_nd<T> &s){
    assert(_size == s.len());
    for (int i=0; i<_size; i++)
      ptr[i] += s[i];
    return *this;
  }

  inline mArray_nd<T>& operator -= (mArray_nd<T> &s){
    assert(_size == s.len());
    for (int i=0; i<_size; i++)
      ptr[i] -= s[i];
    return *this;
  }

  inline mArray_nd<T>& operator /= (const T n){
    for (int i=0; i<_size; i++)
      ptr[i] /= n;
    return *this;
  }


};


// mArray with a destructor.
template <class T>
class mArray : public mArray_nd<T> {
 public:
 mArray(int arg_size=0)  : mArray_nd<T>(arg_size) {}
  ~mArray() {
    if (mArray_nd<T>::ptr) mArray_nd<T>::deallocate();
  }
}; 



template<class T>
inline mArray<T> elementwise_product(const mArray<T> &s1, const mArray<T> &s2){
  assert(s1.len() == s2.len());
  mArray<T> rval;
  rval.resize(s1.len());
  for (int i=0; i<rval.len(); i++)
    rval[i] = s1[i] * s2[i];
  return rval;
}

template<class T>
inline mArray<T> operator + (const mArray<T> &s1, const mArray<T> &s2){
  assert(s1.len() == s2.len());
  mArray<T> rval;
  rval.resize(s1.len());
  for (int i=0; i<rval.len(); i++)
    rval[i] = s1[i] + s1[i];
  return rval;
}


template<class T>
inline mArray<T> operator * (const T a, const mArray<T> &s){
  mArray<T> rval;
  rval.resize(s.len());
  for (int i=0; i<rval.len(); i++)
    rval[i] = a * s[i];
  return rval;
}

template<class T>
inline mArray<T> operator * (const mArray<T> &s, const T a){ return a * s; }


inline mArray<double> operator * (const double a, const mArray<float> &s){
  mArray<double> rval;
  rval.resize(s.len());
  for (int i=0; i<rval.len(); i++)
    rval[i] = a * s[i];
  return rval;
}



#endif
