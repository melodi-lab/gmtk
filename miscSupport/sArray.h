//
//
// sArray: a Simple Array class that is basically a glorified 
// type*.
// 
// This class is NOT meant for protection. It is
// mean to have type* (such as an int*, float*, double*, etc.)
// and where the underlying pointer is easy to access for
// low level loops. But this class also provides
// for convenient managing of lengths, allocation, deallocation,
// and resizing.
//
// ****** NOTE ***** There is no COPY CONSTRUCTOR.

// $Header$
//
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu
//

#ifndef SARRAY_H
#define SARRAY_H

#include <string.h>

#include "error.h"
#include "assert.h"

template <class T>
class sArray {

  int size;

 public:

  T *ptr;

  sArray(int _size=0) {
    size = _size;
    ptr = NULL;
    if (size < 0)
      error("Error: sArray::sArray _size < 0");
    if (size > 0)
      ptr = new T[size];
  }

  ~sArray() {
    delete [] ptr;
  }
  void resize(int _size) {
    if (_size < 0)
      error("Error: Sarray:resize _size < 0");
    delete [] ptr;
    size = _size;
    ptr = new T[size];
  }
  void growIfNeeded(const int _size) {
    if (_size > size)
      resize(_size);
  }

  void growByNIfNeeded(const int n,const int _size) {
    if (_size > size)
      resize(n*_size);
  }

  void resizeAndCopy(const int _size) {
    if (_size < 0)
      error("Error: Sarray:resize _size < 0");
    T* tmp = new T[_size];
    const int nsize = (size<_size?size:_size);
    ::memcpy((void*)tmp,(void*)ptr,sizeof(T)*nsize);
    delete [] ptr;
    size = _size;
    ptr = tmp;
  }

  void growIfNeededAndCopy(int _size) {
    if (_size > size)
      resizeAndCopy(_size);
  }

  void growByNIfNeededAndCopy(const int n,int _size) {
    if (_size > size)
      resizeAndCopy(n*_size);
  }


  void swapPtrs(sArray<T>& sa) {
    if (size != sa.size)
      error("Error: Sarray:swapPtrs, different sizes");
    T *tmp = ptr;
    ptr = sa.ptr;
    sa.ptr = tmp;
  }

  // Append array y to the end of this one
  void concatenate(sArray<T>& y)
  {
    // Resize ourselves
    if (y.size <= 0)
      return;

    int oldsize = size;
    resizeAndCopy(size + y.size);
    ::memcpy((void*)(ptr+oldsize),(void*)y.ptr,sizeof(T)*y.size);

  }
    
  inline void clear() {
    delete [] ptr;
    ptr = NULL;
    size = 0;
  }

  inline bool contains(T& x)
  {
    for (int i=0; i<size; i++)
      if (ptr[i] == x)
        return true;
    return false;
  }
    
  inline int len() const { return size; }
  inline T& operator[](const int i) { 
    assert ( i >= 0 && i < size );
    return ptr[i]; 
  }
  inline T operator[] (const int i) const { 
    assert ( i >= 0 && i < size );
    return ptr[i]; 
  }
  
  void assignAllToValue(const T x) {
    for (int i=0; i<size; i++)
      ptr[i] = x;
  }


};



#endif
