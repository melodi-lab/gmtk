//
//
// sArray: a Simple Array class that is basically a glorified type*.
// 
// This class is NOT meant for protection. It is mean to have type*
// (such as an int*, float*, double*, etc.)  and where the underlying
// pointer is easy to access for low level loops. But this class also
// provides for convenient managing of lengths, allocation,
// deallocation, and resizing.
//
// ******************* NOTE *****************************************
// There is no COPY CONSTRUCTOR, so all copies are bit-wise
// copies. This means that sArrays work *COMPLETELY* differently then
// say STL vector types, etc.. sArrays are meant to be used when you
// wish not to create lots of temporaries during
// construction/destruction in routine return and call.
// ******************* NOTE *****************************************
// 
//
// $Header$
//
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu
//

#ifndef SARRAY_H
#define SARRAY_H

#include <string.h>

#include "error.h"
#include "assert.h"

template <class T>
class sArray {

  int _size;

 public:

  T *ptr;

  sArray(int arg_size=0) {
    _size = arg_size;
    ptr = NULL;
    if (_size < 0)
      error("Error: sArray::sArray arg_size < 0");
    if (_size > 0)
      ptr = new T[_size];
  }

  ~sArray() {
    delete [] ptr;
  }
  void resize(int arg_size) {
    if (arg_size < 0)
      error("Error: Sarray:resize arg_size < 0");
    delete [] ptr;
    _size = arg_size;
    ptr = new T[_size];
  }
  void resizeIfDifferent(int arg_size) {
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

  void resizeAndCopy(const int arg_size) {
    if (arg_size < 0)
      error("Error: Sarray:resize arg_size < 0");
    T* tmp = new T[arg_size];
    const int nsize = (_size<arg_size?_size:arg_size);
    ::memcpy((void*)tmp,(void*)ptr,sizeof(T)*nsize);
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
      error("Error: Sarray:swapPtrs, different sizes");
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


};


#endif
