//
//
// cArray: a Simple Array class that is basically a glorified type*.
//
// This is just like sArray except it uses C malloc/free rather
// than C++'s new/delete. Moreover, it explicitly calls a default constructor
// upon creation of any new objects (like new does) and explicitely
// calls the destructor on deletion of the entire object (like delete) does.
// It does not do any construction/destruction on a resize, the array
// rather is bit-wise copied avoiding creation of temporary objects, etc.
//
// 
// This class is NOT meant for protection. It is mean to have type*
// (such as an int*, float*, double*, etc.)  and where the underlying
// pointer is easy to access for low level loops. But this class also
// provides for convenient managing of lengths, allocation,
// deallocation, and resizing.
//
// ******************* NOTE *****************************************
// There is no COPY CONSTRUCTOR, so all copies are bit-wise
// copies. This means that cArrays work *COMPLETELY* differently then
// say STL vector types, etc.. cArrays are meant to be used when you
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

#ifndef CARRAY_H
#define CARRAY_H

#include <stdlib.h>
#include <string.h>
#include <new>

#include "error.h"
#include "assert.h"

template <class T>
class cArray {

  int _size;

 public:

  T *ptr;

  cArray(int arg_size=0) {
    _size = arg_size;
    ptr = NULL;
    if (_size < 0)
      error("Error: cArray::cArray arg_size < 0");
    if (_size == 0)
      return;
    if (_size > 0)
      ptr = (T*) malloc(sizeof(T)*_size);
    // explicitly call default constructors
    T* ptr_p = ptr;
    T* ptr_endp = ptr + _size;
    while (ptr_p != ptr_endp) {
      new (ptr_p) T();
      ptr_p++;
    }
  }

  ~cArray() { clear(); }
    
  inline void clear() {
    if (ptr == NULL) return;
    T* ptr_p = ptr;
    T* ptr_endp = ptr + _size;
    // explicitly call destructor when cArray object
    // is cleared.
    while (ptr_p != ptr_endp) {
      ptr_p->~T();
      ptr_p++;
    }
    free ((void*)ptr);
    ptr = NULL;
    _size = 0;
  }

  // resize without worrying about anything else.
  void resize(int arg_size) {
    if (arg_size < 0)
      error("Error: cArray:resize arg_size < 0");

    // delete old, calling destructors along the way.
    if (ptr != NULL) {
      T* ptr_p = ptr;
      T* ptr_endp = ptr + _size;
      while (ptr_p != ptr_endp) {
	ptr_p->~T();
	ptr_p++;
      }
      free ((void*)ptr);
    }

    _size = arg_size;
    if (_size == 0) {
      ptr = NULL;
      return;
    }

    ptr = (T*) malloc(sizeof(T)*_size);
    // explicitly call default constructors for new
    T* ptr_p = ptr;
    T* ptr_endp = ptr + _size;
    while (ptr_p != ptr_endp) {
      new (ptr_p) T();
      ptr_p++;
    }

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

  // bit-copy old entries, call default constructors for new entries.
  void resizeAndCopy(const int arg_size) {
    if (arg_size < 0)
      error("Error: cArray:resize arg_size < 0");
    if (arg_size > _size) {
      // allocate new larger size
      T* tmp = (T*) malloc(sizeof(T)*arg_size);
      // bit-wise copy over old portion
      ::memcpy((void*)tmp,(void*)ptr,sizeof(T)*_size);
      // and call default constructors only for new portion
      T* ptr_p = tmp + _size;
      T* ptr_endp = tmp + arg_size;
      while (ptr_p != ptr_endp) {
	new (ptr_p) T();
	ptr_p++;
      }
      // delete old, without calling element-wise destructors
      if (ptr != NULL)
	free ((void*)ptr);
      // finally store info about new
      _size = arg_size;
      ptr = tmp;
    } else if (arg_size < _size) {

      T* tmp;
      if (arg_size > 0) {
	// allocate new smaller size
	tmp = (T*) malloc(sizeof(T)*arg_size);
	// bit-wise copy over portion that will remain
	::memcpy((void*)tmp,(void*)ptr,sizeof(T)*arg_size);
      } else {
	tmp = NULL;
      }

      // and call default destructors for old portion going away
      T* ptr_p = ptr + arg_size;
      T* ptr_endp = ptr + _size;
      while (ptr_p != ptr_endp) {
	ptr_p->~T();
	ptr_p++;
      }
      // delete old, without calling element-wise destructors.
      free ((void*)ptr);
      // finally keep track of new.
      _size = arg_size;
      ptr = tmp;
    } 

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

  void swapPtrs(cArray<T>& sa) {
    if (_size != sa._size)
      error("Error: cArray:swapPtrs, different sizes");
    T *tmp = ptr;
    ptr = sa.ptr;
    sa.ptr = tmp;
  }

  void swap(cArray<T>& sa) {
    T *tmp = ptr;
    ptr = sa.ptr;
    sa.ptr = tmp;
    int tmp_size = _size;
    _size = sa._size;
    sa._size = tmp_size;
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
