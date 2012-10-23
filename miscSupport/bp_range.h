//
// A general streams-like class with iterators for
// parsing and and using bi-polar range specifications of
// the form "-5:-1,0,25,9:16,19:3:45" or "1:3,5:8,12,15,22:35"
// 
// Jeff Bilmes <bilmes@ee.washington.edu>
//
// $Header$
//

#ifndef BP_RANGE_H
#define BP_RANGE_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>

#include "sArray.h"

class BP_Range {


  class sub_range {
    // This class specifies a subrange of
    // the Matlab form [lower:step:upper]
  public:
    // the lowest element of the range
    int lower;
    // the largest value if step==1
    int upper;
    // the step size.
    int step;
    // the actual largest value.
    int max() { 
      const int n = (upper-lower+step)/step; 
      return lower + step*(n-1);
    }
  };

  // the actuall range array
  sArray < sub_range > range;
  // The number of valid elements of the range array
  // (might be less than what is allocated) 
  int range_size;
  // A copy of the range spec string.
  char *range_str;

  // The smallest valid range value
  int lower_limit;
  // The largest valid range value (corresponding to ^0)
  int upper_limit;

  // the total number of entries that a
  // range specifies, i.e., "1,2,4:7" would
  // give a range_set_size of 6.
  int range_set_size;

  // the first element. i.e., (int)begin().
  int _min;
  // the last element. i.e., (int)end().
  int _max;

  // routine to parse a range string.
  void parse_range();

public:


  BP_Range(const char * range_str,  // Range specification.
	const int lower_limit_a, // Limits of a valid range, 
	const int upper_limit_a); // must be in [l:(u-1)]

  // copy constructor and copy.
  BP_Range(const BP_Range& other);
  BP_Range& operator=(const BP_Range& other);


  //////////////////////////////////////////////////////////////////
  // Create an *INVALID* version of this object for re-construction
  // within an array later using C++ placement new&.
  BP_Range() 
    : range_str(NULL),lower_limit(0),upper_limit(0) {}

  ~BP_Range();
  

  // Iterators on range objects.
  class iterator {
    friend class BP_Range;
  public:
    int cur_value;
    int cur_upper;
    int cur_lower;
    int array_pos;
    const BP_Range& myrange;
  public:

    iterator (const BP_Range& rng);
    iterator (const iterator&);

    iterator& operator ++(); // prefix
    iterator operator ++(int) { iterator tmp=*this; ++*this; return tmp; }
    iterator& operator --();  // prefix
    iterator operator --(int) { iterator tmp=*this; --*this; return tmp; }
    
    const int operator *() { return cur_value; }
    operator int() { return cur_value; }

    bool operator ==(const iterator& it){ return (cur_value == it.cur_value);}
    bool operator !=(const iterator& it){ return (cur_value != it.cur_value);}
    bool operator <(const iterator& it) { return (cur_value < it.cur_value);}
    bool operator <=(const iterator& it){ return (cur_value <= it.cur_value);}
    bool operator >(const iterator& it) { return (cur_value > it.cur_value);}
    bool operator >=(const iterator& it){ return (cur_value >= it.cur_value);}
    
    
    // Just a place holder for now, so that other programs start using
    // at_end instead of BP_Ragne::end() -- Karim
    bool at_end() { return (*this > myrange.end()); }

  };
  friend class iterator;

  // The first valid iterator.
  iterator begin () const;
  // The last valid iterator.
  iterator end () const;

  // Returns the number of elements in this range set.
  int length() const { return range_set_size; }
  // Returns the value of the minimum element.
  int min() const { return _min; }
  // Returns the value of the maximum element.
  int max() const { return _max; }

  // Returns true if value is contained in this range set.
  bool contains(const int value) const;

  static bool emptyRangeSpec(const char *str);
  static bool fullRangeSpec(const char *str);

  bool emptyRangeSpec () const { return emptyRangeSpec(range_str); }
  bool fullRangeSpec () const { return fullRangeSpec(range_str); }

  ////////////////////////////////////////////////////////
  // returns a pointer to the range string.
  char*const rangeStr  () const { return range_str; };

  // returns true if r has any overlap with this (i.e., if
  // the intersection is non-null)
  bool overlapP (const BP_Range& r) const;
  bool overlapP (const BP_Range* r) const { return overlapP(*r); }


  bool operator<  (const BP_Range& r)  const { return (max() < r.min()); }
  // operator <=, an "==" range is one that is non-comparable, i.e.,
  // one that has boundaries that overlaping.
  bool operator<=  (const BP_Range& r)  const 
  { return (max() < r.min()) || (min() <= r.max()); }

  bool operator>  (const BP_Range& r)  const { return (min() > r.max()); }
  // operator >=, an "==" range is one that is non-comparable, i.e.,
  // one that has boundaries that overlaping.
  bool operator>=  (const BP_Range& r)  const 
  { return (min() > r.max()) || (max() >= r.min()); }


  // essentially an "==" operator.
  bool boundariesOverlap(const  BP_Range& r) const {
    if (max() < r.min() || min() > r.max())
      return false;
    return true;
  }
  bool operator== (const BP_Range& r)  const 
  { return boundariesOverlap(r); }


  bool operator<  (const int r)  const { return (max() < r); }
  bool operator<= (const int r)  const { return (max() < r || contains(r)); }
  bool operator>  (const int r)  const { return (min() > r); }
  bool operator>= (const int r)  const { return (min() > r || contains(r)); }
  bool operator== (const int r)  const { return contains(r); }

};





#endif

