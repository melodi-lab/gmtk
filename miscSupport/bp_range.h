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
  sub_range* range;
  // Ensures we can add one to the range array
  void make_range_safe_to_add_one();
  // The number of valid elements of the range array
  int range_size;
  // The actuall range array size
  int range_array_length;
  // A copy of the range spec string.
  char *range_str;

  // The smallest valid range value
  const int lower_limit;
  // The largest valid range value (corresponding to ^0)
  const int upper_limit;

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
  };
  friend class iterator;

  // The first valid iterator.
  iterator begin();
  // The last valid iterator.
  iterator end();

  // Returns the number of elements in this range set.
  int length() { return range_set_size; }
  // Returns the value of the minimum element.
  int min() { return _min; }
  // Returns the value of the maximum element.
  int max() { return _max; }

  // Returns true if value is contained in this range set.
  bool contains(const int value);

  static bool emptyRangeSpec(const char *str);
  static bool fullRangeSpec(const char *str);

  bool emptyRangeSpec() { return emptyRangeSpec(range_str); }
  bool fullRangeSpec() { return fullRangeSpec(range_str); }

  ////////////////////////////////////////////////////////
  // returns a pointer to the range string.
  char*const rangeStr() { return range_str; };

  // returns true if r has any overlap with this (i.e., if
  // the intersection is non-null)
  bool overlapP(BP_Range& r);
  bool overlapP(BP_Range* r) { return overlapP(*r); }


};





#endif

