// -*- c++ -*-
// A general streams-like class with iterators for
// parsing and and using range specifications of
// the form "2,5,9-16,19:3:45" or "1-3,5-8,12,15,22-35"
// 
// Jeff Bilmes <bilmes@icsi.berkeley.edu> Sep 1997
//
// $Header$
//

#ifndef RANGE_H
#define RANGE_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>

#ifndef HAVE_BOOL
#ifndef BOOL_DEFINED
enum bool { false = 0, true = 1 };
#define BOOL_DEFINED
#endif
#endif

class Range {


  class sub_range {
    // This class specifies a subrange of
    // the Matlab form [lower:step:upper]
  public:
    // the lowest element of the range
    size_t lower;
    // the largest value if step==1
    size_t upper;
    // the step size.
    size_t step;
    // the actual largest value.
    size_t max() { 
      const size_t n = (upper-lower+step)/step; 
      return lower + step*(n-1);
    }
  };

  // the actual range array
  sub_range* range;
  // Ensures we can add one to the range array
  void make_range_safe_to_add_one();
  // The number of valid elements of the range array
  size_t range_size;
  // The actual range array size
  size_t range_array_length;
  // A copy of the range spec string.
  char *range_str;

  // The smallest valid range value
  const size_t lower_limit;
  // The largest valid range value (corresponding to ^0)
  const size_t upper_limit;

  // the total number of entries that a
  // range specifies, i.e., "1,2,4-7" would
  // give a range_set_size of 6.
  size_t range_set_size;

  // the first element in the un-negated range.
  size_t _min_un;
  // the last element in the un-negated range.
  size_t _max_un;

  // the actual first element
  size_t _min;
  // the actual last element
  size_t _max;


  // "all-but" mode, i.e., negate the meaning.
  bool all_but_mode;

  // routine to parse a range string.
  void parse_range();

  // Internal 'contains' routine
  bool _contains(const size_t value) const;

public:


  Range(const char * range_str,  // Range specification.
	const size_t lower_limit_a, // Limits of a valid range, 
	const size_t upper_limit_a); // must be in [l:(u-1)]
  ~Range();
  

  // Iterators on range objects.
  class iterator {
    friend class Range;
  public:
    size_t cur_value;
    size_t cur_upper;
    size_t cur_lower;
    size_t array_pos;
    const Range& myrange;
  public:

    iterator (const Range& rng);
    iterator (const iterator&);


    iterator& operator ++(); // prefix
    iterator operator ++(int) { iterator tmp=*this; ++*this; return tmp; }
    iterator& operator --();  // prefix
    iterator operator --(int) { iterator tmp=*this; --*this; return tmp; }
    
    const size_t operator *() { return cur_value; }
    operator size_t() { return cur_value; }

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
  size_t length() const { return range_set_size; }
  // Returns the value of the minimum element.
  size_t min() const { return _min; }
  // Returns the value of the maximum element.
  size_t max() const { return _max; }

  // Returns true if value is contained in this range set.
  bool contains(const size_t value) const;

  static bool emptyRangeSpec(const char *str);
  static bool fullRangeSpec(const char *str);

  bool emptyRangeSpec() { return emptyRangeSpec(range_str); }
  bool fullRangeSpec() { return fullRangeSpec(range_str); }

};





#endif

