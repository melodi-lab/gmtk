// 
// A general streams-like class with iterators for
// parsing and and using range specifications of
// the form "2,5,9-16,19:3:45" or "1-3,5-8,12,15,22-35"
// 
// Jeff Bilmes <bilmes@icsi.berkeley.edu>  Sep 1997
//
// $Header$
//


#include <limits.h>
#include <string.h>

#include "rrange.h"
#include "error.h"

#define ALL_BUT_CHAR '/'
// must be the same as ALL_BUT_CHAR but a string.
#define ALL_BUT_STR "/"

#define FROM_END_CHAR '^'


const size_t MAX_SIZE_T = UINT_MAX;

Range::Range(const char *range_str_a,
	     const size_t lower_limit_a,
	     const size_t upper_limit_a)
  : lower_limit(lower_limit_a),upper_limit(upper_limit_a),
    all_but_mode(range_str_a && (*range_str_a == ALL_BUT_CHAR))
  //
  // Note: a range spec must be within the limits 
  // [lower_limit:(upper_limit-1)]
  // The reason for this is that lower_limit == upper_limit is
  // the special case of only a "null" or "full" range being a valid 
  // range spec, in either case the spec designates a null range.
{

  if (range_str_a == NULL) 
    range_str_a = "all";
  range_str = new char[::strlen(range_str_a)+1];
  ::strcpy(range_str,range_str_a);


  if (upper_limit < lower_limit)
    error("Error upper_limit=%u < lower_limit=%u.",upper_limit,lower_limit);

  // if (upper_limit == lower_limit) {
  // empty range spec. 
  // }


  range_array_length = 20; 
  range_size = 0;
  range = new sub_range[range_array_length];

  // parse the current range spec.
  parse_range();

  // The following must be called after parse_range.
  _min_un = range[0].lower;
  _max_un = range[range_size-1].max();
  if (all_but_mode) {
    size_t val;

    // Do dumb linear search. This could take a long
    // time if the user specified a large limit range 
    // and large range.
    val = lower_limit;
    while (_contains(val) && val < upper_limit)
      val++;
    _min = val;

    val = upper_limit-1;
    while (_contains(val) && val >= lower_limit)
      val--;
    _max = val;

  } else {
    _min = _min_un;
    _max = _max_un;
  }


#if 0
  for (size_t i=0;i<range_size;i++) {
    printf ("[%u:%u:%u]\n",range[i].lower,range[i].step,
  	    range[i].upper);
  }
  printf("_max = %u, _min = %u\n",_max,_min);
#endif


}

Range::~Range()
{
  delete [] range;
  delete [] range_str;
}


bool
Range::emptyRangeSpec(const char *str)
{
  if (str == NULL)
    return false;
  return (!strcmp(str,"nil") ||
	  !strcmp(str,"none") ||
	  !strcmp(str,ALL_BUT_STR "all") ||
	  !strcmp(str,ALL_BUT_STR "full") ||
	  !strcmp(str,ALL_BUT_STR "-") ||
	  !strcmp(str,ALL_BUT_STR ":") ||
	  !strcmp(str,ALL_BUT_STR));
}


bool
Range::fullRangeSpec(const char *str)
{
  return (str == NULL 
	|| !strcmp(str,"all")
	|| !strcmp(str,"full")
	|| !strcmp(str,"-")
	|| !strcmp(str,":")
	|| !strcmp(str,ALL_BUT_STR "nil")
	|| !strcmp(str,ALL_BUT_STR "none"));
}




/*
 * Parse a range set specification.
 *
 * A (probably ambiguous) grammer is:
 *  range -> token | list
 *  token -> 'nil' | 'none' | 'all' | 'full' | '-' | ':'
 *  list -> [ brange range_list erange ]
 *  range_list -> range | [ rangelist ',' range ] | nil
 *  range -> num | [num1 '-' num2] | [num1 ':' num3 ':' num2 ]
 *  brange -> [ '-' num ] | [ ':' num ] | [ ':' num3 ':' num2 ] | nil
 *  erange -> [ num '-' ] | [ num ':' ] | [ num1 ':' num3 ':' ] | nil
 *  num1 -> numa   # num1 corresponds to the lower limit of a subrange
 *  num2 -> numa   # num2 corresponds to the upper limit of a subrange
 *  num3 -> num   # num3 corresponds to the step size of a subrange
 *  num -> [0:UINT_MAX]
 *  numa -> [0:UINT_MAX] | ^[0:UINT_MAX]
 *  nil -> ""
 *
 *  Also
 *   num1 < num2
 *   'nil' - select no featuers
 *   'none' - select no featuers
 *   '-' - select all featuers in [0:n_feats-1]
 *   'all' - select all featuers in [0:n_feats-1]
 *   'full' - select all featuers in [0:n_feats-1]
 * 
 * A range may also be preceeded by a '/' character to make it
 * do the negation of the range (within limits) (i.e., selects
 * all but what is specified by the range).
 *   
 * Examples:
 *   0,4-7,14,30-35,40:5:65
 * corresponds to:
 *   0,4,5,6,7,14,30,31,32,32,34,35,40,45,50,55,60,65
 */
void Range::parse_range()
{
  size_t i;
  const char *p = range_str;

  range_set_size = 0;
  if (p == NULL)
    error("Error: No range spec string argument given.");

  if (!strcmp(p,"nil") || 
      !strcmp(p,"none") ||
      !strcmp(p,ALL_BUT_STR "all") ||
      !strcmp(p,ALL_BUT_STR "full") ||
      !strcmp(p,ALL_BUT_STR "-") ||
      !strcmp(p,ALL_BUT_STR ":") ||
      !strcmp(p,ALL_BUT_STR)) {
    // empty range spec, select nothing.
    range[range_size].lower = 1;
    range[range_size].upper = 0;
    range[range_size].step = 1;
    range_size++;
    range_set_size = 0;
    if (all_but_mode)
      all_but_mode = 0;
    return;
  } else if (!strcmp(p,"all") || 
	     !strcmp(p,"full") ||
	     !strcmp(p,"-") || 
	     !strcmp(p,":") ||
	     !strcmp(p,ALL_BUT_STR "nil") ||
	     !strcmp(p,ALL_BUT_STR "none")) {
    // complete subset, select all features.
    if (upper_limit == lower_limit) {
      // empty range spec, select nothing.
      range[range_size].lower = 1;
      range[range_size].upper = 0;
      range[range_size].step = 1;
      range_size++;
      range_set_size = 0;
    } else {
      range[range_size].lower = lower_limit;
      range[range_size].upper = upper_limit-1;
      range[range_size].step = 1;
      range_set_size = upper_limit - lower_limit;
      range_size++;
    }
    if (all_but_mode)
      all_but_mode = 0;
    return;
  } else if (upper_limit == lower_limit) {
    // only empty range is possible but user gave non-empty range.
    error("Error: non-empty range specification given for a nul range, (%s).",range_str);
  }

  const char *next;
  size_t l,u;
  size_t s;

  if (p && (*p == ALL_BUT_CHAR))
    p++;

  while (*p) {
    // set to true if a ~n number is specified.
    bool from_top;

    if (*p == '-' || *p == ':') {
      if (p == range_str) {
	// so we are at the beginning of the range string
	// and we've got a form ":..." or "-n".
	l = lower_limit;
      } else {
	error("Error: invalid range spec string (%s) starting at (%s)",
	      range_str,p);
	l=0; // remove compiler warning.
      }
      next = p;
    } else {
      from_top = false;
      if (*p == FROM_END_CHAR) {
	p++;
	from_top = true;
      }
      l = (size_t) strtol(p,(char**) &next, 0);
      if (next == p) {
	error("Error: invalid range spec string (%s) starting at (%s)",
	      range_str,p);
      }
      if (from_top) {
	l = upper_limit - l - 1;
      }
    }

    if (*next == '-') {
      // must be either "l-u" or "l-" form.
      p = next+1;
      if (!*p) {
	// must be a "l-"
	u = upper_limit - 1;
	next++;
      } else {
	// this is a "l-u"
	from_top = false;
	if (*p == FROM_END_CHAR) {
	  p++;
	  from_top = true;
	}
	u = (size_t) strtol(p,(char**) &next, 0);
	if (next == p) {
	  error("Error: invalid range spec string (%s) starting at (%s)",
		range_str,p);
	}
	if (from_top) {
	  u = upper_limit - u - 1;
	}
      }
      s = 1;
    } else if (*next == ':') {
      // must be either "l:u", "l:s:u", "l:s:", or "l:" form.
      p = next+1;

      if (!*p) {
	// must be a "l:"
	u = upper_limit - 1;
	s = 1;
	next++;
      } else {
	// must be either "l:u", "l:s:u", or "l:s:" form.
	// assume for now "l:s:u or "l:s:"
	from_top = false;
	if (*p == FROM_END_CHAR) {
	  p++;
	  from_top = true;
	}
	s = (size_t) strtol(p,(char**) &next, 0);
	if (next == p) {
	  error("Error: invalid range spec string (%s) starting at (%s)",
		range_str,p);
	}
	if (*next == ':') {
	  // then we indeed have "l:s:u" or "l:s:" form:
	  if (from_top) {
	    // the s can't have a FROM_END_CHAR before it.
	    error("Error: invalid range spec string (%s) starting at (%s)",
		  range_str,p);
	  }
	  p = next+1;
	  if (*p == ',' || !*p) {
	    // we have "l:s:" form
	    u = upper_limit - 1;
	    next = p;
	  } else {
	    // we should have "l:s:u"
	    from_top = false;
	    if (*p == FROM_END_CHAR) {
	      p++;
	      from_top = true;
	    }
	    u = (size_t) strtol(p,(char**)&next,0);
	    if (next == p) {
	      error("Error: invalid range spec string (%s) starting at (%s)",
		    range_str,p);
	    }
	    if (from_top) {
	      u = upper_limit - u - 1;
	    }
	  }
	} else {
	  // we have a "l:u" form, s holds upper, step size 1
	  if (from_top) {
	    u = upper_limit - s - 1;
	  } else 
	    u = s;
	  s = 1;
	}
      }
    } else {
      u = l;
      s = 1;
    }

    make_range_safe_to_add_one();
    range[range_size].lower = l;
    range[range_size].upper = u;    
    range[range_size].step = s;
    range_set_size += 
      ((range[range_size].upper - 
	range[range_size].lower + 
	range[range_size].step)/range[range_size].step);
    range_size++;
    
    p = next;
    if (*p == ',')
      p++;
  }

  // now check that the range specification is:
  //   1) sorted
  //   2) has no overlap
  //   3) within range upper and lower limits.


  i=0;
  if (range[i].lower > range[i].upper)
    error ("Error: lower %u > uupper %u in range string",range[i].lower,range[i].upper);
  if (range[i].lower < lower_limit)
    error ("Error: lower %u < limit %u in range string",range[i].lower,lower_limit);
  if (range[i].upper >= upper_limit)
    error ("Error: upper %u >= limit %u in range string",range[i].upper,upper_limit-1);

  for (i=1;i<range_size;i++) {
    if (range[i].lower > range[i].upper)
      error ("Error: lower %u > upper %u in range string",range[i].lower,range[i].upper);
    if (range[i].lower < lower_limit)
      error ("Error: lower %u < limit %u in range string",range[i].lower,lower_limit);
    if (range[i].upper >= upper_limit)
      error ("Error: upper %u >= limit %u in range string",range[i].upper,upper_limit-1);
    if (range[i-1].max() >= range[i].lower)
      error ("Error: previous upper %u >= current lower %u in range string",
	     range[i-1].upper,range[i].lower);
  }

  if (all_but_mode) {
    range_set_size = upper_limit - lower_limit - range_set_size;
  }

}




void
Range::make_range_safe_to_add_one()
{
  if (range_size+1 > range_array_length) {
    range_array_length *= 2;
    sub_range *tmp = new sub_range[range_array_length];
    ::memcpy((void*)tmp,(void*)range,
	     sizeof(*range)*range_size);
    delete [] range;
    range = tmp;
  }
}

Range::iterator Range::begin()
{
  iterator rc(*this);
  if (all_but_mode) {
    rc.cur_value = lower_limit;
    while (_contains(rc.cur_value) && rc.cur_value < upper_limit)
      rc.cur_value++;
  }
  return rc;
}

Range::iterator Range::end()
{
  iterator rc(*this);
  rc.array_pos = range_size-1;

  rc.cur_lower = range[rc.array_pos].lower;
  rc.cur_upper = range[rc.array_pos].upper;
  if (all_but_mode) {
    rc.cur_value = upper_limit-1;
    while (_contains(rc.cur_value) && rc.cur_value >= lower_limit)
      rc.cur_value--;    
  } else
    rc.cur_value = range[rc.array_pos].max();

  return rc;
}




//********************************************************
//********************************************************
// Range::iterator definitions
//********************************************************
//********************************************************

Range::iterator::iterator(const Range& rng)
  : myrange(rng)
{
  array_pos = 0;
  cur_value = myrange._min;
  cur_lower = myrange.range[array_pos].lower;
  cur_upper = myrange.range[array_pos].upper;
}

Range::iterator::iterator(const Range::iterator& it)
  : myrange(it.myrange)
{
  array_pos = it.array_pos;
  cur_value = it.cur_value;
  cur_upper = it.cur_upper;
  cur_lower = it.cur_lower;
}

Range::iterator&
Range::iterator::operator++()
{
  if (myrange.all_but_mode) {
    cur_value++;
    while (myrange._contains(cur_value) && cur_value < myrange.upper_limit)
      cur_value++;
  } else {
    size_t s = myrange.range[array_pos].step;
    if (cur_value+s <= cur_upper)
      cur_value += s;
    else {
      if (array_pos+1 < myrange.range_size) {
	array_pos ++;
	cur_lower = myrange.range[array_pos].lower;
	cur_upper = myrange.range[array_pos].upper;      
	cur_value = myrange.range[array_pos].lower;
      } else {
	// attempting to iterate off the bounds,
	cur_value = MAX_SIZE_T;
      }
    }
  }
  return *this;
}


Range::iterator&
Range::iterator::operator--()
{
  if (myrange.all_but_mode) {
    cur_value--;
    while (myrange._contains(cur_value) && cur_value >= myrange.lower_limit)
      cur_value--;
  } else {
    size_t s = myrange.range[array_pos].step;
    if (cur_value == MAX_SIZE_T)
      cur_value = myrange._max;
    else if (cur_value > cur_lower+s)
      cur_value -= s;
    else {
      if (array_pos > 0) {
	array_pos --;
	cur_lower = myrange.range[array_pos].lower;
	cur_upper = myrange.range[array_pos].upper;      
	const size_t s = myrange.range[array_pos].step;
	const size_t n = (cur_upper-cur_lower+s)/s;
	cur_value = cur_lower + s*(n-1);
      } else {
	// attempting to iterate off the bounds,
	// do nothing.
	// WARNING: This could lead to unexpected behavior, if
	// doing a loop such as: for (it=max;it>=min;it--).
      }
    }
  }
  return *this;
}


//
// Internal _contains function that tells us if value
// is contained in the none-negated range specification set
// (so if this is an all_but range, it ignores the "all but".
bool Range::_contains(const size_t value) const
{

  if (value < _min_un || value > _max_un || range_set_size == 0)
    return false;

  // Do a binary search that
  // finds the sub-range i such that
  // range[i].lower <= value < range[i+1].lower
  // or if (value == range[range_size-1].lower)
  // set i = range_size-1.
  size_t l,r,i;
  l = 0; r = range_size -1;
  while (l<r) {
    size_t m = (l+r+1)/2;
    if (value == range[m].lower) {
      i = m;
      goto done;
    }
    if (value < range[m].lower)
      r = m-1;
    else if (value > range[m].lower)
      l = m;
  }
  i = l;
done:
  
  // so now, if value is contained in the set, then
  // we must have: range[i].lower <= value
  if (value > range[i].upper)
    return false;

  // We must have range[i].lower <= value <= range[i].upper 
  return (((value-range[i].lower)%range[i].step) == 0);

}

//
// Function that returns true if value is contained
// within the current range.
bool Range::contains(const size_t value) const
{
  if (!all_but_mode) {
    return _contains(value);
  } else {
    return (lower_limit <= value &&
	    value < upper_limit &&
	    !_contains(value));
  }
}


#ifdef MAIN

static char *program_name;

static void
usage(const char* message = 0)
{
  if (message)
    fprintf(stderr, "%s: %s\n", program_name, message);
  fprintf(stderr,"Usage:\n");
  fprintf(stderr," -l # -u # -r rangespec [ -c ] [ -v ] [ -h ] [ -d ]\n");
  fprintf(stderr,"   where 'rangespec' is a range specification\n");
  fprintf(stderr,"   that must lie within the values indicated\n");
  fprintf(stderr,"   by the -l and -u flags.\n");
  fprintf(stderr,"   -v : verbose output\n");
  fprintf(stderr,"   -c : put commas between numbers\n");
  fprintf(stderr,"   -h : print range Horizontally\n");
  fprintf(stderr,"   -d : debug\n");
}


static long
parse_long(const char*const s)
{
    size_t len = strlen(s);
    char *ptr;
    long val;

    val = strtol(s, &ptr, 0);

    if (ptr != (s+len))
        error("Not an integer argument.");

    return val;
}

main(int argc,char *argv[])
{
  char *rng=0;
  size_t l=0,u=0;
  bool c = false;
  bool v = false;
  bool h = false;
  bool d = false;
  program_name = *argv++;
  argc--;  

  if (argc == 0) {
    usage("");
    exit(-1);
  }
  while (argc--) {
    char buf[BUFSIZ];
    const char *argp = *argv++;

    if (strcmp(argp, "-help") == 0)
      {
	usage();
	exit(0);
      }
    else if (strcmp(argp, "-c")==0)
      {
	c = true;
      }
    else if (strcmp(argp, "-v")==0)
      {
	v = true;
      }
    else if (strcmp(argp, "-h")==0)
      {
	h = true;
      }
    else if (strcmp(argp, "-d")==0)
      {
	d = true;
      }
    else if (strcmp(argp, "-l")==0)
      {
	if (argc>0)
	  {
	    l = (size_t) parse_long(*argv++);
	    argc--;
	  }
	else
	  usage("No -l value given.");
      }
    else if (strcmp(argp, "-u")==0)
      {
	if (argc>0)
	  {
	    u = (size_t) 1 + parse_long(*argv++);
	    argc--;
	  }
	else
	  usage("No -u value given.");
      }
    else if (strcmp(argp, "-r")==0)
      {
	if (argc>0)
	  {
	    rng = (*argv++);
	    argc--;
	  }
	else
	  usage("No range value given.");
      }
    else {
      sprintf(buf,"Unrecognized argument (%s).",argp);
      usage(buf);
      exit(-1);
    }
  }
  if (d)
    v = true;

  Range rg(rng,l,u);
  if (v) {
    printf("Range is (%s), a subset of [%u:%u]\n",rng,l,u);
    printf("Total Range Specifies %u Elements\n",rg.length());
    printf("First is at %u, Last is at %u\n",*(rg.begin()), *(rg.end()));
  }
  if (d) {
    printf("Iter\tCur_Val\tUp\tLow\tArr_Pos\n");
    for (Range::iterator it = rg.begin(); it <= rg.max();) {
      printf("%u\t%u\t%u\t%u\t%u\n",
	     *it,
	     it.cur_value,
	     it.cur_upper,it.cur_lower,it.array_pos);
      it++;
    }
    printf("Set defined is: ");
    size_t i = 0;
    for (Range::iterator it = rg.begin();
	 i<=rg.max()+10;i++) {
      if (rg.contains(i)) {
	if (it != i)
	  fprintf(stderr,"Error i = %u, it = %u\n",i,(size_t)it);
	printf("%u,",i);
	it++;
      }
    }
    printf("\n");
  } else {
    for (Range::iterator it = rg.begin(); it <= rg.max();) {
      printf("%u",*it);
      if (c && it < rg.max())
	printf(",");
      if (!h)
	printf("\n");
      it++;
    }
    if (h)
      printf("\n");
  }

}


#endif

