// 
// A general streams-like class with iterators for
// parsing and and using bi-polar range specifications.
//     For example: -3:10 means { -3, -2, -1, 0, ...., 10 }
// 
// Jeff Bilmes <bilmes@ee.washington.edu>
//
// $Header$
//



#include <limits.h>
#include <string.h>

#include "bp_range.h"
#include "error.h"


const int MAX_INT = INT_MAX;
const int MIN_INT = INT_MIN;

BP_Range::BP_Range(const char *range_str_a,
	     const int lower_limit_a,
	     const int upper_limit_a) 
  : lower_limit(lower_limit_a),upper_limit(upper_limit_a)
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
    error("Error upper_limit=%d < lower_limit=%d.",upper_limit,lower_limit);

  // if (upper_limit == lower_limit) {
  // empty range spec. 
  // }


  range_array_length = 20; 
  range_size = 0;
  range = new sub_range[range_array_length];
  parse_range();
  _max = (int)end();
  _min = (int)begin();


#if 0
  for (int i=0;i<range_size;i++) {
    printf ("[%d:%d:%d]\n",range[i].lower,range[i].step,
  	    range[i].upper);
  }
  printf("_max = %d, _min = %d\n",_max,_min);
#endif


}

BP_Range::~BP_Range()
{
  delete [] range;
  delete [] range_str;
}


bool
BP_Range::emptyRangeSpec(const char *str)
{
  if (str == NULL)
    return false;
  return (!strcmp(str,"nil")
	  || !strcmp(str,"none"));
}


bool
BP_Range::fullRangeSpec(const char *str)
{
  return (str == NULL 
	|| !strcmp(str,"all")
	|| !strcmp(str,"full")
	|| !strcmp(str,":"));
}




/*
 * Parse a range set specification.
 *
 * A (probably ambiguous) grammer is:
 *  range -> token | list
 *  token -> 'nil' | 'none' | 'all' | 'full' | '-' | ':'
 *  list -> [ brange range_list erange ]
 *  range_list -> range | [ rangelist ',' range ] | nil
 *  range -> num |  [num1 ':' num3 ':' num2 ]
 *  brange -> [ ':' num ] | [ ':' num3 ':' num2 ] | nil
 *  erange -> [ num ':' ] | [ num1 ':' num3 ':' ] | nil
 *  num1 -> numa   # num1 corresponds to the lower limit of a subrange
 *  num2 -> numa   # num2 corresponds to the upper limit of a subrange
 *  num3 -> num   # num3 corresponds to the step size of a subrange
 *  num -> [1:INT_MAX]
 *  numa -> [INT_MIN:INT_MAX] | ^[INT_MIN:INT_MAX]
 *  nil -> ""
 *
 *  Also
 *   num1 < num2
 *   'nil' - select no elements
 *   'none' - select no elements
 *   ':' - select all elements
 *   'all' - select all  elements
 *   'full' - select all elements
 *   
 * Examples:
 *   0,4:7,14,30:35,40:5:65
 * corresponds to:
 *   -5:4,5,6,7,14,30,31,32,32,34,35,40,45,50,55,60,65
 */
void BP_Range::parse_range()
{
  int i;
  const char *p = range_str;

  range_set_size = 0;
  if (p == NULL)
    error("Error: No range spec string argument given.");

  if (!strcmp(p,"nil") || !strcmp(p,"none")) {
    // empty range spec, select nothing.
    range[range_size].lower = 1;
    range[range_size].upper = 0;
    range[range_size].step = 1;
    range_size++;
    range_set_size = 0;
    return;
  } else if (!strcmp(p,"all") || !strcmp(p,"full") 
	     || !strcmp(p,"-") || !strcmp(p,":") ) {
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
    return;
  } else if (upper_limit == lower_limit) {
    // only empty range is possible but user gave non-empty range.
    error("Error: non-empty range specification given for a nul range.");
  }

  const char *next;
  int l,u;
  int s;

  while (*p) {
    // set to true if a ^n number is specified.
    bool from_top;

    if (*p == ':') {
      if (p == range_str) {
	// so we are at the beginning of the range string
	// and we've got a form ":..."
	l = lower_limit;
      } else {
	error("Error: invalid range spec string (%s) starting at (%s)",
	      range_str,p);
	l=0; // remove compiler warning.
      }
      next = p;
    } else {
      from_top = false;
      if (*p == '^') {
	p++;
	from_top = true;
      }
      l = (int) strtol(p,(char**) &next, 0);
      if (next == p) {
	error("Error: invalid range spec string (%s) starting at (%s)",
	      range_str,p);
      }
      if (from_top) {
	if (l < 0)
	  error("Error: invalid range spec string (%s) starting at (%s), negative ^n",range_str,p);
	l = upper_limit - l - 1;
      }
    }

    if (*next == ':') {
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
	if (*p == '^') {
	  p++;
	  from_top = true;
	}
	s = (int) strtol(p,(char**) &next, 0);
	if (next == p) {
	  error("Error: invalid range spec string (%s) starting at (%s)",
		range_str,p);
	}

	if (*next == ':') {
	  // then we indeed have "l:s:u" or "l:s:" form:

	  if (s < 1) {
	    // the s can't be negative.
	    error("Error: invalid negative range step size (%s) starting at (%s)",
		  range_str,p);
	  }

	  if (from_top) {
	    // the s can't have a '^' before it.
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
	    if (*p == '^') {
	      p++;
	      from_top = true;
	    }
	    u = (int) strtol(p,(char**)&next,0);
	    if (next == p) {
	      error("Error: invalid range spec string (%s) starting at (%s)",
		    range_str,p);
	    }
	    if (from_top) {
	      if (u < 0)
		error("Error: invalid range spec string (%s) starting at (%s), negative ^n",range_str,p);
	      u = upper_limit - u - 1;
	    }
	  }
	} else {
	  // we have a "l:u" form, s holds upper, step size 1
	  if (from_top) {
	    if (s < 0)
		error("Error: invalid range spec string (%s) starting at (%s), negative ^n",range_str,p);
	    u = upper_limit - s - 1;
	  }
	  else 
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
    error ("Error: lower %d > uupper %d in range string",range[i].lower,range[i].upper);
  if (range[i].lower < lower_limit)
    error ("Error: lower %d < limit %d in range string",range[i].lower,lower_limit);
  if (range[i].upper >= upper_limit)
    error ("Error: upper %d >= limit %d in range string",range[i].upper,upper_limit-1);

  for (i=1;i<range_size;i++) {
    if (range[i].lower > range[i].upper)
      error ("Error: lower %d > upper %d in range string",range[i].lower,range[i].upper);
    if (range[i].lower < lower_limit)
      error ("Error: lower %d < limit %d in range string",range[i].lower,lower_limit);
    if (range[i].upper >= upper_limit)
      error ("Error: upper %d >= limit %d in range string",range[i].upper,upper_limit-1);
    if (range[i-1].max() >= range[i].lower)
      error ("Error: previous upper %d >= current lower %d in range string",
	     range[i-1].upper,range[i].lower);
  }

}




void
BP_Range::make_range_safe_to_add_one()
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

BP_Range::iterator BP_Range::begin() const
{
  iterator rc(*this);
  return rc;
}

BP_Range::iterator BP_Range::end() const
{
  iterator rc(*this);
  rc.array_pos = range_size-1;

  rc.cur_lower = range[rc.array_pos].lower;
  rc.cur_upper = range[rc.array_pos].upper;
  rc.cur_value = range[rc.array_pos].max();

  return rc;
}




//********************************************************
//********************************************************
// BP_Range::iterator definitions
//********************************************************
//********************************************************

BP_Range::iterator::iterator(const BP_Range& rng)
  : myrange(rng)
{
  array_pos = 0;
  cur_value = myrange.range[array_pos].lower;
  cur_lower = myrange.range[array_pos].lower;
  cur_upper = myrange.range[array_pos].upper;
}

BP_Range::iterator::iterator(const BP_Range::iterator& it)
  : myrange(it.myrange)
{
  array_pos = it.array_pos;
  cur_value = it.cur_value;
  cur_upper = it.cur_upper;
  cur_lower = it.cur_lower;
}


BP_Range::iterator&
BP_Range::iterator::operator++()
{
  int s = myrange.range[array_pos].step;
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
      cur_value = MAX_INT;
    }
  }
  return *this;
}


BP_Range::iterator&
BP_Range::iterator::operator--()
{
  int s = myrange.range[array_pos].step;
  if (cur_value == MAX_INT)
    cur_value = myrange._max;
  else if (cur_value > cur_lower+s)
    cur_value -= s;
  else {
    if (array_pos > 0) {
      array_pos --;
      cur_lower = myrange.range[array_pos].lower;
      cur_upper = myrange.range[array_pos].upper;      
      const int s = myrange.range[array_pos].step;
      const int n = (cur_upper-cur_lower+s)/s;
      cur_value = cur_lower + s*(n-1);
    } else {
      // attempting to iterate off the bounds,
      // do nothing.
      // WARNING: This could lead to unexpected behavior, if
      // doing a loop such as: for (it=max;it>=min;it--).
    }
  }
  return *this;
}

bool BP_Range::contains(const int value) const
{

  // check for single value range case.
  if (min() == max())
    return (value == min());
  // check for out of bounds, or empty range.
  if (value < min() || value > max() || range_set_size == 0)
    return false;

  // Do a binary search that
  // finds the sub-range i such that
  // range[i].lower <= value < range[i+1].lower
  // or if (value == range[range_size-1].lower)
  // set i = range_size-1.
  int l,r,i;
  l = 0; r = range_size -1;
  while (l<r) {
    int m = (l+r+1)/2;
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
  if (range[i].step == 1)
    return true;
  else
    return (((value-range[i].lower)%range[i].step) == 0);
}




bool BP_Range::overlapP(const BP_Range& r) const
{
  //////////////////////////////////////////////////////////
  // TODO: This is an INEFFICIENT implementation
  // of this operation. This should be REDONE 
  // and optimized to have log time (linear time currently)

  // do simple tests first
  if (range_set_size == 0 || r.range_set_size == 0)
    return false;
  if (max() < r.min() || min() > r.max())
    return false;




  // So, there may be some overlap.
  // Be dumb, and earch over all range elements of r, stopping
  // as soon as we can. Do linear search on the smaller
  // one.
  if (r.length() < length()) {
    for (BP_Range::iterator it = r.begin();
	 it <=r.max();
	 it++) {
      const int val = (*it);
      if (val < min())
	continue;
      if (val > max())
	return false;
      if (contains(val)) {
	return true;
      }
    }
  } else {
    for (BP_Range::iterator it = begin();
	 it <= max();
	 it++) {
      const int val = (*it);
      if (val < r.min())
	continue;
      if (val > r.max())
	return false;
      if (r.contains(val)) {
	return true;
      }
    }    
  }
  return false;
}


bool BP_Range::operator <(const BP_Range& r) const
{
  if (max() < r.min())
    return true;
  if (max() < r.max() && min() < r.min())
    return true;
  return false;
}

bool BP_Range::operator <=(const BP_Range& r) const
{
  return (*this < r) || (*this == r);
}


bool BP_Range::operator >(const BP_Range& r) const
{
  if (min() > r.max())
    return true;
  if (min() > r.min() && max() > r.max())
    return true;
  return false;
}

bool BP_Range::operator >=(const BP_Range& r) const
{
  return (*this > r) || (*this == r);
}


bool BP_Range::operator == (const BP_Range& r) const
{
  if (min() != r.min() || max() != r.max() || length() != r.length())
    return false;
  BP_Range::iterator it1 = begin();
  BP_Range::iterator it2 = r.begin();
  for (;it1 <= max() && it2 <= r.max();it1++,it2++) {
    if ((*it1) != (*it2))
      return false;
  }
  return true;
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
    int len = strlen(s);
    char *ptr;
    long val;

    val = strtol(s, &ptr, 0);

    if (ptr != (s+len))
        error("Not an integer argument.");

    return val;
}

int
main(int argc,char *argv[])
{
  char *rng=0;
  int l=0,u=0;
  bool c = false;
  bool v = false;
  bool h = false;
  bool d = false;
  program_name = *argv++;
  argc--;  


  while (argc--) {
    char buf[BUFSIZ];
    const char *argp = *argv++;

    if (strcmp(argp, "-help")==0)
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
	    l = (int) parse_long(*argv++);
	    argc--;
	  }
	else
	  usage("No -l value given.");
      }
    else if (strcmp(argp, "-u")==0)
      {
	if (argc>0)
	  {
	    u = (int) parse_long(*argv++)+1;
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

  BP_Range rg(rng,l,u);
  if (v) {
    printf("Range is (%s), a subset of [%d:%d]\n",rng,l,u);
    printf("Total Range Specifies %d Elements\n",rg.length());
    printf("First is at %d, Last is at %d\n",*(rg.begin()), *(rg.end()));
  }
  if (d) {
    printf("Iter\tCur_Val\tUp\tLow\tArr_Pos\n");
    for (BP_Range::iterator it = rg.begin(); it <= rg.max();) {
      printf("%d\t%d\t%d\t%d\t%d\n",
	     *it,
	     it.cur_value,
	     it.cur_upper,it.cur_lower,it.array_pos);
      it++;
    }
    printf("Set defined is: ");
    int i = l-10;
    for (BP_Range::iterator it = rg.begin();
	 i<=rg.max()+10;i++) {
      if (rg.contains(i)) {
	if (it != i)
	  fprintf(stderr,"Error i = %d, it = %d\n",i,(int)it);
	printf("%d,",i);
	it++;
      }
    }
    printf("\n");
  } else {
    for (BP_Range::iterator it = rg.begin(); it <= rg.max();) {
      printf("%d",*it);
      if (c)
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

