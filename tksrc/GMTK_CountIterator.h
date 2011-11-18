#ifndef GMTK_COUNTITERATOR_H
#define GMTK_COUNTITERATOR_H

/*
 * simple count iterator that counts the number
 * of insertions made, but doesn't do anything else.
 */
template <typename _Container>
class count_iterator 
  : public iterator<output_iterator_tag, void, void, void, void> {
  unsigned counter;
public:

  count_iterator(_Container& __x) { counter = 0; }
  count_iterator() { counter = 0; }

  // count_iterator(const count_iterator& ci) { counter = ci.counter; }
  // count_iterator& operator=(const count_iterator& ci) { counter = ci.counter; }

  count_iterator& operator=(const typename _Container::const_reference _value) 
  { counter++; return *this; }
  count_iterator& operator*() { return *this; }
  count_iterator& operator++() {  return *this; }
  count_iterator& operator++(int) { return *this; }

  void reset() { counter = 0; }
  unsigned count() { return counter; }
};

#endif
