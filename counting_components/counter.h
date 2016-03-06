/* 

Defines class counter<T, Counting>, for use in measuring
the performance of certain STL generic algorithms.  Objects of this
class behave like those of type T with respect to assignments
and comparison operations, but the class also keeps counts of those
operations, using values of type Counting.

*/

/*
 * Copyright (c) 1997 Rensselaer Polytechnic Institute
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Rensselaer Polytechnic Institute makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */

#include <fstream>

#ifndef COUNTER_H
#define COUNTER_H

template <class T, class Counting>
class counter {
protected:
  T value;
public:
  static Counting assignments;
  static Counting less_comparisons;
  static Counting equal_comparisons;
  static Counting accesses;

  T base() const {++accesses; return value;}

  counter() : value(T()) { ++assignments; }

  counter(const T& v) : value(v) { ++assignments; }

  counter(const counter<T, Counting>& x) : value(x.value) { ++assignments; }

  static Counting total() {
    return assignments + less_comparisons + equal_comparisons + accesses; }
  
  static void reset() {
    assignments = 0;
    less_comparisons = 0;
    equal_comparisons = 0;
    accesses = 0;
  }

  static void report(std::ostream& o) {
      using namespace std;
    o << "Counter report:" << endl 
      << "  Accesses:  "    << accesses << endl
      << "  Assignments:  " << assignments << endl
      << "  Less comps:   " << less_comparisons << endl
      << "  Equality comps:" << equal_comparisons << endl;
  }

  friend bool operator<(const counter& x, 
                        const counter& y) 
  {
    ++less_comparisons;
    return x.value < y.value;
  }

  friend std::ostream& operator<<(std::ostream& o, const counter<T, Counting>& x) {
    return o << x.value;
  }

  counter<T, Counting>& operator=(const counter<T, Counting>& x) {
    ++assignments;
    value = x.value;
    return *this;
  }

  bool operator==(const counter<T, Counting>& x) const {
    ++equal_comparisons;
    return value == x.value;
  }
};

template <class T, class Counting>
Counting counter<T, Counting>::assignments = 0;

template <class T, class Counting>
Counting counter<T, Counting>::less_comparisons = 0;

template <class T, class Counting>
Counting counter<T, Counting>::equal_comparisons = 0;

template <class T, class Counting>
Counting counter<T, Counting>::accesses = 0;

#endif

