/* 

Defines class iteration_counter<RandomAccessIterator, T, Reference,
Distance, DistanceBase, Counting>, for use in measuring the
performance of certain STL generic algorithms.  Objects of this class
behave like those of type RandomAccessIterator, but the class also
keeps counts of all iterator operations, using values of type
Counting.  Type T should be the value type of RandomAccessIterator,
and Reference should be the reference type of T.  Type Distance should
be a distance type for RandomAccessIterator, and Counting is the type
used for the counts.

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
#include <iterator>
#include <fstream>

#ifndef ITERCOUNT_H
#define ITERCOUNT_H

template <class RandomAccessIterator, class T, class Reference,
          class Distance, class DistanceBase, class Counting>
          class iteration_counter : public std::iterator<std::random_access_iterator_tag, T, Distance> {
    typedef iteration_counter<RandomAccessIterator,
                              T, Reference, Distance,
                              DistanceBase, Counting>
            self;
    friend bool operator==(const self& x, const self& y);
    friend bool operator<(const self& x, const self& y);
    friend Distance operator-(const self& x, const self& y);
    friend self operator+(const Distance& n, const self& x);
    friend self operator+(const DistanceBase& n, const self& x);
protected:
    RandomAccessIterator current;
    Counting generation;
public:
    static Counting constructions;
    static Counting assignments;
    static Counting increments;
    static Counting dereferences;
    static Counting bigjumps;
    static Counting equality_comparisons;
    static Counting less_comparisons;
    static Counting differences;
    static Counting max_generation;

    static void reset() {
      constructions = 0;
      assignments = 0;
      increments = 0;
      dereferences = 0;
      bigjumps = 0;
      differences = 0;
      equality_comparisons = 0;
      less_comparisons = 0;
      max_generation = 0;
    }

    static Counting total() {
      return constructions + assignments + increments
              + dereferences + bigjumps + differences
              + equality_comparisons + less_comparisons;
    }

    static void report(std::ostream& o) {
      o << "Iterator stats: \n"
        << "  Constructions:  " << constructions << "\n"
        << "  Assignments:    " << assignments   << "\n"
        << "  Increments:     " << increments    << "\n"
        << "  Dereferences:   " << dereferences  << "\n"
        << "  Big jumps:      " << bigjumps      << "\n"
        << "  Differences     " << differences   << "\n"
        << "  Equality comps: " << equality_comparisons << "\n"
        << "  Less comps:     " << less_comparisons << "\n"
        << "  TOTAL:          " << total() << "\n";
      o << "  Maximum generation: " << max_generation << "\n";
    }
    iteration_counter() : generation(0) { ++constructions; }
    iteration_counter(RandomAccessIterator x) : current(x), generation(0) {
      ++constructions;
    }
    iteration_counter(const self& c) : current(c.current),
      generation(c.generation + 1) {
      ++constructions;
      if (generation > max_generation) {
        max_generation = generation;
      }
    }
    RandomAccessIterator base() { return current; }
    Reference operator*() const {
      ++dereferences;
      return *current;
    }
    self& operator++() {
        ++increments;
        ++current;
        return *this;
    }
    self operator++(int) {
        self tmp = *this;
        ++increments;
        ++current;
        return tmp;
    }
    self& operator--() {
        ++increments;
        --current;
        return *this;
    }
    self operator--(int) {
        self tmp = *this;
        ++increments;
        --current;
        return tmp;
    }
    self operator+(const Distance& n) const {
        ++bigjumps;
        return self(current + n);
    }
    self& operator+=(const Distance& n) {
        ++bigjumps; 
        current += n;
        return *this;
    }
    self operator-(const Distance& n) const {
        ++bigjumps;
        return self(current - n);
    }
    self& operator-=(const Distance& n) {
        ++bigjumps;
        current -= n;
        return *this;
    }
    Reference operator[](const Distance& n) { return *(*this + n); }
    self operator+(const DistanceBase& n) const {
        ++bigjumps;
        return self(current + n);
    }
    self& operator+=(const DistanceBase& n) {
        ++bigjumps;
        current += n;
        return *this;
    }
    self operator-(const DistanceBase& n) const {
        ++bigjumps;
        return self(current - n);
    }
    self& operator-=(const DistanceBase& n) {
        ++bigjumps;
        current -= n;
        return *this;
    }
    Reference operator[](const DistanceBase& n) { return *(*this + n); }
    friend std::ostream& operator << (std::ostream& o, const self& n) {
      return o << n.current;
    }
};

template <class RandomAccessIterator, class T, class Reference,
          class Distance, class DistanceBase, class Counting>
inline bool
  operator==(const iteration_counter<RandomAccessIterator, T,
                                     Reference, Distance, DistanceBase, Counting>& x,
             const iteration_counter<RandomAccessIterator, T,
                                     Reference, Distance, DistanceBase, Counting>& y) {
   ++iteration_counter<RandomAccessIterator,
                       T, Reference, Distance, DistanceBase,
                       Counting>::equality_comparisons;
   return x.current == y.current;
}

template <class RandomAccessIterator, class T, class Reference,
          class Distance, class DistanceBase, class Counting>
inline bool
  operator<(const iteration_counter<RandomAccessIterator, T,
                                    Reference, Distance, DistanceBase, Counting>& x,
            const iteration_counter<RandomAccessIterator, T,
                                    Reference, Distance, DistanceBase, Counting>& y) {
   ++iteration_counter<RandomAccessIterator,
                       T, Reference, Distance, DistanceBase, Counting>::less_comparisons;
   return x.current < y.current;
}

template <class RandomAccessIterator, class T, class Reference,
          class Distance, class DistanceBase, class Counting>
inline Distance
  operator-(const iteration_counter<RandomAccessIterator, T,
                                    Reference, Distance, DistanceBase, Counting>& x,
            const iteration_counter<RandomAccessIterator, T,
                                    Reference, Distance, DistanceBase, Counting>& y) {
   ++iteration_counter<RandomAccessIterator, T,
                       Reference, Distance, DistanceBase, Counting>::differences;
   return x.current - y.current;
}

template <class RandomAccessIterator, class T, class Reference,
          class Distance, class Counting>
inline iteration_counter<RandomAccessIterator, T, Reference,
                         Distance, class DistanceBase, Counting>
  operator+(const Distance& n,
            const iteration_counter<RandomAccessIterator, T, Reference,
                                    Distance, class DistanceBase, Counting>& x) {
   return iteration_counter<RandomAccessIterator, T,
                            Reference, Distance, DistanceBase, Counting>(x.current + n);
}

template <class RandomAccessIterator, class T, class Reference,
          class Distance, class Counting>
inline iteration_counter<RandomAccessIterator, T, Reference,
                         Distance, class DistanceBase, Counting>
  operator+(const DistanceBase& n,
            const iteration_counter<RandomAccessIterator, T, Reference,
                                    Distance, class DistanceBase, Counting>& x) {
   return iteration_counter<RandomAccessIterator, T,
                            Reference, Distance, DistanceBase, Counting>(x.current + n);
}

template <class RandomAccessIterator, class T, class Reference,
          class Distance, class DistanceBase, class Counting>
Counting iteration_counter<RandomAccessIterator, T, Reference,
                           Distance, DistanceBase, Counting>::
  constructions = 0;

template <class RandomAccessIterator, class T,
          class Reference,class Distance, class DistanceBase, class Counting>
Counting iteration_counter<RandomAccessIterator, T, Reference,
                           Distance, DistanceBase, Counting>::
  assignments = 0;

template <class RandomAccessIterator, class T,
          class Reference,class Distance, class DistanceBase, class Counting>
Counting iteration_counter<RandomAccessIterator, T, Reference,
                           Distance, DistanceBase, Counting>::
  increments = 0;

template <class RandomAccessIterator, class T,
          class Reference,class Distance, class DistanceBase, class Counting>
Counting iteration_counter<RandomAccessIterator, T, Reference,
                           Distance, DistanceBase, Counting>::
  dereferences = 0;

template <class RandomAccessIterator, class T,
          class Reference,class Distance, class DistanceBase, class Counting>
Counting iteration_counter<RandomAccessIterator, T, Reference,
                           Distance, DistanceBase, Counting>::
  bigjumps = 0;

template <class RandomAccessIterator, class T,
          class Reference,class Distance, class DistanceBase, class Counting>
Counting iteration_counter<RandomAccessIterator, T, Reference,
                           Distance, DistanceBase, Counting>::
  equality_comparisons = 0;

template <class RandomAccessIterator, class T,
          class Reference,class Distance, class DistanceBase, class Counting>
Counting iteration_counter<RandomAccessIterator, T, Reference,
                           Distance, DistanceBase, Counting>::
  less_comparisons = 0;

template <class RandomAccessIterator, class T,
          class Reference,class Distance, class DistanceBase, class Counting>
Counting iteration_counter<RandomAccessIterator, T, Reference,
                           Distance, DistanceBase, Counting>::
  differences = 0;

template <class RandomAccessIterator, class T,
          class Reference,class Distance, class DistanceBase, class Counting>
Counting iteration_counter<RandomAccessIterator, T, Reference,
                           Distance, DistanceBase, Counting>::
  max_generation = 0;

#endif
