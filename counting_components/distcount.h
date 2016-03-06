/*

Defines class distance_counter<RandomAccessIterator, Distance,
Counting>, for use in measuring the performance of certain STL generic
algorithms.  Objects of this class behave like those of type Distance
(which should be a type capable of representing the difference between
iterators of type RandomAccessIterator) but the class also keeps
counts of all Distance operations, using values of type Counting.
Type RandomAccessIterator is used as the result type of certain
addition and subtraction operators, with the assumption that Distance
is its distance type.

*/

#include <utility>
#include <fstream>
#include <memory>

#ifndef DISTCOUNT_H
#define DISTCOUNT_H

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

template <class RandomAccessIterator, class Distance, class Counting>
class distance_counter {
    typedef distance_counter<RandomAccessIterator, Distance, Counting>
            self;
    friend bool operator==(const self& x, const self& y);
    friend bool operator==(const self& x, const Distance& y);
    friend bool operator<(const self& x, const self& y);
    friend bool operator<(const self& x, const Distance& y);
    friend bool operator<(const Distance& x, const self& y);
    friend self operator+(const Distance& n, const self& x);
    friend self operator*(const Distance& n, const self& x);
    friend RandomAccessIterator
      operator+(RandomAccessIterator i, const self& x);
    friend RandomAccessIterator
      operator+=(RandomAccessIterator& i, const self& x);
    friend RandomAccessIterator
      operator-(RandomAccessIterator i, const self& x);
    friend RandomAccessIterator
      operator-=(RandomAccessIterator& i, const self& x);
protected:
    Distance current;
    Counting generation;
public:
    static Counting constructions;
    static Counting copy_constructions;
    static Counting conversions;
    static Counting assignments;
    static Counting increments;
    static Counting additions;
    static Counting subtractions;
    static Counting multiplications;
    static Counting divisions;
    static Counting equality_comparisons;
    static Counting less_comparisons;
    static Counting max_generation;

    static void reset() {
      constructions = 0;
      copy_constructions = 0;
      conversions = 0;
      assignments = 0;
      increments = 0;
      additions = 0;
      subtractions = 0;
      multiplications = 0;
      divisions = 0;
      equality_comparisons = 0;
      less_comparisons = 0;
      max_generation = 0;
    }

    static Counting total() {
      return constructions + copy_constructions + conversions
              + assignments + increments
              + additions + subtractions
              + multiplications + divisions
              + equality_comparisons + less_comparisons;
    }

    static void report(std::ostream& o) {
      o << "Distance stats: \n"
        << "  Constructions:   " << constructions << "\n"
        << "  Copies:          " << copy_constructions << "\n"
        << "  Conversions:     " << conversions   << "\n"
        << "  Assignments:     " << assignments   << "\n"
        << "  Increments:      " << increments    << "\n"
        << "  Additions:       " << additions     << "\n"
        << "  Subtractions:    " << subtractions  << "\n"
        << "  Multiplications: " << multiplications << "\n"
        << "  Divisions:       " << divisions     << "\n"
        << "  Equality comps:  " << equality_comparisons << "\n"
        << "  Less comps:      " << less_comparisons << "\n"
        << "  TOTAL:           " << total() << "\n";
      o << "  Maximum generation: " << max_generation << "\n";
    }
    distance_counter() : generation(0) { ++constructions; }
    distance_counter(const Distance& x) : current(x), generation(0) {
     ++conversions;
    }
    operator int() const { ++conversions; return current; }
    distance_counter(const self& c) : current(c.current),
      generation(c.generation + 1) {
      ++copy_constructions;
      if (generation > max_generation) {
        max_generation = generation;
      }
    }
    Distance base() const {return current; }
    self& operator=(const self& x) {
      ++assignments;
      current = x.current;
      return *this;
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
    self operator+(const self& n) const {
        self tmp = *this;
        return tmp += n;
    }
    self operator+(const Distance& n) const {
        self tmp = *this;
        return tmp += n;
    }
    self& operator+=(const self& n) {
        ++additions;
        current += n.current;
        return *this;
    }
    self& operator+=(const Distance& n) {
        ++additions;
        current += n;
        return *this;
    }
    self operator-(const self& n) const {
        self tmp = *this;
        return tmp -= n;
    }
    self operator-(const Distance& n) const {
        self tmp = *this;
        return tmp -= n;
    }
    self& operator-=(const self& n) {
        ++subtractions;
        current -= n.current;
        return *this;
    }
    self& operator-=(const Distance& n) {
        ++subtractions;
        current -= n;
        return *this;
    }
    self operator*(const self& n) const {
        self tmp = *this;
        return tmp *= n;
    }
    self operator*(const Distance& n) const {
        self tmp = *this;
        return tmp *= n;
    }
    self& operator*=(const self& n) {
        ++multiplications;
        current *= n.current;
        return *this;
    }
    self& operator*=(const Distance& n) {
        ++multiplications;
        current *= n;
        return *this;
    }
    self operator/(const self& n) const {
        self tmp = *this;
        return tmp /= n;
    }
    self operator/(const Distance& n) const {
        self tmp = *this;
        return tmp /= n;
    }
    self& operator/=(const self& n) {
        ++divisions;
        current /= n.current;
        return *this;
    }
    self& operator/=(const Distance& n) {
        ++divisions;
        current /= n;
        return *this;
    }
};

template <class RandomAccessIterator, class Distance, class Counting>
inline bool
  operator==(const
             distance_counter<RandomAccessIterator, Distance, Counting>& x,
             const
             distance_counter<RandomAccessIterator, Distance, Counting>& y) {
    ++distance_counter<RandomAccessIterator, Distance,
                       Counting>::equality_comparisons;
    return x.current == y.current;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline bool
  operator==(const
             distance_counter<RandomAccessIterator, Distance, Counting>& x,
             const Distance& y) {
    ++distance_counter<RandomAccessIterator, Distance,
                       Counting>::equality_comparisons;
    return x.current == y;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline bool
  operator<(const
              distance_counter<RandomAccessIterator, Distance, Counting>& x,
            const
              distance_counter<RandomAccessIterator, Distance, Counting>& y) {
    ++distance_counter<RandomAccessIterator, Distance,
                       Counting>::less_comparisons;
    return x.current < y.current;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline bool
  operator<(const
              distance_counter<RandomAccessIterator, Distance, Counting>& x,
            const Distance& y) {
    ++distance_counter<RandomAccessIterator, Distance,
                       Counting>::less_comparisons;
    return x.current < y;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline bool
  operator<(const Distance& x,
            const distance_counter<RandomAccessIterator, Distance, Counting>& y) {
    ++distance_counter<RandomAccessIterator, Distance,
                       Counting>::less_comparisons;
    return x < y.current;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline RandomAccessIterator
  operator+(RandomAccessIterator i,
            const
              distance_counter<RandomAccessIterator, Distance, Counting>& x) {
    return i + x.current;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline RandomAccessIterator
  operator+=(RandomAccessIterator& i,
             const
               distance_counter<RandomAccessIterator, Distance, Counting>& x) {
    return i += x.current;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline RandomAccessIterator
  operator-(RandomAccessIterator i,
            const
              distance_counter<RandomAccessIterator, Distance, Counting>& x) {
    return i - x.current;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline RandomAccessIterator
  operator-=(RandomAccessIterator &i,
             const
               distance_counter<RandomAccessIterator, Distance, Counting>& x) {
    return i -= x.current;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline distance_counter<RandomAccessIterator, Distance, Counting>
   operator-(const
               distance_counter<RandomAccessIterator, Distance, Counting>& x,
             const
               distance_counter<RandomAccessIterator, Distance, Counting>& y) {
     ++distance_counter<RandomAccessIterator, Distance,
                        Counting>::subtractions;
    return distance_counter<RandomAccessIterator, Distance,
                            Counting>(x.current - y.current);
}

template <class RandomAccessIterator, class Distance, class Counting>
inline distance_counter<RandomAccessIterator, Distance, Counting>
   operator+(const Distance& n,
             const
               distance_counter<RandomAccessIterator, Distance, Counting>& x) {
     ++distance_counter<RandomAccessIterator, Distance, Counting>::additions;
    return x + n;
}

template <class RandomAccessIterator, class Distance, class Counting>
inline distance_counter<RandomAccessIterator, Distance, Counting>
   operator*(const Distance& n,
             const
               distance_counter<RandomAccessIterator, Distance, Counting>& x) {
    ++distance_counter<RandomAccessIterator, Distance,
                       Counting>::multiplications;
    return x * n;
}

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  constructions = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  copy_constructions = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  conversions = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  assignments = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  increments = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  additions = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  subtractions = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  multiplications = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  divisions = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  equality_comparisons = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  less_comparisons = 0;

template <class RandomAccessIterator, class Distance, class Counting>
Counting distance_counter<RandomAccessIterator, Distance, Counting>::
  max_generation = 0;

/* The purpose of the following version of get_temporarary_buffer is to
   enable the STL stable_sort generic algorithms to work with
   distance_counter objects. */

template <class T, class RandomAccessIterator, class Distance, class Counting>
std::pair<T*,distance_counter<RandomAccessIterator, Distance, Counting> >
  get_temporary_buffer(distance_counter<RandomAccessIterator, Distance,
                       Counting> len, T* p) {
      std::pair<T*, int> tmp = get_temporary_buffer((int)len, p);
      return std::pair<T*, distance_counter<RandomAccessIterator, Distance,
                Counting> >(tmp.first,
                            distance_counter<RandomAccessIterator,
                                             Distance, Counting>(tmp.second));
}

#endif
