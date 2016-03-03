
#line 2027 "gensearch.w"
#ifndef NEW_SEARCH
#  define NEW_SEARCH
#  include <vector.h>
#  include "search_traits.h"
#  include <iterator.h>
#  ifdef __STL_ITERATOR_TRAITS_NEEDED
     
#line 3532 "gensearch.w"
     template <class Iterator>
     struct iterator_traits {
       typedef Iterator::iterator_category iterator_category;
       typedef Iterator::value_type        value_type;
       typedef Iterator::difference_type   difference_type;
       typedef Iterator::pointer           pointer;
       typedef Iterator::reference         reference;
     };
     
     #define ptr_iterator_traits(T)                           \
     struct iterator_traits<T*> {                             \
       typedef random_access_iterator_tag iterator_category;  \
       typedef T                          value_type;         \
       typedef ptrdiff_t                  difference_type;    \
       typedef T*                         pointer;            \
       typedef T&                         reference;          \
     };                                                       \
                                                              \
     struct iterator_traits<const T*> {                       \
       typedef random_access_iterator_tag iterator_category;  \
       typedef T                          value_type;         \
       typedef ptrdiff_t                  difference_type;    \
       typedef const T*                   pointer;            \
       typedef const T&                   reference;          \
     }
     
     ptr_iterator_traits(char);
     ptr_iterator_traits(unsigned char);
     ptr_iterator_traits(int);
     ptr_iterator_traits(unsigned int);
     ptr_iterator_traits(short);
     ptr_iterator_traits(unsigned short);
     
#line 2033 "gensearch.w"

#  endif


#line 2332 "gensearch.w"
template <class RandomAccessIterator, class Distance>
void compute_next(RandomAccessIterator pattern, 
                  RandomAccessIterator patternEnd,
                  vector<Distance>& next)
{
  Distance pattern_size = patternEnd - pattern, j = 0, t = -1;
  next.reserve(32);
  next.push_back(-1);
  while (j < pattern_size - 1) {
    while (t >= 0 && pattern[j] != pattern[t]) 
       t = next[t];
    ++j; ++t;
    if (pattern[j] == pattern[t]) 
      next.push_back(next[t]);
    else
      next.push_back(t);
  }
}

#line 2036 "gensearch.w"


#line 1104 "gensearch.w"
template <class ForwardIterator, class Distance>
void compute_next(ForwardIterator pattern, 
                  ForwardIterator patternEnd,
                  vector<Distance>& next, 
                  vector<ForwardIterator>& pattern_iterator)
{
  Distance t = -1;
  next.reserve(32);
  pattern_iterator.reserve(32);
  next.push_back(-1);
  pattern_iterator.push_back(pattern);
  for (;;) {
    ForwardIterator advance = pattern;
    ++advance;
    if (advance == patternEnd)
      break;
    while (t >= 0 && *pattern != *pattern_iterator[t]) 
       t = next[t];
    ++pattern; ++t;
    if (*pattern == *pattern_iterator[t]) 
      next.push_back(next[t]);
    else
      next.push_back(t);
    pattern_iterator.push_back(pattern);
  }
}

#line 2037 "gensearch.w"


#line 1052 "gensearch.w"
template <class ForwardIterator1, class ForwardIterator2>
inline ForwardIterator1 search(ForwardIterator1 text,
                               ForwardIterator1 textEnd,
                               ForwardIterator2 pattern,
                               ForwardIterator2 patternEnd)
{
  typedef iterator_traits<ForwardIterator1> T;
  return __search(text, textEnd, pattern, patternEnd, T::iterator_category());
}

#line 2038 "gensearch.w"


#line 1066 "gensearch.w"
template <class ForwardIterator1, class ForwardIterator2>
inline ForwardIterator1 __search(ForwardIterator1 text,
                                 ForwardIterator1 textEnd,
                                 ForwardIterator2 pattern,
                                 ForwardIterator2 patternEnd,
                                 forward_iterator_tag)
{
  return __search_L(text, textEnd, pattern, patternEnd);
}

template <class ForwardIterator1, class ForwardIterator2>
ForwardIterator1 __search_L(ForwardIterator1 text,
                            ForwardIterator1 textEnd,
                            ForwardIterator2 pattern,
                            ForwardIterator2 patternEnd)
{
  typedef typename iterator_traits<ForwardIterator2>::difference_type Distance2;
  ForwardIterator1 advance, hold;
  ForwardIterator2 p, p1;
  Distance2 j, m;
  vector<Distance2> next;
  vector<ForwardIterator2> pattern_iterator;
  
#line 1100 "gensearch.w"
  compute_next(pattern, patternEnd, next, pattern_iterator);
  
#line 1088 "gensearch.w"

  m = next.size();
  
#line 1135 "gensearch.w"
  
#line 1149 "gensearch.w"
  if (next.size() == 1)
    return find(text, textEnd, *pattern);
  
#line 1135 "gensearch.w"
  
  p1 = pattern; ++p1;
  while (text != textEnd) {
    
#line 1158 "gensearch.w"
    while (*text != *pattern) 
      if (++text == textEnd)
        return textEnd;
    
#line 1138 "gensearch.w"
  
    
#line 1164 "gensearch.w"
    p = p1; j = 1;
    hold = text;
    if (++text == textEnd)
      return textEnd;
    while (*text == *p) {
      if (++p == patternEnd) 
        return hold;
      if (++text == textEnd)
        return textEnd;
      ++j;
    }
    
#line 1139 "gensearch.w"
  
    
#line 1178 "gensearch.w"
    for (;;) {
      j = next[j];
      if (j < 0) {
        ++text;
        break;
      }
      if (j == 0)
        break;
      p = pattern_iterator[j];
      while (*text == *p) {
        ++text; ++p; ++j;
        if (p == patternEnd) {
          
#line 1202 "gensearch.w"
          advance = hold;
          for (int i = m; --i >= 0;) 
            ++advance;
          while (advance != text) 
            ++advance, ++hold;
          return hold;
          
#line 1190 "gensearch.w"
    
        }
        if (text == textEnd)
          return textEnd;
      }
    }
    
#line 1140 "gensearch.w"
  
  }
  return textEnd;
  
#line 1090 "gensearch.w"

}

#line 2039 "gensearch.w"


#line 2175 "gensearch.w"
template <class BidirectionalIterator1, class BidirectionalIterator2>
inline BidirectionalIterator1 __search(BidirectionalIterator1 text,
                                       BidirectionalIterator1 textEnd,
                                       BidirectionalIterator2 pattern,
                                       BidirectionalIterator2 patternEnd,
                                       bidirectional_iterator_tag)
{
  return __search_L(text, textEnd, pattern, patternEnd);
}

#line 2040 "gensearch.w"


#line 2196 "gensearch.w"
template <class RandomAccessIterator1, class RandomAccessIterator2>
inline RandomAccessIterator1 __search(RandomAccessIterator1 text,
                                      RandomAccessIterator1 textEnd,
                                      RandomAccessIterator2 pattern,
                                      RandomAccessIterator2 patternEnd,
                                      random_access_iterator_tag)
{
  typedef iterator_traits<RandomAccessIterator1>::value_type V;
  typedef search_trait<V> Trait;
  return search_hashed(text, textEnd, pattern, patternEnd, (Trait*)0 ); 
}

#line 2041 "gensearch.w"


#line 2213 "gensearch.w"

template <class RandomAccessIterator1, class RandomAccessIterator2, class Trait>
RandomAccessIterator1 search_hashed(RandomAccessIterator1 text,
                                    RandomAccessIterator1 textEnd,
                                    RandomAccessIterator2 pattern,
                                    RandomAccessIterator2 patternEnd,
                                    Trait*)
{
  typedef typename iterator_traits<RandomAccessIterator1>::difference_type Distance1;
  typedef typename iterator_traits<RandomAccessIterator2>::difference_type Distance2;
  if (pattern == patternEnd) return text;
  Distance2 pattern_size, j, m;
  pattern_size = patternEnd - pattern; 
  if (Trait::suffix_size == 0 || pattern_size < Trait::suffix_size)
    return __search_L(text, textEnd, pattern, patternEnd);
  Distance1 i, k, large, adjustment, mismatch_shift, text_size;
  vector<Distance1> next, skip;
  
#line 2240 "gensearch.w"
  k = 0; 
  text_size = textEnd - text;
  
#line 2353 "gensearch.w"
  compute_next(pattern, patternEnd, next);
  
#line 2242 "gensearch.w"
  
  
#line 1149 "gensearch.w"
  if (next.size() == 1)
    return find(text, textEnd, *pattern);
  
#line 2243 "gensearch.w"
  
  
#line 2315 "gensearch.w"
  m = next.size();
  for (i = 0; i < Trait::hash_range_max; ++i)
    skip.push_back(m - Trait::suffix_size + 1);
  for (j = Trait::suffix_size - 1; j < m - 1; ++j)
    skip[Trait::hash(pattern + j)] = m - 1 - j;
  mismatch_shift = skip[Trait::hash(pattern + m - 1)];
  skip[Trait::hash(pattern + m - 1)] = 0;
  
#line 2244 "gensearch.w"
  
  large = text_size + 1;
  adjustment = large + pattern_size - 1;
  skip[Trait::hash(pattern + pattern_size - 1)] = large;
  k -= text_size;
  for(;;) {
    k += pattern_size - 1;
    if (k >= 0) break;
    
#line 2259 "gensearch.w"
    do {
      k += skip[Trait::hash(textEnd + k)]; 
    } while (k < 0);
    if (k < pattern_size)
      return textEnd;
    k -= adjustment;
    
#line 2252 "gensearch.w"
  
    
#line 2268 "gensearch.w"
    if (textEnd[k] != pattern[0])
      k += mismatch_shift;
    else {
      
#line 2280 "gensearch.w"
      j = 1; 
      for (;;) {
        ++k;
        if (textEnd[k] != pattern[j])
          break;
        ++j;
        if (j == pattern_size)
          return textEnd + k - pattern_size + 1;
      }
      
#line 2271 "gensearch.w"
    
      if (mismatch_shift > j)
        k += mismatch_shift - j;
      else
        
#line 2292 "gensearch.w"
        for (;;) {
          j = next[j];
          if (j < 0) {
            ++k;
            break;
          }
          if (j == 0)
            break;
          while (textEnd[k] == pattern[j]) {
            ++k; ++j;
            if (j == pattern_size) {
              return textEnd + k - pattern_size;
            }
            if (k == 0)
              return textEnd;
          }
        }
        
#line 2275 "gensearch.w"
    
    }
    
#line 2253 "gensearch.w"
  
  }
  return textEnd;
  
#line 2230 "gensearch.w"

}

#line 2042 "gensearch.w"

#endif
