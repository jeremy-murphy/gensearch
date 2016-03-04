

struct large_alphabet_trait {
  typedef unsigned short T;
  enum {suffix_size = 1};
  enum {hash_range_max = (1u << (sizeof(T) * 8)) - 1};
};

#ifdef __STL_MEMBER_TEMPLATES
  template <> struct search_trait<unsigned short> {
    enum {hash_range_max = 256};
    enum {suffix_size = 1};
    template <class RandomAccessIterator>
    inline static unsigned int hash(RandomAccessIterator i) {
      return (unsigned char)(*i);
    }
  };
#else
  struct search_trait<unsigned short> {
    enum {hash_range_max = 256};
    enum {suffix_size = 1};
    inline static unsigned int hash(const unsigned short* i) {
      return (unsigned char)(*i);
    }
  };
#endif

template<class T>
class skewed_value {
  static T skew;
  T value;
public:
  skewed_value() : value(0) {}
  skewed_value(T val) : value(val - skew) {}
  operator T () { return value + skew; }
  static void setSkew(T askew) { skew = askew; }
  void clear() { value = 0; }
};

template<class T> T skewed_value<T>::skew;

template<class T, class RandomAccessIterator, int size>
class skewed_array {
  typedef skewed_value<T> value_type;
  static value_type array[size];
  RandomAccessIterator pattern, patternEnd;
public:
  skewed_array(T skew, RandomAccessIterator pat, RandomAccessIterator patEnd):
    pattern(pat),patternEnd(patEnd){ value_type::setSkew(skew); }
  ~skewed_array() {
    while (pattern != patternEnd) 
      array[*pattern++].clear();
  }
  value_type  operator[] (int index) const { return array[index]; }
  value_type& operator[] (int index)       { return array[index]; }
};

template<class T, class RandomAccessIterator, int size>
skewed_value<T> skewed_array<T,RandomAccessIterator,size>::array[size];

template <class RandomAccessIterator1, class RandomAccessIterator2>
RandomAccessIterator1 search_no_hashing(RandomAccessIterator1 text,
                                        RandomAccessIterator1 textEnd,
                                        RandomAccessIterator2 pattern,
                                        RandomAccessIterator2 patternEnd)
{
  typedef typename iterator_traits<RandomAccessIterator1>::difference_type Distance1;
  typedef typename iterator_traits<RandomAccessIterator2>::difference_type Distance2;
  typedef large_alphabet_trait Trait;
  if (pattern == patternEnd) 
    return text;
  Distance1 k, text_size, large, adjustment, mismatch_shift;
  Distance2 j, m, pattern_size;
  pattern_size = patternEnd - pattern;
  if (pattern_size < Trait::suffix_size)
    return __search_L(text, textEnd, pattern, patternEnd);
  vector<Distance1> next;
  skewed_array<Distance1, RandomAccessIterator2, Trait::hash_range_max+1>
    skip(pattern_size - Trait::suffix_size + 1, pattern, patternEnd);
  
  k = 0; 
  text_size = textEnd - text;
  
  compute_next(pattern, patternEnd, next);
  
  
  
  if (next.size() == 1)
    return find(text, textEnd, *pattern);
  
  
  
  m = next.size();
  for (j = Trait::suffix_size - 1; j < m - 1; ++j)
    skip[*(pattern + j)] = m - 1 - j;
  mismatch_shift = skip[*(pattern + m - 1)];
  skip[*(pattern + m - 1)] = 0;
  
  
  large = text_size + 1;
  adjustment = large + pattern_size - 1;
  skip[*(pattern + m - 1)] = large;
  
  k -= text_size;
  for (;;) {
    k += pattern_size - 1;
    if (k >= 0) break;
    
    do {
      k += skip[*(textEnd + k)]; 
    } while (k < 0);
    if (k < pattern_size)
      return textEnd;
    k -= adjustment;
    
  
    
    if (textEnd[k] != pattern[0])
      k += mismatch_shift;
    else {
      
      j = 1; 
      for (;;) {
        ++k;
        if (textEnd[k] != pattern[j])
          break;
        ++j;
        if (j == pattern_size)
          return textEnd + k - pattern_size + 1;
      }
      
    
      if (mismatch_shift > j)
        k += mismatch_shift - j;
      else
        
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
        
    
    }
    
  
  }
  return textEnd;
  

}


