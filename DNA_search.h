#ifdef __STL_MEMBER_TEMPLATES
  
  struct search_trait_dna2 {
    enum {hash_range_max = 64};
    enum {suffix_size = 2};
    template <class RAI>
    inline static unsigned int hash(RAI i) {
      return (*(i-1) + ((*i) << 3)) & 63;
    }
  };
  
  struct search_trait_dna3 {
    enum {hash_range_max = 512};
    enum {suffix_size = 3};
    template <class RAI>
    inline static unsigned int hash(RAI i) {
      return (*(i-2) + (*(i-1) << 3) + ((*i) << 6)) & 511;
    }
  };
  
  struct search_trait_dna4 {
    enum {hash_range_max = 256};
    enum {suffix_size = 4};
    template <class RAI>
    inline static unsigned int hash(RAI i) {
      return (*(i-3) + (*(i-2) << 2) + (*(i-1) << 4)
             + ((*i) << 6)) & 255;
    }
  };
  
  struct search_trait_dna5 {
    enum {hash_range_max = 256};
    enum {suffix_size = 5};
    template <class RAI>
    inline static unsigned int hash(RAI i) {
      return (*(i-4) + (*(i-3) << 2) + (*(i-2) << 4)
             + (*(i-1) << 6) + ((*i) << 8)) & 255;
    }
  };
  

#else
  
  typedef unsigned char unsigned_char;
  
  #define search_trait_helper_macro(Iterator)     \
    inline static unsigned int hash(Iterator i) {  \
      return (*(i-1) + ((*i) << 3)) & 63; }
  
  struct search_trait_dna2 {
    enum {hash_range_max = 64};
    enum {suffix_size = 2};
    search_trait_helper_macro(unsigned_char*)
    search_trait_helper_macro(const unsigned_char*)
    search_trait_helper_macro(deque<unsigned_char>::iterator)
    search_trait_helper_macro(deque<unsigned_char>::const_iterator)
  };
  #undef search_trait_helper_macro
  
  struct search_trait_dna3 {
    enum {hash_range_max = 512};
    enum {suffix_size = 3};
  # define search_trait_helper_macro(Iterator)             \
    inline static unsigned int hash(Iterator i) {          \
      return (*(i-2) + (*(i-1) << 3) + ((*i) << 6)) & 511; \
    }
    typedef unsigned char unsigned_char;
    search_trait_helper_macro(unsigned_char*)
    search_trait_helper_macro(const unsigned_char*)
    search_trait_helper_macro(deque<unsigned_char>::iterator)
    search_trait_helper_macro(deque<unsigned_char>::const_iterator)
  # undef search_trait_helper_macro
  };
  
  struct search_trait_dna4 {
    enum {hash_range_max = 256};
    enum {suffix_size = 4};
  # define search_trait_helper_macro(Iterator)             \
    inline static unsigned int hash(Iterator i) {          \
      return (*(i-3) + (*(i-2) << 2) + (*(i-1) << 4)       \
             + ((*i) << 6)) & 255;                         \
    }
    typedef unsigned char unsigned_char;
    search_trait_helper_macro(unsigned_char*)
    search_trait_helper_macro(const unsigned_char*)
    search_trait_helper_macro(deque<unsigned_char>::iterator)
    search_trait_helper_macro(deque<unsigned_char>::const_iterator)
  # undef search_trait_helper_macro
  };
  
  struct search_trait_dna5 {
    enum {hash_range_max = 256};
    enum {suffix_size = 5};
  # define search_trait_helper_macro(Iterator)             \
    inline static unsigned int hash(Iterator i) {          \
      return (*(i-4) + (*(i-3) << 2) + (*(i-2) << 4)       \
             + (*(i-1) << 6) + ((*i) << 8)) & 255;         \
    }
    typedef unsigned char unsigned_char;
    search_trait_helper_macro(unsigned_char*)
    search_trait_helper_macro(const unsigned_char*)
    search_trait_helper_macro(deque<unsigned_char>::iterator)
    search_trait_helper_macro(deque<unsigned_char>::const_iterator)
  # undef search_trait_helper_macro
  };
  

#endif

template <class RandomAccessIterator1, class RandomAccessIterator2>
inline RandomAccessIterator1 hal2(RandomAccessIterator1 text, 
                                  RandomAccessIterator1 textEnd,
                                  RandomAccessIterator2 pattern,
                                  RandomAccessIterator2 patternEnd)
{
  return search_hashed(text, textEnd, pattern, patternEnd,
                       (search_trait_dna2*)0);
}

template <class RandomAccessIterator1, class RandomAccessIterator2>
inline RandomAccessIterator1 hal3(RandomAccessIterator1 text, 
                                  RandomAccessIterator1 textEnd,
                                  RandomAccessIterator2 pattern,
                                  RandomAccessIterator2 patternEnd)
{
  return search_hashed(text, textEnd, pattern, patternEnd,
                       (search_trait_dna3*)0);
}
    
template <class RandomAccessIterator1, class RandomAccessIterator2>
inline RandomAccessIterator1 hal4(RandomAccessIterator1 text, 
                                  RandomAccessIterator1 textEnd,
                                  RandomAccessIterator2 pattern,
                                  RandomAccessIterator2 patternEnd)
{
  return search_hashed(text, textEnd, pattern, patternEnd,
                       (search_trait_dna4*)0);
} 

template <class RandomAccessIterator1, class RandomAccessIterator2>
inline RandomAccessIterator1 hal5(RandomAccessIterator1 text, 
                                  RandomAccessIterator1 textEnd,
                                  RandomAccessIterator2 pattern,
                                  RandomAccessIterator2 patternEnd)
{
  return search_hashed(text, textEnd, pattern, patternEnd,
                       (search_trait_dna5*)0);
}
