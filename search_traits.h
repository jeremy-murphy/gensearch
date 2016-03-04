
#ifndef SEARCH_HASH_TRAITS
#  define SEARCH_HASH_TRAITS
#  ifdef __STL_MEMBER_TEMPLATES
     
     template <class T>
     struct search_trait {
       enum {hash_range_max = 0};
       enum {suffix_size = 0};
       template <class RandomAccessIterator>
       inline static unsigned int hash(RandomAccessIterator i) {
         return 0;              
       }
     };
     

     
     template <> struct search_trait<signed char> {
       enum {hash_range_max = 256};
       enum {suffix_size = 1};
       template <class RandomAccessIterator>
       inline static unsigned int hash(RandomAccessIterator i) {
         return *i;              
       }
     };
     
     typedef unsigned char unsigned_char;
     template <> struct search_trait<unsigned_char> {
       enum {hash_range_max = 256};
       enum {suffix_size = 1};
       template <class RandomAccessIterator>
       inline static unsigned int hash(RandomAccessIterator i) {
         return *i;              
       }
     };
     

#  else
     
     #include <vector>
     #include <deque>
     
     template <class T>
     struct search_trait {
       enum {hash_range_max = 0};
       enum {suffix_size = 0};
     # define search_trait_helper_macro(Iterator)     \
       inline static unsigned int hash(Iterator i) {  \
         return 0;                                    \
       }
       search_trait_helper_macro(T*)
       search_trait_helper_macro(const T*)
       search_trait_helper_macro(typename std::deque<T>::iterator)
       search_trait_helper_macro(typename std::deque<T>::const_iterator)
     # undef search_trait_helper_macro
     };
     
     template <>
     struct search_trait<char> {
       enum {hash_range_max = 256};
       enum {suffix_size = 1};
     # define search_trait_helper_macro(Iterator)     \
       inline static unsigned int hash(Iterator i) {  \
         return *i;                                   \
       }
       search_trait_helper_macro(char*)
       search_trait_helper_macro(const char*)
       search_trait_helper_macro(typename std::deque<char>::iterator)
       search_trait_helper_macro(typename std::deque<char>::const_iterator)
     # undef search_trait_helper_macro
     };
     
     typedef unsigned char unsigned_char;
     template <>
     struct search_trait<unsigned_char> {
       enum {hash_range_max = 256};
       enum {suffix_size = 1};
     # define search_trait_helper_macro(Iterator)     \
       inline static unsigned int hash(Iterator i) {  \
         return *i;                                   \
       }
       search_trait_helper_macro(unsigned_char*)
       search_trait_helper_macro(const unsigned_char*)
       search_trait_helper_macro(typename std::deque<unsigned_char>::iterator)
       search_trait_helper_macro(typename std::deque<unsigned_char>::const_iterator)
     # undef search_trait_helper_macro
     };
     

#  endif
#endif
