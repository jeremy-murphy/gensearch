
#ifndef SEARCH_HASH_TRAITS
#  define SEARCH_HASH_TRAITS

// TODO: The return type of the hash function should be the unsigned equivalent
// of the iterator's value type.

     template <class T>
     struct search_trait {
       enum {hash_range_max = 0};
       enum {suffix_size = 0};
       template <class RandomAccessIterator>
       inline static 
       typename std::iterator_traits<RandomAccessIterator>::value_type
       hash(RandomAccessIterator) {
         return 0;              
       }
     };

     template <> struct search_trait<char> {
         enum {hash_range_max = 256};
         enum {suffix_size = 1};
         template <class RandomAccessIterator>
         inline static 
         typename std::iterator_traits<RandomAccessIterator>::value_type
         hash(RandomAccessIterator i) {
             return *i;              
         }
     };

     template <> struct search_trait<signed char> {
       enum {hash_range_max = 256};
       enum {suffix_size = 1};
       template <class RandomAccessIterator>
       inline static 
       typename std::iterator_traits<RandomAccessIterator>::value_type
       hash(RandomAccessIterator i) {
         return *i;              
       }
     };
     
     template <> struct search_trait<unsigned char> {
       enum {hash_range_max = 256};
       enum {suffix_size = 1};
       template <class RandomAccessIterator>
       inline static 
       typename std::iterator_traits<RandomAccessIterator>::value_type
       hash(RandomAccessIterator i) {
         return *i;              
       }
     };
     

#endif

