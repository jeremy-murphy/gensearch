
#line 3383 "gensearch.w"

#line 2661 "gensearch.w"
#define search stl_search
#define __search __stl_search
#include <algo.h>
#undef search
#undef __search

#line 3383 "gensearch.w"

#include "new_search.h"
#include <iterator.h>
#include <vector.h>
#include <map.h>
#include <iostream.h>
#include <fstream.h>
#include <mstring.h>
typedef string data;
typedef vector<string> sequence;
#if __STL_ITERATOR_TRAITS_NEEDED 
  ptr_iterator_traits(data);
#endif
sequence S1, S2;
int Base_Line, Number_Of_Tests, Number_Of_Pattern_Sizes, Increment;

#line 3424 "gensearch.w"
struct search_word_trait {
  typedef vector<string>::const_iterator RAI;
  enum {hash_range_max = 256};
  enum {suffix_size = 1};
  inline static unsigned int hash(RAI i) {
    return (*i)[0];
  }
};

#line 3398 "gensearch.w"


#line 3435 "gensearch.w"
enum algorithm_enumeration {
     Dummy, STL_search, L, HAL
};
const char* algorithm_names[] = {
     "selection code", "SF", "L", "HAL"
};
#ifndef LIST_TEST
const int number_of_algorithms = 4;
#else
const int number_of_algorithms = 3;
#endif

template <class Container, class Container__const_iterator>
inline void
   Algorithm(int k, const Container& x, const Container& y, 
             Container__const_iterator& result)
{
  switch (algorithm_enumeration(k)) {
  case Dummy: 
     result = x.begin(); return;  // does nothing, used for timing overhead of test loop
  case STL_search: 
     result = stl_search(x.begin(), x.end(), y.begin(), y.end()); return;
  case L: 
     result = __search_L(x.begin(), x.end(), y.begin(), y.end()); return;
#ifndef LIST_TEST
  case HAL: 
     result = search_hashed(x.begin(), x.end(), y.begin(), y.end(), 
        (search_word_trait*)0); return;
#endif
  }
  result = x.begin(); return;
}

#line 3399 "gensearch.w"


#line 2805 "gensearch.w"
template <class Container>
void Report(algorithm_enumeration k, const Container& S1, 
            const Container& S2, const char* separator)
{
  typename Container::const_iterator P;
  Algorithm(k, S1, S2, P);
  cout << "  String " << '"';
  copy(S2.begin(), S2.end(), 
       ostream_iterator<typename Container::value_type>(cout, separator));
  if (P == S1.end())
    cout << '"' << " not found" << endl;
  else
    cout << '"' << " found at position " << P - S1.begin() << endl;
  if (Base_Line == 0)
    Base_Line = P - S1.begin();
  else
    if (P - S1.begin() != Base_Line)
      cout << "*****Incorrect result!" << endl;
}

#line 3400 "gensearch.w"

int main() 
{  
  int F, K, j;
  
#line 2870 "gensearch.w"
  cout << "Input number of tests (for each pattern size): " << flush;
  cin >> Number_Of_Tests;
  cout << "Input number of pattern sizes: " << flush;
  cin >> Number_Of_Pattern_Sizes;
  cout << "Input pattern sizes: " << flush;
  vector<int> Pattern_Size(Number_Of_Pattern_Sizes);
  for (j = 0; j < Number_Of_Pattern_Sizes; ++j)
    cin >> Pattern_Size[j];
  cout << "\nNumber of tests: " << Number_Of_Tests << endl;
  cout << "Pattern sizes: ";
  for (j = 0; j < Number_Of_Pattern_Sizes; ++j) 
    cout << Pattern_Size[j] << " ";
  cout << endl;
  
#line 3404 "gensearch.w"

  typedef map<int, vector<sequence >, less<int> > map_type;
  map_type dictionary;
  
#line 3470 "gensearch.w"
  ifstream ifs("long.txt");
  typedef istream_iterator<string, ptrdiff_t> string_input;
  copy(string_input(ifs), string_input(), back_inserter(S1));
  
#line 3407 "gensearch.w"

  cout << S1.size() << " words read." << endl;
  const char* separator = " ";
  for (j = 0; j < Number_Of_Pattern_Sizes; ++j) {
    Increment = (S1.size() - Pattern_Size[j]) / Number_Of_Tests;
    
#line 2922 "gensearch.w"
    cout << "\n\n-----------------------------------------------------------\n"
         << "Searching for patterns of size " << Pattern_Size[j] 
         << "..." << endl;
    cout << "(" << Number_Of_Tests << " patterns from the text, "
         << dictionary[Pattern_Size[j]].size() << "  from the dictionary)" << endl;
    
#line 3412 "gensearch.w"

    
#line 2930 "gensearch.w"
    F = 0;
    for (K = 1; K <= Number_Of_Tests; ++K) {
      
#line 2938 "gensearch.w"
      S2.erase(S2.begin(), S2.end());
      copy(S1.begin() + F, S1.begin() + F + Pattern_Size[j], back_inserter(S2));
      F += Increment;
      
#line 2932 "gensearch.w"
    
      
#line 2796 "gensearch.w"
      Base_Line = 0;
      for (int k = 1; k < number_of_algorithms; ++k) {
        cout << "Using " << algorithm_names[k] << ":" << endl;
        Report(algorithm_enumeration(k), S1, S2, separator);
      }
      cout << endl;
      
#line 2933 "gensearch.w"
    
    }
    
#line 3413 "gensearch.w"

  }
}
