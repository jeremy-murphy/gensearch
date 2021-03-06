

#define search stl_search
#define __search __stl_search
#include <algorithm>
#undef search
#undef __search


#include "new_search.h"
#include <iterator>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
//#include <list>
//#define LIST_TEST
typedef std::string data;
typedef std::vector<data> sequence;
#if __STL_ITERATOR_TRAITS_NEEDED 
  ptr_iterator_traits(data);
#endif

sequence S1, S2;
unsigned Base_Line, Number_Of_Tests, Number_Of_Pattern_Sizes, Increment;
double Base_Time = 0.0;

struct search_word_trait {
    typedef std::vector<std::string>::const_iterator RAI;
  enum {hash_range_max = 256};
  enum {suffix_size = 1};
  inline static unsigned int hash(RAI i) {
    return (*i)[0];
  }
};



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



template <class Container>
void Run(int k, const Container& S1, 
         const std::vector<Container>& dictionary, int Pattern_Size)
{
    using namespace std;
  typename Container::const_iterator P;
  unsigned F = 0, d, K;
  double Start_Time, Finish_Time, Time_Taken;
  long Total_Search = 0;
  Start_Time = clock();
  Container S2;
  for (K = 1; K <= Number_Of_Tests; ++K) {
    typename Container::const_iterator u = S1.begin();
    advance(u, F);
    S2.erase(S2.begin(), S2.end());
    for (int I = 0; I < Pattern_Size; ++I)
      S2.push_back(*u++);
    F += Increment;
    
    Algorithm(k, S1, S2, P);
    d = distance(S1.begin(), P);
    Total_Search += d + Pattern_Size;
    

  }
  for (K = 0; K < dictionary.size(); ++K) {
    S2 = dictionary[K];
    
    Algorithm(k, S1, S2, P);
    d = distance(S1.begin(), P);
    Total_Search += d + Pattern_Size;
    

  }
  Finish_Time = clock();
  
  Time_Taken = (Finish_Time - Start_Time)/CLOCKS_PER_SEC - Base_Time;
  if (k == 0) 
    Base_Time = Time_Taken;  
  else {
    cout << "Total search length: " << Total_Search << " elements" << endl;
    cout << "Time: " << Time_Taken << " seconds." << endl;
    double Speed = Time_Taken == 0.0 ? 0.0 : 
      (double)Total_Search / 1000000 / Time_Taken;
    cout << "Speed: " << Speed << " elements/microsecond." << endl << endl;
  }
  

}


int main(int argc, char **argv) 
{  
    using namespace std;
  unsigned j;
/*  
  cout << "Input number of tests (for each pattern size): " << flush;
  cin >> Number_Of_Tests;
  cout << "Input number of pattern sizes: " << flush;
  cin >> Number_Of_Pattern_Sizes;
  cout << "Input pattern sizes: " << flush;
*/
  if (argc < 4)
      return 1;
  Number_Of_Tests = strtoul(argv[1], NULL, 10);
  Number_Of_Pattern_Sizes = strtoul(argv[2], NULL, 10);
  vector<unsigned> Pattern_Size(Number_Of_Pattern_Sizes);
  for (j = 0; j < Number_Of_Pattern_Sizes; ++j)
      Pattern_Size[j] = strtoul(argv[j + 3], NULL, 10);
  cout << "\nNumber of tests: " << Number_Of_Tests << endl;
  cout << "Pattern sizes: ";
  for (j = 0; j < Number_Of_Pattern_Sizes; ++j) 
    cout << Pattern_Size[j] << " ";
  cout << endl;
  

  typedef map<unsigned, vector<sequence >, less<unsigned> > map_type;
  map_type dictionary;
  
  ifstream ifs("long.txt");
  typedef istream_iterator<string> string_input;
  copy(string_input(ifs), string_input(), back_inserter(S1));
  

  cout << S1.size() << " words read." << endl;
  for (j = 0; j < Number_Of_Pattern_Sizes; ++j) {
    Increment = (S1.size() - Pattern_Size[j]) / Number_Of_Tests;
    
    cout << "\n\n-----------------------------------------------------------\n"
         << "Searching for patterns of size " << Pattern_Size[j] 
         << "..." << endl;
    cout << "(" << Number_Of_Tests << " patterns from the text, "
         << dictionary[Pattern_Size[j]].size() << "  from the dictionary)" << endl;
    

    
    Base_Time = 0.0;
    for (int k = 0; k < number_of_algorithms; ++k) {
      if (k != 0) 
        cout << "Timing " << algorithm_names[k] << ":" << endl;
      Run(k, S1, dictionary[Pattern_Size[j]], Pattern_Size[j]);
    }
    cout << endl;
    

  }
}
