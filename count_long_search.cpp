
#line 3185 "gensearch.w"

#line 2661 "gensearch.w"
#define search stl_search
#define __search __stl_search
#include <algo.h>
#undef search
#undef __search

#line 3185 "gensearch.w"

#include "new_search.h"
#include "hume.hh"
#include <iterator.h>
#include <vector.h>
#include <map.h>
#include <iostream.h>
#include <fstream.h>
#include <time.h>
#include <mstring.h>

#line 3221 "gensearch.w"
#include "counter.h"
#include "itercount.h"
#include "distcount.h"
typedef unsigned char basedata;
typedef long counter_t;
typedef counter<basedata, counter_t> data;
// Should be able to use data in following definitions
// but there is a bug in apCC that prevents it
typedef  distance_counter<vector<counter<basedata, counter_t> >::const_iterator,
                          vector<counter<basedata, counter_t> >::difference_type, 
                                 counter_t>
         cdistance;
typedef iteration_counter<vector<counter<basedata, counter_t> >::const_iterator, 
                          const counter<basedata, counter_t>,
                          const counter<basedata, counter_t> &, cdistance,
                          vector<counter<basedata, counter_t> >::difference_type, 
                          counter_t>
        citer;

struct iterator_traits<citer> {
  typedef random_access_iterator_tag iterator_category;
  typedef data value_type;
  typedef cdistance difference_type;
  typedef data* pointer;
  typedef data&  reference;
};
#if __STL_ITERATOR_TRAITS_NEEDED 
  ptr_iterator_traits(data);
#endif

struct search_trait_for_counting {
  enum {hash_range_max = 256};
  enum {suffix_size = 1};
  inline static unsigned int hash(const citer& i) {return (*i).base();}
};


#line 3195 "gensearch.w"

typedef vector<data> sequence;
sequence S1;
int Base_Line, Number_Of_Tests, Number_Of_Pattern_Sizes, Increment;
double Base_Time = 0.0;

#line 3260 "gensearch.w"
enum algorithm_enumeration {
     Dummy, STL_search, L, HAL, ABM, TBM
};
const char* algorithm_names[] = {
     "selection code", "SF", "L", "HAL", "ABM", "TBM"
};
#ifndef LIST_TEST
const int number_of_algorithms = 6;
#else
const int number_of_algorithms = 3;
#endif

const char textFileName[] = "long.txt";
const char wordFileName[] = "words.txt";

template <class Container, class Container__const_iterator>
void Algorithm(int k, const Container& x, const Container& y,
               Container__const_iterator& result)
{
  switch (algorithm_enumeration(k)) {
  case Dummy: 
     result = x.begin();  // does nothing, used for timing overhead of test loop
     return;
  case STL_search: 
     result = stl_search(citer(x.begin()), citer(x.end()), 
                         citer(y.begin()), citer(y.end())).base();
     return;
  case L: 
     result = __search_L(citer(x.begin()), citer(x.end()), 
                         citer(y.begin()), citer(y.end())).base();
     return;

#ifndef LIST_TEST
  case HAL: 
     result = search_hashed(citer(x.begin()), citer(x.end()), 
                            citer(y.begin()), citer(y.end()), 
                            (search_trait_for_counting*)0).base();
     return;
  case ABM:
     fbmprep((const basedata*)y.begin(), y.size());
     result = (typename Container::const_iterator)
                fbmexec_cnt((const basedata*)x.begin(), x.size());
     data::accesses += ::pat.accs;
     data::equal_comparisons += ::pat.cmps;
     return;
  case TBM:
     humprep((const basedata*)y.begin(), y.size());
     result = (typename Container::const_iterator)
                humexec_cnt((const basedata*)x.begin(), x.size());
     data::accesses += ::pat.accs;
     data::equal_comparisons += ::pat.cmps;
     result = result;
     return;
#endif
  }
  result = x.begin();
  return;
}

#line 3200 "gensearch.w"


#line 3321 "gensearch.w"
template <class Container>
void Run(int k, const Container& S1, 
         const vector<Container>& dictionary, int Pattern_Size)
{
  typename Container::const_iterator P;
  int F = 0, K, d;
  double Start_Time, Finish_Time, Time_Taken;
  long Total_Search = 0;
  data::reset();
  citer::reset();
  cdistance::reset();
  Start_Time = clock();
  Container S2;
  for (K = 1; K <= Number_Of_Tests; ++K) {
    typename Container::const_iterator u = S1.begin();
    advance(u, F);
    S2.erase(S2.begin(), S2.end());
    for (int I = 0; I < Pattern_Size; ++I)
      S2.push_back(*u++);
    F += Increment;
    
#line 3034 "gensearch.w"
    Algorithm(k, S1, S2, P);
    d = 0;
    distance(S1.begin(), P, d);
    Total_Search += d + Pattern_Size;
    
#line 3341 "gensearch.w"

  }
  for (K = 0; K < dictionary.size(); ++K) {
    S2 = dictionary[K];
    
#line 3034 "gensearch.w"
    Algorithm(k, S1, S2, P);
    d = 0;
    distance(S1.begin(), P, d);
    Total_Search += d + Pattern_Size;
    
#line 3345 "gensearch.w"

  }
  Finish_Time = clock();
  
#line 3353 "gensearch.w"
  Time_Taken = (Finish_Time - Start_Time)/CLOCKS_PER_SEC - Base_Time;
  if (k == 0) 
    Base_Time = Time_Taken;  
  else {
    data::report(cout, Total_Search, 4);
    citer::report(cout, Total_Search, 4);
    cdistance::report(cout, Total_Search, 4);
    cout << "Total search length: " << Total_Search << " elements" << endl;
    cout << "Time: " << Time_Taken << " seconds." << endl;
    double Speed = Time_Taken == 0.0 ? 0.0 : 
      (double)Total_Search / 1000000 / Time_Taken;
    cout << "Speed: " << Speed << " elements/microsecond." << endl << endl;
  }
  
#line 3348 "gensearch.w"

}

#line 3201 "gensearch.w"


int main()
{ 
  int j;
  
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
  
#line 3206 "gensearch.w"

  
#line 2914 "gensearch.w"
  ifstream ifs(textFileName);
  char C;
  while (ifs.get(C))
    S1.push_back(C);
  cout << S1.size() << " characters read." << endl;
  
#line 3207 "gensearch.w"

  
#line 2886 "gensearch.w"
  ifstream dictfile(wordFileName);
  typedef istream_iterator<string, ptrdiff_t> string_input;
  typedef map<int, vector<sequence>, less<int> > map_type;
  map_type dictionary;
  sequence S;
  vector<char> S0;
  string_input si(dictfile);
  while (si != string_input()) {
    S0 = *si++;
    S.erase(S.begin(), S.end());
    copy(S0.begin(), S0.end() - 1, back_inserter(S));
    dictionary[S.size()].push_back(S);
  }
  
#line 3208 "gensearch.w"

  for (j = 0; j < Number_Of_Pattern_Sizes; ++j) {
    
#line 2902 "gensearch.w"
    vector<sequence>& diction = dictionary[Pattern_Size[j]];
    if (diction.size() > Number_Of_Tests) {
      vector<sequence> temp;
      int Skip_Amount = diction.size() / Number_Of_Tests;
      for (int T = 0; T < Number_Of_Tests; ++T) {
         temp.push_back(diction[T * Skip_Amount]);
      }
      diction = temp;
    }
    
#line 3210 "gensearch.w"

    Increment = (S1.size() - Pattern_Size[j]) / Number_Of_Tests;
    
#line 2922 "gensearch.w"
    cout << "\n\n-----------------------------------------------------------\n"
         << "Searching for patterns of size " << Pattern_Size[j] 
         << "..." << endl;
    cout << "(" << Number_Of_Tests << " patterns from the text, "
         << dictionary[Pattern_Size[j]].size() << "  from the dictionary)" << endl;
    
#line 3212 "gensearch.w"

    cerr << Pattern_Size[j] << " " << flush;
    
#line 3041 "gensearch.w"
    Base_Time = 0.0;
    for (int k = 0; k < number_of_algorithms; ++k) {
      if (k != 0) 
        cout << "Timing " << algorithm_names[k] << ":" << endl;
      Run(k, S1, dictionary[Pattern_Size[j]], Pattern_Size[j]);
    }
    cout << endl;
    
#line 3214 "gensearch.w"

  }
  cerr << endl;
}
