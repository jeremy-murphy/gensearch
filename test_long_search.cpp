
#line 2834 "gensearch.w"

#line 2661 "gensearch.w"
#define search stl_search
#define __search __stl_search
#include <algo.h>
#undef search
#undef __search

#line 2834 "gensearch.w"

#include "new_search.h"
#include "hume.hh"
#include "DNA_search.h"
#include <iterator.h>
#include <vector.h>
#include <map.h>
#include <iostream.h>
#include <fstream.h>
#include <mstring.h>
typedef unsigned char data;
typedef vector<data> sequence;
sequence S1, S2;

int Base_Line, Number_Of_Tests, Number_Of_Pattern_Sizes, Increment;

#line 2702 "gensearch.w"
enum algorithm_enumeration {
     Dummy, SF, L, HAL, ABM, TBM, GBM, HAL2, HAL3, HAL4, HAL5
};
const char* algorithm_names[] = {
     "selection code", "SF", "L", "HAL", "ABM", "TBM", "GBM", 
     "HAL2", "HAL3", "HAL4", "HAL5"
};

#ifndef DNA_TEST
  algorithm_enumeration alg[] = {Dummy, SF, L, HAL, ABM, TBM};
  const char textFileName[] = "long.txt";
  const char wordFileName[] = "words.txt";
#else
  algorithm_enumeration alg[] = {Dummy, SF, L, HAL, ABM, GBM, 
                                 HAL2, HAL3, HAL4, HAL5};
  const char textFileName[] = "dnatext.txt";
  const char wordFileName[] = "dnaword.txt";
#endif

const int number_of_algorithms = sizeof(alg)/sizeof(alg[0]); 

template <class Container, class Container__const_iterator>
inline void
   Algorithm(int k, const Container& x, const Container& y, 
             Container__const_iterator& result)
{
  switch (alg[k]) {
  case Dummy: 
     // does nothing, used for timing overhead of test loop
     result = x.begin(); return;  
  case SF: 
     result = stl_search(x.begin(), x.end(), y.begin(), y.end()); return;
  case L: 
     result =  __search_L(x.begin(), x.end(), y.begin(), y.end()); return;
  case HAL: 
     result = search(x.begin(), x.end(), y.begin(), y.end()); return;
  case ABM: 
     result = fbm(x.begin(), x.end(), y.begin(), y.end()); return;
  case TBM: 
     result = hume(x.begin(), x.end(), y.begin(), y.end()); return;
  case GBM: 
     result = gdbm(x.begin(), x.end(), y.begin(), y.end()); return;
  case HAL2: 
     result = hal2(x.begin(), x.end(), y.begin(), y.end()); return;
  case HAL3: 
     result = hal3(x.begin(), x.end(), y.begin(), y.end()); return;
  case HAL4: 
     result = hal4(x.begin(), x.end(), y.begin(), y.end()); return;
  case HAL5: 
     result = hal5(x.begin(), x.end(), y.begin(), y.end()); return;
  }
  result = x.begin(); return;
}

#line 2849 "gensearch.w"


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

#line 2850 "gensearch.w"

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
  
#line 2854 "gensearch.w"

  
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
  
#line 2855 "gensearch.w"

  
#line 2914 "gensearch.w"
  ifstream ifs(textFileName);
  char C;
  while (ifs.get(C))
    S1.push_back(C);
  cout << S1.size() << " characters read." << endl;
  
#line 2856 "gensearch.w"

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
    
#line 2858 "gensearch.w"

    Increment = (S1.size() - Pattern_Size[j]) / Number_Of_Tests;
    cerr << Pattern_Size[j] << " " << flush;
    const char* separator = "";
    
#line 2922 "gensearch.w"
    cout << "\n\n-----------------------------------------------------------\n"
         << "Searching for patterns of size " << Pattern_Size[j] 
         << "..." << endl;
    cout << "(" << Number_Of_Tests << " patterns from the text, "
         << dictionary[Pattern_Size[j]].size() << "  from the dictionary)" << endl;
    
#line 2862 "gensearch.w"

    
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
    
#line 2863 "gensearch.w"

    
#line 2944 "gensearch.w"
    for (K = 0; K < dictionary[Pattern_Size[j]].size(); ++K) {
      S2 = dictionary[Pattern_Size[j]][K];
      
#line 2796 "gensearch.w"
      Base_Line = 0;
      for (int k = 1; k < number_of_algorithms; ++k) {
        cout << "Using " << algorithm_names[k] << ":" << endl;
        Report(algorithm_enumeration(k), S1, S2, separator);
      }
      cout << endl;
      
#line 2946 "gensearch.w"
    
    }
    
#line 2864 "gensearch.w"

  }
}
