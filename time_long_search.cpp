

#define search stl_search
#define __search __stl_search
#include <algorithm>
#undef search
#undef __search


#include "new_search.h"
#include "hume.hh"
#include "DNA_search.h"
#include <iterator>
#include <deque>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
typedef unsigned char data;
#ifndef APCC_BUG
typedef std::vector<data> sequence;
#else
# define sequence vector<data>
#endif
sequence S1;
unsigned int Base_Line, Number_Of_Tests, Number_Of_Pattern_Sizes, Increment;
double Base_Time = 0.0;

enum algorithm_enumeration {
     Dummy, SF, L, HAL, ABM, TBM, GBM, HAL2, HAL3, HAL4, HAL5
};
const char* algorithm_names[] = {
     "selection code", "SF", "L", "HAL", "ABM", "TBM", "GBM", 
     "HAL2", "HAL3", "HAL4", "HAL5"
};
  const char textFileName[] = "long.txt";
  const char wordFileName[] = "words.txt";

#define DNA_TEST
#ifndef DNA_TEST
  algorithm_enumeration alg[] = {Dummy, TBM};
#else
  algorithm_enumeration alg[] = {Dummy, SF, L, HAL, ABM, GBM, 
                                 HAL2, HAL3, HAL4, HAL5};
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
  
  ifstream ifs(textFileName);
  char C;
  while (ifs.get(C))
    S1.push_back(C);
  cout << S1.size() << " characters read." << endl;
  

  
  ifstream dictfile(wordFileName);
  typedef istream_iterator<string> string_input;
  typedef map<int, vector<sequence>, less<int> > map_type;
  map_type dictionary;
  sequence S;
  string S0;
  string_input si(dictfile);
  while (si != string_input()) {
    S0 = *si++;
    S.erase(S.begin(), S.end());
    copy(S0.begin(), S0.end() - 1, back_inserter(S));
    dictionary[S.size()].push_back(S);
  }
  

  for (j = 0; j < Number_Of_Pattern_Sizes; ++j) {
    
    vector<sequence>& diction = dictionary[Pattern_Size[j]];
    if (diction.size() > Number_Of_Tests) {
      vector<sequence> temp;
      unsigned Skip_Amount = diction.size() / Number_Of_Tests;
      for (unsigned T = 0; T < Number_Of_Tests; ++T) {
         temp.push_back(diction[T * Skip_Amount]);
      }
      diction = temp;
    }
    

    Increment = (S1.size() - Pattern_Size[j]) / Number_Of_Tests;
    
    cout << "\n\n-----------------------------------------------------------\n"
         << "Searching for patterns of size " << Pattern_Size[j] 
         << "..." << endl;
    cout << "(" << Number_Of_Tests << " patterns from the text, "
         << dictionary[Pattern_Size[j]].size() << "  from the dictionary)" << endl;
    

    cerr << Pattern_Size[j] << " " << flush;
    
    Base_Time = 0.0;
    for (int k = 0; k < number_of_algorithms; ++k) {
      if (k != 0) 
        cout << "Timing " << algorithm_names[k] << ":" << endl;
      Run(k, S1, dictionary[Pattern_Size[j]], Pattern_Size[j]);
    }
    cout << endl;
    

  }
  cerr << endl;
}
