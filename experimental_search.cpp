

#define search stl_search
#define __search __stl_search
#include <algorithm>
#undef search
#undef __search


#include "new_search.h"
#include "experimental_search.h"
#include <iterator>
#include <deque>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <ctime>

typedef unsigned short data;
#ifndef APCC_BUG
typedef std::vector<data> sequence;
#else
# define sequence vector<data>  
#endif
sequence S1;

int Base_Line, Number_Of_Tests, Number_Of_Pattern_Sizes, Increment;
double Base_Time = 0.0;

enum algorithm_enumeration {
     Dummy, SF, L, HAL, NHAL
};
const char* algorithm_names[] = {
     "selection code", "SF", "L", "HAL", "NHAL"
};

const int number_of_algorithms = 5;

template <class Container, class Container__const_iterator>
inline void
   Algorithm(int k, const Container& x, const Container& y, 
             Container__const_iterator& result)
{
  switch (algorithm_enumeration(k)) {
  case Dummy: 
     // does nothing, used for timing overhead of test loop
     result = x.begin(); return; 
  case SF: 
     result = stl_search(x.begin(), x.end(), y.begin(), y.end()); return;
  case L: 
     result =  __search_L(x.begin(), x.end(), y.begin(), y.end() ); return;
  case HAL: 
     result = search(x.begin(), x.end(), y.begin(), y.end() ); return;
  case NHAL: 
     result = search_no_hashing(x.begin(), x.end(), y.begin(), y.end() ); return;
  }
  result = x.begin(); return;
}



template <class Container>
void Run(int k, const Container& S1, 
         const std::vector<Container>& dictionary, int Pattern_Size)
{
    using namespace std;
  typename Container::const_iterator P;
  int F = 0, d, K;
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



int random(int max_value) { return rand() % max_value; }

template <int MAX_VALUE> struct RandomNumberGenerator {
  int operator() () { return random(MAX_VALUE); }
};



int main()
{ 
    using namespace std;
  int j;
  
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
  

  
  generate_n(back_inserter(S1), 100000, RandomNumberGenerator<65535>());
  

  
  typedef map<int, vector<sequence >, less<int> > map_type;
  map_type dictionary;
  
  for(int i = 0; i < Number_Of_Pattern_Sizes; ++i) {
    int pattern_size = Pattern_Size[i];
  
    for(int j = 0; j < Number_Of_Tests; ++j) {
      int position = random(S1.size() - pattern_size);
      dictionary[pattern_size].push_back( sequence() );
      copy(S1.begin() + position, S1.begin() + position + pattern_size, 
           back_inserter( dictionary[pattern_size].back() ) ) ;
    }
  }
  

  for (j = 0; j < Number_Of_Pattern_Sizes; ++j) {
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
