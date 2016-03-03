
#line 3103 "gensearch.w"

#line 2661 "gensearch.w"
#define search stl_search
#define __search __stl_search
#include <algo.h>
#undef search
#undef __search

#line 3103 "gensearch.w"

#include "new_search.h"
#include "experimental_search.h"
#include <iterator.h>
#include <deque.h>
#include <vector.h>
#include <map.h>
#include <iostream.h>
#include <fstream.h>
#include <time.h>

typedef unsigned short data;
#ifndef APCC_BUG
  typedef vector<data> sequence;
#else
# define sequence vector<data>  
#endif
sequence S1;

int Base_Line, Number_Of_Tests, Number_Of_Pattern_Sizes, Increment;
double Base_Time = 0.0;

#line 3071 "gensearch.w"
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

#line 3124 "gensearch.w"


#line 3005 "gensearch.w"
template <class Container>
void Run(int k, const Container& S1, 
         const vector<Container>& dictionary, int Pattern_Size)
{
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
    
#line 3034 "gensearch.w"
    Algorithm(k, S1, S2, P);
    d = 0;
    distance(S1.begin(), P, d);
    Total_Search += d + Pattern_Size;
    
#line 3022 "gensearch.w"

  }
  for (K = 0; K < dictionary.size(); ++K) {
    S2 = dictionary[K];
    
#line 3034 "gensearch.w"
    Algorithm(k, S1, S2, P);
    d = 0;
    distance(S1.begin(), P, d);
    Total_Search += d + Pattern_Size;
    
#line 3026 "gensearch.w"

  }
  Finish_Time = clock();
  
#line 3051 "gensearch.w"
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
  
#line 3029 "gensearch.w"

}

#line 3125 "gensearch.w"


#line 3145 "gensearch.w"
int random(int max_value) { return rand() % max_value; }

template <int MAX_VALUE> struct RandomNumberGenerator {
  int operator() () { return random(MAX_VALUE); }
};

#line 3126 "gensearch.w"


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
  
#line 3131 "gensearch.w"

  
#line 3153 "gensearch.w"
  generate_n(back_inserter(S1), 100000, RandomNumberGenerator<65535>());
  
#line 3132 "gensearch.w"

  
#line 3157 "gensearch.w"
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
  
#line 3133 "gensearch.w"

  for (j = 0; j < Number_Of_Pattern_Sizes; ++j) {
    Increment = (S1.size() - Pattern_Size[j]) / Number_Of_Tests;
    
#line 2922 "gensearch.w"
    cout << "\n\n-----------------------------------------------------------\n"
         << "Searching for patterns of size " << Pattern_Size[j] 
         << "..." << endl;
    cout << "(" << Number_Of_Tests << " patterns from the text, "
         << dictionary[Pattern_Size[j]].size() << "  from the dictionary)" << endl;
    
#line 3136 "gensearch.w"

    cerr << Pattern_Size[j] << " " << flush;
    
#line 3041 "gensearch.w"
    Base_Time = 0.0;
    for (int k = 0; k < number_of_algorithms; ++k) {
      if (k != 0) 
        cout << "Timing " << algorithm_names[k] << ":" << endl;
      Run(k, S1, dictionary[Pattern_Size[j]], Pattern_Size[j]);
    }
    cout << endl;
    
#line 3138 "gensearch.w"

  }
  cerr << endl;
}
