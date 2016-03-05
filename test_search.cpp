

#define search stl_search
#define __search __stl_search
#include <algorithm>
#undef search
#undef __search


#include <iostream>
#include <fstream>
#include "new_search.h"
#include "hume.hh"
#include "DNA_search.h"
int Base_Line;


typedef unsigned char data;

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

template <typename ForwardIterator>
ForwardIterator Algorithm(int k, ForwardIterator first_x, ForwardIterator last_x, ForwardIterator first_y, ForwardIterator last_y)
{
  switch (alg[k]) {
  case Dummy: 
     // does nothing, used for timing overhead of test loop
     return first_x;
  case SF: 
     return stl_search(first_x, last_x, first_y, last_y);
  case L: 
     return  __search_L(first_x, last_x, first_y, last_y);
  case HAL: 
     return search(first_x, last_x, first_y, last_y);
  case ABM: 
     return fbm(first_x, last_x, first_y, last_y);
  case TBM: 
     return hume(first_x, last_x, first_y, last_y);
  case GBM: 
     return gdbm(first_x, last_x, first_y, last_y);
  case HAL2: 
     return hal2(first_x, last_x, first_y, last_y);
  case HAL3: 
     return hal3(first_x, last_x, first_y, last_y);
  case HAL4: 
     return hal4(first_x, last_x, first_y, last_y);
  case HAL5: 
     return hal5(first_x, last_x, first_y, last_y);
  }
  return first_x;
}


template <class Container>
void Report(algorithm_enumeration k, const Container& S1, 
            const Container& S2, const char* separator)
{
    using namespace std;
  typename Container::const_iterator P;
  P = Algorithm(k, S1.begin(), S1.end(), S2.begin(), S2.end());
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


int main() 
{  
    using namespace std;
  ostream_iterator<char> out(cout, "");
  ifstream ifs("small.txt");
  string Comment, S1, S2;
  const char* separator = "";
  for (;;) {
    
    getline(ifs, Comment);
    if (ifs.eof())
      break;
    copy(Comment.begin(), Comment.end(), out); cout << endl;
    
    getline(ifs, S1);
    if (ifs.eof()) {
      cout << "**** Unexpected end of file." << endl;
      exit(1);
    }
    
    cout << "Text string:......";
    copy(S1.begin(), S1.end(), out);
    cout << endl;
    
    getline(ifs, S2);
    
    if (ifs.eof()) {
      cout << "**** Unexpected end of file." << endl;
      exit(1);
    }
    
    cout << "Pattern string:...";
    copy(S2.begin(), S2.end(), out); cout << endl;
    
    Base_Line = 0;
    for (int k = 1; k < number_of_algorithms; ++k) {
      cout << "Using " << algorithm_names[k] << ":" << endl;
      Report(algorithm_enumeration(k), S1, S2, separator);
    }
    cout << endl;
  }
}
