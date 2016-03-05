

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

template <class Container>
void get(std::istream& is, Container& S) {
  S.erase(S.begin(), S.end());
  char ch;
  while (is.get(ch)) {
    if (ch == '\n')
      break;
    S.push_back(ch);
  }
}


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
void Report(algorithm_enumeration k, const Container& S1, 
            const Container& S2, const char* separator)
{
    using namespace std;
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


int main() 
{  
    using namespace std;
  ostream_iterator<char> out(cout, "");
  ifstream ifs("small.txt");
  string Comment, S1, S2;
  const char* separator = "";
  for (;;) {
    
    get(ifs, Comment);
    if (ifs.eof())
      break;
    copy(Comment.begin(), Comment.end(), out); cout << endl;
    
    get(ifs, S1);
    
    if (ifs.eof()) {
      cout << "**** Unexpected end of file." << endl;
      exit(1);
    }
    
    
    cout << "Text string:......";
    copy(S1.begin(), S1.end(), out);
    cout << endl;
    
    get(ifs, S2);
    
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
