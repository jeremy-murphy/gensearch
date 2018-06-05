#include <boost/algorithm/searching/aho_corasick.hpp>
#include <boost/algorithm/searching/boyer_moore.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using boost::iostreams::mapped_file;

namespace std
{
    ostream &operator<<(ostream &out, pair<string, unsigned> const &p)
    {
        out << p.first << ": " << p.second;
        return out;
    }
}

int main(int argc, char **argv)
{
    if (argc <= 2)
    {
        return 1;
    }
    
    using namespace std::placeholders;
    
    ifstream dictionary(argv[1]);
    mapped_file corpus(argv[2], mapped_file::readonly);
    typedef mapped_file::const_iterator corpus_iterator;
    string word;
    unordered_map<string, unsigned> word_count;
    /*
    while (dictionary >> word)
        word_count[word];
    dictionary.seekg(0);
    */
    auto const count_word = [&](corpus_iterator first, corpus_iterator last){ word_count[string(first, last)]++; return true; };
    boost::algorithm::aho_corasick_base<char, boost::container::flat_map> foo;
    boost::algorithm::aho_corasick_search<char>(corpus.const_begin(), corpus.const_end(), istream_iterator<string>(dictionary), istream_iterator<string>(), count_word);
    // cout << "word_count.size(): " << word_count.size() << endl;
    std::vector<pair<string, unsigned>> tmp(begin(word_count), end(word_count));
    // sort(begin(tmp), end(tmp));
    copy(begin(tmp), end(tmp), ostream_iterator<pair<string, unsigned>>(cout, "\n"));
}
