#include <vector>
#include <bitset>
#include <cstring>
#include <fstream>
#include <iostream>
#include "stdint.h"

using namespace std;

class GenomeHash{
private:
    int word_size;
    unsigned long long max_size;
    int subj_size;
    int query_size;
    int num_mm;
    int min_match_size_threshold;
    double cutoff;
    typedef unsigned long long UINT64;
    static const UINT64 ONE=1;
    struct NodeU{char word[13];};//string word;};
    struct NodeS{bool num_occ=0; vector<NodeS*> MM; vector<NodeU*> ptr;};
    struct NodeO{string header; string match; int length; double percent_ident;int total_mm; unsigned long long q_pos; unsigned long long t_pos; char strand;};
    bool remove_submatches();
    NodeU * htU;                // list by appearance
    NodeS * htS;                // list by hash
    vector <NodeO> results;
    vector<string> words;       // vector of sequence as words

public:
    char * file_name;
    char * query_name;
    unsigned long long total_size;
    GenomeHash(char * file_name,  char * query_name, int kmer, int missmath, double cutoff);
    bool get_words(char * file_name);
    bool get_words_host(char * file_name);
    unsigned long long word_to_decimal(string word);
    void Create();
    void Hash();
    vector<uint64_t> is_in_1MM(UINT64 numericSeq);
    bool compare_query_multi(char * query_path);
    bool compare_query(char * query_path);
    /*bool is_in_2MM(UINT64_t numericSeq);
    bool is_in_3MM(UINT64_t numericSeq);*/
    unsigned long long get_match(string header, unsigned long long current_word, uint64_t i, int mm, char strand, double cutoff); // pass in query index, word; return new index
    vector <uint64_t> make_chains(vector <uint64_t> seqs, int mm);
    unsigned long long ptr_to_index(unsigned long long word, int pointer_index);
    void roll_back_mm(int match_num);
    void Get_Results();
    bool check_word(unsigned long long index);
    void get_reverse_complement(); // take query and add reverse compliment kmers

 };
