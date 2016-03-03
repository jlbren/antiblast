#include <vector>
#include <bitset>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cstddef>
#include "GenomeHash.h"
#include "stdint.h"
#include <algorithm>
#include <iterator>

using namespace std;

GenomeHash::GenomeHash(char * file_path, char * query_path, int kmer, int missmatch, double threshold)
{
    query_name = query_path;
    file_name = file_path;
    word_size = kmer;
    num_mm = missmatch;
    cutoff = threshold;

    max_size=1;
    for(int i=0; i<word_size; i++) max_size=max_size*4;

    cout << "Reading: " << file_name << endl;
    if(get_words_host(file_name))
        Create();
    else
        cout<<"Error reading target genome"<<endl;
}


/*Generates words of size k from genome, populates array*/
bool GenomeHash:: get_words_host(char * file_name)
{
    char * word=new char [word_size+1];
    int i,j;
    unsigned int len;
    ifstream in;
    in.open(file_name);
    if(!in.is_open()) {
        cout << "The file could not be opened. Check the location.\n";
        return false;
    }

    //gets the length
    string str_word,line,header;
    getline(in,header);       //gets header
    len=0;
    while(in.peek()!=EOF)
    {
        getline(in,str_word);
        if(str_word[0]!='>')
            len+=str_word.size();
    }
    cout << "len=" << len << endl;
    subj_size=len;
    in.clear();
    in.close();

    htU = new NodeU [len];
    unsigned long long max_size=1;
    for(int i=0; i<word_size; i++)
        max_size=max_size*4;     //calculates = 4^word_size
    htS = new NodeS [max_size];

    int ii=0;

   in.open(file_name);
    getline(in,header);       //gets header
    cout << header << endl;
    len=0;
    char c;
    while(in.peek()!=EOF)
    {
        for(i=0; i<word_size; i++)
        {
            in>>c;
            if(c<97) c+=32;
            htU[ii].word[i]=c;
        }
        ii++;

        str_word="";
        while(!(in.peek()=='>' || in.peek()==EOF))
        {
            getline(in,line);
            str_word+=line;
            for(i=0;i<(str_word.length()-word_size+1);i++)
           {
                for(j=0;j<word_size-1;j++)
                    htU[ii].word[j]=htU[ii-1].word[j+1];
                c=str_word[i];
                if(c<97) c+=32;
                htU[ii].word[word_size-1]=c;
                ii++;
            }
            //adjust str_word
            j=str_word.length()-word_size+1;
            if(!(in.peek()=='>' || in.peek()==EOF || j<(word_size-1))) str_word.erase(0,j);
        }

        if(in.peek()!=EOF)
        {
            getline(in,header);       //gets header
            cout << header << endl;
        }
    }
    subj_size=ii;
    in.clear();
    in.close();
    return true;
}
//sets sizse of tables, fills unsorted array and calls hash
void GenomeHash::Create()
{
    Hash();

    compare_query_multi(query_name);

    cout << results.size() << " results found." << endl;
}



//Hashes host genome to populate hashmap
void GenomeHash::Hash()
{
    unsigned long long hash_val;
    for (int i = 0; i <subj_size; i++)
    {
        hash_val=word_to_decimal(htU[i].word);
        if(!(hash_val<0 || hash_val>max_size))
        {
            htS[hash_val].num_occ=true;
            htS[hash_val].ptr.push_back(&htU[i]);
        }
    }
}


/* Hash function: turn nuc into binary representation
*/
unsigned long long GenomeHash::word_to_decimal(string word)
{
    int word_size=word.length();
    unsigned long long val = 0; // 4 bytes =  32 bits = 16 chars
    char c;
    unsigned char temp;
    int i=0;

    c=word.at(i);

    switch(c)
    {
    case 'a':
        temp=0;
        break; //shifts 00 for a
    case 't':
        temp=1;
        break; //shifts 10 for t
    case 'c':
        temp=2;
        break;//shifts 01 for c
    case 'g':
        temp=3;
        break;//shifts 11 for g
    default: return -1;
    }
    val=temp;

    for (i=1; i<word_size; i++)
    {
        val=val<<2;
        c=word.at(i);
        switch(c)
        {
        case 'a':
            temp=0;
            break; //shifts 00 for a
        case 't':
            temp=1;
            break; //shifts 10 for t
        case 'c':
            temp=2;
            break;//shifts 01 for c
        case 'g':
            temp=3;
            break;//shifts 11 for g
        default: return -1;
        }
        val+=(int)temp;
    }

    if (val>max_size)
        cout <<"Out of range: " + word <<endl;

    return val;
}


//Returns values of all other nodes in hashmap within 1MM of arg sequence
vector<uint64_t> GenomeHash::is_in_1MM(UINT64 numericSeq)   //works
{
               UINT64 t;
               vector<UINT64> oneMM;

               unsigned char j1;
               unsigned char k1;
               unsigned char x1;
               UINT64 t1;
               UINT64 tt1;

               unsigned char value=0;


               for(j1=0;j1<word_size;j1++)
               {
                              t1=numericSeq;
                              t1=t1<<(62-(j1<<1));
                              t1=t1>>62;
                              x1=t1;
                              t1=t1<<(j1<<1);
                              tt1=ONE<<(j1<<1);

                              t=numericSeq-t1;

                              for(k1=0;k1<4;k1++)
                              {
                                             if(k1==x1)continue;                         //if k1= x1, that means this is the original sequence
                                             oneMM.push_back(t+k1*tt1);
                              }
               }
               /*for (int z = 0; z<oneMM.size(); z++)
                    cout<<oneMM[z]<<endl;*/ // test
               return oneMM;
}

//Sets chains between nodes generated by is_in_1mm
vector <uint64_t>  GenomeHash:: make_chains(vector <uint64_t> seqs, int mm)
{
    vector<uint64_t> chains;
    vector<uint64_t> temp;
    for (int i = 0; i < seqs.size(); i++){
        temp=is_in_1MM(seqs[i]);

        chains.insert(chains.end(), temp.begin(), temp.end());
       // cout<<chains.size();

    }
    if ( (mm-1)>0)
        make_chains(chains, (mm-1));

    return chains;
}


//Compare one or more query sequences against HT
bool GenomeHash:: compare_query_multi(char * query_path)
{
    query_name = query_path;
    unsigned long long current_word;
    ifstream in;
    in.open(query_name);
    int i,j,len;
    char * word=new char [word_size+1];

    cout << "Reading: " << query_name << endl;
    string str_word,line,header;
    getline(in,header);       //gets header
    cout << header << endl;

    while(in.peek()!=EOF)
    {
        len=0;
        str_word="";
        while( !(in.peek()=='>' || in.peek()==EOF))
        {
            getline(in,line);
            len+=line.size();
            str_word+=line;
        }

        for(i=0; i<word_size; i++)
        {
            word[i]=str_word[i];
            if(word[i]<97) word[i]+=32;                         //makes it lowercase
        }
        word[word_size]='\0';
        words.push_back(word);

        for(i=1; i<(len-word_size+1); i++)                      //read until the end of the file
        {
            //shift
            for(j=0; j<(word_size-1); j++) word[j]=word[j+1];
            word[(word_size-1)]=str_word[word_size+i-1];
            if(word[word_size-1]<97) word[word_size-1]+=32;     //makes it lowercase
            word[word_size]='\0';
            words.push_back(word);
        }

        query_size=words.size();

        for (uint64_t i=0; i<query_size;i++)                    // forward strand of query
        {
            UINT64 current_word = word_to_decimal(words[i]);
            i = get_match(header, current_word, i, num_mm, '+', cutoff );
        }

        //consider reverse complement of the query sequence
        get_reverse_complement();

        for (uint64_t i =0; i<query_size;i++)                   // reverse strand of query
        {
            UINT64 current_word = word_to_decimal(words[i]);
            i = get_match(header, current_word, i, num_mm, '-', cutoff);
        }
        if(in.peek()!=EOF)
        {
            getline(in,header);       //gets header
            cout << header << endl;
        }
    }

    in.clear();
    in.close();

    Get_Results();

    return true;
 }

/*Get numeric index from pointers used as chains between nodes*/
unsigned long long GenomeHash:: ptr_to_index(uint64_t word, int pointer_index)
{
    ptrdiff_t z;
    NodeU * base = htU;
    NodeU * s=htS[word].ptr[pointer_index];
    z=s-base;
    return (uint64_t)z;
}
/*Build match between query seed and host genome*/
unsigned long long GenomeHash:: get_match(string header, unsigned long long current_word, uint64_t i, int mm, char strand, double cutoff){ // pass in query index, word; return new index
    string match;
    double mismatches;
    double threshold = cutoff;
    double percent_mm=0;
    double match_size = 0;
    uint64_t start;
    uint64_t start_subj;
    uint64_t ii;
    uint64_t min_ii=query_size;

    ofstream out;
    out.open("dump.txt",ios::app);

    if(current_word<0 || current_word>max_size) return 0;

    if (htS[current_word].num_occ)
    {
        //loop through for each instance of query word in the subject sequence
        for(int np=0; np<htS[current_word].ptr.size(); np++)
        {
            NodeO temp;            uint64_t index = ptr_to_index(current_word,np);        // calculate index of current word in unsorted table of target genome
            start_subj=index;
            match = htU[index].word;
            mismatches=0;
            percent_mm=0;
            start=i;
            ii=i;
            //ii is position within the query (words), index is position within the subject

            //jump start => match=2*word_size
            while (((ii-start)<=word_size) && (ii<(query_size-1)) && (index<(subj_size-1)))
            {
                ii++;                                           // increment query position
                index +=1;                                      // increment target position
                match+=htU[index].word[word_size-1];            // add next char to match
                if (words[ii]!=htU[index].word)
                    mismatches++;
            }
            match_size=match.size();
            percent_mm = mismatches/match_size;

            //if jump started match meets the threshold, continue to grow the word
            if(percent_mm<=threshold)
            {                while ((percent_mm<=threshold) && (ii<(query_size-1)) && (index<(subj_size-1)))
                {
                    ii++;                                       // increment query position
                    index +=1;                                  // increment target position
                    match+=htU[index].word[word_size-1];        // add next char to match
                    match_size = match.size();
                    if (words[ii]!=htU[index].word)
                    {
                        mismatches++;
                        percent_mm = mismatches/match_size;
                    }
                }

                //if the match meets the user specified threshold (for length)
                if (match_size>=min_match_size_threshold)
                {
                    temp.q_pos=start;
                    temp.t_pos=start_subj;
                    temp.match=match;
                    temp.length=match_size;
                    temp.total_mm=mismatches;
                    temp.strand=strand;
                    temp.header=header;
                    results.push_back(temp);
                }

                out << match << "\t" << start << "\t" << start_subj << "\t" << mismatches << "\t" << match_size << endl;
                if(ii<min_ii) min_ii=ii;
            }
        }
    }

    out.clear();
    out.close();
    return ii;//min_ii;
}

//remove extraneous mismatches at the end of a match
void GenomeHash::roll_back_mm(int match_num)
{
    uint64_t match_length = results[match_num].length-1;
    uint64_t query_position = (results[match_num].q_pos+match_length)-1;

    uint64_t last_match;
    int mm = results[match_num].total_mm;
    string match = results[match_num].match;

    for(int i=match_length; i>=0; i--)
    {
        string q_word = words[query_position];//query[query_position];
        char match_nuc = match[i];
        char query_nuc = words[query_position][word_size-1];

        if(match[i]==words[query_position][word_size-1])
        {
            results[match_num].match=match.substr(0,i);
            results[match_num].length = results[match_num].match.length();
            return;
        }
        else
        {
            results[match_num].total_mm --;
            query_position--;
        }
    }
}

//write out the results
void GenomeHash:: Get_Results()
{
    ofstream out;
    out.open("results.txt");
    //reduce results to remove "sub-matches"
    remove_submatches();
    for(int i = 0; i < results.size(); i++)
    {
        //roll_back_mm(i);
        double mm = results[i].total_mm;
        double length = results[i].length;
        double total_match = length-mm;
        float percent_id = (total_match / length)*100;
        out<<">Query: "<<results[i].q_pos<<" Target: "<<results[i].t_pos<<" %ID: "<<total_match<<"/"<<length<<" = "<<percent_id<<" Strand: +/"<<results[i].strand<<results[i].header<<endl;
        out<<results[i].match<<endl;
    }
    out.clear();
    out.close();
}

/*
    Detect & remove submatches within greater matches*/

bool GenomeHash::remove_submatches()
{
    for(int i=0; i<results.size(); i++)
    {
        for(int j=0; j<results.size(); j++)
        {
            if(i!=j)
            {
                if( (results[j].q_pos>=results[i].q_pos) && (results[j].q_pos+results[j].length-1<=results[i].q_pos+results[i].length-1) &&
                    (results[j].t_pos>=results[i].t_pos) && (results[j].t_pos+results[j].length-1<=results[i].t_pos+results[i].length-1) &&
                    (results[j].strand==results[i].strand))
                {
                    results.erase(results.begin()+j);
                }
            }
        }
    }
    return true;
}

bool GenomeHash::check_word(unsigned long long index){
    if (htS[index].num_occ)
        return true;
    else return false;
}
/*Get reverse compliment of query words*/
void GenomeHash:: get_reverse_complement()
{
    vector<string> rc;
    string rev_word;
    for(int j=query_size-1; j>=0; j--)
    {
        rev_word=words[j];
        reverse(rev_word.begin(), rev_word.end());
        for(int i=0; i<rev_word.size(); i++)
        {
            switch(rev_word[i]){
            case 'a': rev_word[i]='t'; break;
            case 't': rev_word[i]='a'; break;
            case 'c': rev_word[i]='g'; break;
            case 'g': rev_word[i]='c'; break;
            default: rev_word[i]='n';
            }
        }
        rc.push_back(rev_word);
    }
    words.clear();
    words.swap(rc);
}
