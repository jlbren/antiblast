#include <iostream>
#include <time.h>
#include "GenomeHash.h"

using namespace std;
/*todo:
    GUI
        java applets
        command line arg make file input c++
        python wrapper
    Bugs:
        kkmer sizes > 13



*/

int main (int argc, char *argv[])
{
    clock_t t1, t2;
    t1 = clock();


    char * file_name = argv[0];
    char * query_name = argv[1];


    int word_size = argv[2];// breaking at >13 && <8 5/21
    bool multi_fata = true;
    double cutoff =argv[3];

    int mm = 0;//2;
    GenomeHash * Oddish = new GenomeHash(file_name, query_name, word_size, mm, cutoff);
    t2 = clock();
    float diff = t2-t1;
    float secs = diff/CLOCKS_PER_SEC;

    ofstream out;
    out.open("results.txt",  fstream::in | fstream::out | fstream::app);
    out<<'\n'<<"#RUN TIME: "<<secs;
    //system ("pause");

    return 0;

}
