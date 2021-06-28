#include <iostream>
#include <random>
#include <chrono>

#include "JPMA.hpp"
#include <time.h>

#define InsertSize 107374182

using namespace std;

int main(int argc, char **argv){
    PMA pma;
    //int key;
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);
    //int inserted = 0;
    int *data;
    data =  (int *)malloc((InsertSize + 1) * sizeof(int));
    //int data[1073741825];
    srand (time(NULL));
    /*
    for(int i = 1; i< 100; i++){
        double d = distribution(generator);
        key = (int) (d*100 < 1? 1 : d*100);
        //cout<<"Inserting key "<<key;
        bool success = pma.insert(key, key*10);
        //if(success){ cout<<"inserted "<<key<<endl; inserted++;}
        //else cout<<"Key: "<<key<<" already present" << endl;
        //pma.printAllElements();
    }
    */
    for(int i = 0; i<= InsertSize; i++){
        data[i] = i;
    }
    cout<<"Data created"<<endl;
    
    for(int i = 0; i<=InsertSize; i++){
        int source, dest, buffer;
        source = rand() % InsertSize + 1;
        dest = rand() % InsertSize + 1;
        buffer = data[source];
        data[source] = data[dest];
        data[dest] = buffer;
    }
    cout<<"Data rearranged"<<endl;

    chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    start = chrono::high_resolution_clock::now();
    for(int i = 1; i<=InsertSize; i++){
        pma.insert(data[i], ((int64_t)data[i]) * 10, 0);
    }
    stop = chrono::high_resolution_clock::now();
    int64_t delay = chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    cout << "Total delay: "<<delay<<endl;
    //pma.printAllElements();
    //cout << "total inserted "<<inserted << endl;
    //cout <<"Looking 15: Found? "<< pma.lookup(105) <<endl;
    //cout <<"range query: "<<pma.range_sum(11, 50) <<endl;
    delete[] data;
    return 0;
}
