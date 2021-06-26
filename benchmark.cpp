#include <iostream>
#include <random>

#include "JPMA.hpp"

using namespace std;

int main(int argc, char **argv){
    PMA pma;
    int key;
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);
    
    for(int i = 1; i< 100; i++){
        double d = distribution(generator);
        key = (int) (d*100 < 1? 1 : d*100);
        //cout<<"Inserting key "<<key;
        bool success = pma.insert(key, key*10);
        if(success) cout<<"inserted "<<key<<endl;
        else cout<<"Key: "<<key<<" already present" << endl;
        //pma.printAllElements();
    }
    pma.printAllElements();
    cout <<"Looking 15: Found? "<< pma.lookup(105) <<endl;
    cout <<"range query: "<<pma.range_sum(11, 50) <<endl;
    return 0;
}
