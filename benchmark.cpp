#include <iostream>
#include <random>
#include <chrono>
#include <cstring>

#include<fstream>

#include "JPMA_BT.hpp"
#include <time.h>

#define InsertSize 10737418

using namespace std;

void printArguments(){
    cout<<"USAGE: ./benchmark [options]"<<endl;
    cout<<"Options:"<<endl;
    cout<<"    -i [int]     number of key-value pairs to insert"<<endl;
    cout<<"    -d [int]     number of key-value pairs to delete"<<endl;
    cout<<"    -r [double]  length of range for sacnning expressed as percentage of inserted keys (0 < r <= 1)"<<endl;
    cout<<"    -rr[int]     number of repeating for range scan queries "<<endl;
    cout<<"    -s [int]     number of key-value pairs to search"<<endl;
    cout<<endl;
}

int main(int argc, char **argv){
    //Redirect cout to file out.txt
    //std::ofstream out("out.txt");
    //std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    //std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    chrono::time_point<std::chrono::high_resolution_clock> start, stop;

    if (argc == 1) {
        printArguments();
        return 1;
    }

    type_t totalInsert = 0;
    type_t totalDelete = 0;
    double rangeRatio = -1;
    type_t rangeLength = 0;
    type_t rangeIteration = 0;
    type_t totalSearch = 0;

    for (type_t i = 1; i<argc; i++) {
        if(strcmp(argv[i], "-i") == 0) {
            totalInsert = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-d") == 0) {
            totalDelete = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-r") == 0) {
            rangeRatio = (atof(argv[++i]));
        } else if (strcmp(argv[i], "-rr") == 0) {
            rangeIteration = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-s") == 0) {
            totalSearch = atoi(argv[++i]);
        } else {
            printArguments();
            return 1;
        }
    }

    if(totalInsert == 0){
        totalInsert = InsertSize;
    }
    if(rangeRatio>0){
        rangeLength = (type_t)(totalInsert * rangeRatio);
    }
    if(rangeIteration > 0 && rangeLength<0) {
        cout<<"Range Iteration without specifying range length is ambigous"<<endl;
        exit(1);
    }

    PMA pma(totalInsert);

    if(totalInsert < rangeLength) {
        cout<<"Range length greater than total elements"<<endl;
        exit(1);
    }
    if(totalInsert < totalDelete) {
        cout<<"Number of elements to be deleted greater than total elements"<<endl;
        exit(1);
    }

    int64_t *data;
    if(Insert_type == 1){
        int insertCount = 128, inserted = 0;
        type_t insertDelay = 0;
        while(inserted+insertCount <= totalInsert){
            int64_t *data;
            data =  (int64_t *)malloc((insertCount+1) * sizeof(int64_t));
            srand (time(NULL));
            for(int i = 0; i< insertCount; i++){
                data[i] = i+inserted+1;
                //data[i] = rand() % 1000000;
            }

            for(int i = 0; i<insertCount; i++){
                int64_t source, dest, buffer;
                source = rand() % insertCount;
                dest = rand() % insertCount;
                buffer = data[source];
                data[source] = data[dest];
                data[dest] = buffer;
            }

            start = chrono::high_resolution_clock::now();
            for(int i = 0; i<insertCount; i++){
                pma.insert(data[i], data[i] * 10);
            }
            stop = chrono::high_resolution_clock::now();
            int64_t delay = chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
            if(inserted == 0){
                inserted = insertCount;
            }else{
                inserted += insertCount;
                insertCount *= 2;
            }
            delete data;
            insertDelay += delay;
        }
        cout<<"Inserted "<<(inserted+insertCount)<<" elements in "<<insertDelay<<" microSeconds."<<endl;
    }
    else{
        //int64_t *data;
        data = (int64_t *)malloc((totalInsert+1) * sizeof(int64_t));
        for(int i=0; i<totalInsert; i++) data[i] = i+1;
        srand(time(NULL));
        for(int i = 0; i<totalInsert; i++){
            int64_t source, dest, buffer;
            source = rand() % totalInsert;
            dest = rand() % totalInsert;
            buffer = data[source];
            data[source] = data[dest];
            data[dest] = buffer;
        }
        start = chrono::high_resolution_clock::now();
        for(int i = 0; i<totalInsert; i++){
            if(!pma.insert(data[i], data[i] * 10))
                cout<<"Could not insert: "<<data[i]<<endl;
        }
        stop = chrono::high_resolution_clock::now();
        type_t insertDelay = chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        cout<<"Inserted "<<totalInsert<<" elements in "<<insertDelay<<" microSeconds."<<endl;
        pma.printStat();
    }
    
    //Searching in the PMA
    if(totalSearch>0){
        int notFound = 0;
        start = chrono::high_resolution_clock::now();
        for(type_t i=0; i<totalSearch; i++){
            if(!pma.lookup(data[i])){
                notFound++;
                cout<<"Could not get key: "<<data[i]<<endl;
            }
        }
    
        stop = chrono::high_resolution_clock::now();
        int64_t searchDelay = chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        cout<<"Searched "<<totalSearch<<" elements in "<<searchDelay<<" microSeconds."<<endl;
    }
    
    //Range Scan in the PMA
    if(rangeLength>0){
        start = chrono::high_resolution_clock::now();
        for(int i=0; i<rangeIteration; i++){
            type_t startRange = totalInsert == rangeLength? 0 : rand()%(totalInsert - rangeLength);
            type_t sum_key, sum_value;
            tie(sum_key, sum_value) = pma.range_sum(startRange, startRange+rangeLength);
            if(sum_key*10 != sum_value){
                cout<<"Error in range scan!"<<endl;
                exit(0);
            }
        }
        stop = chrono::high_resolution_clock::now();
        int64_t scanDelay = chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        cout<<"Scanned elements with range "<<rangeLength<<" total "<<rangeIteration<<" times in " <<scanDelay<<" microSeconds."<<endl;
    }
    
    //records = (int64_t *)malloc(totalDelete * sizeof(int64_t));
    if(totalDelete > 0){
        start = chrono::high_resolution_clock::now();
        for(type_t i=0; i<totalDelete; i++){
            //records[i] = numbers(rng);
            //cin>>records[i];
            pma.remove(data[i]);
        }
        stop = chrono::high_resolution_clock::now();
        int64_t deleteDelay = chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        cout<<"Deleted "<<totalDelete<<" elements in "<<deleteDelay<<" microSeconds."<<endl;
    }
    return 0;
}