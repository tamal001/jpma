#ifndef JPMA_HPP_
#define JPMA_HPP_

#include <vector>

#include "defines.hpp"
using namespace std;

class PMA {
private:
    vector<type_t *> key_chunks;
    vector<type_t *> value_chunks;
    vector<type_t> smallest;
    vector<type_t> lastElementPos;
    vector<int> cardinality;
    int totalSegments = 0;
    u_char NonZeroEntries[JacobsonIndexCount][JacobsonIndexSize+1];
    vector<vector<u_short>> bitmap;
    type_t lastValidPos;             //Last accessible slot in each segment
    
public:
    enum redistribution{
        INSERT = 1,
        DELETE = 2,
        REDESIGN = 3
    };
    PMA();
    ~PMA();

    //Library functions
    bool insert(type_t key, type_t value, int insertCount);
    bool remove(type_t key);
    bool lookup(type_t key);
    type_t range_sum(type_t startKey, type_t endKey);

    //Support functions
    void preCalculateJacobson();
    bool insertInPosition(type_t position, int targetSegment, type_t key, type_t value);
    bool deleteInPosition(type_t position, int targetSegment, type_t key);
    void deleteSegment(int targetSegment);
    type_t findLocation(type_t key, int targetSegment);
    void redistribute(int targetSegment, redistribution type);
    void redistributeWithPrev(int targetSegment);
    void redistributeWithNext(int targetSegment);
    void redistributeWithDividing(int targetSegment);
    void swapElements(type_t targetSegment, type_t position, type_t adjust);

    //Testing functions
    void printAllElements(void);
    void printSegElements(int targetSegment);
};
#endif