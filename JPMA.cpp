#include <iostream>
#include <vector>
#include <random>

#include "defines.hpp"
#include "JPMA.hpp"
#include "diymalloc.h"
using namespace std;

PMA::PMA(){
    smallest.push_back(1);                                   //First segment has smallest element 1
    cardinality.push_back(0);                                //First segment contains 0 elements
    totalSegments = 1;                                       //One segment created at the start
    lastElementPos.push_back(0);                             //Position of last element in the segment
    lastValidPos = (SEGMENT_SIZE/sizeof(type_t)) - 1;
    cout << "Last valid position: "<<lastValidPos << endl;

    //Create space for first segment
    type_t *starting_key_chunk, *starting_value_chunk;
    starting_key_chunk = (type_t *) malloc (SEGMENT_SIZE);
    starting_value_chunk = (type_t *) malloc (SEGMENT_SIZE);
    key_chunks.push_back(starting_key_chunk);
    value_chunks.push_back(starting_value_chunk);
    vector <u_short> blocks;
    for(int i = 0; i<BLOCKS_IN_SEGMENT; i++){
        blocks.push_back(0);
    }
    bitmap.push_back(blocks);

    //Create jacobson Index
    preCalculateJacobson();
}

PMA::~PMA(){
    int targetSegment;
    for(targetSegment = 0; targetSegment<totalSegments; targetSegment++){
        delete[] key_chunks[targetSegment];
        delete[] value_chunks[targetSegment];
    }
}

void PMA::preCalculateJacobson(){
    for(int i = 0; i<JacobsonIndexCount; i++){
        int j, k = 1, idx = 0;
        u_char count = 0;
        for(j = 0; j<JacobsonIndexSize; j++){
            if(i & k) count++;
            k = k<<1;
        }
        NonZeroEntries[i][idx++] = count;
        k = 1;
        for(j = 0; j<JacobsonIndexSize; j++){
            if(i & k) NonZeroEntries[i][idx++] = j;
            k = k<<1;
        }
    }
}
bool PMA::insert(type_t key, type_t value){
    //Find the location using Binary Search.
    int targetSegment;
    for(targetSegment = 0; targetSegment < totalSegments; targetSegment++){
        if(smallest[targetSegment]>key) break;
    }
    targetSegment--;

    type_t position = findLocation(key, targetSegment);
    //cout <<" Found position: " << position<< endl;
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t foundKey = *(segmentOffset + position);
    if(foundKey == key) return false;
    else if(foundKey == 0) return insertInPosition(position, targetSegment, key, value);

    int probCount;
    if (position >= lastElementPos[targetSegment]){ //Got the last inserted postion of the current segment
        if(position >= lastValidPos ){ //Got out of the current segment. Traverse backward for vacant space
            probCount = 1;
            type_t *insertPos = segmentOffset + position - probCount;
            while(*insertPos != 0){
                insertPos--;
                probCount++;
                if(probCount > ProbLimit){
                    redistribute(targetSegment);
                    return insert(key, value);
                }
            }
            insertInPosition(position-probCount, targetSegment, key, value);

            while(*insertPos > *(insertPos+1)){
                swapElements(targetSegment, position-probCount, 1);
                probCount--;
                if(probCount == 0) break;
                insertPos++;
            }
            return true;
        }
        else{ //Have some space left in the segment. Go forward in the space max 3 slots
            type_t adjust = min(lastValidPos - lastElementPos[targetSegment], (type_t) MaxGap);
            adjust = min(adjust, abs(*(key_chunks[targetSegment]+lastElementPos[targetSegment])-key));
            if(key > foundKey){
                return insertInPosition(position+adjust, targetSegment, key, value);
            }
            else{
                insertInPosition(position+adjust, targetSegment, key, value);
                swapElements(targetSegment, position, adjust);
                return true;
            }
        }
    }

    //Insert among other inserted elements
    probCount = 1;
    type_t *insertPos = segmentOffset + position + probCount;
    while(*insertPos != 0){
        probCount++;
        insertPos++;
        if(probCount > ProbLimit) {
            redistribute(targetSegment);
            return insert(key, value);
        }
    }
    insertInPosition(position+probCount, targetSegment, key, value);
    while(*insertPos < *(insertPos-1)){
        probCount--;
        swapElements(targetSegment, position+probCount, 1);
        if(probCount == 0) break;
        insertPos--;
    }
    //Adjust position not to be adjacent. Keep some space

    return true;
}

void PMA::swapElements(type_t targetSegment, type_t position, type_t adjust){
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t holdKey = *(segmentOffset + position);
    *(segmentOffset + position) = *(segmentOffset + position + adjust);
    *(segmentOffset + position + adjust) = holdKey;

    segmentOffset = value_chunks[targetSegment];
    type_t holdValue = *(segmentOffset + position);
    *(segmentOffset + position) = *(segmentOffset + position + adjust);
    *(segmentOffset + position + adjust) = holdValue;
}

bool PMA::insertInPosition(type_t position, int targetSegment, type_t key, type_t value){
    //Store key, value and update bitmap, cardinality and last index
    type_t * segmentOffset = key_chunks[targetSegment];
    *(segmentOffset + position) = key;
    segmentOffset = value_chunks[targetSegment];
    *(segmentOffset + position) = value;
    int blockPosition = position/JacobsonIndexSize;
    int bitPosition = position % JacobsonIndexSize;
    u_short mask =  1 << bitPosition;
    bitmap[targetSegment][blockPosition] |= mask;
    cardinality[targetSegment]++;
    if(lastElementPos[targetSegment] < position) lastElementPos[targetSegment] = position;
    if(cardinality[targetSegment] > MaxThreshold) redistribute(targetSegment);
    return true;
}

bool PMA::remove(type_t key){
    int targetSegment;
    for(targetSegment = 0; targetSegment<totalSegments; targetSegment++){
        if(smallest[targetSegment]>key) break;
    }
    targetSegment--;

    type_t position = findLocation(key, targetSegment);
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t foundKey = *(segmentOffset + position);
    if(foundKey != key) return false;    
    return deleteInPosition(position, targetSegment, key);
}

bool PMA::deleteInPosition(type_t position, int targetSegment, type_t key){
    type_t * segmentOffset = key_chunks[targetSegment];
    *(segmentOffset + position) = 0;
    segmentOffset = value_chunks[targetSegment];
    *(segmentOffset + position) = 0;
    int blockPosition = position/JacobsonIndexSize;
    int bitPosition = position % JacobsonIndexSize;
    u_short mask =  1 << bitPosition;
    bitmap[targetSegment][blockPosition] &= (~mask);
    cardinality[targetSegment]--;
    if(lastElementPos[targetSegment] == position){
        position--;
        segmentOffset = key_chunks[targetSegment];
        while(*(segmentOffset + position) == 0){
            position--;
        }
        lastElementPos[targetSegment] = position;
    }
    if(cardinality[targetSegment] > MinThreshold) redistribute(targetSegment);
    return true;
}

void PMA::redistribute(int targetSegment){
    cout << "CALLED FOR REDISTRIBUTION. Cardinality "<<cardinality[targetSegment] <<endl;
    if(cardinality[targetSegment] < MinThreshold){ //Try to delete this segment
        if(targetSegment > 0 && cardinality[targetSegment] + cardinality[targetSegment-1] < MaxThreshold){
            redistributeWithPrev(targetSegment);
        }else if(targetSegment<totalSegments-1 && cardinality[targetSegment] + cardinality[targetSegment+1] < MaxThreshold){
            redistributeWithNext(targetSegment);
        }
    }else if(cardinality[targetSegment] > MaxThreshold){ //Try to add new segment
        redistributeWithDividing(targetSegment);
    }else{ //Redistribute the contents of this Segment
        type_t ar[lastValidPos+1], br[lastValidPos+1];
        int elementCount = 0;
        type_t * moveKeyOffset = key_chunks[targetSegment];
        type_t * moveValOffset = value_chunks[targetSegment];

        for(int i = 0; i <= lastElementPos[targetSegment]; i++){
            if(*(moveKeyOffset + i) != 0){
                ar[elementCount] = *(moveKeyOffset + i);
                br[elementCount] = *(moveValOffset + i);
                *(moveKeyOffset + i) = 0;
                *(moveValOffset + i) = 0;
                elementCount++;
            }
        }
        cout << "copied "<<elementCount<<" elements." <<endl;
        for(int i = 0; i < BLOCKS_IN_SEGMENT; i++){     //Clear existing bitmap
            bitmap[targetSegment][i] = 0;
        }
        cout <<"cleared bitmap blocks"<<endl;
        default_random_engine generator;
        uniform_real_distribution<double> distribution(0.0,1.0);
        double density = cardinality[targetSegment] / (lastValidPos + 1.0);
        type_t i, j=0, gap = 0;
        for(i = 0; i <= lastValidPos; i++){
            if(elementCount == lastValidPos - i){
                //direct insert
                *(moveKeyOffset + i) = ar[j];
                *(moveValOffset + i) = br[j];
                int blockPosition = i / JacobsonIndexSize;
                int bitPosition = i % JacobsonIndexSize;
                u_short mask =  1 << bitPosition;
                bitmap[targetSegment][blockPosition] |= mask;
                j++; elementCount--;
                if(elementCount == 0) {
                    lastElementPos[targetSegment] = i;
                    break;
                }
            }else{
                double d = distribution(generator);
                if(d < density || gap == MaxGap){
                    *(moveKeyOffset + i) = ar[j];
                    *(moveValOffset + i) = br[j];
                    int blockPosition = i / JacobsonIndexSize;
                    int bitPosition = i % JacobsonIndexSize;
                    u_short mask =  1 << bitPosition;
                    bitmap[targetSegment][blockPosition] |= mask;
                    j++; elementCount--; gap = 0;
                    if(elementCount == 0){
                        lastElementPos[targetSegment] = i;
                        break;
                    }
                }else gap++;
            }
        }
        for(i = 0; i <= lastValidPos; i++){
            if(*(moveKeyOffset + i) != 0){
                smallest[targetSegment] = i;
                break;
            }
        }
    }
}

void PMA::redistributeWithPrev(int targetSegment){
    if(lastElementPos[targetSegment-1] + lastElementPos[targetSegment] < lastValidPos){
        //Direct copy with key and spaces
        type_t j = 0;
        type_t * destKeyOffset = key_chunks[targetSegment-1];
        type_t * destValOffset = value_chunks[targetSegment-1];
        type_t * moveKeyOffset = key_chunks[targetSegment];
        type_t * moveValOffset = value_chunks[targetSegment];
        type_t key, value;
        for(type_t i = lastElementPos[targetSegment-1]+1; j <= lastElementPos[targetSegment]; i++,j++){
            key = *(moveKeyOffset + j);
            value = *(moveValOffset + j);
            if(key != 0){
                *(destKeyOffset + i) = key;
                *(destValOffset + i) = value;
                int blockPosition = i / JacobsonIndexSize;
                int bitPosition = i % JacobsonIndexSize;
                u_short mask =  1 << bitPosition;
                bitmap[targetSegment-1][blockPosition] |= mask;
            }
        }
        cardinality[targetSegment-1] += cardinality[targetSegment];
        lastElementPos[targetSegment - 1] += lastElementPos[targetSegment] + 1;
        deleteSegment(targetSegment);
    }else if(lastElementPos[targetSegment-1] + cardinality[targetSegment] < lastValidPos){
        //copy with last part redistribute
        type_t slots = lastValidPos - lastElementPos[targetSegment - 1];
        default_random_engine generator;
        uniform_real_distribution<double> distribution(0.0,1.0);
        double density = cardinality[targetSegment] / (double) slots;
        type_t i, j = 0, elements = cardinality[targetSegment];
        type_t * destKeyOffset = key_chunks[targetSegment-1];
        type_t * destValOffset = value_chunks[targetSegment-1];
        type_t * moveKeyOffset = key_chunks[targetSegment];
        type_t * moveValOffset = value_chunks[targetSegment];

        for(i = lastElementPos[targetSegment-1]+1; i <= lastValidPos; i++){
            if(elements == lastValidPos - i){
                while(*(moveKeyOffset + j) == 0 ){
                    j++;
                }
                *(destKeyOffset + i) = *(moveKeyOffset + j);
                *(destValOffset + i) = *(moveValOffset + j);
                int blockPosition = i / JacobsonIndexSize;
                int bitPosition = i % JacobsonIndexSize;
                u_short mask =  1 << bitPosition;
                bitmap[targetSegment-1][blockPosition] |= mask;
                elements--;
                if(elements == 0){
                    lastElementPos[targetSegment-1] = i;
                    break;
                }
            }else{
                double d = distribution(generator);
                if(d < density){
                    while(*(moveKeyOffset + j) == 0 ){
                        j++;
                    }
                    *(destKeyOffset + i) = *(moveKeyOffset + j);
                    *(destValOffset + i) = *(moveValOffset + j);
                    int blockPosition = i / JacobsonIndexSize;
                    int bitPosition = i % JacobsonIndexSize;
                    u_short mask =  1 << bitPosition;
                    bitmap[targetSegment-1][blockPosition] |= mask;
                    elements--;
                    if(elements == 0){
                        lastElementPos[targetSegment-1] = i;
                        break;
                    }
                }
            }
        }
        cardinality[targetSegment-1] += cardinality[targetSegment];
        //lastElementPos[targetSegment - 1] = i;
        deleteSegment(targetSegment);
    }else{
        type_t *new_key_chunk, *new_value_chunk;
        new_key_chunk = (type_t *) malloc (SEGMENT_SIZE);
        new_value_chunk = (type_t *) malloc (SEGMENT_SIZE);
        type_t * destKeyOffset = new_key_chunk;
        type_t * destValOffset = new_value_chunk;
        type_t * moveKeyOffset = key_chunks[targetSegment-1];
        type_t * moveValOffset = value_chunks[targetSegment-1];
        
        default_random_engine generator;
        uniform_real_distribution<double> distribution(0.0,1.0);
        double density = (cardinality[targetSegment] + cardinality[targetSegment-1]) / (lastValidPos + 1.0);
        type_t i, j = 0, elements = cardinality[targetSegment-1] + cardinality[targetSegment];
        bool whichSegment = true;             //Change the incoming segment only once
        
        //clear previous bitmap
        vector<u_short> blocks;
        for(i = 0; i < BLOCKS_IN_SEGMENT; i++){
            blocks.push_back(0);
        }
        
        for(i = 0; i <= lastValidPos && j <= lastElementPos[targetSegment-1]; i++){
            if(elements == lastValidPos - i){
                while(*(moveKeyOffset + j) == 0 ){
                    j++;
                }
                *(destKeyOffset + i) = *(moveKeyOffset + j);
                *(destValOffset + i) = *(moveValOffset + j);
                int blockPosition = i / JacobsonIndexSize;
                int bitPosition = i % JacobsonIndexSize;
                u_short mask =  1 << bitPosition;
                blocks[blockPosition] |= mask;
                elements--;
            }else{
                double d = distribution(generator);
                if(d < density){
                    while(*(moveKeyOffset + j) == 0 ){
                        j++;
                    }
                    *(destKeyOffset + i) = *(moveKeyOffset + j);
                    *(destValOffset + i) = *(moveValOffset + j);
                    int blockPosition = i / JacobsonIndexSize;
                    int bitPosition = i % JacobsonIndexSize;
                    u_short mask =  1 << bitPosition;
                    blocks[blockPosition] |= mask;
                    elements--;
                }
            }
            if(whichSegment && j == lastElementPos[targetSegment-1]){
                whichSegment = false;
                moveKeyOffset = key_chunks[targetSegment];
                moveValOffset = value_chunks[targetSegment];
                j = 0;
            }
        }
        smallest.erase(smallest.begin() + targetSegment); //Keep smallest element at targetsegment -1
        cardinality[targetSegment - 1] += cardinality [targetSegment];
        totalSegments--;
        lastElementPos[targetSegment - 1] = i;

        //Remove previous key and value chunks. Insert the new one.
        delete[] key_chunks[targetSegment - 1], key_chunks[targetSegment];
        key_chunks.erase(key_chunks.begin() + targetSegment-1, key_chunks.begin() + targetSegment);
        key_chunks.emplace(key_chunks.begin() + targetSegment - 1, new_key_chunk);

        delete[] value_chunks[targetSegment - 1], value_chunks[targetSegment];
        value_chunks.erase(value_chunks.begin() + targetSegment - 1, value_chunks.begin() + targetSegment);
        value_chunks.emplace(value_chunks.begin() + targetSegment - 1, new_value_chunk);

        //Remove previous bitmaps. Insert the new Bitmap.
        bitmap.erase(bitmap.begin() + targetSegment - 1, bitmap.begin() + targetSegment);
        bitmap.emplace(bitmap.begin() + targetSegment - 1, blocks);
    }
}

void PMA::redistributeWithNext(int targetSegment){
    redistributeWithPrev(targetSegment + 1);
}

void PMA::redistributeWithDividing(int targetSegment){
    type_t halfElement = cardinality[targetSegment]/2;
    type_t *new_key_chunk, *new_value_chunk;
    new_key_chunk = (type_t *) malloc (SEGMENT_SIZE);
    new_value_chunk = (type_t *) malloc (SEGMENT_SIZE);

    type_t * moveKeyOffset = key_chunks[targetSegment];
    type_t * moveValOffset = value_chunks[targetSegment];
    type_t * destKeyOffset = new_key_chunk;
    type_t * destValOffset = new_value_chunk;

    type_t i, j, elementCount = 0, lastElement = lastElementPos[targetSegment];
    for(i = 0; i < lastElement; i++){
        if(*(moveKeyOffset + i) != 0) elementCount++;
        if(elementCount == halfElement) break;
    }
    lastElementPos[targetSegment] = i;

    vector<u_short> blocks;
    for(i = 0; i < BLOCKS_IN_SEGMENT; i++){
        blocks.push_back(0);
    }
    //Copy the elements to the new segment and erase from previous segment
    for(++i, j = 0; i <= lastElement; i++, j++){
        *(destKeyOffset + j) = *(moveKeyOffset + i);
        *(destValOffset + j) = *(moveKeyOffset + i);
        if(*(moveKeyOffset + i) != 0){
            *(moveKeyOffset + i) = 0;
            *(moveValOffset + i) = 0;
            int blockPosition = i / JacobsonIndexSize;
            int bitPosition = i % JacobsonIndexSize;
            u_short mask =  1 << bitPosition;
            bitmap[targetSegment][blockPosition] &= (~mask);

            blockPosition = j / JacobsonIndexSize;
            bitPosition = i % JacobsonIndexSize;
            mask = i << bitPosition;
            blocks[blockPosition] |= mask;
        }
    }
    key_chunks.emplace(key_chunks.begin() + targetSegment + 1, new_key_chunk);
    value_chunks.emplace(value_chunks.begin() + targetSegment + 1, new_value_chunk);
    lastElementPos.emplace(lastElementPos.begin() + targetSegment + 1, j-1);
    cardinality.emplace(cardinality.begin() + targetSegment + 1, cardinality[targetSegment] - halfElement);
    cardinality[targetSegment] = halfElement;
    bitmap.emplace(bitmap.begin() + targetSegment + 1, blocks);
    totalSegments++;
    for(i = 0; i < cardinality[targetSegment + 1]; i++){
        if(*(destKeyOffset + i) != 0){
            smallest.emplace(smallest.begin() + targetSegment + 1, *(destKeyOffset + i));
            break;
        }
    }
}

void PMA::deleteSegment(int targetSegment){
    delete[] key_chunks[targetSegment];
    delete[] value_chunks[targetSegment];
    key_chunks.erase(key_chunks.begin() + targetSegment);
    value_chunks.erase(value_chunks.begin() + targetSegment);
    smallest.erase(smallest.begin() + targetSegment);
    lastElementPos.erase(lastElementPos.begin() + targetSegment);
    cardinality.erase(cardinality.begin() + targetSegment);
    bitmap.erase(bitmap.begin() + targetSegment);
    totalSegments--;
}

type_t PMA::range_sum(type_t startKey, type_t endKey){
    int targetSegment;
    for(targetSegment = 0; targetSegment<totalSegments; targetSegment++){
        if(smallest[targetSegment]>startKey) break;
    }
    targetSegment--;

    type_t position = findLocation(startKey, targetSegment);
    type_t sum_key = 0, sum_value = 0;
    type_t blockNo = position/JacobsonIndexSize;
    u_char * ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
    type_t * segmentKeyOffset = key_chunks[targetSegment];
    type_t * segmentValOffset = value_chunks[targetSegment];
    type_t pbase = blockNo * JacobsonIndexSize;

    bool flag = false;
    cout << "start Position "<<position<<", blockNo "<<blockNo<<", ar[0] "<<ar[0]<<endl;
    for(int offset = 1; offset <= ar[0] ; offset++){
        //base = pBase + ar[offset];
        type_t key = *(segmentKeyOffset+pbase+ar[offset]);
        if(flag && key > endKey) flag = false;
        else if(!flag && key >= startKey) flag = true;
        if(flag) {
            sum_key += key;
            sum_value += *(segmentValOffset+pbase+ar[offset]);
        }
    }

    type_t * key_pos, *value_pos;
    int offset;
    while(true){
        pbase += JacobsonIndexSize;
        blockNo++;
        if(blockNo == BLOCKS_IN_SEGMENT){
            blockNo = 0;
            pbase = 0;
            targetSegment++;
            if(targetSegment == totalSegments) return sum_key;
            segmentKeyOffset = key_chunks[targetSegment];
            segmentValOffset = value_chunks[targetSegment];
        }
        ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
        key_pos = segmentKeyOffset + pbase;
        value_pos = segmentValOffset + pbase;
        offset = ar[0];
        switch (offset){
            case 1:
                sum_key += *(key_pos + ar[1]);
                sum_value += *(value_pos + ar[1]);
                break;
            case 2:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                break;
            case 3:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                break;
            case 4:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                break;
            case 5:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                break;
            case 6:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                break;
            case 7:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_key += *(key_pos + ar[7]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                sum_value += *(value_pos + ar[7]);
                break;
            case 8:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_key += *(key_pos + ar[7]);
                sum_key += *(key_pos + ar[8]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                sum_value += *(value_pos + ar[7]);
                sum_value += *(value_pos + ar[8]);
                break;
            case 9:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_key += *(key_pos + ar[7]);
                sum_key += *(key_pos + ar[8]);
                sum_key += *(key_pos + ar[9]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                sum_value += *(value_pos + ar[7]);
                sum_value += *(value_pos + ar[8]);
                sum_value += *(value_pos + ar[9]);
                break;
            case 10:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_key += *(key_pos + ar[7]);
                sum_key += *(key_pos + ar[8]);
                sum_key += *(key_pos + ar[9]);
                sum_key += *(key_pos + ar[10]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                sum_value += *(value_pos + ar[7]);
                sum_value += *(value_pos + ar[8]);
                sum_value += *(value_pos + ar[9]);
                sum_value += *(value_pos + ar[10]);
                break;
            case 11:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_key += *(key_pos + ar[7]);
                sum_key += *(key_pos + ar[8]);
                sum_key += *(key_pos + ar[9]);
                sum_key += *(key_pos + ar[10]);
                sum_key += *(key_pos + ar[11]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                sum_value += *(value_pos + ar[7]);
                sum_value += *(value_pos + ar[8]);
                sum_value += *(value_pos + ar[9]);
                sum_value += *(value_pos + ar[10]);
                sum_value += *(value_pos + ar[11]);
                break;
            case 12:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_key += *(key_pos + ar[7]);
                sum_key += *(key_pos + ar[8]);
                sum_key += *(key_pos + ar[9]);
                sum_key += *(key_pos + ar[10]);
                sum_key += *(key_pos + ar[11]);
                sum_key += *(key_pos + ar[12]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                sum_value += *(value_pos + ar[7]);
                sum_value += *(value_pos + ar[8]);
                sum_value += *(value_pos + ar[9]);
                sum_value += *(value_pos + ar[10]);
                sum_value += *(value_pos + ar[11]);
                sum_value += *(value_pos + ar[12]);
                break;
            case 13:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_key += *(key_pos + ar[7]);
                sum_key += *(key_pos + ar[8]);
                sum_key += *(key_pos + ar[9]);
                sum_key += *(key_pos + ar[10]);
                sum_key += *(key_pos + ar[11]);
                sum_key += *(key_pos + ar[12]);
                sum_key += *(key_pos + ar[13]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                sum_value += *(value_pos + ar[7]);
                sum_value += *(value_pos + ar[8]);
                sum_value += *(value_pos + ar[9]);
                sum_value += *(value_pos + ar[10]);
                sum_value += *(value_pos + ar[11]);
                sum_value += *(value_pos + ar[12]);
                sum_value += *(value_pos + ar[13]);
                break;
            case 14:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_key += *(key_pos + ar[7]);
                sum_key += *(key_pos + ar[8]);
                sum_key += *(key_pos + ar[9]);
                sum_key += *(key_pos + ar[10]);
                sum_key += *(key_pos + ar[11]);
                sum_key += *(key_pos + ar[12]);
                sum_key += *(key_pos + ar[13]);
                sum_key += *(key_pos + ar[14]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                sum_value += *(value_pos + ar[7]);
                sum_value += *(value_pos + ar[8]);
                sum_value += *(value_pos + ar[9]);
                sum_value += *(value_pos + ar[10]);
                sum_value += *(value_pos + ar[11]);
                sum_value += *(value_pos + ar[12]);
                sum_value += *(value_pos + ar[13]);
                sum_value += *(value_pos + ar[14]);
                break;
            case 15:
                sum_key += *(key_pos + ar[1]);
                sum_key += *(key_pos + ar[2]);
                sum_key += *(key_pos + ar[3]);
                sum_key += *(key_pos + ar[4]);
                sum_key += *(key_pos + ar[5]);
                sum_key += *(key_pos + ar[6]);
                sum_key += *(key_pos + ar[7]);
                sum_key += *(key_pos + ar[8]);
                sum_key += *(key_pos + ar[9]);
                sum_key += *(key_pos + ar[10]);
                sum_key += *(key_pos + ar[11]);
                sum_key += *(key_pos + ar[12]);
                sum_key += *(key_pos + ar[13]);
                sum_key += *(key_pos + ar[14]);
                sum_key += *(key_pos + ar[15]);
                sum_value += *(value_pos + ar[1]);
                sum_value += *(value_pos + ar[2]);
                sum_value += *(value_pos + ar[3]);
                sum_value += *(value_pos + ar[4]);
                sum_value += *(value_pos + ar[5]);
                sum_value += *(value_pos + ar[6]);
                sum_value += *(value_pos + ar[7]);
                sum_value += *(value_pos + ar[8]);
                sum_value += *(value_pos + ar[9]);
                sum_value += *(value_pos + ar[10]);
                sum_value += *(value_pos + ar[11]);
                sum_value += *(value_pos + ar[12]);
                sum_value += *(value_pos + ar[13]);
                sum_value += *(value_pos + ar[14]);
                sum_value += *(value_pos + ar[15]);
                break;
            case 16:
                sum_key += *key_pos;
                sum_key += *(key_pos + 1);
                sum_key += *(key_pos + 2);
                sum_key += *(key_pos + 3);
                sum_key += *(key_pos + 4);
                sum_key += *(key_pos + 5);
                sum_key += *(key_pos + 6);
                sum_key += *(key_pos + 7);
                sum_key += *(key_pos + 8);
                sum_key += *(key_pos + 9);
                sum_key += *(key_pos + 10);
                sum_key += *(key_pos + 11);
                sum_key += *(key_pos + 12);
                sum_key += *(key_pos + 13);
                sum_key += *(key_pos + 14);
                sum_key += *(key_pos + 15);
                sum_value += *value_pos;
                sum_value += *(value_pos + 1);
                sum_value += *(value_pos + 2);
                sum_value += *(value_pos + 3);
                sum_value += *(value_pos + 4);
                sum_value += *(value_pos + 5);
                sum_value += *(value_pos + 6);
                sum_value += *(value_pos + 7);
                sum_value += *(value_pos + 8);
                sum_value += *(value_pos + 9);
                sum_value += *(value_pos + 10);
                sum_value += *(value_pos + 11);
                sum_value += *(value_pos + 12);
                sum_value += *(value_pos + 13);
                sum_value += *(value_pos + 14);
                sum_value += *(value_pos + 15);
        }
        if(*(key_pos + ar[offset]) > endKey){
            while(offset > 0 && *(key_pos + ar[offset]) > endKey){
                sum_key -= *(key_pos + ar[offset]);
                sum_value -= *(value_pos + ar[offset]);
                offset--;
            }
		    return sum_key;
	    }
    }
    return sum_key;
}

bool PMA::lookup(type_t key){
    int targetSegment;
    for(targetSegment = 0; targetSegment<totalSegments; targetSegment++){
        if(smallest[targetSegment]>key) break;
    }
    targetSegment--;

    type_t position = findLocation(key, targetSegment);
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t foundKey = *(segmentOffset + position);
    if(foundKey == key) return true;  
    return false;
}

type_t PMA::findLocation(type_t key, int targetSegment){
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t start = 0;
    type_t end = lastElementPos[targetSegment];
    type_t data, mid = 0;
    while(start <= end){
        mid = (start + end) / 2;
        //if(key[mid] == 0){
        data = *(segmentOffset+mid);
        if(data == 0){
            type_t changedMid = mid, offset = -1;
            while(changedMid > start){
                changedMid += offset;
                if(*(segmentOffset + changedMid)!=0) break;
            }
            if(changedMid == start){
                if(*(segmentOffset + start) >= key) {mid = start; break;}
                changedMid = mid;
                offset = 1;
                while(changedMid < end){
                    changedMid += offset;
                    if(*(segmentOffset + changedMid)!=0) break;
                }
                if(changedMid == end){
                    mid = end; 
                    data = *(segmentOffset+mid);
                    break;
                }else {
                    mid = changedMid;
                    data = *(segmentOffset+mid);
                }
            }else {
                mid = changedMid;
                data = *(segmentOffset+mid);
            }
        }
        if(data == key) break;
        else if(data < key) start = mid + 1;
        else end = mid - 1;
    }
    return mid;
}

void PMA::printAllElements(){
    for(type_t i = 0; i<totalSegments; i++){
        type_t * key = key_chunks[i];
        for(type_t j = 0; j<=lastElementPos[i]; j++){
            cout << *(key+j)<< " ";
        }
        cout<<"last offset: "<<lastElementPos[i] <<" Total Segment: "<<totalSegments<< endl;
    }

}
