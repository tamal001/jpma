#include <iostream>
#include <vector>
#include <random>
#include <string.h>
#include <sys/mman.h>
#include <tuple>
#include <cassert>
#include <cmath>

#include "defines.hpp"
#include "JPMA_BT.hpp"
#include "diymalloc.h"

using namespace std;

//uint64_t totalRebalance = 0;
//uint64_t totalShiftingInsert = 0;
//uint64_t totalShiftingReb = 0;
//uint64_t rewiring_count = 0;
//uint64_t totalInserts = 0;
//int maxheight = 0;

PMA::PMA(type_t totalInsert){
    elementsInSegment = SEGMENT_SIZE/sizeof(type_t);
    int estSegment = (int)( totalInsert/(elementsInSegment * 0.75));
    smallest.reserve(estSegment);
    cardinality.reserve(estSegment);
    lastElementPos.reserve(estSegment);
    key_chunks.reserve(estSegment);
    value_chunks.reserve(estSegment);
    freeSegID.reserve(estSegment);
    segCount = -1;
    smallest.push_back(1);                                   //First segment has smallest element 1
    cardinality.push_back(0);                                //First segment contains 0 elements
    totalSegments = 1;                                       //One segment deployed at the start
    lastElementPos.push_back(0);                             //Position of last element in the segment
    lastValidPos = elementsInSegment - 1;
    blocksInSegment = elementsInSegment / JacobsonIndexSize;
    freeSegmentCount = 0;
    int CurSegemnt = getSegment();
    
    tree = new BPlusTree(this);
    tree->insertInTree(CurSegemnt, 0, this); //(segment no, dummy key, current JPMA object)

    //Create jacobson Index
    preCalculateJacobson();
}

PMA::~PMA(){
    u_int size = cleanSegments.size();
    for(u_int i = 0; i<size; i++){
        type_t * temp = cleanSegments.back();
        cleanSegments.pop_back();
        delete temp;
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

int PMA::searchSegment(type_t key){
    return tree->searchSegment(key);
}

int PMA::getSegment(){
    //Get a segment from the pool of free segment IDs in freeSegID vector. If the pool is empty, create a bunch of free segments
    if(UNLIKELY(freeSegmentCount < 1)){
        type_t *new_key_chunk, *new_value_chunk;
        if(Allocation_type == 1){
            new_key_chunk = (type_t *) mmap(ADDR, CHUNK_SIZE, PROTECTION, FLAGS, -1, 0);
            if(new_key_chunk == MAP_FAILED){ 
                cout<<"Cannot allocate the virtual memory: " << CHUNK_SIZE << " bytes. mmap error: " << strerror(errno) << "(" << errno << ")"; 
                exit(0);
            }
            new_value_chunk = (type_t *) mmap(ADDR, CHUNK_SIZE, PROTECTION, FLAGS, -1, 0);    
            if(new_value_chunk == MAP_FAILED){ 
                cout<<"Cannot allocate the virtual memory: " << CHUNK_SIZE*16 << " bytes. mmap error: " << strerror(errno) << "(" << errno << ")"; 
                exit(0);
            }
        }
        else{
            new_key_chunk = (type_t *) malloc (CHUNK_SIZE);
            new_value_chunk = (type_t *) malloc (CHUNK_SIZE);    
        }
        cleanSegments.push_back(new_key_chunk);
        cleanSegments.push_back(new_value_chunk);

        freeSegmentCount = CHUNK_SIZE / SEGMENT_SIZE;
        vector<u_short> musks (blocksInSegment, 0);
        segCount += freeSegmentCount;
        for(int i = 0; i < freeSegmentCount; i++){
            freeSegID.push_back(segCount-i);
            key_chunks.push_back(new_key_chunk + i * elementsInSegment);
            value_chunks.push_back(new_value_chunk + i * elementsInSegment);
            smallest.push_back(0);
            lastElementPos.push_back(0);
            cardinality.push_back(0);
            bitmap.push_back(musks);
        }
    }
    freeSegmentCount--;
    int curSegNo = freeSegID.back();
    freeSegID.pop_back();

    smallest[curSegNo] = 0; lastElementPos[curSegNo] = 0; cardinality[curSegNo] = 0;
    for(int i=0; i<blocksInSegment; i++) bitmap[curSegNo][i] = 0;
    memset(key_chunks[curSegNo], 0, sizeof(int64_t)*elementsInSegment);

    return curSegNo;
}

bool PMA::insert(type_t key, type_t value){
    //Find the location using Binary Search.
    int targetSegment = tree->searchSegment(key);
    type_t position = findLocation1(key, targetSegment);

    int blockNo = position/JacobsonIndexSize;
    int bitPosition = position % JacobsonIndexSize;
    u_short mask =  1 << bitPosition;

    if((bitmap[targetSegment][blockNo] & mask) == 0){
        insertInPosition(position, targetSegment, key, value);
        if(cardinality[targetSegment] > tree->maxLevel[0]) redistributeInsert(targetSegment, smallest[targetSegment]);
        return true;
    }

    type_t * segmentOffset = key_chunks[targetSegment];
    type_t foundKey = *(segmentOffset + position);
    if(foundKey == key) return false;

    //check if need traversing from backside
    if(position >= lastElementPos[targetSegment]){
        if(!insertAfterLast(position, key, value, targetSegment, foundKey)) return false;
        if(cardinality[targetSegment] > tree->maxLevel[0]) redistributeInsert(targetSegment, smallest[targetSegment]);
        return true;
    }


    //Insert among other inserted elements 
    u_char * ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
    //type_t * segmentKeyOffset = key_chunks[targetSegment];
    int pBase = blockNo * JacobsonIndexSize;
    int insertPos = pBase;
    while(true){
        if(UNLIKELY(ar[0] == 0)){//Got a complete empty block. select the first cell. This can not be the starting block
            insertPos = pBase;
            break;
        }
        else if(ar[0] == JacobsonIndexSize){
            blockNo++;
            if(blockNo == blocksInSegment){
                if(!backSearchInsert(position, key, value, targetSegment, lastValidPos+position)) return false;
                if(cardinality[targetSegment] > tree->maxLevel[0]) redistributeInsert(targetSegment, smallest[targetSegment]); 
                return true;
            }
            ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
            pBase += JacobsonIndexSize;
            insertPos = pBase;
        }else{
            //Before first element of the block
            if(ar[1] != 0 && position <= pBase){
                insertPos = pBase;
                break;
            } 
            //In between other inserted elements in the block
            bool got = false;
            for(int j = 2; j <= ar[0]; j++){
                if(ar[j]-ar[j-1] > 1 && pBase+ar[j-1] >= position) {
                    insertPos += ar[j-1] + 1;
                    got = true;
                    break;
                }
            }
            //after last element
            if(UNLIKELY(!got)){
                if(ar[ar[0]] == JacobsonIndexSize-1){
                    blockNo++;
                    if(blockNo == blocksInSegment){
                        if(!backSearchInsert(position, key, value, targetSegment, lastValidPos+position)) return false;
                        if(cardinality[targetSegment] > tree->maxLevel[0]) redistributeInsert(targetSegment, smallest[targetSegment]);
                        return true;
                    }
                    ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
                    pBase += JacobsonIndexSize;
                    insertPos = pBase;
                }else{
                    insertPos += ar[ar[0]] + 1;
                    break;
                }
            }else break;
        }
    }

    if(LIKELY(!backSearchInsert(position, key, value, targetSegment, insertPos))){
        cout<<"error in inserting"<<endl;
        exit(0);
    }
    if(cardinality[targetSegment] > tree->maxLevel[0]) {
        redistributeInsert(targetSegment, smallest[targetSegment]);
    }
    return true;
}


bool PMA::insertForward(type_t position, type_t key, type_t value, int targetSegment, int insertPos){

    type_t * movePosKey = key_chunks[targetSegment] + insertPos;
    type_t * movePosVal = value_chunks[targetSegment] + insertPos;

    insertInPosition(insertPos, targetSegment, *(movePosKey-1), *(movePosVal-1));
    movePosKey--; movePosVal--;    
    insertPos--;
    while(insertPos >= position && key < *(movePosKey-1)){
    //while(insertPos >= position){
        *movePosKey = *(movePosKey-1);
        *movePosVal = *(movePosVal-1);
        movePosKey--; movePosVal--; insertPos--;
    }
    *movePosKey = key;
    *movePosVal = value;
    return true;
}

bool PMA::insertBackward(type_t position, type_t key, type_t value, int targetSegment, int insertPos){
    
    type_t * movePosKey = key_chunks[targetSegment] + insertPos;
    type_t * movePosVal = value_chunks[targetSegment] + insertPos;
    insertInPosition(insertPos, targetSegment, *(movePosKey+1), *(movePosVal+1));
    movePosKey++; movePosVal++;    
    insertPos++;
    while(insertPos < position && key > *(movePosKey+1)){
    //while(insertPos < position){
        *movePosKey = *(movePosKey+1);
        *movePosVal = *(movePosVal+1);
        movePosKey++; movePosVal++; insertPos++;
    }
    *movePosKey = key;
    *movePosVal = value;
    return true;
}

bool PMA::insertAfterLast(type_t position, type_t key, type_t value, int targetSegment, type_t foundKey) {

    if(UNLIKELY(position >= lastValidPos)){ //Got out of the current segment. Traverse backward for vacant space
        for(type_t i = lastValidPos-1; i>=0; i--){
            int blockNo = i/JacobsonIndexSize;
            int bitPosition = i % JacobsonIndexSize;
            u_short mask =  1 << bitPosition;

            if((bitmap[targetSegment][blockNo] & mask) == 0) {
                return insertBackward(position, key, value, targetSegment, i);
            }
        }
        cout<<"Program should never reach hear at InsertAfterLast"<<endl;
        printSegElements(targetSegment);
        exit(0);
    }
    else{ //Have some space left in the segment. Go forward in the space max 3 slots
        type_t adjust = min(abs(*(key_chunks[targetSegment]+lastElementPos[targetSegment])-key), min(lastValidPos - lastElementPos[targetSegment], (type_t) MaxGap));
        //adjust = min(adjust, abs(*(key_chunks[targetSegment]+lastElementPos[targetSegment])-key));
        if(key > foundKey){
            insertInPosition(position+adjust, targetSegment, key, value);
            return true;
        }
        else{
            insertInPosition(position+adjust, targetSegment, key, value);
            swapElements(targetSegment, position, adjust);
            return true;
        }
    }
}

bool PMA::backSearchInsert(type_t position, type_t key, type_t value, int targetSegment, int forwardInsertPos) {
    type_t blockNo = position/JacobsonIndexSize; 
    u_char * ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
    type_t insertPos = blockNo * JacobsonIndexSize;
    
    while(true){
        if(ar[0] == 0){//Got a complete empty block. select the last cell
            insertPos += JacobsonIndexSize - 1;
            break;
        }
        else if(ar[0] == JacobsonIndexSize){ //Fully loaded block
            blockNo--;
            if(blockNo < 0){
                return insertForward(position, key, value, targetSegment, forwardInsertPos);
            }
            ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
            insertPos -= JacobsonIndexSize;
        }else{
            int elements = ar[0];
            bool got = false;
            if(ar[elements] < (JacobsonIndexSize - 1) && insertPos + ar[elements] < position){
                insertPos += JacobsonIndexSize - 1;
                got = true;
            }else{
                for(int j = elements; j > 1; j--){
                    if(ar[j]-ar[j-1] > 1 && insertPos + ar[j] <= position) {
                        insertPos += ar[j] - 1;
                        got = true;
                        break;
                    }
                }
            }
            if(!got){
                if(ar[1] == 0){
                    blockNo--;
                    if(blockNo < 0){
                        return insertForward(position, key, value, targetSegment, forwardInsertPos);
                    }
                    
                    ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
                    insertPos -= JacobsonIndexSize;
                }else{
                    insertPos += ar[1] - 1;
                    break;
                }
            }else break;
        }
    }
    if(abs(insertPos-position)>abs(forwardInsertPos-position)){
        return insertForward(position, key, value, targetSegment, forwardInsertPos);
    }else return insertBackward(position, key, value, targetSegment, insertPos);
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

void PMA::insertInPosition(type_t position, int targetSegment, type_t key, type_t value){
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
    if(lastElementPos[targetSegment] < position) {
        lastElementPos[targetSegment] = position;
    }

}

bool PMA::remove(type_t key){
    int targetSegment = tree->searchSegment(key);
    type_t position = findLocation1(key, targetSegment);
    int blockPosition = position/JacobsonIndexSize;
    int bitPosition = position % JacobsonIndexSize;
    u_short mask =  1 << bitPosition;
    if(!(bitmap[targetSegment][blockPosition] & mask)) return false;

    type_t * segmentOffset = key_chunks[targetSegment];
    type_t foundKey = *(segmentOffset + position);
    if(foundKey != key) return false;
    
    bitmap[targetSegment][blockPosition] &= (~mask);
    cardinality[targetSegment]--;
    //Check if the smallest of the current segment is deleted
    if(key == smallest[targetSegment]){
        u_char * ar = NonZeroEntries[bitmap[targetSegment][0]];
        int blockNo = 0;
        while(ar[0] == 0 && blockNo < blocksInSegment) ar =NonZeroEntries[bitmap[targetSegment][++blockNo]];
        smallest[targetSegment] = *(key_chunks[targetSegment]+blockNo*JacobsonIndexSize+ar[1]);
    }
    //Check if the last element is deleted
    else if(lastElementPos[targetSegment] == position){
        u_char * ar = NonZeroEntries[bitmap[targetSegment][blockPosition]];
        while(ar[0] == 0 && blockPosition >= 0) ar = NonZeroEntries[bitmap[targetSegment][--blockPosition]];
        lastElementPos[targetSegment] = blockPosition * JacobsonIndexSize + ar[ar[0]];
    }
    if(cardinality[targetSegment] < tree->minLevel[0]){
        redistributeRemove(targetSegment);
    }
    return true;
}

void PMA::redistributeRemove(int targetSegment){
    redisUpCount++;
    BPlusTree::leaf *p = tree->findLeaf(smallest[targetSegment]);
    if(UNLIKELY(p->childCount == 1)) return;

    int loc = 0, start, end;
    for(; loc<p->childCount; loc++){
        if(p->segNo[loc] == targetSegment) break;
    }
    if(loc > 0 && cardinality[targetSegment]+cardinality[p->segNo[loc-1]] < tree->maxLevel[1]){
        int totalElements = cardinality[targetSegment] + cardinality[p->segNo[loc-1]];
        p->segNo[loc-1] = mergeTwoSegments(p->segNo[loc-1], targetSegment, totalElements);
        for(int i =loc+1; i < p->childCount; i++){
            p->key[i-2] = p->key[i-1];
            p->segNo[i-1] = p->segNo[i];
        }
        p->childCount--;
        if(p->childCount < Tree_Degree/2) tree->rebalanceOnDelete(p, smallest[p->segNo[0]]);
        return;
    }
    if(loc<p->childCount-1 && cardinality[targetSegment]+cardinality[p->segNo[loc+1]] < tree->maxLevel[1]){
        int totalElements = cardinality[targetSegment] + cardinality[p->segNo[loc+1]];
        p->segNo[loc] = mergeTwoSegments(targetSegment, p->segNo[loc+1], totalElements);
        for(int i =loc+2; i < p->childCount; i++){
            p->key[i-2] = p->key[i-1];
            p->segNo[i-1] = p->segNo[i];
        }
        p->childCount--;
        if(p->childCount < Tree_Degree/2) tree->rebalanceOnDelete(p, smallest[p->segNo[0]]);
        return;
    }
    start = loc == 0 ? loc : loc-1;
    end = loc == start ? loc+1: loc;
    int totalElements = cardinality[p->segNo[start]] + cardinality[p->segNo[end]];
    redistributeTwotoTwo(p, start, end, totalElements);
    cout<<"Called two to two"<<endl;

/*
    if(totalElements < tree->minLevel[1]){
        int startLeft, endRight, curLevel = 1;
        while(true){
            int elementCount = totalElements, startCur = start, endCur = end;
            startLeft = pow(2, curLevel-1), endRight = pow(2, curLevel-1) ;
            if(startCur - startLeft < 0 ){ 
                endRight += (startLeft - startCur); 
                startLeft = startCur; 
            }
            if(endCur + endRight >= p->childCount) {
                startLeft += (endRight-(p->childCount - (endCur + 1)));
                endRight = p->childCount - (endCur + 1);
            }
            for(startCur--, startLeft--; startCur >= 0 && startLeft >= 0; startCur--, startLeft--){
                elementCount += cardinality[startCur];
            }
            startCur = max(startCur,0);
            for(endCur++, endRight--; endCur < p->childCount && endRight >= 0; endCur++, endRight-- ){
                elementCount += cardinality[endCur];
            }
            if(elementCount >= tree->minLevel[curLevel] || startLeft !=0 || endRight != 0) break;
            curLevel++; start = startCur; end = endCur; totalElements = elementCount;
        }
        //redistribute the keys from stat to end
        if(startLeft !=0 || endRight != 0) mergeMultipleSegments(p, 0, p->childCount-1, tree->findCardinality(p, this), true);
        else mergeMultipleSegments(p, start, end, (type_t)totalElements, false);
    }else{
        //Not mergeable in one segments. redistribute from two to two
        redistributeTwotoTWo(p, start, end, totalElements);
    }
*/
}

void PMA::mergeMultipleSegments(BPlusTree::leaf *p, int startLoc, int endLoc, type_t totalElements, bool isAllElements){
    return;
}

void PMA::redistributeTwotoTwo(BPlusTree::leaf *p, int startSeg, int endSeg, type_t totalElements){
    int lSeg = getSegment();
    int rSeg = getSegment();
    int curSegment = lSeg;
    int halfElement = totalElements/2;
    type_t copyBlock = 0, lastInsertkey = 0, j = 0;

    type_t *moveKeyOffset = key_chunks[startSeg];
    type_t *moveValOffset = value_chunks[startSeg];
    type_t *destKeyOffset = key_chunks[curSegment];
    type_t *destValOffset = value_chunks[curSegment];
    u_char * ar = NonZeroEntries[bitmap[startSeg][0]];
    cardinality[lSeg] = totalElements - halfElement;
    cardinality[rSeg] = halfElement;
    while(UNLIKELY(ar[0] == 0)){
        copyBlock++;
        ar = NonZeroEntries[bitmap[startSeg][copyBlock]];
    }
    moveKeyOffset += (copyBlock * JacobsonIndexSize);
    moveValOffset += (copyBlock * JacobsonIndexSize);
    lastInsertkey = smallest[lSeg] = *(moveKeyOffset + ar[1]);

    for(int blockNo = copyBlock; blockNo <blocksInSegment; blockNo++){
        ar = NonZeroEntries[bitmap[startSeg][blockNo]];
        
        for(int i = 1; i<=ar[0]; i++){
            type_t current_element = *(moveKeyOffset + ar[i]);
            j += min((elementsInSegment - j) - totalElements, min(current_element - lastInsertkey, (type_t)MaxGap));
            *(destKeyOffset + j) = lastInsertkey = current_element;
            *(destValOffset + j) = *(moveValOffset + ar[i]);
            int blockPosition = j / JacobsonIndexSize;
            int bitPosition = j % JacobsonIndexSize;
            u_short mask = 1 << bitPosition;
            bitmap[curSegment][blockPosition] |= mask;
            totalElements--;
            if(UNLIKELY(totalElements == halfElement)){
                //load the second segment
                curSegment = rSeg;
                destKeyOffset = key_chunks[rSeg];
                destValOffset = value_chunks[rSeg];
                lastElementPos[lSeg] = j;
                j = 0;
            }
        }
        bitmap[startSeg][blockNo] = 0;
        moveKeyOffset += JacobsonIndexSize;
        moveValOffset += JacobsonIndexSize;
    }

    moveKeyOffset = key_chunks[endSeg];
    moveValOffset = value_chunks[endSeg];

    for(int blockNo = 0; blockNo <blocksInSegment; blockNo++){
        ar = NonZeroEntries[bitmap[endSeg][blockNo]];
        
        for(int i = 1; i<=ar[0]; i++){
            type_t current_element = *(moveKeyOffset + ar[i]);
            j += min((elementsInSegment - j) - totalElements, min(current_element - lastInsertkey, (type_t)MaxGap));
            *(destKeyOffset + j) = lastInsertkey = current_element;
            *(destValOffset + j) = *(moveValOffset + ar[i]);
            int blockPosition = j / JacobsonIndexSize;
            int bitPosition = j % JacobsonIndexSize;
            u_short mask = 1 << bitPosition;
            bitmap[curSegment][blockPosition] |= mask;
            totalElements--;
            if(UNLIKELY(totalElements == halfElement)){
                //load the second segment
                curSegment = rSeg;
                destKeyOffset = key_chunks[rSeg];
                destValOffset = value_chunks[rSeg];
                lastElementPos[lSeg] = j;
                j = 0;
            }
        }
        bitmap[endSeg][blockNo] = 0;
        moveKeyOffset += JacobsonIndexSize;
        moveValOffset += JacobsonIndexSize;
    }
    //update the metrices
    ar = NonZeroEntries[bitmap[rSeg][0]];
    smallest[rSeg] = *(key_chunks[rSeg] + ar[1]);
    lastElementPos[rSeg] = j;
    //Now update the leaf
    for(int i=0; i<p->childCount; i++){
        if(p->segNo[i] == startSeg){
            p->segNo[i] = lSeg;
            p->segNo[i+1] = rSeg;
            p->key[i] = smallest[rSeg];
            break;
        }
    }
    freeSegID.push_back(startSeg);
    freeSegID.push_back(endSeg);
    freeSegmentCount += 2;
}

//merge all elements in one segment and return
int PMA::mergeTwoSegments(int startSeg, int endSeg, int totalElements){
    //copy elements from 1st segment
    int curSegment = getSegment();

    type_t *moveKeyOffset = key_chunks[startSeg];
    type_t *moveValOffset = value_chunks[startSeg];
    type_t *destKeyOffset = key_chunks[curSegment];
    type_t *destValOffset = value_chunks[curSegment];
    u_char * ar = NonZeroEntries[bitmap[startSeg][0]];
    type_t copyBlock = 0, lastInsertkey = 0, j = 0;
    while(UNLIKELY(ar[0] == 0)){
        copyBlock++;
        ar = NonZeroEntries[bitmap[startSeg][copyBlock]];
    }
    moveKeyOffset += (copyBlock * JacobsonIndexSize);
    moveValOffset += (copyBlock * JacobsonIndexSize);
    *destKeyOffset = lastInsertkey = *(moveKeyOffset + ar[1]);
    *destValOffset = *(moveValOffset + ar[1]);
    totalElements--;
    bitmap[curSegment][0] = 1;

    for(int i = 2; i<=ar[0]; i++){
        type_t current_element = *(moveKeyOffset + ar[i]);
          
        j += min((elementsInSegment - j) - totalElements, min((current_element - lastInsertkey), (type_t)MaxGap));
        *(destKeyOffset + j) = lastInsertkey = current_element;
        *(destValOffset + j) = *(moveValOffset + ar[i]);
        int blockPosition = j / JacobsonIndexSize;
        int bitPosition = j % JacobsonIndexSize;
        u_short mask = 1 << bitPosition;
        bitmap[curSegment][blockPosition] |= mask;
        totalElements--;
    }
    bitmap[startSeg][copyBlock] = 0;

    for(int blockno = copyBlock+1; blockno < blocksInSegment; blockno++){
        moveKeyOffset += JacobsonIndexSize;
        moveValOffset += JacobsonIndexSize;
        ar = NonZeroEntries[bitmap[startSeg][blockno]];
        for(int i = 1; i<=ar[0]; i++){
            type_t current_element = *(moveKeyOffset + ar[i]);
            
            j += min((elementsInSegment - j) - totalElements, min(current_element - lastInsertkey, (type_t)MaxGap));
            *(destKeyOffset + j) = lastInsertkey = current_element;
            *(destValOffset + j) = *(moveValOffset + ar[i]);
            int blockPosition = j / JacobsonIndexSize;
            int bitPosition = j % JacobsonIndexSize;
            u_short mask = 1 << bitPosition;
            bitmap[curSegment][blockPosition] |= mask;
            totalElements--;
        }
        bitmap[startSeg][blockno] = 0;
    }

    //Now copy elements from 2nd segment
    moveKeyOffset = key_chunks[endSeg];
    moveValOffset = value_chunks[endSeg];
    for(int blockno = 0; blockno < blocksInSegment; blockno++){
        ar = NonZeroEntries[bitmap[endSeg][blockno]];

        for(int i = 1; i<=ar[0]; i++){
            type_t current_element = *(moveKeyOffset + ar[i]);
            j += min((elementsInSegment - j) - totalElements, min(current_element - lastInsertkey, (type_t)MaxGap));
            *(destKeyOffset + j) = lastInsertkey = current_element;
            *(destValOffset + j) = *(moveValOffset + ar[i]);
            int blockPosition = j / JacobsonIndexSize;
            int bitPosition = j % JacobsonIndexSize;
            u_short mask = 1 << bitPosition;
            bitmap[curSegment][blockPosition] |= mask;
            totalElements--;
        }
        bitmap[endSeg][blockno] = 0;

        moveKeyOffset += JacobsonIndexSize;
        moveValOffset += JacobsonIndexSize;
    }

    lastElementPos[curSegment] = j;
    smallest[curSegment] = *(destKeyOffset);
    cardinality[curSegment] = cardinality[startSeg] + cardinality[endSeg];
    totalSegments--;
    freeSegID.push_back(startSeg);
    freeSegID.push_back(endSeg);
    freeSegmentCount += 2;

    return curSegment;
}

void PMA::deleteSegment(int targetSegment){
    freeSegmentCount++;
    freeSegID.push_back(targetSegment);

    fill(bitmap[targetSegment].begin(), bitmap[targetSegment].end(), 0);
    totalSegments--;
}

bool PMA::lookup(type_t key){
    int targetSegment = tree->searchSegment(key);
    type_t position = findLocation(key, targetSegment);
    int blockNo = position/JacobsonIndexSize;
    int bitPosition = position % JacobsonIndexSize;
    u_short mask =  1 << bitPosition;
    if(!(bitmap[targetSegment][blockNo] & mask)) return false;

    type_t * segmentOffsetKey = key_chunks[targetSegment];
    type_t * segmentOffsetVal = value_chunks[targetSegment];
    type_t foundKey = *(segmentOffsetKey + position);
    type_t foundVal = *(segmentOffsetVal + position);
    assert(foundKey*10 == foundVal);
    return foundKey == key ? true : false;
}

type_t PMA::findLocation(type_t key, int targetSegment){
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t start = 0, end = lastElementPos[targetSegment];
    int blockPosition, bitPosition, mask;
    type_t data, mid = 0;
    while(start <= end){
        mid = (start + end) / 2;
        blockPosition = mid / JacobsonIndexSize;
        bitPosition = mid % JacobsonIndexSize;
        mask = 1 << bitPosition;
        if((bitmap[targetSegment][blockPosition] & mask) == 0){
            int64_t changedMid = mid, offset = -1;
            while(changedMid >= start){
                changedMid += offset;
                blockPosition = changedMid / JacobsonIndexSize;
                bitPosition = changedMid % JacobsonIndexSize;
                mask = 1 << bitPosition;
                if((bitmap[targetSegment][blockPosition] & mask) !=0) break;
            }
            if(changedMid < start){
                changedMid = mid;
                offset = 1;
                while(changedMid <= end){
                    changedMid += offset;
                    blockPosition = changedMid / JacobsonIndexSize;
                    bitPosition = changedMid % JacobsonIndexSize;
                    mask = 1 << bitPosition;
                    if((bitmap[targetSegment][blockPosition] & mask) !=0 ) break;
                }
                if(changedMid > end){
                    return mid;
                }else mid = changedMid;
            }else mid = changedMid;
        }
        data = *(segmentOffset + mid);
        if(data == key) return mid;
        else if(data < key) start = mid + 1;
        else end = mid - 1;
    }
    return mid;
}

type_t PMA::findLocation1(type_t key, int targetSegment){
    type_t * segmentOffset = key_chunks[targetSegment];
    int i = 0;
    type_t start = 0;
    u_char * ar = NonZeroEntries[bitmap[targetSegment][0]];
    while(LIKELY(i < blocksInSegment)){
        int elements = ar[0], lastpos = start + ar[elements];
        if(LIKELY(elements)){
            if(*(segmentOffset+lastpos)>=key){
                //Got the target segment. Find the location in here.
                if(*(segmentOffset+start+ar[1]) >= key){
                    if(*(segmentOffset+start+ar[1]) == key) return start+ar[1];
                    return start;
                }
                // Greater than the first elemenst. Search next elements
                for(int j = 2; j<=elements; j++){
                    if(*(segmentOffset+start+ar[j]) == key) return start + ar[j];
                    else if(*(segmentOffset+start+ar[j]) > key){
                        if(ar[j]-ar[j-1]>1) return start + ar[j] - 1;
                        else return start + ar[j];
                    }
                }
                cout<<"Program must never reach here"<<endl;
                exit(0);
            }else if(UNLIKELY(lastpos == lastElementPos[targetSegment])){ // larger than all elements in current segment!
                return min(lastElementPos[targetSegment] + min((type_t)MaxGap, key-*(segmentOffset+lastElementPos[targetSegment])), lastValidPos);
            }
            ar = NonZeroEntries[bitmap[targetSegment][++i]];
            start += JacobsonIndexSize;
            if(UNLIKELY(ar[0] && *(segmentOffset+start+ar[1])>key )) return lastpos + 1;
        }
        else{
            ar = NonZeroEntries[bitmap[targetSegment][++i]];
            start += JacobsonIndexSize;
        }
    }
    //cout<<"progam should never reach here except the first case"<<endl;
    return 0;
}

type_t PMA::findLocation2(type_t key, int targetSegment){
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t start = 0;
    type_t end = lastElementPos[targetSegment];
    type_t data, mid = 0;
    int blockPosition, bitPosition, mask;
    while(start<=end){
        mid = (start+end)/2;
        blockPosition = mid / JacobsonIndexSize;
        bitPosition = mid % JacobsonIndexSize;
        mask = 1 << bitPosition;
        if((bitmap[targetSegment][blockPosition] & mask) == 0){
            type_t base = blockPosition * JacobsonIndexSize;
            bool over = false, found = false;
            while(true){
                u_char * ar = NonZeroEntries[bitmap[targetSegment][blockPosition]];
                if(UNLIKELY(ar[0] == 0)){
                    blockPosition++;
                    base += JacobsonIndexSize;
                    if(blockPosition >= blocksInSegment || base > end) break;
                    continue;
                }
                for(int i = 1; i<=ar[0]; i++){
                    if(base + ar[i] > mid) {
                        if(base + ar[i] <= end) {mid = base + ar[i]; found = true;}
                        else over = true;
                        break;
                    }
                }
                if(found || over) break;
                blockPosition++;
                base += JacobsonIndexSize;
            }
            if(!found){
                blockPosition = mid / JacobsonIndexSize;
                type_t base = blockPosition * JacobsonIndexSize;
                while(true){
                    u_char * ar = NonZeroEntries[bitmap[targetSegment][blockPosition]];
                    if(UNLIKELY(ar[0] == 0)){
                        blockPosition--;
                        base -= JacobsonIndexSize;
                        if(blockPosition < 0 || base < start) return mid;
                        continue;
                    }
                    for(int i = ar[0]; i>0; i--){
                        if(base + ar[i] < mid){
                            if(base + ar[i] < start) return mid;
                            mid = base + ar[i];
                            found = true;
                            break;
                        }
                    }
                    if(found) break;
                    blockPosition--;
                    base -= JacobsonIndexSize;
                }
            }
        }
        data = *(segmentOffset + mid);
        if(data == key) return mid;
        else if(data < key) start = mid + 1;
        else end = mid - 1;
    }
    return mid;
}

void PMA::printAllElements(){
    tree->printAllElements(this);
}

tuple<type_t, type_t> PMA::range_sum(type_t startKey, type_t endKey){
    BPlusTree::leaf *leaf = tree->findLeaf(startKey);
    int SegNo = tree->findInLeaf(leaf, startKey);
    int targetSegment = leaf->segNo[SegNo++];
    if(SegNo == leaf->childCount){leaf=leaf->nextLeaf; SegNo = 0;}

    type_t position = findLocation(startKey, targetSegment);
    type_t sum_key = 0, sum_value = 0;
    type_t blockNo = position/JacobsonIndexSize;
    u_char * ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
    type_t * segmentKeyOffset = key_chunks[targetSegment];
    type_t * segmentValOffset = value_chunks[targetSegment];
    type_t pbase = blockNo * JacobsonIndexSize;

    //Range starts somewhere within this block
    bool flag = false;
    for(int offset = 1; offset <= ar[0] ; offset++){
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
    while(LIKELY(true)){
        pbase += JacobsonIndexSize;
        blockNo++;
        if(blockNo == blocksInSegment){
            blockNo = 0;
            pbase = 0;
            targetSegment = leaf->segNo[SegNo++];
            if(SegNo == leaf->childCount){
                if(leaf->nextLeaf == NULL)return{sum_key, sum_value};
                leaf=leaf->nextLeaf;
                SegNo = 0;
            }
            //targetSegment++; //Will not work. Use leaf containing the current segment.
            if(UNLIKELY(targetSegment == totalSegments)) return{sum_key, sum_value};
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
		    return {sum_key, sum_value};
	    }
    }
    return {sum_key, sum_value};
}

void PMA::printSegElements(int targetSegment){
    type_t * key = key_chunks[targetSegment];
    type_t pBase = 0;
    for(type_t block = 0; block<blocksInSegment; block++){
        u_short bitpos = 1;
        cout <<" Bitmap: "<<bitmap[targetSegment][block]<<" ";
        for(type_t j = 0; j<JacobsonIndexSize; j++){
            if(bitmap[targetSegment][block] & bitpos)
                cout << *(key+pBase+j) << " ";
            else cout <<"0 ";
            bitpos = bitpos << 1;
        }
        pBase += JacobsonIndexSize;
    }
    cout<<"last offset: "<<lastElementPos[targetSegment] <<" Cardinality: "<<cardinality[targetSegment]<<" Total Segment: "<<totalSegments<< endl;
}

void PMA::printStat(){
    type_t totalElements = 0;
    BPlusTree::leaf * leaf = tree->leftmostLeaf(tree->root);
    while(leaf != NULL){
        for(int i=0; i<leaf->childCount; i++){
            totalElements += cardinality[leaf->segNo[i]];
        }
        leaf = leaf->nextLeaf;
    }
    //for(type_t i = 0; i<totalSegments; i++){
    //    totalElements += cardinality[i];
    //}
    cout<<"Tree Level: "<<tree->totalLevel<<endl;
    cout<<"Total elements: "<<totalElements<<endl;
    cout<<"Total Segment: "<<totalSegments<<", Free Segments: "<<freeSegmentCount<<", Elements in a Segment: "<<elementsInSegment<<endl;
    cout<<"Redistribute with insert: "<<redisInsCount<<", Redistribute with update: "<<redisUpCount<<endl;
    //cout<<"Total Shift for insert: "<<totalShiftingInsert<<", total shift for rebalance: "<<totalShiftingReb<<endl;
}


BPlusTree::BPlusTree(PMA *obj){
    root = NULL;
    calculateThreshold(obj->elementsInSegment);
}

void BPlusTree::insertInTree(int chunkNo, type_t search_key, PMA *obj){
    if(UNLIKELY(root == NULL)){
        root = new Node();
        leaf *leafNode = new leaf();
        root->child_ptr[0] = (node *)leafNode;
        root->key[0] = INT64_MAX;
        root->ptrCount = 1;
        root->nodeLeaf = true;

        leafNode->segNo[0] = chunkNo;
        leafNode->key[0] = INT64_MAX;
        leafNode->childCount = 1;
        return;
    }

    leaf *leaf = findLeaf(search_key);
    if(leaf->childCount < Tree_Degree){ //insert in leaf
        if(leaf->childCount == 1){
            int segNo = leaf->segNo[0];
            if(obj->smallest[segNo] > search_key){
                leaf->segNo[1] = leaf->segNo[0];
                leaf->segNo[0] = chunkNo;
                leaf->key[0] = obj->smallest[segNo];
            }else{
                leaf->segNo[1] = chunkNo;
                leaf->key[0] = search_key;
            }
            leaf->childCount = 2;
            return;
        }
        int position;
        for(position = leaf->childCount-1; position>0; position--){
            if(leaf->key[position-1] > search_key){
                leaf->key[position] = leaf->key[position-1];
                leaf->segNo[position+1] = leaf->segNo[position]; 
            }
            else{
                int segNo = leaf->segNo[position];
                if(obj->smallest[segNo] > search_key){
                    leaf->segNo[position+1] = leaf->segNo[position];
                    leaf->segNo[position] = chunkNo;
                    leaf->key[position] = obj->smallest[segNo];
                }else{
                    leaf->segNo[position+1] = chunkNo;
                    leaf->key[position] = search_key;
                }
                leaf->childCount++;
                return;
            }
        }
        int segNo = leaf->segNo[0];
        if(obj->smallest[segNo] > search_key){
            leaf->segNo[1] = leaf->segNo[0];
            leaf->segNo[0] = chunkNo;
            leaf->key[0] = obj->smallest[segNo];
        }else{
            leaf->segNo[1] = chunkNo;
            leaf->key[0] = search_key;
        }
        leaf->childCount++;
        return;
    }
    //Leaf is not empty. Divide.
    //Copy the key-value pairs
    type_t key_store[Tree_Degree+1];
    int segNo_store[Tree_Degree+1];
    int position;
    bool done = false;
    for(position = leaf->childCount-1; position>0; position--){
        if(leaf->key[position-1] > search_key){
            key_store[position] = leaf->key[position-1];
            segNo_store[position + 1] = leaf->segNo[position];
        }else{
            int segNo = leaf->segNo[position];
            if(obj->smallest[segNo] > search_key){
                segNo_store[position+1] = leaf->segNo[position];
                segNo_store[position] = chunkNo;
                key_store[position] = obj->smallest[segNo];
            }else{
                segNo_store[position+1] = chunkNo;
                segNo_store[position] = leaf->segNo[position];
                key_store[position] = search_key;
            }
            done = true;
            position--;
            break;
        }
    }
    if(done){
        for( ; position>=0; position--){
            segNo_store[position] = leaf->segNo[position];
            key_store[position] = leaf->key[position];
        }
    }else{
        int segNo = leaf->segNo[0];
        if(obj->smallest[segNo] > search_key){
            segNo_store[position+1] = leaf->segNo[position];
            segNo_store[position] = chunkNo;
            key_store[position] = obj->smallest[segNo];
        }else{
            segNo_store[position+1] = chunkNo;
            segNo_store[position] = leaf->segNo[position];
            key_store[position] = search_key;
        }
    }

    //Copying done. Create two nodes
    BPlusTree::leaf *newleaf = new BPlusTree::leaf();
    for(int halfLeaf = 0; halfLeaf <= Tree_Degree/2; halfLeaf++){
        leaf->segNo[halfLeaf] = segNo_store[halfLeaf];
        leaf->key[halfLeaf] = key_store[halfLeaf];
    }
    for(int nextHalf = Tree_Degree/2 + 1, index = 0; nextHalf <= Tree_Degree; nextHalf++, index++){
        newleaf->segNo[index] = segNo_store[nextHalf];
        newleaf->key[index] = key_store[nextHalf];
    }
    leaf->childCount = Tree_Degree/2 + 1;
    newleaf->childCount = Tree_Degree + 1 - leaf->childCount;
    newleaf->nextLeaf = leaf->nextLeaf;
    leaf->nextLeaf = newleaf;
    type_t key_parent = obj->smallest[leaf->segNo[0]];
    insert_in_parent(leaf,key_store[Tree_Degree/2],newleaf, key_parent);
}

void BPlusTree::insert_in_parent(void *left, type_t search_key, void *right, type_t key_parent){
    if(left == root || right == root){
        node *N = new node();
        N->child_ptr[0] = (node *)left;
        N->child_ptr[1] = (node *)right;
        N->key[0] = search_key; 
        N->ptrCount = 2;
        root = N;
        return;
    }
    node *N = findParent(root, left, key_parent);
    if(N->ptrCount < Tree_Degree){
        for(int position = N->ptrCount-1; position >= 0; position--){
            if(N->child_ptr[position] == left){
                N->child_ptr[position+1] = (node *)right;
                N->key[position] = search_key;
                N->ptrCount++;
                return;
            }else{
                N->child_ptr[position+1] = N->child_ptr[position];
                N->key[position] = N->key[position-1];
            }
        }
        cout<<"Program should never reach here. Insert in Parent node of B+ Tree"<<endl;
        exit(0);
    }
    type_t key_store[Tree_Degree+1];
    Node *ptr_store[Tree_Degree+1];
    int position;
    bool done = false;
    for(position = N->ptrCount-1; position >= 0; position--){
        if(N->child_ptr[position] == left){
            ptr_store[position+1] = (node *) right;
            ptr_store[position] = N->child_ptr[position];
            key_store[position] = search_key;
            done = true;
        }else if(done){
            ptr_store[position] = N->child_ptr[position];
            key_store[position] = N->key[position];
        }else{
            ptr_store[position+1] = N->child_ptr[position];
            key_store[position] = N->key[position-1];
        }
    }

    //Spilit records into two nodes
    node *N2 = new node();
    for(int half = 0; half <= Tree_Degree/2; half++){
        N->child_ptr[half] = ptr_store[half];
        N->key[half] = key_store[half];
    }
    for(int nextHalf = Tree_Degree/2 + 1, index = 0; nextHalf <= Tree_Degree; nextHalf++, index++){
        N2->child_ptr[index] = ptr_store[nextHalf];
        N2->key[index] = key_store[nextHalf];
    }
    N->ptrCount = Tree_Degree/2 + 1;
    N2->ptrCount = Tree_Degree + 1 - N->ptrCount;
    N2->nodeLeaf = N->nodeLeaf;
    insert_in_parent(N, key_store[Tree_Degree/2], N2, key_parent);
}

void BPlusTree::rebalanceOnDelete(leaf *p, type_t key){
    node *n = findParent(root, p, key);
    if(UNLIKELY(n->ptrCount==1)) return;
    int loc = 0;
    for( ; loc < n->ptrCount; loc++) if(n->child_ptr[loc] == (node *)p) break;
    //check if leaves can be combined
    leaf * prev = loc > 0 ? (leaf *) n->child_ptr[loc-1] : p;
    leaf * next = prev == p ? (leaf *) n->child_ptr[loc+1]: p;
    if(loc == 0) key = n->key[loc];
    else key = n->key[loc-1];
    if(prev->childCount + next->childCount < Tree_Degree){
        prev->nextLeaf = next->nextLeaf;
        prev->key[prev->childCount-1] = key;
        for(int j = 0; j<next->childCount-1; j++){
            prev->key[prev->childCount+j] = next->key[j];
            prev->segNo[prev->childCount+j] = next->segNo[j];
        }
        prev->childCount += next->childCount;
        prev->segNo[prev->childCount-1] = next->segNo[next->childCount-1];
        
        for(int pos = loc>0 ? loc+1: loc+2 ; pos<n->ptrCount; pos++){
            n->child_ptr[pos-1] = n->child_ptr[pos];
            n->key[pos-2] = n->key[pos-1];
        }
        delete next;
        n->ptrCount--;
        if(n->ptrCount < Tree_Degree/2) rebalanceOnDelete(n, key);
        return;
    }
    //Not mergeable: redistribute with neighbor
    if(loc>0){ //get a key-value pair from left segment
        for(int i=next->childCount-1; i>0; i--){
            next->segNo[i+1] = next->segNo[i];
            next->key[i] = next->key[i-1];
        }next->segNo[1] = next->segNo[0];
        next->key[0] = key;
        n->key[loc-1] = prev->key[prev->childCount-2];
        next->segNo[0] = prev->segNo[prev->childCount-1];
        prev->childCount--; next->childCount++;
        return;
    }
    //Get a key-value pair from the right segment

    prev->segNo[prev->childCount] = next->segNo[0];
    prev->key[prev->childCount-1] = n->key[loc];
    prev->childCount++;
    n->key[loc] = next->key[0];
    for(int i=1; i < next->childCount-1; i++){
        next->segNo[i-1] = next->segNo[i];
        next->key[i-1] = next->key[i];
    }
    next->childCount--;
    next->segNo[next->childCount-1] = next->segNo[next->childCount];
}

void BPlusTree::rebalanceOnDelete(node *p, type_t key){
    //CHECK THE VERY FIRST CONDITION
    node *n = findParent(root, p, key);
    if(UNLIKELY(n == NULL || n->ptrCount == 1)){
        if(n == root && !n->nodeLeaf) {
            root = p;
            delete n;
        }
        return;
    }
    int loc = 0;
    for( ; loc < n->ptrCount; loc++) if(n->child_ptr[loc] == p) break;
    //If the current node is leaf-parent and the siblings are not... then do nothing
    node *prev = NULL, *next = NULL;
    if(loc>0 && n->child_ptr[loc-1]->nodeLeaf == p->nodeLeaf){        
        prev = n->child_ptr[loc-1];
        next = p;
        key = n->key[loc-1];
    }else if(loc<n->ptrCount-1 && n->child_ptr[loc+1]->nodeLeaf == p->nodeLeaf){
        prev = p;
        next = n->child_ptr[loc+1];
        key = n->key[loc];
    }else return; //Got no same status sibling

    if(prev->ptrCount + next->ptrCount < Tree_Degree){
        prev->key[prev->ptrCount-1] = key;
        for(int j = 0; j<next->ptrCount-1; j++){
            prev->key[prev->ptrCount+j] = next->key[j];
            prev->child_ptr[prev->ptrCount+j] = next->child_ptr[j];
        }
        prev->ptrCount += next->ptrCount;
        prev->child_ptr[prev->ptrCount-1] = next->child_ptr[next->ptrCount-1];
        for(int pos = loc>0? loc+1: loc+2; pos<n->ptrCount; pos++){
            n->child_ptr[pos-1] = n->child_ptr[pos];
            n->key[pos-2] = n->key[pos-1];
        }
        delete next;
        n->ptrCount--;
        if(n->ptrCount < Tree_Degree/2) rebalanceOnDelete(n, key);
        return;
    }
    //Not mergeable: redistribute with neighbor
    if(loc>0){ //get a key-value pair from left segment
        for(int i=next->ptrCount-1; i>0; i--){
            next->child_ptr[i+1] = next->child_ptr[i];
            next->key[i] = next->key[i-1];
        }next->child_ptr[1] = next->child_ptr[0];
        next->key[0] = key;
        n->key[loc-1] = prev->key[prev->ptrCount-2];
        next->child_ptr[0] = prev->child_ptr[prev->ptrCount-1];
        prev->ptrCount--; next->ptrCount++;
        return;
    }
    //Get a key-value pair from the right segment

    prev->child_ptr[prev->ptrCount] = next->child_ptr[0];
    prev->key[prev->ptrCount-1] = key;
    prev->ptrCount++;
    n->key[loc] = next->key[0];
    for(int i=1; i < next->ptrCount-1; i++){
        next->child_ptr[i-1] = next->child_ptr[i];
        next->key[i-1] = next->key[i];
    }
    next->ptrCount--;
    next->child_ptr[next->ptrCount-1] = next->child_ptr[next->ptrCount];
}

BPlusTree::node* BPlusTree::findParent(node *r, void *n, type_t search_key){
    if(n == r) return NULL; //Root will return null as parent
    for(int position = r->ptrCount-1; position>0; position--){
        if(r->key[position-1] <= search_key){
            if(n == r->child_ptr[position]) return r;
            return findParent(r->child_ptr[position], n, search_key);
        }
    }
    if(n == r->child_ptr[0]) return r;
    return findParent(r->child_ptr[0], n, search_key);
}

BPlusTree::leaf* BPlusTree::findLeaf(type_t search_key){
    node *temp = root;
    while(!temp->nodeLeaf){
        int smallest;
        for(smallest = temp->ptrCount - 1; smallest > 0; smallest--){
            if(temp->key[smallest-1] <= search_key) break;
        }
        temp = temp->child_ptr[smallest];
    }
    for(int smallest = temp->ptrCount - 1; smallest > 0; smallest--){
        if(temp->key[smallest-1] <= search_key) return (leaf *)temp->child_ptr[smallest];
    }
    return (leaf *)temp->child_ptr[0];
}

BPlusTree::leaf* BPlusTree::findLeftSiblingLeaf(void *cur, int key){
    node * temp = findParent(root, cur, key);
    if(temp->child_ptr[0] == cur){
        while(true){
            node * temp2 = findParent(root, temp, key);
            if(temp2 == NULL) return NULL;
            if(temp2->child_ptr[0] != temp){
                for(int i=1; i<temp2->ptrCount; i++){
                    if(temp2->child_ptr[i] == temp) return rightmostLeaf(temp2->child_ptr[i-1]);
                }
            }
            temp = temp2;
        }
    }
    for(int i = 1; i<temp->ptrCount; i++){
        if(temp->child_ptr[i] == cur) return (leaf *)temp->child_ptr[i-1];
    }
    return NULL;
}

int BPlusTree::searchSegment(type_t search_key){
    leaf *leaf = findLeaf(search_key);
    if(leaf->childCount == 1) return leaf->segNo[0];
    for(int smallest = leaf->childCount - 1; smallest > 0; smallest--){
        if(leaf->key[smallest-1] <= search_key) return leaf->segNo[smallest];
    }
    return leaf->segNo[0];
}

void BPlusTree::calculateThreshold(int elements){
    int estSegments = INT32_MAX/(SEGMENT_SIZE/sizeof(type_t));
    totalLevel = (int)(ceil((log(estSegments)/log(Tree_Degree))));
    maxLevel = (double *) malloc(sizeof(double) * (totalLevel+1));
    minLevel = (double *) malloc(sizeof(double) * (totalLevel+1));

    minLevel[0] = RHO_L * elements; maxLevel[0] = TOU_L * elements;
    minLevel[totalLevel] = RHO_H * elements; maxLevel[totalLevel] = TOU_H * elements;

    for(int i = 1; i<totalLevel; i++){
        maxLevel[i] = (TOU_L - ((TOU_L-TOU_H)/(totalLevel-i))) * pow(2,i) * elements;
        minLevel[i] = (RHO_L + ((RHO_H-RHO_L)/(totalLevel-i))) * pow(2,i) * elements;
    }
}

void PMA::redistributeInsert(int segment, type_t SKey){ 
    int segNo = redistributeWithDividing(segment);
    tree->insertInTree(segNo, smallest[segNo], this);
    //return;
    /*
    BPlusTree::leaf * leaf=tree->findLeaf(SKey);
    if(UNLIKELY(leaf->childCount == 1)){
        redisInsCount++;
        int segNo = redistributeWithDividing(segment);
        tree->insertInTree(segNo, smallest[segNo], this);
        return;
    }
    int pos = tree->findInLeaf(leaf, SKey);
    int start, end, totalElements;
    start = pos>0 ? pos-1 : pos;
    end = start == pos? pos+1 : pos;
    totalElements = cardinality[leaf->segNo[start]]+cardinality[leaf->segNo[end]];
    if(totalElements > tree->maxLevel[1]){
        redisUpCount++;
        int startLeft, endRight, curLevel = 1;
        while(true){
            int elementCount = totalElements, startCur = start, endCur = end;
            startLeft = endRight = pow(2, curLevel-1);
            if(startCur - startLeft < 0 ){ 
                endRight += (startLeft - startCur); 
                startLeft = startCur; 
            }
            if(endCur + endRight >= leaf->childCount) {
                startLeft += (endRight-(leaf->childCount - (endCur + 1)));
                endRight = leaf->childCount - (endCur + 1);
            }

            for(startCur--, startLeft--; startCur>=0 && startLeft >= 0; startCur--, startLeft--){
                elementCount += cardinality[leaf->segNo[startCur]];
            }
            startCur++;
            for(endCur++, endRight--; endCur < leaf->childCount && endRight >= 0; endCur++, endRight-- ){
                elementCount += cardinality[leaf->segNo[endCur]];
            }
            endCur--;
            //cout<<"startCur: "<<startCur<<" endCur: "<<endCur<<" totalEle: "<<elementCount<<endl;
            if(elementCount <= tree->maxLevel[curLevel+1] || startLeft >= 0 || endRight >= 0) break;
            curLevel++; start = startCur; end = endCur; totalElements = elementCount;
        }
        //cout<<"start: "<<start<<" end: "<<end<<" totalEle: "<<totalElements<<" startLeft: "<<startLeft<<" endRight: "<<endRight<<endl;
        if(startLeft >=0 || endRight >= 0){
            //cout<<"Got inside Full Leaf! Leaf Child: "<<leaf->childCount<<endl;
            type_t allElements = 0;
            for(int i=0; i<leaf->childCount; i++) allElements += cardinality[leaf->segNo[i]];
            int limitElement = (tree->maxLevel[curLevel]*leaf->childCount) / pow(2, curLevel);
            if(limitElement < allElements) return redistributeNToM(leaf, 0, leaf->childCount-1, allElements);
            return redistributeNToM(leaf, start, end, (type_t)totalElements);
        }
        else redistributeNToM(leaf, start, end, (type_t)totalElements);
    }else{
        redisInsCount++;
        int segNo = redistributeWithDividing(segment);
        tree->insertInTree(segNo, smallest[segNo], this);
    }
    */
}

void PMA::redistributeNToM(BPlusTree::leaf *p, int start, int end, type_t totalElements){
    double capacity = elementsInSegment * 0.5;
    vector<int> newSegments;
    int curSegment = getSegment();
    newSegments.push_back(curSegment);
    int curCount = min((type_t)capacity, totalElements);
    totalElements -= capacity;
    cardinality[curSegment] = capacity;
    type_t *moveKeyOffset, *moveValOffset;
    type_t *destKeyOffset = key_chunks[curSegment];
    type_t *destValOffset = value_chunks[curSegment];
    u_char * ar = NonZeroEntries[bitmap[p->segNo[start]][0]];
    type_t j=0, lastInsertKey = 0;

    for(int segPos = start; segPos<=end; segPos++){
        int movSegment = p->segNo[segPos];

        moveKeyOffset = key_chunks[movSegment];
        moveValOffset = value_chunks[movSegment];
        for(int blockNo = 0; blockNo < blocksInSegment; blockNo++){
            ar = NonZeroEntries[bitmap[movSegment][blockNo]];
            for(int i=1; i<=ar[0]; i++){
                type_t current_element = *(moveKeyOffset + ar[i]);
                j += min((lastValidPos - j) - curCount, min(current_element - lastInsertKey, (type_t)MaxGap));
                *(destKeyOffset + j) = lastInsertKey = current_element;
                *(destValOffset + j) = *(moveValOffset + ar[i]);
                int blockPosition = j / JacobsonIndexSize;
                int bitPosition = j % JacobsonIndexSize;
                u_short mask = 1 << bitPosition;
                bitmap[curSegment][blockPosition] |= mask;
                curCount--;
                if(curCount == 0){
                    lastElementPos[curSegment] = j;
                    //Load a new segment
                    if(totalElements > 0){
                        curCount = min((type_t)capacity, totalElements);
                        totalElements -= curCount;
                        j = 0;
                        curSegment = getSegment();
                        newSegments.push_back(curSegment);
                        cardinality[curSegment] = curCount;
                        destKeyOffset = key_chunks[curSegment];
                        destValOffset = value_chunks[curSegment];
                    }
                }
            }
            //bitmap[movSegment][blockNo] = 0;
            moveKeyOffset += JacobsonIndexSize;
            moveValOffset += JacobsonIndexSize;
        }
        //freeSegID.push_back(movSegment);
        //freeSegmentCount++;
    }
    lastElementPos[curSegment] = j;

    for(int i=start; i<=end; i++){
        for(int j=0; j<blocksInSegment; j++)
            bitmap[p->segNo[i]][j] = 0;
    }
    for(u_int i=0; i<newSegments.size(); i++){
        ar = NonZeroEntries[bitmap[newSegments[i]][0]];
        smallest[newSegments[i]] = *(key_chunks[newSegments[i]] + ar[1]);
    }
    //Update the tree with the new segments. First copy all the segments in a buffer
    vector<int> segs;
    for(int i=0; i<start; i++) segs.push_back(p->segNo[i]);

    segs.insert(segs.end(), newSegments.begin(), newSegments.end());

    for(int i=end+1; i<p->childCount; i++) segs.push_back(p->segNo[i]);
    //Now redistribute to p and possibly a new child
    totalSegments -= p->childCount;
    totalSegments += segs.size();
    if(segs.size()>Tree_Degree){
        BPlusTree::leaf *nleaf = new BPlusTree::leaf();
        int halfSegs = segs.size()/2;
        p->segNo[0] = segs[0];
        nleaf->segNo[0] = segs[halfSegs];
        for(int i=1; i<halfSegs; i++){
            p->segNo[i] = segs[i];
            p->key[i-1] = smallest[segs[i]];
        }
        for(u_int i=1; i+halfSegs<segs.size(); i++){
            nleaf->segNo[i] = segs[i+halfSegs];
            nleaf->key[i-1] = smallest[segs[i+halfSegs]];
        }
        nleaf->nextLeaf = p->nextLeaf;
        p->nextLeaf = nleaf;
        p->childCount = halfSegs;
        nleaf->childCount = segs.size() - halfSegs;
        cout<<"prev leaf:"<<p<<" new leaf: "<<nleaf<<endl;
        cout<<"previous tree"<<endl;
        tree->showTreeStat();
        tree->insert_in_parent(p, smallest[nleaf->segNo[0]], nleaf, smallest[p->segNo[0]]);
        cout<<"after insert: tree---"<<endl;
        tree->showTreeStat();
    }else{
        p->segNo[0] = segs[0];
        for(u_int i=1; i<segs.size(); i++){
            p->segNo[i] = segs[i];
            p->key[i-1] = smallest[segs[i]];
        }
        p->childCount = segs.size();
    }
}
void PMA::checkElementWise(vector<int> usedSegment, vector<int>newSegments){
    int xblock = 0, j=1;
    u_int x=0;
    u_char * ar, *xar;
    xar = NonZeroEntries[bitmap[newSegments[x]][xblock]];
    type_t *chcekKeyOffset1, *chcekKeyOffset2;
    chcekKeyOffset2 = key_chunks[newSegments[x]];
    int compMade = 0;
    for(u_int a=0; a<usedSegment.size(); a++){
        chcekKeyOffset1 = key_chunks[usedSegment[a]];
        for(int block = 0; block<blocksInSegment; block++){
            ar = NonZeroEntries[bitmap[usedSegment[a]][block]];
            for(int i=1; i<=ar[0]; i++){
                type_t element1 = *(chcekKeyOffset1+ar[i]);
                type_t element2 = *(chcekKeyOffset2+xar[j++]);
                if(element1 != element2){
                    cout<<"got the case in mismatch3! comparison done: "<<compMade<< endl;
                    cout<<"element1: "<<element1<<" element2: "<<element2<<" newSeg: "<<newSegments[x]<<" oldSeg: "<<usedSegment[a]<<" nblock "<<xblock<<" oblock "<<block<<endl;
                    cout<<"old pos: "<<(int)ar[i]<<" new pos: "<<(int)xar[j-1]<<endl;
                    for(u_int ii = 0; ii<newSegments.size(); ii++){cout<<newSegments[ii]<<" "; printSegElements(newSegments[ii]);}
                    cout<<endl;
                    for(u_int ii = 0; ii<usedSegment.size(); ii++) {cout<<usedSegment[ii]<<" "; printSegElements(usedSegment[ii]);}
                    exit(0);
                }
                compMade++;
                //cout<<"j: "<<(j-1)<<" xar[j] "<<(int)(xar[0])<<" .. ";
                if(j > xar[0]){
                    j = 1;
                    xblock++;
                    if(xblock < blocksInSegment){ //Avoid the trailing zeros of a segment
                        xar = NonZeroEntries[bitmap[newSegments[x]][xblock]];
                        while(xar[0] == 0 && xblock < blocksInSegment) {
                            xblock++;
                            xar = NonZeroEntries[bitmap[newSegments[x]][xblock]];
                        }
                    }
                    if(xblock == blocksInSegment){
                        xblock = 0;
                        x++;
                        if(x == newSegments.size()){
                            //check for left over elements in usedSegments
                            bool flag = false;
                            if(i<ar[0]) flag = true;
                            for(int b = block+1; !flag && b<blocksInSegment; b++){
                                ar = NonZeroEntries[bitmap[usedSegment[a]][b]];
                                if(ar[0]>0) {cout<<"Got more non-zero blocks"<<endl; flag = true;}
                            }
                            if((a < usedSegment.size()-1) || flag){
                                cout<<"got the case in mismatch1! comparison done: "<<compMade<<endl;
                                if(i<ar[0]) cout<<"ar[0] has more elements left"<<endl;
                                if(a<usedSegment.size()-1) cout<<"more usedSegments left"<<endl;
                                for(u_int ii = 0; ii<newSegments.size(); ii++){cout<<newSegments[ii]<<" "; printSegElements(newSegments[ii]);}
                                cout<<endl;
                                for(u_int ii = 0; ii<usedSegment.size(); ii++){cout<<usedSegment[ii]<<" "; printSegElements(usedSegment[ii]);}
                                exit(0);
                            }
                            //New segments are over and previous segments have no extra elements
                            return;
                        }
                        else{
                            chcekKeyOffset2 = key_chunks[newSegments[x]];
                        }
                    }
                    else{
                        chcekKeyOffset2 += JacobsonIndexSize;
                    }
                    xar = NonZeroEntries[bitmap[newSegments[x]][xblock]];
                }
            }
            chcekKeyOffset1 += JacobsonIndexSize;
        }
    }
    //Check for extra element in NewSegments
    bool flag = false;

    if(j<xar[0]){
        flag = true;
        for(int ii = j; ii<=xar[0]; ii++) cout<< *(chcekKeyOffset2+xar[ii])<<" ";
    }
    for(int b = xblock+1; !flag && b<blocksInSegment; b++){
        xar = NonZeroEntries[bitmap[newSegments[x]][b]];
        chcekKeyOffset2 += JacobsonIndexSize;
        if(xar[0]>0) {
            cout<<"Got more non-zero blocks"<<endl;
            for(int iii = 1; iii<=xar[0]; iii++) cout<< *(chcekKeyOffset2+xar[iii])<<" ";
            flag = true;
        }
        cout<<"b: "<<b<<" xblock: "<<xblock<<endl;
    }
    if((x < newSegments.size()-1) || flag){
        cout<<"got the case in mismatch2! comparison done"<<compMade<<endl;
        if(j<xar[0]) cout<<"xar[0] has more elements left"<<endl;
        cout<<"J: "<<j<<" xar[0]: "<<(int)xar[0]<<" xblock "<<xblock<<" x: "<<x<<" newSegment Size: "<<newSegments.size()<<endl;
        if(x<newSegments.size()-1) cout<<"more newSegments left"<<endl;
        for(u_int ii = 0; ii<newSegments.size(); ii++) {cout<<newSegments[ii]<<" "; printSegElements(newSegments[ii]);}
            cout<<endl;
        for(u_int ii = 0; ii<usedSegment.size(); ii++) { cout<<usedSegment[ii]<<" "; printSegElements(usedSegment[ii]);}
        exit(0);
    }
}

type_t PMA::findSegmentElements(int segno){
    type_t total = 0;
    u_char * ar;
    for(int i=0; i<blocksInSegment; i++){
        ar = NonZeroEntries[bitmap[segno][i]];
        total += ar[0];
    }
    return total;
}

/*
    Returns new segment nubmer. Unsed in cases only one new segment needs to be created
 */
int PMA::redistributeWithDividing(int targetSegment){
    type_t halfElement = cardinality[targetSegment]/2;
    int curSegment = getSegment();

    type_t * moveKeyOffset = key_chunks[targetSegment];
    type_t * moveValOffset = value_chunks[targetSegment];
    type_t * destKeyOffset = key_chunks[curSegment];
    type_t * destValOffset = value_chunks[curSegment];

    int copyBlock, i, elementCount = 0;
    type_t j;

    for(copyBlock=0; copyBlock<blocksInSegment; copyBlock++){
        u_char * ar = NonZeroEntries[bitmap[targetSegment][copyBlock]];
        if((elementCount + ar[0]) >= halfElement){
            halfElement = elementCount + ar[0];
            lastElementPos[targetSegment] = copyBlock * JacobsonIndexSize + ar[ar[0]];
            copyBlock++;
            break;
        }
        elementCount += ar[0];
    }        

    //Copy the elements of current block
    elementCount = cardinality[targetSegment]-halfElement;
    u_char * ar = NonZeroEntries[bitmap[targetSegment][copyBlock]];
    while(UNLIKELY(ar[0] == 0)){
        copyBlock++;
        ar = NonZeroEntries[bitmap[targetSegment][copyBlock]];
    }
    type_t * pKeyBase = moveKeyOffset + copyBlock * JacobsonIndexSize;
    type_t * pValBase = moveValOffset + copyBlock * JacobsonIndexSize;
    type_t lastInsertkey = 0;

    *destKeyOffset = lastInsertkey = *(pKeyBase + ar[1]);
    *destValOffset = *(pValBase + ar[1]);
    bitmap[curSegment][0] = 1;
    for(i = 2, j = 0; i<=ar[0]; i++){
        type_t current_element = *(pKeyBase + ar[i]);
        j += min(lastValidPos - (j + elementCount), min(current_element - lastInsertkey, (type_t)MaxGap));
        *(destKeyOffset + j) = lastInsertkey = current_element;
        *(destValOffset + j) = *(pValBase + ar[i]);
        int blockPosition = j / JacobsonIndexSize;
        int bitPosition = j % JacobsonIndexSize;
        u_short mask = 1 << bitPosition;
        bitmap[curSegment][blockPosition] |= mask;
        elementCount--;
    }
    //lastInput = j;
    bitmap[targetSegment][copyBlock] = 0;

    //type_t lastAccessPos = lastValidPos;
    for(int blockno = copyBlock+1; blockno < blocksInSegment; blockno++){
        pKeyBase += JacobsonIndexSize;
        pValBase += JacobsonIndexSize;
        ar = NonZeroEntries[bitmap[targetSegment][blockno]];
        if(UNLIKELY(ar[0] == 0)){ 
            continue;
        }
        for(i = 1; i<=ar[0]; i++){
            type_t current_element = *(pKeyBase + ar[i]);
            
            j += min(lastValidPos - (j + elementCount), min(current_element - lastInsertkey, (type_t)MaxGap));   
            *(destKeyOffset + j) = lastInsertkey = current_element;
            *(destValOffset + j) = *(pValBase + ar[i]);
            int blockPosition = j / JacobsonIndexSize;
            int bitPosition = j % JacobsonIndexSize;
            u_short mask = 1 << bitPosition;
            bitmap[curSegment][blockPosition] |= mask;
            elementCount--;
        }
        bitmap[targetSegment][blockno] = 0;
    }

    lastElementPos[curSegment] = j;
    smallest[curSegment] = *(destKeyOffset);
    cardinality[curSegment] = cardinality[targetSegment] - halfElement;
    cardinality[targetSegment] = halfElement;
    totalSegments++;
    return curSegment;
}

type_t BPlusTree::findCardinality(leaf *l, PMA *obj){
    type_t total = 0;
    for(int i=0; i<l->childCount; i++){
        total += obj->cardinality[l->segNo[i]];
    }
    return total;
}

//Use binary search to find the position of key and segment in leaf
int BPlusTree::findInLeaf(leaf * leaf, int key){
    if(leaf->key[0]>key || leaf->childCount == 1) return 0;
    else if(leaf->key[leaf->childCount-2] <= key) return leaf->childCount-1;
    int start = 0, end = leaf->childCount-2, mid;
    while(start<end-1){
        mid = (start+end+1)/2;
        if(leaf->key[mid] < key) start = mid;
        else if(leaf->key[mid] > key) end = mid;
        else return mid+1;
    }
    return start+1;
}

BPlusTree::leaf* BPlusTree::rightmostLeaf(node *parent){
    if(parent->nodeLeaf) return (leaf *)parent->child_ptr[parent->ptrCount-1];
    while(!parent->nodeLeaf){
        parent = parent->child_ptr[root->ptrCount-1];
    }
    return (leaf *)parent->child_ptr[parent->ptrCount-1];
}

BPlusTree::leaf* BPlusTree::leftmostLeaf(node *parent){
    if(parent->nodeLeaf) return (leaf *)root->child_ptr[0];
    while(!parent->nodeLeaf){
        parent = parent->child_ptr[0];
    }
    return (leaf *)parent->child_ptr[0];
}

void BPlusTree::listSegments(vector<int> &segments, node *parent){
    if(parent->nodeLeaf){
        for(int i = 0; i<parent->ptrCount; i++){
            leaf *l = (leaf *)parent->child_ptr[i];
            for(int j = 0; j<l->childCount; j++){
                segments.push_back(l->segNo[j]);
            }
            delete l;
        }
        return;
    }
    for(int i = 0; i<parent->ptrCount; i++){
        listSegments(segments, parent->child_ptr[i]);
        delete parent->child_ptr[i];
    }
}

void BPlusTree::deleteNode(node *parent){
    if(parent->nodeLeaf){
        for(int i=0; i<parent->ptrCount; i++){
            delete (leaf *)parent->child_ptr[i];
        }
        delete parent;
    }else{
        for(int i = 0; i<parent->ptrCount; i++){
            deleteNode(parent->child_ptr[i]);
            delete parent;
        }
    }
}

void BPlusTree::deleteLeaf(leaf *l, type_t SKey){
    leaf *prev = findLeaf(SKey);
    prev->nextLeaf = l->nextLeaf;

    node *par = findParent(root, l, SKey);
    while(par->ptrCount == 1){
        if(par == root) return;
        node *par2 = findParent(root, par, SKey);
        par2->child_ptr[0] = par->child_ptr[0];
        delete par;
        par = par2;
    }

}

type_t BPlusTree::findCardinality(node *n, PMA *obj){
    type_t total = 0;
    for(int i=0; i < n->ptrCount; i++){
        if(n->nodeLeaf) total += findCardinality((leaf *)n->child_ptr[i], obj);
        else total += findCardinality(n->child_ptr[i], obj);
    }
    return total;
}

void BPlusTree::printAllElements(PMA *obj){
    int totalElements = 0;
    leaf *leaf;
    for(leaf = leftmostLeaf(root); leaf != NULL; leaf = leaf->nextLeaf){
        for(int i = 0; i<leaf->childCount; i++){
            int segNo = leaf->segNo[i];
            type_t *key = obj->key_chunks[segNo];
            type_t pBase = 0;
            for(type_t block = 0; block<obj->blocksInSegment; block++){
                cout <<" Bitmap: "<<obj->bitmap[segNo][block]<<" ";
                u_short bitpos = 1;
                for(int j = 0; j<JacobsonIndexSize; j++){
                    if(obj->bitmap[segNo][block] & bitpos)
                        cout << *(key+pBase+j) << " ";
                    else cout <<"0 ";
                    bitpos = bitpos << 1;
                }
                pBase += JacobsonIndexSize;
                if(pBase>obj->lastElementPos[segNo]) break;
            }
            totalElements += obj->cardinality[segNo];
            cout<<"SIGMENT NO: "<<segNo<<" CARDINALITY: "<<obj->cardinality[segNo]<<" LAST Position: "<<obj->lastElementPos[segNo]<<endl;
        }
    }
    cout<<"Total element inserted in the PMA: "<<totalElements<<" Total Segments: "<<obj->totalSegments<<endl;
}

void BPlusTree::printTree(vector<Node *> nodes, int level){
    if(nodes[0] == NULL) return;
    if(nodes[0]->nodeLeaf){ //Create list of leaf child nodes from the current list
        cout<<"Printing level: "<<level<<endl;
        vector<Leaf *> list;
        for(u_int i=0; i<nodes.size(); i++){
            Node *temp = nodes[i];
            cout<<"Child:- "<<(int)temp->ptrCount<<"::";
            for(int j=0; j<temp->ptrCount; j++){
                cout<<" "<<temp->child_ptr[j];
                if(j<temp->ptrCount-1) cout<<" "<<temp->key[j];
                list.push_back((Leaf *)temp->child_ptr[j]);
            }
            cout<<" || ";
        }
        cout<<endl;
        printTree(list, level+1);
    }else{
        cout<<"Printing level: "<<level<<endl;
        vector<Node *> list;
        for(u_int i=0; i<nodes.size(); i++){
            Node *temp = nodes[i];
            cout<<"Child:- "<<(int) temp->ptrCount<<"::";
            for(int j=0; j<temp->ptrCount; j++){
                cout<<" "<<temp->child_ptr[j];
                if(j<temp->ptrCount-1) cout<<" "<<temp->key[j];
                list.push_back(temp->child_ptr[j]);
            }
            cout<<" || ";
        }
        cout<<endl;
        printTree(list, level+1);
    }
}

void BPlusTree::printTree(vector<Leaf *> nodes, int level){
    cout<<"Printing level: "<<level<<endl;
    
    for(u_int i=0; i<nodes.size(); i++){
        Leaf *temp = nodes[i];
        cout<<"child:- "<<(int)temp->childCount<<"::";
        for(int j=0; j<temp->childCount; j++){
            cout<<" "<<temp->segNo[j];
            if(j<temp->childCount-1) cout<<" "<<temp->key[j];
        }
        cout<<" || ";
    }
    cout<<endl;
}

void BPlusTree::showTreeStat(){
    vector<Node *> temp;
    temp.push_back(root);
    printTree(temp, 0);
}