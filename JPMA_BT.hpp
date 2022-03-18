#ifndef JPMA_HPP_
#define JPMA_HPP_

#include <vector>
#include <tuple>
#include <array>

#include "defines.hpp"
using namespace std;

class PMA;

class BPlusTree{
public:
    typedef struct Leaf{
        type_t key[Tree_Degree-1];
        int segNo[Tree_Degree];
        short childCount = 0;
        Leaf *nextLeaf = NULL;
        //Leaf() : childCount(0), nextLeaf(NULL) {}
    }leaf;

    //No node should have a combination of child of leaf and node
    typedef struct Node{
        type_t key[Tree_Degree-1];
        Node *child_ptr[Tree_Degree];
        bool nodeLeaf = false; //Last non-leaf node has value true
        short ptrCount = 0;
        //Node() : ptrCount(0), nodeLeaf(false){}
    }node;
    node *root;
    double *minLevel, *maxLevel;
    int totalLevel;
    //int maxElementInSegment;

    BPlusTree(PMA *obj);
    leaf* findLeaf(type_t search_key);
    int findInLeaf(leaf *leaf, int SKey);
    leaf* findLeftSiblingLeaf(void *p, int key);
    inline int searchSegment(type_t search_key);
    void insertInTree(int chunkNo, type_t search_key, PMA *obj);
    void insert_in_parent(void *left, type_t search_key, void *right, type_t key_for_leaf);
    void rebalanceOnDelete(leaf *p, type_t key);
    void rebalanceOnDelete(node *n, type_t key);
    node * findParent(node *root, void *n, type_t key_parent);

    void calculateThreshold(int elements);
    void listSegments(vector<int> &segments, node *parent);
    type_t findCardinality(BPlusTree::leaf *l, PMA *obj);
    type_t findCardinality(BPlusTree::node *n, PMA *obj);

    leaf* leftmostLeaf(node *root);
    leaf* rightmostLeaf(node *root);
    void deleteNode(node *parent);
    void deleteLeaf(leaf *l, type_t SKey);
    void printAllElements(PMA *obj);
    void printTree(vector<Node *> nodes, int level);
    void printTree(vector<Leaf *> nodes, int level);
    void showTreeStat();
};

class PMA{
public:
    vector<type_t *> key_chunks;
    vector<type_t *> value_chunks;
    vector<type_t> smallest;
    vector<type_t> lastElementPos;
    vector<int> cardinality;
    int totalSegments;
    int elementsInSegment;
    u_char NonZeroEntries[JacobsonIndexCount][JacobsonIndexSize+1];
    vector<vector<u_short>> bitmap;
    BPlusTree *tree;
    type_t lastValidPos;             //Last accessible slot in each segment
    int freeSegmentCount;
    type_t blocksInSegment;
    vector<type_t *> cleanSegments;
    vector<int> freeSegID;
    int segCount;
    int redisInsCount = 0, redisUpCount = 0;

    PMA(type_t totalInsert);
    ~PMA();

    //Library functions
    bool insert(type_t key, type_t value);
    bool remove(type_t key);
    bool lookup(type_t key);
    tuple<type_t, type_t> range_sum(type_t startKey, type_t endKey);
    type_t range_sum2(type_t startKey, type_t endKey);

    //Support functions
    inline int searchSegment(type_t key);
    //tuple<type_t *, type_t *> getSegment();
    int getSegment();
    void preCalculateJacobson();
    void insertInPosition(type_t position, int targetSegment, type_t key, type_t value);
    bool backSearchInsert(type_t position, type_t key, type_t value, int targetSegment, int count);
    bool insertForward(type_t position, type_t key, type_t value, int targetSegment, int count); //Extra
    bool insertBackward(type_t position, type_t key, type_t value, int targetSegment, int count); //Extra
    bool insertAfterLast(type_t position, type_t key, type_t value, int targetSegment, type_t foundKey);
    void deleteSegment(int targetSegment);
    type_t findLocation(type_t key, int targetSegment);
    type_t findLocation1(type_t key, int targetSegment);
    type_t findLocation2(type_t key, int targetSegment);
    void redistributeInsert(int segment, type_t Skey);
    int redistributeWithDividing(int targetSegment);
    int mergeTwoSegments(int startSeg, int endSeg, int totalElements);
    void mergeMultipleSegments(BPlusTree::leaf *p, int startLoc, int endLoc, type_t totalElements, bool isAllElements);
    void redistributeTwotoTwo(BPlusTree::leaf *p, int startSeg, int endSeg, type_t totalElements);
    void redistributeNToM(BPlusTree::leaf *p, int start, int end, type_t totalElements);
    void redistributeRemove(int targetSegment);
    inline void swapElements(type_t targetSegment, type_t position, type_t adjust);

    //Testing functions
    void printStat();
    void printAllElements();
    void printSegElements(int targetSegment);
    type_t findSegmentElements(int segno);

    void checkElementWise(vector<int> usedSegment, vector<int>newSegments);
};

#endif