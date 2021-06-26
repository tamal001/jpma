#ifndef DEFINES_HPP_
#define DEFINES_HPP_

#define SEGMENT_SIZE 2097152

#define MAX_SEGMENTS 32768

#define BLOCKS_IN_SEGMENT 16384

#define type_t int64_t

#define JacobsonIndexSize 16
#define JacobsonIndexCount 65536

#define ProbLimit JacobsonIndexSize // * 4
#define MaxGap 3

#define MinThreshold 16 //Change
#define MaxThreshold 1000000 //Change to 95% of total key in a segment
#endif