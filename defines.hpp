#ifndef DEFINES_HPP_
#define DEFINES_HPP_

//#define SEGMENT_SIZE 2097152
#define CHUNK_SIZE 2097152UL
//#define CHUNK_SIZE 262144
//#define SEGMENT_SIZE 16384//1024//16384//32768

#define SEGMENT_SIZE 1024
//#define TOTAL_SEGMENTS 65536

#define type_t int64_t

#define LIKELY(x) __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)

#define PROTECTION (PROT_READ | PROT_WRITE)

#ifndef MAP_HUGETLB
#define MAP_HUGETLB 0x40000 /* arch specific */
#endif

//1 for mmap, 2 for malloc
#ifndef Allocation_type
#define Allocation_type 2
#endif

#define Tree_Degree 16

//Follows 0 <= rho(l) <= rho(h) <= tou(h) <= tou(l) <= 1
#define RHO_L 0.20
#define RHO_H 0.50
#define TOU_L 0.95
#define TOU_H 0.75

#ifndef Insert_type
#define Insert_type 2
#endif

#ifdef __ia64__
#define ADDR (void *)(0x8000000000000000UL)
#define FLAGS (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_FIXED)
#else
#define ADDR (void *)(0x0UL)
//#define FLAGS (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB)
#define FLAGS (MAP_SHARED |MAP_ANONYMOUS | MAP_HUGETLB)
#endif

#define JacobsonIndexSize 16
#define JacobsonIndexCount 65536
#define MaxGap 3

#endif