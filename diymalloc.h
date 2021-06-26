#ifndef _DIYMALLOC_H_
#define _DIYMALLOC_H_

#include "jemalloc.h"

//define to jemalloc
#define malloc(size) je_malloc(size)
#define calloc(count,size) je_calloc(count,size)
#define realloc(ptr,size) je_realloc(ptr,size)
#define free(ptr) je_free(ptr)

#endif