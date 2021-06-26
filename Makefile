CC=g++
CFLAGS=-Wall -g -O3 -std=c++17
INCLUDES=-I ./include/
ALLOC_DEP=./lib/libjemalloc.a
ALLOC_LINK=$(ALLOC_DEP) -lpthread -ldl

PROGRAMS = benchmark

all: $(PROGRAMS)

jpma: 
	$(CC) $(INCLUDES) $(CFLAGS) -c JPMA.cpp -o jpma.o 
# $(ALLOC_LINK)

benchmark: jpma
	$(CC) $(INCLUDES) $(CFLAGS) jpma.o benchmark.cpp -o benchmark $(ALLOC_LINK)

clean:
	rm -f benchmark jpma.o