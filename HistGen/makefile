CC = gcc
CPP = g++
CFLAGS = -g
RCFLAGS = `root-config --cflags`
RLDFLAGS = `root-config --libs`

all: test

clean:
	find . -type f -name '*.o' -delete
	rm -f test
	
test: test.cpp histGen.o ../SFMT-src-1.4.1/SFMT.o
	$(CPP) $(CFLAGS) $(RCFLAGS) -o test test.cpp ./histGen.o ../SFMT-src-1.4.1/SFMT.o -DSFMT_MEXP=19937 $(RLDFLAGS)
	
histGen: histGen.c
	$(CC) $(CFLAGS) -c -o histGen.o histGen.c
	
../SFMT-src-1.4.1/SFMT.o: ../SFMT-src-1.4.1/SFMT.c
	$(CC) $(CFLAGS) -c -DSFMT_MEXP=19937 -o ../SFMT-src-1.4.1/SFMT.o ../SFMT-src-1.4.1/SFMT.c
