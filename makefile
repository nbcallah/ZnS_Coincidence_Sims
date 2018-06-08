#CFLAGS = -g -DDEBUG
CFLAGS = -O3
SFMTFLAGS = -msse2 -DHAVE_SSE2 -DSFMT_MEXP=19937

PARALLEL := $(shell command -v mpic++ 2> /dev/null)

debug: CFLAGS = -g -DDEBUG
debug: all

all:
ifndef PARALLEL
all: single
else
all: parallel
endif

single: CPP = g++
single: CC = gcc
single: LOC = -DLOCAL
single: randomGenerator-2016-2017_Efficiency

parallel: CPP = mpic++
parallel: CC = mpicc
parallel: LOC = LOCAL
parallel: PAR = -DPARALLEL
parallel: randomGenerator-2016-2017_Efficiency

#all: randomGenerator-2016-2017_BRII_dip1 randomGenerator-2016-2017_BRII_dip3 randomGenerator-2016-2017_BRII_dip9 efficiency-2016-2017_BRII randomGenerator-2016-2017_Efficiency randomGenerator-2016-2017_Efficiency_NoPileup randomGenerator-2016-2017_Efficiency_NoPileup_FixedWin

randomGenerator-2016-2017_Efficiency: randomGenerator-2016-2017_Efficiency.o SFMT-src-1.4.1/SFMT.o SFMT-src-1.4.1/jump/SFMT-jump.o HistGen/HistGen.o
	$(CPP) -o randomGenerator-2016-2017_Efficiency randomGenerator-2016-2017_Efficiency.o SFMT-src-1.4.1/SFMT.o SFMT-src-1.4.1/jump/SFMT-jump.o HistGen/HistGen.o
	
randomGenerator-2016-2017_Efficiency.o: randomGenerator-2016-2017_Efficiency.cpp
	$(CPP) -std=c++11 $(CFLAGS) $(SFMTFLAGS) $(PAR) -D$(LOC) -ISFMT-src-1.4.1/ -c -o randomGenerator-2016-2017_Efficiency.o randomGenerator-2016-2017_Efficiency.cpp

HistGen/HistGen.o: HistGen/HistGen.c
	$(CC) $(CFLAGS) -c -o HistGen/HistGen.o HistGen/HistGen.c

SFMT-src-1.4.1/SFMT.o: SFMT-src-1.4.1/SFMT.c
	$(CC) $(CFLAGS) $(SFMTFLAGS) -c -o SFMT-src-1.4.1/SFMT.o SFMT-src-1.4.1/SFMT.c

SFMT-src-1.4.1/jump/SFMT-jump.o: SFMT-src-1.4.1/jump/SFMT-jump.c
	$(CC) $(CFLAGS) $(SFMTFLAGS) -ISFMT-src-1.4.1/ -c -o SFMT-src-1.4.1/jump/SFMT-jump.o SFMT-src-1.4.1/jump/SFMT-jump.c
	
clean:
	rm -rf randomGenerator-2016-2017_BRII_dip1 randomGenerator-2016-2017_BRII_dip1.o randomGenerator-2016-2017_BRII_dip3 randomGenerator-2016-2017_BRII_dip3.o randomGenerator-2016-2017_BRII_dip9 randomGenerator-2016-2017_BRII_dip9.o SFMT-src-1.4.1/jump/SFMT-jump.o SFMT-src-1.4.1/SFMT.o HistGen/HistGen.o efficiency-2016-2017_BRII.o randomGenerator-2016-2017_Efficiency randomGenerator-2016-2017_Efficiency.o randomGenerator-2016-2017_Efficiency_NoPileup randomGenerator-2016-2017_Efficiency_NoPileup.o randomGenerator-2016-2017_Efficiency_NoPileup_FixedWin randomGenerator-2016-2017_Efficiency_NoPileup_FixedWin.o
