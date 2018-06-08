#ifndef HISTGEN_H
#define HISTGEN_H

#include<stdint.h>

struct histGen {
	uint32_t nEntries;
	uint32_t cutoff;  //Will be used for random generation
	uint32_t *alias; //Will hold array of indices
	uint32_t *prob;   //Will hold array of probabilities to determine whether alias1 or 2
	uint32_t maxIndex;
	uint32_t maxProb;
};

struct histGen* createGen(uint64_t* histogram, uint32_t n);

uint32_t genIndex(struct histGen* gen, uint64_t u);

void freeGen(struct histGen* gen);

#endif