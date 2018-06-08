#include "HistGen.h"
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>

struct histGen* createGen(uint64_t* histogram, uint32_t n) {
	if(n == UINT32_MAX) {
		return NULL;
	}
	
	struct histGen* gen = (struct histGen *)malloc(sizeof(struct histGen));

	uint64_t sum = 0;
	for(uint32_t i = 0; i < n; i++) {
		sum += histogram[i];
	}

	gen->nEntries = n+1;
	uint64_t cutoff = (sum / (n+1)) + 1;
	if(cutoff > UINT32_MAX) {
		return NULL;
	}
	gen->cutoff = cutoff;
	
	uint64_t extras = cutoff * (n+1) - sum;
	
	gen->maxIndex = (n+1) * (UINT32_MAX / (n+1));
	gen->maxProb = cutoff * (UINT32_MAX / cutoff);
	
	gen->alias = (uint32_t *)malloc(sizeof(uint32_t) * (n+1));
	gen->prob = (uint32_t *)malloc(sizeof(uint32_t) * (n+1));
	
	uint64_t* smallProb = (uint64_t *)malloc(sizeof(uint64_t) * (n+1));
	uint32_t* smallIndex = (uint32_t *)malloc(sizeof(uint32_t) * (n+1));
	uint64_t* largeProb = (uint64_t *)malloc(sizeof(uint64_t) * (n+1));
	uint32_t* largeIndex = (uint32_t *)malloc(sizeof(uint32_t) * (n+1));
	uint32_t nSmall = 0;
	uint32_t nLarge = 0;
	
	for(uint32_t i = 0; i < n; i++) {
		if(histogram[i] < cutoff) {
			smallProb[nSmall] = histogram[i];
			smallIndex[nSmall] = i;
			nSmall++;
		}
		else {
			largeProb[nLarge] = histogram[i];
			largeIndex[nLarge] = i;
			nLarge++;
		}
	}

	if(extras < cutoff) {
		smallProb[nSmall] = extras;
		smallIndex[nSmall] = n;
		nSmall++;
	}
	else {
		largeProb[nLarge] = extras;
		largeIndex[nLarge] = n;
		nLarge++;
	}
	
	while(nSmall > 0) {
		gen->alias[smallIndex[nSmall-1]] = largeIndex[nLarge-1];
		gen->prob[smallIndex[nSmall-1]] = smallProb[nSmall-1];
		
		largeProb[nLarge-1] -= (cutoff - smallProb[nSmall-1]);
		
		nSmall--;
		
		if(largeProb[nLarge-1] < cutoff) {
			smallIndex[nSmall] = largeIndex[nLarge-1];
			smallProb[nSmall] = largeProb[nLarge-1];
			nSmall++;
			nLarge--;
		}
	}

	while(nLarge > 0) {
		gen->alias[largeIndex[nLarge-1]] = largeIndex[nLarge-1];
		gen->prob[largeIndex[nLarge-1]] = cutoff;
		nLarge -= 1;
	}
		
	free(smallProb);
	free(smallIndex);
	free(largeProb);
	free(largeIndex);
	
	return gen;
}

uint32_t genIndex(struct histGen* gen, uint64_t u) {
	uint32_t a;
	uint32_t b;
	a = u >> 32;
	b = (u << 32) >> 32;
	
	if(a > gen->maxIndex || b > gen->maxProb) {
		return gen->nEntries-1;
	}
	
	uint32_t index = a % (gen->nEntries);
	uint32_t prob = b % (gen->cutoff) + 1;
	
	return prob <= gen->prob[index] ? index : gen->alias[index];
}

void freeGen(struct histGen* gen) {
	free(gen->alias);
	free(gen->prob);
	free(gen);
}