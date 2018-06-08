#include <stdio.h>
#include <stdint.h>
#include "TFile.h"
#include "TH1D.h"

extern "C" {
	#include "../SFMT-src-1.4.1/SFMT.h"
	#include "./HistGen.h"
}

int main(int argc, char** argv) {
	TFile* f = new TFile("./waveformTest.root", "READ");
	TH1D* rHist = (TH1D *)f->Get("hist");
	uint64_t* hist = new uint64_t[rHist->GetNbinsX()];
	for(int i = 0; i < rHist->GetNbinsX(); i++) {
		hist[i] = rHist->GetBinContent(i+1);
	}
	
	TH1D* myHist = new TH1D("name", "title", 40000, 0, 40000);
	//uint64_t hist[6] = {10,10,5,7,8,8};
	struct histGen* gen = createGen(hist, rHist->GetNbinsX());
	
//	unsigned int numBoxes = gen->prob[7];
//	for(int j = 0; j < gen->nEntries; j++) {
//		if(gen->alias[j] == 7) {
//			numBoxes += gen->cutoff - gen->prob[j];
//		}
//	}
//	printf("NumBoxes: %u\n", numBoxes);
//	return 0;
	
	sfmt_t sfmt;
	sfmt_init_gen_rand(&sfmt, 4321);
	
	uint64_t* us;
	us = (uint64_t *)malloc(sizeof(uint64_t)*102400);
	sfmt_fill_array64(&sfmt, us, 102400);
	
	int i = 0;
	unsigned int n = 0;
	unsigned int r = 0;;
	while(n < 10000000) {
		if(i >= 102400) {sfmt_fill_array64(&sfmt, us, 102400); i = 0;}
		uint32_t f = genIndex(gen, us[i]);
		i++;
		while(f == gen->nEntries-1) {
			r++;
			if(i >= 102400) {sfmt_fill_array64(&sfmt, us, 102400); i = 0;}
			f = genIndex(gen, us[i]);
			i++;
		}
		myHist->Fill(f);
		n++;
	}
	
	printf("%u,%u,%f\n", r, n, ((double)r)/((double)n));
	
	myHist->SaveAs("./TestHist.root");
	
	freeGen(gen);
	delete[] hist;
	
	return 1;
}