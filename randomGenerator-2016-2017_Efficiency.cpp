#include <stdio.h>
#include <stdexcept>
#include <cmath>
#include <getopt.h>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <mpi.h>
#include <limits>
#include <functional>

extern "C" {
	#include "SFMT-src-1.4.1/SFMT.h"
	#include "SFMT-src-1.4.1/jump/SFMT-jump.h"
	//#include "SFMT-src-1.4.1/jump/SFMT-calc-jump.hpp"
	#include "HistGen/HistGen.h"
}

#include "jumpPoly.h"

#ifdef LOCAL
	#define LOCATION "./2016_2017_Waveform_Data"
#endif

#define LENGTH 100.0

#define BKGRATEA 560
#define BKGRATEB 545

#define COINCBKGRATE 4.89

#define PMTAPROB 0.584975

#define NANOSECOND .000000001
#define CLOCKTONS 0.8

#define COINCWINDOW 50.0

struct event {
	//long time;
	double realtime;
	int ch;
	//int full;
	//int edge;
	//int tag;
};

struct ph {
	int a;
	int b;
};

#ifdef DEBUG
	template class std::vector<double>;
	template class std::vector<event>;
#endif

std::vector<ph> phsVect;
std::vector<ph> phsVect_bkg;

histGen* ch1_1;
histGen* ch1_2;
histGen* ch2_1;
histGen* ch2_2;
histGen* phs;
histGen* ch1_1_bkg;
histGen* ch1_2_bkg;
histGen* ch2_1_bkg;
histGen* ch2_2_bkg;
histGen* phs_bkg;

template<typename T, typename F>
void loadTextHistogram(histGen* &gen, std::string fName, std::vector<T> &vec, F f);

template<typename T, typename F>
void fillRandomArray(sfmt_t* sfmt, T* ary, int num, F f);

std::vector<event> createSynthRun(sfmt_t* sfmt, double rate, std::function<double(sfmt_t*)> pop, std::function<uint64_t(sfmt_t*)>pop64, long &numCreated);

double sumCtsCorrected(std::vector<event> data, int numPH, double window) {
    int i;
	int cur;
	int tailIt;
    
    double numCoinc;
    
	for(i = 0; i < data.size(); i++) { //Look at each entry of the data vector
		int ch1PESum = 0;
		int ch2PESum = 0;
		event prevEvt;
				
		if(data.at(i).ch != 1 && data.at(i).ch != 2) { continue; } //Skip all but dagger for coincidences
		
		//Put these 05/06/2017 - Should get correct order for pmtAHits
		if(data.at(i).ch == 1) {
			ch1PESum += 1;
		}
		if(data.at(i).ch == 2) {
			ch2PESum += 1;
		}
		
		for(cur = i+1; cur < data.size(); cur++) { //Now we'll search forward for a coincidence
			if(data.at(cur).realtime - data.at(i).realtime > COINCWINDOW*NANOSECOND) {
				break; //If we exceed the coincidence window, continue at the event after i
			}
			if(data.at(cur).ch != 1 && data.at(cur).ch != 2) { continue; } //ignore non-dagger
			if(data.at(cur).ch == data.at(i).ch) { //If we find a second event in channel i before coincidence, count it
				if(data.at(i).ch == 1) {
					ch1PESum += 1;
				}
				if(data.at(i).ch == 2) {
					ch2PESum += 1;
				}
			}
			if(data.at(cur).ch != data.at(i).ch) { //Only trigger coincidence on different channel

				if(data.at(cur).ch == 1) {
					ch1PESum += 1;
				}
				if(data.at(cur).ch == 2) {
					ch2PESum += 1;
				}
				
				prevEvt = data.at(cur);
				
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) { //Now integrate the tail
					double diff = (data.at(tailIt).realtime - prevEvt.realtime);
					if((data.at(tailIt).realtime - prevEvt.realtime) > window*NANOSECOND) {
						break;
					}
					if(data.at(tailIt).ch != 1 && data.at(tailIt).ch != 2) { continue; }
					if(data.at(tailIt).ch == 1) {
						ch1PESum += 1;
					}
					if(data.at(tailIt).ch == 2) {
						ch2PESum += 1;
					}
					prevEvt = data.at(tailIt);
				}
				if(ch1PESum + ch2PESum >= numPH) {
					numCoinc += 1;
					i = tailIt-1; //Only impose deadtime on counted neutrons; Otherwise continue looking.
					cur = tailIt-1;
					break;
				}
				break;
			}
//			break;
		}
	}
    return numCoinc;
}

int main(int argc, char** argv) {
	int c;

//	double rate = 0.0;
	double precision = 0.0;

    //Read in command line options
	while (1)
	{
		static struct option long_options[] =
		{
			/* These options don't set a flag.
			 We distinguish them by their indices. */
//			{"rate", required_argument, 0, 'r'},
			{"precision", required_argument, 0, 'p'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "p", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1) { break; }

		switch (c)
		{
			case 0:
				break;
			
//			case 'r':
//				rate = atof(optarg);
//				break;
				
			case 'p':
				precision = atof(optarg);
				break;

			default:
				abort();
		}
	}

	if(precision == 0.0) {
		fprintf(stderr, "Error! Usage: ./randomGenerator --precision=final_rel_precision\n");
		exit(1);
	}
	
    //Set up PRNG depending on if parallelizing or not
    #ifdef PARALLEL
        MPI_Init(NULL, NULL);
         // Get the number of processes
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        // Get the rank of the process
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
        sfmt_t sfmt;
        uint32_t seed[4] = {2736687128, 234302120, 3355772407, 657836083};
        sfmt_init_by_array(&sfmt, seed, 4);
        //sfmt_init_gen_rand(&sfmt, (uint32_t)4321);

        for(int i = 0; i < world_rank; i++) {
            SFMT_jump(&sfmt, jumpPolynomial);
        }
    #else
        sfmt_t sfmt;
        uint32_t seed[4] = {2736687128, 234302120, 3355772407, 657836083};
        sfmt_init_by_array(&sfmt, seed, 4);
        //sfmt_init_gen_rand(&sfmt, (uint32_t)4321);
    #endif
	
    //Read in ZnS waveform data to resample
	std::vector<unsigned long> dump;

	auto waveformFunc = [](std::string s)->unsigned long{return std::stoul(s);};
	auto phsFunc = [](std::string s)->ph{
		std::stringstream ss(s);
		std::string a, b;
		std::getline(ss, a, ':');
		std::getline(ss, b, ':');
		return ph {std::stoi(a), std::stoi(b)};
	};
	
	loadTextHistogram(ch1_1, LOCATION"/WF-ch1-1start.csv", dump, waveformFunc);
	dump.clear();
	loadTextHistogram(ch1_2, LOCATION"/WF-ch1-2start.csv", dump, waveformFunc);
	dump.clear();
	loadTextHistogram(ch2_1, LOCATION"/WF-ch2-1start.csv", dump, waveformFunc);
	dump.clear();
	loadTextHistogram(ch2_2, LOCATION"/WF-ch2-2start.csv", dump, waveformFunc);
	dump.clear();
	loadTextHistogram(phs, LOCATION"/PHSTextHist.csv", phsVect, phsFunc);

	loadTextHistogram(ch1_1_bkg, LOCATION"/WF-ch1-1start_bkg.csv", dump, waveformFunc);
	dump.clear();
	loadTextHistogram(ch1_2_bkg, LOCATION"/WF-ch1-2start_bkg.csv", dump, waveformFunc);
	dump.clear();
	loadTextHistogram(ch2_1_bkg, LOCATION"/WF-ch2-1start_bkg.csv", dump, waveformFunc);
	dump.clear();
	loadTextHistogram(ch2_2_bkg, LOCATION"/WF-ch2-2start_bkg.csv", dump, waveformFunc);
	dump.clear();
	loadTextHistogram(phs_bkg, LOCATION"/PHSTextHist_bkg.csv", phsVect_bkg, phsFunc);
	
	int i = 0;
	
	//Create a little object that acts as a buffer of random numbers.
	int numPops = 0;
	double *src = new double[1000];
	fillRandomArray(&sfmt, src, 1000,
					[](uint64_t val){return sfmt_to_res53(val);});
	auto pop = [&numPops, src](sfmt_t* sfmt)->double{
		//printf("%d\n", numPops);
		if(numPops >= 1000) {
			fillRandomArray(sfmt, src, 1000, [](uint64_t val){return sfmt_to_res53(val);});
			numPops = 0;
		}
		double ret = src[numPops];
		numPops++;
		return ret;
	};
	
	int numPops64 = 0;
	uint64_t *src64 = new uint64_t[1000];
	fillRandomArray(&sfmt, src64, 1000,
					[](uint64_t val){return val;});
	auto pop64 = [&numPops64, src64](sfmt_t* sfmt)->uint64_t{
		//printf("%d\n", numPops);
		if(numPops64 >= 1000) {
			fillRandomArray(sfmt, src64, 1000, [](uint64_t val){return val;});
			numPops64 = 0;
		}
		uint64_t ret = src64[numPops64];
		numPops64++;
		return ret;
	};
    
    
    
    //Start the simulation
    int numPh = 8;
    double window = 500.0;
    
    double mean = 0.0;
	double stdDev = 0.0;
    std::vector<double> bkgCts;
    int numRuns = 0;
    //First measure the background rate by generating lots of runs with 1 nanoHz rate
    do {
        long createdCts;
        std::vector<event> vec = createSynthRun(&sfmt, 1e-9, pop, pop64, createdCts);
        double cts = sumCtsCorrected(vec, numPh, window);
        bkgCts.push_back(cts);
        numRuns += 1;

        mean = std::accumulate(bkgCts.begin(), bkgCts.end(), 0.0, [](double sum, double e){return sum+e;})/bkgCts.size();
        stdDev = sqrt(std::accumulate(bkgCts.begin(), bkgCts.end(), 0.0, [mean](double sum, double e){return sum+pow((e-mean), 2);})/(bkgCts.size()-1));
        stdDev = numRuns >= 10 ? stdDev : std::numeric_limits<double>::infinity();
//        printf("mean: %f stdDev: %f precision: %f\n", mean/100, stdDev/100, (stdDev/sqrt(bkgCts.size()))/mean);
    } while((stdDev/sqrt(bkgCts.size()))/mean > 0.005); //relative precision of 0.5% on bkg. rate
    
    double bkgRate = mean/100;
//    printf("bkg: %f (%d runs)\n", bkgRate, numRuns);
    
    std::vector<double> efficiencies;
    std::vector<double> efficiencies_err;
    std::vector<double> rates;
    
    //Check efficiency for 3 different rates
    for(int i = 100; i < 14000; i+= 5000) {
        double sumCts = 0.0;
        int ltNum = 0;
        std::vector<double> eff;
        
        double rate = (double)i;
        //We'll loop until the desired precision is reached
        do {
            long createdCts;
            ltNum += 1;

            //Create a run with a constant rate. Will report back the total number of UCN created during the simulation in createdCts
            std::vector<event> vec = createSynthRun(&sfmt, rate, pop, pop64, createdCts);

            ///!\/!\/!\/!\/!\/!\/!\
            //Put your code here!
            ///!\/!\/!\/!\/!\/!\/!\
            //Function takes vec and returns deadtime and pileup corrected sum of counts in the run
            double cts = sumCtsCorrected(vec, numPh, window);
            ///!\/!\/!\/!\/!\/!\/!\
            //Put your code here!
            ///!\/!\/!\/!\/!\/!\/!\
            

            cts -= bkgRate * LENGTH; //Subtrack bkg previously measured
            
            sumCts += cts; //keep track of sum

            eff.push_back(cts/((double)createdCts)); //Add most recent measurement to list

            //Calc mean and stddev. stddev/sqrt(n) is the error
            mean = std::accumulate(eff.begin(), eff.end(), 0.0, [](double sum, double e){return sum+e;})/eff.size();
            stdDev = sqrt(std::accumulate(eff.begin(), eff.end(), 0.0, [mean](double sum, double e){return sum+pow((e-mean), 2);})/(eff.size()-1));

            stdDev = ltNum >= 10 ? stdDev : std::numeric_limits<double>::infinity();

        } while((stdDev/sqrt(eff.size()))/mean > precision);
        
        efficiencies.push_back(mean);
        efficiencies_err.push_back((stdDev/sqrt(eff.size())));
        rates.push_back(rate);
        
        //Print for Python program
        printf("%.12f %.12f %.12f\n", rate, mean, (stdDev/sqrt(eff.size())));

//        printf("Final : %f,%.12f,%.12f,%d,%d,%f\n", rate, mean, (stdDev/sqrt(eff.size())), ltNum, 8, window);
//        printf("Sum: %f,%d,%d,%f\n", sumCts, ltNum, numPh, window);
    }
    
    

	delete[] src;
	delete[] src64;
	freeGen(ch1_1);
	freeGen(ch1_2);
	freeGen(ch2_1);
	freeGen(ch2_2);
	freeGen(phs);
	freeGen(ch1_1_bkg);
	freeGen(ch1_2_bkg);
	freeGen(ch2_1_bkg);
	freeGen(ch2_2_bkg);
	freeGen(phs_bkg);
	
    #ifdef PARALLEL
        MPI_Finalize();
    #endif
}

template<typename T, typename F>
void loadTextHistogram(histGen* &gen, std::string fName, std::vector<T> &vec, F f) {
	std::vector<uint64_t> hist;
	std::ifstream infile(fName);
	std::string line;
	while(std::getline(infile, line)) {
		std::istringstream iss(line);
		std::vector<std::string> strings;
		std::string s;
		while(std::getline(iss, s, ',')) {
			strings.push_back(s);
		}
		if(strings.size() != 2) {
			gen = NULL;
			vec.empty();
		}
		hist.push_back(std::stoull(strings[0]));
		vec.push_back(f(strings[1]));
	}
	assert(hist.size() > 0);
	gen = createGen(&(hist[0]), hist.size());
}

template<typename T, typename F>
void fillRandomArray(sfmt_t* sfmt, T* ary, int num, F f) {
	//numInts += num;
	int i = 0;
	int sourceNum = num;
	sourceNum = sourceNum > 312 ? sourceNum : 312;
	sourceNum = sourceNum % 2? sourceNum + 1 : sourceNum;
	uint64_t* sourceAry = new uint64_t[sourceNum];
	sfmt_fill_array64(sfmt, sourceAry, sourceNum);
	for(i=0; i < num; i++) {
		ary[i] = f(sourceAry[i]);
	}
	delete[] sourceAry;
}

std::vector<event> createSynthRun(sfmt_t* sfmt, double rate, std::function<double(sfmt_t*)> pop, std::function<uint64_t(sfmt_t*)>pop64, long &numCreated) {
//	double rate = num / LENGTH;
	//printf("%f\n", rate);
	
	std::vector<event> evts;
	
	event evt;
	
	numCreated = 0;
	
	double coincHitTime = 0.0;
	while(true) {
		coincHitTime += -log(pop(sfmt))/rate;
		if(coincHitTime > LENGTH) {
			break;
		}
		numCreated += 1;
		uint32_t index = genIndex(phs, pop64(sfmt));
		while(index == phs->nEntries - 1) {
			index = genIndex(phs, pop64(sfmt));
		}
		int numA = phsVect[index].a;
		int numB = phsVect[index].b;
		//printf("%f - %d:%d\n", coincHitTime, numA, numB);
		int srcChan = pop(sfmt) < PMTAPROB ? 1 : 2;
		histGen* ch1 = srcChan == 1 ? ch1_1 : ch1_2;
		histGen* ch2 = srcChan == 1 ? ch2_1 : ch2_2;
		if(srcChan == 1) {numA--; evts.push_back(event {coincHitTime, 1});} //Automatically put hit at t=0
		if(srcChan == 2) {numB--; evts.push_back(event {coincHitTime, 2});} //Automatically put hit at t=0
		evt.realtime = 0.0;
		evt.ch = 0;
		for(int i = 0; i < numA; i++) {
			uint32_t index = genIndex(ch1, pop64(sfmt));
			while(index == ch1->nEntries - 1) {
				index = genIndex(ch1, pop64(sfmt));
			}
			evt.ch = 1;
			evt.realtime = coincHitTime + index * CLOCKTONS * NANOSECOND;// * 0.01;
			evts.push_back(evt);
		}
		for(int i = 0; i < numB; i++) {
			uint32_t index = genIndex(ch2, pop64(sfmt));
			while(index == ch2->nEntries - 1) {
				index = genIndex(ch2, pop64(sfmt));
			}
			evt.ch = 2;
			evt.realtime = coincHitTime + index * CLOCKTONS * NANOSECOND;// * 0.01;
			evts.push_back(evt);
		}
	}
	
	
	double bkgHitTime = 0.0;
	evt.ch = 1;
	while(true) {
		bkgHitTime += -log(pop(sfmt))/BKGRATEA;
		if(bkgHitTime > LENGTH) {
			break;
		}
		evt.realtime = bkgHitTime;
		evts.push_back(evt);
	}
	bkgHitTime = 0.0;
	evt.ch = 2;
	while(true) {
		bkgHitTime += -log(pop(sfmt))/BKGRATEB;
		if(bkgHitTime > LENGTH) {
			break;
		}
		evt.realtime = bkgHitTime;
		evts.push_back(evt);
	}
	
	bkgHitTime = 0.0;
	while(true) {
		bkgHitTime += -log(pop(sfmt))/COINCBKGRATE;
		if(bkgHitTime > LENGTH) {
			break;
		}
		uint32_t index = genIndex(phs_bkg, pop64(sfmt));
		while(index == phs_bkg->nEntries - 1) {
			index = genIndex(phs_bkg, pop64(sfmt));
		}
		int numA = phsVect_bkg[index].a;
		int numB = phsVect_bkg[index].b;
		//printf("%f - %d:%d\n", bkgHitTime, numA, numB);
		int srcChan = pop(sfmt) < PMTAPROB ? 1 : 2;
		histGen* ch1 = srcChan == 1 ? ch1_1_bkg : ch1_2_bkg;
		histGen* ch2 = srcChan == 1 ? ch2_1_bkg : ch2_2_bkg;
		if(srcChan == 1) {numA--; evts.push_back(event {bkgHitTime, 1});} //Automatically put hit at t=0
		if(srcChan == 2) {numB--; evts.push_back(event {bkgHitTime, 2});} //Automatically put hit at t=0
		evt.realtime = 0.0;
		evt.ch = 0;
		for(int i = 0; i < numA; i++) {
			uint32_t index = genIndex(ch1, pop64(sfmt));
			while(index == ch1->nEntries - 1) {
				index = genIndex(ch1, pop64(sfmt));
			}
			evt.ch = 1;
			evt.realtime = bkgHitTime + index * CLOCKTONS * NANOSECOND;// * 0.01;
			evts.push_back(evt);
		}
		for(int i = 0; i < numB; i++) {
			uint32_t index = genIndex(ch2, pop64(sfmt));
			while(index == ch2->nEntries - 1) {
				index = genIndex(ch2, pop64(sfmt));
			}
			evt.ch = 2;
			evt.realtime = bkgHitTime + index * CLOCKTONS * NANOSECOND;// * 0.01;
			evts.push_back(evt);
		}
	}
	
	std::sort(evts.begin(), evts.end(), [](event x, event y)->bool{return (x.realtime < y.realtime);});
	
	return evts;
}
