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

#ifdef BRII
	#define LOCATION "/N/u/nbcallah/BigRed2/TauPileup/BRII/2016_2017_Waveform_Data"
#endif

#ifdef KRST
	#define LOCATION "/N/u/nbcallah/Karst/TauPileup/BRII/2016_2017_Waveform_Data"
#endif

#ifdef LOCAL
	#define LOCATION "/Users/nate/Ultra_Cold_Neutrons/2016-2017/NatesAnalyzer-Clone/synthetic_data/BRII/2016_2017_Waveform_Data"
#endif

//#define LENGTH 100.0

#define BKGRATEA 560
#define BKGRATEB 545

#define COINCBKGRATE 4.89

#define PMTAPROB 0.584975

#define NANOSECOND .000000001
#define CLOCKTONS 0.8

//#define COINCBKG 0.106

//#define COINCBKG 0.117058656786906 //For 50_1000_8
//#define COINCBKG 0.0458956
//#define COINCBKG 0.0

#define COINCWINDOW 50.0
//#define PESUMWINDOW 500.0
//#define PESUM 8

//#define COINCWINDOW 50.0
//#define PESUMWINDOW 1000.0
//#define PESUM 8

#define TAU_N 877.6

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

struct dip {
	double tau;
	double time;
	double cts;
};

struct measurement {
	double val;
	double err;
	measurement operator/(const measurement& rhs) {
		return measurement {val/rhs.val, val/rhs.val*sqrt(pow(err/val, 2.0)+pow(rhs.err/rhs.val, 2.0))};
	}
	measurement operator*(const measurement& rhs) {
		return measurement {val*rhs.val, val*rhs.val*sqrt(pow(err/val, 2.0)+pow(rhs.err/rhs.val, 2.0))};
	}
	measurement operator+(const measurement& rhs) {
		return measurement {val+rhs.val, sqrt(err*err+rhs.err*rhs.err)};
	}
	measurement operator-(const measurement& rhs) {
		return measurement {val-rhs.val, sqrt(err*err+rhs.err*rhs.err)};
	}
};

#ifdef DEBUG
	template class std::vector<double>;
	template class std::vector<event>;
#endif

//struct iat { //Interarrival time
//	double t;
//	double dt;
//	double dtDiff;
//};

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

//double coincbkg[4][12] = {	{1.22210389613692,1.30066258889395,1.34907820243417,1.38275826497236,1.40687940813643,1.42791495843683,1.46129400702337, 0.0, 0.0, 0.0, 0.0, 0.0},
//{0.282826323932059,0.313503829545169,0.333364064754208,0.347147436915004,0.357783022240608,0.366888149486427,0.38136399369093, 0.0, 0.0, 0.0, 0.0, 0.0},
//{0.0966203728098595,0.106320118034727,0.11259287261949,0.117058656786906,0.120474121319817,0.123423919623387,0.128081970694195, 0.0941322,0.0989980,0.101074,0.103058,0.104825},
//{0.0549676624345168,0.0588065508982811,0.0612549496510915,0.0629702802928673,0.0643182657562066,0.0654809907683621,0.0672582478010832, 0.0, 0.0, 0.0, 0.0, 0.0}
//};

int nph[7] = {4, 5, 6, 7, 8, 9, 10};
double window[8] = {200, 250, 300, 350, 400, 500, 750, 1000};
//double coincbkg[7][8] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//						 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//						 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//						 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//						 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//						 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//						 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
double coincbkg[7][8] ={{1.198134730539,1.223737003058,1.241981424149,1.258849056604,1.271412698413,1.297482200647,1.346610738255,1.381534482759,},
{0.562505617978,0.576792507205,0.589269513991,0.600774774775,0.612515290520,0.629973228346,0.663006622517,0.689825862069,},
{0.274744505495,0.283271054494,0.290568627451,0.297097253155,0.303093181818,0.314630503145,0.333898998331,0.346852686308,},
{0.150076519130,0.154173025048,0.159034578696,0.162370941558,0.165102352456,0.171304366438,0.183113501144,0.190070308789,},
{0.094266258247,0.096429260063,0.098691662556,0.101049002273,0.102924620530,0.105916070956,0.112568936410,0.117171060340,},
{0.067602568435,0.069441416421,0.070615710503,0.071992262012,0.073076909024,0.075151606237,0.078647336348,0.081723391216,},
{0.054136419001,0.055077103125,0.056018764879,0.056862402274,0.057636075493,0.058907377411,0.061307126437,0.063176563487,},
};
//double window[12] = {250, 500, 750, 1000, 1250, 1500, 2000, 200, 300, 350, 400, 450};

template<typename T, typename F>
void loadTextHistogram(histGen* &gen, std::string fName, std::vector<T> &vec, F f);

template<typename T, typename F>
void fillRandomArray(sfmt_t* sfmt, T* ary, int num, F f);

std::vector<event> createSynthRun(sfmt_t* sfmt, double rate, std::function<double(sfmt_t*)> pop, std::function<uint64_t(sfmt_t*)>pop64, long &numCreated);

double getDeadTimeCounts(std::vector<double> dtHist, std::vector<double> coincHist) {
	double deadTimeCounts = 0.0;
	auto dtIt = dtHist.begin();
	for(auto cIt = coincHist.begin(); cIt < coincHist.end(); cIt++) {
//		deadTimeCounts += (*cIt/(1.0-*dtIt/(1.0/10.0)) - *cIt);
		//deadTimeCounts += (*cIt/((1.0/10.0 - *dtIt)/(1.0/10.0)) - *cIt);
		printf("%.12e\n", *dtIt/(*cIt));
		dtIt++;
	}
	return deadTimeCounts;
}

//void getCoincs2(std::vector<event> data, std::vector<double> &coincs, std::vector<double> &dtHist, std::vector<double> &coincHist) {
//	for(int i = 0; i < 10*ceil(data.back().realtime); i++) {
//		dtHist.push_back(0.0);
//		coincHist.push_back(0.0);
//	}
//	std::vector<iat> ttne;
//	for(auto it = data.begin(); it < data.end(); it++) {
//		if(it->ch != 1 && it->ch != 2) { continue; }
//
//		auto nextDiff = std::find_if(it+1, data.end(), [it](event x){return (x.ch == 1 && it->ch == 2) || (x.ch == 2 && it->ch == 1);});
//		auto next = std::find_if(it+1, data.end(), [](event x){return (x.ch == 1) || (x.ch == 2);});
//
//		//if(next == data.end()) { break; }
//		
//		iat evt;
//		evt.t = it->realtime;
//		evt.dt =  next < data.end() ? next->realtime - it->realtime : std::numeric_limits<double>::infinity();
//		//evt.dt =  next->realtime - it->realtime;
//		evt.dtDiff = nextDiff < data.end() ? nextDiff->realtime - it->realtime : std::numeric_limits<double>::infinity();
//		
//		ttne.push_back(evt);
//	}
//	for(auto it = ttne.begin(); it < ttne.end(); it++) {
//		if(it->dtDiff < COINCWINDOW*NANOSECOND) { //Have a coincidence!
//			auto cIt = std::find_if(it+1, ttne.end(), [](iat x){return x.dt > PESUMWINDOW*NANOSECOND;});
//			int numPE = std::distance(it, cIt) + 1;
//			if(numPE >= PESUM) {
//				dtHist[(int)(10*it->t)] += ((cIt)->t - it->t);
//				coincHist[(int)(10*it->t)] += 1;
//				it = cIt;
//			}
//		}
//	}
//}
void getCoincs(std::vector<event> data, std::vector<double> &coincs, std::vector<double> &dtHist, std::vector<double> &coincHist, int inumPH, int ipesumWindow) {
	int i;
	int cur;
	int tailIt;

	for(i = 0; i < 10*ceil(data.back().realtime); i++) {
		dtHist.push_back(0.0);
		coincHist.push_back(0.0);
	}
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
//					if((data.at(tailIt).realtime - prevEvt.realtime) > PESUMWINDOW*NANOSECOND) {
					if((data.at(tailIt).realtime - prevEvt.realtime) > window[ipesumWindow]*NANOSECOND) {
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
				//if(ch1PESum >= peSum && ch2PESum >= peSum) { //Found Neutron!
//				if(ch1PESum + ch2PESum >= PESUM) {
				if(ch1PESum + ch2PESum >= nph[inumPH]) {
					//coincs.push_back(data.at(i).realtime);
					dtHist[(int)(10*data[i].realtime)] += (prevEvt.realtime - data.at(i).realtime) + window[ipesumWindow]*NANOSECOND;
					//printf("Boo!! %d, %f, %f\n", (int)data[i].realtime, (prevEvt.realtime - data.at(i).realtime), dtHist[(int)data[i].realtime]);
					coincHist[(int)(10*data[i].realtime)] += 1;
					i = tailIt-1; //Only impose deadtime on counted neutrons; Otherwise continue looking.
					cur = tailIt-1;
					break;
				}
				break;
			}
//			break;
		}
	}
}

measurement normNByDip(std::vector<event> run, int inumPH, int ipesumWindow) {
	measurement num;
	
	std::vector<double> coincs;
	std::vector<double> dtHist;
	std::vector<double> coincHist;
	getCoincs(run, coincs, dtHist, coincHist, inumPH, ipesumWindow);
	
	//double dagCts = std::accumulate(coincs.begin(), coincs.end(), 0.0, [startTime, stopTime](double sum, double time){return time > startTime && time < stopTime ? sum + 1.0 : sum;});
	
	double dagCts = std::accumulate(coincHist.begin(), coincHist.end(), 0.0, [](double sum, double cts){return sum + cts;});
	//printf("numCts: %f\n", dagCts / 120.0);

	double deadTimeCounts = getDeadTimeCounts(dtHist, coincHist);
	//printf("dt: %f\n", deadTimeCounts);
	
	double bkgrate = coincbkg[inumPH][ipesumWindow];
	
	//num = {dagCts+deadTimeCounts-bkgrate*(coincHist.size()/10.0), sqrt(dagCts+deadTimeCounts+bkgrate*(coincHist.size()/10.0))};
	num = {dagCts-bkgrate*(coincHist.size()/10.0), sqrt(dagCts+bkgrate*(coincHist.size()/10.0))};
	//num = {dagCts, sqrt(dagCts)};
	
	return num;
}

int main(int argc, char** argv) {
	int c;
	
//	double tau = -1.0;
//	double hits = 0;
	double rate = 0.0;
	double precision = 0.0;
	int inumPH = -1;
	int ipesumWindow = -1;
//	double start = -1.0;
//	double length = -1.0;
//	char* pattern;
//	int startNum = -1;
//	int endNum = -1;
//	char fName[256];
	
	while (1)
	{
		static struct option long_options[] =
		{
			/* These options don't set a flag.
			 We distinguish them by their indices. */
			{"rate", required_argument, 0, 'r'},
			{"precision", required_argument, 0, 'p'},
			{"numpe", required_argument, 0, 'n'},
			{"window", required_argument, 0, 'w'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "r:p:n:w", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1) { break; }

		switch (c)
		{
			case 0:
				break;
			
			case 'r':
				rate = atof(optarg);
				break;
				
			case 'p':
				precision = atof(optarg);
				break;
				
			case 'n':
				inumPH = atoi(optarg);
				break;
				
			case 'w':
				ipesumWindow = atoi(optarg);
				break;

			default:
				abort();
		}
	}

	if(rate == 0.0 || precision == 0.0 || inumPH < 0 || ipesumWindow < 0) {
		fprintf(stderr, "Error! Usage: ./randomGenerator --hits=num_hits_100s --precision=final_rel_precision --numpe=numpe --window=integration_window\n");
		exit(1);
	}
	
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
	
//	histGen* gen;
	
	std::vector<unsigned long> dump;
//	std::vector<ph> phsVect;
//	std::vector<ph> phsVect_bkg;
//	
//	histGen* gen_ch1_1;
//	histGen* gen_ch1_2;
//	histGen* gen_ch2_1;
//	histGen* gen_ch2_2;
//	histGen* gen_phs;
//	histGen* gen_ch1_1_bkg;
//	histGen* gen_ch1_2_bkg;
//	histGen* gen_ch2_1_bkg;
//	histGen* gen_ch2_2_bkg;
//	histGen* gen_phs_bkg;
	
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
	
	measurement cts= {0.0, 0.0};
//	measurement eff;
	
	long createdCts;
	std::vector<double> eff;
	double mean = 0.0;
	double stdDev = 0.0;
	double sumCts = 0.0;
	//measurement sumCreated = {0.0, 0.0};
	//measurement sum = {0.0, 0.0};
	
	int ltNum = 0;
	
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
	//
//	printf("{");
//	for(inumPH = 0; inumPH < 7; inumPH++) {
//		printf("{");
//		for(ipesumWindow = 0; ipesumWindow < 8; ipesumWindow++) {
			sumCts = 0.0;
			ltNum = 0;
			eff.clear();
			do {
		//		measurement cts= {0.0, 0.0};
		//		long createdCts;
		//		measurement sumCreated = {0.0, 0.0};
		//		measurement sum = {0.0, 0.0};
				ltNum += 1;

				std::vector<event> vec = createSynthRun(&sfmt, rate, pop, pop64, createdCts);

				cts = normNByDip(vec, inumPH, ipesumWindow);
				sumCts += cts.val;

				eff.push_back(cts.val/((double)createdCts));

				mean = std::accumulate(eff.begin(), eff.end(), 0.0, [](double sum, double e){return sum+e;})/eff.size();
				stdDev = sqrt(std::accumulate(eff.begin(), eff.end(), 0.0, [mean](double sum, double e){return sum+pow((e-mean), 2);})/(eff.size()-1));

				stdDev = ltNum >= 10 ? stdDev : std::numeric_limits<double>::infinity();

		//		
		//		
		//		sumCreated = sumCreated + measurement{(double)createdCts, sqrt((double)createdCts)};
		//		
		//		cts = normNByDip(vec, inumPH, ipesumWindow);
		//		sum = sum + cts;
		//		
		//		eff = sum / sumCreated;
		//		eff.err = sqrt(eff.val*(1-eff.val)/sumCreated.val);
		//	} while(sumCts < 400000);
			} while((stdDev/sqrt(eff.size()))/mean > precision);
			printf("Final : %f,%.12f,%.12f,%d,%d,%f\n", rate, mean, (stdDev/sqrt(eff.size())), ltNum, nph[inumPH], window[ipesumWindow]);
			printf("Sum: %f,%d,%d,%f\n", sumCts, ltNum, nph[inumPH], window[ipesumWindow]);
//			printf("%.12f,", sumCts/(ltNum*1000.0));
//		}
//		printf("},\n");
//	}
//	printf("}\n");
	
//	printf("Final : %f,%.12f,%.12f,%d,%d,%f\n", rate, mean, (stdDev/sqrt(eff.size())), ltNum, nph[inumPH], window[ipesumWindow]);
//	printf("Sum: %f,%d,%d,%f\n", sumCts, ltNum, nph[inumPH], window[ipesumWindow]);

//	printf("Bkg: %f,%d,%f\n", sumCts/(ltNum*1000.0), nph[inumPH], window[ipesumWindow]);
		//fflush(stdout);
//	} while(sqrt(sqErrSum)/ltNum > sqrt((double)world_size)*precision);
//	} while(err > sqrt((double)world_size)*precision);
//	} while(ltNum < 19999);
	
//	printf("%f\n", sumShort); printf("%f\n", sumLong); printf("%d\n", ltNum);
	
//	printf("Final : %f,%f\n", ltSum/ltNum, sqrt(sqErrSum)/ltNum);
//	printf("Final : %f,%f\n", lt, err);
//	printf("Sums: %f, %f\n", sumShort, sumLong);
//	printf("num: %d\n", ltNum);
	
	//printf("Num Ints: %llu\n", numInts);
	
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
	
	MPI_Finalize();
	
	
	
	/*for(i = startNum; i <= endNum; i++) {
		sprintf(fName, pattern, i);
		if(shortRun) {
			createSynthRun(&sfmt, fName, tau, hits, start, length, hitTimes, phsAnum, phsBnum);
		}
		else {
			createSynthRun(&sfmt, fName, tau, lround(hits*exp(-HOLD_T/TAU_N)), start+HOLD_T, length+HOLD_T, hitTimes, phsAnum, phsBnum);
		}
		shortRun = !shortRun;
	}*/
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
	
	double length = 10.0;
	
	std::vector<event> evts;
	
	event evt;
	
	numCreated = 0;
	
	double coincHitTime = 0.0;
	while(true) {
		coincHitTime += -log(pop(sfmt))/rate;
		if(coincHitTime > length) {
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
		if(bkgHitTime > length) {
			break;
		}
		evt.realtime = bkgHitTime;
		evts.push_back(evt);
	}
	bkgHitTime = 0.0;
	evt.ch = 2;
	while(true) {
		bkgHitTime += -log(pop(sfmt))/BKGRATEB;
		if(bkgHitTime > length) {
			break;
		}
		evt.realtime = bkgHitTime;
		evts.push_back(evt);
	}
	
	bkgHitTime = 0.0;
	while(true) {
		bkgHitTime += -log(pop(sfmt))/COINCBKGRATE;
		if(bkgHitTime > length) {
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
