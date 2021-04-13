/*
SOMA CLUSTER version with linear adaptation of PRT
SOMA-CLP
	-	CEC 2021
	-	multithread
	-	advanced log system
	-   log only error values
	-   use mersen twister with seed specified by CEC21

	#parameters:
		-	number of test function
		-	dimension size
		-	population size
		-	FES size
		-	run size (runMax)
		-   CEC21 C parameter       - which type to use, <0-7>
                                    - 0 :   000   :   Basic
                                    - 1 :   100   :   Bias
                                    - 2 :   010   :   Shift
                                    - 3 :   001   :   Rotation
                                    - 4 :   110   :   Bias & Shift
                                    - 5 :   101   :   Bias & Rotation
                                    - 6 :   011   :   Shift & Rotation
                                    - 7 :   111   :   Bias, Shift & Rotation


	#optional parameters (#-name #value):
		-prt #
		-step #
		-path #

	#descriptions:
		Two switching strategies (All-To-All, All-To-One)
		Save all points from All-To-All
		For All-To-ONE, sort by OFV, cluster & kill, rank select from rest ONE
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <thread>
//#include <Windows.h>
#include <random>
#include <filesystem>
#include <algorithm>
#include <chrono>
#include <iostream>

#include <cstring>

#include "kmeans.h"

//CEC 2021 by definition
//------------------------------------------------------
void cec21_bias_shift_rot_func(double *, double *,int,int,int);
void cec21_bias_shift_func(double *, double *,int,int,int);
void cec21_bias_rot_func(double *, double *,int,int,int);
void cec21_shift_rot_func(double *, double *,int,int,int);
void cec21_rot_func(double *, double *,int,int,int);
void cec21_shift_func(double *, double *,int,int,int);
void cec21_bias_func(double *, double *,int,int,int);
void cec21_basic_func(double *, double *,int,int,int);

//CEC21 Binary parameter switch
int CEC_C = 0;

//minimal values of CEC2021 functions
//used only for certain conbination
double cecMin[] = {100, 1100,700,1900,1700,1600,2100,2200,2400,2500};
//mask mins values of CEC based on binary parameter
double cecMask[] = {0,1,0,0,1,1,0,1};

//load CEC Rand Seed MAX ITEMS
#define CEC_RAND_SEED_MAX 1000
//load CEC Rand Seed
double Rand_Seeds[CEC_RAND_SEED_MAX];

/*vars MUST be changed in benchmark file

	static thread_local double *OShift, *M, *y, *z, *x_bound;
	static thread_local int ini_flag, n_flag, func_flag, *SS;

*/
//double *OShift, *M, *y, *z, *x_bound;
//int ini_flag = 0, n_flag, func_flag, *SS;
//------------------------------------------------------



/*
Global variables
-------------------------------------------------------
*/
//number of selected test
int test = 0;
//sizes (dimension,iteration,population)
int sizeDim = 0;
int sizeIt = 0;
int sizePop = 0;
int maxFES = 0;
int FES = 0;
//bounds (defined by CEC 2020)
double static limMin = -100;
double static limMax = 100;
//number of algoritm runs (also number of threads)
int runMax = 0;
//FILE INFO pointer
FILE* fileInfo = NULL;
//FILE INFO name
const char* FILE_INFO_NAME = "SOMA_CLP_CEC21";
//DIR name (for runs)
char DIR_NAME[200];

//SOMA control parameters
double prt = 0.3; /* {0.1, 0.3, 0.5, 0.9} */
double step = 0.11; /* {0.11, 0.31} */
double path = 2;

//SOMA special control params for All-To-Random
double stepA = 0.33;
double pathA = 3;
double prtA = 0.5;
int maxLeaders = 10;	//should be sizePop/10 

//Shared memory for results from runs
double* gBestCFs = NULL;
double** gBests = NULL;

// generator
static thread_local std::mt19937* generator = nullptr;

/*
Common Functions
-------------------------------------------------------
*/
//Main function
void run(int runID);
//Argument parsing
void argumentParse(int argc, char **argv);
//Create and fill INFO FILE
void createInfoFile();
//Random generator (uniform - rovnomerne)
int intRand(const int & min, const int & max);
double doubleRand(const double & min, const double & max);
//Fit function ([input],[dest])
void fit(double* input, double* destination);

/*
Special Functions
-------------------------------------------------------
*/
Point leaderRoulette(std::vector<Point>*);
Point leaderRank(std::vector<Point>*);

/*
	MAIN ENTRY
*/
int main(int argc, char **argv) {
	
	auto startT = std::chrono::high_resolution_clock::now();

	//parse arguments
	argumentParse(argc, argv);

	//print name of algorithm
	printf("--- %s ---\n\n", FILE_INFO_NAME);

	//create and fill INFO FILE
	createInfoFile();

    //preload random seed data for CEC21
    FILE* tempFile = fopen("input_data/Rand_Seeds.txt", "r");
    if(tempFile == NULL) {
        printf("\t ERROR file not found input_data/Rand_Seeds.txt");
        return 1;
    }
    int count = 0;
    char tempBuffer[1000];
    while(fgets(tempBuffer,1000, tempFile) != NULL && count < CEC_RAND_SEED_MAX) {
        sscanf(tempBuffer, "%lf", &Rand_Seeds[count++]);
    }

	//create threads
	std::thread* holder = new std::thread[runMax];

	//create global arrays for results
	gBestCFs = new double[runMax];
	gBests = new double*[runMax];
	for (int i = 0; i < runMax; i++) {
		gBests[i] = new double[sizeDim];
	}

	//create dir for runs
	sprintf(DIR_NAME, "%s_d%d_t%d_c%d", FILE_INFO_NAME, sizeDim, test, CEC_C);
	std::filesystem::remove_all(DIR_NAME);
	std::filesystem::create_directory(DIR_NAME);

	//time spend on algoritm
	clock_t timeStart, timeEnd;
	timeStart = clock();
	auto start = std::chrono::high_resolution_clock::now();


	//run threads
	for (int i = 0; i < runMax; i++) {
		//TODO odkomentovat
		holder[i] = std::thread(run, i);
	}

	//wait for end of all threads
	for (int i = 0; i < runMax; i++) {
		//TODO odkomentovat
		holder[i].join();
	}

	//log time
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	std::cout << "Time taken by function: "
		<< duration.count() << " microseconds" << std::endl;
	timeEnd = clock();
	double time = (double)(timeEnd - timeStart) / CLOCKS_PER_SEC;
	//fprintf(fileInfo, "time = %f\n", time);
	fprintf(fileInfo, "time = %lf\n", duration.count());
	fprintf(fileInfo, "(*---------------*)\n");

	//statistical evaluation and final log INFO FILE
	//still need to find index of this minimum :(
	//double* min = std::min_element(gBestCFs, gBestCFs + runMax);
	double* max = std::max_element(gBestCFs, gBestCFs + runMax);
	double min = gBestCFs[0];
	int minIndex = 0;
	for (int i = 1; i < runMax; i++) {
		if (gBestCFs[i] < min) {
			min = gBestCFs[i];
			minIndex = i;
		}
	}

	double sum = 0.0, standardDeviation = 0.0, mean = 0.0;

	for (int i = 0; i < runMax; ++i) {
		sum += gBestCFs[i];
	}
	mean = sum / runMax;
	for (int i = 0; i < runMax; ++i) {
		standardDeviation += pow(gBestCFs[i] - mean, 2);
	}
	standardDeviation = sqrt(standardDeviation / runMax);

	fprintf(fileInfo, "min = %.8e\n", min);
	fprintf(fileInfo, "max = %.8e\n", *max);
	fprintf(fileInfo, "std = %.8e\n", standardDeviation);
	fprintf(fileInfo, "mean = %.8e\n", mean);
	fprintf(fileInfo, "(*---------------*)\n");

	//find best and print to INFO FILE
	fprintf(fileInfo, "bestCF = %.8e\n", min);
	fprintf(fileInfo, "best = {");
	for (int i = 0; i < sizeDim - 1; i++) {
		fprintf(fileInfo, "%.8e,", gBests[minIndex][i]);
	}
	fprintf(fileInfo, "%.8e}\n", gBests[minIndex][sizeDim - 1]);
	fprintf(fileInfo, "(*---------------*)\n");

	//log all final runs
	fprintf(fileInfo, "lastRuns = {");
	for (int i = 0; i < runMax - 1; i++) {
		fprintf(fileInfo, "%.8e,", gBestCFs[i]);
	}
	fprintf(fileInfo, "%.8e}\n", gBestCFs[runMax - 1]);
	fprintf(fileInfo, "(*---------------*)\n");

	//release threads and other dynamic arrays
	delete[] holder;
	delete[] gBestCFs;
	for (int i = 0; i < runMax; i++) {
		delete[] gBests[i];
	}
	delete[] gBests;

	//close INFO FILL
	fclose(fileInfo);

	return 0;
}

void argumentParse(int argc, char **argv) {
	/*
	required arguments(!number & position!):
	number of test
	dimension size
	population size
	FES size
	run size (runMax)
	*/

	if (argc < 7) {
		printf("BAD ARGUMENTS (number of arguments < 6)");
		exit(1);
	}
	test = strtol(argv[1], NULL, 10);
	if (test < 1 || test>10) {
		printf("BAD ARGUMENTS (selected test = %s)", argv[1]);
		exit(1);
	}
	sizeDim = strtol(argv[2], NULL, 10);
	if (sizeDim != 5 && sizeDim != 10 && sizeDim != 15 && sizeDim != 20) {
		printf("BAD ARGUMENTS (dimension size = %s)", argv[2]);
		exit(1);
	}
	sizePop = strtol(argv[3], NULL, 10);
	if (sizePop < 0) {
		printf("BAD ARGUMENTS (population size = %s)", argv[3]);
		exit(1);
	}
	maxFES = strtol(argv[4], NULL, 10);
	if (maxFES < 0) {
		printf("BAD ARGUMENTS (number of FES = %s)", argv[4]);
		exit(1);
	}
	runMax = strtol(argv[5], NULL, 10);
	if (runMax < 1) {
		printf("BAD ARGUMENTS (number of runs of program = %s)", argv[5]);
		exit(1);
	}
    CEC_C = strtol(argv[6], NULL, 10);
	if (CEC_C < 0 || CEC_C > 7) {
		printf("BAD ARGUMENTS (CEC binary switch = %s)", argv[6]);
		exit(1);
	}

	sizeIt = maxFES / sizePop;

	for (int i = 8; i < argc; i += 2) {

		if (strcmp(argv[i - 1], "-prt") == 0) {
			prt = strtod(argv[i], NULL);
			printf("\n  -prt: %f", prt);
		}

		if (strcmp(argv[i - 1], "-path") == 0) {
			path = strtod(argv[i], NULL);
			printf("\n  -path: %f", path);
		}

		if (strcmp(argv[i - 1], "-step") == 0) {
			step = strtod(argv[i], NULL);
			printf("\n  -step: %f", step);
		}

	}
}

void createInfoFile() {
	//create full info name file
	//print also on stdout
	char buffer[200];
	sprintf(buffer, "info_%s_d%d_t%d_c%d.txt", FILE_INFO_NAME, sizeDim, test, CEC_C);
	fileInfo = stdout;
	FILE* tempFile = fopen(buffer, "w");
	while (1) {
		fprintf(fileInfo, "testFunction = %d\n", test);
        fprintf(fileInfo, "CEC_C = %d\n", CEC_C);
		fprintf(fileInfo, "dim = %d\n", sizeDim);
		fprintf(fileInfo, "pop = %d\n", sizePop);
		fprintf(fileInfo, "FE = %d\n", maxFES);
		fprintf(fileInfo, "It = %d\n", sizeIt);
		fprintf(fileInfo, "runs = %d\n", runMax);
		fprintf(fileInfo, "(*---------------*)\n");
		fprintf(fileInfo, "prt = %.8e\n", prt);
		fprintf(fileInfo, "step = %.8e\n", step);
		fprintf(fileInfo, "path = %.8e\n", path);
		fprintf(fileInfo, "(*---------------*)\n");
		if (fileInfo == stdout) {
			fileInfo = tempFile;
		}
		else {
			break;
		}
	}
}

/* Thread-safe function that returns a random number between min and max (inclusive).
This function takes ~142% the time that calling rand() would take. For this extra
cost you get a better uniform distribution and thread-safety. */
int intRand(const int & min, const int & max) {
	static thread_local std::mt19937* generator = nullptr;
	if (!generator) generator = new std::mt19937(clock() +
		std::hash<std::thread::id>()(std::this_thread::get_id()));
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(*generator);
}

double doubleRand(const double & min, const double & max) {
	//static thread_local std::mt19937* generator = nullptr;
	//if (!generator) generator = new std::mt19937(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
   // if (!generator) generator = new std::mt19937(sizeDim/10*test*runMax + std::this_thread::get_id())-runMax;
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(*generator);
}

//change due to new CEC21
void fit(double* input, double* destination) {
    /*
    void cec21_bias_shift_rot_func(double *, double *,int,int,int);
    void cec21_bias_shift_func(double *, double *,int,int,int);
    void cec21_bias_rot_func(double *, double *,int,int,int);
    void cec21_shift_rot_func(double *, double *,int,int,int);
    void cec21_rot_func(double *, double *,int,int,int);
    void cec21_shift_func(double *, double *,int,int,int);
    void cec21_bias_func(double *, double *,int,int,int);
    void cec21_basic_func(double *, double *,int,int,int);
    */
    switch(CEC_C) {
        case 0: {
            cec21_basic_func(input, destination, sizeDim, 1, test);
        };
        break;
        case 1: {
            cec21_bias_func(input, destination, sizeDim, 1, test);
        };
        break;
        case 2: {
            cec21_shift_func(input, destination, sizeDim, 1, test);
        };
        break;
        case 3: {
            cec21_rot_func(input, destination, sizeDim, 1, test);
        };
        break;
        case 4: {
            cec21_bias_shift_func(input, destination, sizeDim, 1, test);
        };
        break;
        case 5: {
            cec21_bias_rot_func(input, destination, sizeDim, 1, test);
        };
        break;
        case 6: {
            cec21_shift_rot_func(input, destination, sizeDim, 1, test);
        };
        break;
        case 7: {
            cec21_bias_shift_rot_func(input, destination, sizeDim, 1, test);
        };
        break;
    }

    //update min value and mask it if needed
    *destination -= (cecMin[test - 1] * cecMask[CEC_C]);

}

//custom function to sort by CLUSTERID
bool sortByClusterID(Point i, Point j) {
	return (i.getCluster() < j.getCluster());
}

//custom function to sort by CF
bool sortByCF(Point i, Point j) {
	return (i.getCluster() < j.getCluster());
}

//rank selection of leader
//todo fix rank
Point leaderRoulette(std::vector<Point>* leaders) {

	double max = 0;
	for (std::vector<Point>::iterator i = leaders->begin(); i != leaders->end(); i++) {
		//max += i->CF;
		max += 1/i->getCF();
	}

	double rand = doubleRand(0, max);
	std::vector<Point>::iterator selected;

	for (std::vector<Point>::iterator i = leaders->begin(); i != leaders->end(); i++) {
		if (rand <  (1/i->getCF())) {
			selected = i;
			//*leader = i;
			return *i;
			//break;
		}
		rand -= (1/i->getCF());
	}
}

Point leaderRank(std::vector<Point>* leaders) {

	int max = ((leaders->size() / 2) * (2 + leaders->size() - 1));

	//int rand = intRand(0, max);
    int rand = (int)doubleRand(0, max);
	std::vector<Point>::iterator selected;
	int start = leaders->size();
	//int start = 1;
	for (std::vector<Point>::iterator i = leaders->begin(); i != leaders->end(); i++) {
		//i->CF = start--;
		i->setCF(start--);
	}

	start = 0;
	for (std::vector<Point>::iterator i = leaders->begin(); i != leaders->end(); i++) {
		start += i->getCF();
		if (rand < start) {
			return *i;
		}
	}

	return leaders->at(leaders->size() - 1);
}

//Main algorithm
void run(int runID) {


    //init random generator
    int seed = (sizeDim/10*test*runMax + runID + 1)-runMax;
    seed %= CEC_RAND_SEED_MAX;
    if (!generator) generator = new std::mt19937(Rand_Seeds[seed]);


	//create file in DIR for this run
	FILE* logFile = NULL;
	{
		char buffer[200];
		sprintf(buffer, "%s/%d.txt", DIR_NAME, runID);
		logFile = fopen(buffer, "w");
		fprintf(logFile, "{");
	}

	//FES counter
	int FES = 0;

	//gBest & gBestCF
	double gBestCF = 0;
	double* gBest = new double[sizeDim];

	//initial population
	double** pop = new double*[sizePop];
	for (int i = 0; i < sizePop; i++) {
		pop[i] = new double[sizeDim];
		for (int j = 0; j < sizeDim; j++) {
			pop[i][j] = doubleRand(limMin, limMax);
		}
	}
	double* popCF = new double[sizePop];
	fit(pop[0], &popCF[0]);
	gBestCF = popCF[0];
	for (int i = 0; i < sizeDim; i++) {
		gBest[i] = pop[0][i];
	}
	FES++;
	//log first entry
	fprintf(logFile, "{%d,%.8e},", FES, gBestCF);
	for (int i = 1; i < sizePop; i++) {
		fit(pop[i], &popCF[i]);
		FES++;
		if (popCF[i] < gBestCF) {
			gBestCF = popCF[i];
			for (int j = 0; j < sizeDim; j++) {
				gBest[j] = pop[i][j];
			}
			//log if gBest is changed
			fprintf(logFile, "{%d,%.8e},", FES, gBestCF);
		}
	}

	//for III. ALL-TO-LEADER
	double** popTemp = new double*[sizePop];
	for (int i = 0; i < sizePop; i++) {
		popTemp[i] = new double[sizeDim];
	}
	double* popCFTemp = new double[sizePop];

	double tempCF = 0;
	double* temp = new double[sizeDim];

	//prepare prtVector
	double* prtVector = new double[sizeDim];

	//prepare field of points
	std::vector<Point> points;
	int pointID = 0;

	//main run
	char on = 1;
	if (FES >= maxFES) {
		on = 0;
	}
	while (on) {
		
		//clear all points
		points.clear();
		pointID = 0;


		//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
		//WWWWWWWWW		I. ALL-TO-RANDOM
		//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
		for (int i = 0; i < sizePop; i++) {
			//copy individual 
			popCFTemp[i] = popCF[i];
			tempCF = popCF[i];
			for (int j = 0; j < sizeDim; j++) {
				popTemp[i][j] = pop[i][j];
				temp[j] = pop[i][j];
			}
			
			//select leader
			int leader = -1;
			do {
				//leader = intRand(0, sizePop - 1);
                leader = (int)doubleRand(0, sizePop - 1);
			} while (leader == i);

			//move over path
			for (double t = 0; t <= pathA; t += stepA) {

				//generate PRT vector for current leader
				for (int j = 0; j < sizeDim; j++) {
                    prtA = 0.08 + 0.9 * ((double)FES/(double)maxFES);
					if (doubleRand(0, 1) > prtA) {
						prtVector[j] = 0;
					}
					else {
						prtVector[j] = 1;
					}
				}

				// point of interest

				for (int j = 0; j < sizeDim; j++) {
					temp[j] = pop[i][j] + (pop[leader][j] - pop[i][j])*t*prtVector[j];

					//simple random strategy
					if (temp[j] > limMax) {
						temp[j] = doubleRand(limMin, limMax);
					}
					if (temp[j] < limMin) {
						temp[j] = doubleRand(limMin, limMax);
					}

				}
				//evaluate new substep (actualy fitness popTemp)
				fit(temp, &tempCF);
				FES++;
				if (tempCF < gBestCF) {
					gBestCF = tempCF;
					for (int j = 0; j < sizeDim; j++) {
						gBest[j] = temp[j];
					}
					//log if gBest is changed
					fprintf(logFile, "{%d,%.8e},", FES, gBestCF);
				}
				if (tempCF < popCFTemp[i]) {
					popCFTemp[i] = tempCF;
					for (int j = 0; j < sizeDim; j++) {
						popTemp[i][j] = temp[j];
					}
				}

				//added point to all points
				Point p(pointID++, temp, tempCF, sizeDim);
				points.push_back(p);

				if (FES >= maxFES) {
					on = 0;
					leader = sizePop;
					i = sizePop;
					break;
				}
			}
		}

		//COPY TEMPpopulation to actual POPULATION
		for (int i = 0; i < sizePop; i++) {
			popCF[i] = popCFTemp[i];
			for (int j = 0; j < sizeDim; j++) {
				pop[i][j] = popTemp[i][j];
			}
		}


		//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
		//WWWWWWWWW		II. cluster (k-means)
		//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
		
		KMeans kmeans(maxLeaders, 100);
		kmeans.run(points);

		//prepare possible leaders for III. All-To-Leader
		std::vector<Point> leaders;

		std::sort(points.begin(), points.end(), sortByClusterID);

		int id = points.at(0).getCluster();
		int sum = 0;
		leaders.push_back(points.at(0));
		int leaderIndex = 0;
		for (std::vector<Point>::iterator it = points.begin(); it != points.end(); it++) {
			if (it->getCluster() == id) {
				sum++;

				if (leaders.at(leaderIndex).getCF() > it->getCF()) {
					leaders.pop_back();
					leaders.push_back(*it);
				}
			}
			else {
				sum = 0;
				id = it->getCluster();
				leaders.push_back(*it);
				leaderIndex++;
			}
		}
		

		//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
		//WWWWWWWWW		III. ALL-TO-LEADER
		//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
		for (int i = 0; i < sizePop; i++) {
			//copy individual 
			popCFTemp[i] = popCF[i];
			tempCF = popCF[i];
			for (int j = 0; j < sizeDim; j++) {
				popTemp[i][j] = pop[i][j];
				temp[j] = pop[i][j];
			}
			//select leader
			Point leader = leaderRank(&leaders);

            if (FES >= maxFES) {
					on = 0;
					i = sizePop;
					break;
			}

			//move over path
			for (double t = step; t <= path; t += step) {

				//generate PRT vector for current leader
				for (int j = 0; j < sizeDim; j++) {
                    prt = 0.08 + 0.9 * ((double)FES/(double)maxFES);
					if (doubleRand(0, 1) > prt) {
						prtVector[j] = 0;
					}
					else {
						prtVector[j] = 1;
					}
				}
				
				for (int j = 0; j < sizeDim; j++) {
					temp[j] = pop[i][j] + (leader.getVal(j) - pop[i][j])*t*prtVector[j];

					//simple random strategy
					if (temp[j] > limMax) {
						temp[j] = doubleRand(limMin, limMax);
					}
					if (temp[j] < limMin) {
						temp[j] = doubleRand(limMin, limMax);
					}

				}
				//evaluate new substep (actualy fitness popTemp)
				fit(temp, &tempCF);
				FES++;
				if (tempCF < gBestCF) {
					gBestCF = tempCF;
					for (int j = 0; j < sizeDim; j++) {
						gBest[j] = temp[j];
					}
					//log if gBest is changed
					fprintf(logFile, "{%d,%.8e},", FES, gBestCF);
				}
				if (tempCF < popCFTemp[i]) {
					popCFTemp[i] = tempCF;
					for (int j = 0; j < sizeDim; j++) {
						popTemp[i][j] = temp[j];
					}
				}

				if (FES >= maxFES) {
					on = 0;
					i = sizePop;
					break;
				}
			}
		}

		//copy become new pop in next generation
		for (int i = 0; i < sizePop; i++) {
			for (int j = 0; j < sizeDim; j++) {
				pop[i][j] = popTemp[i][j];
			}
			popCF[i] = popCFTemp[i];
		}

		//Testing
		//#################################################################
		if (FES >= maxFES) {
			break;
		}
		//#################################################################
	}




	//save gBest to shared memory
	gBestCFs[runID] = gBestCF;
	for (int i = 0; i < sizeDim; i++) {
		gBests[runID][i] = gBest[i];
	}

	//close log file & log last gBest
	fprintf(logFile, "{%d,%.8e}", FES, gBestCF);
	fprintf(logFile, "}");
	fclose(logFile);

	//print gBest on console
	printf(" - id: %d | FES: %d | gBestCF: %.8e\n", runID, FES, gBestCF);

	//clear all points
	points.clear();
	

	//free memory
	delete[] gBest;
	for (int i = 0; i < sizePop; i++) {
		delete[] pop[i];
		delete[] popTemp[i];
	}
	delete[] popTemp;
	delete[] popCFTemp;
	delete[] pop;
	delete[] popCF;
	delete[] prtVector;
	delete[] temp;
}