#ifndef generateHRD
#define generateHRD

void* parallelThreadWorker(void* arg);
void storeData(int offset, double* arrayAddress, double* curStar);

double flux(double* starVars);
double* bisectStar(double h, double tempCore);

void createMainSequence(int numThreads, int numStars, double minCoreTemp, double maxCoreTemp);

// arguments for threads
typedef struct tData {
    double* coreTemps;
    double*  HRDataAddress;
    int start;
    int end;
} tData;


#endif