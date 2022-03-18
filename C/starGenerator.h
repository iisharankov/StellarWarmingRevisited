#include <stdio.h>
#ifndef starGenerator
#define starGenerator


#define KROWS 5
#define KCOLS 5
#define TEMPDATASIZE 5000


void writeIterationToOpenFile(FILE* fPtr, int index, double new_rk4[], char* dataLinee);

// RK4 functions
void func(double* dep_var, double radius, double h, double* outputVec);
void rk4(double* y, double h); 

int opticalDepthLimit(double* newVars); 
int radiativeStar(double* taus, int tauSize);


double* createStar(double radius, double density, double temp, int writeData); // 

#endif