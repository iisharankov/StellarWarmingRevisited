#include <stdio.h>
#ifndef starGenerator
#define starGenerator


#define VECTORSIZE 5


double* func(double* density, double temp, double h); // initial implementation of Eulers timestep
double* rk4(double* y, double h); // initial implementation of Eulers timestep

int opticalDepthLimit(double* newVars); // 

double* mainLoop(double radius, double density, double temp, int writeData); // 
void writeIterationToOpenFile(FILE* fPtr, double* new_rk4);

#endif