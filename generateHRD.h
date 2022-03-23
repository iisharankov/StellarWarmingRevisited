#ifndef generateHRD
#define generateHRD

double flux(double* starVars);
double* fixDensity(double h, double tempCore);
void createMainSequence(double numStars, double minCoreTemp, double maxCoreTemp); // initial implementation of Eulers timestep


#endif