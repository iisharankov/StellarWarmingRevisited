#include "const.h"
#include "equations.h"
#include "starGenerator.h"
#include "generateHRD.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    
    // double* firstVec = malloc(5*sizeof(double));
    // firstVec[0] = 1;
    // firstVec[1] = 2;
    // firstVec[2] = 3;
    // firstVec[3] = 4;
    // firstVec[4] = 5;

    // double radius = 4;
    // double h = 4;
    // double* result = rk4(firstVec, radius, h);
    // double newRadius = radius + h;

    // for (int i = 0; i < 5; i++) { printf("%e\n",result[i]); }
    // printf("%e\n", newRadius);


    // writeIterationToFile(result, newRadius);

    // mainLoop(1000.0, 0.5e3, 1.5e7, 0);

    createMainSequence(100, pow(10.0, 6.6), pow(10.0, 7.4));
    // double outputVec[5];
    // for (int i = 0; i<5; i++) { printf("%e\n",firstVec[i]); }
    // *outputVec = func(firstVec, radius);
    // for (int i = 0; i<5; i++) { printf("%e\n",firstVec[i]); }
    return 0;
}