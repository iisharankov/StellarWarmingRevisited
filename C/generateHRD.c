

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "starGenerator.h"
#include "const.h"


double flux(double* starVars) {
    // TODO: no self.Rstar index, just using last

    double Lum1 = starVars[3];  // Note: Lum = 4pi BOLT * r^2 * T^4
    double Lum2 = (stefanBoltRadii * (starVars[5]*starVars[5]) * (starVars[1]*starVars[1]*starVars[1]*starVars[1]));

    return (Lum1 - Lum2) / sqrt(Lum1 * Lum2);
}

double* fixDensity(double h, double tempCore) {
    double tol = 0.01;
    double starADensity = 0.3e3;
    double starBDensity = 500.0e3;
    double starCDensity = 500.3e3/2.0;

    double* starA = mainLoop(h, starADensity, tempCore, 0);
    double* starB = mainLoop(h, starBDensity, tempCore, 0);
    double* starC = mainLoop(h, starCDensity, tempCore, 0);
    
    // Bisection start
    double starCrho = (starADensity + starBDensity) / 2.0;

    while ((starBDensity - starADensity) / 2.0 > tol) {
        if (flux(starA) * flux(starC) < 0) {
            for (int i = 0; i < 6; i++) { starB[i] = starC[i]; }
            starBDensity = starCDensity;
        } else {
            for (int i = 0; i < 6; i++) { starA[i] = starC[i]; }
            starADensity = starCDensity;

        }
    
        starCDensity = (starADensity + starBDensity) / 2.0;
        starC = mainLoop(h, starCDensity, tempCore, 0);
    }

    starCrho = MAX(starADensity, starBDensity);
    starC = mainLoop(h, starCrho, tempCore, 0);
    // Bisection end

    free(starA);
    free(starB);
    return starC;
}


void createMainSequence(double numStars, double minCoreTemp, double maxCoreTemp) {
    
    FILE * fPtr;  // File pointer to hold reference to our file 
    fPtr = fopen("hrData.txt", "a");  // Open file in append mode
    if(fPtr == NULL) { // If NULL, filefopen unsuccessful
        printf("Unable to create file.\n");
        exit(EXIT_FAILURE); // File not created, exit
    }

    double* coreTemps = malloc(numStars*sizeof(double));

    // Create stepsize and make list or ranged temperatures
    double stepSize = (maxCoreTemp - minCoreTemp) / numStars;
    for (int i=0; i < numStars; i++) coreTemps[i] = minCoreTemp + (i * stepSize); 

    
    
    for (int i=0; i<numStars; i++) {
        double* curStar = fixDensity(1000.0, coreTemps[i]);
        
        double denom = (stefanBoltRadii * curStar[5]*curStar[5]);
        double surfaceTemp = sqrt(sqrt(curStar[3] / denom));

        
        // Create the line of vectors to save into the file
        char *dataLine = (char*)malloc( 8 * sizeof(double) ); // Allocate storage
        sprintf(dataLine, "%e\t%e\t%e\t%e", surfaceTemp, curStar[3], curStar[5], curStar[2]);
        fputs(dataLine, fPtr); // Write data to file
        fputs("\n", fPtr); // Write new line

        free(dataLine);
    }


    // close file
    fclose(fPtr); 
}

