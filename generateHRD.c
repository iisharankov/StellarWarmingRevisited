

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#include "starGenerator.h"
#include "generateHRD.h"
#include "const.h"


//     double (*tempData)[TEMPDATASIZE][6]     = ((tData*) arg)->tempData;

void* parallelThreadWorker(void* arg) {
    double* coreTemps   = ((tData*) arg)->coreTemps;
    double* pHRData     = ((tData*) arg)->HRDataAddress;
    long int start      = ((tData*) arg)->start;
    long int end        = ((tData*) arg)->end;

    for(int starNo = start; starNo < end; starNo++) {

        // call bisectStar on star, store result into HRData, and free curStar
        double* curStar = bisectStar(1000.0, coreTemps[starNo]);
        storeData(starNo, pHRData, curStar);\
        free(curStar);
    }
    return NULL;
}


void storeData(int offset, double HRData[], double* restrict curStar) {
    double surfaceTemp = sqrt(sqrt(curStar[3] /  (stefanBoltRadii * curStar[5]*curStar[5])  ));
    
    // store into array at correct location
    int loc = offset * 4; //* sizeof(double);
    HRData[loc + 0] = surfaceTemp;
    HRData[loc + 1] = curStar[3];
    HRData[loc + 2] = curStar[5];
    HRData[loc + 3] = curStar[2];
}


double flux(double* restrict starVars) {  // TODO: no self.Rstar index, just using last
    double Temp = starVars[1];
    double Lum = starVars[3];
    double radius = starVars[5];

    // Note: Lum = 4pi BOLT * r^2 * T^4
    double LumCorrected = (stefanBoltRadii * (radius*radius) * (Temp*Temp*Temp*Temp));
    return (Lum - LumCorrected) / sqrt(Lum * LumCorrected);
}


double* bisectStar(double h, double coreTemp) {


    double tol = 0.01;
    double starADensity = 0.3e3;
    double starBDensity = 500.0e3;
    double starCDensity = 500.3e3/2.0;

    double* starA = createStar(h, starADensity, coreTemp, 0);
    double* starB = createStar(h, starBDensity, coreTemp, 0);
    double* starC = createStar(h, starCDensity, coreTemp, 0);
    
    // Bisection start
    while ((starBDensity - starADensity) / 2.0 > tol) {

        if (flux(starA) * flux(starC) < 0) {
            for (int i = 0; i < 6; i++) starB[i] = starC[i];
            starBDensity = starCDensity;
        } else {
            for (int i = 0; i < 6; i++) starA[i] = starC[i];
            starADensity = starCDensity;
        }
        free(starC);

        starCDensity = (starADensity + starBDensity) / 2.0;
        starC = createStar(h, starCDensity, coreTemp, 0);
    }
    free(starC);

    starCDensity = MAX(starADensity, starBDensity);
    starC = createStar(h, starCDensity, coreTemp, 0);
    // Bisection end

    free(starA);
    free(starB);
    return starC;
}


void createMainSequence(int numThreads, int numStars, double minCoreTemp, double maxCoreTemp) {
    double HRData[numStars * (6 * sizeof(double))];
    double* coreTemps = malloc(numStars*sizeof(double));
    // printf("%d %f %f\n", numStars, minCoreTemp, maxCoreTemp);

    // Create stepsize and make list or ranged temperatures
    double stepSize = (maxCoreTemp - minCoreTemp) / numStars;
    for (int i=0; i < numStars; i++) coreTemps[i] = minCoreTemp + (i * stepSize); 

    if (numThreads > 1) {
        int splits = (int)((double)numStars/(double)numThreads);

        // The data to pass to each thread. 
        tData threadData[numThreads];
        for (int i = 0; i < numThreads; i++) {
            threadData[i].coreTemps = coreTemps;
            threadData[i].HRDataAddress = HRData;
            threadData[i].start = i*splits;  // start
            threadData[i].end   = (i+1)*splits; // end

            // If splits is not an even cut, make sure last thread goes to real end. i.e. consider
            // 10 elem / 3 threads -> (0-3), (3-6), and (6-9). So make last (6-10) to cover domain
            if (numThreads - i == 1) threadData[i].end = numStars;
        }

        pthread_t threads[numThreads];
        for (int i = 0; i < numThreads; i++) {
            pthread_create(&(threads[i]), NULL, parallelThreadWorker, (void *) &threadData[i]);
        }

        // Join threads when complete. Nothing needs to be returned from threads
        for (int i = 0; i < numThreads; i++) pthread_join(threads[i], NULL);

    } else {  // run in serial
        for (int starNo=0; starNo<numStars; starNo++) {
            // printf("\nStar %d - temp=%f\n", starNo, coreTemps[starNo]);

            double* curStar = bisectStar(1000.0, coreTemps[starNo]);
            storeData(starNo, HRData, curStar);
            free(curStar);
        }
    }
    free(coreTemps);

    char *dataLine = (char*)malloc( 8 * sizeof(double) ); // Allocate storage
    FILE * fPtr;  // File pointer to hold reference to our file 
    fPtr = fopen("hrData.txt", "a");  // Open file in append mode
    if(fPtr == NULL) { // If NULL, file fopen unsuccessful
        printf("Unable to create file.\n");
        exit(EXIT_FAILURE); // File not created, exit
    }

    // Store the data into the file
    for (int i=0; i<numStars; i++) {
        int j = i*4; // offset

        // Create the line of vectors to save into the file
        sprintf(dataLine, "%e\t%e\t%e\t%e\n", HRData[j], HRData[j+1], HRData[j+2], HRData[j+3]);
        fputs(dataLine, fPtr); // Write data to file
    }

    // Free mallocs and close file
    free(dataLine);
    fclose(fPtr);

}

