#include "const.h"
#include "equations.h"
#include "starGenerator.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>


void writeIterationToOpenFile(FILE* fPtr, int index, double new_rk4[], char *dataLine) {
    
    // These next lines append to the dLdr list, kappa list, pressure list, and dlogPdlogT lists.
    double density  = new_rk4[(index*6) + 0];
    double temp     = new_rk4[(index*6) + 1];
    double mass     = new_rk4[(index*6) + 2];
    double lum      = new_rk4[(index*6) + 3];
    double tau      = new_rk4[(index*6) + 4];
    double radius   = new_rk4[(index*6) + 5];
    

    double new_dLdr =  dLdr(radius, density, temp);
    double kappa = Kappa(density, temp);
    double pressure = P(density, temp);

    const double dPdr = -(G * mass * density / (radius * radius));
    double new_dlogPdlogT = (temp / pressure) * (dPdr / dTdr(radius, mass, density, temp, lum, kappa));
    double new_dTau = dtau(radius, mass, density, temp, lum, kappa);

    // Create the line of vectors to save into the file
    sprintf(dataLine, "%e\t%e\t%e\t%e\t%f\t%e\t%e\t%e\t%e\t%e\t%e\t\n", 
    density, temp, mass, lum, tau, radius, new_dLdr, kappa, pressure, new_dlogPdlogT, new_dTau);

    fputs(dataLine, fPtr); // Write data to file
}


void func(double* restrict inputVec, double radius, double h, double* restrict outputVec) {

    // Extract array to individual variable names. Divide by divisor if needed by rk4
    double density = inputVec[0];
    double temp  = inputVec[1];
    double mass = inputVec[2];
    double lum = inputVec[3];

    double curKappa = Kappa(density, temp);
    double cur_dMdr = 4.0 * pi * radius * radius * density;

    // Output vec is of order "Density, Temp, Mass, Lum, Tau"
    outputVec[0] = h * dpdr(radius, mass, density, temp, lum, curKappa);
    outputVec[1] = h * dTdr(radius, mass, density, temp, lum, curKappa);
    outputVec[2] = h * cur_dMdr; //dMdr(radius, density);
    outputVec[3] = h * cur_dMdr * epsilon(density, temp); // dLdr(radius, density, temp);
    outputVec[4] = h * curKappa * density;
}


void rk4(double* restrict y, double h) {
    // This is the runge kutta method. It runs the function func on the current 
    // dependant variable values and the radius
    // y: current values for dependent variables
    // r: radius, the independent variable
    // h: step-size

    double array[KROWS][KCOLS];

    // calculate k1
    func(y, y[5], h, array[1]);

    // Update tempVector (array[0]) with modified k1 (array[1]) from func
    for(int i = 0; i < 5; i++) array[0][i] = y[i] + (array[1][i]/2.0); 
    
    // calculate k2
    func(array[0], y[5] + h/2.0, h, array[2]);

    // Update tempVector (array[0]) with modified k2 (array[2]) from func
    for(int i = 0; i < 5; i++) array[0][i] = y[i] + (array[2][i]/2.0);
    
    // calculate k3
    func(array[0], y[5] + h/2.0, h, array[3]);

    // Update tempVector (array[0]) with modified k3 (array[3]) from func
    for(int i = 0; i < 5; i++) array[0][i] = y[i] + array[3][i];

    // calculate k4
    func(array[0], y[5] + h, h, array[4]);

    // Update vars using the runge kutta formula
    for(int i = 0; i < 5; i++) { 
        // This is the runge kutta formula 
        double rkSum = (array[1][i] + (2.0*array[2][i]) + (2.0*array[3][i]) + array[4][i]);
        y[i] += (rkSum/6.0); 
    }
    y[5] += h;
}

int opticalDepthLimit(double* restrict newVars) {
    // tests the optical depth to make sure that d_tau is not very small,
    // or that the mass hasn't blown up (1e3 Mass sun). While those conditions 
    // have not been met, returns False
    double newKappa = Kappa(newVars[0], newVars[1]);
    double new_dTau = dtau(newVars[5], newVars[2], newVars[0], newVars[1], newVars[3], newKappa);

    if (new_dTau < 0.001 || newVars[2] > OneEThreeM_SUN) return 1;
    return 0;
}


int radiativeStar( double* restrict taus, int tauSize){

    double lastTau = taus[tauSize-1];
    if (lastTau < 0) lastTau *= -1.0;

    double reducedTaus[tauSize];
    for (int i=0; i<tauSize; i++) reducedTaus[i] = lastTau - taus[i];

    int index = 0;
    for(int i = 1; i < tauSize; i++) {
        // printf("%d - %f %f %d\n", i, reducedTaus[i], reducedTaus[index], reducedTaus[i] < reducedTaus[index]);
        if (reducedTaus[i] < reducedTaus[index]) index = i;
    }

    // printf("%f, %d, %d, %d\n", lastTau, index, tauSize, tauSize-index);
    // int gap = tauSize - index + 20;
    // for (int i = 0; i < gap; i++) printf("%f ", taus[tauSize-gap+i]);
    // printf("\n");
    // printf("\n");

    // if (tauSize - index == 499) exit(EXIT_FAILURE);
    return index;
}


double* createStar(double radius, double density, double temp, int writeData) {

    // Create arrays needed for storing data
    double starData[TEMPDATASIZE * (6 * sizeof(double))];
    double taus[TEMPDATASIZE];

    // File pointer to hold reference to our file and allocate storage line for writing to
    FILE * fPtr; // Below: size adds (6 * sizeof(double)) for the whitespace
    char *dataLine = (char*)malloc( (12 * sizeof(double)) + (6 * sizeof(double)) ); 

    if (writeData) {
        fPtr = fopen("file1.txt", "a");  // Open file in append mode
        if(fPtr == NULL) { // If NULL, filefopen unsuccessful
            printf("Unable to create file.\n");
            exit(EXIT_FAILURE); // File not created, exit
        }
    }  

    // allocate the var vector we'll update
    double* vars = malloc(6*sizeof(double));
    vars[0] = density;
    vars[1] = temp;
    vars[2] = (4.0 / 3.0) * pi * (radius*radius*radius) * density;
    vars[3] = vars[2] * epsilon(density, temp);
    vars[4] = Kappa(density, temp) * density;
    vars[5] = radius;

    // Write initial variables if writing to file
    if (writeData) writeIterationToOpenFile(fPtr, 0, vars, dataLine);

    int counter = 0;
    double stepSize = r0;

    do {
        for (int j = 0; j < LOOP_UNROLL_FACTOR; j++) {
            rk4(vars, stepSize);
            for(int i = 0; i < 6; i++) starData[counter*6 + i] = vars[i];
            counter++;

            // Finds the new dr (increment in radius) by checking if the temp is < 5e4
            stepSize = vars[1] < 5e4 ? (0.00001 * vars[5]) + 1000 : (0.001 * vars[5]) + 1000;
        }

        if (counter >= TEMPDATASIZE) {
                /*
                    If writeData, then write all of starData to file so we can reset counter and overwrite. If
                    !writeData, then we don't care about data so flush starData. Due to some mathematical 
                    properties of the DEs in question, we only need the last few thousand to properly calculate 
                    radiativeStar. 
                */
                if (writeData) for (int i=0; i<counter; i++) writeIterationToOpenFile(fPtr, i, starData, dataLine);
                else for (int j = 0; j < TEMPDATASIZE * 6; j++) starData[j] = 0; 
                counter = 0;
        }
    } while ( !opticalDepthLimit(vars) );
    free (vars);    

    // Write last data, delete dataLine and close file
    if (writeData) {
        for (int i=0; i<counter; i++) writeIterationToOpenFile(fPtr, i, starData, dataLine);
        fclose(fPtr);
    }
    free (dataLine);

    // printf("%d - %d\n", TEMPDATASIZE, counter*6);
    // for (int i=0; i<counter; i++) taus[i] = tempData[i][4];
    for (int i=0; i<counter; i++) taus[i] = starData[i*6 + 4];
    int tauMinIndex = radiativeStar(taus, counter);
    
    // printf("%d - %d\n", tauMinIndex, counter);
    // printf("%e\n", vars[3]);

    double* output = malloc(6*sizeof(double));
    for (int i=0; i<6; i++) output[i] = starData[tauMinIndex*6 + i];

    return output;
}   



