#include "const.h"
#include "equations.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define MAX(x, y) (((x) > (y)) ? (x) : (y))


void writeIterationToOpenFile(FILE* fPtr, double* new_rk4) {
    
    //////////////
    // These next lines append to the dLdr list, kappa list, pressure list, and dlogPdlogT lists.
    double newDist = new_rk4[0];
    double newTemp  = new_rk4[1];
    double newMass = new_rk4[2];
    double newLum = new_rk4[3];
    double newTau = new_rk4[4];
    double newRadius = new_rk4[5];
    

    double new_dLdr =  dLdr(newRadius, newDist, newTemp);
    double newKappa = Kappa(newDist, newTemp);
    double newDensity = P(newDist, newTemp);
    double new_dlogPdlogT = (newTemp / newDensity) * (dPdr(newRadius, newMass, newDist) / dTdr(newRadius, newMass, newDist, newTemp, newLum));

    double new_dTau = dtau(newRadius, newMass, newDist, newTemp, newLum);
    /////////////

    
    // Allocates storage
    char *dataLine = (char*)malloc( (12 * sizeof(double)) + (6 * sizeof(double))); // (6 * sizeof(double)) for the whitespace
    
    // Create the line of vectores to save into the file
    sprintf(dataLine, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t", 
    newDist, newTemp, newMass, newLum, newTau, newRadius, new_dLdr, newKappa, newDensity, new_dlogPdlogT, new_dTau);


    fputs(dataLine, fPtr); // Write data to file
    fputs("\n", fPtr); // Write new line
}


double* func(double* dep_var, double radius, double h) {

    // Extract array to individual variable names. Divide by divisor if needed by rk4
    double density = dep_var[0];
    double temp  = dep_var[1];
    double mass = dep_var[2];
    double lum = dep_var[3];

    double newDensity = dpdr(radius, mass, density, temp, lum);
    double newTemp = dTdr(radius, mass, density, temp, lum);
    double newMass = dMdr(radius, density);
    double newLum = dLdr(radius, density, temp);
    double newTau = dtaudr(density, temp);

   
    double* outputVec = malloc(5*sizeof(double));
    outputVec[0] = h * newDensity;
    outputVec[1] = h * newTemp;
    outputVec[2] = h * newMass;
    outputVec[3] = h * newLum;
    outputVec[4] = h * newTau;

    return outputVec;
}


double* rk4(double* y, double h) {
    // This is the runge kutta method. It runs the function func on the current 
    // dependant variable values and the radius
    // y: current values for dependent variables
    // r: radius, the independent variable
    // h: step-size

    double* inputVector = malloc(5*sizeof(double));

    // calculate k1
    double r = y[5];
    double* k1 = func(y, r, h);

    // calculate k2
    for(int i = 0; i < 5; i++) inputVector[i] = y[i] + (k1[i]/2.0); 
    double* k2 = func(inputVector, r + h/2.0, h);

    // calculate k3
    for(int i = 0; i < 5; i++) inputVector[i] = y[i] + (k2[i]/2.0);
    double* k3 = func(inputVector, r + h/2.0, h);
    
    // calculate k4
    for(int i = 0; i < 5; i++) inputVector[i] = y[i] + k3[i];
    double* k4 = func(inputVector, r + h, h);


    double* result = malloc(6*sizeof(double));

    // This is the runge kutta formula 
    for(int i = 0; i < 5; i++) { 
        double rkSum = (k1[i] + (2.0*k2[i]) + (2.0*k3[i]) + k4[i]);
        result[i] = y[i] + (rkSum/6.0); 
    }
    result[5] = r + h;

    // Free everything not needed
    free(inputVector);
    free(k1);
    free(k2);
    free(k3);
    free(k4);

    return result;
}

int opticalDepthLimit(double* newVars) {
    // tests the optical depth to make sure that d_tau is not very small,
    // or that the mass hasn't blown up (1e3 Mass sun). While those conditions 
    // have not been met, returns False

    double new_dTau = dtau(newVars[5], newVars[2], newVars[0], newVars[1], newVars[3]);

    if (new_dTau < 0.001 || newVars[2] > 1e3 * M_SUN) return 1;
    return 0;
}


double* mainLoop(double radius, double density, double temp, int writeData) {

    FILE * fPtr;  // File pointer to hold reference to our file 
    if (writeData) {
        fPtr = fopen("file1.txt", "a");  // Open file in append mode
        if(fPtr == NULL) { // If NULL, filefopen unsuccessful
            printf("Unable to create file.\n");
            exit(EXIT_FAILURE); // File not created, exit
        }
    }  

    double* newVars = malloc(6*sizeof(double));
    double* oldVars = malloc(6*sizeof(double));
    newVars[0] = density;
    newVars[1] = temp;
    newVars[2] = (4.0 / 3.0) * pi * (radius*radius*radius) * density;
    newVars[3] = newVars[2] * epsilon(density, temp);
    newVars[4] = Kappa(density, temp) * density;
    newVars[5] = radius;

    if (writeData) writeIterationToOpenFile(fPtr, newVars);  // Write initial variables

    double stepSize = r0;
    do {
        // Overwrite
        for (int i = 0; i < 6; i++) { oldVars[i] = newVars[i]; }
    
        newVars = rk4(oldVars, stepSize); 
        if (writeData)  writeIterationToOpenFile(fPtr, newVars);

        // Finds the new dr (increment in radius) by checking if the temp is < 5e4
        if (newVars[1] < 5e4) { 
            stepSize = (0.00001 * newVars[5]) + 1000;
        } else {
            stepSize = (0.001 * newVars[5]) + 1000;
        }

    } while ( !opticalDepthLimit(newVars));

    if (writeData)  fclose(fPtr); // Close file to save file data

    free (oldVars);
    return newVars;
}   




