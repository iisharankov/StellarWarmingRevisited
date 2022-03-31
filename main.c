#include "const.h"
#include "equations.h"
#include "starGenerator.h"
#include "generateHRD.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main() {
    // createStar(1000.0, 0.5e3, 1.5e7, 0);
    // double* test = createStar(1000.0, 31531.25, 3981071.7055349695, 0);
    // free(test);
    // double* t = createStar(1000.0, 47146.875, 3981071.7055349695, 1);
    // createStar(1000.0, 54954.6875, 3981071.7055349695, 0);
    
    // int numThreads = 1;
    // int stars = 500;
    // struct timespec start, finish;
    // double elapsed;
    // clock_gettime(CLOCK_MONOTONIC, &start);

    // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    // // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    // // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    // // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    // // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    // // createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));


    // clock_gettime(CLOCK_MONOTONIC, &finish);
    // elapsed = (finish.tv_sec - start.tv_sec);
    // elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    // printf("Total time was %.3f - Average time for %d star was %.3f seconds\n", elapsed,  stars, elapsed / 10.0);
    createMainSequence(1, 100, pow(10.0, 6.6), pow(10.0, 7.4));

    // createStar(1000.0, 1000.0, pow(10, 6.6), 1);
    // bisectStar(1000.0, pow(10, 6.6));

    return 0;
}