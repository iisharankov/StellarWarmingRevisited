#include "const.h"
#include "equations.h"
#include "starGenerator.h"
#include "generateHRD.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


void testTimings() {
    int numThreads = 1;
    int stars = 500;
    struct timespec start, finish;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // This was in a block and not a for loop as I noticed inconsistencies 
    // with vectorization enabled, even with a conditional if added.
    createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));
    createMainSequence(numThreads, stars, pow(10.0, 6.6), pow(10.0, 7.4));

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Total time was %.3f - Average time for %d star was %.3f seconds\n", elapsed,  stars, elapsed / 5.0);
}

int main() {
    // double* test1 = createStar(1000.0, 31531., 3981071., 0);
    // double* test2 = createStar(1000.0, 54954., 3981071., 0);
    // free(test1);
    // free(test2);
    
    
    createMainSequence(2, 2, pow(10.0, 6.6), pow(10.0, 7.4));
    return 0;
}