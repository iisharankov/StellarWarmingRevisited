#include "const.h"
#include "equations.h"
#include "starGenerator.h"
#include "generateHRD.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    mainLoop(1000.0, 0.5e3, 1.5e7, 0);

    //createMainSequence(100, pow(10.0, 6.6), pow(10.0, 7.4));
    return 0;
}