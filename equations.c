#include "const.h"
#include <math.h>
#include <stdlib.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// Pressures
double P(double density, double temp) {

    // Pressure degenerate - Eqn 5 in the project_description
    double P_deg = (pow(3.0 * pi2, 2.0/3.0) * (HBAR * HBAR) * pow(density / MASS_P, 5.0/3.0)) / (5.0 * MASS_E);

    // Pressure Ideal Gas
    double P_ig = (K * temp * density) / (mu * MASS_P);

    // Pressure Radiative
    double P_rad = (1.0 / 3.0) * RAD_CONST * pow(temp, 4.0);

    return P_deg + P_ig + P_rad;
}

// Pressure differentials
double dPdp(double density, double temp) {

    // Degenerate
    double dPdp_deg = (pow(3.0 * pi2, 2.0/3.0) * pow(HBAR, 2.0) * pow(density / MASS_P, 2.0/3.0)  ) / (3.0 * MASS_P * MASS_E);

    // Ideal Gas
    double dPdp_ig = (K * temp) / (mu * MASS_P);

    return dPdp_deg + dPdp_ig;
}


double dPdT(double density, double temp) {
    double dPdT_ig = (density * K) / (mu * MASS_P);
    double dPdT_rad = (4.0 / 3.0) * RAD_CONST * pow(temp, 3.0);
    return dPdT_ig + dPdT_rad;
}

// Energy generation
double epsilon(double density, double temp) {
    // Value for epsilon, used often below
    double epp = (1.07e-7) * (density / 1.0e5) * pow(FRAC_X, 2) * pow(temp / 1.0e6, 4.0);
    double ecno = (8.24e-26) * (density / 1.0e5) * 0.03 * pow(FRAC_X, 2.0) * pow(temp / 1.0e6, 19.9);

    return epp + ecno;
}


// Opacity
double Kappa(double density, double temp) {
    double Kes = 0.02 * (1.0 +  FRAC_X);
    double Kff = 1.0e24 * ( FRAC_Z + 0.0001) * pow(density / 1.0e3, 0.7) * pow(temp, -3.5);
    double Khminus = 2.5e-32 * ( FRAC_Z / 0.02) * pow(density / 1.0e3, 0.5) * pow(temp, 9.0);

    return 1/ ( (1.0 / Khminus) + (1.0 / MAX(Kes, Kff)) );
}

double dTdr(double radius, double mass, double density, double temp, double lum) {
    // second equation in the set of 5 equations in the project description file
    double dTdr_rad = (3.0 * Kappa(density, temp) * density * lum) / (
                16.0 * pi * RAD_CONST * C * (temp * temp * temp) * (radius * radius));

    double dTdr_conv = (1.0 - (1.0 / GAMMA)) * (temp / P(density, temp)) * (
                (G * mass * density) / (radius * radius));

    return - MIN(dTdr_rad, dTdr_conv);
}

// Stellar Structure ODEs
double dpdr(double radius, double mass, double density, double temp, double lum) {
    // First last equation in the set of 5 equations in the project description file
    return -((G * mass * density / pow(radius, 2.0)) +
             dPdT(density, temp) * dTdr(radius, mass, density, temp, lum)) / (dPdp(density, temp));
}


double dMdr(double radius, double density) {
    // middle last equation in the set of 5 equations in the project description file
    return 4.0 * pi * (radius * radius) * density;
}

double dLdr(double radius, double density, double temp) {
    // Second last equation in the set of 5 equations in the project description file
    return dMdr(radius, density) * epsilon(density, temp);
}

double dtaudr(double density, double temp) {
    // Last equation in the set of 5 equations in the project description file
    return Kappa(density, temp) * density;
}

double dPdr(double radius, double mass, double density) {
    return -(G * mass * density / (radius * radius));
}

// delta(tau) for optical depth limit
double dtau(double radius, double mass, double density, double temp, double lum) {
    return (Kappa(density, temp) * (density * density)) / (double)(-(dpdr(radius, mass, density, temp, lum)));
}