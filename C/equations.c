#include "const.h"
#include <math.h>
#include <stdlib.h>


// Pressures
double P(const double density, const double temp) {

    // Pressure degenerate - Eqn 5 in the project_description
    double P_deg = P_deg_const * pow(density / MASS_P, 5.0/3.0);

    // Pressure Ideal Gas
    double P_ig = (temp * density) * IdealGasFactor;

    // Pressure Radiative
    double P_rad = (RAD_CONST / 3.0) * temp*temp*temp*temp;

    return P_deg + P_ig + P_rad;
}

// Pressure differentials
double dPdp(const double density, const double temp) {

    // Degenerate
    double dPdp_deg = dPdp_deg_const * pow(density / MASS_P, 2.0/3.0);
    double dPdp_ig = IdealGasFactor * temp;  // Ideal Gas

    return dPdp_deg + dPdp_ig;
}


double dPdT(const double density, const double temp) {
    double dPdT_ig = density * IdealGasFactor;
    double dPdT_rad = (4.0 / 3.0) * RAD_CONST * temp*temp*temp;
    return dPdT_ig + dPdT_rad;
}

// Energy generation
double epsilon(const double density, const double temp) {
    double rTemp = temp / 1.0e6; // reducedTemp

    double epp = epsilon_epp * density * rTemp*rTemp*rTemp*rTemp;
    double ecno = epsilon_ecno * density * pow(rTemp, 19.9);

    return epp + ecno;
}


// Opacity
double Kappa(const double density, const double temp) {
    // double Kes = 0.02 * (1.0 +  FRAC_X);
    double Kff = Kappa_Kff_const * pow(density / 1.0e3, 0.7) * pow(temp, -3.5);
    double Khmin = Kappa_Khmin_const * pow(density / 1.0e3, 0.5) * pow(temp, 9.0);

    // return 1/ ( (1.0 / Khmin) + (1.0 / MAX(Kappa_Kes_const, Kff)));
    double max_result = MAX(Kappa_Kes_const, Kff);
    return (Khmin * max_result)/(Khmin + max_result);
}

// Stellar Structure ODEs

// dTdr = change in temperature over radius
double dTdr(const double radius, const double mass, const double density, const double temp, const double lum) {
    // second equation in the set of 5 equations in the project description file
    double radiative = (3.0 * Kappa(density, temp) * density * lum) / (
        16.0 * pi * RAD_CONST * C * (temp * temp * temp) * (radius * radius));

    double convection = (1.0 - (1.0 / GAMMA)) * (temp / P(density, temp)) * (
                (G * mass * density) / (radius * radius));

    return - MIN(radiative, convection);
}


// dpdr = change in density over radius
double dpdr(const double radius, const double mass, const double density, const double temp, const double lum) {
    // First last equation in the set of 5 equations in the project description file
    return -((G * mass * density / (radius*radius)) +
             dPdT(density, temp) * dTdr(radius, mass, density, temp, lum)) / (dPdp(density, temp));
}

// dMdr = change in mass over radius
double dMdr(const double radius, const double density) {
    // middle last equation in the set of 5 equations in the project description file
    return 4.0 * pi * (radius * radius) * density;
}

// DLdr = change in luminosity over radius
double dLdr(const double radius, const double density, const double temp) {
    // Second last equation in the set of 5 equations in the project description file
    return dMdr(radius, density) * epsilon(density, temp);
}

double dtaudr(const double density, const double temp) {
    // Last equation in the set of 5 equations in the project description file
    return Kappa(density, temp) * density;
}

// dPdr = change in pressure over radius
double dPdr(const double radius, const double mass, const double density) {
    return -(G * mass * density / (radius * radius));
}

// delta(tau) for optical depth limit
double dtau(const double radius, const double mass, const double density, const double temp, const double lum) {
    return (Kappa(density, temp) * (density * density)) / (double)(-(dpdr(radius, mass, density, temp, lum)));
}