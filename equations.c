#include "const.h"
#include <math.h>
#include <stdlib.h>


// Pressures
double P(const double density, const double temp) {

    // Degenerate Pressure
    double dDivMp = density / MASS_P;
    double P_deg = P_deg_const * cbrt(dDivMp*dDivMp*dDivMp*dDivMp*dDivMp);

    // Ideal Gas Pressure
    double P_ig = (temp * density) * IdealGasFactor;

    // Radiative Pressure
    double P_rad = (RAD_CONST / 3.0) * temp*temp*temp*temp;

    return P_deg + P_ig + P_rad;
}

// Pressure differentials
static double dPdp(const double density, const double temp) {

    // Degenerate   
    double dDivMp = density / MASS_P;
    double dPdp_deg = dPdp_deg_const * cbrt(dDivMp*dDivMp);
    double dPdp_ig = IdealGasFactor * temp;  // Ideal Gas

    return dPdp_deg + dPdp_ig;
}


double dPdT(const double density, const double temp) {
    double dPdT_ig = density * IdealGasFactor;
    double dPdT_rad = (4.0/3.0) * RAD_CONST * temp*temp*temp;
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
double Kappa(double density, double temp) {
    double tempToPowThreeAndHalf =  temp*temp*temp*sqrt(temp);
    double reducedDen = density / 1.0e3;
    double tempToPowNine = temp*temp*temp*temp*temp*temp*temp*temp*temp;

    double Kff = Kappa_Kff_const * pow(reducedDen, 0.7) / tempToPowThreeAndHalf; 
    double Khmin = Kappa_Khmin_const * sqrt(reducedDen) * tempToPowNine;

    double max_result = MAX(Kappa_Kes_const, Kff);
    return (Khmin * max_result)/(Khmin + max_result);
}


// Stellar Structure ODEs //

// dTdr = change in temperature over radius
double dTdr(const double radius, const double mass, const double density, const double temp, const double lum, double curKappa) {
    // second equation in the set of 5 equations in the project description file
    double radiative = (3.0 * curKappa * density * lum) / (
        16.0 * pi * RAD_CONST * C * (temp * temp * temp) * (radius * radius));

    double convection = (1.0 - (1.0 / GAMMA)) * (temp / P(density, temp)) * (
                (G * mass * density) / (radius * radius));

    return - MIN(radiative, convection);
}


// dpdr = change in density over radius
double dpdr(const double radius, const double mass, const double density, const double temp, const double lum, double curKappa) {
    // First last equation in the set of 5 equations in the project description file
    return -((G * mass * density / (radius*radius)) +
             dPdT(density, temp) * dTdr(radius, mass, density, temp, lum, curKappa)) / (dPdp(density, temp));
}

// dMdr = change in mass over radius
// double dMdr(const double radius, const double density) {
//     // middle last equation in the set of 5 equations in the project description file
//     return 4.0 * pi * radius * radius * density;
// }

// DLdr = change in luminosity over radius
double dLdr(const double radius, const double density, const double temp) {
    // Second last equation in the set of 5 equations in the project description file
    const double dMdr = 4.0 * pi * radius * radius * density;
    return dMdr * epsilon(density, temp);
}

// double dtaudr(const double density, const double temp) {
//     // Last equation in the set of 5 equations in the project description file
//     return Kappa(density, temp) * density;
// }

// dPdr = change in pressure over radius
// double dPdr(const double radius, const double mass, const double density) {
//     return -(G * mass * density / (radius * radius));
// }

// delta(tau) for optical depth limit
double dtau(const double radius, const double mass, const double density, const double temp, const double lum, double curKappa) {
    double dpdrSoln = dpdr(radius, mass, density, temp, lum, curKappa);
    if (dpdrSoln < 0) dpdrSoln *= -1.0; // Take absolute

    return (Kappa(density, temp) * (density * density)) / dpdrSoln;
}