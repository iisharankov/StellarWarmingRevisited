#ifndef const_h
#define const_h

#define MASS_E 9.10938291e-31 // electron mass
#define MASS_P 1.67262178e-27 // proton mass
#define STEF_BOLT 5.670373e-8 // Stefan-Boltzmann constant
#define K_BOLT 1.381e-23 // Boltzmann constant
#define C 299792458 // speed of light
#define G 6.673e-11 // Gravitation constant
#define HBAR 1.054571817e-34

// Astronomy constants
#define RAD_CONST (4*STEF_BOLT)/C // radiation constant
#define FRAC_X 0.73 // Hydrogen mass fraction
#define FRAC_Y 0.25 // Helium mass fraction
#define FRAC_Z 0.02 // Metals mass fraction)
#define GAMMA 5/3 // adiabatic constant
#define M_SUN 1.989e30 // Kg
#define R_SUN 6.955e8 // meters
#define L_SUN 3.827e26 // watts

// Others
#define K 1.38e-23
#define pi  3.14159265
#define pi2 9.86960440

#define mu 1/(2.0 * FRAC_X + 0.75 * FRAC_Y + 0.5 * FRAC_Z)
#define r0 0.001  // m
#define S 1.0  // error tolerance


#endif