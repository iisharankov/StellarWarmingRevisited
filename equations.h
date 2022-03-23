#ifndef equations_h
#define equations_h



double P(double density, double temp); // initial implementation of Eulers timestep
double dPdp(double density, double temp); // Pressure differentials
double dPdT(double density, double temp);
double epsilon(double density, double temp); //Energy generation
double Kappa(double density, double temp); // Opacity

// Stellar Structure ODEs
double dpdr(double radius, double mass, double density, double temp, double lum);
double dTdr(double radius, double mass, double density, double temp, double lum);
double dMdr(double radius, double density);

double dLdr(double radius, double density, double temp);
double dtaudr(double density, double temp);
double dPdr(double radius, double mass, double density);
double dtau(double radius, double mass, double density, double temp, double lum);  // delta(tau) for optical depth limit

#endif