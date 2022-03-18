#ifndef equations_h
#define equations_h


double P(const double density, const double temp); // initial implementation of Eulers timestep
// static double dPdp(const double density, const double temp); // Pressure differentials
double dPdT(const double density, const double temp);
double epsilon(const double density, const double temp); //Energy generation
double Kappa(const double density, const double temp); // Opacity

// Stellar Structure ODEs
double dpdr(const double radius, const double mass, const double density, const double temp, const double lum, double curKappa);
double dTdr(const double radius, const double mass, const double density, const double temp, const double lum, double curKappa);
double dMdr(const double radius, const double density);

double dLdr(const double radius, const double density, const double temp);
double dtaudr(const double density, const double temp);
// double dPdr(const double radius, const double mass, const double density);
double dtau(const double radius, const double mass, const double density, const double temp, const double lum, double curKappa);  // delta(tau) for optical depth limit

#endif