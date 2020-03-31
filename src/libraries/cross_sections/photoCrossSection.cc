#include <iostream>
#include <cmath>
#include <cstdlib>

#include "photoCrossSection.hh"
#include "helpers.hh"
#include "constants.hh"

photoCrossSection::photoCrossSection()
{
}

photoCrossSection::~photoCrossSection(){}

double photoCrossSection::sigma_pi0_p(double s, double t)
{
  double a = 14.79;
  double b = 7.37;
  double c = 7.49;
  double d = 5.36;

  double m_m = mpi0;
  double m_B = mN;

  double A = s;
  double B = sq(m_B) - sq(m_m) - s;
  double C = sq(m_m);
  double D = sq(B) - 4*A*C;

  double k = (s - sq(mN))/(2.*mN);
  double tmax = sq(m_m) - k*mN/s*(- B - sqrt(D));
  double tmin = sq(m_m) - k*mN/s*(- B + sqrt(D));
  double x = (tmax - t)/(tmax - tmin);
  
  return pow(a/s,b)*pow(x,c)*exp(d*sq(log(x)));
}

double photoCrossSection::sigma_pi0_n(double s, double t)
{
  return sigma_pi0_p(s,t);
}

double photoCrossSection::sigma_pip_n(double s, double t)
{
  double a = 4.03;
  double b = 8.52;
  double cplusd = 10.58;
  double cminusd = 0.67;
  double c = 0.5*(cplusd + cminusd);
  double d = 0.5*(cplusd - cminusd);

  double m_m = mpip;
  double m_B = mN;

  double A = s;
  double B = sq(m_B) - sq(m_m) - s;
  double C = sq(m_m);
  double D = sq(B) - 4*A*C;

  double k = (s - sq(mN))/(2.*mN);
  double tmax = sq(m_m) - k*mN/s*(- B - sqrt(D));
  double tmin = sq(m_m) - k*mN/s*(- B + sqrt(D));
  double x = (tmax - t)/(tmax - tmin);
  
  return pow(a/s,b)*pow(x,-c)*pow(1.-x,-d);
}

double photoCrossSection::sigma_pim_p(double s, double t)
{
  return sigma_pip_n(s,t);
}

