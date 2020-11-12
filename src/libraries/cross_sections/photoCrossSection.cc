#include <iostream>
#include <cmath>
#include "photoCrossSection.hh"
#include "helpers.hh"
#include "constants.hh"

using namespace std;

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

double photoCrossSection::sigma_pip_n(double s, double cosThetaCM)
{

  double A = 9.490;
  double b = 5.329;
  double c = 4.638;
      
  return pow(A/s,7)*pow(1-cosThetaCM,-b)*pow(1.+cosThetaCM,-c);
}

double photoCrossSection::sigma_pim_p(double s, double cosThetaCM)
{

  double A = 9.490;
  double b = 5.329;
  double c = 4.638;
      
  return pow(A/s,7)*pow(1-cosThetaCM,-b)*pow(1.+cosThetaCM,-c);
}


double photoCrossSection::sigma_rho0_p(double s, double cosThetaCM)
{
	const double b=-3.7;
	const double c=-2.2;
	const double a=5.82005e7;

return 0.75*(pow(s,-7)*a*pow(1-cosThetaCM,b)*pow(1+cosThetaCM,c));
}

double photoCrossSection::sigma_omega_p(double s, double cosThetaCM)
{
	const double b=-3.7;
	const double c=-2.2;
	const double a=5.82005e7;

return 0.25*(pow(s,-7)*a*pow(1-cosThetaCM,b)*pow(1+cosThetaCM,c));
} 

double photoCrossSection::sigma_psi_p(double s, double t)
{
  double sig0 = 11.3; //nb
  double beta = 1.3;
  double Lambda = 1.14;
  
  double x = (sq(mpsi) + 2*mN*mpsi)/(s - sq(mN));
  if ((x > 1.) || (x < 0.))
    return 0.;
  
  double sig = sig0 * pow(1. - x,beta);

  double pi = (sq(mN) - s)/(2.*sqrt(s));
  double pf = 0.5*sqrt(sq(sq(mN) - sq(mpsi))/s - 2.*(sq(mN) + sq(mpsi)) + s);
  double tmin = sq(sq(mpsi))/(4.*s) - sq(pi + pf);
  double tmax = sq(sq(mpsi))/(4.*s) - sq(pi - pf);
  
  return sig * dipole_F(t, tmin, tmax);
}

double photoCrossSection::sigma_psi_p(double s, double t, double QSq)
{
  double sig0 = 11.3; //nb
  double beta = 1.3;
  double Lambda = 1.14;
  double n = 2.44;
  
  double x = (sq(mpsi) + 2*mN*mpsi)/(s - sq(mN));
  if ((x > 1.) || (x < 0.))
    return 0.;
  
  double sig = sig0 * pow(1. - x,beta);

  double pi = (sq(mN) - s)/(2.*sqrt(s));
  double pf = 0.5*sqrt(sq(sq(mN) - sq(mpsi))/s - 2.*(sq(mN) + sq(mpsi)) + s);
  double tmin = sq(sq(mpsi) + QSq)/(4.*s) - sq(pi + pf);
  double tmax = sq(sq(mpsi) + QSq)/(4.*s) - sq(pi - pf);
  
  return sig * dipole_F(t, tmin, tmax) * pow(sq(mpsi)/(QSq + sq(mpsi)),n);
}

double photoCrossSection::R_psi_p(double QSq)
{
  double a = 2.164;
  double n = 2.131;
  return pow((a*sq(mpsi) + QSq)/(a*sq(mpsi)),n) - 1.;
}

double photoCrossSection::dipole_F(double t, double tmin, double tmax)
{
  double Lambda = 1.14;

  return pow(sq(Lambda) - t,-4.) * 3./(pow(sq(Lambda) - tmin,-3.) - pow(sq(Lambda) - tmax,-3.));
}
