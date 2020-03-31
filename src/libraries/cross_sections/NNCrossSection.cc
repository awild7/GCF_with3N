#include <iostream>
#include <cmath>
#include <cstdlib>

#include "NNCrossSection.hh"
#include "helpers.hh"
#include "constants.hh"

NNCrossSection::NNCrossSection()
{
  // Set defaults
  myParam=Panin;
}

NNCrossSection::NNCrossSection(csParam thisParam)
{
  std::cerr << "Cross_Sections: you have selected parameterization: " << thisParam <<"\n";
  myParam=thisParam;
}

NNCrossSection::~NNCrossSection(){}

double NNCrossSection::sigma_pp(double s, double t)
{
  switch (myParam)
    {
    case Panin:
      return sigma_pp_Panin(s,t);
    default:
      {
	std::cerr << "Invalid cross section parameterization! Double check and fix!\n";
        exit(-1);
      }
    }
  return 0;
}

double NNCrossSection::sigma_pp_Panin(double s, double t)
{

  double tmin = 4 * sq(mN) - s;
  double A = 0., B = 0.;

  if (t < tmin / 2.)
    t = tmin - t;

  if (s < 4.79 && s > 3.8)
  {
    A = -1799.02 + 937.563 * s - 19.1253 * pow(s, 2) - 52.883 * pow(s, 3) + 6.88961 * pow(s, 4);
    //B = -776.822 + 586.016 * s - 175.347 * pow(s, 2) + 26.1823 * pow(s, 3) - 1.94889 * pow(s, 4) + 0.0578352 * pow(s, 5);
    B = 3.69683 -0.221464 * s -0.185574 * pow(s, 2) -0.0469019 * pow(s, 3)  -0.00587074* pow(s, 4) +0.000975612 * pow(s, 5) + 0.000981429*pow(s,6);
  }
  else if (s >= 4.79 && s < 8.3)
  {
    A = -16885.8 + 10792.3 * s - 2264.77 * pow(s, 2) + 72.9756 * pow(s, 3) + 37.4671 * pow(s, 4) - 5.15719 * pow(s, 5) + 0.20543 * pow(s, 6);
    //B = -3283.75 + 3064.11*s - 1068.44*pow(s,2) + 164.844*pow(s,3) - 9.48152*pow(s,4);
    B =  -27.2777 +  8.03803* s +0.342604  * pow(s, 2) -0.0886274 * pow(s, 3) -0.01421 * pow(s, 4) -0.000251407 * pow(s, 5) + 0.000213174*pow(s,6) +2.9559e-05*pow(s,7) -3.45896e-06*pow(s,8);
  }
  else
  {
    return 0.;
  }

  double xs = A * exp(B * t) * (1 + 0.02 * exp(-6 * t));

  return xs;
}
