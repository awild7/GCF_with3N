#include <iostream>
#include <cmath>
#include <cstdlib>

#include "DISCrossSection.hh"
#include "helpers.hh"
#include "constants.hh"

DISCrossSection::DISCrossSection(string pdsFile)
{
  myCTEQ.setct11(pdsFile);
}

DISCrossSection::~DISCrossSection(){}

double DISCrossSection::sigma_xQSq_DIS(double x, double y, double QSq, int nucleon_type)
{
  
  return nbGeVSq*(4*M_PI*sq(alpha))/(sq(QSq))*((1.-y-sq(mN*y)/QSq)*F2(x,QSq,nucleon_type)/x + sq(y)*F1(x,QSq,nucleon_type));

}

int DISCrossSection::getParton(double x, double QSq, int nucleon_type, double r)
{
  if (nucleon_type == pCode)
    return getPartonp(x,QSq,r);
  else if (nucleon_type == nCode)
    return getPartonn(x,QSq,r);
  else
    {
      std::cerr << "Invalid nucleon type in getParton call. Check and fix.\n\n\n";
      return 0.;
    }
}

double DISCrossSection::F1(double x, double QSq, int nucleon_type)
{
  
  return F2(x,QSq,nucleon_type)/(2.*x);

}

double DISCrossSection::F2(double x, double QSq, int nucleon_type)
{
  if (nucleon_type == pCode)
    return F2p(x,QSq);
  else if (nucleon_type == nCode)
    return F2n(x,QSq);
  else
    {
      std::cerr << "Invalid nucleon type in Structure Function call. Check and fix.\n\n\n";
      return 0.;
    }

}

double DISCrossSection::F2p(double x, double QSq)
{
  
  return x*(sq(2./3.)*(u(x,QSq) + ubar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	    sq(1./3.)*(d(x,QSq) + dbar(x,QSq) + s(x,QSq) + sbar(x,QSq) + b(x,QSq) + bbar(x,QSq)));

}

double DISCrossSection::F2n(double x, double QSq)
{
  
  return x*(sq(2./3.)*(d(x,QSq) + dbar(x,QSq) + c(x,QSq) + cbar(x,QSq)) + 
	    sq(1./3.)*(u(x,QSq) + ubar(x,QSq) + s(x,QSq) + sbar(x,QSq) + b(x,QSq) + bbar(x,QSq)));


}

double DISCrossSection::u(double x, double QSq)
{
  return myCTEQ.parton(1,x,sqrt(QSq));
}

double DISCrossSection::d(double x, double QSq)
{
  return myCTEQ.parton(2,x,sqrt(QSq));
}

double DISCrossSection::s(double x, double QSq)
{
  return myCTEQ.parton(3,x,sqrt(QSq));
}

double DISCrossSection::c(double x, double QSq)
{
  return myCTEQ.parton(4,x,sqrt(QSq));
}

double DISCrossSection::b(double x, double QSq)
{
  return myCTEQ.parton(5,x,sqrt(QSq));
}

double DISCrossSection::ubar(double x, double QSq)
{
  return myCTEQ.parton(-1,x,sqrt(QSq));
}

double DISCrossSection::dbar(double x, double QSq)
{
  return myCTEQ.parton(-2,x,sqrt(QSq));
}

double DISCrossSection::sbar(double x, double QSq)
{
  return myCTEQ.parton(-3,x,sqrt(QSq));
}

double DISCrossSection::cbar(double x, double QSq)
{
  return myCTEQ.parton(-4,x,sqrt(QSq));
}

double DISCrossSection::bbar(double x, double QSq)
{
  return myCTEQ.parton(-5,x,sqrt(QSq));
}

int DISCrossSection::getPartonp(double x, double QSq, double r)
{
  double R = r*(sq(2./3.)*(u(x,QSq) + ubar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
		sq(1./3.)*(d(x,QSq) + dbar(x,QSq) + s(x,QSq) + sbar(x,QSq) + b(x,QSq) + bbar(x,QSq)));

  if (R < sq(2./3.)*(u(x,QSq)))
    return 1;
  else if (R < sq(2./3.)*(u(x,QSq) + ubar(x,QSq)))
    return -1;
  else if (R < sq(2./3.)*(u(x,QSq) + ubar(x,QSq) + c(x,QSq)))
    return 4;
  else if (R < sq(2./3.)*(u(x,QSq) + ubar(x,QSq) + c(x,QSq) + cbar(x,QSq)))
    return -4;
  else if (R < sq(2./3.)*(u(x,QSq) + ubar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(d(x,QSq)))
    return 2;
  else if (R < sq(2./3.)*(u(x,QSq) + ubar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(d(x,QSq) + dbar(x,QSq)))
    return -2;
  else if (R < sq(2./3.)*(u(x,QSq) + ubar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(d(x,QSq) + dbar(x,QSq) + s(x,QSq)))
    return 3;
  else if (R < sq(2./3.)*(u(x,QSq) + ubar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(d(x,QSq) + dbar(x,QSq) + s(x,QSq) + sbar(x,QSq)))
    return -3;
  else if (R < sq(2./3.)*(u(x,QSq) + ubar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(d(x,QSq) + dbar(x,QSq) + s(x,QSq) + sbar(x,QSq) + b(x,QSq)))
    return 5;
  else
    return -5;

  return 0;
}

int DISCrossSection::getPartonn(double x, double QSq, double r)
{
  double R = r*(sq(2./3.)*(d(x,QSq) + dbar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
		sq(1./3.)*(u(x,QSq) + ubar(x,QSq) + s(x,QSq) + sbar(x,QSq) + b(x,QSq) + bbar(x,QSq)));

  if (R < sq(2./3.)*(d(x,QSq)))
    return 1;
  else if (R < sq(2./3.)*(d(x,QSq) + dbar(x,QSq)))
    return -1;
  else if (R < sq(2./3.)*(d(x,QSq) + dbar(x,QSq) + c(x,QSq)))
    return 4;
  else if (R < sq(2./3.)*(d(x,QSq) + dbar(x,QSq) + c(x,QSq) + cbar(x,QSq)))
    return -4;
  else if (R < sq(2./3.)*(d(x,QSq) + dbar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(u(x,QSq)))
    return 2;
  else if (R < sq(2./3.)*(d(x,QSq) + dbar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(u(x,QSq) + ubar(x,QSq)))
    return -2;
  else if (R < sq(2./3.)*(d(x,QSq) + dbar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(u(x,QSq) + ubar(x,QSq) + s(x,QSq)))
    return 3;
  else if (R < sq(2./3.)*(d(x,QSq) + dbar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(u(x,QSq) + ubar(x,QSq) + s(x,QSq) + sbar(x,QSq)))
    return -3;
  else if (R < sq(2./3.)*(d(x,QSq) + dbar(x,QSq) + c(x,QSq) + cbar(x,QSq)) +
	   sq(1./3.)*(u(x,QSq) + ubar(x,QSq) + s(x,QSq) + sbar(x,QSq) + b(x,QSq)))
    return 5;
  else
    return -5;

  return 0;
}
