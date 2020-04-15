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

  double Q = sqrt(QSq);
  
  return x*(sq(2./3.)*(myCTEQ.parton(1,x,Q) + myCTEQ.parton(-1,x,Q) + myCTEQ.parton(4,x,Q) + myCTEQ.parton(-4,x,Q)) +
	    sq(1./3.)*(myCTEQ.parton(2,x,Q) + myCTEQ.parton(-2,x,Q) + myCTEQ.parton(3,x,Q) + myCTEQ.parton(-3,x,Q) + myCTEQ.parton(5,x,Q) + myCTEQ.parton(-5,x,Q)));

}

double DISCrossSection::F2n(double x, double QSq)
{
  
  double Q = sqrt(QSq);
  
  return x*(sq(2./3.)*(myCTEQ.parton(2,x,Q) + myCTEQ.parton(-2,x,Q) + myCTEQ.parton(4,x,Q) + myCTEQ.parton(-4,x,Q)) +
	    sq(1./3.)*(myCTEQ.parton(1,x,Q) + myCTEQ.parton(-1,x,Q) + myCTEQ.parton(3,x,Q) + myCTEQ.parton(-3,x,Q) + myCTEQ.parton(5,x,Q) + myCTEQ.parton(-5,x,Q)));

}
