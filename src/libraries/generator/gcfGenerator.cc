#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "gcfGenerator.hh"

#include "TVector3.h"
#include "TRandom3.h"

#include "constants.hh"
#include "helpers.hh"
#include "nucleus/gcfNucleus.hh"

using namespace std;

gcfGenerator::gcfGenerator(gcfNucleus * thisInfo, TRandom3 * thisRand)
{

  myInfo = thisInfo;
  myRand = thisRand;
  
  Anum = myInfo->get_A();
  mA = myInfo->get_mA();
  mbar = myInfo->get_mbar();
  mAmpp = myInfo->get_mAmpp(); // this includes the effect of Estar
  mAmpn = myInfo->get_mAmpn();
  mAmnn = myInfo->get_mAmnn();
  sigCM = myInfo->get_sigmaCM();

  pRelmin = 0.2;
  pRelmax = 1.05;
  phiRelmin = 0.;
  phiRelmax = 2.*M_PI;
  thetaRelmin = 0.;
  thetaRelmax = M_PI;  
  cosThetaRelmin = -1.;
  cosThetaRelmax = 1.;

}

gcfGenerator::~gcfGenerator()
{
}

void gcfGenerator::set_phiRel_range(double low, double high)
{
  phiRelmin = low;
  phiRelmax = high;
}

void gcfGenerator::set_phiRel_range_deg(double low, double high)
{
  phiRelmin = low*M_PI/180.;
  phiRelmax = high*M_PI/180.;
}

void gcfGenerator::set_thetaRel_range(double low, double high)
{
  thetaRelmin = low;
  thetaRelmax = high;
  cosThetaRelmin = cos(thetaRelmax);
  cosThetaRelmax = cos(thetaRelmin);
}

void gcfGenerator::set_thetaRel_range_deg(double low, double high)
{
  thetaRelmin = low*M_PI/180.;
  thetaRelmax = high*M_PI/180.;
  cosThetaRelmin = cos(thetaRelmax);
  cosThetaRelmax = cos(thetaRelmin);
}

void gcfGenerator::set_pRel_range(double low, double high)
{
  pRelmin = low;
  pRelmax = high;
}

void gcfGenerator::set_pRel_cut(double new_cutoff)
{
  pRel_cut = new_cutoff;
}

void gcfGenerator::decay_function(double &weight, double &lcweight, int lead_type, int rec_type, TVector3 &vi, TVector3 &vRec)
{
  
  // Pick random CM motion
  TVector3 vCM(myRand->Gaus(0.,sigCM),myRand->Gaus(0.,sigCM),myRand->Gaus(0.,sigCM));
  
  // Pick random relative motion
  TVector3 vRel;
  double phiRel = phiRelmin + (phiRelmax-phiRelmin)*myRand->Rndm();
  double cosThetaRel = cosThetaRelmin + (cosThetaRelmax - cosThetaRelmin)*myRand->Rndm();
  double thetaRel = acos(cosThetaRel);
  double pRel_Mag = pRelmin + (pRelmax - pRelmin)*myRand->Rndm();
  vRel.SetMagThetaPhi(pRel_Mag,thetaRel,phiRel);

  // Factor universal functions into the weights
  weight *= (phiRelmax-phiRelmin) * (cosThetaRelmax - cosThetaRelmin) * (pRelmax - pRelmin) // Phase Space
    * sq(pRel_Mag)*myInfo->get_S(pRel_Mag,lead_type,rec_type)/pow(2.*M_PI,3); // Contacts
  lcweight *= (phiRelmax-phiRelmin) * (cosThetaRelmax - cosThetaRelmin) * (pRelmax - pRelmin)
    * 0.; // Lightcone weight not yet implemented

  // Do a safeguard cut
  if (pRel_Mag < pRel_cut)
    weight=0.;

  // Determine initial nucleon momenta
  vi = 0.5*vCM + vRel;
  vRec = 0.5*vCM - vRel;

}
