#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "gcfGenerator.hh"

#include "TVector3.h"
#include "TVector2.h"
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
  alphaRelmin = 0.;
  alphaRelmax = 2.;

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

void gcfGenerator::decay_function(double &weight, int lead_type, int rec_type, TVector3 &vi, TVector3 &vRec)
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

  // Do a safeguard cut
  if (pRel_Mag < pRel_cut)
    weight=0.;

  // Determine initial nucleon momenta
  vi = 0.5*vCM + vRel;
  vRec = 0.5*vCM - vRel;

}

void gcfGenerator::decay_function_lc(double &weight, int lead_type, int rec_type, double &alphai, TVector2 &vi_perp, double &alphaRec, TVector2 &vRec_perp)
{
  
  // Pick random CM motion
  double alphaCM = myRand->Gaus(2.,sigCM/mbar);
  TVector2 vCM_perp(myRand->Gaus(0.,sigCM),myRand->Gaus(0.,sigCM));
  
  // Pick random relative motion
  double alphaRel = alphaRelmin + (alphaRelmax-alphaRelmin)*myRand->Rndm();
  double kmin = abs(1.-alphaRel)/sqrt(alphaRel*(2.-alphaRel));
  double k = max(pRelmin,kmin) + (pRelmax - max(pRelmin,kmin))*myRand->Rndm();
  double phiRel = phiRelmin + (phiRelmax-phiRelmin)*myRand->Rndm();
  double kperpSq = alphaRel*(2.-alphaRel)*(sq(k)+sq(mN)) - sq(mN);
  TVector2 vRel_perp;
  vRel_perp.SetMagPhi(sqrt(kperpSq),phiRel);

  // Factor universal functions into the weights
  weight *= (alphaRelmax-alphaRelmin) * (pRelmax - max(pRelmin,kmin)) * (phiRelmax-phiRelmin) // Phase Space
    *k*sqrt(sq(k)+sq(mN))*myInfo->get_S(k,lead_type,rec_type)/pow(2.*M_PI,3); // Contacts

  // Do a safeguard cut
  if (k < pRel_cut)
    weight=0.;

  // Determine initial nucleon momenta
  alphai = alphaRel*alphaCM/2.;
  vRec_perp = 0.5*vCM_perp + vRel_perp;
  alphaRec = alphaCM - alphai;
  vRec_perp = 0.5*vCM_perp - vRel_perp;

}


void gcfGenerator::t_scatter(double &weight, double m3, double m4, TLorentzVector v1, TLorentzVector v2, TLorentzVector &v3, TLorentzVector &v4)
  {

    TLorentzVector Z_lab = v1 + v2;
    double s = Z_lab.M2();
    if (s < sq(m3 + m4))
      {
	weight=0.;
	return;
      }
    double W_cm = sqrt(s);
    
    // Boost vectors to scattering CM frame
    TVector3 m = Z_lab.BoostVector();
    TLorentzVector v1_cm = v1;
    TLorentzVector v2_cm = v2;
    v1_cm.Boost(-m);
    v2_cm.Boost(-m);

    // Rotate to scattering along z-axis
    double rot_phi = v1_cm.Vect().Phi();
    double rot_theta = v1_cm.Vect().Theta();

    // Determine scattered energy and momentum
    double E3_cm = (s + sq(m3) - sq(m4))/(2.*W_cm);
    double E4_cm = W_cm - E3_cm;
    double p_cm = sqrt(sq(E3_cm) - sq(m3));

    // Define Jacobian between cm scattering angle and t
    double J = 2.*p_cm*v1_cm.Vect().Mag()/(2.*M_PI);

    // Pick random CM scattering angle
    double phi_cm = 2.*M_PI*myRand->Rndm();
    double cosTheta_cm = -1. + 2.*myRand->Rndm();
    double theta_cm = acos(cosTheta_cm);

    weight *= J * 4.*M_PI;

    // Set outgoind particle vetors
    TVector3 v_cm;
    v_cm.SetMagThetaPhi(p_cm,theta_cm,phi_cm);

    // Rotate back
    v_cm.RotateY(rot_theta);
    v_cm.RotateZ(rot_phi);
    TLorentzVector v3_cm(v_cm,E3_cm);
    TLorentzVector v4_cm(-v_cm,E4_cm);
    
    // Boost to lab system
    v3 = v3_cm;
    v4 = v4_cm;
    v3.Boost(m);
    v4.Boost(m);
    
  }
