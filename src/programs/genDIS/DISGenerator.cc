#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "DISGenerator.hh"

#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "constants.hh"
#include "helpers.hh"
#include "cross_sections/DISCrossSection.hh"
#include "nucleus/gcfNucleus.hh"

using namespace std;

DISGenerator::DISGenerator(double E, gcfNucleus * thisInfo, DISCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{

  myCS = thisCS;

  Ebeam = E;
  vbeam.SetXYZ(0.,0.,Ebeam);
  vbeam_target.SetXYZT(0.,0.,Ebeam,Ebeam);

  xBmin = 0.3;
  xBmax = 0.7;
  QSqmin = 2.5;
  QSqmax = 100.0;
  phikmin = 0.;
  phikmax = 2.*M_PI;
  
}

DISGenerator::~DISGenerator()
{
}

void DISGenerator::generate_event(double &weight, int &lead_type, int &rec_type, int &ipart, TLorentzVector& vk_target, TLorentzVector &vq_target, TLorentzVector &v1_target_onshell, TLorentzVector &vHadron_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what kind of proton or neutron pair we are dealing with
  lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  // Determine mass of A-2 system
  double mAm2 = get_mAm2(lead_type, rec_type);

  double alpha1, alphaRec;
  TVector2 v1_perp, vRec_perp;
  decay_function_lc(weight, lead_type, rec_type, alpha1, v1_perp, alphaRec, vRec_perp);
  
  if (weight <= 0.)
    return;

  double alphaCM = alpha1 + alphaRec;
  double alphaAm2 = Anum - alphaCM;
  TVector2 vAm2_perp = -1.*v1_perp - vRec_perp;
  
  if (xBmin > alpha1)
    {
      weight=0.;
      return;
    }
  
  // Pick random scattering quantities
  double xB = xBmin + (min(xBmax,alpha1) - xBmin)*myRand->Rndm();
  double QSq = QSqmin + (QSqmax - QSqmin)*myRand->Rndm();
  double phik = phikmin + (phikmax - phikmin)*myRand->Rndm();
  weight *= (min(xBmax,alpha1) - xBmin) * (QSqmax - QSqmin) * (phikmax - phikmin);
  
  double nu = QSq/(2.*mN*xB);
  double pe_Mag = Ebeam - nu;
  double cosThetak = 1. - QSq/(2.*Ebeam*pe_Mag);
  
  if (fabs(cosThetak) > 1.)
    {
      weight=0.;
      return;
    }
  
  double thetak = acos(cosThetak);
  TVector3 vk;
  vk.SetMagThetaPhi(pe_Mag,thetak,phik);
  vk_target.SetVect(vk);
  vk_target.SetT(pe_Mag);
  
  vq_target = vbeam_target - vk_target;
  double rot_phi = vq_target.Vect().Phi();
  double rot_theta = vq_target.Vect().Theta();
  
  // Constructing 4-momenta
  double pAm2_plus = mbar*alphaAm2;
  double pAm2_minus = (vAm2_perp.Mod2() + sq(mAm2))/pAm2_plus;
  double EAm2 = 0.5*(pAm2_plus + pAm2_minus);
  double pAm2z = 0.5*(pAm2_plus - pAm2_minus);
  vAm2_target.SetXYZT(vAm2_perp.X(),vAm2_perp.Y(),pAm2z,EAm2);
  vAm2_target.RotateY(rot_theta);
  vAm2_target.RotateZ(rot_phi);
  
  double pRec_plus = mbar*alphaRec;
  double pRec_minus = (vRec_perp.Mod2() + sq(mN))/pRec_plus;
  double ERec = 0.5*(pRec_plus + pRec_minus);
  double pRecz = 0.5*(pRec_plus - pRec_minus);
  vRec_target.SetXYZT(vRec_perp.X(),vRec_perp.Y(),pRecz,ERec);
  vRec_target.RotateY(rot_theta);
  vRec_target.RotateZ(rot_phi);
  
  double p1_plus = mbar*alpha1;
  double p1_minus_onshell = (v1_perp.Mod2() + sq(mN))/p1_plus;
  double E1_onshell = 0.5*(p1_plus + p1_minus_onshell);
  double p1z_onshell = 0.5*(p1_plus - p1_minus_onshell);
  v1_target_onshell.SetXYZT(v1_perp.X(),v1_perp.Y(),p1z_onshell,E1_onshell);
  v1_target_onshell.RotateY(rot_theta);
  v1_target_onshell.RotateZ(rot_phi);
  
  double p1_minus = mA - pAm2_plus - pRec_plus;
  double E1 = 0.5*(p1_plus + p1_minus);
  double p1z = 0.5*(p1_plus - p1_minus);
  TLorentzVector v1_target(v1_perp.X(),v1_perp.Y(),p1z,E1);
  v1_target.RotateY(rot_theta);
  v1_target.RotateZ(rot_phi);
  
  vHadron_target = v1_target + vq_target;
  
  double y = (Ebeam - pe_Mag)/Ebeam;
  
  // Calculate the weight
  weight *= myCS->sigma_xQSq_DIS(xB/alpha1,y,QSq,lead_type)/(2.*M_PI); // DIS cross section
  
  // Determine struck parton
  ipart = myCS->getParton(xB/alpha1,QSq,lead_type,gRandom->Rndm());
	  
  
}
