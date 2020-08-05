#include "TVector3.h"
#include "photoGenerator.hh"
#include "spectra/defaultSpectrum.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

photoGenerator::photoGenerator(gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{
  
  myCS = thisCS;
  myReaction = pim;

  photonSpectrum = new TH1D("photonEnergy","photonEnergy",280,5.,12.);
  for (int i = 0; i < 280; i++) {
    photonSpectrum->SetBinContent(i+1,defaultSpectrum[i]);
  }

}

photoGenerator::photoGenerator(gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand, reaction thisReaction) : gcfGenerator(thisInfo, thisRand)
{
  
  myCS = thisCS;
  myReaction = thisReaction;

  photonSpectrum = new TH1D("photonEnergy","photonEnergy",280,5.,12.);
  for (int i = 0; i < 280; i++) {
    photonSpectrum->SetBinContent(i+1,defaultSpectrum[i]);
  }

}

photoGenerator::~photoGenerator()
{
}

void photoGenerator::generate_event(double &weight, double &Ephoton, int &meson_type, int &baryon_type, int &rec_type, TLorentzVector &vMeson_target, TLorentzVector &vBaryon_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{

  Ephoton = photonSpectrum->GetRandom();
  
  TLorentzVector vphoton_target(0.,0.,Ephoton,Ephoton);
  
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what branching ratio to generate.
  int lead_type;
  double mMesonMean;
  double gammaMeson;
  double mBaryonMean;
  double gammaBaryon;
  if (myReaction==pim)
    {
      lead_type = nCode;
      meson_type = pimCode;
      baryon_type = pCode;
      mMesonMean = mpip;
      gammaMeson = 0.;
      mBaryonMean = mN;
      gammaBaryon = 0.;
    }
  else if (myReaction==rho0)
    {
      lead_type = pCode; 
      meson_type = rho0Code;
      baryon_type = pCode;
      mMesonMean = mrho0;
      gammaMeson = gammarho0;
      mBaryonMean = mN;
      gammaBaryon = 0.;
    }

  double mMeson = fabs(myRand->BreitWigner(mMesonMean,gammaMeson));
  double mBaryon = fabs(myRand->BreitWigner(mBaryonMean,gammaBaryon));
  
  // Decide recoil nucleon type
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 2.;

  // Determine mass of A-2 system
  double mAm2 = get_mAm2(lead_type, rec_type);

  TVector3 v1, vRec;

  decay_function(weight, lead_type, rec_type, v1, vRec);
  
  if (weight <= 0.)
    return;
      
  TVector3 vAm2 = - v1 - vRec;
  double EAm2 = sqrt(vAm2.Mag2() + sq(mAm2));
  vAm2_target.SetVect(vAm2);
  vAm2_target.SetT(EAm2);
  
  double Erec = sqrt(sq(mN) + vRec.Mag2());
  vRec_target.SetVect(vRec);
  vRec_target.SetT(Erec);
  
  double E1 = mA - EAm2 - Erec;
  TLorentzVector v1_target(v1,E1);

  double cosThetaCM;
  t_scatter(weight, mMeson, mBaryon, vphoton_target, v1_target, vMeson_target, vBaryon_target, cosThetaCM);
  
  if (weight <= 0.)
    return;
  
  double s = sq(mMeson) + sq(mBaryon) + 2.*vMeson_target.Dot(vBaryon_target);
  
  // Calculate the flux factor on the cross section
  double vgamma1 = vphoton_target.Dot(v1_target)/(Ephoton*E1);
  
  // Calculate the weight
  double thisCS=0;
  if (myReaction==pim)
    thisCS=myCS->sigma_pip_n(s,cosThetaCM);
  else if (myReaction==rho0)
    thisCS=myCS->sigma_rho0_p(s,cosThetaCM);  

  weight *= vgamma1*thisCS; // Photoproduction cross section
  
}
