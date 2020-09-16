#include <iostream>
#include <fstream>
#include "TVector2.h"
#include "gcfGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

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
  random_Estar = myInfo->get_Estar_randomization();
  doRad = false;

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

void gcfGenerator::set_doRad(bool rad)
{
  doRad = rad;
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

void gcfGenerator::set_nu_range(double low, double high)
{
  numin = low;
  numax = high;
}

void gcfGenerator::set_xB_range(double low, double high)
{
  xBmin = low;
  xBmax = high;
}

void gcfGenerator::set_QSq_range(double low, double high)
{
  QSqmin = low;
  QSqmax = high;
}

void gcfGenerator::set_phik_range(double low, double high)
{
  phikmin = low;
  phikmax = high;
}

void gcfGenerator::set_phik_range_deg(double low, double high)
{
  phikmin = low*M_PI/180.;
  phikmax = high*M_PI/180.;
}

void gcfGenerator::randomize_cutoff()
{
  pRel_cut = pRel_cut + (myRand->Uniform() - 0.5)*pRel_cut_range;
}

bool gcfGenerator::parse_phase_space_file(char* phase_space)
{

  ifstream ps_file(phase_space);
  string param;
  double low, high;
  while (ps_file >> param >> low >> high)
    {
      if (param == "phiRel" || param == "phirel")
	{
	  set_phiRel_range(low,high);
	}
      else if (param == "phiRel_deg" || param == "phirel_deg")
	{
	  set_phiRel_range_deg(low,high);
	}
      else if (param == "thetaRel" || param == "thetarel")
	{
	  set_thetaRel_range(low,high);
	}
      else if (param == "thetaRel_deg" || param == "thetarel_deg")
	{
	  set_thetaRel_range_deg(low,high);
	}
      else if (param == "pRel" || param == "prel")
	{
	  set_pRel_range(low,high);
	}
      else if (param == "nu" || param == "omega")
	{
	  set_nu_range(low,high);
	}
      else if (param == "x" || param == "xB")
	{
	  set_xB_range(low,high);
	}
      else if (param == "QSq" || param == "Qsq")
	{
	  set_QSq_range(low,high);
	}
      else if (param == "phik" || param == "phie")
	{
	  set_phik_range(low,high);
	}
      else if (param == "phik_deg" || param == "phie_deg")
	{
	  set_phik_range_deg(low,high);
	}
      else
	{
	  cerr << "Invalid phase space parameter provided. Aborting...\n";
	  return false;
	}
    }
  
  return true;
  
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

  if ((phiRelmax < phiRelmin) or (cosThetaRelmax < cosThetaRelmin) or (pRelmax < pRelmin))
    {
      weight = 0.;
      return;
    }
  
  // Factor universal functions into the weights
  weight *= (phiRelmax-phiRelmin) * (cosThetaRelmax - cosThetaRelmin) * (pRelmax - pRelmin) // Phase Space
    * sq(pRel_Mag)*myInfo->get_S(pRel_Mag,lead_type,rec_type)/pow(2.*M_PI,3); // Contacts

  // Do a safeguard cut
  if (pRel_Mag < pRel_cut)
    {
      weight=0.;
      return;
    }

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

  if ((alphaRelmax < alphaRelmin) or (pRelmax < max(pRelmin,kmin)) or (phiRelmax < phiRelmin))
    {
      weight = 0.;
      return;
    }

  // Factor universal functions into the weights
  weight *= (alphaRelmax - alphaRelmin) * (pRelmax - max(pRelmin,kmin)) * (phiRelmax-phiRelmin) // Phase Space
    *k*sqrt(sq(k)+sq(mN))*myInfo->get_S(k,lead_type,rec_type)/pow(2.*M_PI,3); // Contacts
  
  // Do a safeguard cut
  if (k < pRel_cut)
    {
      weight=0.;
      return;
    }

  // Determine initial nucleon momenta
  alphai = alphaRel*alphaCM/2.;
  vi_perp = alphai/alphaCM*vCM_perp + vRel_perp;
  alphaRec = alphaCM - alphai;
  vRec_perp = alphaRec/alphaCM*vCM_perp - vRel_perp;
  
}

double gcfGenerator::get_mAm2(int lead_type, int rec_type)
{
  double Estar;
  return get_mAm2(lead_type,rec_type,Estar);
}

double gcfGenerator::get_mAm2(int lead_type, int rec_type, double &Estar)
{
  if (random_Estar)
    {
      if (lead_type == pCode and rec_type == pCode)
	return myInfo->get_mAmpp_random(Estar,myRand);
      else if (lead_type == nCode and rec_type == nCode)
	return myInfo->get_mAmnn_random(Estar,myRand);
      else
	return myInfo->get_mAmpn_random(Estar,myRand);
    }
  else
    {
      if (lead_type == pCode and rec_type == pCode)
	return mAmpp;
      else if (lead_type == nCode and rec_type == nCode)
	return mAmnn;
      else
	return mAmpn;
    }
}

void gcfGenerator::t_scatter(double &weight, double m3, double m4, TLorentzVector v1, TLorentzVector v2, TLorentzVector &v3, TLorentzVector &v4)
{
  double cosThetaCM;
  t_scatter(weight, m3, m4, v1, v2, v3, v4, cosThetaCM);
}

void gcfGenerator::t_scatter(double &weight, double m3, double m4, TLorentzVector v1, TLorentzVector v2, TLorentzVector &v3, TLorentzVector &v4, double &cosThetaCM)
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
  cosThetaCM = -1. + 2.*myRand->Rndm();
  double theta_cm = acos(cosThetaCM);

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

TVector3 gcfGenerator::radiateElectron(TVector3 ve)
{
  double Ee = ve.Mag();
  double lambda_e = alpha/M_PI*(log(4*sq(Ee)/sq(me)) -1);
  double DeltaE = pow(myRand->Rndm(),1./lambda_e)*Ee;
  double Ee_rad = Ee - DeltaE;
  TVector3 ve_rad = ve;
  ve_rad.SetMag(Ee_rad);
  return ve_rad;
}

double gcfGenerator::deltaHard(double QSq)
{
  return 2.*alpha/M_PI * ( -13./12.*log(QSq/sq(me)) + 8./3.);
}

double gcfGenerator::radiationFactor(double Ebeam, double Ek, double QSq)
{
  double lambda_ei = alpha/M_PI*(log(4*sq(Ebeam)/sq(me)) -1);
  double lambda_ef = alpha/M_PI*(log(4*sq(Ek)/sq(me)) -1);

  return (1 - deltaHard(QSq)) * pow(Ebeam/sqrt(Ebeam*Ek),lambda_ei) * pow(Ek/sqrt(Ebeam*Ek),lambda_ef);
}
