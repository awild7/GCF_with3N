#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "photoGenerator.hh"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "constants.hh"
#include "helpers.hh"
#include "cross_sections/photoCrossSection.hh"
#include "nucleus/gcfNucleus.hh"

using namespace std;

photoGenerator::photoGenerator(double E, gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{

  myCS = thisCS;
  
  Ebeam = E;
  vbeam.SetXYZ(0.,0.,Ebeam);
  vbeam_target.SetXYZT(0.,0.,Ebeam,Ebeam);

}

photoGenerator::~photoGenerator()
{
}

void photoGenerator::generate_event(double &weight, int &meson_type, int &baryon_type, int &rec_type, TLorentzVector &vMeson_target, TLorentzVector &vBaryon_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what kind of proton or neutron pair we are dealing with
  int lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  // Decide what branching ratio to generate (Currently only implemented gamma p -> pi+ n and equivalent neutron channel)
  double mMeson;
  double mBaryon;
  if (lead_type == pCode)
    {
      meson_type = pipCode;
      baryon_type = nCode;
    }
  else
    {
      meson_type = pimCode;
      baryon_type = pCode;
    }
  mMeson = mpip;
  mBaryon = mN;
  weight *= 1;
  
  // Determine mass of A-2 system
  double mAm2;
  if (lead_type == pCode and rec_type == pCode)
    mAm2 = mAmpp;
  else if (lead_type == nCode and rec_type == nCode)
    mAm2 = mAmnn;
  else
    mAm2 = mAmpn;

  TVector3 v1, vRec;

  decay_function(weight, lead_type, rec_type, v1, vRec);
  
  if (weight > 0.)
    {
      
      TVector3 vAm2 = - v1 - vRec;
      double EAm2 = sqrt(vAm2.Mag2() + sq(mAm2));
      vAm2_target.SetVect(vAm2);
      vAm2_target.SetT(EAm2);
      
      double Erec = sqrt(sq(mN) + vRec.Mag2());
      vRec_target.SetVect(vRec);
      vRec_target.SetT(Erec);
      
      double E1 = mA - EAm2 - Erec;
      TLorentzVector v1_target(v1,E1);

      t_scatter(weight, mMeson, mBaryon, vbeam_target, v1_target, vMeson_target, vBaryon_target);

      if (weight > 0.)
	{

	  double s = sq(mMeson) + sq(mBaryon) + 2.*vMeson_target.Dot(vBaryon_target);
	  double t = sq(mMeson) - 2.*vbeam_target.Dot(vMeson_target);
	  
	  // Calculate the flux factor on the cross section
	  double vgamma1 = vbeam_target.Dot(v1_target)/(Ebeam*E1);
	  
	  // Calculate the weight
	  weight *= vgamma1*myCS->sigma_pip_n(s,t); // Photoproduction cross section
	    
	}
  
    }
  
}

bool photoGenerator::parse_phase_space_file(char* phase_space)
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
      else
	{
	  cerr << "Invalid phase space parameter provided. Aborting...\n";
	  return false;
	}
    }
  
  return true;
  
}
