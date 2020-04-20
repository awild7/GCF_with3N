#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "QEGenerator.hh"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "constants.hh"
#include "helpers.hh"
#include "cross_sections/eNCrossSection.hh"
#include "nucleus/gcfNucleus.hh"

using namespace std;

QEGenerator::QEGenerator(double E, gcfNucleus * thisInfo, eNCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{

  myCS = thisCS;

  Ebeam = E;
  vbeam.SetXYZ(0.,0.,Ebeam);
  vbeam_target.SetXYZT(0.,0.,Ebeam,Ebeam);

  QSqmin = 1.0;
  QSqmax = 5.0;
  phikmin = 0.;
  phikmax = 2.*M_PI;
  
}

QEGenerator::~QEGenerator()
{
}

void QEGenerator::generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what kind of proton or neutron pair we are dealing with
  lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  // Determine mass of A-2 system
  double mAm2;
  if (lead_type == pCode and rec_type == pCode)
    mAm2 = mAmpp;
  else if (lead_type == nCode and rec_type == nCode)
    mAm2 = mAmnn;
  else
    mAm2 = mAmpn;

  // Sample decay function
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
  double p1_minus = E1 - v1.Z();

  // Pick random electron scattering
  double QSqmax_kine = 2*Ebeam*p1_minus;
  if (QSqmax_kine < QSqmin)
    {
      weight=0.;
      return;
    }
  
  double QSq = QSqmin + (min(QSqmax,QSqmax_kine) - QSqmin)*myRand->Rndm();
  double phik = phikmin + (phikmax - phikmin)*myRand->Rndm();
  weight *= (min(QSqmax,QSqmax_kine) - QSqmin) * (phikmax - phikmin);
  
  // Calculate outgoing electron kinematics
  double k_minus = QSq/(2*Ebeam);
  double plead_minus = p1_minus - k_minus;
  if (plead_minus < 0.)
    {
      weight=0.;
      return;
    }
  
  double p1_plus = E1 + v1.Z();
  double virt = v1_target.Mag2() - sq(mN);
  
  double A = k_minus/p1_minus;
  double c = A*(p1_plus*k_minus - 2*Ebeam*plead_minus - virt);
  
  double delta_phi = phik - v1.Phi();
  double p1_perp = v1.Perp();
  double y = p1_perp*cos(delta_phi);
  double b = -2*A*y;
  double D = sq(b) - 4*c;
  if (D < 0.)
    {
      weight=0.;
      return;
    }
  
  double k_perp = (- b + sqrt(D))/2.;

  // May be two possible solutions
  if (sqrt(D) < -b)
    {
      weight *= 2.;
      if (myRand->Rndm() > 0.5)
	k_perp = (- b - sqrt(D))/2.;
    }

  double k_plus = sq(k_perp)/k_minus;
  double kz = (k_plus - k_minus)/2.;
  TVector3 vk(k_perp*cos(phik),k_perp*sin(phik),kz);
  double Ek = vk.Mag();
  vk_target.SetVect(vk);
  vk_target.SetT(Ek);

  vLead_target = v1_target + vbeam_target - vk_target;
  TVector3 vLead = vLead_target.Vect();
  double Elead = vLead_target.T();
  
  // Jacobian for delta function
  double J = 2.*Ebeam*Ek*fabs(1. - (v1.Z() + y * tan(vk.Theta()/2.) + Ebeam - Ek)/Elead);
  weight *= 1./J;

  // Calculate the weight
  weight *= myCS->sigma_eN(Ebeam, vk, vLead, (lead_type==pCode)); // eN cross section
  
}
