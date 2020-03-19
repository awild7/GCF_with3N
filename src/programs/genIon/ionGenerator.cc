#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "ionGenerator.hh"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "constants.hh"
#include "helpers.hh"
#include "cross_sections/NNCrossSection.hh"
#include "nucleus/gcfNucleus.hh"

using namespace std;

ionGenerator::ionGenerator(double E, gcfNucleus * thisInfo, NNCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{

  myCS = thisCS;
  
  Ebeam = E;
  E1 = Ebeam + mN;
  p1 = sqrt(sq(E1) - sq(mN));
  v1.SetXYZ(0.,0.,p1);
  v1_lab.SetXYZT(0.,0.,p1,E1);

}

ionGenerator::~ionGenerator()
{
}

void ionGenerator::generate_event(double &weight, int &lead_type, int &rec_type, TVector3 &v3, TVector3 &v4, TVector3 &vRec, TVector3 &vAm2)
{
  
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what kind of proton or neutron pair we are dealing with
  lead_type = pCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 2.;
  
  // Determine mass of A-2 system
  double mAm2;
  if (lead_type == pCode and rec_type == pCode)
    mAm2 = mAmpp;
  else if (lead_type == nCode and rec_type == nCode)
    mAm2 = mAmnn;
  else
    mAm2 = mAmpn;

  TVector3 v2;

  decay_function(weight, lead_type, rec_type, v2, vRec);
  
  if (weight > 0.)
    {
      
      vAm2 = - v2 - vRec;
      double EAm2 = sqrt(vAm2.Mag2() + sq(mAm2));
      
      double Erec = sqrt(sq(mN) + vRec.Mag2());  
      double E2 = mA - EAm2 - Erec;
      
      TLorentzVector v2_lab(v2,E2);
      TLorentzVector v3_lab, v4_lab;

      t_scatter(weight, mN, mN, v1_lab, v2_lab, v3_lab, v4_lab);

      if (weight > 0.)
	{
	  
	  v3 = v3_lab.Vect();
	  v4 = v4_lab.Vect();

	  double s = 2.*sq(mN) + 2.*v3_lab.Dot(v4_lab);
	  double t = 2.*sq(mN) - 2.*v1_lab.Dot(v3_lab);
	  double u = 2.*sq(mN) - 2.*v1_lab.Dot(v4_lab);
	  
	  // Calculate the flux factor on the cross section
	  double v12 = sqrt(v1_lab.Dot(v2_lab) - sq(mN*mN))/(E1*E2);
	  double v1A = p1/E1;
	  
	  // Calculate the weight
	  weight *= v12/v1A*myCS->sigma_pp(s,max(t,u)); // Scattering cross section
	    
	}
  
    }
  
}
