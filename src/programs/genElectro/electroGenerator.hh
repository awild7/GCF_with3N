#ifndef __ELECTRO_GENERATOR_H__
#define __ELECTRO_GENERATOR_H__

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "generator/gcfGenerator.hh"
#include "nucleus/gcfNucleus.hh"
#include "cross_sections/photoCrossSection.hh"

class electroGenerator: public gcfGenerator
{
 public:
  electroGenerator(double E, gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand);
  ~electroGenerator();
  void generate_event(double &weight, int &meson_type, int &baryon_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vMeson_target, TLorentzVector &vBaryon_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, double &r);
  
 private:
  double photon_flux(double nu, double QSq);
  
  photoCrossSection * myCS;
    
  double Ebeam;
  TVector3 vbeam;
  TLorentzVector vbeam_target;

};

#endif
