#ifndef __PHOTO_GENERATOR_H__
#define __PHOTO_GENERATOR_H__

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "generator/gcfGenerator.hh"
#include "nucleus/gcfNucleus.hh"
#include "cross_sections/photoCrossSection.hh"

class photoGenerator: public gcfGenerator
{
 public:
  photoGenerator(double E, gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand);
  ~photoGenerator();
  void generate_event(double &weight, int &meson_type, int &baryon_type, int &rec_type, TLorentzVector &vMeson_target, TLorentzVector &vBaryon_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target);
  
 private:
  photoCrossSection * myCS;
    
  double Ebeam;
  TVector3 vbeam;
  TLorentzVector vbeam_target;
  
};

#endif
