#ifndef __PHOTO_GENERATOR_H__
#define __PHOTO_GENERATOR_H__

#include "TH1D.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "generator/gcfGenerator.hh"
#include "nucleus/gcfNucleus.hh"
#include "cross_sections/photoCrossSection.hh"


enum Reaction {pi,rho};
class photoGenerator: public gcfGenerator
{
public:
  photoGenerator(gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand, Reaction thisreaction);
  ~photoGenerator();
  void generate_event(double &weight, double &Ephoton, int &meson_type, int &baryon_type, int &rec_type, TLorentzVector &vMeson_target, TLorentzVector &vBaryon_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target);
  
private:
  photoCrossSection * myCS;

  TH1D * photonSpectrum;

Reaction myreacion

};

#endif
