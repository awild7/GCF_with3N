#ifndef __PHOTO_GENERATOR_H__
#define __PHOTO_GENERATOR_H__

#include "TH1D.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "generator/gcfGenerator.hh"
#include "nucleus/gcfNucleus.hh"
#include "cross_sections/photoCrossSection.hh"

enum spectrum {diamond,amorphous};
enum reaction {pim,proton_rho0,proton_rhoMinus,proton_omega,proton_phi,neutron_phi,DeltaPlusPlus_piMinus,DeltaPlus_piMinus};

class photoGenerator: public gcfGenerator
{
public:
  photoGenerator(gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand);
  photoGenerator(gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand, spectrum thisSpectrum, reaction thisReaction);
  ~photoGenerator();
  void generate_event(double &weight, double &Ephoton, int &meson_type, double &mMeson, int &baryon_type, double &mBaryon, int &rec_type, TLorentzVector &vMeson_target, TLorentzVector &vBaryon_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target);
  void setfixedE(double newfixedE);
  void print_beam_info();

private:
  photoCrossSection * myCS;
  bool usingfixedE;
  double fixedE;
  TH1D * photonSpectrum;

  void fillDiamond();
  void fillAmorphous();

  reaction myReaction;

};

#endif
