#ifndef __CCQE_GENERATOR_H__
#define __CCQE_GENERATOR_H__

#include "TH1D.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "generator/gcfGenerator.hh"
#include "nucleus/gcfNucleus.hh"
#include "cross_sections/eNCrossSection.hh"

class CCQEGenerator: public gcfGenerator
{
 public:
  CCQEGenerator(gcfNucleus * thisInfo, eNCrossSection * thisCS, TRandom3 * thisRand);
  ~CCQEGenerator();
  void generate_event(double &weight, double &Eneutrino, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target);
  
 private:
  eNCrossSection * myCS;

  TH1D * neutrinoSpectrum;
  
};

#endif
