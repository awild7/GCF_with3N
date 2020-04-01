#ifndef __ELECTRO_GENERATOR_H__
#define __ELECTRO_GENERATOR_H__

#include "TVector3.h"
#include "TLorentzVector.h"

#include "generator/gcfGenerator.hh"

class gcfNucleus;
class photoCrossSection;
class TRandom3;

class electroGenerator: public gcfGenerator
{
 public:
  electroGenerator(double E, gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand);
  ~electroGenerator();
  void generate_event(double &weight, int &meson_type, int &baryon_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vMeson_target, TLorentzVector &vBaryon_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, double &r);
  bool parse_phase_space_file(char* phase_space);
  void set_nu_range(double low, double high);
  void set_QSq_range(double low, double high);
  void set_phik_range(double low, double high);
  void set_phik_range_deg(double low, double high);
  
 private:
  double photon_flux(double nu, double QSq);
  
  photoCrossSection * myCS;
    
  double Ebeam;
  TVector3 vbeam;
  TLorentzVector vbeam_target;

  double numin;
  double numax;
  double QSqmin;
  double QSqmax;
  double phikmin;
  double phikmax;

};

#endif
