#ifndef __QE_GENERATOR_H__
#define __QE_GENERATOR_H__

#include "TVector3.h"
#include "TLorentzVector.h"

#include "generator/gcfGenerator.hh"

class gcfNucleus;
class eNCrossSection;
class TRandom3;

class QEGenerator: public gcfGenerator
{
 public:
  QEGenerator(double E, gcfNucleus * thisInfo, eNCrossSection * thisCS, TRandom3 * thisRand);
  ~QEGenerator();
  void generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target);
  void generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, double &Estar);
  
 private:
  eNCrossSection * myCS;
    
  double Ebeam;
  TVector3 vbeam;
  TLorentzVector vbeam_target;

};

#endif
