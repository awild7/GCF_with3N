#ifndef __DIS_GENERATOR_H__
#define __DIS_GENERATOR_H__

#include "TVector3.h"
#include "TLorentzVector.h"

#include "generator/gcfGenerator.hh"

class gcfNucleus;
class DISCrossSection;
class TRandom3;

class DISGenerator: public gcfGenerator
{
 public:
  DISGenerator(double E, gcfNucleus * thisInfo, DISCrossSection * thisCS, TRandom3 * thisRand);
  ~DISGenerator();
  void generate_event(double &weight, int &lead_type, int &rec_type, int &ipart, TLorentzVector& vk_target, TLorentzVector &vq_target, TLorentzVector &v1_target_onshell, TLorentzVector &vHadron_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target);
  
 private:
  DISCrossSection * myCS;
    
  double Ebeam;
  TVector3 vbeam;
  TLorentzVector vbeam_target;

};

#endif
