#ifndef __ION_GENERATOR_H__
#define __ION_GENERATOR_H__

#include "TVector3.h"
#include "TLorentzVector.h"

#include "generator/gcfGenerator.hh"

class gcfNucleus;
class NNCrossSection;
class TRandom3;

class ionGenerator: public gcfGenerator
{
 public:
  ionGenerator(double E, gcfNucleus * thisInfo, NNCrossSection * thisCS, TRandom3 * thisRand);
  ~ionGenerator();
  void generate_event(double &weight, int &lead_type, int &rec_type, TVector3 &v3, TVector3 &v4, TVector3 &vRec, TVector3 &vAm2);
  
 private:
  NNCrossSection * myCS;
    
  double Ebeam;
  double E1;
  double p1;
  TVector3 v1;
  TLorentzVector v1_lab;
  
};

#endif
