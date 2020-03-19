#ifndef __GCF_GENERATOR_H__
#define __GCF_GENERATOR_H__

#include "TVector3.h"

class gcfNucleus;
class TRandom3;

class gcfGenerator
{
 public:
  gcfGenerator(gcfNucleus * thisInfo, TRandom3 * thisRand);
  ~gcfGenerator();
  void set_phiRel_range(double low, double high);
  void set_phiRel_range_deg(double low, double high);
  void set_thetaRel_range(double low, double high);
  void set_thetaRel_range_deg(double low, double high);
  void set_pRel_range(double low, double high);
  void set_pRel_cut(double new_cutoff);
  
 private:
  void decay_function(double &weight, double &lcweight, int lead_type, int rec_type, TVector3 &vi, TVector3 &vRec);
  gcfNucleus * myInfo;
  TRandom3 * myRand;  

  int Anum;
  double mA;
  double mbar;
  double mAmpp;
  double mAmpn;
  double mAmnn;
  double sigCM;

  double pRel_cut = 0.25;
  
  double pRelmin;
  double pRelmax;
  double phiRelmin;
  double phiRelmax;
  double thetaRelmin;
  double thetaRelmax; 
  double cosThetaRelmin;
  double cosThetaRelmax;
  
};

#endif
