#ifndef __GCF_GENERATOR_H__
#define __GCF_GENERATOR_H__

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "nucleus/gcfNucleus.hh"

class gcfGenerator
{
 public:
  gcfGenerator(gcfNucleus * thisInfo, TRandom3 * thisRand);
  ~gcfGenerator();
  void set_doRad(bool rad);
  void set_doCoul(bool coul);
  void set_phiRel_range(double low, double high);
  void set_phiRel_range_deg(double low, double high);
  void set_thetaRel_range(double low, double high);
  void set_thetaRel_range_deg(double low, double high);
  void set_pRel_range(double low, double high);
  void set_pRel_cut(double new_cutoff);
  void set_nu_range(double low, double high);
  void set_xB_range(double low, double high);
  void set_QSq_range(double low, double high);
  void set_phik_range(double low, double high);
  void set_phik_range_deg(double low, double high);

  void randomize_cutoff();

  bool parse_phase_space_file(char* phase_space);
  
  void decay_function(double &weight, int lead_type, int rec_type, TVector3 &vi, TVector3 &vRec);
  void decay_function_lc(double &weight, int lead_type, int rec_type, double &alphai, TVector2 &vi_perp, double &alphaRec, TVector2 &vRec_perp);
  
  protected:
  double get_mAm2(int lead_type, int rec_type);
  double get_mAm2(int lead_type, int rec_type, double &Estar);

  void t_scatter(double &weight, double m3, double m4, TLorentzVector v1, TLorentzVector v2, TLorentzVector &v3, TLorentzVector &v4);
  void t_scatter(double &weight, double m3, double m4, TLorentzVector v1, TLorentzVector v2, TLorentzVector &v3, TLorentzVector &v4, double &cosThetaCM);

  TVector3 radiateElectron(TVector3 ve);
  double deltaHard(double QSq);
  double radiationFactor(double Ebeam, double Ek, double QSq);

  double calcCoulombEnergy();
  void coulombCorrection(TLorentzVector &p, double deltaE);
  
  gcfNucleus * myInfo;
  TRandom3 * myRand;  

  int Anum;
  double mA;
  double mbar;
  double mAmpp;
  double mAmpn;
  double mAmnn;
  double sigCM;
  bool random_Estar;
  bool doRad;
  bool doCoul;

  double pRel_cut = 0.25;
  double pRel_cut_range = 0.05;
  
  double pRelmin;
  double pRelmax;
  double phiRelmin;
  double phiRelmax;
  double thetaRelmin;
  double thetaRelmax; 
  double cosThetaRelmin;
  double cosThetaRelmax;
  double alphaRelmin;
  double alphaRelmax;
  
  double numin;
  double numax;
  double xBmin;
  double xBmax;
  double QSqmin;
  double QSqmax;
  double phikmin;
  double phikmax;
  
};

#endif
