#ifndef __QE_GENERATOR_3N_H__
#define __QE_GENERATOR_3N_H__

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "generator/gcfGenerator.hh"
#include "nucleus/gcfNucleus.hh"
#include "cross_sections/eNCrossSection.hh"
#include <fstream>
#include <iostream>

#define ADIR ARRAY_DIR

class QEGenerator_3N
{
 public:
  QEGenerator_3N(double E, eNCrossSection * thisCS, int thisU, TRandom3 * thisRand);
  ~QEGenerator_3N();
  void generate_event(double &weight, int &N1_type, int &N2_type, int &N3_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &v2_target, TLorentzVector &v3_target, TLorentzVector &vAm3_target);
  void generate_event(double &weight, int &N1_type, int &N2_type, int &N3_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &v2_target, TLorentzVector &v3_target, TLorentzVector &vAm3_target, double &Estar);
  void set_theta_k_maxmin(double min, double max);
  double get_rho(double N2_Type, double N3_Type, double p_a, double p_b, double theta_ab);
  double get_rho_ptot_f1f2f3(double k_1, double k_2, double k_3);
  
 private:
  eNCrossSection * myCS;
    
  double Ebeam;
  TVector3 vbeam;
  TLorentzVector vbeam_target;
  TRandom3 * myRand;

  int u;
  TString uType;

  double sigCM;
  double phi_a_max;
  double phi_a_min;
  double cos_theta_a_max;
  double cos_theta_a_min;
  double phi_baz_max;
  double phi_baz_min;
  double mom_a_max;
  double mom_a_min;
  double mom_b_max;
  double mom_b_min;
  double theta_ab_max;
  double theta_ab_min;
  double cos_theta_k_max;
  double cos_theta_k_min;
  double phi_k_max;
  double phi_k_min;

  //Functions to get the 3D wavefunction
  static const int theta_bins = 61;
  const double theta_min = 0;
  const double theta_max = 180;
  const double theta_width = (theta_max - theta_min)/((double)theta_bins-1);

  static const int k_cm_bins = 101;
  const double k_cm_min = 0;
  const double k_cm_max = 10;
  const double k_cm_width = (k_cm_max - k_cm_min)/((double)k_cm_bins-1);

  static const int k_rel_bins = 101;
  const double k_rel_min = 0;
  const double k_rel_max = 10;
  const double k_rel_width = (k_rel_max - k_rel_min)/((double)k_rel_bins-1);

  double density_matrix_pp[theta_bins][k_cm_bins][k_rel_bins];
  //double density_matrix_pn[theta_bins][k_cm_bins][k_rel_bins];

  void fill_array();
  double interpolate(double matrix[k_rel_bins], double k_rel);
  double interpolate(double matrix[k_cm_bins][k_rel_bins], double k_cm, double k_rel);
  double interpolate(double matrix[theta_bins][k_cm_bins][k_rel_bins], double theta, double k_cm, double k_rel);
  
};

#endif
