#ifndef __PHOTO_CROSS_SECTIONS_H__
#define __PHOTO_CROSS_SECTIONS_H__

class photoCrossSection
{
 public:
  photoCrossSection();
  ~photoCrossSection();
  double sigma_pi0_p(double s, double t);
  double sigma_pip_n(double s, double cosThetaCM);
  double sigma_pim_p(double s, double cosThetaCM);
  double sigma_rho0_p_old(double s, double cosThetaCM);
  double sigma_rho0_p(double s, double t, double cosThetaCM);
  double sigma_omega_p_old(double s, double cosThetaCM);
  double sigma_omega_p(double s, double t, double cosThetaCM);
  double sigma_phi_p(double s, double t, double cosThetaCM);
  double sigma_phi_n(double s, double t, double cosThetaCM);
  double sigma_Jpsi_p(double s, double t);
  double sigma_Jpsi_p(double s, double t, double QSq);
  double R_Jpsi_p(double QSq);
  double sigma_deltapp_pim(double s, double cosThetaCM);

private:
  double dipole_F(double t, double tmin, double tmax);

};

#endif
