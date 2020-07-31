#ifndef __PHOTO_CROSS_SECTIONS_H__
#define __PHOTO_CROSS_SECTIONS_H__

class photoCrossSection
{
 public:
  photoCrossSection();
  ~photoCrossSection();
  double sigma_pi0_p(double s, double t);
  double sigma_pi0_n(double s, double t);
  double sigma_pip_n(double s, double cosThetaCM);
  double sigma_pim_p(double s, double cosThetaCM);
  double sigma_psi_p(double s, double t);
  double sigma_psi_p(double s, double t, double QSq);
  double R_psi_p(double QSq);

private:
  double dipole_F(double t, double tmin, double tmax);
  
};

#endif
