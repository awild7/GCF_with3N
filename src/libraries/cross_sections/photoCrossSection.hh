#ifndef __PHOTO_CROSS_SECTIONS_H__
#define __PHOTO_CROSS_SECTIONS_H__

#include "TVector3.h"

class photoCrossSection
{
 public:
  photoCrossSection();
  ~photoCrossSection();
  double sigma_pi0_p(double s, double t);
  double sigma_pi0_n(double s, double t);
  double sigma_pip_n(double s, double t);
  double sigma_pim_p(double s, double t);
  
};

#endif
