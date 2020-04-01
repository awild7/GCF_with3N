#ifndef __NN_CROSS_SECTIONS_H__
#define __NN_CROSS_SECTIONS_H__

enum csParam {Panin};

class NNCrossSection
{
 public:
  NNCrossSection();
  NNCrossSection(csParam thisParam);
  ~NNCrossSection();
  double sigma_pp(double s, double t);
  double sigma_pp_Panin(double s, double t);

 private:
  csParam myParam;

};

#endif
