#ifndef __DIS_CROSS_SECTIONS_H__
#define __DIS_CROSS_SECTIONS_H__

#include "CTEQ/ct11pdf.h"

class DISCrossSection
{
 public:
  DISCrossSection(string pdsFile);
  ~DISCrossSection();
  double F1(double x, double QSq, int nucleon_type);
  double F2(double x, double QSq, int nucleon_type);
  int getParton(double x, double QSq, int nucleon_type, double r);

 private:
  cteqpdf myCTEQ;
  double F2p(double x, double QSq);
  double F2n(double x, double QSq);

  double u(double x, double QSq);
  double d(double x, double QSq);
  double s(double x, double QSq);
  double c(double x, double QSq);
  double b(double x, double QSq);
  double ubar(double x, double QSq);
  double dbar(double x, double QSq);
  double sbar(double x, double QSq);
  double cbar(double x, double QSq);
  double bbar(double x, double QSq);
  int getPartonp(double x, double QSq, double r);
  int getPartonn(double x, double QSq, double r);
  
};

#endif
