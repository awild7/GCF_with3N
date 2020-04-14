#ifndef __DIS_CROSS_SECTIONS_H__
#define __DIS_CROSS_SECTIONS_H__

#include "CTEQ/ct11pdf.h"

class DISCrossSection
{
 public:
  DISCrossSection(string pdsFile);
  ~DISCrossSection();
  double sigma_xQSq_DIS(double x, double y, double QSq, int nucleon_type);

 private:
  cteqpdf myCTEQ;
  double F1(double x, double QSq, int nucleon_type);
  double F2(double x, double QSq, int nucleon_type);
  double F2p(double x, double QSq);
  double F2n(double x, double QSq);
  
};

#endif
