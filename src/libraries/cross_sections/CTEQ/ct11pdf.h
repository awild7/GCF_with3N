#ifndef CTPDF_H
#define CTPDF_H

//--------------------------------------------------------------
// C++ version of the CTEQ PDF (***only CTEQ6.6 or later***)
// interface by Jun Gao and Pavel Nadolsky on Nov 2013.
//        <jung@smu.edu or nadolsky@physics.smu.edu>
//
// There is a new class "cteqpdf" for the CTEQ PDFs. Each PDF is
// defined as, e.g., "cteqpdf ct10". Then user can read the PDF
// table by "ct10.setct11(pdsname)", where "pdsname" is the name
// of the corresponding table file (.pds). After that the PDFs can
// be called as usual "ct10.parton(IP, XX, QQ)", as well as the
// QCD coupling constant "ct10.alphas(QQ)". User can also call
// "ct10.pdfexit()" to release the memory of the large talbe files.
//
// The flavor assignment is as usual for CTEQ, meaning   
// Ip is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
//                    for (b, c, s, d, u, g, u_bar, ..., b_bar).
// 
// This inteface works for the CTEQ6.6 format of table files, as 
// well as the newer ct10, ct11 formats. But only for the latter
// ones the alphas function should be used, for which the table
// files include additional column for the alphas interpolation
// values. As for CTEQ6.6 format, user can read the alphas value
// at a specific scale "ct10.Qalfa" as "ct10.AalfQ". Then may
// use external subroutine for running.
//
// More information could be found in the demo file.
//-------------------------------------------------------------

#include <vector>

//--------------------------------------------

using namespace std;

class cteqpdf {

 public:

// name of pds file  
  string filepds;

// QCD parameters
  double AlfaQ, Qalfa, amass[6];
  double QINI, QMAX, XMIN;
  int ipk, Iorder, Nfl;

// initializing of the table
  void setct11(string);
  double parton (int, double, double);
  double alphas (double);
  void pdfexit () {  UPD.clear(); };

 private:

// commons
  static const int MXX=201, MXQ=40, MXF=5, MaxVal=4;
  static const int MXPQX = (MXF+1+MaxVal) * MXQ * MXX;
  int ipdsformat, N0, Nfmx, MxVal;
  int NX, NT, NG, Npts, Nblk, Isetch, ipdsset;
  double qv[MXQ+1], TV[MXQ+1], AlsCTEQ[MXQ+1], XV[MXX+1];
  double Dr, fl, aimass, fswitch, xvpow[MXX+1];
  double Alambda, dummy, qbase, qbase1, qbase2, aa;

  static const int nqvec = 4;
  static const int ientry = 0;
  const double OneP = 1.00001;
  const double xpow = 0.3;

  int JX, JQ;
  double X, Q;
  int J1, JLX, JLQ, JU, JM, Ip, jtmp;
  double fvec[5], fij[5], fx;
  double ss, const1, const2, const3, const4, const5, const6;
  double sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3;
  double ff, tf2, tf3, sf2, sf3, s12, s13, s34, s24, s1213, s2434;
  double h00, tmp1, tmp2, tdet, sdet, tmp, g1, g4;
  double svec1, svec2, svec3, svec4, tvec1, tvec2, tvec3, tvec4;

// vector of PDF tables
  vector<double> UPD;
  void POLINT4F(double *, double *, double &, double &);


};

#endif
