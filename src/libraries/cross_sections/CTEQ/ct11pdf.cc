#include "ct11pdf.h"
#include <cstdlib>

//------------------------------------------------------------------------------------
void cteqpdf::setct11(string fname)
{
     JX=0;
     JQ=0;
     X=-1;
     Q=-1;
     filepds=fname;
     string aline;
     Isetch=1;
     ipdsset=1;

     ifstream infile;
     infile.open(fname.c_str());

     if (!infile) {
	     cerr << "error: unable to open input file: "
	     << fname <<endl;
     }

     getline(infile, aline);
     getline(infile, aline);

// check the format of the table file 

     if (aline.substr(0, 11) == "  ipk, Ordr") {
     	    cout << "pds file " << fname << " read with Post-CT10 format..." <<endl;
	    ipdsformat=10;
     	    infile >> ipk >> Dr >> Qalfa >> AlfaQ;
	    for(int i=0; i<=5; i++) infile >> amass[i];
	    Iorder=(int)Dr;
            getline(infile, aline);
            getline(infile, aline);
	    	if (aline.substr(0, 7) == "  IMASS") {
			ipdsformat=11;
			infile >> aimass >> fswitch >> N0 >> N0 >> N0 >> Nfmx >> MxVal;
			Nfl=Nfmx;
            		getline(infile, aline);
		}
		else {
			infile >> N0 >> N0 >> N0 >> Nfmx >> MxVal;
            		getline(infile, aline);
		}

     }
     else {
     	    cout << "pds file " << fname << " read with CTEQ6.6 format..." <<endl;
	    ipdsformat=6;
            AlfaQ=0.118;
	    Qalfa=91.187;
     	    infile >> Dr >> fl >> Alambda;
	    for(int i=0; i<=5; i++) infile >> amass[i];
            Iorder = (int)Dr;
	    Nfl = (int)fl;
            getline(infile, aline);
            getline(infile, aline);
	    infile >> dummy >> dummy >> dummy >> Nfmx >> MxVal >>N0;
            getline(infile, aline);
     }

     getline(infile, aline);
     infile >> NX >> NT >> N0 >> NG >> N0;
     getline(infile, aline);

     if (NG > 0) {
	     for(int i=1; i<=(NG+1); i++) getline(infile, aline); 
     }
     getline(infile, aline);

     if (ipdsformat >= 11) {
	     infile >> QINI >> QMAX;
	     for(int i=0; i<=NT; i++) infile >> qv[i] >> TV[i] >> AlsCTEQ[i];
     }
     else {
	     infile >> QINI >> QMAX;
	     for(int i=0; i<=NT; i++) infile >> qv[i] >> TV[i];
     }
     getline(infile, aline);

// check that qBase is consistent with the definition of Tv(0:nQ) for 2 values of Qv

     qbase1=qv[1]/exp(exp(TV[1]));
     qbase2=qv[NT]/exp(exp(TV[NT]));
     if (fabs(qbase1-qbase2) > 1e-5) {
	     cout << "Readpds0: something wrong with qbase" <<endl;
	     cout << "qbase1, qbase2= " << qbase1 << " " <<qbase2 <<endl;
	     exit(0);
     }
     else{
	     qbase=(qbase1+qbase2)/2;
     }

     getline(infile, aline);
     infile >> XMIN >> aa;
     XV[0]=0;
     for (int i=1; i<=NX; i++) infile >> XV[i];
     Nblk=(NX+1)*(NT+1);
     Npts=Nblk*(Nfmx+1+MxVal);
     getline(infile, aline);
     getline(infile, aline);

// read the PDFs into a vector

     UPD.push_back(0);
     double num;
     int rch=0;

     while (infile >> num) {
    
     UPD.push_back(num);
     rch++;

     }

     if (rch != Npts) {
	     cout << "Wrong in the table!" << "\n" << "length not match!" << endl;
	     exit(0);
     }

     int j=Npts-1;

     infile.clear();
     infile.close();

}

//------------------------------------------------------------------------------------

void cteqpdf::POLINT4F(double *XA, double *YA, double &X, double &Y) {
     
     double H1, H2, H3, H4, W, DEN, C1, C2, C3, D1, D2, D3;
     double CC1, CC2, CD1, CD2, DD1, DD2, DC1, DC2;  

     H1=XA[0]-X;
     H2=XA[1]-X;
     H3=XA[2]-X;
     H4=XA[3]-X;

     W=YA[1]-YA[0];
     DEN=W/(H1-H2);
     D1=H2*DEN;
     C1=H1*DEN;

     W=YA[2]-YA[1];
     DEN=W/(H2-H3);
     D2=H3*DEN;
     C2=H2*DEN;

     W=YA[3]-YA[2];
     DEN=W/(H3-H4);
     D3=H4*DEN;
     C3=H3*DEN;

     W=C2-D1;
     DEN=W/(H1-H3);
     CD1=H3*DEN;
     CC1=H1*DEN;

     W=C3-D2;
     DEN=W/(H2-H4);
     CD2=H4*DEN;
     CC2=H2*DEN;

     W=CC2-CD1;
     DEN=W/(H1-H4);
     DD1=H4*DEN;
     DC1=H1*DEN;

     if ((H3+H4) < 0) {
        Y=YA[3]+D3+CD2+DD1;
     }
     else if ((H2+H3) < 0) {
        Y=YA[2]+D2+CD1+DC1;
     }
     else if ((H1+H2) < 0) {
        Y=YA[1]+C2+CD1+DC1;
     }
     else {
        Y=YA[0]+C1+CC1+DC1;
     }
}

//------------------------------------------------------------------------------------
double cteqpdf::parton (int IPRTN, double XX, double QQ)  {

//  Given the parton distribution function in the array U in
//  COMMON / PEVLDT / , this routine interpolates to find
//  the parton distribution at an arbitray point in x and q.

//check the x and Q range

      if (ipdsset != 1) { 
         cout<< "CT11Pdf: the PDF table was not initialized"<<endl;
	 exit(0);
      }

      if ((XX < 0.) || (XX > 1.)) {
         cout<< "X out of range in CT11Pdf: " << XX <<endl;
	 return 0.;
      }

      if (QQ < 0.3) {
         cout<< "Q out of range in CT11Pdf: "<<QQ<<endl;
	 return 0.;
      }

      if (abs(IPRTN) > Nfmx) {
//        print a warning for calling extra flavor
         cout<<"Warning: Iparton out of range in CT11Pdf! "<<endl;
         cout<<"Iparton, MxFlvN0: "<< IPRTN <<" "<<Nfmx<<endl;
	 return 0.;
      }

// store the pows used for interpolation on first call...
      if (Isetch == 1) {
         Isetch = 0;

         xvpow[0] = 0e0;
         for (int i=1; i<=NX; i++) xvpow[i] = pow(XV[i],xpow);
      }
      if ((Isetch == 1) || (XX != X) || (QQ != Q)) {

      X = XX;
      Q = QQ;
      tt = log(log(Q/qbase));

//      -------------    find lower end of interval containing x, i.e.,
//                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLX = -1;
      JU = NX+1;
      while (JU-JLX > 1) {
         JM = (JU+JLX) / 2;

         if (X >= XV[JM]) {
            JLX = JM;
	 }
         else {
            JU = JM;
	 }
      }
// working to here
//                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
//                           |---|---|---|...|---|-x-|---|...|---|---|
//                     x     0  Xmin               x                 1

      if (JLX <= -1) {
        cout << "Severe error: x <= 0 in PartonX11! x = " << X <<endl;
        exit(0);
      }
      else if (JLX == 0) {
        JX = 0;
      }
      else if (JLX <= (NX-2)) {

//                For interrior points, keep x in the middle, as shown above
        JX = JLX - 1;
      }
      else if ((JLX == NX-1) || (X < OneP)) {

//                  We tolerate a slight over-shoot of one (OneP=1.00001),
//              perhaps due to roundoff or whatever, but not more than that.
//                                      Keep at least 4 points >= Jx
        JX = JLX - 2;
      }
      else {
        cout << "Severe error: x > 1 in PartonX11! x = " << X <<endl;
        exit(0);
      }
//          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

//                       This is the variable to be interpolated in
      ss = pow(X, xpow);

      if ((JLX >= 2) && (JLX <= (NX-2))) {

//     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow[JX];
      svec2 = xvpow[JX+1];
      svec3 = xvpow[JX+2];
      svec4 = xvpow[JX+3];

      s12 = svec1 - svec2;
      s13 = svec1 - svec3;
      s23 = svec2 - svec3;
      s24 = svec2 - svec4;
      s34 = svec3 - svec4;

      sy2 = ss - svec2;
      sy3 = ss - svec3;

//     constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23;
      const2 = s12/s23;
      const3 = s34/s23;
      const4 = s24/s23;
      s1213 = s12 + s13;
      s2434 = s24 + s34;
      sdet = s12*s34 - s1213*s2434;
      tmp = sy2*sy3/sdet;
      const5 = (s34*sy2-s2434*sy3)*tmp/s12;
      const6 = (s1213*sy2-s12*sy3)*tmp/s34;

      }

//         --------------Now find lower end of interval containing Q, i.e.,
//                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLQ = -1;
      JU = NT+1;
      while (JU-JLQ > 1) {
         JM = (JU+JLQ) / 2;
         if (tt >= TV[JM]) {
            JLQ = JM;
	 }
         else {
            JU = JM;
         }
       }

      if (JLQ <= 0) {
         JQ = 0;
      }
      else if (JLQ <= NT-2) {
//                                  keep q in the middle, as shown above
         JQ = JLQ - 1;
      }
      else {
//                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        JQ = NT - 3;
      }

//                                   This is the interpolation variable in Q

      if ((JLQ >= 1) && (JLQ <= NT-2)) {
//                                        store the lattice points in t...
      tvec1 = TV[JQ];
      tvec2 = TV[JQ+1];
      tvec3 = TV[JQ+2];
      tvec4 = TV[JQ+3];

      t12 = tvec1 - tvec2;
      t13 = tvec1 - tvec3;
      t23 = tvec2 - tvec3;
      t24 = tvec2 - tvec4;
      t34 = tvec3 - tvec4;

      ty2 = tt - tvec2;
      ty3 = tt - tvec3;

      tmp1 = t12 + t13;
      tmp2 = t24 + t34;

      tdet = t12*t34 - tmp1*tmp2;

        }
      } // end of recalculations...

// get the pdf function values at the lattice points...

      if (IPRTN > MxVal) {
         Ip = - IPRTN;
      }
      else {
         Ip = IPRTN;
      }
      jtmp = ((Ip + Nfmx)*(NT+1)+(JQ-1))*(NX+1)+JX+1;

      for (int it = 1; it <= nqvec; it++) {
         J1  = jtmp + it*(NX+1);

         if (JX == 0) {
//                      For the first 4 x points, interpolate x^2*f(x,Q)
//                      This applies to the two lowest bins JLx = 0, 1
//            We can not put the JLx.eq.1 bin into the "interrior" section
//                           (as we do for q), since Upd(J1) is undefined.
         fij[1] = 0;
         fij[2] = UPD[J1+1] * XV[1]*XV[1];
         fij[3] = UPD[J1+2] * XV[2]*XV[2];
         fij[4] = UPD[J1+3] * XV[3]*XV[3];

//                 Use Polint which allows x to be anywhere w.r.t. the grid

         POLINT4F (&xvpow[0], &fij[1], ss, fx);

         if (X > 0)  {
		 fvec[it] =  fx / X/X;
	 }
	 }
//                                              Pdf is undefined for x.eq.0
         else if  (JLX == NX-1) {
//                                                This is the highest x bin:

         POLINT4F (&xvpow[NX-3], &UPD[J1], ss, fx);

         fvec[it] = fx;
	 }

         else {
//                       for all interior points, use Jon's in-line function
//                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
         sf2 = UPD[J1+1];
         sf3 = UPD[J1+2];

         g1 =  sf2*const1 - sf3*const2;
         g4 = -sf2*const3 + sf3*const4;

         fvec[it] = (const5*(UPD[J1]-g1)+const6*(UPD[J1+3]
	          -g4)+sf2*sy3 - sf3*sy2) / s23;

           }
      } //end Q loop
//                                   We now have the four values Fvec(1:4)
//     interpolate in t...

      if (JLQ <= 0) {
//                         1st Q-bin, as well as extrapolation to lower Q
        POLINT4F (&TV[0], &fvec[1], tt, ff);
      }

      else if (JLQ >= NT-1) {
//                         Last Q-bin, as well as extrapolation to higher Q
        POLINT4F (&TV[NT-3], &fvec[1], tt, ff);
      }

      else {
//                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
//       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
//                         the full range QV(0:Nt)  (in contrast to XV)
      tf2 = fvec[2];
      tf3 = fvec[3];

      g1 = ( tf2*t13 - tf3*t12) / t23;
      g4 = (-tf2*t34 + tf3*t24) / t23;

      h00 = ((t34*ty2-tmp2*ty3)*(fvec[1]-g1)/t12
         +  (tmp1*ty2-t12*ty3)*(fvec[4]-g4)/t34);

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23;
      }

      return ff;

} //end of parton

//------------------------------------------------------------------------------------
double  cteqpdf::alphas (double QQ) {

      int JLQ, JU, JM, JQ;
      double Q, tt, Alsout;
      double tvec1, tvec2, tvec3, tvec4;
      double tdet, ty2, ty3, tmp1, tmp2, t12, t13, t23, t24, t34;
       
      if (ipdsset != 1) { 
          cout << "CT11Alphas: the PDF table was not initialized" <<endl;
	  exit(0);
      }
      
      if (ipdsformat < 11) {
         cout << "-------------------------Warning!----------------------"<<endl; 
         cout << "*  CT11alphas: the table of QCD coupling values"<<endl;
         cout << "*  is not included with this version of the table file."<<endl;
         cout << "*  You can still compute the PDFs, but do not call"<<endl;
         cout << "*  the CT11alphas function for the interpolation of"<<endl;
	 cout << "*  alpha_s."<<endl;
         cout << "-------------------------------------------------------"<<endl; 
         exit(0);
         return AlfaQ;
      }
      Q = QQ;
      tt = log(log(Q/qbase));

//         --------------   Find lower end of interval containing Q, i.e.,
//                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLQ = -1;
      JU = NT+1;
      while ((JU-JLQ) > 1) {
        JM = (JU+JLQ) / 2;
        if (tt >= TV[JM]) {
            JLQ = JM;
	}
        else {
            JU = JM;
	}
       }

      if (JLQ <= 0) {
         JQ = 0;
      }
      else if (JLQ <= (NT-2)) {
//                                  keep q in the middle, as shown above
         JQ = JLQ - 1;
      }
      else {
//                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        JQ = NT - 3;
      }

//                                   This is the interpolation variable in Q

      if ((JLQ >= 1) && (JLQ <= (NT-2))) {
//                                        store the lattice points in t...

        tvec1 = TV[JQ];
        tvec2 = TV[JQ+1];
        tvec3 = TV[JQ+2];
        tvec4 = TV[JQ+3];
        
        t12 = tvec1 - tvec2;
        t13 = tvec1 - tvec3;
        t23 = tvec2 - tvec3;
        t24 = tvec2 - tvec4;
        t34 = tvec3 - tvec4;
        
        ty2 = tt - tvec2;
        ty3 = tt - tvec3;
        
        tmp1 = t12 + t13;
        tmp2 = t24 + t34;
      
        tdet = t12*t34 - tmp1*tmp2;
        
      }
      
      POLINT4F (&TV[JQ], &AlsCTEQ[JQ], tt, Alsout);
      
      return Alsout;

}// end of alphas
