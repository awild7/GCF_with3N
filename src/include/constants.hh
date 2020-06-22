#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

// constants
const double mU = 0.9314941024;
const double GeVfm  =0.1973;
const double alpha = 0.0072973525664;
const double cmSqGeVSq = GeVfm*GeVfm*1.E-26;
const double nbGeVSq = cmSqGeVSq*1.E33;
const double Vud = 0.97417;

// particle masses
const double me = 0.000511;
const double mmu = 0.10566;

const double mN = 0.93892;
const double mDelta = 1.232;
const double mLambda = 1.115683;
const double mSigmap = 1.18937;
const double mSigma0 = 1.192642;
const double mSigmam = 1.197449;

const double mpip = 0.139570;
const double mpi0 = 0.1349766;
const double mrhop = 0.7754;
const double mrho0 = 0.77549;
const double mKp = 0.493677;
const double mK0 = 0.497648;
const double momega = 0.78265;
const double mphi = 1.019461;
const double mpsi = 3.096916;

const double mW = 80.379;

// nuclear masses
const double m_1H = mN;
const double m_2H = 2.01410178 * mU - me;
const double m_3H = 3.01604928199 * mU - me;
const double m_3He = 3.0160293 * mU - 2*me;
const double m_4He = 4.00260325415 * mU - 2*me;
const double m_6Li = 6.015122795 * mU - 3*me;
const double m_8Be = 8.00530510 * mU - 4*me;
const double m_10Be = 10.013534 * mU - 4*me;
const double m_10B = 10.0129370 * mU - 5*me;
const double m_11B = 11.0093054 * mU - 5*me;
const double m_10C = 10.016853 * mU - 6*me;
const double m_12C = 12. * mU - 6*me;
const double m_14N = 14.0030740048 * mU - 7*me;
const double m_16O = 15.99491461956 * mU - 8*me;

const double m_25Na = 24.9899540 * mU - 11*me;
const double m_25Mg = 24.98583696 * mU - 12*me;
const double m_25Al = 24.99042831 * mU - 13*me;
const double m_27Al = 26.98153841 * mU - 13*me;

const double m_54Cr = 53.9388804 * mU - 24*me;
const double m_54Mn = 53.9403589 * mU - 25*me;
const double m_54Fe = 53.9396090 * mU - 26*me;
const double m_56Fe = 55.9349363 * mU - 26*me;

const double m_206Hg = 205.977514 * mU - 80*me + 0.5682;
const double m_206Tl = 205.9761103 * mU - 81*me + 0.5682;
const double m_206Pb = 205.9744653 * mU - 82*me + 0.5682;
const double m_208Pb = 207.9766521 * mU - 82*me + 0.5682;

// nucleon codes
const int pCode = 2212;
const int nCode = 2112;
const int DelatppCode = 2224;
const int DeltapCode = 2214;
const int Delta0Code = 2114;
const int DeltamCode = 1114;
const int LambdaCode = 3122;
const int SigmapCode = 3222;
const int Sigma0Code = 3212;
const int SigmamCode = 3112;

const int pipCode = 211;
const int pi0Code = 111;
const int pimCode = -pipCode;
const int rhopCode = 213;
const int rho0Code = 113;
const int rhomCode = - rhopCode;
const int KpCode = 321;
const int K0Code = 311;
const int KmCode = -KpCode;
const int omegaCode = 223;
const int phiCode = 333;
const int psiCode = 443;

#endif
