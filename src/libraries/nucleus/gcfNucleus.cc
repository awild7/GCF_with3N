#include <iostream>
#include "gcfNucleus.hh"
#include "constants.hh"

gcfNucleus::gcfNucleus(int thisZ, int thisN, char* uType)
{
  Estar = 0;
  sigmaE = 0;
  random_Estar = false;
  set_Nucleus(thisZ, thisN);
  mySRC = new gcfSRC(thisZ,thisN,uType);
}  

gcfNucleus::gcfNucleus(int thisZ, int thisN, NNModel uType)
{
  Estar = 0;
  sigmaE = 0;
  random_Estar = false;
  set_Nucleus(thisZ, thisN);
  mySRC = new gcfSRC(thisZ,thisN,uType);
}  

gcfNucleus::gcfNucleus(int thisZ, int thisN, gcfSRC * thisSRC)
{
  Estar = 0;
  sigmaE = 0;
  random_Estar = false;
  set_Nucleus(thisZ, thisN);
  mySRC = thisSRC;
}  

gcfNucleus::~gcfNucleus()
{
}

void gcfNucleus::set_Nucleus(int thisZ, int thisN){
  Z = thisZ;
  N = thisN;
  A = Z + N;
  
  if ((Z==1) && (N==1))
    {
      sigmaCM=0.;
      d_sigmaCM=0.;
      mA=m_2H;
      mAmpp=0.;
      mAmpn=0.;
      mAmnn=0.;
      Estar_max = 0.0;
    }
  else if ((Z==2) && (N==1))
    {
      sigmaCM=0.1;
      d_sigmaCM=0.02;
      mA=m_3He;
      mAmpp=mN;
      mAmpn=mN;
      mAmnn=mN;
      Estar_max = 0.0;
    }  
 else if ((Z==1) && (N==2))
    {
      sigmaCM=0.1;
      d_sigmaCM=0.02;
      mA=m_3H;
      mAmpp=mN;
      mAmpn=mN;
      mAmnn=mN;
      Estar_max = 0.0;
    }
  else if ((Z==2) && (N==2))
    {
      sigmaCM=0.1;
      d_sigmaCM=0.02;
      mA=m_4He;
      mAmpp=2*mN;
      mAmpn=m_2H;
      mAmnn=2*mN;
      Estar_max = 0.010;      
    }
  else if ((Z==6) && (N==6))
    {
      sigmaCM=0.15;
      d_sigmaCM=0.02;
      mA=m_12C;
      mAmpp=m_10Be;
      mAmpn=m_10B;
      mAmnn=m_10C;
      Estar_max = 0.030;
    }
  else if ((Z==13) && (N==14))
    {
      sigmaCM=0.15;
      d_sigmaCM=0.02;
      mA=m_27Al;
      mAmpp=m_25Na;
      mAmpn=m_25Mg;
      mAmnn=m_25Al;
      Estar_max = 0.040;
    }
  else if ((Z==18) && (N==22))
    {
      sigmaCM=0.15;
      d_sigmaCM=0.02;
      mA=m_40Ar;
      mAmpp=m_38S;
      mAmpn=m_38Cl;
      mAmnn=m_38Ar;
      Estar_max = 0.050;
    }
  else if ((Z==20) && (N==20))
    {
      sigmaCM=0.15;
      d_sigmaCM=0.02;
      mA=m_40Ca;
      mAmpp=m_38Ar;
      mAmpn=m_38K;
      mAmnn=m_38Ca;
      Estar_max = 0.050;
    }
  else if ((Z==26) && (N==30))
    {
      sigmaCM=0.15;
      d_sigmaCM=0.02;
      mA=m_56Fe;
      mAmpp=m_54Cr;
      mAmpn=m_54Mn;
      mAmnn=m_54Fe;
      Estar_max = 0.050;
    }
  else if ((Z==82) && (N==126))
    {
      sigmaCM=0.15;
      d_sigmaCM=0.02;
      mA=m_208Pb;
      mAmpp=m_206Hg;
      mAmpn=m_206Tl;
      mAmnn=m_206Pb;
      Estar_max = 0.100;
    }
  else
    {
      std::cerr << "You selected a nucleus with Z=" << Z << " and with N=" << N << "\n"
	   << " which is not in the library. Aborting...\n";
      exit(-2);
    }

}

void gcfNucleus::randomize(TRandom3* myRand)
{
  randomize_sigmaCM(myRand);
  randomize_Estar(myRand);
  mySRC->randomize_Contacts(myRand);
}

void gcfNucleus::randomize_sigmaCM(TRandom3* myRand)
{
  sigmaCM += myRand->Gaus(0.,d_sigmaCM);
}

void gcfNucleus::randomize_Estar(TRandom3* myRand)
{
  Estar = myRand->Uniform()*Estar_max;
}

void gcfNucleus::setCustomValues(double newSigma, double newEstar, double newCpp0 ,double newCnn0 ,double newCpn0, double newCpn1){

  sigmaCM = newSigma;
  Estar = newEstar;
  mySRC->set_Cpp0(newCpp0);
  mySRC->set_Cnn0(newCnn0);
  mySRC->set_Cpn0(newCpn0);
  mySRC->set_Cpn1(newCpn1);

}

void gcfNucleus::set_sigmaCM(double newSigma){

  sigmaCM = newSigma;
  
}

void gcfNucleus::set_Estar(double newEstar){

  Estar = newEstar;

}

void gcfNucleus::set_sigmaE(double newSigE){

  sigmaE = newSigE;
  random_Estar = true;
  
}

int gcfNucleus::get_Z()
{
  return Z;
}

int gcfNucleus::get_N()
{
  return N;
}

int gcfNucleus::get_A()
{
  return A;
}

double gcfNucleus::get_mA()
{
  return mA;
}

double gcfNucleus::get_mbar()
{
  return mA/A;
}

double gcfNucleus::get_mAmpp()
{
  return mAmpp + Estar;
}

double gcfNucleus::get_mAmpn()
{
  return mAmpn + Estar;
}

double gcfNucleus::get_mAmnn()
{
  return mAmnn + Estar;
}

double gcfNucleus::get_mAmpp_random(double &Estar, TRandom3* myRand)
{
  Estar = get_Estar_random(myRand);
  return mAmpp + Estar;
}

double gcfNucleus::get_mAmpn_random(double &Estar, TRandom3* myRand)
{
  Estar = get_Estar_random(myRand);
  return mAmpn + Estar;
}

double gcfNucleus::get_mAmnn_random(double &Estar, TRandom3* myRand)
{
  Estar = get_Estar_random(myRand);
  return mAmnn + Estar;
}
double gcfNucleus::get_sigmaCM()
{
  return sigmaCM;
}

double gcfNucleus::get_Estar()
{
  return Estar;
}

double gcfNucleus::get_sigmaE()
{
  return sigmaE;
}

double gcfNucleus::get_Estar_random(TRandom3* myRand)
{
  while(true)
    {
      double Estar_rand = myRand->Gaus(Estar,sigmaE);
      if (Estar_rand > 0.)
	return Estar_rand;
    }
}

gcfSRC * gcfNucleus::get_SRC()
{
  return mySRC;
}

bool gcfNucleus::get_Estar_randomization()
{
  return random_Estar;
}

double gcfNucleus::get_S(double k_rel, int l_type, int r_type)
{
  return mySRC->get_S(k_rel, l_type, r_type);
}
