#include <iostream>
#include <unistd.h>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "QEGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

int nEvents;
TFile * outfile;
bool verbose = false;
TVector3 boost_vector;
gcfNucleus * myInfo;
TRandom3 * myRand;
eNCrossSection * myCS;
QEGenerator * myGen;
TTree * outtree;
bool doLC;
bool doCoul = false;
double deltaECoul = 0;

// Tree variables
Double_t pe[3], pLead[3], pRec[3], pAm2[3];
Double_t weight;
Int_t lead_type, rec_type;

void Usage()
{
  cerr << "Usage: ./genQE <Z> <N> <Beam energy (GeV)> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-u: Specify NN interaction (default AV18)\n"
       << "-s: Specify sigma_CM [GeV/c]\n"
       << "-E: Specify E* [GeV]\n"
       << "-M: Use randomized E* according to Barack's values\n"
       << "-O: Turn on peaking radiation\n"
       << "-C: Turn on coulomb correction\n"
       << "-l: Use Lightcone cross section\n"
       << "-h: Print this message and exit\n\n\n";
}

bool init(int argc, char ** argv)
{
  int numargs = 6;
 
  if (argc < numargs)
    {
      Usage();
      return false;
    }

  // Read in the arguments
  int Z = atoi(argv[1]);
  int N = atoi(argv[2]);
  double Ebeam = atof(argv[3]);
  outfile = new TFile(argv[4],"RECREATE");
  nEvents = atoi(argv[5]);

  csMethod csMeth=cc1;
  ffModel ffMod=kelly;
  
  // Optional flags
  bool custom_ps = false;
  char * phase_space;
  char * uType = "AV18";
  double sigCM = 0.;
  bool do_sigCM = false;
  double Estar = 0.;
  bool do_Estar = false;
  double sigmaE = 0.;
  bool do_sigmaE = false;
  bool doRad = false;
  doLC = false;
  
  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:u:s:E:MOClh")) != -1)
    switch(c)
      {
	
      case 'v':
	verbose = true;
	break;
      case 'P':
	custom_ps = true;
	phase_space = optarg;
	break;
      case 'u':
	uType = optarg;
	break;
      case 's':
	do_sigCM = true;
	sigCM = atof(optarg);
	break;
      case 'E':
	do_Estar = true;
	Estar = atof(optarg);
	break;
      case 'M':
        do_Estar = true;
	Estar = 0.01732;
        do_sigmaE = true;
	sigmaE = 0.009571;
	break;
      case 'O':
	doRad = true;
	break;
      case 'C':
  doCoul = true;
  break;
      case 'l':
	doLC = true;
	break;
      case 'h':
	Usage();
	return false;
      default:
	abort();
	
      }

  // Initialize objects
  myInfo = new gcfNucleus(Z,N,uType);
  myRand = new TRandom3(0);
  myCS = new eNCrossSection(csMeth,ffMod);

  if (do_sigCM)
    myInfo->set_sigmaCM(sigCM);
  if (do_Estar)
    myInfo->set_Estar(Estar);
  if (do_sigmaE)
    myInfo->set_sigmaE(sigmaE);

  // Coulomb correction only available on Carbon 12
  if (doCoul && Z == 6 && N == 6)
    deltaECoul = carbonDeltaECoulombGeV;
  
  // Initialize generator
  myGen = new QEGenerator(Ebeam + deltaECoul, myInfo, myCS, myRand);
  if (custom_ps)
    myGen->parse_phase_space_file(phase_space);
  if (doRad)
    myGen->set_doRad(true);
  
  // Set up the tree
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("pAm2",pAm2,"pAm2[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  
  return true;
  
}

void coulombCorrection(TLorentzVector &p, double deltaE) 
{
  // If p is 0, don't correct
  if(p.E() == 0) {
    return;
  }

  TVector3 p3 = p.Vect();
  double pMag = p3.Mag();

  if(deltaE < 0 && pMag < abs(deltaE)) {
    p.SetPxPyPzE(0, 0, 0, 0);
    return;
  }

  double deltaP = - pMag + sqrt(pow(pMag, 2) + 2 * p.E() * deltaE + pow(deltaE, 2));
  p3.SetMag(pMag + deltaP);

  p.SetPxPyPzE(p3.X(), p3.Y(), p3.Z(), p.E() + deltaE);
}

void evnt(int event)
{

  TLorentzVector vk;
  TLorentzVector vLead;
  TLorentzVector vRec;
  TLorentzVector vAm2;

  if (doLC)
    myGen->generate_event_lightcone(weight, lead_type, rec_type, vk, vLead, vRec, vAm2);
  else
    myGen->generate_event(weight, lead_type, rec_type, vk, vLead, vRec, vAm2);

  // cout << "(" << vk.X() << ", " << vk.Y() << ", " << vk.Z() << ", " << vk.E() << ")" << endl;
  // cout << "(" << vLead.X() << ", " << vLead.Y() << ", " << vLead.Z() << ", " << vLead.E() << ")" << endl;

  // if (doCoul) {
  //   vk = coulombCorrection(vk, -deltaECoul);

  //   if (lead_type == pCode) {
  //     vLead = coulombCorrection(vLead, deltaECoul); 
  //   }
  //   if (rec_type == pCode) {
  //     vRec = coulombCorrection(vRec, deltaECoul);
  //   }
  // }

  // cout << "(" << vk.X() << ", " << vk.Y() << ", " << vk.Z() << ", " << vk.E() << ")" << endl;
  // cout << "(" << vLead.X() << ", " << vLead.Y() << ", " << vLead.Z() << ", " << vLead.E() << ")" << endl;

  // string testing;
  // cin >> testing;
  // cout << endl;
  
  pe[0] = vk.X();
  pe[1] = vk.Y();
  pe[2] = vk.Z();
  pLead[0] = vLead.X();
  pLead[1] = vLead.Y();
  pLead[2] = vLead.Z();
  pRec[0] = vRec.X();
  pRec[1] = vRec.Y();
  pRec[2] = vRec.Z();
  pAm2[0] = vAm2.X();
  pAm2[1] = vAm2.Y();
  pAm2[2] = vAm2.Z();

  if (weight > 0.)
    outtree->Fill();
  
}

void fini()
{
  outtree->SetName("genT");
  outtree->Write();
  outfile->Delete("genTbuffer;*");
  outfile->Close();
}

int main(int argc, char ** argv)
{

  if (not init(argc,argv))
    return -1;

  for (int event=0; event < nEvents; event++)
    {
      if ((event %100000==0) && (verbose))
	cout << "Working on event " << event << "\n";

      evnt(event);
    }

  fini();
  
  return 0;
  
}
