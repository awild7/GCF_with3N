#include <iostream>
#include <unistd.h>
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

double E;

// Tree variables
Double_t pei[3], pef[3], pLead[3], pRec[3], pAm2[3], q[3], pMiss[3], pCM[3], pRel[3];
Double_t QSq, xB, nu, pef_Mag, q_Mag, pLead_Mag, pRec_Mag, pMiss_Mag, pCM_Mag, pRel_Mag, theta_pmq, theta_prq, mass_Aminus2, estar, mAm2, weight, lcweight;
Int_t lead_type, rec_type;

void Usage()
{
  cerr << "Usage: ./genQE <Z> <N> <Beam energy (GeV)> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-M: Use randomized E* according to Barack's values\n"
       << "-O: Turn on peaking radiation\n"
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
  E = Ebeam;
  outfile = new TFile(argv[4],"RECREATE");
  nEvents = atoi(argv[5]);

  csMethod csMeth=cc1;
  ffModel ffMod=kelly;
  
  // Optional flags
  bool custom_ps = false;
  char * phase_space;
  double Estar = 0.;
  bool do_Estar = false;
  double sigmaE = 0.;
  bool do_sigmaE = false;
  bool doRad = false;
  
  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:MOh")) != -1)
    switch(c)
      {
	
      case 'v':
	verbose = true;
	break;
      case 'P':
	custom_ps = true;
	phase_space = optarg;
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
      case 'h':
	Usage();
	return false;
      default:
	abort();
	
      }

  // Initialize objects
  myInfo = new gcfNucleus(Z,N,AV18);
  myRand = new TRandom3(0);
  myCS = new eNCrossSection(csMeth,ffMod);
  
  if (do_Estar)
    myInfo->set_Estar(Estar);
  if (do_sigmaE)
    myInfo->set_sigmaE(sigmaE);
  
  // Initialize generator
  myGen = new QEGenerator(Ebeam, myInfo, myCS, myRand);
  if (custom_ps)
    myGen->parse_phase_space_file(phase_space);
  if (doRad)
    myGen->set_doRad(true);

  // Set up the tree
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pei",pei,"pei[3]/D");
  outtree->Branch("pef",pef,"pef[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("pAm2",pAm2,"pAm2[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  outtree->Branch("lcweight",&lcweight,"lcweight/D");
  outtree->Branch("QSq",&QSq,"QSq/D");
  outtree->Branch("xB",&xB,"xB/D");
  outtree->Branch("q",q,"q[3]/D");
  outtree->Branch("pMiss",pMiss,"pMiss[3]/D");
  outtree->Branch("pCM",pCM,"pCM[3]/D");
  outtree->Branch("pRel",pRel,"pRel[3]/D");
  outtree->Branch("nu",&nu,"nu/D");
  outtree->Branch("pef_Mag",&pef_Mag,"pef_Mag/D");
  outtree->Branch("q_Mag",&q_Mag,"q_Mag/D");
  outtree->Branch("pLead_Mag",&pLead_Mag,"pLead_Mag/D");
  outtree->Branch("pRec_Mag",&pRec_Mag,"pRec_Mag/D");
  outtree->Branch("pMiss_Mag",&pMiss_Mag,"pMiss_Mag/D");
  outtree->Branch("pCM_Mag",&pCM_Mag,"pCM_Mag/D");
  outtree->Branch("pRel_Mag",&pRel_Mag,"pRel_Mag/D");
  outtree->Branch("theta_pmq",&theta_pmq,"theta_pmq/D");
  outtree->Branch("theta_prq",&theta_prq,"theta_prq/D");
  outtree->Branch("Estar",&estar,"Estar/D");
  outtree->Branch("mAm2",&mAm2,"mAm2/D");
  
  return true;
  
}

void evnt(int event)
{

  TLorentzVector vbeam(0.,0.,E,E);
  
  TLorentzVector vk;
  TLorentzVector vLead;
  TLorentzVector vRec;
  TLorentzVector vAm2;

  myGen->generate_event(weight, lead_type, rec_type, vk, vLead, vRec, vAm2, estar);

  TLorentzVector vq = vbeam - vk;
  TLorentzVector vMiss = vLead - vq;
  TLorentzVector vCM = vMiss + vRec;
  TLorentzVector vRel = 0.5*(vMiss - vRec);

  pei[0] = vbeam.X();
  pei[1] = vbeam.Y();
  pei[2] = vbeam.Z();
  pef[0] = vk.X();
  pef[1] = vk.Y();
  pef[2] = vk.Z();
  pLead[0] = vLead.X();
  pLead[1] = vLead.Y();
  pLead[2] = vLead.Z();
  pRec[0] = vRec.X();
  pRec[1] = vRec.Y();
  pRec[2] = vRec.Z();
  pAm2[0] = vAm2.X();
  pAm2[1] = vAm2.Y();
  pAm2[2] = vAm2.Z();
  q[0] = vq.X();
  q[1] = vq.Y();
  q[2] = vq.Z();
  pMiss[0] = vMiss.X();
  pMiss[1] = vMiss.Y();
  pMiss[2] = vMiss.Z();
  pCM[0] = vCM.X();
  pCM[1] = vCM.Y();
  pCM[2] = vCM.Z();
  pRel[0] = vRel.X();
  pRel[1] = vRel.Y();
  pRel[2] = vRel.Z();
  nu = vq.T();
  QSq = -vq.Mag2();
  xB = QSq/(2.*mN*nu);
  pef_Mag = vk.Vect().Mag();
  q_Mag = vq.Vect().Mag();
  pLead_Mag = vLead.Vect().Mag();
  pRec_Mag = vRec.Vect().Mag();
  pMiss_Mag = vMiss.Vect().Mag();
  pCM_Mag = vCM.Vect().Mag();
  pRel_Mag = vRel.Vect().Mag();
  theta_pmq = vq.Vect().Angle(vMiss.Vect());
  theta_prq = vq.Vect().Angle(vRec.Vect());
  lcweight = weight;
  mAm2 = vAm2.Mag();
  
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
