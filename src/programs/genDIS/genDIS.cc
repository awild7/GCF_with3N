#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "constants.hh"
#include "helpers.hh"
#include "cross_sections/DISCrossSection.hh"
#include "nucleus/gcfNucleus.hh"
#include "DISGenerator.hh"

using namespace std;

int nEvents;
TFile * outfile;
bool verbose = false;
TVector3 boost_vector;
gcfNucleus * myInfo;
TRandom3 * myRand;
DISCrossSection * myCS;
DISGenerator * myGen;
TTree * outtree;

// Tree variables
Double_t pe[3], q[3], p1_onshell[3], pHadron[3], pRec[3], pAm2[3];
Double_t nu, E1_onshell, EHadron, ERec, EAm2, weight;
Int_t lead_type, rec_type, ipart;

void Usage()
{
  cerr << "Usage: ./genDIS <Z> <N> <Beam energy (GeV)> <path/to/CTEQ.pds> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-h: Print this message and exit\n\n\n";
}

bool init(int argc, char ** argv)
{
  int numargs = 7;
 
  if (argc < numargs)
    {
      Usage();
      return false;
    }

  // Read in the arguments
  int Z = atoi(argv[1]);
  int N = atoi(argv[2]);
  double Ebeam = atof(argv[3]);
  char * pdsFile = argv[4];
  outfile = new TFile(argv[5],"RECREATE");
  nEvents = atoi(argv[6]);

  // Optional flags
  bool custom_ps = false;
  char * phase_space;
  double Estar = 0.;
  bool do_Estar = false;
  double sigmaE = 0.;
  bool do_sigmaE = false;
  
  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:E:Mh")) != -1)
    switch(c)
      {
	
      case 'v':
	verbose = true;
	break;
      case 'P':
	custom_ps = true;
	phase_space = optarg;
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
      case 'h':
	Usage();
	return false;
      default:
	abort();
	
      }

  // Initialize objects
  myInfo = new gcfNucleus(Z,N,AV18);
  myRand = new TRandom3(0);
  myCS = new DISCrossSection(pdsFile);

  if (do_Estar)
    myInfo->set_Estar(Estar);
  if (do_sigmaE)
    myInfo->set_sigmaE(sigmaE);
  
  // Initialize generator
  myGen = new DISGenerator(Ebeam, myInfo, myCS, myRand);
  if (custom_ps)
    myGen->parse_phase_space_file(phase_space);
  
  // Set up the tree
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("ipart",&ipart,"ipart/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("q",q,"q[3]/D");
  outtree->Branch("nu",&nu,"nu/D");
  outtree->Branch("p1_onshell",p1_onshell,"p1_onshell[3]/D");
  outtree->Branch("E1_onshell",&E1_onshell,"E1_onshell/D");
  outtree->Branch("pHadron",pHadron,"pHadron[3]/D");
  outtree->Branch("EHadron",&EHadron,"EHadron/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("ERec",&ERec,"ERec/D");
  outtree->Branch("pAm2",pAm2,"pAm2[3]/D");
  outtree->Branch("EAm2",&EAm2,"EAm2/D");
  outtree->Branch("weight",&weight,"weight/D");
  
  return true;
  
}

void evnt(int event)
{
  
  TLorentzVector vk;
  TLorentzVector vq;
  TLorentzVector v1_onshell;
  TLorentzVector vHadron;
  TLorentzVector vRec;
  TLorentzVector vAm2;

  myGen->generate_event(weight, lead_type, rec_type, ipart, vk, vq, v1_onshell, vHadron, vRec, vAm2);

  pe[0] = vk.X();
  pe[1] = vk.Y();
  pe[2] = vk.Z();
  q[0] = vq.X();
  q[1] = vq.Y();
  q[2] = vq.Z();
  nu = vq.T();
  p1_onshell[0] = v1_onshell.X();
  p1_onshell[1] = v1_onshell.Y();
  p1_onshell[2] = v1_onshell.Z();
  E1_onshell = v1_onshell.T();
  pHadron[0] = vHadron.X();
  pHadron[1] = vHadron.Y();
  pHadron[2] = vHadron.Z();
  EHadron = vHadron.T();
  pRec[0] = vRec.X();
  pRec[1] = vRec.Y();
  pRec[2] = vRec.Z();
  ERec = vRec.T();
  pAm2[0] = vAm2.X();
  pAm2[1] = vAm2.Y();
  pAm2[2] = vAm2.Z();
  EAm2 = vAm2.T();
  
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
