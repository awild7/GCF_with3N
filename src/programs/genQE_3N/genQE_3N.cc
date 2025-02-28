#include <iostream>
#include <unistd.h>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "QEGenerator_3N.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

int nEvents;
TFile * outfile;
bool verbose = false;
gcfNucleus * myInfo;
TRandom3 * myRand;
eNCrossSection * myCS;
QEGenerator_3N * myGen;
TTree * outtree;

// Tree variables
Double_t pe[3], pLead[3], p2[3], p3[3], pAm3[3];
Double_t weight;
Int_t N1_type, N2_type, N3_type;

void Usage()
{
  cerr << "Usage: ./code <Beam energy (GeV)> [Interaction #] <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-c: Specify eN cross section model (default cc1)\n"
       << "-h: Print this message and exit\n\n\n";
}

bool init(int argc, char ** argv)
{
  int numargs = 5;
 
  if (argc < numargs)
    {
      Usage();
      return false;
    }

  // Read in the arguments
  double Ebeam = atof(argv[1]);
  double u = atoi(argv[2]);
  outfile = new TFile(argv[3],"RECREATE");
  nEvents = atoi(argv[4]);

  ffModel ffMod=kelly;
  csMethod csMeth=cc1;

  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vc:h")) != -1)
    switch(c)
      {
	
      case 'v':
	verbose = true;
	break;
      case 'c':
	if (strcmp(optarg,"onshell")==0)
	  csMeth=onshell;
	else if ((strcmp(optarg,"cc1")==0) or (atoi(optarg)==1))
	  csMeth=cc1;
	else if ((strcmp(optarg,"cc2")==0) or (atoi(optarg)==2))
	  csMeth=cc2;
	else
	  {
	    cerr << "Invalid cross section designation. Allowed values are onshell, cc1 and cc2. Aborting...\n";
	    return -1;
	  }
	break;
      case 'h':
	Usage();
	return false;
      default:
	abort();	
      }

  // Initialize objects
  myRand = new TRandom3(0);
  myCS = new eNCrossSection(csMeth,ffMod);
  
  // Initialize generator
  myGen = new QEGenerator_3N(Ebeam, myCS, u, myRand);
  myGen->set_theta_k_maxmin(5,40);
  // Set up the tree
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("N1_type",&N1_type,"N1_type/I");
  outtree->Branch("N2_type",&N2_type,"N2_type/I");
  outtree->Branch("N3_type",&N3_type,"N3_type/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("p2",p2,"p2[3]/D");
  outtree->Branch("p3",p3,"p3[3]/D");
  outtree->Branch("pAm3",pAm3,"pAm3[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  
  return true;
  
}

void evnt(int event)
{

  TLorentzVector vk;
  TLorentzVector vLead;
  TLorentzVector v2;
  TLorentzVector v3;
  TLorentzVector vAm3;

  myGen->generate_event(weight, N1_type, N2_type, N3_type, vk, vLead, v2, v3, vAm3);
  
  pe[0] = vk.X();
  pe[1] = vk.Y();
  pe[2] = vk.Z();
  pLead[0] = vLead.X();
  pLead[1] = vLead.Y();
  pLead[2] = vLead.Z();
  p2[0] = v2.X();
  p2[1] = v2.Y();
  p2[2] = v2.Z();
  p3[0] = v3.X();
  p3[1] = v3.Y();
  p3[2] = v3.Z();
  pAm3[0] = vAm3.X();
  pAm3[1] = vAm3.Y();
  pAm3[2] = vAm3.Z();

  if(weight>0.){
    outtree->Fill();
  }
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
