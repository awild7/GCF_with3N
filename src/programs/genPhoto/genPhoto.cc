#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "constants.hh"
#include "helpers.hh"
#include "cross_sections/photoCrossSection.hh"
#include "nucleus/gcfNucleus.hh"
#include "photoGenerator.hh"

using namespace std;

int nEvents;
TFile * outfile;
bool verbose = false;
TVector3 boost_vector;
gcfNucleus * myInfo;
TRandom3 * myRand;
photoCrossSection * myCS;
photoGenerator * myGen;
TTree * outtree;

// Tree variables
Double_t pMeson[3], pBaryon[3], pRec[3], pAm2[3];
Double_t weight;
Int_t meson_type, baryon_type, rec_type;

void Usage()
{
  cerr << "Usage: ./genPhoto <Z> <N> <Beam energy (GeV)> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-T: Boost to inverse kinematics. In this case, the given beam energy is the kinetic energy per nucleon"
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

  // Optional flags
  bool custom_ps = false;
  char * phase_space;
  
  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:h")) != -1)
    switch(c)
      {
	
      case 'v':
	verbose = true;
	break;
      case 'P':
	custom_ps = true;
	phase_space = optarg;
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
  myCS = new photoCrossSection();
  
  // Initialize generator
  myGen = new photoGenerator(Ebeam, myInfo, myCS, myRand);
  if (custom_ps)
    myGen->parse_phase_space_file(phase_space);

  // Set up the tree
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("meson_type",&meson_type,"meson_type/I");
  outtree->Branch("baryon_type",&baryon_type,"baryon_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pMeson",pMeson,"pMeson[3]/D");
  outtree->Branch("pBaryon",pBaryon,"pBaryon[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("pAm2",pAm2,"pAm2[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  
  return true;
  
}

void evnt(int event)
{
  
  TLorentzVector vMeson;
  TLorentzVector vBaryon;
  TLorentzVector vRec;
  TLorentzVector vAm2;

  myGen->generate_event(weight, meson_type, baryon_type, rec_type, vMeson, vBaryon, vRec, vAm2);

  pMeson[0] = vMeson.X();
  pMeson[1] = vMeson.Y();
  pMeson[2] = vMeson.Z();
  pBaryon[0] = vBaryon.X();
  pBaryon[1] = vBaryon.Y();
  pBaryon[2] = vBaryon.Z();
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
