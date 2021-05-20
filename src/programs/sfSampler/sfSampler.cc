#include <iostream>
#include <unistd.h>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "generator/gcfGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

int nEvents;
TFile * outfile;
bool verbose = false;
TVector3 boost_vector;
gcfNucleus * myInfo;
TRandom3 * myRand;
gcfGenerator * myGen;
TTree * outtree;

// Tree variables
Double_t pI[3], pRec[3];
Double_t weight;
Int_t lead_type, rec_type;

void Usage()
{
  cerr << "Usage: ./genQE <Z> <N> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-u: Specify NN interaction (default AV18)\n"
       << "-s: Specify sigma_CM [GeV/c]\n"
       << "-E: Specify E* [GeV]\n"
       << "-k: Specify pRel hard cutoff [GeV/c]\n"
       << "-r: Randomize nuclear properties\n"
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
  int Z = atoi(argv[1]);
  int N = atoi(argv[2]);
  outfile = new TFile(argv[3],"RECREATE");
  nEvents = atoi(argv[4]);
  
  // Optional flags
  bool custom_ps = false;
  char * phase_space;
  char * uType = (char *)"AV18";
  double sigCM = 0.;
  bool do_sigCM = false;
  double Estar = 0.;
  bool do_Estar = false;
  double sigmaE = 0.;
  bool do_sigmaE = false;
  double kCut = 0.25;
  bool do_kCut = false;
  bool rand_flag = false;
  if ((Z == 1) and (N == 1))
    uType = (char *)"AV18_deut";
  
  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:u:s:E:k:rh")) != -1)
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
      case 'k':
	do_kCut = true;
	kCut = atof(optarg);
	break;
      case 'r':
	rand_flag = true;
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

  if (rand_flag)
    myInfo->randomize(myRand);
  if (do_sigCM)
    myInfo->set_sigmaCM(sigCM);
  if (do_Estar)
    myInfo->set_Estar(Estar);
  if (do_sigmaE)
    myInfo->set_sigmaE(sigmaE);

  // Initialize generator
  myGen = new gcfGenerator(myInfo, myRand);
  if ((Z == 1) and (N == 1))
    myGen->set_deuteron();
  if (custom_ps)
    myGen->parse_phase_space_file(phase_space);
  if (do_kCut)
    myGen->set_pRel_cut(kCut);

  // Set up the tree
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pI",pI,"pI[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  
  return true;
  
}

void evnt(int event)
{

  weight = 1;

  // Decide what kind of proton or neutron pair we are dealing with
  lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  // Call decay function
  TVector3 vi, vr;  
  myGen->decay_function(weight,lead_type,rec_type,vi,vr);

  pI[0] = vi.X();
  pI[1] = vi.Y();
  pI[2] = vi.Z();
  pRec[0] = vr.X();
  pRec[1] = vr.Y();
  pRec[2] = vr.Z();
  
  if ((weight > 0.) or true)
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
