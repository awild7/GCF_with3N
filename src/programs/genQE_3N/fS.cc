#include <fstream>
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

const int theta_bins = 61;
const double theta_min = 0;
const double theta_max = 180;
const double theta_width = (theta_max - theta_min)/((double)theta_bins-1);

const int k_cm_bins = 51;
const double k_cm_min = 0;
const double k_cm_max = 5;
const double k_cm_width = (k_cm_max - k_cm_min)/((double)k_cm_bins-1);

const int k_rel_bins = 51;
const double k_rel_min = 0;
const double k_rel_max = 5;
const double k_rel_width = (k_rel_max - k_rel_min)/((double)k_rel_bins-1);

double interpolate(double matrix[k_rel_bins], double k_rel);
double interpolate(double matrix[k_cm_bins][k_rel_bins], double k_cm, double k_rel);
double interpolate(double matrix[theta_bins][k_cm_bins][k_rel_bins], double theta, double k_cm, double k_rel);

void Usage()
{
  cerr << "Usage: ./code <output file>\n\n";
}

int main(int argc, char ** argv)
{

  if(argc != 1)
    {
      std::cerr<<"Wrong number of arguments.\n";
      Usage();
      return -1;
    }

  ifstream densityFile("3N_Distribution.dnst");
  if(!densityFile.is_open()){
    cout<<"3 nucleon distribution file failed to load.\n"
	<<"Aborting...\n\n\n";
    exit(-2);
  }
  double density_matrix_pp[theta_bins][k_cm_bins][k_rel_bins];
  double density_matrix_pn[theta_bins][k_cm_bins][k_rel_bins];


  densityFile.ignore(1000,'\n');  
  for(int i = 0; i < theta_bins; i++){
    for(int j = 0; j < k_cm_bins; j++){
      for(int k = 0; k < k_rel_bins; k++){
	std::string theta;
	std::string k_cm;
	std::string k_rel;
	std::string C_tot;
	std::string C_pp;
	std::string C_nn;
	std::string C_pn;
	
	densityFile >> theta >> k_cm >> k_rel >> C_tot >> C_pp >> C_nn >> C_pn;
	density_matrix_pp[i][j][k] = std::stod(C_pp);
	density_matrix_pn[i][j][k] = std::stod(C_pn);
	densityFile.ignore(1000,'\n');
      }
    }
  }

  
  double p_all_min = 0.35;
  double S_min = 0.0;
  double S_max = 15.0;
  const int array_bins = 500;
  cout<<"S: "<<S_min<<","<<S_max<<","<<array_bins<<endl;
       
  double f1_arr[4]={0.3333,0.4998,0.4998,0.4998};
  double f2_arr[4]={0.3333,0.0004,0.2501,0.4998};
  for(int k = 0; k<4; k++){
    double f1 = f1_arr[k];
    double f2 = f2_arr[k];
    char temp_name[100];

    sprintf(temp_name,"1D_new/Shape_3N_f1_%i_f2_%i.txt",(int)round(f1*100),(int)round(f2*100));
    ofstream ppFile(temp_name);
    
    for(int i = 0; i < array_bins+1; i++){
      double S = (((double)i/array_bins)*(S_max-S_min)) + S_min;
      double k_tot = S;
      double f3 = 1 - f1 - f2;
      double k_1 = f1 * k_tot;
      double k_2 = f2 * k_tot;
      double k_3 = f3 * k_tot;
      double theta_12 = M_PI - acos((sq(k_1) + sq(k_2) - sq(k_3))/(2*k_1*k_2));

      TVector3 v_k_1(0.0,0.0,k_1);
      TVector3 v_k_2;
      v_k_2.SetMagThetaPhi(k_2,theta_12,0.0);
      TVector3 v_k_3 = - v_k_1 - v_k_2;
      if(fabs(v_k_3.Mag()-k_3)>0.000001){
	cout<<"Geometry Problem!"<<endl;
      }
      double sin_theta_12 = sin(theta_12);

      TVector3 v_k_rel = (v_k_1 - v_k_2)*0.5;
      double k_rel = v_k_rel.Mag();
      TVector3 v_k_cm = -v_k_3;
      double k_cm = v_k_cm.Mag();
      double sin_theta_prime = sin(v_k_cm.Angle(v_k_rel));
      //Units of Degrees
      double theta_prime = asin(sin_theta_prime)*180/M_PI;
      
      double p_1 = k_1 * GeVfm;
      double p_2 = k_2 * GeVfm;
      double p_3 = k_3 * GeVfm;
      //Define the Jacobian for this set of variables
      double J = (sq(k_rel)*sq(k_cm)*sin_theta_prime)/(sq(k_1)*sq(k_2)*sin_theta_12);
      //Now swtich variables again to S,f1,f2
      J /= k_1*k_2*sin_theta_12 / (k_tot*k_tot*k_tot);
      J=1;
      //if((k_cm>5) || (k_rel>5) || (p_1>1.0) || (p_2>1.0) || (p_3>1.0) || (p_1<p_all_min) || (p_2<p_all_min) || (p_3<p_all_min)){
      if((k_cm>5) || (k_rel>5) || (p_1>1.0) || (p_2>1.0) || (p_3>1.0)){
	ppFile << 0 << "\n";
      }
      else{	  	  
	ppFile << J * interpolate(density_matrix_pp,theta_prime,k_cm,k_rel) << "\n";	
      }      
    }
    ppFile.close();
  }  
  
  return 0;  
}
  
double interpolate(double matrix[k_rel_bins], double k_rel){
  //Get the "bin" for the value
  double float_bin = (k_rel - k_rel_min)/k_rel_width;
  int lower_bin = float_bin;
  int upper_bin = lower_bin+1;
  double delta = float_bin - (double)lower_bin;

  //Check to see if bins are out of bounds.
  //If so, return boundary
  if(lower_bin<0){
    return matrix[0];
  }
  if(upper_bin>k_rel_bins-1){
    return matrix[k_rel_bins-1];
  }
  //Otherwise return interpolation
  return delta*matrix[upper_bin] + (1-delta)*matrix[lower_bin];
}

double interpolate(double matrix[k_cm_bins][k_rel_bins], double k_cm, double k_rel){
  //Get the "bin" for the value
  double float_bin = (k_cm - k_cm_min)/k_cm_width;
  int lower_bin = float_bin;
  int upper_bin = lower_bin+1;
  double delta = float_bin - (double)lower_bin;

  //Check to see if bins are out of bounds.
  //If so, return boundary
  if(lower_bin<0){
    return interpolate(matrix[0],k_rel);
  }
  if(upper_bin>k_cm_bins-1){
    return interpolate(matrix[k_cm_bins-1],k_rel);
  }
  //Otherwise return interpolation
  return delta*interpolate(matrix[upper_bin],k_rel) + (1-delta)*interpolate(matrix[lower_bin],k_rel);
}

double interpolate(double matrix[theta_bins][k_cm_bins][k_rel_bins], double theta, double k_cm, double k_rel){
  //Get the "bin" for the value
  double float_bin = (theta - theta_min)/theta_width;
  int lower_bin = float_bin;
  int upper_bin = lower_bin+1;
  double delta = float_bin - (double)lower_bin;

  //Check to see if bins are out of bounds.
  //If so, return boundary
  if(lower_bin<0){
    return interpolate(matrix[0],k_cm,k_rel);
  }
  if(upper_bin>theta_bins-1){
    return interpolate(matrix[theta_bins-1],k_cm,k_rel);
  }  
  //Otherwise return interpolation
  return delta*interpolate(matrix[upper_bin],k_cm,k_rel) + (1-delta)*interpolate(matrix[lower_bin],k_cm,k_rel);
}
