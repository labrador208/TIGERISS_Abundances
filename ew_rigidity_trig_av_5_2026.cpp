#include <TROOT.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <iostream>
#include "Riostream.h"
#include <string.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TApplication.h>
#include <math.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <string.h>
#include "stdio.h"
#include <stdlib.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TVirtualFitter.h>
#include <TSpline.h>
#include <sstream>
//#include <TPaletteAxis.h>
#include <TGraphAsymmErrors.h>
#include <TColor.h>


#define ISS_inc 51.6445 // Inclination of ISS in degrees
#define amu_MeV 931.49432 // unified mass unit in MeV/c^2
#define c_vac 2.99792458e8 // speed of light in vacuum in m/s (exact)
#define E_CALET_hi 10e8 // 100 TeV maximum energy for CALET

#define long_bins 25 // number of longitude bins in geomagnetic model: 15 degree increments from 0 to 360 (inclusive)
#define long_bins_calc 24 // limit for actual calculation, as 0 360 degrees are the same
#define lat_bins 37  // number of latitude bins in geomagnetic model: 5 degree increments from 90 to -90 degrees

#define long_bins_1 360 // number of longitude bins for 1 degree pitch
#define lat_bins_1 180 // number of latitude bins for 1 degree pitch

#define geo_factor_m 0.4397  // 0.12   // Geometric Factor CALET m^2-sr
#define geo_factor_cm 4.397e3 // 1.2e3 // Geometric Factor of CALET cm^2-sr
#define geo_factor_m_s  0.12   // Full Instrument Geometric Factor CALET m^2-sr
#define geo_factor_cm_s  1.2e3 // Full Instrument Geometric Factor of CALET cm^2-sr



// Interaction correction parameters
#define S_fac 1.277 // scaling factor: Rn = SRe
#define NA 6.02214179e23 // Avogadro's Number (particles/mole)

//#define E_gf_thresh_Fe 600.0 // (MeV/nuc) Energy threshold between high and low acceptance geometry factor regimes for Fe

//#define E_gf_thresh_Fe 350.0 // (MeV/nuc) Energy threshold between high and low acceptance geometry factor regimes for Fe


FILE *outfile;

// KE/A = - amu_MeV * c_vac^2 + sqrt(amu_MeV^2 * c_vac^4 + (R^2 * Z^2 *c_vac^2)/A^2)

// R = sqrt(KE^2 + 2*KE*A*m*c^2)/(z*c)
/*
 To do:
 Generate plots of the geomagnetic cutoff potential and elemental cutoff energies.
 Generate histograms of geomagnetic cutoff energies for elements.
 Plot the latitude as a function of time.
 */

Double_t Integralfunc(TGraph *gr, Double_t a, Double_t b, Option_t *option, Double_t epsilon);

TGraph2D* Thresholds();

using namespace std;


int main(int argc, char* argv[])
{
    
    // ======================================================================
    // NICE ROOT STYLE: Put this at the top of your root script
    // or program to produce nice publication-quality output.  You
    // can also add it to your root.login script so all root plots
    // will use it automatically.
    // Settings courtesy of Karl Kosack for quality plot generation.
    int nicefont = 132;
    gROOT->SetStyle("Plain");               // Start with the plain style
    gStyle->SetLabelFont(nicefont,"x,y,z"); // nicer font for axes titles
    gStyle->SetTitleFont(nicefont,"");      // nicer font for pad title
    gStyle->SetTitleFont(nicefont,"");      // nicer font for pad title
    gStyle->SetTitleSize(0.06,"");          // bigger title font
    gStyle->SetTitleFont(nicefont,"x,y,z"); // nicer font for axes title
    gStyle->SetStatBorderSize(0);           // remove ugly border from stats
    gStyle->SetTitleBorderSize(0);          // remove ugly border from title
    gStyle->SetPadTickX(1);                 // Tics on opposite X axis
    gStyle->SetPadTickY(1);                 // Tics on opposite Y axis
    gStyle->SetTitleX(0.15f);                // Center the title
    gStyle->SetTitleW(0.7f);                // Center the title
    gStyle->SetOptStat(0);                  // turn off stats (1 to enable)
    gStyle->SetOptFit(0);                   // turn off fit (1 to enable)
    gStyle->SetOptTitle(1);                 // 0 to disable title
    // gStyle->SetOptStat(110010);
    
    
    
    
    
//  TApplication theApp("App", &argc, argv);
    
    ifstream in_1, in_2;
    ofstream out_1;
  
    
    int Z_iter=92;
    double E_gf_thresh_Fe = 350.0;
    char name_ext1[20] = "_tigeriss1";
    char name_ext2[20] = "_tigeriss2";
    
    if(argc>=2) Z_iter=atoi(argv[1]);
    
    if(argc>=4) {
        strcpy(name_ext1, argv[2]);
        strcpy(name_ext2, argv[3]);
    }
    
    // epsilon variable required for the Integralfunc
	double epsilon=1e-6;
    
    int i=0; // index
    
    // Geometry factors
    
    // source directory
    char source_dir[200] = "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/geofac_dir/";
    


    
    // East-West geometry factor
    //char dgf_ew[200] = "/data/calet/abundance_estimate/rigidity/code/geofac_dir/dgf_g_bin1_UH_geom.txt";
    char dgf_ew[200] = "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/TIGERISS_Abundances/angle/dgf_1d_bin1_";

    if(argc>=5) {
        strcpy(dgf_ew, source_dir);
        strcat(dgf_ew, "dgf_g_bin1");
        strcat(dgf_ew, name_ext1);
        strcat(dgf_ew, ".txt");
    }
    fprintf(stderr, "dgf_ew = %s\n",dgf_ew);
    
    int dgf_ew_b1_index=0;
	int dgf_ew_b1_count=1;
	
	double angle_dgf_ew_b1[181];
	double geo_factor_dgf_ew_b1[181];
	double geo_factor_dgf_ew_b1_adj[181];
    
    double blah[181];

  // Input files.
  char S_max_rel_in[200] = "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/TIGERISS_Abundances/S_max_integral_rel.txt";
  char S_min_rel_in[200] = "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/TIGERISS_Abundances/S_min_integral_rel.txt";

  char S_max_top[200] = "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/S_max_tigeriss_Z_theta_top_PS_Al_thick.txt";
  char S_min_top[200] = "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/S_min_tigeriss_Z_theta_top_PS_Al_thick.txt";

  int S_max_Z[100];
  double S_max_A[100];
  double S_max_R[100];
  double S_max_abund_in[100];
  double S_max_abund_out[100][90];
  int S_min_Z[100];
  double S_min_A[100];
  double S_min_R[100];
  double S_min_abund_in[100];
  double S_min_abund_out[100][90];

  int ang_ind=0;
  int index = 0;

  for(index=0; index<100; index++) {
    S_max_Z[index] = 0;
    S_max_A[index] = 0.0;
    S_max_R[index] = 0.0;
    S_max_abund_in[index] = 0.0;
    S_min_Z[index] = 0;
    S_min_A[index] = 0.0;
    S_min_R[index] = 0.0;
    S_min_abund_in[index] = 0.0;
    for(ang_ind=0; ang_ind<90; ang_ind++) {
      S_max_abund_out[index][ang_ind] = 0.0;  
      S_min_abund_out[index][ang_ind] = 0.0;
    }
  }

  int index_S_max_rel_in = 0;
  int index_S_min_rel_in = 0;
  int index_S_max_rel_in_fe = 0;
  int index_S_min_rel_in_fe = 0;


  index=0;
  in_1.open(S_max_rel_in);

  while(1) {
    in_1 >> S_max_Z[index] >> S_max_A[index] >> S_max_abund_in[index];
    if(in_1.eof() || !in_1.good()) break;
    S_max_R[index] = 1.58*pow(S_max_A[index],0.281);
    //fprintf(stderr,"Z: %2d A: %f     %E\n",S_max_Z[index],S_max_A[index],S_max_abund_in[index]);  

    if(S_max_Z[index] == Z_iter) {
      index_S_max_rel_in = index;
    }
    if(S_max_Z[index] == 26) {
      index_S_max_rel_in_fe = index;
    }

    index++;
  }
  in_1.close();
  in_1.clear();

  in_1.open(S_min_rel_in);
  index=0;

  while(1) {
    in_1 >> S_min_Z[index] >> S_min_A[index] >> S_min_abund_in[index];
    if(in_1.eof() || !in_1.good()) break;
    S_min_R[index] = 1.58*pow(S_min_A[index],0.281);

    if(S_min_Z[index] == Z_iter) {
      index_S_min_rel_in = index;
    }
    if(S_min_Z[index] == 26) {
      index_S_min_rel_in_fe = index;
    }

    index++;
  }
  in_1.close();
  in_1.clear();

  // reading in TIGERISS top interaction spectra modification factors

    char char_Z[2];
    char char_blah[10];

    int angle_deg[90];
    int Z_el[100];
    double S_max_delta_abund[100][90];
    double S_min_delta_abund[100][90];
    int el_ind=0;

    int index_S_max_top=0;
    int index_S_min_top=0;


    for(i=0; i < 100; i++) {
      Z_el[i] = 0;
      for(el_ind=0; el_ind<100; el_ind++) {
        S_max_delta_abund[el_ind][i] = 0.0;
        S_min_delta_abund[el_ind][i] = 0.0;
      }
    }



    el_ind=0;

    in_1.open(S_max_top);
    fprintf(stderr,"S_max top spectra modification file\n");

    
    in_1 >> char_blah >> char_Z;
    
    //fprintf(stderr,"%s %s %s\n",Z_char,A_char,X_char);
    
    for (i=0; i<90; i++) {
        in_1 >> angle_deg[i];
    }
    
    while (1) {
        in_1 >> Z_el[el_ind];
        //fprintf(stderr,"Z: %2d A: %.2f %s\n",Z_S_min,A_S_min,cs_S_min);
         for (i=0; i<90; i++) {
            in_1 >> S_max_delta_abund[el_ind][i];
        }
        if(in_1.eof() || !in_1.good()) break;

        if(Z_el[el_ind] == Z_iter) {
          index_S_max_top = el_ind;
        }


        el_ind++;
    }
    
    
    in_1.close();
    in_1.clear();

    el_ind=0;

    in_1.open(S_min_top);
    fprintf(stderr,"S_min top spectra modification file\n");

    
    in_1 >> char_blah >> char_Z;
    
    //fprintf(stderr,"%s %s %s\n",Z_char,A_char,X_char);
    
    for (i=0; i<90; i++) {
        in_1 >> angle_deg[i];
    }
    
    while (1) {
        in_1 >> Z_el[el_ind];
        fprintf(stderr,"Z_el: %2d \n",Z_el[el_ind]);
         for (i=0; i<90; i++) {
            in_1 >> S_min_delta_abund[el_ind][i];
            fprintf(stderr,"S_min_delta_abund[%d][%d]: %e\n",el_ind,i,S_min_delta_abund[el_ind][i]);
        }
        if(in_1.eof() || !in_1.good()) break;

        if(Z_el[el_ind] == Z_iter) {
          index_S_min_top = el_ind;
        }


        el_ind++;
    }
    
    
    in_1.close();
    in_1.clear();


    
    
    for (i=0; i<181; i++) {
        angle_dgf_ew_b1[i]=0.0;
        geo_factor_dgf_ew_b1[i]=0.0;
        geo_factor_dgf_ew_b1_adj[i]=0.0;
        blah[i]=0.0;
    }

	
    
	angle_dgf_ew_b1[0]=0.0;
	geo_factor_dgf_ew_b1[0]=0.0;
	
	in_1.open(dgf_ew);
	while(1) {
		in_1 >> angle_dgf_ew_b1[dgf_ew_b1_count] >> geo_factor_dgf_ew_b1[dgf_ew_b1_count];
        if(angle_dgf_ew_b1[dgf_ew_b1_count] < 45 || angle_dgf_ew_b1[dgf_ew_b1_count] > 135) {
            geo_factor_dgf_ew_b1_adj[dgf_ew_b1_count] = 0.75*geo_factor_dgf_ew_b1[dgf_ew_b1_count]/1E4;
        }
        else {
            geo_factor_dgf_ew_b1_adj[dgf_ew_b1_count] = geo_factor_dgf_ew_b1[dgf_ew_b1_count]/1E4;
        }
        blah[dgf_ew_b1_count]=geo_factor_dgf_ew_b1_adj[dgf_ew_b1_count];
        fprintf(stderr,"geo_factor_dgf_ew_b1_adj[%d]: %f\n",dgf_ew_b1_count,geo_factor_dgf_ew_b1_adj[dgf_ew_b1_count]);
        fprintf(stderr,"blah[%d]: %f\n",dgf_ew_b1_count,blah[dgf_ew_b1_count]);
        
		if(in_1.eof() || !in_1.good()) break;
		
		
		dgf_ew_b1_count++;
	}
	in_1.close();
	in_1.clear();
    
    
    // Getting the energy threshold function
    TGraph2D *gt = Thresholds();
    //get->Interpolate(Z,zenith);   // Command to get the threshold for Z and incidence angle

    
    
    // Limit Differential Geometry factor
    
    
    
    
    TGraph *gr_dgf_ew_b1 = new TGraph(dgf_ew_b1_count,angle_dgf_ew_b1,geo_factor_dgf_ew_b1);
	
	gr_dgf_ew_b1->SetMarkerStyle(1);
	gr_dgf_ew_b1->SetMarkerSize(1);
    //	gr_dgf_ew_b1->SetTitle("CALET Differential Geometry Factor");
	gr_dgf_ew_b1->SetTitle("");
	gr_dgf_ew_b1->GetXaxis()->SetTitle("East-West Incidence Angle (degrees)");
	gr_dgf_ew_b1->GetXaxis()->CenterTitle();
	gr_dgf_ew_b1->GetYaxis()->SetTitle("Geometry Factor (cm^{2}sr/(1 degree bins))");
	gr_dgf_ew_b1->GetYaxis()->CenterTitle();
	gr_dgf_ew_b1->GetXaxis()->SetTitleOffset(1.2);
	gr_dgf_ew_b1->GetYaxis()->SetTitleOffset(1.3);
	gr_dgf_ew_b1->GetYaxis()->SetLabelFont(nicefont);
	gr_dgf_ew_b1->GetYaxis()->SetTitleFont(nicefont);
	gr_dgf_ew_b1->GetXaxis()->SetLabelFont(nicefont);
	gr_dgf_ew_b1->GetXaxis()->SetTitleFont(nicefont);
	gr_dgf_ew_b1->GetXaxis()->SetRangeUser(0,180);
	gr_dgf_ew_b1->GetYaxis()->SetRangeUser(0,60);

    TGraph *gr_dgf_ew_b1_adj = new TGraph(dgf_ew_b1_count,angle_dgf_ew_b1,geo_factor_dgf_ew_b1_adj);
	
	gr_dgf_ew_b1_adj->SetMarkerStyle(1);
	gr_dgf_ew_b1_adj->SetMarkerSize(1);
    //	gr_dgf_ew_b1_adj->SetTitle("CALET Differential Geometry Factor");
	gr_dgf_ew_b1_adj->SetTitle("");
	gr_dgf_ew_b1_adj->GetXaxis()->SetTitle("East-West Incidence Angle (degrees)");
	gr_dgf_ew_b1_adj->GetXaxis()->CenterTitle();
	gr_dgf_ew_b1_adj->GetYaxis()->SetTitle("Geometry Factor (cm^{2}sr/(1 degree bins))");
	gr_dgf_ew_b1_adj->GetYaxis()->CenterTitle();
	gr_dgf_ew_b1_adj->GetXaxis()->SetTitleOffset(1.2);
	gr_dgf_ew_b1_adj->GetYaxis()->SetTitleOffset(1.3);
	gr_dgf_ew_b1_adj->GetYaxis()->SetLabelFont(nicefont);
	gr_dgf_ew_b1_adj->GetYaxis()->SetTitleFont(nicefont);
	gr_dgf_ew_b1_adj->GetXaxis()->SetLabelFont(nicefont);
	gr_dgf_ew_b1_adj->GetXaxis()->SetTitleFont(nicefont);
	gr_dgf_ew_b1_adj->GetXaxis()->SetRangeUser(0,180);
	gr_dgf_ew_b1_adj->GetYaxis()->SetRangeUser(0,60);
	
	
	TCanvas *C_dgf_ew_b1 = new TCanvas("C_dgf_ew_b1","CALET Differental Geometry Factor",0,0,800,800);
	C_dgf_ew_b1->SetLeftMargin(0.139594);
	C_dgf_ew_b1->SetRightMargin(0.0596447);
	C_dgf_ew_b1->SetBottomMargin(0.11911);
	C_dgf_ew_b1->SetTopMargin(0.0811518);
	
	gr_dgf_ew_b1->SetLineWidth(2);
	gr_dgf_ew_b1->Draw("AL");

    gr_dgf_ew_b1_adj->SetLineColor(2);
	gr_dgf_ew_b1_adj->SetLineWidth(2);
	gr_dgf_ew_b1_adj->Draw("L");
    
   /*
    // 2-D differential geometry factor for East-West and the transverse angle (delta)
    char dgf_2d_ew[200] = "/data/calet/abundance_estimate/rigidity/code/geofac_dir/dgf_2d_bin1_UH_geom.txt";

    double gamma_2d_g_b1[180][180];
    double delta_2d_g_b1[180][180];
    double dgf_2d_g_b1[180][180];
    double dgf_2d_g_b1_adj[180][180];
    
    double theta=0;
    
    int gamma_2d_g_b1_count=0;
    int delta_2d_g_b1_count=0;
    
    in_1.open(dgf_2d_ew);
	while(1) {
		in_1 >> gamma_2d_g_b1[gamma_2d_g_b1_count][delta_2d_g_b1_count]
        >> delta_2d_g_b1[gamma_2d_g_b1_count][delta_2d_g_b1_count]
        >> dgf_2d_g_b1[gamma_2d_g_b1_count][delta_2d_g_b1_count];
        if(angle_dgf_ew_b1[dgf_ew_b1_count] < 45 || angle_dgf_ew_b1[dgf_ew_b1_count] > 135) {
            geo_factor_dgf_ew_b1_adj[dgf_ew_b1_count] = 0.75*geo_factor_dgf_ew_b1[dgf_ew_b1_count];
        }
        else {
            geo_factor_dgf_ew_b1_adj[dgf_ew_b1_count] = geo_factor_dgf_ew_b1[dgf_ew_b1_count];
        }
		if(in_1.eof() || !in_1.good()) break;
		
		
		dgf_ew_b1_count++;
	}
	in_1.close();
	in_1.clear();
    
    */
    
    // East-West angle (gamma) versus theta 
    
    //char geom_file_g_theta_b1[200]="/data/calet/abundance_estimate/rigidity/code/geofac_dir/dgf_2d_g_theta_bin1_UH_test.txt";
    
    char geom_file_g_theta_b1_1[200]="/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/TIGERISS_Abundances/angle/dgf_2d_g_theta_bin1";
    
    char geom_file_g_theta_b1_2[200]="/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/TIGERISS_Abundances/angle/dgf_2d_g_theta_bin1";
    
    
    if(argc>=5) {
        strcpy(geom_file_g_theta_b1_1, source_dir);
        strcat(geom_file_g_theta_b1_1, "dgf_2d_g_theta_bin1");
        strcat(geom_file_g_theta_b1_1, name_ext1);
        strcat(geom_file_g_theta_b1_1, ".txt");
        strcpy(geom_file_g_theta_b1_2, source_dir);
        strcat(geom_file_g_theta_b1_2, "dgf_2d_g_theta_bin1");
        strcat(geom_file_g_theta_b1_2, name_ext2);
        strcat(geom_file_g_theta_b1_2, ".txt");
    }
    fprintf(stderr, "geom_file_g_theta_b1_1 = %s\n",geom_file_g_theta_b1_1);
    fprintf(stderr, "geom_file_g_theta_b1_2 = %s\n",geom_file_g_theta_b1_2);
    
	int g_dgf_count_theta_b1=0;
	
	double gamma_g_dgf_theta_b1;
	double theta_g_dgf_theta_b1;
	double geo_factor_g_dgf_theta_b1;
	double geo_factor_g_dgf_theta_b1_adj;
    
    int g_index=0;
    int t_index=0;
    double dgf_g_t[180][90];
    double dgf_g_t_adj[180][90];
	
    for (g_index=0; g_index < 180; g_index++) {
        for (t_index=0; t_index < 90; t_index++) {
            dgf_g_t[g_index][t_index] = 0.0;
            dgf_g_t_adj[g_index][t_index] = 0.0;
        }
    }
    
    //sprintf(geom_file,"%sdgf_2d_g_theta_bin1%s.txt",dir,extension);
    
	
	in_1.open(geom_file_g_theta_b1_1);
	while(1) {
		in_1 >> gamma_g_dgf_theta_b1 >> theta_g_dgf_theta_b1 >> geo_factor_g_dgf_theta_b1;
        if(theta_g_dgf_theta_b1 > 45) {
            geo_factor_g_dgf_theta_b1_adj = 0.75*geo_factor_g_dgf_theta_b1;
        }
        else {
            geo_factor_g_dgf_theta_b1_adj = geo_factor_g_dgf_theta_b1;
        }
        g_index=(int)(gamma_g_dgf_theta_b1 - 0.5);
        t_index=(int)(theta_g_dgf_theta_b1 - 0.5);
        
        dgf_g_t[g_index][t_index]+=0.5*geo_factor_g_dgf_theta_b1;
        dgf_g_t_adj[g_index][t_index]+=0.5*geo_factor_g_dgf_theta_b1_adj;
        
        
		if(in_1.eof() || !in_1.good()) break;
		g_dgf_count_theta_b1++;
	}
	in_1.close();
	in_1.clear();
	
    g_dgf_count_theta_b1 = 0;
    
    in_1.open(geom_file_g_theta_b1_2);
    while(1) {
        in_1 >> gamma_g_dgf_theta_b1 >> theta_g_dgf_theta_b1 >> geo_factor_g_dgf_theta_b1;
        if(theta_g_dgf_theta_b1 > 45) {
            geo_factor_g_dgf_theta_b1_adj = 0.75*geo_factor_g_dgf_theta_b1;
        }
        else {
            geo_factor_g_dgf_theta_b1_adj = geo_factor_g_dgf_theta_b1;
        }
        g_index=(int)(gamma_g_dgf_theta_b1 - 0.5);
        t_index=(int)(theta_g_dgf_theta_b1 - 0.5);
        
        dgf_g_t[g_index][t_index]+=0.5*geo_factor_g_dgf_theta_b1;
        dgf_g_t_adj[g_index][t_index]+=0.5*geo_factor_g_dgf_theta_b1_adj;
        
        
        if(in_1.eof() || !in_1.good()) break;
        g_dgf_count_theta_b1++;
    }
    in_1.close();
    in_1.clear();
    
    
    // Geomagnetic rigidity analysis
    
    char cut_450[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/TBLV450K.txt";
    
    char header_mag[100];
	int long_index=0;
	int lat_index=0;
    int el_index=0;

    double longitude[long_bins];
    
    struct Mag_Cut
	{
		double latitude;
		double cut_off[long_bins];
	} alt_450[lat_bins];
    
    
    // Reading geomagnetic cutoff data and calculating elemental minimum energies
	in_1.open(cut_450);
    
	in_1.getline(header_mag,100);
    
	fprintf(stderr,"%s\n",header_mag);
    
	for(long_index=0;long_index<25;long_index++)
    {
	  	in_1 >> longitude[long_index];
	  	// printf("%f\n",longitude[long_index]);
    }
	
	lat_index=0;
    
	while(1) {
        
		in_1 >> alt_450[lat_index].latitude;
		// fprintf(stderr,"\n%f ",alt_450[lat_index].latitude);
        if(in_1.eof() || !in_1.good()) break;
        
        for(long_index=0;long_index<long_bins;long_index++)
		{
			in_1 >> alt_450[lat_index].cut_off[long_index];
            // fprintf(stderr,"%f ",alt_450[lat_index].cut_off[long_index]);
		}
		lat_index++;
    }
    
    
	in_1.close();
	in_1.clear();
    

	
	
  
    
    
    const Int_t NCont = 99;
    
    /*
    const Int_t NRGBs = 5;
    //const Int_t NCont = 255;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
     */
    gStyle->SetNumberContours(NCont);
    
	
    
    
    
    // Read in geomagnetic latitude residence times
	char geo_lat_file[200] = "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/geomagnetic_latitude_orbit_fraction_SAA.txt";
	int g_m_l_count=0;
	double angle_g_m_l[90];
	double orbit_fraction_g_m_l[90];
	
	char tstring_h_g_m_l[100];
	sprintf(tstring_h_g_m_l,"ISS Geomagnetic Latitude Residence");
	
    //	TH1F *h_g_m_l = new TH1F("h_g_m_l_p",tstring_h_g_m_l_p,100,0,20);
    TH1F *h_g_m_l= new TH1F("h_g_m_l",tstring_h_g_m_l,90,0,90);
    //TH1* h_g_m_l = new TH1F("h_g_m_l",tstring_h_g_m_l,90,0,90);
	
    // h_g_m_l
	h_g_m_l->GetXaxis()->SetTitle("Geomagnetic Latitude (#lambda in degrees)");
	h_g_m_l->GetXaxis()->CenterTitle();
	h_g_m_l->GetYaxis()->SetTitle("Orbit Fraction");
	h_g_m_l->GetYaxis()->CenterTitle();
	h_g_m_l->GetYaxis()->SetTitleOffset(1.7);
	h_g_m_l->GetYaxis()->SetLabelFont(nicefont);
	h_g_m_l->GetYaxis()->SetTitleFont(nicefont);
	h_g_m_l->GetXaxis()->SetLabelFont(nicefont);
	h_g_m_l->GetXaxis()->SetTitleFont(nicefont);
	h_g_m_l->SetTitleOffset(1.3);
	h_g_m_l->SetLineWidth(2);
	
	in_1.open(geo_lat_file);
	while(1) {
		in_1 >> angle_g_m_l[g_m_l_count] >> orbit_fraction_g_m_l[g_m_l_count];
		if(in_1.eof() || !in_1.good()) break;
		fprintf(stderr,"%4.1f %.6f\n",angle_g_m_l[g_m_l_count],orbit_fraction_g_m_l[g_m_l_count]);
		h_g_m_l->Fill(angle_g_m_l[g_m_l_count],orbit_fraction_g_m_l[g_m_l_count]);
        
		g_m_l_count++;
	}
	in_1.close();
	in_1.clear();
	
	TCanvas *C_g_m_l = new TCanvas("C_g_m_l","ISS Geomagnetic Latitudes",0,0,800,800);
	C_g_m_l->SetLeftMargin(0.148477);
	C_g_m_l->SetRightMargin(0.0507614);
	C_g_m_l->SetBottomMargin(0.125654);
	C_g_m_l->SetTopMargin(0.0746073);
	h_g_m_l->GetXaxis()->SetRangeUser(0,90);
	h_g_m_l->Draw("HIST");
    

    
    // Interaction Correction Factors for TIGER-ISS Material stack
    
    // 2 Silicon strip detectors for top hodoscope and dE/dx charge
    
    // Aluminum top
    //
    
    // Aerogel Cherenkov light collection box
    
    // Acrylic Cherenkov light collection box
    
    // 2 Silicon strip detectors for top hodoscope and dE/dx charge
    
    

    // material densities
    double den_Al = 2.700; // g/cm^3
    double den_PET = 1.390; // g/cm^3
    double den_PVT = 1.030; // g/cm^3
    double den_PS = 0.040; // g/cm^3
    double den_Si = 2.329; // g/cm^3
    double den_SiO2 = 0.205; // g/cm^3
    double den_Epoxy = 1.135; // g/cm^3
    double den_PMI = 0.032; // g/cm^3
    double den_PTFE = 0.600; // g/cm^3
    double den_PE = 0.919; // g/cm^3
    double den_Kapton = 1.42; // g/cm^3
    double den_PU = 0.08;  // g/cm^3 polyurathane foam
    
    double width_PVT = 0.0;
    double width_Al = 0.0;
    double width_PET = 0.0;
    double width_PS = 0.0;
    double width_Si = 0.0;
    double width_SiO2 = 0.0;
    double width_Epoxy = 0.0;
    double width_PMI = 0.0;
    double width_PTFE = 0.0;
    double width_PE = 0.0;
    double width_Kapton = 0.0;
    
    
    // material areal densities to top of active instrument
    double x_PVT_top=0.68581;
    double x_PS_top = 0.08300;

    // material areal densities to bottom of instrument
    double x_PVT=0.0;
    double x_Al = 0.4536;
    double x_PET = 0.1022;
    double x_PS = 0.0;
    double x_Si = 0.9316;
    double x_SiO2 = 0.41;
    double x_Epoxy = 0.0;
    double x_PMI = 0.2794;
    double x_PTFE = 0.06;
    double x_PE = 0.00552;
    double x_Kapton = 0.01803;
    double x_PU = 0.0508;
    double x_PMMA = 1.4986;

    // material areal densities to bottom of second Si layer
    double x_PVT_Si2=0.0;
    double x_Al_Si2 = 0.14580;
    double x_PET_Si2 = 0.03405;
    double x_PS_Si2 = 0.0;
    double x_Si_Si2 = 0.46580;
    double x_SiO2_Si2 = 0.00000;
    double x_Epoxy_Si2 = 0.0;
    double x_PMI_Si2 = 0.08128;
    double x_PTFE_Si2 = 0.00000;
    double x_PE_Si2 = 0.00000;
    double x_Kapton_Si2 = 0.00000;
    double x_PU_Si2 = 0.00000;
    double x_PMMA_Si2 = 0.00000;


	
    // material properties
	
    // PVT: C9H10
	int n_PVT_H = 10;
	int n_PVT_C = 9;
	
    // PS: C8H8
	int n_PS_H = 8;
	int n_PS_C = 8;
	
    // PMMA: C5H8O2
	int n_PMMA_H = 8;
	int n_PMMA_C = 5;
	int n_PMMA_O = 2;
    
    // PMI: C8H11O2N
    int n_PMI_H = 11;
    int n_PMI_C = 8;
    int n_PMI_O = 2;
    int n_PMI_N = 1;
    
    // PTFE: C2F4
    int n_PTFE_C = 2;
    int n_PTFE_F = 4;
    
    // PE: C2H4
    int n_PE_C = 2;
    int n_PE_H = 4;
    
    // PU: C25H42N2O6
    int n_PU_C = 25;
    int n_PU_H = 42;
    int n_PU_N = 2;
    int n_PU_O = 6;
    
    // Kapton
    int n_Kapton_C = 22;
    int n_Kapton_H = 10;
    int n_Kapton_N = 2;
    int n_Kapton_O = 5;
    
    // PET: C10H8O4 Polyethylene terepthalate (Mylar)
    int n_PET_C = 10;
    int n_PET_H = 8;
    int n_PET_O = 4;
    
    // Silicon
    int n_Si = 1;
	
    // Glass: SiO2
	int n_Glass_Si = 1;
	int n_Glass_O = 2;
	
    // Aluminum: Al
	int n_Al = 1;
    
	
    // Correction Factors
	
	
    // Target Radii
    
    
    // Atmoic masses for interaction cross section and interaction path calculations
    
	double A_N=14.0067;
	double A_O=15.9994;
    double A_F=17.4228;
	double A_H=1.00794;
	double A_C=12.0107;
	double A_Si=28.0855;
	double A_Al=26.981538;
	
	double R_N=0;
	double R_O=0;
	double R_H=0;
	double R_C=0;
    double R_F=0;
	double R_Si=0;
	double R_Al=0;
    
    double R_Z=0; // Radius of interaction 

	
	R_N = 1.58*pow(A_N,0.281);
	R_O = 1.58*pow(A_O,0.281);
    //  R_H = 1.58*pow(A_H,0.281);  // forrmula does NOT work for H
	R_H = 1.034;
	R_C = 1.58*pow(A_C,0.281);
	R_Si = 1.58*pow(A_Si,0.281);
	R_Al = 1.58*pow(A_Al,0.281);
    R_F = 1.58*pow(A_F,0.281);
	
    
	double angle_fac=0.0;
	

    
    double sig_tot_H_eq_9 = 0.0;
    double sig_tot_C_eq_9 = 0.0;
    double sig_tot_N_eq_9 = 0.0;
    double sig_tot_O_eq_9 = 0.0;
    double sig_tot_F_eq_9 = 0.0;
    double sig_tot_Si_eq_9 = 0.0;
    double sig_tot_Al_eq_9 = 0.0;
    
    double lambda_PVT = 0.0;
    double lambda_PS = 0.0;
    double lambda_PMMA = 0.0;
    double lambda_SiO2 = 0.0;
    double lambda_Al = 0.0;
    double lambda_PMI = 0.0;
    double lambda_PTFE = 0.0;
    double lambda_Kapton = 0.0;
    double lambda_PE = 0.0;
    double lambda_PU = 0.0;
    double lambda_Si = 0.0;
    double lambda_PET = 0.0;


    // Material mean free paths
	
    // initial[index].R = 1.58*pow(initial[index].A,0.281);
	double interaction_Z[100];
	double interaction_survival_fraction[100];
    
    
    
    // Integral spectra files from write_integral_spectra.cpp
    // Derived from ACE spectra courtesy of Kelly Lave for 5 <= Z <= 32, and approximated by scaling the
    // Fe-spectrum with TOA abundances from HEAO and TIGER
    // Format from email from Kelly Lave: /data/tiger/calet/abundance_estimate/rigidity/tables/CRIS stuff.html
    char S_max_in[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/TIGERISS_Abundances/S_max_integral_spectra.txt";
    
    char S_min_in[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/TIGERISS_Abundances/S_min_integral_spectra.txt";
    
    
    char Z_char[2]; // file header dummy value - will read in "Z"
    char A_char[2]; // file header dummy value - will read in "A"
    char X_char[2]; // file header dummy value - will read in "X"
    
    int Z_S_max; // Nuclear charge
    double A_S_max; // Atomic number
    char cs_S_max[2]; // chemical symbol

    int Z_S_min; // Nuclear charge
    double A_S_min; // Atomic number
    char cs_S_min[2]; // chemical symbol
    
    double E_int_S_max[201]; // array of integral spectrum energy steps
    double F_int_S_max[201]; // array of integral spectrum fluxes
    // - last entry for extrapolated integral of spectrum tail out to CALET maximum energy
    double E_int_S_min[201]; // array of integral spectrum energy steps
    double F_int_S_min[201]; // array of integral spectrum fluxes
    // - last entry for extrapolated integral of spectrum tail out to CALET maximum energy
 
    // need integral iron spectra for elements assumed to have no flux in GCR that are 
    // produced by fragmentation in non-active top layers of detector
    double F_int_S_max_fe[201]; // array of integral spectrum fluxes   
    double F_int_S_min_fe[201]; // array of integral spectrum fluxes
    
    for(i=0; i<201; i++) {
      F_int_S_max_fe[i] = 0.0;
      F_int_S_min_fe[i] = 0.0;
    }


    
    
  
    in_1.open(S_max_in);
    fprintf(stderr,"S_max integral file\n");
    
    in_1 >> Z_char >> A_char >> X_char;
    
    fprintf(stderr,"%s %s %s\n",Z_char,A_char,X_char);

    
    for (i=0; i<201; i++) {
        in_1 >> E_int_S_max[i];
    }
    
    while (1) {
        in_1 >> Z_S_max >> A_S_max >> cs_S_max;
        fprintf(stderr,"Z: %2d A: %.2f %s\n",Z_S_max,A_S_max,cs_S_max);
        for (i=0; i<201; i++) {
            in_1 >> F_int_S_max[i];
        }
        if (Z_S_max==Z_iter) {
            break;
        }
        if (Z_S_max==26) {
          for (i=0; i<201; i++) {
            F_int_S_max_fe[i] = F_int_S_max[i];
          }
        }
        if(in_1.eof() || !in_1.good()) break;
    }
    
    
    in_1.close();
    in_1.clear();
    
    in_1.open(S_min_in);
    fprintf(stderr,"S_min integral file\n");

    
    in_1 >> Z_char >> A_char >> X_char;
    
    fprintf(stderr,"%s %s %s\n",Z_char,A_char,X_char);
    
    for (i=0; i<201; i++) {
        in_1 >> E_int_S_min[i];
    }
    
    while (1) {
        in_1 >> Z_S_min >> A_S_min >> cs_S_min;
        fprintf(stderr,"Z: %2d A: %.2f %s\n",Z_S_min,A_S_min,cs_S_min);
        for (i=0; i<201; i++) {
            in_1 >> F_int_S_min[i];
        }
        if (Z_S_min==Z_iter) {
            break;
        }
        if (Z_S_min==26) {
          for (i=0; i<201; i++) {
            F_int_S_min_fe[i] = F_int_S_min[i];
          }
        }

        if(in_1.eof() || !in_1.good()) break;
    }
    
    
    in_1.close();
    in_1.clear();
    
    
    double A=0; // Atomic mass of elment for spectral integration
    int Z=0; // Atomic number (nuclear charge) for spectral integration
    
    if (A_S_min == A_S_max) {
        A=A_S_min;
    }
    else {
        fprintf(stderr,"Ahhhhhh! The atomic masses for the S_min and S_max files do not agree!\n");
        fprintf(stderr,"A_S_min: %.2f A_S_max: %.2f\n",A_S_min,A_S_max);
    }

    
    if (Z_S_min == Z_S_max) {
        Z=Z_S_min;
    }
    else {
        fprintf(stderr,"Ahhhhhh! The atomic numbers for the S_min and S_max files do not agree!\n");
        fprintf(stderr,"Z_S_min: %d Z_S_max: %d\n",Z_S_min,Z_S_max);
    }
    
    
    // Interaction factor calculation
    if (Z_iter == 1) {
        R_Z = R_H;
    }
    else {
        R_Z = 1.58*pow(A,0.281);
    }
    
    
    // Generate Cross Sections
		
    sig_tot_H_eq_9 = 0.01*TMath::Pi()*pow((S_fac*R_H + S_fac*R_Z - 3.2),2);
    sig_tot_C_eq_9 = 0.01*TMath::Pi()*pow((S_fac*R_C + S_fac*R_Z - 3.2),2);
    sig_tot_N_eq_9 = 0.01*TMath::Pi()*pow((S_fac*R_N + S_fac*R_Z - 3.2),2);
    sig_tot_O_eq_9 = 0.01*TMath::Pi()*pow((S_fac*R_O + S_fac*R_Z - 3.2),2);
    sig_tot_F_eq_9 = 0.01*TMath::Pi()*pow((S_fac*R_F + S_fac*R_Z - 3.2),2);
    sig_tot_Si_eq_9 = 0.01*TMath::Pi()*pow((S_fac*R_Si + S_fac*R_Z - 3.2),2);
    sig_tot_Al_eq_9 = 0.01*TMath::Pi()*pow((S_fac*R_Al + S_fac*R_Z - 3.2),2);


    lambda_PVT = 1.6624*(n_PVT_H*A_H + n_PVT_C*A_C)/(n_PVT_H*sig_tot_H_eq_9 + n_PVT_C*sig_tot_C_eq_9);

    lambda_PS = 1.6624*(n_PS_H*A_H + n_PS_C*A_C)/(n_PS_H*sig_tot_H_eq_9 + n_PS_C*sig_tot_C_eq_9);

    lambda_PMMA = 1.6624*(n_PMMA_H*A_H + n_PMMA_C*A_C + n_PMMA_O*A_O)/
    (n_PMMA_H*sig_tot_H_eq_9 + n_PMMA_C*sig_tot_C_eq_9 + n_PMMA_O*sig_tot_O_eq_9);

    lambda_SiO2 = 1.6624*(n_Glass_Si*A_Si + n_Glass_O*A_O)/
    (n_Glass_Si*sig_tot_Si_eq_9 + n_Glass_O*sig_tot_O_eq_9);
    
    lambda_Al = 1.6624*(n_Al*A_Al)/
    (n_Al*sig_tot_Al_eq_9);
    
    lambda_PMI = 1.6624*(n_PMI_H*A_H + n_PMI_C*A_C + n_PMI_O*A_O + n_PMI_N*A_N)/
    (n_PMI_H*sig_tot_H_eq_9 + n_PMI_C*sig_tot_C_eq_9 + n_PMI_O*sig_tot_O_eq_9 + n_PMI_N*sig_tot_N_eq_9);
    
    lambda_PTFE = 1.6624*(n_PTFE_C*A_C + n_PTFE_F*A_F)/
    (n_PTFE_C*sig_tot_C_eq_9 + n_PTFE_F*sig_tot_F_eq_9);

    lambda_PE = 1.6624*(n_PE_H*A_H + n_PE_C*A_C)/
    (n_PE_H*sig_tot_H_eq_9 + n_PE_C*sig_tot_C_eq_9);

    lambda_PET = 1.6624*(n_PET_H*A_H + n_PET_C*A_C + n_PET_O*A_O)/
    (n_PET_H*sig_tot_H_eq_9 + n_PET_C*sig_tot_C_eq_9 + n_PET_O*sig_tot_O_eq_9);
    
    lambda_Kapton = 1.6624*(n_Kapton_H*A_H + n_Kapton_C*A_C + n_Kapton_N*A_N + n_Kapton_O*A_O)/
    (n_Kapton_H*sig_tot_H_eq_9 + n_Kapton_C*sig_tot_C_eq_9 +
     n_Kapton_N*sig_tot_N_eq_9 + n_Kapton_O*sig_tot_O_eq_9);
    
    lambda_PU = 1.6624*(n_PU_H*A_H + n_PU_C*A_C + n_PU_N*A_N + n_PU_O*A_O)/
    (n_PU_H*sig_tot_H_eq_9 + n_PU_C*sig_tot_C_eq_9 + n_PU_N*sig_tot_N_eq_9 + n_PU_O*sig_tot_O_eq_9);
    
    lambda_Si = 1.6624*(n_Si*A_Si)/(n_Si*sig_tot_Si_eq_9);
    



		
    fprintf(stderr,"A: %f lambda_PVT: %f\n",A,lambda_PVT);
    fprintf(stderr,"sig_tot_H_eq_9: %f\n",sig_tot_H_eq_9);
    fprintf(stderr,"sig_tot_C_eq_9: %f\n",sig_tot_C_eq_9);


    
    
    
    // Graphs of integral solar max and solar min to be evaluated at critical energies
    
    char t_string[100];
    sprintf(t_string,"Critical Energy for _{%d}%s",Z,cs_S_min);
    
    
    TGraph *gr_int_S_max = new TGraph(200,E_int_S_max,F_int_S_max);
    gr_int_S_max->SetMarkerStyle(1);
    gr_int_S_max->SetMarkerSize(1);
    // gr_int_S_max->SetTitle("Solar Maximum Integral Flux");
    gr_int_S_max->SetTitle(t_string);
    gr_int_S_max->GetXaxis()->SetTitle("MeV/nuc");
    gr_int_S_max->GetXaxis()->CenterTitle();
    gr_int_S_max->GetYaxis()->SetTitle("1/(cm^{2} s sr)");
    gr_int_S_max->GetYaxis()->CenterTitle();
    gr_int_S_max->GetXaxis()->SetTitleOffset(1.2);
    gr_int_S_max->GetYaxis()->SetTitleOffset(1.4);
    gr_int_S_max->GetYaxis()->SetLabelFont(nicefont);
    gr_int_S_max->GetYaxis()->SetTitleFont(nicefont);
    gr_int_S_max->GetXaxis()->SetLabelFont(nicefont);
    gr_int_S_max->GetXaxis()->SetTitleFont(nicefont);
    gr_int_S_max->GetXaxis()->SetRangeUser(1e1,1e5);
    // gr_int_S_max->GetYaxis()->SetRangeUser(1e-13,1e-2);
    
    
    
    TGraph *gr_int_S_min = new TGraph(200,E_int_S_min,F_int_S_min);
    gr_int_S_min->SetMarkerStyle(1);
    gr_int_S_min->SetMarkerSize(1);
    // gr_int_S_min->SetTitle("Solar Minimum Integral Flux");
    gr_int_S_min->SetTitle(t_string);
    gr_int_S_min->GetXaxis()->SetTitle("MeV/nuc");
    gr_int_S_min->GetXaxis()->CenterTitle();
    gr_int_S_min->GetYaxis()->SetTitle("1/(cm^{2} s sr)");
    gr_int_S_min->GetYaxis()->CenterTitle();
    gr_int_S_min->GetXaxis()->SetTitleOffset(1.2);
    gr_int_S_min->GetYaxis()->SetTitleOffset(1.4);
    gr_int_S_min->GetYaxis()->SetLabelFont(nicefont);
    gr_int_S_min->GetYaxis()->SetTitleFont(nicefont);
    gr_int_S_min->GetXaxis()->SetLabelFont(nicefont);
    gr_int_S_min->GetXaxis()->SetTitleFont(nicefont);
    gr_int_S_min->GetXaxis()->SetRangeUser(1e1,1e5);
    // gr_int_S_min->GetYaxis()->SetRangeUser(1e-16,1e-2);


    TGraph *gr_int_S_max_fe = new TGraph(200,E_int_S_max,F_int_S_max_fe);
    gr_int_S_max_fe->SetMarkerStyle(1);
    gr_int_S_max_fe->SetMarkerSize(1);
    // gr_int_S_max_fe->SetTitle("Solar Maximum Integral Flux");
    gr_int_S_max_fe->SetTitle(t_string);
    gr_int_S_max_fe->GetXaxis()->SetTitle("MeV/nuc");
    gr_int_S_max_fe->GetXaxis()->CenterTitle();
    gr_int_S_max_fe->GetYaxis()->SetTitle("1/(cm^{2} s sr)");
    gr_int_S_max_fe->GetYaxis()->CenterTitle();
    gr_int_S_max_fe->GetXaxis()->SetTitleOffset(1.2);
    gr_int_S_max_fe->GetYaxis()->SetTitleOffset(1.4);
    gr_int_S_max_fe->GetYaxis()->SetLabelFont(nicefont);
    gr_int_S_max_fe->GetYaxis()->SetTitleFont(nicefont);
    gr_int_S_max_fe->GetXaxis()->SetLabelFont(nicefont);
    gr_int_S_max_fe->GetXaxis()->SetTitleFont(nicefont);
    gr_int_S_max_fe->GetXaxis()->SetRangeUser(1e1,1e5);
    // gr_int_S_max_fe->GetYaxis()->SetRangeUser(1e-13,1e-2);
    
    
    
    TGraph *gr_int_S_min_fe = new TGraph(200,E_int_S_min,F_int_S_min_fe);
    gr_int_S_min_fe->SetMarkerStyle(1);
    gr_int_S_min_fe->SetMarkerSize(1);
    // gr_int_S_min_fe->SetTitle("Solar Minimum Integral Flux");
    gr_int_S_min_fe->SetTitle(t_string);
    gr_int_S_min_fe->GetXaxis()->SetTitle("MeV/nuc");
    gr_int_S_min_fe->GetXaxis()->CenterTitle();
    gr_int_S_min_fe->GetYaxis()->SetTitle("1/(cm^{2} s sr)");
    gr_int_S_min_fe->GetYaxis()->CenterTitle();
    gr_int_S_min_fe->GetXaxis()->SetTitleOffset(1.2);
    gr_int_S_min_fe->GetYaxis()->SetTitleOffset(1.4);
    gr_int_S_min_fe->GetYaxis()->SetLabelFont(nicefont);
    gr_int_S_min_fe->GetYaxis()->SetTitleFont(nicefont);
    gr_int_S_min_fe->GetXaxis()->SetLabelFont(nicefont);
    gr_int_S_min_fe->GetXaxis()->SetTitleFont(nicefont);
    gr_int_S_min_fe->GetXaxis()->SetRangeUser(1e1,1e5);
    // gr_int_S_min_fe->GetYaxis()->SetRangeUser(1e-16,1e-2);
    
    
    int g_m_l_index=0; // geomagneitc latitude index
	int a_index=0; // west-east azimuth angle index
	
	double lambda=0.0; // geomagnetic latitude
	double gamma=0.0; // west-east azimith angle
    double p_p[90][180];
    
    for (g_m_l_index=0; g_m_l_index<90; g_m_l_index++) {
		for (a_index=0; a_index<180; a_index++) {
			p_p[g_m_l_index][a_index]=0.0;
		}
	}
    
    
    for (g_m_l_index=0; g_m_l_index<90; g_m_l_index++) {
		for (a_index=0; a_index<180; a_index++) {
			gamma=a_index+0.5;
			lambda=angle_g_m_l[g_m_l_index];
            
            p_p[g_m_l_index][a_index]=60*Z*pow((1-sqrt(1-cos(gamma*TMath::Pi()/180.0)*pow(cos(lambda*TMath::Pi()/180.0),3.0)))/(cos(gamma*TMath::Pi()/180)*cos(lambda*TMath::Pi()/180)),2.0);
        }
    }
    double E_crit[90][180];
    int n_E_crit=0;
    double x_E_crit[90*180];
    double y_E_crit[90*180];
    double z_E_crit[90*180];
    
    double time_in_orbit=365*24*60*60*5;  // time in orbit in seconds for 5 years
    
    double S_max_events=0.0;
    double S_min_events=0.0;

    double S_max_events_tot=0.0;
    double S_min_events_tot=0.0;
    
    double S_max_events_int=0.0;
    double S_min_events_int=0.0;

    double S_max_events_int_Si2=0.0;
    double S_min_events_int_Si2=0.0;
    
    double S_max_events_int_top=0.0;
    double S_min_events_int_top=0.0;    
    
    double S_max_events_int_Si2_top=0.0;
    double S_min_events_int_Si2_top=0.0;

    double S_max_top_factor=0.0;
    double S_min_top_factor=0.0;

    
    double zenith; // incidence angle for trigger threshold energy calculation
    double E_Thresh; // Threshold energy for triggering
    Double_t *E_Thresh_p; // Threshold energy for triggering
    
    int theta_index=0;
    double theta_val=0;
    double I_factor[90];
    double I_factor_Si2[90];
    
    
    for (theta_index=0; theta_index<90; theta_index++) {
        
        theta_val=theta_index+0.5;
        
        angle_fac=1/cos((theta_val)*(TMath::Pi()/180.0));
        
        // x_PVT = angle_fac*(den_PVT*width_PVT);
        // fprintf(stderr,"x_PVT: %.2f\n",x_PVT);
        
        // I_factor[theta_index] = exp(-x_PVT/lambda_PVT);
        
        I_factor[theta_index] =
        exp(-x_PS*angle_fac/lambda_PS)*
        exp(-x_Al*angle_fac/lambda_Al)*
        exp(-x_Si*angle_fac/lambda_Si)*
        exp(-x_PET*angle_fac/lambda_PET)*
        exp(-x_PMI*angle_fac/lambda_PMI)*
        exp(-x_PMMA*angle_fac/lambda_PMMA)*
        exp(-x_Kapton*angle_fac/lambda_Kapton)*
        exp(-x_SiO2*angle_fac/lambda_SiO2)*
        exp(-x_PE*angle_fac/lambda_PE)*
        exp(-x_PU*angle_fac/lambda_PU)*
        exp(-x_PTFE*angle_fac/lambda_PTFE);

        I_factor_Si2[theta_index] =
        exp(-x_PS_Si2*angle_fac/lambda_PS)*
        exp(-x_Al_Si2*angle_fac/lambda_Al)*
        exp(-x_Si_Si2*angle_fac/lambda_Si)*
        exp(-x_PET_Si2*angle_fac/lambda_PET)*
        exp(-x_PMI_Si2*angle_fac/lambda_PMI)*
        exp(-x_PMMA_Si2*angle_fac/lambda_PMMA)*
        exp(-x_Kapton_Si2*angle_fac/lambda_Kapton)*
        exp(-x_SiO2_Si2*angle_fac/lambda_SiO2)*
        exp(-x_PE_Si2*angle_fac/lambda_PE)*
        exp(-x_PU_Si2*angle_fac/lambda_PU)*
        exp(-x_PTFE_Si2*angle_fac/lambda_PTFE);
        
        fprintf(stderr,"theta_val: %.1f I_factor[%d]: %.6f I_factor_Si2[%d]: %.6f\n",
                        theta_val,theta_index,I_factor[theta_index],theta_index,I_factor_Si2[theta_index]);
    }
    

    
    for (g_m_l_index=0; g_m_l_index<90; g_m_l_index++) {
		for (a_index=0; a_index<180; a_index++) {
            E_crit[g_m_l_index][a_index]=0;
        }
    }
    
    
    for (g_m_l_index=0; g_m_l_index<90; g_m_l_index++) {
		for (a_index=0; a_index<180; a_index++) {
			gamma=a_index+0.5;
			lambda=angle_g_m_l[g_m_l_index];
            
            
            // Calculate the critical energy for element in MeV/nucleon
            E_crit[g_m_l_index][a_index] = -amu_MeV +
            sqrt(pow(amu_MeV,2.0) + pow(p_p[g_m_l_index][a_index]*1e3,2.0)/pow(A,2.0));
            
            
            // Get the threshold energy for element at zenith angle to trigger TIGERISS SOX
            if(a_index <= 90) {
                zenith = 90.5 - a_index*1.0;
            }     
            if(a_index > 90) {
                zenith = a_index*1.0 - 89.5;
            }      
            E_Thresh=gt->Interpolate(Z_iter,zenith);
            //fprintf(stderr,"Z: %2d zenith: %2.2f E_Thresh: %f\n",Z_iter,zenith,E_Thresh);

                
                
                 
                S_max_events_tot+=(gr_int_S_max->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_max[200])*
                orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*blah[a_index];
                
                S_min_events_tot+=(gr_int_S_min->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_min[200])*
                orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*blah[a_index];
            
            
            
            //if (E_crit[g_m_l_index][a_index] >= E_gf_thresh_Fe) { 
            if (E_crit[g_m_l_index][a_index] >= E_Thresh) {
                // Wherever the critical energy is above the threshold energy add the integral flux
                // times the East-West geometry factor times the total time in orbit * interaction survival
                // probability


                
                S_max_events+=(gr_int_S_max->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_max[200])*
                orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*blah[a_index];
                
                S_min_events+=(gr_int_S_min->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_min[200])*
                orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*blah[a_index];
                
                // Got to include interactions
                for (theta_index=0; theta_index<90; theta_index++) {
                    

                    S_max_events_int+=(gr_int_S_max->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_max[200])*
                        orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                    dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index];
                
                 
                    S_min_events_int+=(gr_int_S_min->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_min[200])*
                    orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                    dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index];   

                    // Tracking for interaction losses only to bottom of the second top Si detector
                    S_max_events_int_Si2+=(gr_int_S_max->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_max[200])*
                        orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                    dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index];
                
                 
                    S_min_events_int_Si2+=(gr_int_S_min->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_min[200])*
                    orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                    dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]; 
 



                    // check if integral spectra are nonzero
                    if(S_max_abund_in[index_S_max_rel_in] > 0) {

                      // Find average scaling factors for angle bins centered on half integers
                      S_max_top_factor = 0.5*((S_max_delta_abund[index_S_max_rel_in][theta_index] + S_max_abund_in[index_S_max_rel_in])/
                        S_max_abund_in[index_S_max_rel_in] + 
                        (S_max_delta_abund[index_S_max_rel_in][theta_index+1] + S_max_abund_in[index_S_max_rel_in])/
                        S_max_abund_in[index_S_max_rel_in]);

                      S_min_top_factor = 0.5*((S_min_delta_abund[index_S_min_rel_in][theta_index] + S_min_abund_in[index_S_min_rel_in])/
                        S_min_abund_in[index_S_min_rel_in] + 
                        (S_min_delta_abund[index_S_min_rel_in][theta_index+1] + S_min_abund_in[index_S_min_rel_in])/
                        S_min_abund_in[index_S_min_rel_in]);

                      S_max_events_int_top+=(gr_int_S_max->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_max[200])*
                          orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index]*S_max_top_factor;
                
                 
                      S_min_events_int_top+=(gr_int_S_min->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_min[200])*
                      orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index]*S_min_top_factor;  

                      // Tracking for interaction losses only to bottom of the second top Si detector
                      S_max_events_int_Si2_top+=(gr_int_S_max->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_max[200])*
                          orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]*S_max_top_factor;
                
                 
                      S_min_events_int_Si2_top+=(gr_int_S_min->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_min[200])*
                      orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]*S_min_top_factor;
              
                    }
                    // Here I would add code to generate scaled Fe spectra for the elements assumed to have no 
                    // flux in the GCR, but which are generated by interations in the TIGERISS non-active top layers
                    if(S_max_abund_in[index_S_max_rel_in] == 0 && S_max_Z[index_S_max_rel_in] > 26) {

                      // Find average scaling factors for angle bins centered on half integers
                      S_max_top_factor = 0.5*((S_max_delta_abund[index_S_max_rel_in][theta_index]/S_max_abund_in[index_S_max_rel_in_fe]) + 
                        (S_max_delta_abund[index_S_max_rel_in][theta_index+1]/S_max_abund_in[index_S_max_rel_in_fe]));

                      S_min_top_factor = 0.5*((S_min_delta_abund[index_S_min_rel_in][theta_index]/S_min_abund_in[index_S_min_rel_in_fe]) + 
                        (S_min_delta_abund[index_S_min_rel_in][theta_index+1]/S_min_abund_in[index_S_min_rel_in_fe]));



                      S_max_events_int_top+=(gr_int_S_max_fe->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_max[200])*
                          orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index]*S_max_top_factor;
                
                 
                      S_min_events_int_top+=(gr_int_S_min_fe->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_min[200])*
                      orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index]*S_min_top_factor;  

                      // Tracking for interaction losses only to bottom of the second top Si detector
                      S_max_events_int_Si2_top+=(gr_int_S_max_fe->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_max[200])*
                          orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]*S_max_top_factor;
                
                 
                      S_min_events_int_Si2_top+=(gr_int_S_min_fe->Eval(E_crit[g_m_l_index][a_index]) + F_int_S_min[200])*
                      orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]*S_min_top_factor;


                    }
                    
                }
                
            }
  
                // ADD the ELSE!!! If e_crit < E_thresh evaluate the integrals for E_thresh
            if (E_crit[g_m_l_index][a_index] < E_Thresh) {
                // Wherever the critical energy is above the threshold energy add the integral flux
                // times the East-West geometry factor times the total time in orbit * interaction survival
                // probability

                
                S_max_events+=(gr_int_S_max->Eval(E_Thresh) + F_int_S_max[200])*
                orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*blah[a_index];
                
                S_min_events+=(gr_int_S_min->Eval(E_Thresh) + F_int_S_min[200])*
                orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*blah[a_index];
                
                // Got to include interactions
                for (theta_index=0; theta_index<90; theta_index++) {
                    

                    S_max_events_int+=(gr_int_S_max->Eval(E_Thresh) + F_int_S_max[200])*
                        orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                    dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index];
                
                 
                    S_min_events_int+=(gr_int_S_min->Eval(E_Thresh) + F_int_S_min[200])*
                    orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                    dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index];   

                    // Tracking for interaction losses only to bottom of the second top Si detector
                    S_max_events_int_Si2+=(gr_int_S_max->Eval(E_Thresh) + F_int_S_max[200])*
                        orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                    dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index];
                
                 
                    S_min_events_int_Si2+=(gr_int_S_min->Eval(E_Thresh) + F_int_S_min[200])*
                    orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                    dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]; 
 



                    // check if integral spectra are nonzero
                    if(S_max_abund_in[index_S_max_rel_in] > 0) {

                      // Find average scaling factors for angle bins centered on half integers
                      S_max_top_factor = 0.5*((S_max_delta_abund[index_S_max_rel_in][theta_index] + S_max_abund_in[index_S_max_rel_in])/
                        S_max_abund_in[index_S_max_rel_in] + 
                        (S_max_delta_abund[index_S_max_rel_in][theta_index+1] + S_max_abund_in[index_S_max_rel_in])/
                        S_max_abund_in[index_S_max_rel_in]);

                      S_min_top_factor = 0.5*((S_min_delta_abund[index_S_min_rel_in][theta_index] + S_min_abund_in[index_S_min_rel_in])/
                        S_min_abund_in[index_S_min_rel_in] + 
                        (S_min_delta_abund[index_S_min_rel_in][theta_index+1] + S_min_abund_in[index_S_min_rel_in])/
                        S_min_abund_in[index_S_min_rel_in]);

                      S_max_events_int_top+=(gr_int_S_max->Eval(E_Thresh) + F_int_S_max[200])*
                          orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index]*S_max_top_factor;
                
                 
                      S_min_events_int_top+=(gr_int_S_min->Eval(E_Thresh) + F_int_S_min[200])*
                      orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index]*S_min_top_factor;  

                      // Tracking for interaction losses only to bottom of the second top Si detector
                      S_max_events_int_Si2_top+=(gr_int_S_max->Eval(E_Thresh) + F_int_S_max[200])*
                          orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]*S_max_top_factor;
                
                 
                      S_min_events_int_Si2_top+=(gr_int_S_min->Eval(E_Thresh) + F_int_S_min[200])*
                      orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]*S_min_top_factor;
              
                    }
                    // Here I would add code to generate scaled Fe spectra for the elements assumed to have no 
                    // flux in the GCR, but which are generated by interations in the TIGERISS non-active top layers
                    if(S_max_abund_in[index_S_max_rel_in] == 0 && S_max_Z[index_S_max_rel_in] > 26) {

                      // Find average scaling factors for angle bins centered on half integers
                      S_max_top_factor = 0.5*((S_max_delta_abund[index_S_max_rel_in][theta_index]/S_max_abund_in[index_S_max_rel_in_fe]) + 
                        (S_max_delta_abund[index_S_max_rel_in][theta_index+1]/S_max_abund_in[index_S_max_rel_in_fe]));

                      S_min_top_factor = 0.5*((S_min_delta_abund[index_S_min_rel_in][theta_index]/S_min_abund_in[index_S_min_rel_in_fe]) + 
                        (S_min_delta_abund[index_S_min_rel_in][theta_index+1]/S_min_abund_in[index_S_min_rel_in_fe]));



                      S_max_events_int_top+=(gr_int_S_max_fe->Eval(E_Thresh) + F_int_S_max[200])*
                          orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index]*S_max_top_factor;
                
                 
                      S_min_events_int_top+=(gr_int_S_min_fe->Eval(E_Thresh) + F_int_S_min[200])*
                      orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor[theta_index]*S_min_top_factor;  

                      // Tracking for interaction losses only to bottom of the second top Si detector
                      S_max_events_int_Si2_top+=(gr_int_S_max_fe->Eval(E_Thresh) + F_int_S_max[200])*
                          orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]*S_max_top_factor;
                
                 
                      S_min_events_int_Si2_top+=(gr_int_S_min_fe->Eval(E_Thresh) + F_int_S_min[200])*
                      orbit_fraction_g_m_l[g_m_l_index]*time_in_orbit*
                      dgf_g_t_adj[a_index][theta_index]*I_factor_Si2[theta_index]*S_min_top_factor;


                    }
                    
                }
                
            }
                
/*
 
                fprintf(stderr,"gr_int_S_max->Eval(E_crit[%d][%d]: %e\n",g_m_l_index,a_index,
                        gr_int_S_max->Eval(E_crit[g_m_l_index][a_index]));
                fprintf(stderr,"F_int_S_max[200]: %e\n",F_int_S_max[200]);
                fprintf(stderr,"orbit_fraction_g_m_l[%d]: %f\n",g_m_l_index,orbit_fraction_g_m_l[g_m_l_index]);
 */
//                fprintf(stderr,"geo_factor_dgf_ew_b1[%d]: %f\n",a_index,blah[a_index]);
                
            
            
         /*
            x_E_crit[n_E_crit]=gamma;
            y_E_crit[n_E_crit]=lambda;
            z_E_crit[n_E_crit]=E_crit[g_m_l_index][a_index];
            
            fprintf(stderr,"x_E_crit[%d]: %f y_E_crit[%d]: %f z_E_crit[%d]: %f\n",
                    n_E_crit,x_E_crit[n_E_crit],n_E_crit,y_E_crit[n_E_crit],n_E_crit,z_E_crit[n_E_crit]);
            
            n_E_crit++;
          */
        }
    }
    
	fprintf(stderr,"n_E_crit: %d\n",n_E_crit);
    fprintf(stderr,"S_max_events: %f\n",S_max_events);
    fprintf(stderr,"S_min_events: %f\n",S_min_events);
    fprintf(stderr,"S_max_events_int: %f\n",S_max_events_int);
    fprintf(stderr,"S_min_events_int: %f\n",S_min_events_int);
    fprintf(stderr,"S_max_events_int_top: %f\n",S_max_events_int_top);
    fprintf(stderr,"S_min_events_int_top: %f\n",S_min_events_int_top);
    fprintf(stdout,"%2d %e %e %e ",Z,0.5*(S_max_events_int + S_min_events_int),S_max_events_int,S_min_events_int);
    fprintf(stdout,"%e %e %e ",0.5*(S_max_events_int_top + S_min_events_int_top),S_max_events_int_top,S_min_events_int_top);
    fprintf(stdout,"%e %e %e ",0.5*(S_max_events_int_Si2 + S_min_events_int_Si2),S_max_events_int_Si2,S_min_events_int_Si2);
    fprintf(stdout,"%e %e %e\n",0.5*(S_max_events_int_Si2_top + S_min_events_int_Si2_top),S_max_events_int_Si2_top,S_min_events_int_Si2_top);
    fprintf(stderr,"%d %f %f %f\n",Z,0.5*(S_max_events_tot + S_min_events_tot),S_max_events_tot,S_min_events_tot);


    
/*
    
	TGraph2D *gr_E_crit = new TGraph2D(n_E_crit,x_E_crit,y_E_crit,z_E_crit);
	
	//gr_E_crit->SetMarkerStyle(6);
	//gr_E_crit->SetMarkerSize(1);
	gr_E_crit->SetTitle(t_string);
	gr_E_crit->GetXaxis()->SetTitle("East-West Angle (degrees)");
	gr_E_crit->GetXaxis()->CenterTitle();
	gr_E_crit->GetYaxis()->SetTitle("Geomagnetic Latitude (degrees)");
	gr_E_crit->GetYaxis()->CenterTitle();
	gr_E_crit->GetXaxis()->SetTitleOffset(1.2);
	gr_E_crit->GetYaxis()->SetTitleOffset(1.2);
	gr_E_crit->GetYaxis()->SetLabelFont(nicefont);
	gr_E_crit->GetYaxis()->SetTitleFont(nicefont);
	gr_E_crit->GetXaxis()->SetLabelFont(nicefont);
	gr_E_crit->GetXaxis()->SetTitleFont(nicefont);
	gr_E_crit->GetXaxis()->SetRangeUser(0,180);
	gr_E_crit->GetYaxis()->SetRangeUser(0,90);
    
     //gr_E_crit->SetContour(NCont);

    

    
    TCanvas *c_E_crit = new TCanvas("c_E_crit",t_string,1200,800);
    gr_E_crit->Draw("COLZ");
    
*/    
    
/*    
	//gr_E_crit->SetMarkerStyle(6);
	//gr_E_crit->SetMarkerSize(1);
	gt->SetTitle("; Z; Zenith Angle (degrees); Energy (MeV/nuc)");
	gt->GetXaxis()->SetTitle("Z");
	gt->GetXaxis()->CenterTitle();
	gt->GetYaxis()->SetTitle("Zenith Angle)");
	gt->GetYaxis()->CenterTitle();
	gt->GetZaxis()->SetTitle("Theshold Energy (MeV/nuc)");
	gt->GetZaxis()->CenterTitle();
	gt->GetXaxis()->SetTitleOffset(1.4);
	gt->GetYaxis()->SetTitleOffset(1.4);
	gt->GetZaxis()->SetTitleOffset(1.4);
	gt->GetXaxis()->SetLabelFont(nicefont);
	gt->GetXaxis()->SetTitleFont(nicefont);
	gt->GetYaxis()->SetLabelFont(nicefont);
	gt->GetYaxis()->SetTitleFont(nicefont);
	gt->GetZaxis()->SetLabelFont(nicefont);
	gt->GetZaxis()->SetTitleFont(nicefont);
	gt->GetXaxis()->SetRangeUser(0,90);
	gt->GetYaxis()->SetRangeUser(0,90);
    
     //gt->SetContour(NCont);

    

    
    TCanvas *c_E_Thresh = new TCanvas("c_E_Thresh",t_string,1200,800);
   gt->Draw("LEGO");
    
*/    

    
//	theApp.Run();
    
}

    


Double_t Integralfunc(TGraph *gr, Double_t a, Double_t b, Option_t
                      *option, Double_t epsilon)
{
    //Compute the integral of TGraph *gr from a to b.
    //if option "5 is specified a quintic spline is used to approximate
    //the graph, otherwise a cubic spline is used.
    //The method used is described in TF1::Integral
    //This function assumes that
    //  -the points in the graph are sorted along X
    //  -there is only one value of Y for a given value of X in [a,b].
    //if ga = gr(a) and gb=gr(b), the value returned excludes the area
    //outside the line connecting (a,ga) and (b,gb) by default.
    //To also include this area, specify the option "out"
    
    const Double_t kHF = 0.5;
    const Double_t kCST = 5./1000;
    
    Double_t x[12] = { 0.96028985649753623,  0.79666647741362674,
        0.52553240991632899,  0.18343464249564980,
        0.98940093499164993,  0.94457502307323258,
        0.86563120238783174,  0.75540440835500303,
        0.61787624440264375,  0.45801677765722739,
        0.28160355077925891,  0.09501250983763744};
    
    Double_t w[12] = { 0.10122853629037626,  0.22238103445337447,
        0.31370664587788729,  0.36268378337836198,
        0.02715245941175409,  0.06225352393864789,
        0.09515851168249278,  0.12462897125553387,
        0.14959598881657673,  0.16915651939500254,
        0.18260341504492359,  0.18945061045506850};
    
    Double_t h, aconst, bb, aa, c1, c2, u, s8, s16, f1, f2;
    Double_t xx;
    Int_t i;
    h = 0;
    if (b == a) return h;
    //Evaluate a spline smoothing this graph
    TSpline *spline = 0;
    TString opt = option;
    opt.ToLower();
    if (opt.Contains("5")) spline= new TSpline3("integral",gr);
    else                   spline= new TSpline5("integral",gr);
    Bool_t absValue = kFALSE;
    if (opt.Contains("abs")) absValue = kTRUE;
    //Evaluate the integral
    Double_t ga = gr->Eval(a,spline);
    Double_t gb = gr->Eval(b,spline);
    aconst = kCST/TMath::Abs(b-a);
    bb = a;
CASE1:
    aa = bb;
    bb = b;
CASE2:
    c1 = kHF*(bb+aa);
    c2 = kHF*(bb-aa);
    s8 = 0;
    for (i=0;i<4;i++) {
        u     = c2*x[i];
        xx = c1+u;
        f1    = gr->Eval(xx,spline);
        if (absValue) f1 = TMath::Abs(f1);
        xx = c1-u;
        f2    = gr->Eval(xx,spline);
        if (absValue) f2 = TMath::Abs(f2);
        s8   += w[i]*(f1 + f2);
    }
    s16 = 0;
    for (i=4;i<12;i++) {
        for (i=4;i<12;i++) {
            u     = c2*x[i];
            xx = c1+u;
            f1    = gr->Eval(xx,spline);
            if (absValue) f1 = TMath::Abs(f1);
            xx = c1-u;
            f2    = gr->Eval(xx,spline);
            if (absValue) f2 = TMath::Abs(f2);
            s16  += w[i]*(f1 + f2);
        }
        s16 = c2*s16;
        if (TMath::Abs(s16-c2*s8) <= epsilon*(1. + TMath::Abs(s16))) {
            h += s16;
            if(bb != b) goto CASE1;
        } else {
            bb = c1;
            if(1. + aconst*TMath::Abs(c2) != 1) goto CASE2;
            h = s8;
        }
        delete spline;
        if (!opt.Contains("out")) return h;
        else                      return h - 0.5*(ga+gb)/(b-a);
    }
}

TGraph2D* Thresholds()
{
//Read thresholds txt file and return TGraph2D for interpolation
  //Use like this:
  //TGraph2D *gt = Thresholds();
  //get->Interpolate(Z,zenith);
  
//  string fname = "Thresholds_SOX_20230307.txt";
  string fname = "input/EThresh.txt";
  const int MAXZ=80;
  const int MAXZEN = 70;
  double ch[MAXZ], p0[MAXZ], p1[MAXZ], p2[MAXZ], p3[MAXZ];
  string line;
  
  ifstream ifile(fname, ios::in);

  TGraph2D* dd;
  if (!ifile.is_open()) {
    cout << "unable to open file " << fname << endl;
    return dd;
  } 

  int iknt = 0;
  while (!ifile.eof()) {
      getline(ifile, line);
      if (line[0] == '#'  || line == "") continue; 
      stringstream ss(line);
      ss >> ch[iknt] >> p0[iknt] >> p1[iknt] >> p2[iknt] >> p3[iknt];
      cout << ch[iknt] << "\t" << p0[iknt] << "\t" <<  p1[iknt] << "\t" <<  p2[iknt] << "\t" <<  p3[iknt] << "\t" << endl;
      iknt++;
  }
  
  //prepare arrays of (Z,zenith,Ethr)
  // Z=charge
  // zenith angle in degrees
  // Ethr in MeV/nuc
  
  const int MAXN = MAXZ*MAXZEN;
  double Z[MAXN], zen[MAXN], Ethr[MAXN];
  
  int npt = 0;
  TF1 *fp = new TF1("fp","pol3",0,90);
  for (int iZ=0; iZ<iknt; iZ++) //charge
  {
    fp->SetParameters(p0[iZ],p1[iZ],p2[iZ],p3[iZ]);
    for (int izen = 0; izen<=MAXZEN; izen++)
    {
      Z[npt] = ch[iZ];
      zen[npt] = izen;
      Ethr[npt] = fp->Eval(zen[npt]);
      cout << npt << "\t" << Z[npt] << "\t" << zen[npt] << "\t" << Ethr[npt] << endl; 
      npt++;
      if (npt > MAXN) cout << "DANGER! npt, MAXN: " << npt << "\t" << MAXN << endl;
    }
  }
    
  TGraph2D *g0 = new TGraph2D(npt,Z,zen,Ethr); 
  g0->SetTitle("E_{Thresh} for SOX short stack;Charge;Zenith angle (deg);Threshold energy (MeV/nuc)");
//  g0->Draw("surf");

  return g0;

} 
