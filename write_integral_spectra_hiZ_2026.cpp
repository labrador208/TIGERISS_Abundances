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
#include <TPaletteAxis.h>
#include <TGraphAsymmErrors.h>


#define ISS_inc 51.6445 // Inclination of ISS in degrees
#define amu_MeV 931.49432 // unified mass unit in MeV/c^2
#define c_vac 2.99792458e8 // speed of light in vacuum in m/s (exact)
#define E_CALET_hi 10e8 // 100 TeV maximum energy for CALET

#define tail_int 20 // number of points at the end of the ACE spectra to use in fitting power law tail

#define long_bins 25 // number of longitude bins in geomagnetic model: 15 degree increments from 0 to 360 (inclusive)
#define long_bins_calc 24 // limit for actual calculation, as 0 360 degrees are the same
#define lat_bins 37  // number of latitude bins in geomagnetic model: 5 degree increments from 90 to -90 degrees 

#define long_bins_1 360 // number of longitude bins for 1 degree pitch
#define lat_bins_1 180 // number of latitude bins for 1 degree pitch

#define geo_factor_m 0.4397  // 0.12   // Geometric Factor CALET m^2-sr
#define geo_factor_cm 4.397e3 // 1.2e3 // Geometric Factor of CALET cm^2-sr
#define geo_factor_m_s  0.12   // Full Instrument Geometric Factor CALET m^2-sr
#define geo_factor_cm_s  1.2e3 // Full Instrument Geometric Factor of CALET cm^2-sr

// Here I'm just making a WAG
// 2016-03-16 On Bob's suggestion I am going from 0.8/0.2 to 0.75/0.25
#define even_frac 0.75 // fractional abundance of HEAO for even nuclei
#define odd_frac 0.25  // fractional abundance of HEAO for odd nuclei

	// Interaction correction parameters
#define S_fac 1.277 // scaling factor: Rn = SRe
#define NA 6.02214179e23 // Avogadro's Number (particles/mole)

#define E_gf_thresh_Fe 600.0 // (MeV/nuc) Energy threshold between high and low acceptance geometry factor regimes for Fe


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





  TApplication theApp("App", &argc, argv);

  ifstream in_1, in_2;
  ofstream out_1;

  // Input files.
  char ab_input[200] =
		//	"/data/tiger/calet/abundance_estimate/rigidity/code/ab_input";
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/sim_abunds/ab_input_HiZ.txt";
  //  "/home/brian/brian/propagation/code/ab_input";

  char supertiger_toa[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/sim_abunds/supertiger_toa.txt";

  char uhcre[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/sim_abunds/UHCRE.txt";

  char spec_smax[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/spectra.1au.phi900.txt";

  char spec_smin[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/spectra.1au.phi325.txt";    

  char cut_450[200] =
		"/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/TBLV450K.txt"; 

  char abund_smin[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/solmin_relabund97.160.dat";

	char abund_smax[200] =
//    "/data/calet/abundance_estimate/rigidity/tables/test.txt";
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/solmax_relabund.160.dat";

	char chem_s[3];
	
	int plot_Z=0;
	if(argc==2) plot_Z=atoi(argv[1]);
	
	switch (plot_Z)
    {
		case 1:
			sprintf(chem_s,"H");
			break;
		case 2:
			sprintf(chem_s,"He");
			break;
		case 3:
			sprintf(chem_s,"Li");
			break;
		case 4:
			sprintf(chem_s,"Be");
			break;
		case 5:
			sprintf(chem_s,"B");
			break;
		case 6:
			sprintf(chem_s,"C");
			break;
		case 7:
			sprintf(chem_s,"N");
			break;
		case 8:
			sprintf(chem_s,"O");
			break;
		case 9:
			sprintf(chem_s,"F");
			break;
		case 10:
			sprintf(chem_s,"Ne");
			break;
		case 11:
			sprintf(chem_s,"Na");
			break;
		case 12:
			sprintf(chem_s,"Mg");
			break;
		case 13:
			sprintf(chem_s,"Al");
			break;
		case 14:
			sprintf(chem_s,"Si");
			break;
		case 15:
			sprintf(chem_s,"P");
			break;
		case 16:
			sprintf(chem_s,"S");
			break;
		case 17:
			sprintf(chem_s,"Cl");
			break;
		case 18:
			sprintf(chem_s,"Ar");
			break;
		case 19:
			sprintf(chem_s,"K");
			break;
		case 20:
			sprintf(chem_s,"Ca");
			break;
		case 21:
			sprintf(chem_s,"Sc");
			break;
		case 22:
			sprintf(chem_s,"Ti");
			break;
		case 23:
			sprintf(chem_s,"V");
			break;
		case 24:
			sprintf(chem_s,"Cr");
			break;
		case 25:
			sprintf(chem_s,"Mn");
			break;
		case 26:
			sprintf(chem_s,"Fe");
			break;
		case 27:
			sprintf(chem_s,"Co");
			break;
		case 28:
			sprintf(chem_s,"Ni");
			break;
		case 29:
			sprintf(chem_s,"Cu");
			break;
		case 30:
			sprintf(chem_s,"Zn");
			break;
		case 31:
			sprintf(chem_s,"Ga");
			break;
		case 32:
			sprintf(chem_s,"Ge");
			break;
		case 33:
			sprintf(chem_s,"As");
			break;
		case 34:
			sprintf(chem_s,"Se");
			break;
		case 35:
			sprintf(chem_s,"Br");
			break;
		case 36:
			sprintf(chem_s,"Kr");
			break;
		case 37:
			sprintf(chem_s,"Rb");
			break;
		case 38:
			sprintf(chem_s,"Sr");
			break;
		case 39:
			sprintf(chem_s,"Y");
			break;
		case 40:
			sprintf(chem_s,"Zr");
			break;
		case 41:
			sprintf(chem_s,"Nb");
			break;
		case 42:
			sprintf(chem_s,"Mo");
			break;
		case 43:
			sprintf(chem_s,"Tc");
			break;
		case 44:
			sprintf(chem_s,"Ru");
			break;
		case 45:
			sprintf(chem_s,"Rh");
			break;
		case 46:
			sprintf(chem_s,"Pd");
			break;
		case 47:
			sprintf(chem_s,"Ag");
			break;
		case 48:
			sprintf(chem_s,"Cd");
			break;
		case 49:
			sprintf(chem_s,"In");
			break;
		case 50:
			sprintf(chem_s,"Sn");
			break;
		case 51:
			sprintf(chem_s,"Sb");
			break;
		case 52:
			sprintf(chem_s,"Te");
			break;
		case 53:
			sprintf(chem_s,"I ");
			break;
		case 54:
			sprintf(chem_s,"Xe");
			break;
		case 55:
			sprintf(chem_s,"Cs");
			break;
		case 56:
			sprintf(chem_s,"Ba");
			break;
		case 57:
			sprintf(chem_s,"La");
			break;
		case 58:
			sprintf(chem_s,"Ce");
			break;
		case 59:
			sprintf(chem_s,"Pr");
			break;
		case 60:
			sprintf(chem_s,"Nd");
			break;
		case 61:
			sprintf(chem_s,"Pm");
			break;
		case 62:
			sprintf(chem_s,"Sm");
			break;
		case 63:
			sprintf(chem_s,"Eu");
			break;
		case 64:
			sprintf(chem_s,"Gd");
			break;
		case 65:
			sprintf(chem_s,"Tb");
			break;
		case 66:
			sprintf(chem_s,"Dy");
			break;
		case 67:
			sprintf(chem_s,"Ho");
			break;
		case 68:
			sprintf(chem_s,"Er");
			break;
		case 69:
			sprintf(chem_s,"Tm");
			break;
		case 70:
			sprintf(chem_s,"Yb");
			break;
		case 71:
			sprintf(chem_s,"Lu");
			break;
		case 72:
			sprintf(chem_s,"Hf");
			break;
		case 73:
			sprintf(chem_s,"Ta");
			break;
		case 74:
			sprintf(chem_s,"W ");
			break;
		case 75:
			sprintf(chem_s,"Re");
			break;
		case 76:
			sprintf(chem_s,"Os");
			break;
		case 77:
			sprintf(chem_s,"Ir");
			break;
		case 78:
			sprintf(chem_s,"Pt");
			break;
		case 79:
			sprintf(chem_s,"Au");
			break;
		case 80:
			sprintf(chem_s,"Hg");
			break;
		case 81:
			sprintf(chem_s,"Tl");
			break;			
		case 82:
			sprintf(chem_s,"Pb");
			break;
		case 83:
			sprintf(chem_s,"Bi");
			break;
		case 84:
			sprintf(chem_s,"Po");
			break;
		case 85:
			sprintf(chem_s,"At");
			break;
		case 86:
			sprintf(chem_s,"Rn");
			break;
		case 87:
			sprintf(chem_s,"Fr");
			break;
		case 88:
			sprintf(chem_s,"Ra");
			break;
		case 89:
			sprintf(chem_s,"Ac");
			break;
		case 90:
			sprintf(chem_s,"Th");
			break;
		case 91:
			sprintf(chem_s,"Pa");
			break;
		case 92:
			sprintf(chem_s,"U ");
			break;
        case 93:
            sprintf(chem_s,"Np");
            break;
        case 94:
            sprintf(chem_s,"Pu");
            break;
        case 95:
            sprintf(chem_s,"Am");
            break;
        case 96:
            sprintf(chem_s,"Cm");
            break;
        case 97:
            sprintf(chem_s,"Bk");
            break;
        case 98:
            sprintf(chem_s,"Cf");
            break;
        case 99:
            sprintf(chem_s,"Es");
            break;
        case 100:
            sprintf(chem_s,"Fm");
            break;
    }
	
	
	
  struct element
  {
    int Z;
    double A;
    double SS_abund;
    double tiger;
    double HEAO_HNE;
    double HEAO_C2;
		double TOA_abund;
		double E_min_mag[long_bins][lat_bins];  // Element geomagnetic cutoff energy for longitude-latitude bins
		double int_S_max[long_bins][lat_bins]; 	// Element integrated spectra for energy geomagnetic energy cutoffs
		double int_S_min[long_bins][lat_bins]; 	// Element integrated spectra for energy geomagnetic energy cutoffs
	  double int_S_max_gf[long_bins][lat_bins]; 	// Element integrated spectra for energy geomagnetic energy cutoffs
	  double int_S_min_gf[long_bins][lat_bins]; 	// Element integrated spectra for energy geomagnetic energy cutoffs
	  double int_S_max_gf_s[long_bins][lat_bins]; 	// Element integrated spectra for energy geomagnetic energy cutoffs
	  double int_S_min_gf_s[long_bins][lat_bins]; 	// Element integrated spectra for energy geomagnetic energy cutoffs
	  double abund_S_max;											// Expected abundance of element from solar max
	  double abund_S_min;											// Expected abundance of element from solar min
	  double abund_S_max_gf;											// Expected abundance of element from solar max
	  double abund_S_min_gf;											// Expected abundance of element from solar min
	  double abund_S_max_dgf;											// Expected abundance of element from solar max
	  double abund_S_min_dgf;											// Expected abundance of element from solar min
	  double abund_S_max_dgf_noic;											// Expected abundance of element from solar max
	  double abund_S_min_dgf_noic;											// Expected abundance of element from solar min
		double abund_S_max_tot;									// Expected abundance of element from solar max with integrated power law tail
		double abund_S_min_tot;									// Expected abundance of element from solar min with integrated power law tail
		double abund_S_max_no_mag;
		double abund_S_min_no_mag;
	  double R; // radius for cross section model
	  double sig_tot_C_eq_9; // total interaction cross section on C
	  double sig_tot_H_eq_9; // total interaction cross section on H
	  double lambda_PVT; // mean free path in PVT
	  double I_factor; // Interaction factor
	  double E_gf_thresh;  // Energy threshold between high and low acceptance geometry factor regimes
	  double R_cut; // Rigidity cutoff at energy threshold for high and low acceptance geometry factor regimes
      char chem_s[3]; // Chemical symbol
	} el_dat[100];
	
    int j;
    for (j=0; j<=99; j++) {
        switch (j+1)
        {
            case 1:
                sprintf(el_dat[j].chem_s,"H");
                break;
            case 2:
                sprintf(el_dat[j].chem_s,"He");
                break;
            case 3:
                sprintf(el_dat[j].chem_s,"Li");
                break;
            case 4:
                sprintf(el_dat[j].chem_s,"Be");
                break;
            case 5:
                sprintf(el_dat[j].chem_s,"B");
                break;
            case 6:
                sprintf(el_dat[j].chem_s,"C");
                break;
            case 7:
                sprintf(el_dat[j].chem_s,"N");
                break;
            case 8:
                sprintf(el_dat[j].chem_s,"O");
                break;
            case 9:
                sprintf(el_dat[j].chem_s,"F");
                break;
            case 10:
                sprintf(el_dat[j].chem_s,"Ne");
                break;
            case 11:
                sprintf(el_dat[j].chem_s,"Na");
                break;
            case 12:
                sprintf(el_dat[j].chem_s,"Mg");
                break;
            case 13:
                sprintf(el_dat[j].chem_s,"Al");
                break;
            case 14:
                sprintf(el_dat[j].chem_s,"Si");
                break;
            case 15:
                sprintf(el_dat[j].chem_s,"P");
                break;
            case 16:
                sprintf(el_dat[j].chem_s,"S");
                break;
            case 17:
                sprintf(el_dat[j].chem_s,"Cl");
                break;
            case 18:
                sprintf(el_dat[j].chem_s,"Ar");
                break;
            case 19:
                sprintf(el_dat[j].chem_s,"K");
                break;
            case 20:
                sprintf(el_dat[j].chem_s,"Ca");
                break;
            case 21:
                sprintf(el_dat[j].chem_s,"Sc");
                break;
            case 22:
                sprintf(el_dat[j].chem_s,"Ti");
                break;
            case 23:
                sprintf(el_dat[j].chem_s,"V");
                break;
            case 24:
                sprintf(el_dat[j].chem_s,"Cr");
                break;
            case 25:
                sprintf(el_dat[j].chem_s,"Mn");
                break;
            case 26:
                sprintf(el_dat[j].chem_s,"Fe");
                break;
            case 27:
                sprintf(el_dat[j].chem_s,"Co");
                break;
            case 28:
                sprintf(el_dat[j].chem_s,"Ni");
                break;
            case 29:
                sprintf(el_dat[j].chem_s,"Cu");
                break;
            case 30:
                sprintf(el_dat[j].chem_s,"Zn");
                break;
            case 31:
                sprintf(el_dat[j].chem_s,"Ga");
                break;
            case 32:
                sprintf(el_dat[j].chem_s,"Ge");
                break;
            case 33:
                sprintf(el_dat[j].chem_s,"As");
                break;
            case 34:
                sprintf(el_dat[j].chem_s,"Se");
                break;
            case 35:
                sprintf(el_dat[j].chem_s,"Br");
                break;
            case 36:
                sprintf(el_dat[j].chem_s,"Kr");
                break;
            case 37:
                sprintf(el_dat[j].chem_s,"Rb");
                break;
            case 38:
                sprintf(el_dat[j].chem_s,"Sr");
                break;
            case 39:
                sprintf(el_dat[j].chem_s,"Y");
                break;
            case 40:
                sprintf(el_dat[j].chem_s,"Zr");
                break;
            case 41:
                sprintf(el_dat[j].chem_s,"Nb");
                break;
            case 42:
                sprintf(el_dat[j].chem_s,"Mo");
                break;
            case 43:
                sprintf(el_dat[j].chem_s,"Tc");
                break;
            case 44:
                sprintf(el_dat[j].chem_s,"Ru");
                break;
            case 45:
                sprintf(el_dat[j].chem_s,"Rh");
                break;
            case 46:
                sprintf(el_dat[j].chem_s,"Pd");
                break;
            case 47:
                sprintf(el_dat[j].chem_s,"Ag");
                break;
            case 48:
                sprintf(el_dat[j].chem_s,"Cd");
                break;
            case 49:
                sprintf(el_dat[j].chem_s,"In");
                break;
            case 50:
                sprintf(el_dat[j].chem_s,"Sn");
                break;
            case 51:
                sprintf(el_dat[j].chem_s,"Sb");
                break;
            case 52:
                sprintf(el_dat[j].chem_s,"Te");
                break;
            case 53:
                sprintf(el_dat[j].chem_s,"I ");
                break;
            case 54:
                sprintf(el_dat[j].chem_s,"Xe");
                break;
            case 55:
                sprintf(el_dat[j].chem_s,"Cs");
                break;
            case 56:
                sprintf(el_dat[j].chem_s,"Ba");
                break;
            case 57:
                sprintf(el_dat[j].chem_s,"La");
                break;
            case 58:
                sprintf(el_dat[j].chem_s,"Ce");
                break;
            case 59:
                sprintf(el_dat[j].chem_s,"Pr");
                break;
            case 60:
                sprintf(el_dat[j].chem_s,"Nd");
                break;
            case 61:
                sprintf(el_dat[j].chem_s,"Pm");
                break;
            case 62:
                sprintf(el_dat[j].chem_s,"Sm");
                break;
            case 63:
                sprintf(el_dat[j].chem_s,"Eu");
                break;
            case 64:
                sprintf(el_dat[j].chem_s,"Gd");
                break;
            case 65:
                sprintf(el_dat[j].chem_s,"Tb");
                break;
            case 66:
                sprintf(el_dat[j].chem_s,"Dy");
                break;
            case 67:
                sprintf(el_dat[j].chem_s,"Ho");
                break;
            case 68:
                sprintf(el_dat[j].chem_s,"Er");
                break;
            case 69:
                sprintf(el_dat[j].chem_s,"Tm");
                break;
            case 70:
                sprintf(el_dat[j].chem_s,"Yb");
                break;
            case 71:
                sprintf(el_dat[j].chem_s,"Lu");
                break;
            case 72:
                sprintf(el_dat[j].chem_s,"Hf");
                break;
            case 73:
                sprintf(el_dat[j].chem_s,"Ta");
                break;
            case 74:
                sprintf(el_dat[j].chem_s,"W ");
                break;
            case 75:
                sprintf(el_dat[j].chem_s,"Re");
                break;
            case 76:
                sprintf(el_dat[j].chem_s,"Os");
                break;
            case 77:
                sprintf(el_dat[j].chem_s,"Ir");
                break;
            case 78:
                sprintf(el_dat[j].chem_s,"Pt");
                break;
            case 79:
                sprintf(el_dat[j].chem_s,"Au");
                break;
            case 80:
                sprintf(el_dat[j].chem_s,"Hg");
                break;
            case 81:
                sprintf(el_dat[j].chem_s,"Tl");
                break;			
            case 82:
                sprintf(el_dat[j].chem_s,"Pb");
                break;
            case 83:
                sprintf(el_dat[j].chem_s,"Bi");
                break;
            case 84:
                sprintf(el_dat[j].chem_s,"Po");
                break;
            case 85:
                sprintf(el_dat[j].chem_s,"At");
                break;
            case 86:
                sprintf(el_dat[j].chem_s,"Rn");
                break;
            case 87:
                sprintf(el_dat[j].chem_s,"Fr");
                break;
            case 88:
                sprintf(el_dat[j].chem_s,"Ra");
                break;
            case 89:
                sprintf(el_dat[j].chem_s,"Ac");
                break;
            case 90:
                sprintf(el_dat[j].chem_s,"Th");
                break;
            case 91:
                sprintf(el_dat[j].chem_s,"Pa");
                break;
            case 92:
                sprintf(el_dat[j].chem_s,"U ");
                break;
            case 93:
                sprintf(el_dat[j].chem_s,"Np");
                break;
            case 94:
                sprintf(el_dat[j].chem_s,"Pu");
                break;
            case 95:
                sprintf(el_dat[j].chem_s,"Am");
                break;
            case 96:
                sprintf(el_dat[j].chem_s,"Cm");
                break;
            case 97:
                sprintf(el_dat[j].chem_s,"Bk");
                break;
            case 98:
                sprintf(el_dat[j].chem_s,"Cf");
                break;
            case 99:
                sprintf(el_dat[j].chem_s,"Es");
                break;
            case 100:
                sprintf(el_dat[j].chem_s,"Fm");
                break;
        }
    }
    
    

  struct ACE_spectra
  {
    char S[10];						// string character for structure index
    double EI[201];				// energy for structure index [0], flux otherwise
		double EI_int[201];		// integrated flux spectra
		double I_spec;				// integrated spectrum from 10 to 10e5 MeV
		double I_E_hi;				// integrated power law tail to CALET high energy limit
		double I_E_max;				// integrated power law tail to infinity
		double I_tot;					// total integrated spectrum with added power law tail
		double rel_ab_I;			// relative abundance due to integral of spectum
		double c_pl;					// power law constat
		double i_pl;					// power law index
  } S_max[100], S_min[100];

	int Z_index=0;
	int E_index=0;
	int long_index=0;
	int lat_index=0;
	int el_index=0;


	int i=0;
	int k=0;
	
	
	
	for (i=0; i<100; i++) {
		el_dat[i].abund_S_max=0;
		el_dat[i].abund_S_min=0;
		el_dat[i].abund_S_max_gf=0;
		el_dat[i].abund_S_min_gf=0;
		el_dat[i].abund_S_max_dgf=0;
		el_dat[i].abund_S_min_dgf=0;
		el_dat[i].abund_S_max_dgf_noic=0;
		el_dat[i].abund_S_min_dgf_noic=0;
		el_dat[i].abund_S_max_tot=0;
		el_dat[i].abund_S_min_tot=0;
		el_dat[i].abund_S_max_no_mag=0;
		el_dat[i].abund_S_min_no_mag=0;
	}
	
	
	double x_plot[250];
	double y_plot[250];
	
	int plot_count=0;

	double longitude[long_bins];

	struct Mag_Cut
	{
		double latitude;
		double cut_off[long_bins];
	} alt_450[lat_bins];

	
	double x_cont[long_bins_calc*lat_bins];
	double y_cont[long_bins_calc*lat_bins];
	double z_cont[long_bins_calc*lat_bins];
	
	double ISS_angle_rad=0.0;
  ISS_angle_rad=ISS_inc*TMath::Pi()/180.0;

	double lat_max_rad;
	double lat_max_deg;
	lat_max_rad=ISS_angle_rad;
	lat_max_deg=180.0*lat_max_rad/TMath::Pi();
	printf("lat_max: %f %f\n",lat_max_rad,lat_max_deg);

	// Format from email from Kelly Lave: /data/tiger/calet/abundance_estimate/rigidity/tables/CRIS stuff.html
  struct ACE_abund
  {
	int Z;					// charge
	double I;				// flux ([cm^2 s sr MeV/n]^-1)
	double pu; 			// percent uncertainty in flux
	double s_stat;	// squared statistical uncertainty
	double s_sys;		// squared systematic uncertainty
	double rel_a;		// relative abundance
  } A_max[24], A_min[24];

	int ACE_Z_index_smin=0;
	int ACE_Z_index_smax=0;

	int ref_index=0; // index for element that relative abundances are in reference to (usually Fe)

	double E_spec_hi=0.0;

	double time_in_orbit=365*24*60*60;

    double SS_LS=0.0; // Solar System abundances of the light secondary group 62 <= Z <= 69
    double SS_HS=0.0; // Solar System abundances of the heavy secondary group 70 <= Z <= 73
    double SS_Pt=0.0; // Solar System abundances of the Platinum group 74 <= Z <= 80
    double SS_Pb=0.0; // Solar System abundances of the Lead group 81 <= Z <= 83
    double SS_Ac=0.0; // Solar System abundances of the Actinide group 88 <= Z <= 100
    
	char header[100];
	char header_mag[100];

	// Reading in element information: Z, A, SS abudnaces..
	int el_dat_Fe_index=0;

  in_1.open(ab_input);

	el_index=0;
  while(1) {
    in_1 >> el_dat[el_index].Z >> el_dat[el_index].A >> el_dat[el_index].SS_abund  >> el_dat[el_index].HEAO_HNE >> el_dat[el_index].HEAO_C2;
    if(in_1.eof() || !in_1.good()) break;
	  el_dat[el_index].R = 1.58*pow(el_dat[el_index].A,0.281);

		el_dat[el_index].tiger = 0.0;
		if(el_dat[el_index].Z == 26) el_dat_Fe_index=el_index;
      if (el_dat[el_index].Z >= 26 && el_dat[el_index].Z <= 29) {
          el_dat[el_index].TOA_abund = el_dat[el_index].HEAO_C2/
          el_dat[el_dat_Fe_index].HEAO_C2;
      }

    if(el_dat[el_index].Z > 29 && el_dat[el_index].Z <=60) {
      if(el_dat[el_index].Z%2==0) {
				el_dat[el_index].TOA_abund = even_frac*el_dat[el_index].HEAO_HNE/
					el_dat[el_dat_Fe_index].HEAO_HNE;
       //   if(el_dat[el_index].Z==)
      }
      if(el_dat[el_index].Z%2==1) {
				el_dat[el_index].TOA_abund = odd_frac*el_dat[el_index].HEAO_HNE/
					el_dat[el_dat_Fe_index].HEAO_HNE;
      }
    }
    if(el_dat[el_index].Z == 61) el_dat[el_index].TOA_abund = 0.0;
    // Find total SS relative abundance in the light secondary group
    if(el_dat[el_index].Z > 61 && el_dat[el_index].Z <= 69) {
        SS_LS+=el_dat[el_index].SS_abund;
    }
    // Find total SS relative abundance in the heavy secondary group
    if(el_dat[el_index].Z > 69 && el_dat[el_index].Z <= 73) {
        SS_HS+=el_dat[el_index].SS_abund;
    }
    // Find total SS relative abundance in the Platinum group
    if(el_dat[el_index].Z > 73 && el_dat[el_index].Z <= 80) {
        SS_Pt+=el_dat[el_index].SS_abund;
    }
    // Find total SS relative abundance in the Lead group
    if(el_dat[el_index].Z > 80 && el_dat[el_index].Z <= 83) {
        SS_Pb+=el_dat[el_index].SS_abund;
    }
    // No particles in this range
    if(el_dat[el_index].Z > 83 && el_dat[el_index].Z <= 87) {
        el_dat[i].TOA_abund=0.0;
    }
    // Find total SS relative abundance in the Actinide group
    if(el_dat[el_index].Z > 87 && el_dat[el_index].Z <= 92) {
        SS_Ac+=el_dat[el_index].SS_abund;
    }
      
    el_index++;
  }
  in_1.close();
  in_1.clear();
    
    

 
    

	
	
	for(i=0;i<el_index;i++)
	{
		
			//el_dat[i].E_gf_thresh=E_gf_thresh_Fe*pow((el_dat[i].Z/26.0),1.71);
		el_dat[i].E_gf_thresh=E_gf_thresh_Fe;
		el_dat[i].R_cut = (sqrt(pow(el_dat[i].E_gf_thresh,2.0) + 2.0*el_dat[i].E_gf_thresh*amu_MeV)*el_dat[i].A/el_dat[i].Z)/1000;

		fprintf(stderr,"%2d E_gf_thresh_Fe: %.2f R_cut: %.2f\n",el_dat[i].Z,E_gf_thresh_Fe,el_dat[i].R_cut);

        
        // Assign TOA relative abundances in HEAO-HNE charge group ranges based on SS relative abundnaces
        
        // Weight the LS-group (Light Secondaries) 62 <= Z <= 69 by group fractional SS abundance
        if(el_dat[i].Z > 61 && el_dat[i].Z <= 69) {
            el_dat[i].TOA_abund = (el_dat[i].SS_abund/SS_LS)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        // Weight the HS-group (Heavy Secondaries) 70 <= Z <= 73 by group fractional SS abundance
        if(el_dat[i].Z > 69 && el_dat[i].Z <= 73) {
            el_dat[i].TOA_abund = (el_dat[i].SS_abund/SS_HS)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        // Weight the Pt-group 74 <= Z <= 80 by group fractional SS abundance
        // for W and Re
        if(el_dat[i].Z > 73 && el_dat[i].Z <= 75) {
            el_dat[i].TOA_abund = (el_dat[i].SS_abund/SS_Pt)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        // Os is 0.38 times the "Pt-group"
        // 75 <= Z <=79 in Westphal and 74 <= Z <= 80 in HEAO-HNE
        // Westphal et al., Nature, 396, 50, 1998
        // Multiply HEAO_HNE Pt-group by weighted Westphal fractional
        // abundances relative to the 76 <= Z <= 80 group
        // 0.38 + 0.21 + 0.48 + 0.05 + 0.12 = 1.24
        // (0.38 / 1.24) * HEAO-HNE Pt-group / HEAO-HNE Fe
        if(el_dat[i].Z == 76) {
            el_dat[i].TOA_abund = (0.38/1.24)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        // Ir is 0.21 times the "Pt-group"
        // 75 <= Z <=79 in Westphal and 74 <= Z <= 80 in HEAO-HNE
        // Westphal et al., Nature, 396, 50, 1998
        // Multiply HEAO_HNE Pt-group by weighted Westphal fractional
        // abundances relative to the 76 <= Z <= 80 group
        // 0.38 + 0.21 + 0.48 + 0.05 + 0.12 = 1.24
        // (0.21 / 1.24) * HEAO-HNE Pt-group / HEAO-HNE Fe
        if(el_dat[i].Z == 77) {
            el_dat[i].TOA_abund = (0.21/1.24)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        // Pt is 0.48 times the "Pt-group"
        // 75 <= Z <=79 in Westphal and 74 <= Z <= 80 in HEAO-HNE
        // Westphal et al., Nature, 396, 50, 1998
        // Multiply HEAO_HNE Pt-group by weighted Westphal fractional
        // abundances relative to the 76 <= Z <= 80 group
        // 0.38 + 0.21 + 0.48 + 0.05 + 0.12 = 1.24
        // (0.48 / 1.24) * HEAO-HNE Pt-group / HEAO-HNE Fe
        if(el_dat[i].Z == 78) {
            el_dat[i].TOA_abund = (0.48/1.24)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        // Au is 0.05 times the "Pt-group"
        // 75 <= Z <=79 in Westphal and 74 <= Z <= 80 in HEAO-HNE
        // Westphal et al., Nature, 396, 50, 1998
        // Multiply HEAO_HNE Pt-group by weighted Westphal fractional
        // abundances relative to the 76 <= Z <= 80 group
        // 0.38 + 0.21 + 0.48 + 0.05 + 0.12 = 1.24
        // (0.05 / 1.24) * HEAO-HNE Pt-group / HEAO-HNE Fe
        if(el_dat[i].Z == 79) {
            el_dat[i].TOA_abund = (0.05/1.24)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        // Hg is 0.12 times the "Pt-group"
        // 75 <= Z <=79 in Westphal and 74 <= Z <= 80 in HEAO-HNE
        // Westphal et al., Nature, 396, 50, 1998
        // Multiply HEAO_HNE Pt-group by weighted Westphal fractional
        // abundances relative to the 76 <= Z <= 80 group
        // 0.38 + 0.21 + 0.48 + 0.05 + 0.12 = 1.24
        // (0.12 / 1.24) * HEAO-HNE Pt-group / HEAO-HNE Fe
        if(el_dat[i].Z == 80) {
            el_dat[i].TOA_abund = (0.12/1.24)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        
        // Weight the Pb-group 81 <= Z <= 83 by group fractional SS abundance
        // for Tl and Bi
        if(el_dat[i].Z > 80 && el_dat[i].Z <= 83 && el_dat[i].Z != 82) {
            el_dat[i].TOA_abund = (el_dat[i].SS_abund/SS_Pb)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        // Pb is 0.29 times the "Pt-group"
        // 75 <= Z <=79 in Westphal and 74 <= Z <= 80 in HEAO-HNE
        // Westphal et al., Nature, 396, 50, 1998
        // Multiply HEAO_HNE Pt-group by weighted Westphal fractional
        // abundances relative to the 76 <= Z <= 80 group
        // 0.38 + 0.21 + 0.48 + 0.05 + 0.12 = 1.24
        // (0.29 / 1.24) * HEAO-HNE Pt-group / HEAO-HNE Fe
        if(el_dat[i].Z == 82) {
            el_dat[i].TOA_abund = 2.34E-01*el_dat[i-2].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        // Weight the Actinide-group 88 <= Z <= 92 by group fractional SS abundance
        if(el_dat[i].Z > 87 && el_dat[i].Z <= 92) {
            el_dat[i].TOA_abund = (el_dat[i].SS_abund/SS_Ac)*el_dat[i].HEAO_HNE/
            el_dat[el_dat_Fe_index].HEAO_HNE;
        }
        
	}
	
	// Reading in TIGER TOA values

	in_1.open(supertiger_toa);

	int tiger_index=25;
	int Z_temp=0;
	while(1) {
		in_1 >> Z_temp >> el_dat[tiger_index].tiger;
    if(in_1.eof() || !in_1.good()) break;
		if(el_dat[tiger_index].Z != Z_temp)
			{
			fprintf(stderr,"Data synchronization error: %d != %d\n",
							el_dat[tiger_index].Z,Z_temp);
			}
		fprintf(stderr,"%2d %6.2f %.3e %.3e %.3e %.3e\n",el_dat[tiger_index].Z,
						el_dat[tiger_index].A,
						el_dat[tiger_index].SS_abund/el_dat[el_dat_Fe_index].SS_abund,
						el_dat[tiger_index].HEAO_C2/el_dat[el_dat_Fe_index].HEAO_C2,
						el_dat[tiger_index].tiger,el_dat[tiger_index].TOA_abund);
		if(Z_temp>=30)
			{
			el_dat[tiger_index].TOA_abund = el_dat[tiger_index].tiger;
			}
    tiger_index++;
  }
  in_1.close();
  in_1.clear();
	fprintf(stderr,"tiger_index: %d el_index: %d\n",tiger_index,el_index);
	
    int i_ind=0;
    
    fprintf(stdout,"GCR TOA\n");
    for (i_ind=0; i_ind<el_index; i_ind++) {
        fprintf(stdout,"%d %.6e\n",el_dat[i_ind].Z,el_dat[i_ind].TOA_abund);
    }
    
    

    
    

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

			for(i=0;i<el_index;i++)
			{
				el_dat[i].E_min_mag[long_index][lat_index] = -amu_MeV + sqrt(pow(amu_MeV,2.0) +\
						 (pow(alt_450[lat_index].cut_off[long_index]*1e3,2.0)*pow(el_dat[i].Z,2.0))/pow(el_dat[i].A,2.0));

				/*
 			 	fprintf(stderr,"%2d el_dat[%2d].E_min_mag[%2d][%2d]: %.4e Lat: %2.0f Long: %3.0f\n",el_dat[i].Z,i,long_index,lat_index,\
				el_dat[i].E_min_mag[long_index][lat_index],alt_450[lat_index].latitude,longitude[long_index]);
				*/
                 
			// KE/A = -amu_MeV * c_vac^2 + sqrt(amu_MeV^2 * c_vac^4 + (R^2 * Z^2 *c_vac^2)/A^2)
			}
		}


		lat_index++;		

  }


	in_1.close();
	in_1.clear();
	
	double northern_hemisphere_cut[long_bins_calc];
	double southern_hemisphere_cut[long_bins_calc];

	int c_count=0;
	
	for (i=0; i<long_bins_calc; i++) {
		for (j=0; j<lat_bins; j++) {
			x_cont[c_count]=longitude[i];
			if (longitude[i]>= 180) {
				x_cont[c_count]=longitude[i]-360;
			}
			y_cont[c_count]=alt_450[j].latitude;
			z_cont[c_count]=alt_450[j].cut_off[i];
			c_count++;

		}
	}


	// Constructing latitude bin ranges

	double lat_bin[lat_index+1];

	lat_bin[0]=alt_450[0].latitude; // setting the upper bound for the lattitude range to 90
	lat_bin[lat_index]=alt_450[lat_index-1].latitude; // setting the lower bound for the lattitude range to -90

	for(i=0;i<lat_index-1;i++)
	{
		lat_bin[i+1]=alt_450[i].latitude-(alt_450[i].latitude-alt_450[i+1].latitude)/2;
	}

	/* // Checking the latitude bin ranges
	for(i=0;i<=lat_index;i++)
	{
		printf("%6.2f\n",lat_bin[i]);
	}
	*/

	double dt[lat_index+1];
	double t_tot=0.0;

	for(i=0;i<lat_index;i++)
	{
		if(lat_bin[i] < lat_max_deg && -lat_max_deg < lat_bin[i+1])
		{
			printf("%6.2f %6.2f %f\n",lat_bin[i],lat_bin[i+1],lat_max_deg);
			dt[i]=(asin(sin(lat_bin[i]*TMath::Pi()/180.0)/sin(lat_max_rad)) - asin(sin(lat_bin[i+1]*TMath::Pi()/180.0)/sin(lat_max_rad)))/TMath::Pi();
			printf("%f\n\n",dt[i]);
			t_tot+=dt[i];
		}
		else if(lat_bin[i] > lat_max_deg && lat_bin[i+1] < lat_max_deg)
		{
			printf("%6.2f %6.2f %f !!\n",lat_bin[i],lat_bin[i+1],lat_max_deg);
			dt[i]=(asin(sin(lat_max_deg*TMath::Pi()/180.0)/sin(lat_max_rad)) - asin(sin(lat_bin[i+1]*TMath::Pi()/180.0)/sin(lat_max_rad)))/TMath::Pi();
			printf("%f\n\n",dt[i]);
			t_tot+=dt[i];
		}
		else if(lat_bin[i] > -lat_max_deg && lat_bin[i+1] < -lat_max_deg)
		{
			printf("%6.2f %6.2f %f !\n",lat_bin[i],lat_bin[i+1],lat_max_deg);
			dt[i]=(asin(sin(lat_bin[i]*TMath::Pi()/180.0)/sin(lat_max_rad)) - asin(sin(lat_max_deg*TMath::Pi()/180.0)/sin(-lat_max_rad)))/TMath::Pi();
			printf("%f\n\n",dt[i]);
			t_tot+=dt[i];
		}
		else
		{
			printf("%6.2f %6.2f\n",lat_bin[i],lat_bin[i+1]);
			dt[i]=0;
			t_tot+=dt[i];
		}
	}

	printf("t_tot: %f\n",t_tot);






in_1.open(spec_smax);
  
//    in_1.open("test");
	in_1.getline(header,100);
	printf("%s\n",header);
//n_1 >> header;

	E_index=0;

	for(Z_index=0;Z_index<33;Z_index++)
		{
	  	in_1 >> S_max[Z_index].S;
	  	// printf("%s\n",S_max[Z_index].S);
		}
		
  while(1) {
  
  	for(Z_index=0;Z_index<33;Z_index++)
		{
			in_1 >> S_max[Z_index].EI[E_index];
			// printf("S_max[%d][%d] = %e\n",Z_index,E_index,S_max[Z_index].EI[E_index]);
		}

    if(in_1.eof() || !in_1.good()) break;
		E_index++;		
  }
  in_1.close();
  in_1.clear();
  
	E_spec_hi=S_max[0].EI[E_index-1];
	fprintf(stderr,"E_spec_hi: %e\n",E_spec_hi);
  
  // Solar Max Plot Label
  
  double label_x[2];
  double label_y[2];
  label_x[0] = 1e1;
  label_x[1] = 1e5;
  label_y[0] = 1e-17;
  label_y[1] = 1e-6;

 
  TGraph *gr_label_S_max_spectra = new TGraph(2,label_x,label_y);

  gr_label_S_max_spectra->SetMarkerStyle(1);
  gr_label_S_max_spectra->SetMarkerSize(1);
  gr_label_S_max_spectra->SetTitle("ACE-CRIS Solar Maximum Spectra");
  gr_label_S_max_spectra->GetXaxis()->SetTitle("MeV/nuc");
  gr_label_S_max_spectra->GetXaxis()->CenterTitle();
  gr_label_S_max_spectra->GetYaxis()->SetTitle("1/(cm^{2} s sr MeV/nuc)");
  gr_label_S_max_spectra->GetYaxis()->CenterTitle();
  gr_label_S_max_spectra->GetXaxis()->SetTitleOffset(1.2);
  gr_label_S_max_spectra->GetYaxis()->SetTitleOffset(1.4);
  gr_label_S_max_spectra->GetYaxis()->SetLabelFont(nicefont);
  gr_label_S_max_spectra->GetYaxis()->SetTitleFont(nicefont);
  gr_label_S_max_spectra->GetXaxis()->SetLabelFont(nicefont);
  gr_label_S_max_spectra->GetXaxis()->SetTitleFont(nicefont);
  gr_label_S_max_spectra->GetXaxis()->SetRangeUser(1e1,1e5);
  gr_label_S_max_spectra->GetYaxis()->SetRangeUser(1e-17,1e-6);

  // Solar Max Integral Plot Label
  
  double label_max_int_x[2];
  double label_max_int_y[2];
  label_max_int_x[0] = 1e1;
  label_max_int_x[1] = 1e5;
  label_max_int_y[0] = 1e-13;
  label_max_int_y[1] = 1e-3;

 
  TGraph *gr_label_S_max_int_spectra = new TGraph(2,label_max_int_x,label_max_int_y);

  gr_label_S_max_int_spectra->SetMarkerStyle(1);
  gr_label_S_max_int_spectra->SetMarkerSize(1);
  gr_label_S_max_int_spectra->SetTitle("ACE-CRIS Solar Maximum Integral Spectra");
  gr_label_S_max_int_spectra->GetXaxis()->SetTitle("MeV/nuc");
  gr_label_S_max_int_spectra->GetXaxis()->CenterTitle();
  gr_label_S_max_int_spectra->GetYaxis()->SetTitle("1/(cm^{2} s sr)");
  gr_label_S_max_int_spectra->GetYaxis()->CenterTitle();
  gr_label_S_max_int_spectra->GetXaxis()->SetTitleOffset(1.2);
  gr_label_S_max_int_spectra->GetYaxis()->SetTitleOffset(1.4);
  gr_label_S_max_int_spectra->GetYaxis()->SetLabelFont(nicefont);
  gr_label_S_max_int_spectra->GetYaxis()->SetTitleFont(nicefont);
  gr_label_S_max_int_spectra->GetXaxis()->SetLabelFont(nicefont);
  gr_label_S_max_int_spectra->GetXaxis()->SetTitleFont(nicefont);
  gr_label_S_max_int_spectra->GetXaxis()->SetRangeUser(1e1,1e5);
  gr_label_S_max_int_spectra->GetYaxis()->SetRangeUser(1e-13,1e-3);

	// Declaring the linear fit to be used to find the integral of the high energy tail of the cosmic ray spectra
	TF1 *fp1 = new TF1("fp1","pol1");

	double parp1[2];

	// epsilon variable required for the Integralfunc
	double epsilon=1e-6;

  // Solar Max Plot Spectra

  // Z = 5, Boron 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[5].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_5 = new TGraph(plot_count,x_plot,y_plot);
  S_max[5].I_spec = Integralfunc(gr_S_max_Z_5,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[5].I_spec: %e\n",S_max[5].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[5].EI_int[i]=Integralfunc(gr_S_max_Z_5,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_max[5].EI_int[%d]: %e\n",i,S_max[5].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[5].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_5 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

//	fprintf(stderr,"E_index: %d\n",E_index);

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[5].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_5_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_5_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[5].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[5].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[5].I_tot=S_max[5].I_spec + S_max[5].I_E_max;
	S_max[5].c_pl=pow(10.0,parp1[0]);
	S_max[5].i_pl=parp1[1];

//	fprintf(stderr,"%e\n%f\n",pow(10.0,parp1[0]),parp1[1]);
//	fprintf(stderr,"S_max[5].I_E_max: %e\n",S_max[5].I_E_max);
//	fprintf(stderr,"S_max[5].I_E_hi: %e\n",S_max[5].I_E_hi);
//	fprintf(stderr,"S_max[5].I_tot: %e\n",S_max[5].I_tot);
  
  TPaveText *p_Z_5_S_max = new TPaveText(0.826,0.340,0.866,0.380,"NDC");
  p_Z_5_S_max->SetBorderSize(0);
  p_Z_5_S_max->SetFillColor(kWhite);
  TText *t_Z_5_S_max = p_Z_5_S_max->AddText(S_max[5].S);
  t_Z_5_S_max->SetTextFont(nicefont);
  t_Z_5_S_max->SetTextColor(1);
  
  
  // Z = 6, Carbon 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[6].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_6 = new TGraph(plot_count,x_plot,y_plot);
  S_max[6].I_spec = Integralfunc(gr_S_max_Z_6,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[6].I_spec: %e\n",S_max[5].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[6].EI_int[i]=Integralfunc(gr_S_max_Z_6,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[6].EI_int[%d]: %e\n",i,S_max[6].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[6].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_6 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum


  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[6].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_6_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_6_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[6].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[6].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[6].I_tot=S_max[6].I_spec + S_max[6].I_E_max;
	S_max[6].c_pl=pow(10.0,parp1[0]);
	S_max[6].i_pl=parp1[1];
  
  TPaveText *p_Z_6_S_max = new TPaveText(0.826,0.340,0.866,0.380,"NDC");
  p_Z_6_S_max->SetBorderSize(0);
  p_Z_6_S_max->SetFillColor(kWhite);
  TText *t_Z_6_S_max = p_Z_6_S_max->AddText(S_max[6].S);
  t_Z_6_S_max->SetTextFont(nicefont);
  t_Z_6_S_max->SetTextColor(1);
   
  // Z = 7, Nitrogen 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[7].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_7 = new TGraph(plot_count,x_plot,y_plot);
  S_max[7].I_spec = Integralfunc(gr_S_max_Z_7,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[7].I_spec: %e\n",S_max[7].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[7].EI_int[i]=Integralfunc(gr_S_max_Z_7,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[7].EI_int[%d]: %e\n",i,S_max[7].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[7].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_7 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum


  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[7].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_7_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_7_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[7].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[7].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[7].I_tot=S_max[7].I_spec + S_max[7].I_E_max;
	S_max[7].c_pl=pow(10.0,parp1[0]);
	S_max[7].i_pl=parp1[1];
  
  TPaveText *p_Z_7_S_max = new TPaveText(0.827,0.340,0.877,0.380,"NDC");
  p_Z_7_S_max->SetBorderSize(0);
  p_Z_7_S_max->SetFillColor(kWhite);
  TText *t_Z_7_S_max = p_Z_7_S_max->AddText(S_max[7].S);
  t_Z_7_S_max->SetTextFont(nicefont);
  t_Z_7_S_max->SetTextColor(1);
  
    // Z = 8, Oxygen 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[8].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_8 = new TGraph(plot_count,x_plot,y_plot);
  S_max[8].I_spec = Integralfunc(gr_S_max_Z_8,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[8].I_spec: %e\n",S_max[8].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[8].EI_int[i]=Integralfunc(gr_S_max_Z_8,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[8].EI_int[%d]: %e\n",i,S_max[8].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[8].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_8 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[8].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_8_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_8_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[8].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[8].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[8].I_tot=S_max[8].I_spec + S_max[8].I_E_max;
	S_max[8].c_pl=pow(10.0,parp1[0]);
	S_max[8].i_pl=parp1[1];

    // Z = 9, Flourine
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[9].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_9 = new TGraph(plot_count,x_plot,y_plot);
  S_max[9].I_spec = Integralfunc(gr_S_max_Z_9,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[9].I_spec: %e\n",S_max[9].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[9].EI_int[i]=Integralfunc(gr_S_max_Z_9,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[9].EI_int[%d]: %e\n",i,S_max[9].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[9].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_9 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[9].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_9_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_9_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[9].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[9].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[9].I_tot=S_max[9].I_spec + S_max[9].I_E_max;
	S_max[9].c_pl=pow(10.0,parp1[0]);
	S_max[9].i_pl=parp1[1];
  
    // Z = 10, Neon 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[10].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_10 = new TGraph(plot_count,x_plot,y_plot);  
  S_max[10].I_spec = Integralfunc(gr_S_max_Z_10,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[10].I_spec: %e\n",S_max[10].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[10].EI_int[i]=Integralfunc(gr_S_max_Z_10,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[10].EI_int[%d]: %e\n",i,S_max[10].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[10].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_10 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[10].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_10_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_10_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[10].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[10].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[10].I_tot=S_max[10].I_spec + S_max[10].I_E_max;
	S_max[10].c_pl=pow(10.0,parp1[0]);
	S_max[10].i_pl=parp1[1];
        
    // Z = 11, Sodium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[11].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_11 = new TGraph(plot_count,x_plot,y_plot);
  S_max[11].I_spec = Integralfunc(gr_S_max_Z_11,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[11].I_spec: %e\n",S_max[11].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[11].EI_int[i]=Integralfunc(gr_S_max_Z_11,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[11].EI_int[%d]: %e\n",i,S_max[11].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[11].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_11 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[11].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_11_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_11_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[11].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[11].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[11].I_tot=S_max[11].I_spec + S_max[11].I_E_max;
	S_max[11].c_pl=pow(10.0,parp1[0]);
	S_max[11].i_pl=parp1[1];
  
    // Z = 12, Magnesium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[12].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_12 = new TGraph(plot_count,x_plot,y_plot);  
  S_max[12].I_spec = Integralfunc(gr_S_max_Z_12,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[12].I_spec: %e\n",S_max[12].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[12].EI_int[i]=Integralfunc(gr_S_max_Z_12,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[12].EI_int[%d]: %e\n",i,S_max[12].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[12].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_12 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[12].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_12_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_12_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[12].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[12].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[12].I_tot=S_max[12].I_spec + S_max[12].I_E_max;
	S_max[12].c_pl=pow(10.0,parp1[0]);
	S_max[12].i_pl=parp1[1];
  
  
    // Z = 13, Aluminum 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[13].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_13 = new TGraph(plot_count,x_plot,y_plot);  
  S_max[13].I_spec = Integralfunc(gr_S_max_Z_13,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[13].I_spec: %e\n",S_max[13].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[13].EI_int[i]=Integralfunc(gr_S_max_Z_13,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[13].EI_int[%d]: %e\n",i,S_max[13].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[13].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_13 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[13].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_13_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_13_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[13].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[13].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[13].I_tot=S_max[13].I_spec + S_max[13].I_E_max;
	S_max[13].c_pl=pow(10.0,parp1[0]);
	S_max[13].i_pl=parp1[1];
    
    // Z = 14, Silicon 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[14].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_14 = new TGraph(plot_count,x_plot,y_plot);
  S_max[14].I_spec = Integralfunc(gr_S_max_Z_14,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[14].I_spec: %e\n",S_max[14].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[14].EI_int[i]=Integralfunc(gr_S_max_Z_14,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[14].EI_int[%d]: %e\n",i,S_max[14].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[14].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_14 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[14].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_14_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_14_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[14].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[14].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[14].I_tot=S_max[14].I_spec + S_max[14].I_E_max;
	S_max[14].c_pl=pow(10.0,parp1[0]);
	S_max[14].i_pl=parp1[1];
     
    // Z = 15, Phosphorous
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[15].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_15 = new TGraph(plot_count,x_plot,y_plot);
  S_max[15].I_spec = Integralfunc(gr_S_max_Z_15,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[15].I_spec: %e\n",S_max[15].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[15].EI_int[i]=Integralfunc(gr_S_max_Z_15,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[15].EI_int[%d]: %e\n",i,S_max[15].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[15].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_15 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[15].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_15_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_15_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[15].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[15].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[15].I_tot=S_max[15].I_spec + S_max[15].I_E_max;
	S_max[15].c_pl=pow(10.0,parp1[0]);
	S_max[15].i_pl=parp1[1];
 
    // Z = 16, Sulfur
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[16].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_16 = new TGraph(plot_count,x_plot,y_plot);
  S_max[16].I_spec = Integralfunc(gr_S_max_Z_16,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[16].I_spec: %e\n",S_max[16].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[16].EI_int[i]=Integralfunc(gr_S_max_Z_16,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[16].EI_int[%d]: %e\n",i,S_max[16].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[16].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_16 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[16].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_16_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_16_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[16].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[16].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[16].I_tot=S_max[16].I_spec + S_max[16].I_E_max;
	S_max[16].c_pl=pow(10.0,parp1[0]);
	S_max[16].i_pl=parp1[1];
   
    // Z = 17, Chlorine 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[17].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_17 = new TGraph(plot_count,x_plot,y_plot);
  S_max[17].I_spec = Integralfunc(gr_S_max_Z_17,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[17].I_spec: %e\n",S_max[17].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[17].EI_int[i]=Integralfunc(gr_S_max_Z_17,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[17].EI_int[%d]: %e\n",i,S_max[17].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[17].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_17 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[17].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_17_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_17_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[17].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[17].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[17].I_tot=S_max[17].I_spec + S_max[17].I_E_max;
	S_max[17].c_pl=pow(10.0,parp1[0]);
	S_max[17].i_pl=parp1[1];
     
    // Z = 18, Argon 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[18].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_18 = new TGraph(plot_count,x_plot,y_plot);
  S_max[18].I_spec = Integralfunc(gr_S_max_Z_18,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[18].I_spec: %e\n",S_max[18].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[18].EI_int[i]=Integralfunc(gr_S_max_Z_18,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[18].EI_int[%d]: %e\n",i,S_max[18].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[18].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_18 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[18].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_18_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_18_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[18].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[18].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[18].I_tot=S_max[18].I_spec + S_max[18].I_E_max;
	S_max[18].c_pl=pow(10.0,parp1[0]);
	S_max[18].i_pl=parp1[1];
    
    // Z = 19, Potassium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[19].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_19 = new TGraph(plot_count,x_plot,y_plot);
  S_max[19].I_spec = Integralfunc(gr_S_max_Z_19,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[19].I_spec: %e\n",S_max[19].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[19].EI_int[i]=Integralfunc(gr_S_max_Z_19,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[19].EI_int[%d]: %e\n",i,S_max[19].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[19].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_19 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[19].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_19_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_19_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[19].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[19].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[19].I_tot=S_max[19].I_spec + S_max[19].I_E_max;
	S_max[19].c_pl=pow(10.0,parp1[0]);
	S_max[19].i_pl=parp1[1];
     
    // Z = 20, Calcium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[20].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_20 = new TGraph(plot_count,x_plot,y_plot);
  S_max[20].I_spec = Integralfunc(gr_S_max_Z_20,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[20].I_spec: %e\n",S_max[20].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[20].EI_int[i]=Integralfunc(gr_S_max_Z_20,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[20].EI_int[%d]: %e\n",i,S_max[20].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[20].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_20 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[20].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_20_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_20_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[20].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[20].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[20].I_tot=S_max[20].I_spec + S_max[20].I_E_max;
	S_max[20].c_pl=pow(10.0,parp1[0]);
	S_max[20].i_pl=parp1[1];
    
    // Z = 21, Scandium
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[21].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_21 = new TGraph(plot_count,x_plot,y_plot);
  S_max[21].I_spec = Integralfunc(gr_S_max_Z_21,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[21].I_spec: %e\n",S_max[21].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[21].EI_int[i]=Integralfunc(gr_S_max_Z_21,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[21].EI_int[%d]: %e\n",i,S_max[21].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[21].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_21 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum


  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[21].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_21_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_21_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[21].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[21].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[21].I_tot=S_max[21].I_spec + S_max[21].I_E_max;
	S_max[21].c_pl=pow(10.0,parp1[0]);
	S_max[21].i_pl=parp1[1];
    
    // Z = 22, Titanium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[22].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_22 = new TGraph(plot_count,x_plot,y_plot);
  S_max[22].I_spec = Integralfunc(gr_S_max_Z_22,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[22].I_spec: %e\n",S_max[22].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[22].EI_int[i]=Integralfunc(gr_S_max_Z_22,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[22].EI_int[%d]: %e\n",i,S_max[22].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[22].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_22 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[22].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_22_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_22_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[22].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[22].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[22].I_tot=S_max[22].I_spec + S_max[22].I_E_max;
	S_max[22].c_pl=pow(10.0,parp1[0]);
	S_max[22].i_pl=parp1[1];
    
    // Z = 23, Vanadium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[23].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_23 = new TGraph(plot_count,x_plot,y_plot);
  S_max[23].I_spec = Integralfunc(gr_S_max_Z_23,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[23].I_spec: %e\n",S_max[23].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[23].EI_int[i]=Integralfunc(gr_S_max_Z_23,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[23].EI_int[%d]: %e\n",i,S_max[23].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[23].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_23 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[23].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_23_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_23_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[23].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[23].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[23].I_tot=S_max[23].I_spec + S_max[23].I_E_max;
	S_max[23].c_pl=pow(10.0,parp1[0]);
	S_max[23].i_pl=parp1[1];
    
    // Z = 24, Chromium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[24].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_24 = new TGraph(plot_count,x_plot,y_plot);
  S_max[24].I_spec = Integralfunc(gr_S_max_Z_24,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[24].I_spec: %e\n",S_max[24].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[24].EI_int[i]=Integralfunc(gr_S_max_Z_24,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[24].EI_int[%d]: %e\n",i,S_max[24].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[24].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_24 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[24].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_24_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_24_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[24].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[24].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[24].I_tot=S_max[24].I_spec + S_max[24].I_E_max;
	S_max[24].c_pl=pow(10.0,parp1[0]);
	S_max[24].i_pl=parp1[1];
   
    // Z = 25, Manganese 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[25].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_25 = new TGraph(plot_count,x_plot,y_plot);  
  S_max[25].I_spec = Integralfunc(gr_S_max_Z_25,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[25].I_spec: %e\n",S_max[25].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[25].EI_int[i]=Integralfunc(gr_S_max_Z_25,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[25].EI_int[%d]: %e\n",i,S_max[25].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[25].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_25 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[25].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_25_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_25_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[25].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[25].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[25].I_tot=S_max[25].I_spec + S_max[25].I_E_max;
	S_max[25].c_pl=pow(10.0,parp1[0]);
	S_max[25].i_pl=parp1[1];
   
    // Z = 26, Iron 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[26].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_26 = new TGraph(plot_count,x_plot,y_plot);    
  S_max[26].I_spec = Integralfunc(gr_S_max_Z_26,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[26].I_spec: %e\n",S_max[26].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[26].EI_int[i]=Integralfunc(gr_S_max_Z_26,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[26].EI_int[%d]: %e\n",i,S_max[26].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[26].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_26 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[26].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_26_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_26_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[26].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[26].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[26].I_tot=S_max[26].I_spec + S_max[26].I_E_max;
	S_max[26].c_pl=pow(10.0,parp1[0]);
	S_max[26].i_pl=parp1[1];
 
    // Z = 27, Cobalt 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[27].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_27 = new TGraph(plot_count,x_plot,y_plot);    
  S_max[27].I_spec = Integralfunc(gr_S_max_Z_27,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[27].I_spec: %e\n",S_max[27].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[27].EI_int[i]=Integralfunc(gr_S_max_Z_27,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[27].EI_int[%d]: %e\n",i,S_max[27].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[27].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_27 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[27].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_27_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_27_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[27].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[27].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[27].I_tot=S_max[27].I_spec + S_max[27].I_E_max;
	S_max[27].c_pl=pow(10.0,parp1[0]);
	S_max[27].i_pl=parp1[1];
   
    // Z = 28, Nickel 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[28].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_28 = new TGraph(plot_count,x_plot,y_plot);    
  S_max[28].I_spec = Integralfunc(gr_S_max_Z_28,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[28].I_spec: %e\n",S_max[28].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[28].EI_int[i]=Integralfunc(gr_S_max_Z_28,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[28].EI_int[%d]: %e\n",i,S_max[28].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[28].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_28 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[28].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_28_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_28_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[28].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[28].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[28].I_tot=S_max[28].I_spec + S_max[28].I_E_max;
	S_max[28].c_pl=pow(10.0,parp1[0]);
	S_max[28].i_pl=parp1[1];
    
    // Z = 29, Copper 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[29].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_29 = new TGraph(plot_count,x_plot,y_plot);    
  S_max[29].I_spec = Integralfunc(gr_S_max_Z_29,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[29].I_spec: %e\n",S_max[29].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[29].EI_int[i]=Integralfunc(gr_S_max_Z_29,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[29].EI_int[%d]: %e\n",i,S_max[29].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[29].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_29 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[29].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_29_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_29_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[29].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[29].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[29].I_tot=S_max[29].I_spec + S_max[29].I_E_max;
	S_max[29].c_pl=pow(10.0,parp1[0]);
	S_max[29].i_pl=parp1[1];
    
    // Z = 30, Zink 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[30].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_30 = new TGraph(plot_count,x_plot,y_plot);    
  S_max[30].I_spec = Integralfunc(gr_S_max_Z_30,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[30].I_spec: %e\n",S_max[30].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[30].EI_int[i]=Integralfunc(gr_S_max_Z_30,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[30].EI_int[%d]: %e\n",i,S_max[30].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[30].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_30 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum	

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[30].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_30_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_30_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[30].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[30].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[30].I_tot=S_max[30].I_spec + S_max[30].I_E_max;
	S_max[30].c_pl=pow(10.0,parp1[0]);
	S_max[30].i_pl=parp1[1];
  
    // Z = 31, Gallium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[31].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_31 = new TGraph(plot_count,x_plot,y_plot);      
  S_max[31].I_spec = Integralfunc(gr_S_max_Z_31,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[31].I_spec: %e\n",S_max[31].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[31].EI_int[i]=Integralfunc(gr_S_max_Z_31,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[31].EI_int[%d]: %e\n",i,S_max[31].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[31].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_31 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum	
 
  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[31].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_31_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_31_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[31].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[31].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[31].I_tot=S_max[31].I_spec + S_max[31].I_E_max;
	S_max[31].c_pl=pow(10.0,parp1[0]);
	S_max[31].i_pl=parp1[1];
    
    // Z = 32, Germanium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[32].EI[i];
		plot_count++;
		}
  TGraph *gr_S_max_Z_32 = new TGraph(plot_count,x_plot,y_plot);    
  S_max[32].I_spec = Integralfunc(gr_S_max_Z_32,S_max[0].EI[0],S_max[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_max[32].I_spec: %e\n",S_max[32].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_max[32].EI_int[i]=Integralfunc(gr_S_max_Z_32,S_max[0].EI[i],S_max[0].EI[E_index-1],"5",epsilon); 
		// fprintf(stderr,"S_max[32].EI_int[%d]: %e\n",i,S_max[32].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_max[0].EI[i];
    y_plot[i]=S_max[32].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_max_int_Z_32 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum	
 
  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_max[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_max[32].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_max_Z_32_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_max_Z_32_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_max[32].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_max[32].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_max[32].I_tot=S_max[32].I_spec + S_max[32].I_E_max;
	S_max[32].c_pl=pow(10.0,parp1[0]);
	S_max[32].i_pl=parp1[1];
 

  
  TCanvas *C_S_max_spectra = new TCanvas("C_S_max_spectra","ACE-CRIS S_Max Spectra",0,0,800,800);
  
  C_S_max_spectra->SetLeftMargin(0.124365);
  C_S_max_spectra->SetRightMargin(0.0748731);
  C_S_max_spectra->SetBottomMargin(0.117801);
  C_S_max_spectra->SetTopMargin(0.0824607);

  C_S_max_spectra->SetGrid();
  
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  
  gr_label_S_max_spectra->Draw("AP");

  gr_S_max_Z_5->Draw("L");
  gr_S_max_Z_6->Draw("L");
  gr_S_max_Z_7->Draw("L");
  gr_S_max_Z_8->Draw("L");
  gr_S_max_Z_9->Draw("L");
  gr_S_max_Z_10->Draw("L");
  gr_S_max_Z_11->Draw("L");
  gr_S_max_Z_12->Draw("L");
  gr_S_max_Z_13->Draw("L");
  gr_S_max_Z_14->Draw("L");
  gr_S_max_Z_15->Draw("L");
  gr_S_max_Z_16->Draw("L");
  gr_S_max_Z_17->Draw("L");
  gr_S_max_Z_18->Draw("L");
  gr_S_max_Z_19->Draw("L");
  gr_S_max_Z_20->Draw("L");
  gr_S_max_Z_21->Draw("L");
  gr_S_max_Z_22->Draw("L");
  gr_S_max_Z_23->Draw("L");
  gr_S_max_Z_24->Draw("L");
  gr_S_max_Z_25->Draw("L");
  gr_S_max_Z_26->Draw("L");
  gr_S_max_Z_27->Draw("L");
  gr_S_max_Z_28->Draw("L");
  gr_S_max_Z_29->Draw("L");
  gr_S_max_Z_30->Draw("L");
  gr_S_max_Z_31->Draw("L");
  gr_S_max_Z_32->Draw("L");


  TCanvas *C_S_max_int_spectra = new TCanvas("C_S_max_int_spectra","ACE-CRIS S_Max Spectra",0,0,800,800);
  
  C_S_max_int_spectra->SetLeftMargin(0.124365);
  C_S_max_int_spectra->SetRightMargin(0.0748731);
  C_S_max_int_spectra->SetBottomMargin(0.117801);
  C_S_max_int_spectra->SetTopMargin(0.0824607);

  C_S_max_int_spectra->SetGrid();
  
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  
  gr_label_S_max_int_spectra->Draw("AP");

  gr_S_max_int_Z_5->Draw("L");
  gr_S_max_int_Z_6->Draw("L");
  gr_S_max_int_Z_7->Draw("L");
  gr_S_max_int_Z_8->Draw("L");
  gr_S_max_int_Z_9->Draw("L");
  gr_S_max_int_Z_10->Draw("L");
  gr_S_max_int_Z_11->Draw("L");
  gr_S_max_int_Z_12->Draw("L");
  gr_S_max_int_Z_13->Draw("L");
  gr_S_max_int_Z_14->Draw("L");
  gr_S_max_int_Z_15->Draw("L");
  gr_S_max_int_Z_16->Draw("L");
  gr_S_max_int_Z_17->Draw("L");
  gr_S_max_int_Z_18->Draw("L");
  gr_S_max_int_Z_19->Draw("L");
  gr_S_max_int_Z_20->Draw("L");
  gr_S_max_int_Z_21->Draw("L");
  gr_S_max_int_Z_22->Draw("L");
  gr_S_max_int_Z_23->Draw("L");
  gr_S_max_int_Z_24->Draw("L");
  gr_S_max_int_Z_25->Draw("L");
  gr_S_max_int_Z_26->Draw("L");
  gr_S_max_int_Z_27->Draw("L");
  gr_S_max_int_Z_28->Draw("L");
  gr_S_max_int_Z_29->Draw("L");
  gr_S_max_int_Z_30->Draw("L");
  gr_S_max_int_Z_31->Draw("L");
  gr_S_max_int_Z_32->Draw("L");

                                   
//	p_Z_5_S_max->Draw();


// Getting Solar Minimum
	
  in_1.open(spec_smin);
  
	in_1.getline(header,100);
//	printf("%s\n",header);
//n_1 >> header;

	E_index=0;

	for(Z_index=0;Z_index<33;Z_index++)
		{
	  	in_1 >> S_min[Z_index].S;
	  	// printf("%s\n",S_min[Z_index].S);
		}
		
  while(1) {
  
  	for(Z_index=0;Z_index<33;Z_index++)
		{
			in_1 >> S_min[Z_index].EI[E_index];
			// printf("S_min[%d][%d] = %e\n",Z_index,E_index,S_min[Z_index].EI[E_index]);
		}

    if(in_1.eof() || !in_1.good()) break;
		E_index++;		
  }
  in_1.close();
  in_1.clear();
  

  
  // Solar min Plot Label
  
  double label_x_min[2];
  double label_y_min[2];
  label_x_min[0] = 1e1;
  label_x_min[1] = 1e5;
  label_y_min[0] = 1e-17;
  label_y_min[1] = 1e-6;

 
  TGraph *gr_label_S_min_spectra = new TGraph(2,label_x_min,label_y_min);

  gr_label_S_min_spectra->SetMarkerStyle(1);
  gr_label_S_min_spectra->SetMarkerSize(1);
  gr_label_S_min_spectra->SetTitle("ACE-CRIS Solar Minimum Spectra");
  gr_label_S_min_spectra->GetXaxis()->SetTitle("MeV/nuc");
  gr_label_S_min_spectra->GetXaxis()->CenterTitle();
  gr_label_S_min_spectra->GetYaxis()->SetTitle("1/(cm^{2} s sr MeV/nuc)");
  gr_label_S_min_spectra->GetYaxis()->CenterTitle();
  gr_label_S_min_spectra->GetXaxis()->SetTitleOffset(1.2);
  gr_label_S_min_spectra->GetYaxis()->SetTitleOffset(1.4);
  gr_label_S_min_spectra->GetYaxis()->SetLabelFont(nicefont);
  gr_label_S_min_spectra->GetYaxis()->SetTitleFont(nicefont);
  gr_label_S_min_spectra->GetXaxis()->SetLabelFont(nicefont);
  gr_label_S_min_spectra->GetXaxis()->SetTitleFont(nicefont);
  gr_label_S_min_spectra->GetXaxis()->SetRangeUser(1e1,1e5);
  gr_label_S_min_spectra->GetYaxis()->SetRangeUser(1e-17,1e-6);

  // Solar Min Integral Plot Label
  
  double label_min_int_x[2];
  double label_min_int_y[2];
  label_min_int_x[0] = 1e1;
  label_min_int_x[1] = 1e5;
  label_min_int_y[0] = 1e-13;
  label_min_int_y[1] = 1e-2;

 
  TGraph *gr_label_S_min_int_spectra = new TGraph(2,label_min_int_x,label_min_int_y);

  gr_label_S_min_int_spectra->SetMarkerStyle(1);
  gr_label_S_min_int_spectra->SetMarkerSize(1);
  gr_label_S_min_int_spectra->SetTitle("ACE-CRIS Solar Minimum Integral Spectra");
  gr_label_S_min_int_spectra->GetXaxis()->SetTitle("MeV/nuc");
  gr_label_S_min_int_spectra->GetXaxis()->CenterTitle();
  gr_label_S_min_int_spectra->GetYaxis()->SetTitle("1/(cm^{2} s sr)");
  gr_label_S_min_int_spectra->GetYaxis()->CenterTitle();
  gr_label_S_min_int_spectra->GetXaxis()->SetTitleOffset(1.2);
  gr_label_S_min_int_spectra->GetYaxis()->SetTitleOffset(1.4);
  gr_label_S_min_int_spectra->GetYaxis()->SetLabelFont(nicefont);
  gr_label_S_min_int_spectra->GetYaxis()->SetTitleFont(nicefont);
  gr_label_S_min_int_spectra->GetXaxis()->SetLabelFont(nicefont);
  gr_label_S_min_int_spectra->GetXaxis()->SetTitleFont(nicefont);
  gr_label_S_min_int_spectra->GetXaxis()->SetRangeUser(1e1,1e5);
  gr_label_S_min_int_spectra->GetYaxis()->SetRangeUser(1e-13,1e-2);
  
  // Solar min Plot Spectra

  // Z = 5, Boron 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[5].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_5 = new TGraph(plot_count,x_plot,y_plot);
  S_min[5].I_spec = Integralfunc(gr_S_min_Z_5,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[5].I_spec: %e\n",S_min[5].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[5].EI_int[i]=Integralfunc(gr_S_min_Z_5,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[5].EI_int[%d]: %e\n",i,S_min[5].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[5].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_5 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[5].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_5_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_5_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[5].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[5].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[5].I_tot=S_min[5].I_spec + S_min[5].I_E_max;
	S_min[5].c_pl=pow(10.0,parp1[0]);
	S_min[5].i_pl=parp1[1];
       
  
  TPaveText *p_Z_5_S_min = new TPaveText(0.826,0.340,0.866,0.380,"NDC");
  p_Z_5_S_min->SetBorderSize(0);
  p_Z_5_S_min->SetFillColor(kWhite);
  TText *t_Z_5_S_min = p_Z_5_S_min->AddText(S_min[5].S);
  t_Z_5_S_min->SetTextFont(nicefont);
  t_Z_5_S_min->SetTextColor(1);
  
  
  // Z = 6, Carbon 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[6].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_6 = new TGraph(plot_count,x_plot,y_plot);
  S_min[6].I_spec = Integralfunc(gr_S_min_Z_6,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[6].I_spec: %e\n",S_min[6].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[6].EI_int[i]=Integralfunc(gr_S_min_Z_6,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[6].EI_int[%d]: %e\n",i,S_min[6].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[6].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_6 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[6].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_6_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_6_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[6].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[6].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[6].I_tot=S_min[6].I_spec + S_min[6].I_E_max;
	S_min[6].c_pl=pow(10.0,parp1[0]);
	S_min[6].i_pl=parp1[1];
      
  TPaveText *p_Z_6_S_min = new TPaveText(0.826,0.340,0.866,0.380,"NDC");
  p_Z_6_S_min->SetBorderSize(0);
  p_Z_6_S_min->SetFillColor(kWhite);
  TText *t_Z_6_S_min = p_Z_6_S_min->AddText(S_min[6].S);
  t_Z_6_S_min->SetTextFont(nicefont);
  t_Z_6_S_min->SetTextColor(1);
   
  // Z = 7, Nitrogen 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[7].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_7 = new TGraph(plot_count,x_plot,y_plot);
  S_min[7].I_spec = Integralfunc(gr_S_min_Z_7,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[7].I_spec: %e\n",S_min[7].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[7].EI_int[i]=Integralfunc(gr_S_min_Z_7,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[7].EI_int[%d]: %e\n",i,S_min[7].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[7].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_7 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[7].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_7_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_7_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[7].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[7].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[7].I_tot=S_min[7].I_spec + S_min[7].I_E_max;
	S_min[7].c_pl=pow(10.0,parp1[0]);
	S_min[7].i_pl=parp1[1];
    
  TPaveText *p_Z_7_S_min = new TPaveText(0.827,0.340,0.877,0.380,"NDC");
  p_Z_7_S_min->SetBorderSize(0);
  p_Z_7_S_min->SetFillColor(kWhite);
  TText *t_Z_7_S_min = p_Z_7_S_min->AddText(S_min[7].S);
  t_Z_7_S_min->SetTextFont(nicefont);
  t_Z_7_S_min->SetTextColor(1);
  
    // Z = 8, Oxygen 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[8].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_8 = new TGraph(plot_count,x_plot,y_plot);
  S_min[8].I_spec = Integralfunc(gr_S_min_Z_8,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[8].I_spec: %e\n",S_min[8].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[8].EI_int[i]=Integralfunc(gr_S_min_Z_8,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[8].EI_int[%d]: %e\n",i,S_min[8].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[8].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_8 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[8].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_8_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_8_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[8].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[8].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[8].I_tot=S_min[8].I_spec + S_min[8].I_E_max;
	S_min[8].c_pl=pow(10.0,parp1[0]);
	S_min[8].i_pl=parp1[1];
   
    // Z = 9, Flourine
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[9].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_9 = new TGraph(plot_count,x_plot,y_plot);
  S_min[9].I_spec = Integralfunc(gr_S_min_Z_9,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[9].I_spec: %e\n",S_min[9].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[9].EI_int[i]=Integralfunc(gr_S_min_Z_9,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[9].EI_int[%d]: %e\n",i,S_min[9].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[9].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_9 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[9].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_9_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_9_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[9].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[9].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[9].I_tot=S_min[9].I_spec + S_min[9].I_E_max;
	S_min[9].c_pl=pow(10.0,parp1[0]);
	S_min[9].i_pl=parp1[1];
    
    // Z = 10, Neon 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[10].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_10 = new TGraph(plot_count,x_plot,y_plot);  
  S_min[10].I_spec = Integralfunc(gr_S_min_Z_10,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[10].I_spec: %e\n",S_min[10].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[10].EI_int[i]=Integralfunc(gr_S_min_Z_10,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[10].EI_int[%d]: %e\n",i,S_min[10].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[10].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_10 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[10].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_10_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_10_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[10].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[10].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[10].I_tot=S_min[10].I_spec + S_min[10].I_E_max;
	S_min[10].c_pl=pow(10.0,parp1[0]);
	S_min[10].i_pl=parp1[1];
    
    // Z = 11, Sodium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[11].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_11 = new TGraph(plot_count,x_plot,y_plot);
  S_min[11].I_spec = Integralfunc(gr_S_min_Z_11,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[11].I_spec: %e\n",S_min[11].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[11].EI_int[i]=Integralfunc(gr_S_min_Z_11,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[11].EI_int[%d]: %e\n",i,S_min[11].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[11].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_11 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[11].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_11_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_11_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[11].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[11].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[11].I_tot=S_min[11].I_spec + S_min[11].I_E_max;
	S_min[11].c_pl=pow(10.0,parp1[0]);
	S_min[11].i_pl=parp1[1];
    
    // Z = 12, Magnesium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[12].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_12 = new TGraph(plot_count,x_plot,y_plot);  
  S_min[12].I_spec = Integralfunc(gr_S_min_Z_12,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[12].I_spec: %e\n",S_min[12].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[12].EI_int[i]=Integralfunc(gr_S_min_Z_12,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[12].EI_int[%d]: %e\n",i,S_min[12].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[12].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_12 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[12].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_12_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_12_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[12].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[12].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[12].I_tot=S_min[12].I_spec + S_min[12].I_E_max;
	S_min[12].c_pl=pow(10.0,parp1[0]);
	S_min[12].i_pl=parp1[1];
   
    // Z = 13, Aluminum 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[13].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_13 = new TGraph(plot_count,x_plot,y_plot);  
  S_min[13].I_spec = Integralfunc(gr_S_min_Z_13,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[13].I_spec: %e\n",S_min[13].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[13].EI_int[i]=Integralfunc(gr_S_min_Z_13,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[13].EI_int[%d]: %e\n",i,S_min[13].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[13].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_13 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[13].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_13_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_13_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[13].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[13].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[13].I_tot=S_min[13].I_spec + S_min[13].I_E_max;
	S_min[13].c_pl=pow(10.0,parp1[0]);
	S_min[13].i_pl=parp1[1];
   
    // Z = 14, Silicon 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[14].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_14 = new TGraph(plot_count,x_plot,y_plot);
  S_min[14].I_spec = Integralfunc(gr_S_min_Z_14,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[14].I_spec: %e\n",S_min[14].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[14].EI_int[i]=Integralfunc(gr_S_min_Z_14,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[14].EI_int[%d]: %e\n",i,S_min[14].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[14].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_14 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[14].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_14_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_14_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[14].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[14].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[14].I_tot=S_min[14].I_spec + S_min[14].I_E_max;
	S_min[14].c_pl=pow(10.0,parp1[0]);
	S_min[14].i_pl=parp1[1];
    
    // Z = 15, Phosphorous
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[15].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_15 = new TGraph(plot_count,x_plot,y_plot);
  S_min[15].I_spec = Integralfunc(gr_S_min_Z_15,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[15].I_spec: %e\n",S_min[15].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[15].EI_int[i]=Integralfunc(gr_S_min_Z_15,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[15].EI_int[%d]: %e\n",i,S_min[15].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[15].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_15 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[15].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_15_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_15_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[15].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[15].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[15].I_tot=S_min[15].I_spec + S_min[15].I_E_max;
	S_min[15].c_pl=pow(10.0,parp1[0]);
	S_min[15].i_pl=parp1[1];
   
    // Z = 16, Sulfur
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[16].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_16 = new TGraph(plot_count,x_plot,y_plot);
  S_min[16].I_spec = Integralfunc(gr_S_min_Z_16,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[16].I_spec: %e\n",S_min[16].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[16].EI_int[i]=Integralfunc(gr_S_min_Z_16,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[16].EI_int[%d]: %e\n",i,S_min[16].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[16].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_16 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[16].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_16_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_16_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[16].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[16].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[16].I_tot=S_min[16].I_spec + S_min[16].I_E_max;
	S_min[16].c_pl=pow(10.0,parp1[0]);
	S_min[16].i_pl=parp1[1];
   
    // Z = 17, Chlorine 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[17].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_17 = new TGraph(plot_count,x_plot,y_plot);
  S_min[17].I_spec = Integralfunc(gr_S_min_Z_17,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[17].I_spec: %e\n",S_min[17].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[17].EI_int[i]=Integralfunc(gr_S_min_Z_17,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[17].EI_int[%d]: %e\n",i,S_min[17].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[17].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_17 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[17].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_17_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_17_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[17].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[17].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[17].I_tot=S_min[17].I_spec + S_min[17].I_E_max;
	S_min[17].c_pl=pow(10.0,parp1[0]);
	S_min[17].i_pl=parp1[1];
   
    // Z = 18, Argon 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[18].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_18 = new TGraph(plot_count,x_plot,y_plot);
  S_min[18].I_spec = Integralfunc(gr_S_min_Z_18,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[18].I_spec: %e\n",S_min[18].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[18].EI_int[i]=Integralfunc(gr_S_min_Z_18,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[18].EI_int[%d]: %e\n",i,S_min[18].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[18].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_18 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[18].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_18_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_18_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[18].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[18].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[18].I_tot=S_min[18].I_spec + S_min[18].I_E_max;
	S_min[18].c_pl=pow(10.0,parp1[0]);
	S_min[18].i_pl=parp1[1];
 
    // Z = 19, Potassium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[19].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_19 = new TGraph(plot_count,x_plot,y_plot);
  S_min[19].I_spec = Integralfunc(gr_S_min_Z_19,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[19].I_spec: %e\n",S_min[19].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[19].EI_int[i]=Integralfunc(gr_S_min_Z_19,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[19].EI_int[%d]: %e\n",i,S_min[19].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[19].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_19 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[19].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_19_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_19_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[19].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[19].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[19].I_tot=S_min[19].I_spec + S_min[19].I_E_max;
	S_min[19].c_pl=pow(10.0,parp1[0]);
	S_min[19].i_pl=parp1[1];
   
    // Z = 20, Calcium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[20].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_20 = new TGraph(plot_count,x_plot,y_plot);
  S_min[20].I_spec = Integralfunc(gr_S_min_Z_20,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[20].I_spec: %e\n",S_min[20].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[20].EI_int[i]=Integralfunc(gr_S_min_Z_20,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[20].EI_int[%d]: %e\n",i,S_min[20].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[20].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_20 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[20].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_20_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_20_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[20].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[20].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[20].I_tot=S_min[20].I_spec + S_min[20].I_E_max;
	S_min[20].c_pl=pow(10.0,parp1[0]);
	S_min[20].i_pl=parp1[1];
   
    // Z = 21, Scandium
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[21].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_21 = new TGraph(plot_count,x_plot,y_plot);
  S_min[21].I_spec = Integralfunc(gr_S_min_Z_21,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[21].I_spec: %e\n",S_min[21].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[21].EI_int[i]=Integralfunc(gr_S_min_Z_21,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[21].EI_int[%d]: %e\n",i,S_min[21].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[21].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_21 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[21].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_21_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_21_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[21].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[21].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[21].I_tot=S_min[21].I_spec + S_min[21].I_E_max;
	S_min[21].c_pl=pow(10.0,parp1[0]);
	S_min[21].i_pl=parp1[1];
   
    // Z = 22, Titanium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[22].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_22 = new TGraph(plot_count,x_plot,y_plot);
  S_min[22].I_spec = Integralfunc(gr_S_min_Z_22,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[22].I_spec: %e\n",S_min[22].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[22].EI_int[i]=Integralfunc(gr_S_min_Z_22,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[22].EI_int[%d]: %e\n",i,S_min[22].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[22].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_22 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[22].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_22_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_22_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[22].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[22].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[22].I_tot=S_min[22].I_spec + S_min[22].I_E_max;
	S_min[22].c_pl=pow(10.0,parp1[0]);
	S_min[22].i_pl=parp1[1];
   
    // Z = 23, Vanadium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[23].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_23 = new TGraph(plot_count,x_plot,y_plot);
  S_min[23].I_spec = Integralfunc(gr_S_min_Z_23,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[23].I_spec: %e\n",S_min[23].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[23].EI_int[i]=Integralfunc(gr_S_min_Z_23,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[23].EI_int[%d]: %e\n",i,S_min[23].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[23].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_23 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[23].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_23_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_23_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[23].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[23].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[23].I_tot=S_min[23].I_spec + S_min[23].I_E_max;
	S_min[23].c_pl=pow(10.0,parp1[0]);
	S_min[23].i_pl=parp1[1];
   
    // Z = 24, Chromium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[24].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_24 = new TGraph(plot_count,x_plot,y_plot);
  S_min[24].I_spec = Integralfunc(gr_S_min_Z_24,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[24].I_spec: %e\n",S_min[24].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[24].EI_int[i]=Integralfunc(gr_S_min_Z_24,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[24].EI_int[%d]: %e\n",i,S_min[24].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[24].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_24 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[24].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_24_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_24_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[24].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[24].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[24].I_tot=S_min[24].I_spec + S_min[24].I_E_max;
	S_min[24].c_pl=pow(10.0,parp1[0]);
	S_min[24].i_pl=parp1[1];
   
    // Z = 25, Manganese 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[25].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_25 = new TGraph(plot_count,x_plot,y_plot);  
  S_min[25].I_spec = Integralfunc(gr_S_min_Z_25,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[25].I_spec: %e\n",S_min[25].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[25].EI_int[i]=Integralfunc(gr_S_min_Z_25,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[25].EI_int[%d]: %e\n",i,S_min[25].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[25].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_25 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[25].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_25_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_25_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[25].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[25].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[25].I_tot=S_min[25].I_spec + S_min[25].I_E_max;
	S_min[25].c_pl=pow(10.0,parp1[0]);
	S_min[25].i_pl=parp1[1];
   
    // Z = 26, Iron 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[26].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_26 = new TGraph(plot_count,x_plot,y_plot);    
  S_min[26].I_spec = Integralfunc(gr_S_min_Z_26,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[26].I_spec: %e\n",S_min[26].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[26].EI_int[i]=Integralfunc(gr_S_min_Z_26,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[26].EI_int[%d]: %e\n",i,S_min[26].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[26].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_26 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[26].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_26_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_26_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[26].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[26].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[26].I_tot=S_min[26].I_spec + S_min[26].I_E_max;
	S_min[26].c_pl=pow(10.0,parp1[0]);
	S_min[26].i_pl=parp1[1];
   
    // Z = 27, Cobalt 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[27].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_27 = new TGraph(plot_count,x_plot,y_plot);    
  S_min[27].I_spec = Integralfunc(gr_S_min_Z_27,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[27].I_spec: %e\n",S_min[27].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[27].EI_int[i]=Integralfunc(gr_S_min_Z_27,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[27].EI_int[%d]: %e\n",i,S_min[27].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[27].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_27 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[27].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_27_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_27_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[27].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[27].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[27].I_tot=S_min[27].I_spec + S_min[27].I_E_max;
	S_min[27].c_pl=pow(10.0,parp1[0]);
	S_min[27].i_pl=parp1[1];
   
    // Z = 28, Nickel 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[28].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_28 = new TGraph(plot_count,x_plot,y_plot);    
  S_min[28].I_spec = Integralfunc(gr_S_min_Z_28,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[28].I_spec: %e\n",S_min[28].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[28].EI_int[i]=Integralfunc(gr_S_min_Z_28,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[28].EI_int[%d]: %e\n",i,S_min[28].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[28].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_28 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[28].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_28_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_28_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[28].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[28].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[28].I_tot=S_min[28].I_spec + S_min[28].I_E_max;
	S_min[28].c_pl=pow(10.0,parp1[0]);
	S_min[28].i_pl=parp1[1];
   
    // Z = 29, Copper 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[29].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_29 = new TGraph(plot_count,x_plot,y_plot);    
  S_min[29].I_spec = Integralfunc(gr_S_min_Z_29,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[29].I_spec: %e\n",S_min[29].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[29].EI_int[i]=Integralfunc(gr_S_min_Z_29,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[29].EI_int[%d]: %e\n",i,S_min[29].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[29].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_29 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[29].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_29_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_29_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[29].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[29].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[29].I_tot=S_min[29].I_spec + S_min[29].I_E_max;
	S_min[29].c_pl=pow(10.0,parp1[0]);
	S_min[29].i_pl=parp1[1];
   
    // Z = 30, Zink 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[30].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_30 = new TGraph(plot_count,x_plot,y_plot);    
  S_min[30].I_spec = Integralfunc(gr_S_min_Z_30,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[30].I_spec: %e\n",S_min[30].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[30].EI_int[i]=Integralfunc(gr_S_min_Z_30,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[30].EI_int[%d]: %e\n",i,S_min[30].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[30].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_30 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[30].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_30_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_30_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[30].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[30].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[30].I_tot=S_min[30].I_spec + S_min[30].I_E_max;
	S_min[30].c_pl=pow(10.0,parp1[0]);
	S_min[30].i_pl=parp1[1];
   
    // Z = 31, Gallium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[31].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_31 = new TGraph(plot_count,x_plot,y_plot);      
  S_min[31].I_spec = Integralfunc(gr_S_min_Z_31,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[31].I_spec: %e\n",S_min[31].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[31].EI_int[i]=Integralfunc(gr_S_min_Z_31,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[31].EI_int[%d]: %e\n",i,S_min[31].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[31].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_31 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[31].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_31_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_31_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[31].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[31].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[31].I_tot=S_min[31].I_spec + S_min[31].I_E_max;
	S_min[31].c_pl=pow(10.0,parp1[0]);
	S_min[31].i_pl=parp1[1];
   
    // Z = 32, Germanium 
  
  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[32].EI[i];
		plot_count++;
		}
  TGraph *gr_S_min_Z_32 = new TGraph(plot_count,x_plot,y_plot);    
  S_min[32].I_spec = Integralfunc(gr_S_min_Z_32,S_min[0].EI[0],S_min[0].EI[E_index-1],"5",epsilon);

	// Generating Integral Spectrum
	// fprintf(stderr,"S_min[32].I_spec: %e\n",S_min[32].I_spec);
	for(i=0;i<E_index;i++)
		{
		S_min[32].EI_int[i]=Integralfunc(gr_S_min_Z_32,S_min[0].EI[i],S_min[0].EI[E_index-1],"5",epsilon);
		// fprintf(stderr,"S_min[32].EI_int[%d]: %e\n",i,S_min[32].EI_int[i]);
		}

  plot_count=0;
  for(i=0;i<E_index;i++) {
    x_plot[i]=S_min[0].EI[i];
    y_plot[i]=S_min[32].EI_int[i];
		plot_count++;
		}
  TGraph *gr_S_min_int_Z_32 = new TGraph(plot_count,x_plot,y_plot);
	// End Integral Spectrum

  plot_count=0;
  for(i=0;i<tail_int;i++) {
 	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
 	  y_plot[i]=log10(S_min[32].EI[E_index-tail_int+i]);
		// fprintf(stderr,"x_plot[%d]: %e y_plot[%d]: %e\n",i,x_plot[i],i,y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_Z_32_tail = new TGraph(plot_count,x_plot,y_plot);  
	
	gr_S_min_Z_32_tail->Fit("fp1","Q");
	fp1->GetParameters(parp1);
	S_min[32].I_E_max=(-pow(10.0,parp1[0])/(parp1[1]+1))*pow(E_spec_hi,parp1[1]+1.0);
	S_min[32].I_E_hi=(-pow(10.0,parp1[0])/(parp1[1]+1))*(pow(E_spec_hi,(parp1[1]+1.0)) - pow(E_CALET_hi,(parp1[1]+1.0)));
	S_min[32].I_tot=S_min[32].I_spec + S_min[32].I_E_max;
	S_min[32].c_pl=pow(10.0,parp1[0]);
	S_min[32].i_pl=parp1[1];
   
  
  TCanvas *C_S_min_spectra = new TCanvas("C_S_min_spectra","ACE-CRIS S_Min Spectra",0,0,800,800);
  
  C_S_min_spectra->SetLeftMargin(0.124365);
  C_S_min_spectra->SetRightMargin(0.0748731);
  C_S_min_spectra->SetBottomMargin(0.117801);
  C_S_min_spectra->SetTopMargin(0.0824607);

  C_S_min_spectra->SetGrid();
  
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  
  gr_label_S_min_spectra->Draw("AP");

  gr_S_min_Z_5->Draw("PL");
  gr_S_min_Z_6->Draw("PL");
  gr_S_min_Z_7->Draw("PL");
  gr_S_min_Z_8->Draw("PL");
  gr_S_min_Z_9->Draw("PL");
  gr_S_min_Z_10->Draw("PL");
  gr_S_min_Z_11->Draw("PL");
  gr_S_min_Z_12->Draw("PL");
  gr_S_min_Z_13->Draw("PL");
  gr_S_min_Z_14->Draw("PL");
  gr_S_min_Z_15->Draw("PL");
  gr_S_min_Z_16->Draw("PL");
  gr_S_min_Z_17->Draw("PL");
  gr_S_min_Z_18->Draw("PL");
  gr_S_min_Z_19->Draw("PL");
  gr_S_min_Z_20->Draw("PL");
  gr_S_min_Z_21->Draw("PL");
  gr_S_min_Z_22->Draw("PL");
  gr_S_min_Z_23->Draw("PL");
  gr_S_min_Z_24->Draw("PL");
  gr_S_min_Z_25->Draw("PL");
  gr_S_min_Z_26->Draw("PL");
  gr_S_min_Z_27->Draw("PL");
  gr_S_min_Z_28->Draw("PL");
  gr_S_min_Z_29->Draw("PL");
  gr_S_min_Z_30->Draw("PL");
  gr_S_min_Z_31->Draw("PL");
  gr_S_min_Z_32->Draw("PL");


  TCanvas *C_S_min_int_spectra = new TCanvas("C_S_min_int_spectra","ACE-CRIS S_min Spectra",0,0,800,800);
  
  C_S_min_int_spectra->SetLeftMargin(0.124365);
  C_S_min_int_spectra->SetRightMargin(0.0748731);
  C_S_min_int_spectra->SetBottomMargin(0.117801);
  C_S_min_int_spectra->SetTopMargin(0.0824607);

  C_S_min_int_spectra->SetGrid();
  
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  
  gr_label_S_min_int_spectra->Draw("AP");

  gr_S_min_int_Z_5->Draw("L");
  gr_S_min_int_Z_6->Draw("L");
  gr_S_min_int_Z_7->Draw("L");
  gr_S_min_int_Z_8->Draw("L");
  gr_S_min_int_Z_9->Draw("L");
  gr_S_min_int_Z_10->Draw("L");
  gr_S_min_int_Z_11->Draw("L");
  gr_S_min_int_Z_12->Draw("L");
  gr_S_min_int_Z_13->Draw("L");
  gr_S_min_int_Z_14->Draw("L");
  gr_S_min_int_Z_15->Draw("L");
  gr_S_min_int_Z_16->Draw("L");
  gr_S_min_int_Z_17->Draw("L");
  gr_S_min_int_Z_18->Draw("L");
  gr_S_min_int_Z_19->Draw("L");
  gr_S_min_int_Z_20->Draw("L");
  gr_S_min_int_Z_21->Draw("L");
  gr_S_min_int_Z_22->Draw("L");
  gr_S_min_int_Z_23->Draw("L");
  gr_S_min_int_Z_24->Draw("L");
  gr_S_min_int_Z_25->Draw("L");
  gr_S_min_int_Z_26->Draw("L");
  gr_S_min_int_Z_27->Draw("L");
  gr_S_min_int_Z_28->Draw("L");
  gr_S_min_int_Z_29->Draw("L");
  gr_S_min_int_Z_30->Draw("L");
  gr_S_min_int_Z_31->Draw("L");
  gr_S_min_int_Z_32->Draw("L");

                                  
	
	
		// Interaction Correction Factors for Scintillator Material
	
	
		// material properties
		// PVT: C9H10
		// PS: C8H8
		// PMMA: C5H8O2   
		// Glass: SiO2
		// Aluminum: Al
	
	int n_PVT_H = 10;
	int n_PVT_C = 9;
	
	int n_PS_H = 8;
	int n_PS_C = 8;
	
	int n_PMMA_H = 8;
	int n_PMMA_C = 5;
	int n_PMMA_O = 2;
	
	int n_Glass_Si = 1;
	int n_Glass_O = 2;
	
	int n_Al = 1;
	
		// Correction Factors
	
	
		// Target Radii
	
	double R_N=0;
	double R_O=0;
	double R_H=0;
	double R_C=0;
	double R_Si=0;
	double R_Al=0;
	
	R_N = 1.58*pow(el_dat[6].A,0.281);
	R_O = 1.58*pow(el_dat[7].A,0.281);
		//  R_H = 1.58*pow(A_H,0.281);  // forrmula does NOT work for H
	R_H = 1.034;
	R_C = 1.58*pow(el_dat[5].A,0.281);
	R_Si = 1.58*pow(el_dat[13].A,0.281);
	R_Al = 1.58*pow(el_dat[12].A,0.281);
	
	double angle_fac=0.0;
	angle_fac=1/cos(30.0*(TMath::Pi()/180.0));
	
	double density_PVT=1.032; // g/cm^3
	double width_PVT=2.0; // Two layers of EJ-204 Plastic Scintillator, here called PVT
	
	double x_PVT=0.0;
	x_PVT = angle_fac*(density_PVT*width_PVT);
	fprintf(stderr,"x_PVT: %.2f\n",x_PVT);
		// Material mean free paths
	
		// initial[index].R = 1.58*pow(initial[index].A,0.281);
	double interaction_Z[100];
	double interaction_survival_fraction[100];
	
	for(i=0;i<el_index;i++) {
		
		el_dat[i].sig_tot_H_eq_9 = 0.01*TMath::Pi()*pow((S_fac*R_H + S_fac*el_dat[i].R - 3.2),2);
		el_dat[i].sig_tot_C_eq_9 = 0.01*TMath::Pi()*pow((S_fac*R_C + S_fac*el_dat[i].R - 3.2),2);
		
		el_dat[i].lambda_PVT = 1.6624*(n_PVT_H*el_dat[0].A + n_PVT_C*el_dat[11].A)/
		(n_PVT_H*el_dat[i].sig_tot_H_eq_9 + n_PVT_C*el_dat[i].sig_tot_C_eq_9);
		el_dat[i].I_factor = exp(-x_PVT/el_dat[i].lambda_PVT);
		interaction_Z[i]=el_dat[i].Z;
		interaction_survival_fraction[i]=el_dat[i].I_factor;
		fprintf(stderr,"el_dat[%2d] Z %2d A: %6.2f lambda_PVT: %6.2f I_factor: %.6f\n",i,el_dat[i].Z,el_dat[i].A,el_dat[i].lambda_PVT,el_dat[i].I_factor);
		
	}
	
		// Plot Interaction Survival Fractions as a Function of Z
	
	
	
	TGraph *gr_survival_fraction = new TGraph(el_index,interaction_Z,interaction_survival_fraction);
	
	gr_survival_fraction->SetMarkerStyle(1);
	gr_survival_fraction->SetMarkerSize(1);
	gr_survival_fraction->SetTitle("CHD Interaction Survival Fractions");
	gr_survival_fraction->GetXaxis()->SetTitle("Z");
	gr_survival_fraction->GetXaxis()->CenterTitle();
	gr_survival_fraction->GetYaxis()->SetTitle("Survival Fraction");
	gr_survival_fraction->GetYaxis()->CenterTitle();
	gr_survival_fraction->GetXaxis()->SetTitleOffset(1.2);
	gr_survival_fraction->GetYaxis()->SetTitleOffset(1.2);
	gr_survival_fraction->GetYaxis()->SetLabelFont(nicefont);
	gr_survival_fraction->GetYaxis()->SetTitleFont(nicefont);
	gr_survival_fraction->GetXaxis()->SetLabelFont(nicefont);
	gr_survival_fraction->GetXaxis()->SetTitleFont(nicefont);
	gr_survival_fraction->GetXaxis()->SetRangeUser(0,100);
	gr_survival_fraction->GetYaxis()->SetRangeUser(0.7,1);
	

	
	TCanvas *C_survival_fraction = new TCanvas("C_survival_fraction","Interaction Survival Fraction",0,0,800,800);
	
		//	C_survival_fraction->SetLeftMargin(0.124365);
		//	C_survival_fraction->SetRightMargin(0.0748731);
		//	C_survival_fraction->SetBottomMargin(0.117801);
		//	C_survival_fraction->SetTopMargin(0.0824607);
	
  	C_survival_fraction->SetGrid();
	gr_survival_fraction->SetLineWidth(2);
	gr_survival_fraction->Draw("AL");
	
	
	

// Plot specific Comparison Graphs

	
	char tstring[100];
	char tstring_h[100];


  if(argc==2 && plot_Z>0 && plot_Z<33)
		{


		sprintf(tstring,"ACE-CRIS Solar Minimum and Maximum Spectra for _{%d}%s",plot_Z,S_min[plot_Z].S);


	  TGraph *gr_label_comp = new TGraph(2,label_x_min,label_y_min);

		gr_label_comp->SetMarkerStyle(1);
		gr_label_comp->SetMarkerSize(1);
		gr_label_comp->SetTitle(tstring);
		gr_label_comp->GetXaxis()->SetTitle("MeV/nuc");
		gr_label_comp->GetXaxis()->CenterTitle();
		gr_label_comp->GetYaxis()->SetTitle("1/(cm^{2} s sr MeV/nuc)");
		gr_label_comp->GetYaxis()->CenterTitle();
		gr_label_comp->GetXaxis()->SetTitleOffset(1.2);
		gr_label_comp->GetYaxis()->SetTitleOffset(1.4);
		gr_label_comp->GetYaxis()->SetLabelFont(nicefont);
		gr_label_comp->GetYaxis()->SetTitleFont(nicefont);
		gr_label_comp->GetXaxis()->SetLabelFont(nicefont);
		gr_label_comp->GetXaxis()->SetTitleFont(nicefont);
		gr_label_comp->GetXaxis()->SetRangeUser(1e1,1e5);
		gr_label_S_min_spectra->GetYaxis()->SetRangeUser(1e-17,1e-6);

  
	  plot_count=0;
	  for(i=0;i<E_index;i++) {
  	  x_plot[i]=S_min[0].EI[i];
  	  y_plot[i]=S_min[plot_Z].EI[i];
			plot_count++;
			}
	  TGraph *gr_S_min_plot_Z = new TGraph(plot_count,x_plot,y_plot);    

	  plot_count=0;
	  for(i=0;i<tail_int;i++) {
  	  x_plot[i]=log10(S_min[0].EI[E_index-tail_int+i]);
  	  y_plot[i]=log10(S_min[plot_Z].EI[E_index-tail_int+i]);
			plot_count++;
			}
	  TGraph *gr_S_min_plot_Z_tail = new TGraph(plot_count,x_plot,y_plot);  


		gr_S_min_plot_Z_tail->Fit("fp1","Q");
		fp1->GetParameters(parp1);
		fprintf(stderr,"%e\n%f\n",pow(10.0,parp1[0]),parp1[1]);


	  plot_count=0;
	  for(i=0;i<E_index;i++) {
  	  x_plot[i]=S_max[0].EI[i];
  	  y_plot[i]=S_max[plot_Z].EI[i];
			plot_count++;
			}
	  TGraph *gr_S_max_plot_Z = new TGraph(plot_count,x_plot,y_plot);  

  
	  TCanvas *C_plot_Z = new TCanvas("C_plot_Z","ACE-CRIS Comparison Spectra",0,0,800,800);
  
	  C_plot_Z->SetLeftMargin(0.124365);
	  C_plot_Z->SetRightMargin(0.0748731);
	  C_plot_Z->SetBottomMargin(0.117801);
	  C_plot_Z->SetTopMargin(0.0824607);

  	C_plot_Z->SetGrid();
  
	  gPad->SetLogx(1);
	  gPad->SetLogy(1);
  
	  gr_label_comp->Draw("AP");
			gr_S_max_plot_Z->SetLineColor(2);
			gr_S_max_plot_Z->SetLineWidth(2);
			gr_S_max_plot_Z->Draw("PL");
			gr_S_min_plot_Z->SetLineWidth(2);
			gr_S_min_plot_Z->Draw("PL");
			
			
			TLegend *leg_S_comp = new TLegend(0.59,0.15,0.89,0.27);
			leg_S_comp->AddEntry(gr_S_max_plot_Z,"Solar Maximum","l");
			leg_S_comp->AddEntry(gr_S_min_plot_Z,"Solar Minimum","l");

			
			leg_S_comp->SetFillColor(10);
			leg_S_comp->SetTextFont(nicefont);
			leg_S_comp->SetTextSize(0.03);
			leg_S_comp->SetBorderSize(1);       // remove ugly border from stats 
			leg_S_comp->Draw();
			
			
			
			
				// Integral Spectra for selected charge
			char tstring_int[100];
			
			sprintf(tstring_int,"ACE-CRIS Solar Minimum and Maximum Integral Spectra for _{%d}%s",plot_Z,S_min[plot_Z].S);
			
			
			
			double label_x_int[2];
			double label_y_int[2];
			
			label_x_int[0]=1e1;
			label_x_int[1]=1e5;
			label_y_int[0]=1e-12;
			label_y_int[1]=1e-2;
			
			TGraph *gr_label_comp_int = new TGraph(2,label_x_int,label_y_int);
			
			gr_label_comp_int->SetMarkerStyle(1);
			gr_label_comp_int->SetMarkerSize(1);
			gr_label_comp_int->SetTitle(tstring);
			gr_label_comp_int->GetXaxis()->SetTitle("MeV/nuc");
			gr_label_comp_int->GetXaxis()->CenterTitle();
			gr_label_comp_int->GetYaxis()->SetTitle("1/(cm^{2} s sr)");
			gr_label_comp_int->GetYaxis()->CenterTitle();
			gr_label_comp_int->GetXaxis()->SetTitleOffset(1.2);
			gr_label_comp_int->GetYaxis()->SetTitleOffset(1.4);
			gr_label_comp_int->GetYaxis()->SetLabelFont(nicefont);
			gr_label_comp_int->GetYaxis()->SetTitleFont(nicefont);
			gr_label_comp_int->GetXaxis()->SetLabelFont(nicefont);
			gr_label_comp_int->GetXaxis()->SetTitleFont(nicefont);
			gr_label_comp_int->GetXaxis()->SetRangeUser(1e1,1e5);
			gr_label_S_min_spectra->GetYaxis()->SetRangeUser(1e-13,1e-3);
			
			plot_count=0;
			for(i=0;i<E_index;i++) {
				x_plot[i]=S_min[0].EI[i];
				y_plot[i]=S_min[plot_Z].EI_int[i];
				plot_count++;
			}
			TGraph *gr_S_min_plot_Z_int = new TGraph(plot_count,x_plot,y_plot);    

			
			plot_count=0;
			for(i=0;i<E_index;i++) {
				x_plot[i]=S_max[0].EI[i];
				y_plot[i]=S_max[plot_Z].EI_int[i];
				plot_count++;
			}
			TGraph *gr_S_max_plot_Z_int = new TGraph(plot_count,x_plot,y_plot);  
			
			
			TCanvas *C_plot_Z_int = new TCanvas("C_plot_Z_int","ACE-CRIS Comparison Integral Spectra",0,0,800,800);
			
			C_plot_Z_int->SetLeftMargin(0.124365);
			C_plot_Z_int->SetRightMargin(0.0748731);
			C_plot_Z_int->SetBottomMargin(0.117801);
			C_plot_Z_int->SetTopMargin(0.0824607);
			
			C_plot_Z_int->SetGrid();
			
			gPad->SetLogx(1);
			gPad->SetLogy(1);
			
			gr_label_comp_int->Draw("AP");
			gr_S_max_plot_Z_int->SetLineColor(2);
			gr_S_max_plot_Z_int->SetLineWidth(2);
			gr_S_max_plot_Z_int->Draw("PL");
			gr_S_min_plot_Z_int->SetLineWidth(2);
			gr_S_min_plot_Z_int->Draw("PL");
			
			
			TLegend *leg_S_comp_int = new TLegend(0.59,0.15,0.89,0.27);
			leg_S_comp_int->AddEntry(gr_S_max_plot_Z_int,"Solar Maximum","l");
			leg_S_comp_int->AddEntry(gr_S_min_plot_Z_int,"Solar Minimum","l");
			
			
			leg_S_comp_int->SetFillColor(10);
			leg_S_comp_int->SetTextFont(nicefont);
			leg_S_comp_int->SetTextSize(0.03);
			leg_S_comp_int->SetBorderSize(1);       // remove ugly border from stats 
			leg_S_comp_int->Draw();
			
			
				// Geomagnetic Cutoff Energy Distribution for selected element
			
			
			sprintf(tstring_h,"Geomagnetic Cutoff Energy for _{%d}%s",plot_Z,chem_s);
			
			TH1F *h_Z_plot_E = new TH1F("h_Z_plot_E",tstring_h,70,0,7000);
			
				// h_Z_plot_E
			h_Z_plot_E->GetXaxis()->SetTitle("Minimum Energy (MeV/nuc)");
			h_Z_plot_E->GetXaxis()->CenterTitle();
			h_Z_plot_E->GetYaxis()->SetTitle("Orbit Fraction");
			h_Z_plot_E->GetYaxis()->CenterTitle();
			h_Z_plot_E->GetYaxis()->SetTitleOffset(1.7);
			h_Z_plot_E->GetYaxis()->SetLabelFont(nicefont);
			h_Z_plot_E->GetYaxis()->SetTitleFont(nicefont);
			h_Z_plot_E->GetXaxis()->SetLabelFont(nicefont);
			h_Z_plot_E->GetXaxis()->SetTitleFont(nicefont);
			
			

		for(j=0;j<long_bins;j++)
			{
			for(k=0;k<lat_bins;k++)
				{
				if(dt[k] > 0.0)
					{
					h_Z_plot_E->Fill(el_dat[plot_Z-1].E_min_mag[j][k],dt[k]/long_bins);
					}
				}
			}
			
	  TCanvas *C_Z_plot_E = new TCanvas("C_Z_plot_E","Minimum Energies",0,0,800,800);
			C_Z_plot_E->SetLeftMargin(0.148477);
			C_Z_plot_E->SetRightMargin(0.0507614);
			C_Z_plot_E->SetBottomMargin(0.125654);
			C_Z_plot_E->SetTopMargin(0.0746073);
		h_Z_plot_E->Draw("");

		}
	
	
	
	sprintf(tstring_h,"Weighted Cutoff Rigidities for ISS Orbit");
	
		//	TH1F *h_Z_plot_R = new TH1F("h_Z_plot_R",tstring_h,100,0,20);
	TH1F *h_Z_plot_R = new TH1F("h_Z_plot_R","",100,0,20);
	
		// h_Z_plot_R
	h_Z_plot_R->GetXaxis()->SetTitle("Rigidity Cutoff at 450 km (GV)");
	h_Z_plot_R->GetXaxis()->CenterTitle();
	h_Z_plot_R->GetYaxis()->SetTitle("Orbit Fraction");
	h_Z_plot_R->GetYaxis()->CenterTitle();
	h_Z_plot_R->GetYaxis()->SetTitleOffset(1.7);
	h_Z_plot_R->GetYaxis()->SetLabelFont(nicefont);
	h_Z_plot_R->GetYaxis()->SetTitleFont(nicefont);
	h_Z_plot_R->GetXaxis()->SetLabelFont(nicefont);
	h_Z_plot_R->GetXaxis()->SetTitleFont(nicefont);
	h_Z_plot_R->SetTitleOffset(1.3);

	double weight_check=0.0;
	
	for(j=0;j<long_bins;j++)
	{
		for(k=0;k<lat_bins;k++)
		{
			if(dt[k] > 0.0)
			{
				h_Z_plot_R->Fill(alt_450[k].cut_off[j],dt[k]/(long_bins));
				weight_check+=dt[k]/(long_bins);
			}
		}
	}
	fprintf(stderr,"weight_check: %.4f\n",weight_check);
	
	
		// Fe Rigidity Cutoff for Threshold beween Geometric Area Regimes (500 MeV/nuc)
	
	
	double x_R_cut[2];
	double y_R_cut[2];
	x_R_cut[0] = el_dat[plot_Z-1].R_cut;
	x_R_cut[1] = el_dat[plot_Z-1].R_cut;
	y_R_cut[0] = -0.1;;
	y_R_cut[1] = 1.1;
	
	
	TGraph *gr_R_cut = new TGraph(2,x_R_cut,y_R_cut);
	
	gr_R_cut->SetMarkerStyle(1);
	gr_R_cut->SetMarkerSize(1);
	gr_R_cut->SetTitle("");
	gr_R_cut->GetXaxis()->SetTitle("GV");
	gr_R_cut->GetXaxis()->CenterTitle();
	gr_R_cut->GetYaxis()->SetTitle("counts");
	gr_R_cut->GetYaxis()->CenterTitle();
	gr_R_cut->GetXaxis()->SetTitleOffset(1.2);
	gr_R_cut->GetYaxis()->SetTitleOffset(1.4);
	gr_R_cut->GetYaxis()->SetLabelFont(nicefont);
	gr_R_cut->GetYaxis()->SetTitleFont(nicefont);
	gr_R_cut->GetXaxis()->SetLabelFont(nicefont);
	gr_R_cut->GetXaxis()->SetTitleFont(nicefont);
	gr_R_cut->GetXaxis()->SetRangeUser(0,5);
	gr_R_cut->GetYaxis()->SetRangeUser(0,0.4);
	
	
	TCanvas *C_Z_plot_R = new TCanvas("C_Z_plot_R","ISS Rigidity Thresholds",0,0,800,800);
	C_Z_plot_R->SetLeftMargin(0.148477);
	C_Z_plot_R->SetRightMargin(0.0507614);
	C_Z_plot_R->SetBottomMargin(0.125654);
	C_Z_plot_R->SetTopMargin(0.0746073);
	h_Z_plot_R->GetXaxis()->SetRangeUser(0,16);
	h_Z_plot_R->Draw("");
	gr_R_cut->SetLineWidth(2);
	gr_R_cut->SetLineColor(2);
	gr_R_cut->Draw("L");

	char tstring_h_R_cut[100];
	sprintf(tstring_h_R_cut,"Rigidity Threshold for _{%d}%s",plot_Z,chem_s);

	
	
	TLegend *leg_R_cut = new TLegend(0.3530,0.8290,0.6231,0.9197);
	leg_R_cut->AddEntry(gr_R_cut,tstring_h_R_cut,"l");	
	leg_R_cut->SetFillColor(10);
	leg_R_cut->SetTextFont(nicefont);
	leg_R_cut->SetBorderSize(0); // remove ugly border from stats 
	leg_R_cut->SetTextSize(0.03);

	
	leg_R_cut->Draw();
	
	double E_gf_thresh=0;
	double E_int=0;
	

		//E_gf_thresh=E_gf_thresh_Fe*pow((el_dat[4].Z/26.0),1.71);
	E_gf_thresh=E_gf_thresh_Fe;
	fprintf(stderr,"%2d E_gf_thresh_Fe: %.2f E_gf_thresh: %.2f\n",el_dat[4].Z,E_gf_thresh_Fe,E_gf_thresh);
	
	//  Actually find the modulated integral spectra!


	for(j=0;j<long_bins;j++)
		{
		for(k=0;k<lat_bins;k++)
			{
			el_dat[4].int_S_max[j][k]=gr_S_max_int_Z_5->Eval(el_dat[4].E_min_mag[j][k]);
			el_dat[4].int_S_min[j][k]=gr_S_min_int_Z_5->Eval(el_dat[4].E_min_mag[j][k]);
			// fprintf(stderr,"el_dat[4].int_S_max[%d][%d]: %e el_dat[4].int_S_min[%d][%d]: %e\n",j,k,el_dat[5].int_S_max[j][k],j,k,el_dat[5].int_S_min[j][k]);
			el_dat[5].int_S_max[j][k]=gr_S_max_int_Z_6->Eval(el_dat[5].E_min_mag[j][k]);
			el_dat[5].int_S_min[j][k]=gr_S_min_int_Z_6->Eval(el_dat[5].E_min_mag[j][k]);
			el_dat[6].int_S_max[j][k]=gr_S_max_int_Z_7->Eval(el_dat[6].E_min_mag[j][k]);
			el_dat[6].int_S_min[j][k]=gr_S_min_int_Z_7->Eval(el_dat[6].E_min_mag[j][k]);
			el_dat[7].int_S_max[j][k]=gr_S_max_int_Z_8->Eval(el_dat[7].E_min_mag[j][k]);
			el_dat[7].int_S_min[j][k]=gr_S_min_int_Z_8->Eval(el_dat[7].E_min_mag[j][k]);
			el_dat[8].int_S_max[j][k]=gr_S_max_int_Z_9->Eval(el_dat[8].E_min_mag[j][k]);
			el_dat[8].int_S_min[j][k]=gr_S_min_int_Z_9->Eval(el_dat[8].E_min_mag[j][k]);
			el_dat[9].int_S_max[j][k]=gr_S_max_int_Z_10->Eval(el_dat[9].E_min_mag[j][k]);
			el_dat[9].int_S_min[j][k]=gr_S_min_int_Z_10->Eval(el_dat[9].E_min_mag[j][k]);
			el_dat[10].int_S_max[j][k]=gr_S_max_int_Z_11->Eval(el_dat[10].E_min_mag[j][k]);
			el_dat[10].int_S_min[j][k]=gr_S_min_int_Z_11->Eval(el_dat[10].E_min_mag[j][k]);
			el_dat[11].int_S_max[j][k]=gr_S_max_int_Z_12->Eval(el_dat[11].E_min_mag[j][k]);
			el_dat[11].int_S_min[j][k]=gr_S_min_int_Z_12->Eval(el_dat[11].E_min_mag[j][k]);
			el_dat[12].int_S_max[j][k]=gr_S_max_int_Z_13->Eval(el_dat[12].E_min_mag[j][k]);
			el_dat[12].int_S_min[j][k]=gr_S_min_int_Z_13->Eval(el_dat[12].E_min_mag[j][k]);
			el_dat[13].int_S_max[j][k]=gr_S_max_int_Z_14->Eval(el_dat[13].E_min_mag[j][k]);
			el_dat[13].int_S_min[j][k]=gr_S_min_int_Z_14->Eval(el_dat[13].E_min_mag[j][k]);
			el_dat[14].int_S_max[j][k]=gr_S_max_int_Z_15->Eval(el_dat[14].E_min_mag[j][k]);
			el_dat[14].int_S_min[j][k]=gr_S_min_int_Z_15->Eval(el_dat[14].E_min_mag[j][k]);
			el_dat[15].int_S_max[j][k]=gr_S_max_int_Z_16->Eval(el_dat[15].E_min_mag[j][k]);
			el_dat[15].int_S_min[j][k]=gr_S_min_int_Z_16->Eval(el_dat[15].E_min_mag[j][k]);
			el_dat[16].int_S_max[j][k]=gr_S_max_int_Z_17->Eval(el_dat[16].E_min_mag[j][k]);
			el_dat[16].int_S_min[j][k]=gr_S_min_int_Z_17->Eval(el_dat[16].E_min_mag[j][k]);
			el_dat[17].int_S_max[j][k]=gr_S_max_int_Z_18->Eval(el_dat[17].E_min_mag[j][k]);
			el_dat[17].int_S_min[j][k]=gr_S_min_int_Z_18->Eval(el_dat[17].E_min_mag[j][k]);
			el_dat[18].int_S_max[j][k]=gr_S_max_int_Z_19->Eval(el_dat[18].E_min_mag[j][k]);
			el_dat[18].int_S_min[j][k]=gr_S_min_int_Z_19->Eval(el_dat[18].E_min_mag[j][k]);
			el_dat[19].int_S_max[j][k]=gr_S_max_int_Z_20->Eval(el_dat[19].E_min_mag[j][k]);
			el_dat[19].int_S_min[j][k]=gr_S_min_int_Z_20->Eval(el_dat[19].E_min_mag[j][k]);
			el_dat[20].int_S_max[j][k]=gr_S_max_int_Z_21->Eval(el_dat[20].E_min_mag[j][k]);
			el_dat[20].int_S_min[j][k]=gr_S_min_int_Z_21->Eval(el_dat[20].E_min_mag[j][k]);
			el_dat[21].int_S_max[j][k]=gr_S_max_int_Z_22->Eval(el_dat[21].E_min_mag[j][k]);
			el_dat[21].int_S_min[j][k]=gr_S_min_int_Z_22->Eval(el_dat[21].E_min_mag[j][k]);
			el_dat[22].int_S_max[j][k]=gr_S_max_int_Z_23->Eval(el_dat[22].E_min_mag[j][k]);
			el_dat[22].int_S_min[j][k]=gr_S_min_int_Z_23->Eval(el_dat[22].E_min_mag[j][k]);
			el_dat[23].int_S_max[j][k]=gr_S_max_int_Z_24->Eval(el_dat[23].E_min_mag[j][k]);
			el_dat[23].int_S_min[j][k]=gr_S_min_int_Z_24->Eval(el_dat[23].E_min_mag[j][k]);
			el_dat[24].int_S_max[j][k]=gr_S_max_int_Z_25->Eval(el_dat[24].E_min_mag[j][k]);
			el_dat[24].int_S_min[j][k]=gr_S_min_int_Z_25->Eval(el_dat[24].E_min_mag[j][k]);
			el_dat[25].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[25].E_min_mag[j][k]);
			el_dat[25].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[25].E_min_mag[j][k]);
			//	el_dat[26].int_S_max[j][k]=gr_S_max_int_Z_27->Eval(el_dat[26].E_min_mag[j][k]);
			//	el_dat[26].int_S_min[j][k]=gr_S_min_int_Z_27->Eval(el_dat[26].E_min_mag[j][k]);
			el_dat[26].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[26].E_min_mag[j][k])*el_dat[26].TOA_abund;
			el_dat[26].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[26].E_min_mag[j][k])*el_dat[26].TOA_abund;
			el_dat[27].int_S_max[j][k]=gr_S_max_int_Z_28->Eval(el_dat[27].E_min_mag[j][k]);
			el_dat[27].int_S_min[j][k]=gr_S_min_int_Z_28->Eval(el_dat[27].E_min_mag[j][k]);
			//	el_dat[28].int_S_max[j][k]=gr_S_max_int_Z_29->Eval(el_dat[28].E_min_mag[j][k]);
			//	el_dat[28].int_S_min[j][k]=gr_S_min_int_Z_29->Eval(el_dat[28].E_min_mag[j][k]);
			el_dat[28].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[28].E_min_mag[j][k])*el_dat[28].TOA_abund;
			el_dat[28].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[28].E_min_mag[j][k])*el_dat[28].TOA_abund;
			//	el_dat[29].int_S_max[j][k]=gr_S_max_int_Z_30->Eval(el_dat[29].E_min_mag[j][k]);
			//	el_dat[29].int_S_min[j][k]=gr_S_min_int_Z_30->Eval(el_dat[29].E_min_mag[j][k]);
			el_dat[29].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[29].E_min_mag[j][k])*el_dat[29].TOA_abund;
			el_dat[29].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[29].E_min_mag[j][k])*el_dat[29].TOA_abund;
			//	el_dat[30].int_S_max[j][k]=gr_S_max_int_Z_31->Eval(el_dat[30].E_min_mag[j][k]);
			//	el_dat[30].int_S_min[j][k]=gr_S_min_int_Z_31->Eval(el_dat[30].E_min_mag[j][k]);
			el_dat[30].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[30].E_min_mag[j][k])*el_dat[30].TOA_abund;
			el_dat[30].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[30].E_min_mag[j][k])*el_dat[30].TOA_abund;
			//	el_dat[31].int_S_max[j][k]=gr_S_max_int_Z_32->Eval(el_dat[31].E_min_mag[j][k]);
			//	el_dat[31].int_S_min[j][k]=gr_S_min_int_Z_32->Eval(el_dat[31].E_min_mag[j][k]);
			el_dat[31].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[31].E_min_mag[j][k])*el_dat[31].TOA_abund;
			el_dat[31].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[31].E_min_mag[j][k])*el_dat[31].TOA_abund;
			el_dat[32].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[32].E_min_mag[j][k])*el_dat[32].TOA_abund;
			el_dat[32].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[32].E_min_mag[j][k])*el_dat[32].TOA_abund;
			el_dat[33].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[33].E_min_mag[j][k])*el_dat[33].TOA_abund;
			el_dat[33].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[33].E_min_mag[j][k])*el_dat[33].TOA_abund;
			el_dat[34].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[34].E_min_mag[j][k])*el_dat[34].TOA_abund;
			el_dat[34].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[34].E_min_mag[j][k])*el_dat[34].TOA_abund;
			el_dat[35].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[35].E_min_mag[j][k])*el_dat[35].TOA_abund;
			el_dat[35].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[35].E_min_mag[j][k])*el_dat[35].TOA_abund;
			el_dat[36].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[36].E_min_mag[j][k])*el_dat[36].TOA_abund;
			el_dat[36].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[36].E_min_mag[j][k])*el_dat[36].TOA_abund;
			el_dat[37].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[37].E_min_mag[j][k])*el_dat[37].TOA_abund;
			el_dat[37].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[37].E_min_mag[j][k])*el_dat[37].TOA_abund;
			el_dat[38].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[38].E_min_mag[j][k])*el_dat[38].TOA_abund;
			el_dat[38].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[38].E_min_mag[j][k])*el_dat[38].TOA_abund;
			el_dat[39].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[39].E_min_mag[j][k])*el_dat[39].TOA_abund;
			el_dat[39].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[39].E_min_mag[j][k])*el_dat[39].TOA_abund;
			// Beyond TIGER
			el_dat[40].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[40].E_min_mag[j][k])*el_dat[40].TOA_abund;
			el_dat[40].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[40].E_min_mag[j][k])*el_dat[40].TOA_abund;
			el_dat[41].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[41].E_min_mag[j][k])*el_dat[41].TOA_abund;
			el_dat[41].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[41].E_min_mag[j][k])*el_dat[41].TOA_abund;
			el_dat[42].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[42].E_min_mag[j][k])*el_dat[42].TOA_abund;
			el_dat[42].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[42].E_min_mag[j][k])*el_dat[42].TOA_abund;
			el_dat[43].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[43].E_min_mag[j][k])*el_dat[43].TOA_abund;
			el_dat[43].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[43].E_min_mag[j][k])*el_dat[43].TOA_abund;
			el_dat[44].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[44].E_min_mag[j][k])*el_dat[44].TOA_abund;
			el_dat[44].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[44].E_min_mag[j][k])*el_dat[44].TOA_abund;
			el_dat[45].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[45].E_min_mag[j][k])*el_dat[45].TOA_abund;
			el_dat[45].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[45].E_min_mag[j][k])*el_dat[45].TOA_abund;
			el_dat[46].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[46].E_min_mag[j][k])*el_dat[46].TOA_abund;
			el_dat[46].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[46].E_min_mag[j][k])*el_dat[46].TOA_abund;
			el_dat[47].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[47].E_min_mag[j][k])*el_dat[47].TOA_abund;
			el_dat[47].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[47].E_min_mag[j][k])*el_dat[47].TOA_abund;
			el_dat[48].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[48].E_min_mag[j][k])*el_dat[48].TOA_abund;
			el_dat[48].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[48].E_min_mag[j][k])*el_dat[48].TOA_abund;
			el_dat[49].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[49].E_min_mag[j][k])*el_dat[49].TOA_abund;
			el_dat[49].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[49].E_min_mag[j][k])*el_dat[49].TOA_abund;
			el_dat[50].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[50].E_min_mag[j][k])*el_dat[50].TOA_abund;
			el_dat[50].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[50].E_min_mag[j][k])*el_dat[50].TOA_abund;
			el_dat[51].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[51].E_min_mag[j][k])*el_dat[51].TOA_abund;
			el_dat[51].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[51].E_min_mag[j][k])*el_dat[51].TOA_abund;
			el_dat[52].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[52].E_min_mag[j][k])*el_dat[52].TOA_abund;
			el_dat[52].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[52].E_min_mag[j][k])*el_dat[52].TOA_abund;
			el_dat[53].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[53].E_min_mag[j][k])*el_dat[53].TOA_abund;
			el_dat[53].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[53].E_min_mag[j][k])*el_dat[53].TOA_abund;
			el_dat[54].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[54].E_min_mag[j][k])*el_dat[54].TOA_abund;
			el_dat[54].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[54].E_min_mag[j][k])*el_dat[54].TOA_abund;
			el_dat[55].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[55].E_min_mag[j][k])*el_dat[55].TOA_abund;
			el_dat[55].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[55].E_min_mag[j][k])*el_dat[55].TOA_abund;
			el_dat[56].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[56].E_min_mag[j][k])*el_dat[56].TOA_abund;
			el_dat[56].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[56].E_min_mag[j][k])*el_dat[56].TOA_abund;
			el_dat[57].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[57].E_min_mag[j][k])*el_dat[57].TOA_abund;
			el_dat[57].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[57].E_min_mag[j][k])*el_dat[57].TOA_abund;
			el_dat[58].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[58].E_min_mag[j][k])*el_dat[58].TOA_abund;
			el_dat[58].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[58].E_min_mag[j][k])*el_dat[58].TOA_abund;
			el_dat[59].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[59].E_min_mag[j][k])*el_dat[59].TOA_abund;
			el_dat[59].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[59].E_min_mag[j][k])*el_dat[59].TOA_abund;
			el_dat[60].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[60].E_min_mag[j][k])*el_dat[60].TOA_abund;
			el_dat[60].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[60].E_min_mag[j][k])*el_dat[60].TOA_abund;
			el_dat[61].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[61].E_min_mag[j][k])*el_dat[61].TOA_abund;
			el_dat[61].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[61].E_min_mag[j][k])*el_dat[61].TOA_abund;
			el_dat[62].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[62].E_min_mag[j][k])*el_dat[62].TOA_abund;
			el_dat[62].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[62].E_min_mag[j][k])*el_dat[62].TOA_abund;
			el_dat[63].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[63].E_min_mag[j][k])*el_dat[63].TOA_abund;
			el_dat[63].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[63].E_min_mag[j][k])*el_dat[63].TOA_abund;
			el_dat[64].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[64].E_min_mag[j][k])*el_dat[64].TOA_abund;
			el_dat[64].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[64].E_min_mag[j][k])*el_dat[64].TOA_abund;
			el_dat[65].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[65].E_min_mag[j][k])*el_dat[65].TOA_abund;
			el_dat[65].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[65].E_min_mag[j][k])*el_dat[65].TOA_abund;
			el_dat[66].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[66].E_min_mag[j][k])*el_dat[66].TOA_abund;
			el_dat[66].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[66].E_min_mag[j][k])*el_dat[66].TOA_abund;
			el_dat[67].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[67].E_min_mag[j][k])*el_dat[67].TOA_abund;
			el_dat[67].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[67].E_min_mag[j][k])*el_dat[67].TOA_abund;
			el_dat[68].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[68].E_min_mag[j][k])*el_dat[68].TOA_abund;
			el_dat[68].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[68].E_min_mag[j][k])*el_dat[68].TOA_abund;
			el_dat[69].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[69].E_min_mag[j][k])*el_dat[69].TOA_abund;
			el_dat[69].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[69].E_min_mag[j][k])*el_dat[69].TOA_abund;
			el_dat[70].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[70].E_min_mag[j][k])*el_dat[70].TOA_abund;
			el_dat[70].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[70].E_min_mag[j][k])*el_dat[70].TOA_abund;
			el_dat[71].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[71].E_min_mag[j][k])*el_dat[71].TOA_abund;
			el_dat[71].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[71].E_min_mag[j][k])*el_dat[71].TOA_abund;
			el_dat[72].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[72].E_min_mag[j][k])*el_dat[72].TOA_abund;
			el_dat[72].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[72].E_min_mag[j][k])*el_dat[72].TOA_abund;
			el_dat[73].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[73].E_min_mag[j][k])*el_dat[73].TOA_abund;
			el_dat[73].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[73].E_min_mag[j][k])*el_dat[73].TOA_abund;
			el_dat[74].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[74].E_min_mag[j][k])*el_dat[74].TOA_abund;
			el_dat[74].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[74].E_min_mag[j][k])*el_dat[74].TOA_abund;
			el_dat[75].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[75].E_min_mag[j][k])*el_dat[75].TOA_abund;
			el_dat[75].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[75].E_min_mag[j][k])*el_dat[75].TOA_abund;
			el_dat[76].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[76].E_min_mag[j][k])*el_dat[76].TOA_abund;
			el_dat[76].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[76].E_min_mag[j][k])*el_dat[76].TOA_abund;
			el_dat[77].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[77].E_min_mag[j][k])*el_dat[77].TOA_abund;
			el_dat[77].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[77].E_min_mag[j][k])*el_dat[77].TOA_abund;
			el_dat[78].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[78].E_min_mag[j][k])*el_dat[78].TOA_abund;
			el_dat[78].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[78].E_min_mag[j][k])*el_dat[78].TOA_abund;
			el_dat[79].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[79].E_min_mag[j][k])*el_dat[79].TOA_abund;
			el_dat[79].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[79].E_min_mag[j][k])*el_dat[79].TOA_abund;
			el_dat[80].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[80].E_min_mag[j][k])*el_dat[80].TOA_abund;
			el_dat[80].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[80].E_min_mag[j][k])*el_dat[80].TOA_abund;
			el_dat[81].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[81].E_min_mag[j][k])*el_dat[81].TOA_abund;
			el_dat[81].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[81].E_min_mag[j][k])*el_dat[81].TOA_abund;
			el_dat[82].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[82].E_min_mag[j][k])*el_dat[82].TOA_abund;
			el_dat[82].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[82].E_min_mag[j][k])*el_dat[82].TOA_abund;
			el_dat[83].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[83].E_min_mag[j][k])*el_dat[83].TOA_abund;
			el_dat[83].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[83].E_min_mag[j][k])*el_dat[83].TOA_abund;
			el_dat[84].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[84].E_min_mag[j][k])*el_dat[84].TOA_abund;
			el_dat[84].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[84].E_min_mag[j][k])*el_dat[84].TOA_abund;
			el_dat[85].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[85].E_min_mag[j][k])*el_dat[85].TOA_abund;
			el_dat[85].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[85].E_min_mag[j][k])*el_dat[85].TOA_abund;
			el_dat[86].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[86].E_min_mag[j][k])*el_dat[86].TOA_abund;
			el_dat[86].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[86].E_min_mag[j][k])*el_dat[86].TOA_abund;
			el_dat[87].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[87].E_min_mag[j][k])*el_dat[87].TOA_abund;
			el_dat[87].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[87].E_min_mag[j][k])*el_dat[87].TOA_abund;
			el_dat[88].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[88].E_min_mag[j][k])*el_dat[88].TOA_abund;
			el_dat[88].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[88].E_min_mag[j][k])*el_dat[88].TOA_abund;
			el_dat[89].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[89].E_min_mag[j][k])*el_dat[89].TOA_abund;
			el_dat[89].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[89].E_min_mag[j][k])*el_dat[89].TOA_abund;
			el_dat[90].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[90].E_min_mag[j][k])*el_dat[90].TOA_abund;
			el_dat[90].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[90].E_min_mag[j][k])*el_dat[90].TOA_abund;
			el_dat[91].int_S_max[j][k]=gr_S_max_int_Z_26->Eval(el_dat[91].E_min_mag[j][k])*el_dat[91].TOA_abund;
			el_dat[91].int_S_min[j][k]=gr_S_min_int_Z_26->Eval(el_dat[91].E_min_mag[j][k])*el_dat[91].TOA_abund;
				
					// Geometry Factor Energy Threshold Based Calculations
					//				E_gf_thresh=E_gf_thresh_Fe*pow((el_dat[4].Z/26.0),1.71);
				E_gf_thresh=E_gf_thresh_Fe;
					//				fprintf(stderr,"%2d E_gf_thresh_Fe: %.2f E_gf_thresh: %.2f\n",el_dat[4].Z,E_gf_thresh_Fe,E_gf_thresh);
				
				
				if (el_dat[4].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[4].E_min_mag[j][k];
					el_dat[4].int_S_max_gf[j][k]=gr_S_max_int_Z_5->Eval(E_int);
					el_dat[4].int_S_min_gf[j][k]=gr_S_min_int_Z_5->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[4].int_S_max_gf[j][k]=gr_S_max_int_Z_5->Eval(E_int);
					el_dat[4].int_S_min_gf[j][k]=gr_S_min_int_Z_5->Eval(E_int);
					el_dat[4].int_S_max_gf_s[j][k]=gr_S_max_int_Z_5->Eval(el_dat[4].E_min_mag[j][k]) - gr_S_max_int_Z_5->Eval(E_int);
					el_dat[4].int_S_min_gf_s[j][k]=gr_S_min_int_Z_5->Eval(el_dat[4].E_min_mag[j][k]) - gr_S_min_int_Z_5->Eval(E_int);
				}
					//fprintf(stderr,"E_int: %.2f el_dat[4].E_min_mag[%d][%d]: %.2f E_gf_thresh: %.2f\n",E_int,j,k,el_dat[4].E_min_mag[j][k],E_gf_thresh);
				
				if (el_dat[5].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[5].E_min_mag[j][k];
					el_dat[5].int_S_max_gf[j][k]=gr_S_max_int_Z_6->Eval(E_int);
					el_dat[5].int_S_min_gf[j][k]=gr_S_min_int_Z_6->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[5].int_S_max_gf[j][k]=gr_S_max_int_Z_6->Eval(E_int);
					el_dat[5].int_S_min_gf[j][k]=gr_S_min_int_Z_6->Eval(E_int);
					el_dat[5].int_S_max_gf_s[j][k]=gr_S_max_int_Z_6->Eval(el_dat[5].E_min_mag[j][k]) - gr_S_max_int_Z_6->Eval(E_int);
					el_dat[5].int_S_min_gf_s[j][k]=gr_S_min_int_Z_6->Eval(el_dat[5].E_min_mag[j][k]) - gr_S_min_int_Z_6->Eval(E_int);
				}
				
				if (el_dat[6].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[6].E_min_mag[j][k];
					el_dat[6].int_S_max_gf[j][k]=gr_S_max_int_Z_7->Eval(E_int);
					el_dat[6].int_S_min_gf[j][k]=gr_S_min_int_Z_7->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[6].int_S_max_gf[j][k]=gr_S_max_int_Z_7->Eval(E_int);
					el_dat[6].int_S_min_gf[j][k]=gr_S_min_int_Z_7->Eval(E_int);
					el_dat[6].int_S_max_gf_s[j][k]=gr_S_max_int_Z_7->Eval(el_dat[6].E_min_mag[j][k]) - gr_S_max_int_Z_7->Eval(E_int);
					el_dat[6].int_S_min_gf_s[j][k]=gr_S_min_int_Z_7->Eval(el_dat[6].E_min_mag[j][k]) - gr_S_min_int_Z_7->Eval(E_int);
				}
				
				
				if (el_dat[7].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[7].E_min_mag[j][k];
					el_dat[7].int_S_max_gf[j][k]=gr_S_max_int_Z_8->Eval(E_int);
					el_dat[7].int_S_min_gf[j][k]=gr_S_min_int_Z_8->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[7].int_S_max_gf[j][k]=gr_S_max_int_Z_8->Eval(E_int);
					el_dat[7].int_S_min_gf[j][k]=gr_S_min_int_Z_8->Eval(E_int);
					el_dat[7].int_S_max_gf_s[j][k]=gr_S_max_int_Z_8->Eval(el_dat[7].E_min_mag[j][k]) - gr_S_max_int_Z_8->Eval(E_int);
					el_dat[7].int_S_min_gf_s[j][k]=gr_S_min_int_Z_8->Eval(el_dat[7].E_min_mag[j][k]) - gr_S_min_int_Z_8->Eval(E_int);
				}
				
				if (el_dat[8].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[8].E_min_mag[j][k];
					el_dat[8].int_S_max_gf[j][k]=gr_S_max_int_Z_9->Eval(E_int);
					el_dat[8].int_S_min_gf[j][k]=gr_S_min_int_Z_9->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[8].int_S_max_gf[j][k]=gr_S_max_int_Z_9->Eval(E_int);
					el_dat[8].int_S_min_gf[j][k]=gr_S_min_int_Z_9->Eval(E_int);
					el_dat[8].int_S_max_gf_s[j][k]=gr_S_max_int_Z_9->Eval(el_dat[8].E_min_mag[j][k]) - gr_S_max_int_Z_9->Eval(E_int);
					el_dat[8].int_S_min_gf_s[j][k]=gr_S_min_int_Z_9->Eval(el_dat[8].E_min_mag[j][k]) - gr_S_min_int_Z_9->Eval(E_int);
				}
				
				
				if (el_dat[9].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[9].E_min_mag[j][k];
					el_dat[9].int_S_max_gf[j][k]=gr_S_max_int_Z_10->Eval(E_int);
					el_dat[9].int_S_min_gf[j][k]=gr_S_min_int_Z_10->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[9].int_S_max_gf[j][k]=gr_S_max_int_Z_10->Eval(E_int);
					el_dat[9].int_S_min_gf[j][k]=gr_S_min_int_Z_10->Eval(E_int);
					el_dat[9].int_S_max_gf_s[j][k]=gr_S_max_int_Z_10->Eval(el_dat[9].E_min_mag[j][k]) - gr_S_max_int_Z_10->Eval(E_int);
					el_dat[9].int_S_min_gf_s[j][k]=gr_S_min_int_Z_10->Eval(el_dat[9].E_min_mag[j][k]) - gr_S_min_int_Z_10->Eval(E_int);
				}
				
				
				if (el_dat[10].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[10].E_min_mag[j][k];
					el_dat[10].int_S_max_gf[j][k]=gr_S_max_int_Z_11->Eval(E_int);
					el_dat[10].int_S_min_gf[j][k]=gr_S_min_int_Z_11->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[10].int_S_max_gf[j][k]=gr_S_max_int_Z_11->Eval(E_int);
					el_dat[10].int_S_min_gf[j][k]=gr_S_min_int_Z_11->Eval(E_int);
					el_dat[10].int_S_max_gf_s[j][k]=gr_S_max_int_Z_11->Eval(el_dat[10].E_min_mag[j][k]) - gr_S_max_int_Z_11->Eval(E_int);
					el_dat[10].int_S_min_gf_s[j][k]=gr_S_min_int_Z_11->Eval(el_dat[10].E_min_mag[j][k]) - gr_S_min_int_Z_11->Eval(E_int);
				}
				
				
				if (el_dat[11].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[11].E_min_mag[j][k];
					el_dat[11].int_S_max_gf[j][k]=gr_S_max_int_Z_12->Eval(E_int);
					el_dat[11].int_S_min_gf[j][k]=gr_S_min_int_Z_12->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[11].int_S_max_gf[j][k]=gr_S_max_int_Z_12->Eval(E_int);
					el_dat[11].int_S_min_gf[j][k]=gr_S_min_int_Z_12->Eval(E_int);
					el_dat[11].int_S_max_gf_s[j][k]=gr_S_max_int_Z_12->Eval(el_dat[11].E_min_mag[j][k]) - gr_S_max_int_Z_12->Eval(E_int);
					el_dat[11].int_S_min_gf_s[j][k]=gr_S_min_int_Z_12->Eval(el_dat[11].E_min_mag[j][k]) - gr_S_min_int_Z_12->Eval(E_int);
				}
				
				
				if (el_dat[12].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[12].E_min_mag[j][k];
					el_dat[12].int_S_max_gf[j][k]=gr_S_max_int_Z_13->Eval(E_int);
					el_dat[12].int_S_min_gf[j][k]=gr_S_min_int_Z_13->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[12].int_S_max_gf[j][k]=gr_S_max_int_Z_13->Eval(E_int);
					el_dat[12].int_S_min_gf[j][k]=gr_S_min_int_Z_13->Eval(E_int);
					el_dat[12].int_S_max_gf_s[j][k]=gr_S_max_int_Z_13->Eval(el_dat[12].E_min_mag[j][k]) - gr_S_max_int_Z_13->Eval(E_int);
					el_dat[12].int_S_min_gf_s[j][k]=gr_S_min_int_Z_13->Eval(el_dat[12].E_min_mag[j][k]) - gr_S_min_int_Z_13->Eval(E_int);
				}
				
				
				if (el_dat[13].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[13].E_min_mag[j][k];
					el_dat[13].int_S_max_gf[j][k]=gr_S_max_int_Z_14->Eval(E_int);
					el_dat[13].int_S_min_gf[j][k]=gr_S_min_int_Z_14->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[13].int_S_max_gf[j][k]=gr_S_max_int_Z_14->Eval(E_int);
					el_dat[13].int_S_min_gf[j][k]=gr_S_min_int_Z_14->Eval(E_int);
					el_dat[13].int_S_max_gf_s[j][k]=gr_S_max_int_Z_14->Eval(el_dat[13].E_min_mag[j][k]) - gr_S_max_int_Z_14->Eval(E_int);
					el_dat[13].int_S_min_gf_s[j][k]=gr_S_min_int_Z_14->Eval(el_dat[13].E_min_mag[j][k]) - gr_S_min_int_Z_14->Eval(E_int);
				}
				
				
				if (el_dat[14].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[14].E_min_mag[j][k];
					el_dat[14].int_S_max_gf[j][k]=gr_S_max_int_Z_15->Eval(E_int);
					el_dat[14].int_S_min_gf[j][k]=gr_S_min_int_Z_15->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[14].int_S_max_gf[j][k]=gr_S_max_int_Z_15->Eval(E_int);
					el_dat[14].int_S_min_gf[j][k]=gr_S_min_int_Z_15->Eval(E_int);
					el_dat[14].int_S_max_gf_s[j][k]=gr_S_max_int_Z_15->Eval(el_dat[14].E_min_mag[j][k]) - gr_S_max_int_Z_15->Eval(E_int);
					el_dat[14].int_S_min_gf_s[j][k]=gr_S_min_int_Z_15->Eval(el_dat[14].E_min_mag[j][k]) - gr_S_min_int_Z_15->Eval(E_int);
				}
				
				
				if (el_dat[15].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[15].E_min_mag[j][k];
					el_dat[15].int_S_max_gf[j][k]=gr_S_max_int_Z_16->Eval(E_int);
					el_dat[15].int_S_min_gf[j][k]=gr_S_min_int_Z_16->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[15].int_S_max_gf[j][k]=gr_S_max_int_Z_16->Eval(E_int);
					el_dat[15].int_S_min_gf[j][k]=gr_S_min_int_Z_16->Eval(E_int);
					el_dat[15].int_S_max_gf_s[j][k]=gr_S_max_int_Z_16->Eval(el_dat[15].E_min_mag[j][k]) - gr_S_max_int_Z_16->Eval(E_int);
					el_dat[15].int_S_min_gf_s[j][k]=gr_S_min_int_Z_16->Eval(el_dat[15].E_min_mag[j][k]) - gr_S_min_int_Z_16->Eval(E_int);
				}
				
				
				if (el_dat[16].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[16].E_min_mag[j][k];
					el_dat[16].int_S_max_gf[j][k]=gr_S_max_int_Z_17->Eval(E_int);
					el_dat[16].int_S_min_gf[j][k]=gr_S_min_int_Z_17->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[16].int_S_max_gf[j][k]=gr_S_max_int_Z_17->Eval(E_int);
					el_dat[16].int_S_min_gf[j][k]=gr_S_min_int_Z_17->Eval(E_int);
					el_dat[16].int_S_max_gf_s[j][k]=gr_S_max_int_Z_17->Eval(el_dat[16].E_min_mag[j][k]) - gr_S_max_int_Z_17->Eval(E_int);
					el_dat[16].int_S_min_gf_s[j][k]=gr_S_min_int_Z_17->Eval(el_dat[16].E_min_mag[j][k]) - gr_S_min_int_Z_17->Eval(E_int);
				}
				
				
				if (el_dat[17].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[17].E_min_mag[j][k];
					el_dat[17].int_S_max_gf[j][k]=gr_S_max_int_Z_18->Eval(E_int);
					el_dat[17].int_S_min_gf[j][k]=gr_S_min_int_Z_18->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[17].int_S_max_gf[j][k]=gr_S_max_int_Z_18->Eval(E_int);
					el_dat[17].int_S_min_gf[j][k]=gr_S_min_int_Z_18->Eval(E_int);
					el_dat[17].int_S_max_gf_s[j][k]=gr_S_max_int_Z_18->Eval(el_dat[17].E_min_mag[j][k]) - gr_S_max_int_Z_18->Eval(E_int);
					el_dat[17].int_S_min_gf_s[j][k]=gr_S_min_int_Z_18->Eval(el_dat[17].E_min_mag[j][k]) - gr_S_min_int_Z_18->Eval(E_int);
				}
				
				
				if (el_dat[18].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[18].E_min_mag[j][k];
					el_dat[18].int_S_max_gf[j][k]=gr_S_max_int_Z_19->Eval(E_int);
					el_dat[18].int_S_min_gf[j][k]=gr_S_min_int_Z_19->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[18].int_S_max_gf[j][k]=gr_S_max_int_Z_19->Eval(E_int);
					el_dat[18].int_S_min_gf[j][k]=gr_S_min_int_Z_19->Eval(E_int);
					el_dat[18].int_S_max_gf_s[j][k]=gr_S_max_int_Z_19->Eval(el_dat[18].E_min_mag[j][k]) - gr_S_max_int_Z_19->Eval(E_int);
					el_dat[18].int_S_min_gf_s[j][k]=gr_S_min_int_Z_19->Eval(el_dat[18].E_min_mag[j][k]) - gr_S_min_int_Z_19->Eval(E_int);
				}
				
				
				if (el_dat[19].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[19].E_min_mag[j][k];
					el_dat[19].int_S_max_gf[j][k]=gr_S_max_int_Z_20->Eval(E_int);
					el_dat[19].int_S_min_gf[j][k]=gr_S_min_int_Z_20->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[19].int_S_max_gf[j][k]=gr_S_max_int_Z_20->Eval(E_int);
					el_dat[19].int_S_min_gf[j][k]=gr_S_min_int_Z_20->Eval(E_int);
					el_dat[19].int_S_max_gf_s[j][k]=gr_S_max_int_Z_20->Eval(el_dat[19].E_min_mag[j][k]) - gr_S_max_int_Z_20->Eval(E_int);
					el_dat[19].int_S_min_gf_s[j][k]=gr_S_min_int_Z_20->Eval(el_dat[19].E_min_mag[j][k]) - gr_S_min_int_Z_20->Eval(E_int);
				}
				
				
				if (el_dat[20].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[20].E_min_mag[j][k];
					el_dat[20].int_S_max_gf[j][k]=gr_S_max_int_Z_21->Eval(E_int);
					el_dat[20].int_S_min_gf[j][k]=gr_S_min_int_Z_21->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[20].int_S_max_gf[j][k]=gr_S_max_int_Z_21->Eval(E_int);
					el_dat[20].int_S_min_gf[j][k]=gr_S_min_int_Z_21->Eval(E_int);
					el_dat[20].int_S_max_gf_s[j][k]=gr_S_max_int_Z_21->Eval(el_dat[20].E_min_mag[j][k]) - gr_S_max_int_Z_21->Eval(E_int);
					el_dat[20].int_S_min_gf_s[j][k]=gr_S_min_int_Z_21->Eval(el_dat[20].E_min_mag[j][k]) - gr_S_min_int_Z_21->Eval(E_int);
				}
				
				
				if (el_dat[21].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[21].E_min_mag[j][k];
					el_dat[21].int_S_max_gf[j][k]=gr_S_max_int_Z_22->Eval(E_int);
					el_dat[21].int_S_min_gf[j][k]=gr_S_min_int_Z_22->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[21].int_S_max_gf[j][k]=gr_S_max_int_Z_22->Eval(E_int);
					el_dat[21].int_S_min_gf[j][k]=gr_S_min_int_Z_22->Eval(E_int);
					el_dat[21].int_S_max_gf_s[j][k]=gr_S_max_int_Z_22->Eval(el_dat[21].E_min_mag[j][k]) - gr_S_max_int_Z_22->Eval(E_int);
					el_dat[21].int_S_min_gf_s[j][k]=gr_S_min_int_Z_22->Eval(el_dat[21].E_min_mag[j][k]) - gr_S_min_int_Z_22->Eval(E_int);
				}
				
				
				if (el_dat[22].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[22].E_min_mag[j][k];
					el_dat[22].int_S_max_gf[j][k]=gr_S_max_int_Z_23->Eval(E_int);
					el_dat[22].int_S_min_gf[j][k]=gr_S_min_int_Z_23->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[22].int_S_max_gf[j][k]=gr_S_max_int_Z_23->Eval(E_int);
					el_dat[22].int_S_min_gf[j][k]=gr_S_min_int_Z_23->Eval(E_int);
					el_dat[22].int_S_max_gf_s[j][k]=gr_S_max_int_Z_23->Eval(el_dat[22].E_min_mag[j][k]) - gr_S_max_int_Z_23->Eval(E_int);
					el_dat[22].int_S_min_gf_s[j][k]=gr_S_min_int_Z_23->Eval(el_dat[22].E_min_mag[j][k]) - gr_S_min_int_Z_23->Eval(E_int);
				}
				
				
				if (el_dat[23].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[23].E_min_mag[j][k];
					el_dat[23].int_S_max_gf[j][k]=gr_S_max_int_Z_24->Eval(E_int);
					el_dat[23].int_S_min_gf[j][k]=gr_S_min_int_Z_24->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[23].int_S_max_gf[j][k]=gr_S_max_int_Z_24->Eval(E_int);
					el_dat[23].int_S_min_gf[j][k]=gr_S_min_int_Z_24->Eval(E_int);
					el_dat[23].int_S_max_gf_s[j][k]=gr_S_max_int_Z_24->Eval(el_dat[23].E_min_mag[j][k]) - gr_S_max_int_Z_24->Eval(E_int);
					el_dat[23].int_S_min_gf_s[j][k]=gr_S_min_int_Z_24->Eval(el_dat[23].E_min_mag[j][k]) - gr_S_min_int_Z_24->Eval(E_int);
				}
				
				
				if (el_dat[24].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[24].E_min_mag[j][k];
					el_dat[24].int_S_max_gf[j][k]=gr_S_max_int_Z_25->Eval(E_int);
					el_dat[24].int_S_min_gf[j][k]=gr_S_min_int_Z_25->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[24].int_S_max_gf[j][k]=gr_S_max_int_Z_25->Eval(E_int);
					el_dat[24].int_S_min_gf[j][k]=gr_S_min_int_Z_25->Eval(E_int);
					el_dat[24].int_S_max_gf_s[j][k]=gr_S_max_int_Z_25->Eval(el_dat[24].E_min_mag[j][k]) - gr_S_max_int_Z_25->Eval(E_int);
					el_dat[24].int_S_min_gf_s[j][k]=gr_S_min_int_Z_25->Eval(el_dat[24].E_min_mag[j][k]) - gr_S_min_int_Z_25->Eval(E_int);
				}
				
				
				if (el_dat[25].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[25].E_min_mag[j][k];
					el_dat[25].int_S_max_gf[j][k]=gr_S_max_int_Z_26->Eval(E_int);
					el_dat[25].int_S_min_gf[j][k]=gr_S_min_int_Z_26->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[25].int_S_max_gf[j][k]=gr_S_max_int_Z_26->Eval(E_int);
					el_dat[25].int_S_min_gf[j][k]=gr_S_min_int_Z_26->Eval(E_int);
					el_dat[25].int_S_max_gf_s[j][k]=gr_S_max_int_Z_26->Eval(el_dat[25].E_min_mag[j][k]) - gr_S_max_int_Z_26->Eval(E_int);
					el_dat[25].int_S_min_gf_s[j][k]=gr_S_min_int_Z_26->Eval(el_dat[25].E_min_mag[j][k]) - gr_S_min_int_Z_26->Eval(E_int);
				}


				if (el_dat[26].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[26].E_min_mag[j][k];
					el_dat[26].int_S_max_gf[j][k]=gr_S_max_int_Z_26->Eval(E_int)*el_dat[26].TOA_abund;
					el_dat[26].int_S_min_gf[j][k]=gr_S_min_int_Z_26->Eval(E_int)*el_dat[26].TOA_abund;
				}
				else {
					E_int=E_gf_thresh;
					el_dat[26].int_S_max_gf[j][k]=gr_S_max_int_Z_26->Eval(E_int)*el_dat[26].TOA_abund;
					el_dat[26].int_S_min_gf[j][k]=gr_S_min_int_Z_26->Eval(E_int)*el_dat[26].TOA_abund;
					el_dat[26].int_S_max_gf_s[j][k]=(gr_S_max_int_Z_26->Eval(el_dat[26].E_min_mag[j][k]) - gr_S_max_int_Z_26->Eval(E_int))*el_dat[26].TOA_abund;
					el_dat[26].int_S_min_gf_s[j][k]=(gr_S_min_int_Z_26->Eval(el_dat[26].E_min_mag[j][k]) - gr_S_min_int_Z_26->Eval(E_int))*el_dat[26].TOA_abund;
				}
				
				if (el_dat[27].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[27].E_min_mag[j][k];
					el_dat[27].int_S_max_gf[j][k]=gr_S_max_int_Z_28->Eval(E_int);
					el_dat[27].int_S_min_gf[j][k]=gr_S_min_int_Z_28->Eval(E_int);
				}
				else {
					E_int=E_gf_thresh;
					el_dat[27].int_S_max_gf[j][k]=gr_S_max_int_Z_28->Eval(E_int);
					el_dat[27].int_S_min_gf[j][k]=gr_S_min_int_Z_28->Eval(E_int);
					el_dat[27].int_S_max_gf_s[j][k]=gr_S_max_int_Z_28->Eval(el_dat[27].E_min_mag[j][k]) - gr_S_max_int_Z_28->Eval(E_int);
					el_dat[27].int_S_min_gf_s[j][k]=gr_S_min_int_Z_28->Eval(el_dat[27].E_min_mag[j][k]) - gr_S_min_int_Z_28->Eval(E_int);
				}
				
				for (i=28; i<92; i++) {
					
				
				if (el_dat[i].E_min_mag[j][k] > E_gf_thresh) {
					E_int=el_dat[i].E_min_mag[j][k];
					el_dat[i].int_S_max_gf[j][k]=gr_S_max_int_Z_26->Eval(E_int)*el_dat[i].TOA_abund;
					el_dat[i].int_S_min_gf[j][k]=gr_S_min_int_Z_26->Eval(E_int)*el_dat[i].TOA_abund;
				}
				else {
					E_int=E_gf_thresh;
					el_dat[i].int_S_max_gf[j][k]=gr_S_max_int_Z_26->Eval(E_int)*el_dat[i].TOA_abund;
					el_dat[i].int_S_min_gf[j][k]=gr_S_min_int_Z_26->Eval(E_int)*el_dat[i].TOA_abund;
					el_dat[i].int_S_max_gf_s[j][k]=(gr_S_max_int_Z_26->Eval(el_dat[i].E_min_mag[j][k]) - gr_S_max_int_Z_26->Eval(E_int))*el_dat[i].TOA_abund;
					el_dat[i].int_S_min_gf_s[j][k]=(gr_S_min_int_Z_26->Eval(el_dat[i].E_min_mag[j][k]) - gr_S_min_int_Z_26->Eval(E_int))*el_dat[i].TOA_abund;
				}
				}
			}
		}

		for(i=0;i<el_index;i++)
			{
				el_dat[i].abund_S_max=0.0;
				el_dat[i].abund_S_min=0.0;
				el_dat[i].abund_S_max_gf=0.0;
				el_dat[i].abund_S_min_gf=0.0;
			el_dat[i].abund_S_max_tot=0.0;
			el_dat[i].abund_S_min_tot=0.0;
			el_dat[i].abund_S_max_no_mag=0.0;
			el_dat[i].abund_S_min_no_mag=0.0;
			}

	el_dat[4].abund_S_max_no_mag=gr_S_max_int_Z_5->Eval(S_max[4].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[4].abund_S_min_no_mag=gr_S_min_int_Z_5->Eval(S_min[4].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[5].abund_S_max_no_mag=gr_S_max_int_Z_6->Eval(S_max[5].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[5].abund_S_min_no_mag=gr_S_min_int_Z_6->Eval(S_min[5].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[6].abund_S_max_no_mag=gr_S_max_int_Z_7->Eval(S_max[6].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[6].abund_S_min_no_mag=gr_S_min_int_Z_7->Eval(S_min[6].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[7].abund_S_max_no_mag=gr_S_max_int_Z_8->Eval(S_max[7].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[7].abund_S_min_no_mag=gr_S_min_int_Z_8->Eval(S_min[7].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[8].abund_S_max_no_mag=gr_S_max_int_Z_9->Eval(S_max[8].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[8].abund_S_min_no_mag=gr_S_min_int_Z_9->Eval(S_min[8].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[9].abund_S_max_no_mag=gr_S_max_int_Z_10->Eval(S_max[9].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[9].abund_S_min_no_mag=gr_S_min_int_Z_10->Eval(S_min[9].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[10].abund_S_max_no_mag=gr_S_max_int_Z_11->Eval(S_max[10].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[10].abund_S_min_no_mag=gr_S_min_int_Z_11->Eval(S_min[10].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[11].abund_S_max_no_mag=gr_S_max_int_Z_12->Eval(S_max[11].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[11].abund_S_min_no_mag=gr_S_min_int_Z_12->Eval(S_min[11].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[12].abund_S_max_no_mag=gr_S_max_int_Z_13->Eval(S_max[12].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[12].abund_S_min_no_mag=gr_S_min_int_Z_13->Eval(S_min[12].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[13].abund_S_max_no_mag=gr_S_max_int_Z_14->Eval(S_max[13].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[13].abund_S_min_no_mag=gr_S_min_int_Z_14->Eval(S_min[13].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[14].abund_S_max_no_mag=gr_S_max_int_Z_15->Eval(S_max[14].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[14].abund_S_min_no_mag=gr_S_min_int_Z_15->Eval(S_min[14].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[15].abund_S_max_no_mag=gr_S_max_int_Z_16->Eval(S_max[15].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[15].abund_S_min_no_mag=gr_S_min_int_Z_16->Eval(S_min[15].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[16].abund_S_max_no_mag=gr_S_max_int_Z_17->Eval(S_max[16].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[16].abund_S_min_no_mag=gr_S_min_int_Z_17->Eval(S_min[16].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[17].abund_S_max_no_mag=gr_S_max_int_Z_18->Eval(S_max[17].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[17].abund_S_min_no_mag=gr_S_min_int_Z_18->Eval(S_min[17].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[18].abund_S_max_no_mag=gr_S_max_int_Z_19->Eval(S_max[18].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[18].abund_S_min_no_mag=gr_S_min_int_Z_19->Eval(S_min[18].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[19].abund_S_max_no_mag=gr_S_max_int_Z_20->Eval(S_max[19].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[19].abund_S_min_no_mag=gr_S_min_int_Z_20->Eval(S_min[19].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[20].abund_S_max_no_mag=gr_S_max_int_Z_21->Eval(S_max[20].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[20].abund_S_min_no_mag=gr_S_min_int_Z_21->Eval(S_min[20].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[21].abund_S_max_no_mag=gr_S_max_int_Z_22->Eval(S_max[21].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[21].abund_S_min_no_mag=gr_S_min_int_Z_22->Eval(S_min[21].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[22].abund_S_max_no_mag=gr_S_max_int_Z_23->Eval(S_max[22].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[22].abund_S_min_no_mag=gr_S_min_int_Z_23->Eval(S_min[22].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[23].abund_S_max_no_mag=gr_S_max_int_Z_24->Eval(S_max[23].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[23].abund_S_min_no_mag=gr_S_min_int_Z_24->Eval(S_min[23].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[24].abund_S_max_no_mag=gr_S_max_int_Z_25->Eval(S_max[24].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[24].abund_S_min_no_mag=gr_S_min_int_Z_25->Eval(S_min[24].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[25].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[25].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[25].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[25].EI[0])*geo_factor_cm*time_in_orbit;
//	el_dat[26].abund_S_max_no_mag=gr_S_max_int_Z_27->Eval(S_max[26].EI[0])*geo_factor_cm*time_in_orbit;
//	el_dat[26].abund_S_min_no_mag=gr_S_min_int_Z_27->Eval(S_min[26].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[26].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[26].EI[0])*geo_factor_cm*time_in_orbit*el_dat[26].TOA_abund;
	el_dat[26].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[26].EI[0])*geo_factor_cm*time_in_orbit*el_dat[26].TOA_abund;
	el_dat[27].abund_S_max_no_mag=gr_S_max_int_Z_27->Eval(S_max[27].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[27].abund_S_min_no_mag=gr_S_min_int_Z_27->Eval(S_min[27].EI[0])*geo_factor_cm*time_in_orbit;
//	el_dat[28].abund_S_max_no_mag=gr_S_max_int_Z_29->Eval(S_max[28].EI[0])*geo_factor_cm*time_in_orbit;
//	el_dat[28].abund_S_min_no_mag=gr_S_min_int_Z_29->Eval(S_min[28].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[28].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[28].EI[0])*geo_factor_cm*time_in_orbit*el_dat[28].TOA_abund;
	el_dat[28].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[28].EI[0])*geo_factor_cm*time_in_orbit*el_dat[28].TOA_abund;
//	el_dat[29].abund_S_max_no_mag=gr_S_max_int_Z_30->Eval(S_max[29].EI[0])*geo_factor_cm*time_in_orbit;
//	el_dat[29].abund_S_min_no_mag=gr_S_min_int_Z_30->Eval(S_min[29].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[29].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[29].EI[0])*geo_factor_cm*time_in_orbit*el_dat[29].TOA_abund;
	el_dat[29].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[29].EI[0])*geo_factor_cm*time_in_orbit*el_dat[29].TOA_abund;
//	el_dat[30].abund_S_max_no_mag=gr_S_max_int_Z_31->Eval(S_max[30].EI[0])*geo_factor_cm*time_in_orbit;
//	el_dat[30].abund_S_min_no_mag=gr_S_min_int_Z_31->Eval(S_min[30].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[30].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[30].EI[0])*geo_factor_cm*time_in_orbit*el_dat[30].TOA_abund;
	el_dat[30].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[30].EI[0])*geo_factor_cm*time_in_orbit*el_dat[30].TOA_abund;
//	el_dat[31].abund_S_max_no_mag=gr_S_max_int_Z_32->Eval(S_max[31].EI[0])*geo_factor_cm*time_in_orbit;
//	el_dat[31].abund_S_min_no_mag=gr_S_min_int_Z_32->Eval(S_min[31].EI[0])*geo_factor_cm*time_in_orbit;
	el_dat[31].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[31].EI[0])*geo_factor_cm*time_in_orbit*el_dat[31].TOA_abund;
	el_dat[31].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[31].EI[0])*geo_factor_cm*time_in_orbit*el_dat[31].TOA_abund;
	el_dat[32].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[32].EI[0])*geo_factor_cm*time_in_orbit*el_dat[32].TOA_abund;
	el_dat[32].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[32].EI[0])*geo_factor_cm*time_in_orbit*el_dat[32].TOA_abund;
	el_dat[33].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[33].EI[0])*geo_factor_cm*time_in_orbit*el_dat[33].TOA_abund;
	el_dat[33].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[33].EI[0])*geo_factor_cm*time_in_orbit*el_dat[33].TOA_abund;
	el_dat[34].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[34].EI[0])*geo_factor_cm*time_in_orbit*el_dat[34].TOA_abund;
	el_dat[34].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[34].EI[0])*geo_factor_cm*time_in_orbit*el_dat[34].TOA_abund;
	el_dat[35].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[35].EI[0])*geo_factor_cm*time_in_orbit*el_dat[35].TOA_abund;
	el_dat[35].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[35].EI[0])*geo_factor_cm*time_in_orbit*el_dat[35].TOA_abund;
	el_dat[36].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[36].EI[0])*geo_factor_cm*time_in_orbit*el_dat[36].TOA_abund;
	el_dat[36].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[36].EI[0])*geo_factor_cm*time_in_orbit*el_dat[36].TOA_abund;
	el_dat[37].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[37].EI[0])*geo_factor_cm*time_in_orbit*el_dat[37].TOA_abund;
	el_dat[37].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[37].EI[0])*geo_factor_cm*time_in_orbit*el_dat[37].TOA_abund;
	el_dat[38].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[38].EI[0])*geo_factor_cm*time_in_orbit*el_dat[38].TOA_abund;
	el_dat[38].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[38].EI[0])*geo_factor_cm*time_in_orbit*el_dat[38].TOA_abund;
	el_dat[39].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[39].EI[0])*geo_factor_cm*time_in_orbit*el_dat[39].TOA_abund;
	el_dat[39].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[39].EI[0])*geo_factor_cm*time_in_orbit*el_dat[39].TOA_abund;
	// Beyond TIGER
	el_dat[40].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[40].EI[0])*geo_factor_cm*time_in_orbit*el_dat[40].TOA_abund;
	el_dat[40].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[40].EI[0])*geo_factor_cm*time_in_orbit*el_dat[40].TOA_abund;
	el_dat[41].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[41].EI[0])*geo_factor_cm*time_in_orbit*el_dat[41].TOA_abund;
	el_dat[41].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[41].EI[0])*geo_factor_cm*time_in_orbit*el_dat[41].TOA_abund;
	el_dat[42].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[42].EI[0])*geo_factor_cm*time_in_orbit*el_dat[42].TOA_abund;
	el_dat[42].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[42].EI[0])*geo_factor_cm*time_in_orbit*el_dat[42].TOA_abund;
	el_dat[43].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[43].EI[0])*geo_factor_cm*time_in_orbit*el_dat[43].TOA_abund;
	el_dat[43].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[43].EI[0])*geo_factor_cm*time_in_orbit*el_dat[43].TOA_abund;
	el_dat[44].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[44].EI[0])*geo_factor_cm*time_in_orbit*el_dat[44].TOA_abund;
	el_dat[44].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[44].EI[0])*geo_factor_cm*time_in_orbit*el_dat[44].TOA_abund;
	el_dat[45].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[45].EI[0])*geo_factor_cm*time_in_orbit*el_dat[45].TOA_abund;
	el_dat[45].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[45].EI[0])*geo_factor_cm*time_in_orbit*el_dat[45].TOA_abund;
	el_dat[46].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[46].EI[0])*geo_factor_cm*time_in_orbit*el_dat[46].TOA_abund;
	el_dat[46].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[46].EI[0])*geo_factor_cm*time_in_orbit*el_dat[46].TOA_abund;
	el_dat[47].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[47].EI[0])*geo_factor_cm*time_in_orbit*el_dat[47].TOA_abund;
	el_dat[47].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[47].EI[0])*geo_factor_cm*time_in_orbit*el_dat[47].TOA_abund;
	el_dat[48].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[48].EI[0])*geo_factor_cm*time_in_orbit*el_dat[48].TOA_abund;
	el_dat[48].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[48].EI[0])*geo_factor_cm*time_in_orbit*el_dat[48].TOA_abund;
	el_dat[49].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[49].EI[0])*geo_factor_cm*time_in_orbit*el_dat[49].TOA_abund;
	el_dat[49].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[49].EI[0])*geo_factor_cm*time_in_orbit*el_dat[49].TOA_abund;
	el_dat[50].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[50].EI[0])*geo_factor_cm*time_in_orbit*el_dat[50].TOA_abund;
	el_dat[50].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[50].EI[0])*geo_factor_cm*time_in_orbit*el_dat[50].TOA_abund;
	el_dat[51].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[51].EI[0])*geo_factor_cm*time_in_orbit*el_dat[51].TOA_abund;
	el_dat[51].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[51].EI[0])*geo_factor_cm*time_in_orbit*el_dat[51].TOA_abund;
	el_dat[52].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[52].EI[0])*geo_factor_cm*time_in_orbit*el_dat[52].TOA_abund;
	el_dat[52].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[52].EI[0])*geo_factor_cm*time_in_orbit*el_dat[52].TOA_abund;
	el_dat[53].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[53].EI[0])*geo_factor_cm*time_in_orbit*el_dat[53].TOA_abund;
	el_dat[53].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[53].EI[0])*geo_factor_cm*time_in_orbit*el_dat[53].TOA_abund;
	el_dat[54].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[54].EI[0])*geo_factor_cm*time_in_orbit*el_dat[54].TOA_abund;
	el_dat[54].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[54].EI[0])*geo_factor_cm*time_in_orbit*el_dat[54].TOA_abund;
	el_dat[55].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[55].EI[0])*geo_factor_cm*time_in_orbit*el_dat[55].TOA_abund;
	el_dat[55].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[55].EI[0])*geo_factor_cm*time_in_orbit*el_dat[55].TOA_abund;
	el_dat[56].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[56].EI[0])*geo_factor_cm*time_in_orbit*el_dat[56].TOA_abund;
	el_dat[56].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[56].EI[0])*geo_factor_cm*time_in_orbit*el_dat[56].TOA_abund;
	el_dat[57].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[57].EI[0])*geo_factor_cm*time_in_orbit*el_dat[57].TOA_abund;
	el_dat[57].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[57].EI[0])*geo_factor_cm*time_in_orbit*el_dat[57].TOA_abund;
	el_dat[58].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[58].EI[0])*geo_factor_cm*time_in_orbit*el_dat[58].TOA_abund;
	el_dat[58].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[58].EI[0])*geo_factor_cm*time_in_orbit*el_dat[58].TOA_abund;
	el_dat[59].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[59].EI[0])*geo_factor_cm*time_in_orbit*el_dat[59].TOA_abund;
	el_dat[59].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[59].EI[0])*geo_factor_cm*time_in_orbit*el_dat[59].TOA_abund;
	el_dat[60].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[60].EI[0])*geo_factor_cm*time_in_orbit*el_dat[60].TOA_abund;
	el_dat[60].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[60].EI[0])*geo_factor_cm*time_in_orbit*el_dat[60].TOA_abund;
	el_dat[61].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[61].EI[0])*geo_factor_cm*time_in_orbit*el_dat[61].TOA_abund;
	el_dat[61].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[61].EI[0])*geo_factor_cm*time_in_orbit*el_dat[61].TOA_abund;
	el_dat[62].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[62].EI[0])*geo_factor_cm*time_in_orbit*el_dat[62].TOA_abund;
	el_dat[62].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[62].EI[0])*geo_factor_cm*time_in_orbit*el_dat[62].TOA_abund;
	el_dat[63].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[63].EI[0])*geo_factor_cm*time_in_orbit*el_dat[63].TOA_abund;
	el_dat[63].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[63].EI[0])*geo_factor_cm*time_in_orbit*el_dat[63].TOA_abund;
	el_dat[64].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[64].EI[0])*geo_factor_cm*time_in_orbit*el_dat[64].TOA_abund;
	el_dat[64].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[64].EI[0])*geo_factor_cm*time_in_orbit*el_dat[64].TOA_abund;
	el_dat[65].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[65].EI[0])*geo_factor_cm*time_in_orbit*el_dat[65].TOA_abund;
	el_dat[65].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[65].EI[0])*geo_factor_cm*time_in_orbit*el_dat[65].TOA_abund;
	el_dat[66].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[66].EI[0])*geo_factor_cm*time_in_orbit*el_dat[66].TOA_abund;
	el_dat[66].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[66].EI[0])*geo_factor_cm*time_in_orbit*el_dat[66].TOA_abund;
	el_dat[67].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[67].EI[0])*geo_factor_cm*time_in_orbit*el_dat[67].TOA_abund;
	el_dat[67].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[67].EI[0])*geo_factor_cm*time_in_orbit*el_dat[67].TOA_abund;
	el_dat[68].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[68].EI[0])*geo_factor_cm*time_in_orbit*el_dat[68].TOA_abund;
	el_dat[68].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[68].EI[0])*geo_factor_cm*time_in_orbit*el_dat[68].TOA_abund;
	el_dat[69].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[69].EI[0])*geo_factor_cm*time_in_orbit*el_dat[69].TOA_abund;
	el_dat[69].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[69].EI[0])*geo_factor_cm*time_in_orbit*el_dat[69].TOA_abund;
	el_dat[70].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[70].EI[0])*geo_factor_cm*time_in_orbit*el_dat[70].TOA_abund;
	el_dat[70].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[70].EI[0])*geo_factor_cm*time_in_orbit*el_dat[70].TOA_abund;
	el_dat[71].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[71].EI[0])*geo_factor_cm*time_in_orbit*el_dat[71].TOA_abund;
	el_dat[71].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[71].EI[0])*geo_factor_cm*time_in_orbit*el_dat[71].TOA_abund;
	el_dat[72].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[72].EI[0])*geo_factor_cm*time_in_orbit*el_dat[72].TOA_abund;
	el_dat[72].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[72].EI[0])*geo_factor_cm*time_in_orbit*el_dat[72].TOA_abund;
	el_dat[73].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[73].EI[0])*geo_factor_cm*time_in_orbit*el_dat[73].TOA_abund;
	el_dat[73].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[73].EI[0])*geo_factor_cm*time_in_orbit*el_dat[73].TOA_abund;
	el_dat[74].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[74].EI[0])*geo_factor_cm*time_in_orbit*el_dat[74].TOA_abund;
	el_dat[74].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[74].EI[0])*geo_factor_cm*time_in_orbit*el_dat[74].TOA_abund;
	el_dat[75].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[75].EI[0])*geo_factor_cm*time_in_orbit*el_dat[75].TOA_abund;
	el_dat[75].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[75].EI[0])*geo_factor_cm*time_in_orbit*el_dat[75].TOA_abund;
	el_dat[76].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[76].EI[0])*geo_factor_cm*time_in_orbit*el_dat[76].TOA_abund;
	el_dat[76].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[76].EI[0])*geo_factor_cm*time_in_orbit*el_dat[76].TOA_abund;
	el_dat[77].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[77].EI[0])*geo_factor_cm*time_in_orbit*el_dat[77].TOA_abund;
	el_dat[77].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[77].EI[0])*geo_factor_cm*time_in_orbit*el_dat[77].TOA_abund;
	el_dat[78].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[78].EI[0])*geo_factor_cm*time_in_orbit*el_dat[78].TOA_abund;
	el_dat[78].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[78].EI[0])*geo_factor_cm*time_in_orbit*el_dat[78].TOA_abund;
	el_dat[79].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[79].EI[0])*geo_factor_cm*time_in_orbit*el_dat[79].TOA_abund;
	el_dat[79].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[79].EI[0])*geo_factor_cm*time_in_orbit*el_dat[79].TOA_abund;
	el_dat[80].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[80].EI[0])*geo_factor_cm*time_in_orbit*el_dat[80].TOA_abund;
	el_dat[80].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[80].EI[0])*geo_factor_cm*time_in_orbit*el_dat[80].TOA_abund;
	el_dat[81].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[81].EI[0])*geo_factor_cm*time_in_orbit*el_dat[81].TOA_abund;
	el_dat[81].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[81].EI[0])*geo_factor_cm*time_in_orbit*el_dat[81].TOA_abund;
	el_dat[82].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[82].EI[0])*geo_factor_cm*time_in_orbit*el_dat[82].TOA_abund;
	el_dat[82].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[82].EI[0])*geo_factor_cm*time_in_orbit*el_dat[82].TOA_abund;
	el_dat[83].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[83].EI[0])*geo_factor_cm*time_in_orbit*el_dat[83].TOA_abund;
	el_dat[83].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[83].EI[0])*geo_factor_cm*time_in_orbit*el_dat[83].TOA_abund;
	el_dat[84].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[84].EI[0])*geo_factor_cm*time_in_orbit*el_dat[84].TOA_abund;
	el_dat[84].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[84].EI[0])*geo_factor_cm*time_in_orbit*el_dat[84].TOA_abund;
	el_dat[85].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[85].EI[0])*geo_factor_cm*time_in_orbit*el_dat[85].TOA_abund;
	el_dat[85].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[85].EI[0])*geo_factor_cm*time_in_orbit*el_dat[85].TOA_abund;
	el_dat[86].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[86].EI[0])*geo_factor_cm*time_in_orbit*el_dat[86].TOA_abund;
	el_dat[86].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[86].EI[0])*geo_factor_cm*time_in_orbit*el_dat[86].TOA_abund;
	el_dat[87].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[87].EI[0])*geo_factor_cm*time_in_orbit*el_dat[87].TOA_abund;
	el_dat[87].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[87].EI[0])*geo_factor_cm*time_in_orbit*el_dat[87].TOA_abund;
	el_dat[88].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[88].EI[0])*geo_factor_cm*time_in_orbit*el_dat[88].TOA_abund;
	el_dat[88].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[88].EI[0])*geo_factor_cm*time_in_orbit*el_dat[88].TOA_abund;
	el_dat[89].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[89].EI[0])*geo_factor_cm*time_in_orbit*el_dat[89].TOA_abund;
	el_dat[89].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[89].EI[0])*geo_factor_cm*time_in_orbit*el_dat[89].TOA_abund;
	el_dat[90].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[90].EI[0])*geo_factor_cm*time_in_orbit*el_dat[90].TOA_abund;
	el_dat[90].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[90].EI[0])*geo_factor_cm*time_in_orbit*el_dat[90].TOA_abund;
	el_dat[91].abund_S_max_no_mag=gr_S_max_int_Z_26->Eval(S_max[91].EI[0])*geo_factor_cm*time_in_orbit*el_dat[91].TOA_abund;
	el_dat[91].abund_S_min_no_mag=gr_S_min_int_Z_26->Eval(S_min[91].EI[0])*geo_factor_cm*time_in_orbit*el_dat[91].TOA_abund;

		// Differential Geometry Factor interlude:
	
	char geom_file[200]="/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/TIGERISS_Abundances/angle/dgf_1d_bin1.txt";

	int dgf_index=0;
	int dgf_count=1;
	
	double angle_dgf[20];
	double geo_factor_dgf[20];
	double geo_factor_dgf_w[20];
	
	angle_dgf[0]=0.0;
	geo_factor_dgf[0]=0.0;
	geo_factor_dgf_w[0]=0.0;
	
	in_1.open(geom_file);
	while(1) {
		in_1 >> angle_dgf[dgf_count] >> geo_factor_dgf[dgf_count];
		if(in_1.eof() || !in_1.good()) break;
		geo_factor_dgf_w[dgf_count]=geo_factor_dgf[dgf_count];
		
		if (angle_dgf[dgf_count]>=45.0) {
			geo_factor_dgf_w[dgf_count]=(3.0/4.0)*geo_factor_dgf[dgf_count];
		}
		dgf_count++;
	}
	in_1.close();
	in_1.clear();
	
	TGraph *gr_dgf_w = new TGraph(dgf_count,angle_dgf,geo_factor_dgf_w);

	TGraph *gr_dgf = new TGraph(dgf_count,angle_dgf,geo_factor_dgf);
	
	gr_dgf->SetMarkerStyle(1);
	gr_dgf->SetMarkerSize(1);
		//	gr_dgf->SetTitle("CALET Differential Geometry Factor");
	gr_dgf->SetTitle("");
	gr_dgf->GetXaxis()->SetTitle("Incidence Angle (degrees)");
	gr_dgf->GetXaxis()->CenterTitle();
	gr_dgf->GetYaxis()->SetTitle("Geometry Factor (cm^{2}sr)");
	gr_dgf->GetYaxis()->CenterTitle();
	gr_dgf->GetXaxis()->SetTitleOffset(1.2);
	gr_dgf->GetYaxis()->SetTitleOffset(1.3);
	gr_dgf->GetYaxis()->SetLabelFont(nicefont);
	gr_dgf->GetYaxis()->SetTitleFont(nicefont);
	gr_dgf->GetXaxis()->SetLabelFont(nicefont);
	gr_dgf->GetXaxis()->SetTitleFont(nicefont);
	gr_dgf->GetXaxis()->SetRangeUser(0,90);
	gr_dgf->GetYaxis()->SetRangeUser(0,450);
	
	
	
	TCanvas *C_dgf = new TCanvas("C_dgf","CALET Differental Geometry Factor",0,0,800,800);
	C_dgf->SetLeftMargin(0.139594);
	C_dgf->SetRightMargin(0.0596447);
	C_dgf->SetBottomMargin(0.11911);
	C_dgf->SetTopMargin(0.0811518);

	gr_dgf->SetLineWidth(2);
	gr_dgf->Draw("AL");
	gr_dgf_w->SetLineWidth(2);
	gr_dgf_w->SetLineColor(2);
	gr_dgf_w->Draw("L");
	
	C_dgf->SetGrid();

	
	TLegend *leg_dgf = new TLegend(0.251,0.229,0.501,0.380);
	leg_dgf->AddEntry(gr_dgf,"CHD and 4 IMC layers","l");
	leg_dgf->AddEntry(gr_dgf_w,"45^{o} limit on one side","l");	
	
	leg_dgf->SetFillColor(10);
	leg_dgf->SetTextFont(nicefont);
	leg_dgf->SetBorderSize(0); // remove ugly border from stats 
	leg_dgf->SetTextSize(0.04);
	
	leg_dgf->Draw();
	

	double geo_factor_calc=0.0;
	
	for(i=0;i<el_index;i++)
		{
		for(j=0;j<long_bins_calc;j++)
			{
			for(k=0;k<lat_bins;k++)
				{
					el_dat[i].abund_S_max+=el_dat[i].int_S_max[j][k]*dt[k]*geo_factor_cm*time_in_orbit/long_bins_calc;
					el_dat[i].abund_S_min+=el_dat[i].int_S_min[j][k]*dt[k]*geo_factor_cm*time_in_orbit/long_bins_calc;
					
					el_dat[i].abund_S_max_gf+=el_dat[i].int_S_max_gf[j][k]*dt[k]*geo_factor_cm*time_in_orbit/long_bins_calc;
					el_dat[i].abund_S_min_gf+=el_dat[i].int_S_min_gf[j][k]*dt[k]*geo_factor_cm*time_in_orbit/long_bins_calc;
					
					for (dgf_index=0; dgf_index<dgf_count; dgf_index++) {
						
						geo_factor_calc=geo_factor_dgf[dgf_index];

						if (angle_dgf[dgf_index]>=45.0) {
							geo_factor_calc=(3.0/4.0)*geo_factor_dgf[dgf_index];
						}
						if (geo_factor_calc != geo_factor_dgf_w[dgf_index]) {
							fprintf(stderr,"oops!!!!!!!!!!!!!!!!!!!\n");
						}
						
						el_dat[i].abund_S_max_dgf+=(el_dat[i].int_S_max_gf[j][k]*dt[k]*time_in_orbit/long_bins_calc)*
						exp(-width_PVT/(el_dat[i].lambda_PVT*cos(TMath::Pi()*angle_dgf[dgf_index]/180.0)))*geo_factor_calc;
						
						el_dat[i].abund_S_min_dgf+=(el_dat[i].int_S_min_gf[j][k]*dt[k]*time_in_orbit/long_bins_calc)*
						exp(-width_PVT/(el_dat[i].lambda_PVT*cos(TMath::Pi()*angle_dgf[dgf_index]/180.0)))*geo_factor_calc;
						
						el_dat[i].abund_S_max_dgf_noic+=(el_dat[i].int_S_max_gf[j][k]*dt[k]*time_in_orbit/long_bins_calc)*geo_factor_calc;
						
						el_dat[i].abund_S_min_dgf_noic+=(el_dat[i].int_S_min_gf[j][k]*dt[k]*time_in_orbit/long_bins_calc)*geo_factor_calc;
						
					}
					
					if (el_dat[i].E_min_mag[j][k] < E_gf_thresh) {
						el_dat[i].abund_S_max_gf+=el_dat[i].int_S_max_gf_s[j][k]*dt[k]*geo_factor_cm_s*time_in_orbit/long_bins_calc;
						el_dat[i].abund_S_min_gf+=el_dat[i].int_S_min_gf_s[j][k]*dt[k]*geo_factor_cm_s*time_in_orbit/long_bins_calc;
							// The lower energy events requiring passage through the TASC use the full instrument geometry factor and average interaction corrections
						el_dat[i].abund_S_max_dgf+=(el_dat[i].int_S_max_gf_s[j][k]*dt[k]*geo_factor_cm_s*time_in_orbit/long_bins_calc)*el_dat[i].I_factor;
						el_dat[i].abund_S_min_dgf+=(el_dat[i].int_S_min_gf_s[j][k]*dt[k]*geo_factor_cm_s*time_in_orbit/long_bins_calc)*el_dat[i].I_factor;
					}

				}
			}
		}

	char t_string_abund[100]="CALET 5 Year Abundances";
	char t_string_S_min_abund[100]="CALET 5 Year Abundances for Solar Minimum Modulation";
	
	TH1F *h_S_min_abund = new TH1F("h_S_min_abund",t_string_S_min_abund,99,0.5,99.5);
	
		// h_S_min_abund
	h_S_min_abund->GetXaxis()->SetTitle("Z");
	h_S_min_abund->GetXaxis()->CenterTitle();
	h_S_min_abund->GetYaxis()->SetTitle("counts");
	h_S_min_abund->GetYaxis()->CenterTitle();
		//	h_S_min_abund->GetYaxis()->SetTitleOffset(1.4);
	h_S_min_abund->GetYaxis()->SetLabelFont(nicefont);
	h_S_min_abund->GetYaxis()->SetTitleFont(nicefont);
	h_S_min_abund->GetXaxis()->SetLabelFont(nicefont);
	h_S_min_abund->GetXaxis()->SetTitleFont(nicefont);

	char t_string_S_min_abund_gf[100]="CALET 5 Year Abundances for Solar Minimum Modulation";
	char l_string_S_min_abund_gf[100]="Solar Minimum";
	
	TH1F *h_S_min_abund_gf = new TH1F("h_S_min_abund_gf",t_string_S_min_abund_gf,99,0.5,99.5);
	
		// h_S_min_abund_gf
	h_S_min_abund_gf->GetXaxis()->SetTitle("Z");
	h_S_min_abund_gf->GetXaxis()->CenterTitle();
	h_S_min_abund_gf->GetYaxis()->SetTitle("counts");
	h_S_min_abund_gf->GetYaxis()->CenterTitle();
		//	h_S_min_abund_gf->GetYaxis()->SetTitleOffset(1.4);
	h_S_min_abund_gf->GetYaxis()->SetLabelFont(nicefont);
	h_S_min_abund_gf->GetYaxis()->SetTitleFont(nicefont);
	h_S_min_abund_gf->GetXaxis()->SetLabelFont(nicefont);
	h_S_min_abund_gf->GetXaxis()->SetTitleFont(nicefont);

	
	char t_string_S_min_abund_dgf[100]="CALET 5 Year Abundances for Solar Minimum Modulation";
	char l_string_S_min_abund_dgf[100]="Solar Minimum";
	
	TH1F *h_S_min_abund_dgf = new TH1F("h_S_min_abund_dgf","CALET 5 Year Abundances",99,0.5,99.5);
	
		// h_S_min_abund_dgf
	h_S_min_abund_dgf->GetXaxis()->SetTitle("Z");
	h_S_min_abund_dgf->GetXaxis()->CenterTitle();
	h_S_min_abund_dgf->GetYaxis()->SetTitle("counts");
	h_S_min_abund_dgf->GetYaxis()->CenterTitle();
		//	h_S_min_abund_dgf->GetYaxis()->SetTitleOffset(1.4);
	h_S_min_abund_dgf->GetYaxis()->SetLabelFont(nicefont);
	h_S_min_abund_dgf->GetYaxis()->SetTitleFont(nicefont);
	h_S_min_abund_dgf->GetXaxis()->SetLabelFont(nicefont);
	h_S_min_abund_dgf->GetXaxis()->SetTitleFont(nicefont);
	
	char t_string_S_max_abund[100]="CALET 5 Year Abundances for Solar Maximum Modulation";
	
	TH1F *h_S_max_abund = new TH1F("h_S_max_abund",t_string_S_max_abund,99,0.5,99.5);
	
		// h_S_max_abund
	h_S_max_abund->GetXaxis()->SetTitle("Z");
	h_S_max_abund->GetXaxis()->CenterTitle();
	h_S_max_abund->GetYaxis()->SetTitle("counts");
	h_S_max_abund->GetYaxis()->CenterTitle();
		//	h_S_max_abund->GetYaxis()->SetTitleOffset(1.4);
	h_S_max_abund->GetYaxis()->SetLabelFont(nicefont);
	h_S_max_abund->GetYaxis()->SetTitleFont(nicefont);
	h_S_max_abund->GetXaxis()->SetLabelFont(nicefont);
	h_S_max_abund->GetXaxis()->SetTitleFont(nicefont);
	
	char t_string_S_max_abund_gf[100]="CALET 5 Year Abundances for Solar Maximum Modulation";
	char l_string_S_max_abund_gf[100]="Solar Maximum";

	TH1F *h_S_max_abund_gf = new TH1F("h_S_max_abund_gf",t_string_S_max_abund_gf,99,0.5,99.5);
	
		// h_S_max_abund_gf
	h_S_max_abund_gf->GetXaxis()->SetTitle("Z");
	h_S_max_abund_gf->GetXaxis()->CenterTitle();
	h_S_max_abund_gf->GetYaxis()->SetTitle("counts");
	h_S_max_abund_gf->GetYaxis()->CenterTitle();
		//	h_S_max_abund_gf->GetYaxis()->SetTitleOffset(1.4);
	h_S_max_abund_gf->GetYaxis()->SetLabelFont(nicefont);
	h_S_max_abund_gf->GetYaxis()->SetTitleFont(nicefont);
	h_S_max_abund_gf->GetXaxis()->SetLabelFont(nicefont);
	h_S_max_abund_gf->GetXaxis()->SetTitleFont(nicefont);

	
	char t_string_S_max_abund_dgf[100]="CALET 5 Year Abundances for Solar Maximum Modulation";
	char l_string_S_max_abund_dgf[100]="Solar Maximum";
	
	TH1F *h_S_max_abund_dgf = new TH1F("h_S_max_abund_dgf",t_string_S_max_abund_dgf,99,0.5,99.5);
	
		// h_S_max_abund_dgf
	h_S_max_abund_dgf->GetXaxis()->SetTitle("Z");
	h_S_max_abund_dgf->GetXaxis()->CenterTitle();
	h_S_max_abund_dgf->GetYaxis()->SetTitle("counts");
	h_S_max_abund_dgf->GetYaxis()->CenterTitle();
		//	h_S_max_abund_dgf->GetYaxis()->SetTitleOffset(1.4);
	h_S_max_abund_dgf->GetYaxis()->SetLabelFont(nicefont);
	h_S_max_abund_dgf->GetYaxis()->SetTitleFont(nicefont);
	h_S_max_abund_dgf->GetXaxis()->SetLabelFont(nicefont);
	h_S_max_abund_dgf->GetXaxis()->SetTitleFont(nicefont);
	
	char t_string_S_ave_abund[100]="CALET 5 Year Abundances for Average Solar Modulation";
	
	TH1F *h_S_ave_abund = new TH1F("h_S_ave_abund",t_string_S_ave_abund,99,0.5,99.5);
	
		// h_S_ave_abund
	h_S_ave_abund->GetXaxis()->SetTitle("Z");
	h_S_ave_abund->GetXaxis()->CenterTitle();
	h_S_ave_abund->GetYaxis()->SetTitle("counts");
	h_S_ave_abund->GetYaxis()->CenterTitle();
		//	h_S_ave_abund->GetYaxis()->SetTitleOffset(1.4);
	h_S_ave_abund->GetYaxis()->SetLabelFont(nicefont);
	h_S_ave_abund->GetYaxis()->SetTitleFont(nicefont);
	h_S_ave_abund->GetXaxis()->SetLabelFont(nicefont);
	h_S_ave_abund->GetXaxis()->SetTitleFont(nicefont);
	
	char t_string_S_ave_abund_gf[100]="CALET 5 Year Abundances for Average Solar Modulation";
	char l_string_S_ave_abund_gf[100]="With Interactions";
	
	TH1F *h_S_ave_abund_gf = new TH1F("h_S_ave_abund_gf",t_string_S_ave_abund_gf,99,0.5,99.5);
	
		// h_S_ave_abund_gf
	h_S_ave_abund_gf->GetXaxis()->SetTitle("Z");
	h_S_ave_abund_gf->GetXaxis()->CenterTitle();
	h_S_ave_abund_gf->GetYaxis()->SetTitle("counts");
	h_S_ave_abund_gf->GetYaxis()->CenterTitle();
		//	h_S_ave_abund_gf->GetYaxis()->SetTitleOffset(1.4);
	h_S_ave_abund_gf->GetYaxis()->SetLabelFont(nicefont);
	h_S_ave_abund_gf->GetYaxis()->SetTitleFont(nicefont);
	h_S_ave_abund_gf->GetXaxis()->SetLabelFont(nicefont);
	h_S_ave_abund_gf->GetXaxis()->SetTitleFont(nicefont);
	
	char t_string_S_ave_abund_dgf[100]="CALET 5 Year Abundances for Differential Geometry Factors";
	char l_string_S_ave_abund_dgf[100]="With Interactions";
	
	TH1F *h_S_ave_abund_dgf = new TH1F("h_S_ave_abund_dgf",t_string_S_ave_abund_dgf,99,0.5,99.5);
	
		// h_S_ave_abund_dgf
	h_S_ave_abund_dgf->GetXaxis()->SetTitle("Z");
	h_S_ave_abund_dgf->GetXaxis()->CenterTitle();
	h_S_ave_abund_dgf->GetYaxis()->SetTitle("counts");
	h_S_ave_abund_dgf->GetYaxis()->CenterTitle();
		//	h_S_ave_abund_dgf->GetYaxis()->SetTitleOffset(1.4);
	h_S_ave_abund_dgf->GetYaxis()->SetLabelFont(nicefont);
	h_S_ave_abund_dgf->GetYaxis()->SetTitleFont(nicefont);
	h_S_ave_abund_dgf->GetXaxis()->SetLabelFont(nicefont);
	h_S_ave_abund_dgf->GetXaxis()->SetTitleFont(nicefont);	

	
	char t_string_S_ave_abund_dgf_noic[100]="CALET 5 Year Abundances for Differential Geometry Factors";
	char l_string_S_ave_abund_dgf_noic[100]="Without Interactions";
	
	TH1F *h_S_ave_abund_dgf_noic = new TH1F("h_S_ave_abund_dgf_noic",t_string_S_ave_abund_dgf_noic,99,0.5,99.5);
	
		// h_S_ave_abund_dgf_noic
	h_S_ave_abund_dgf_noic->GetXaxis()->SetTitle("Z");
	h_S_ave_abund_dgf_noic->GetXaxis()->CenterTitle();
	h_S_ave_abund_dgf_noic->GetYaxis()->SetTitle("counts");
	h_S_ave_abund_dgf_noic->GetYaxis()->CenterTitle();
		//	h_S_ave_abund_dgf_noic->GetYaxis()->SetTitleOffset(1.4);
	h_S_ave_abund_dgf_noic->GetYaxis()->SetLabelFont(nicefont);
	h_S_ave_abund_dgf_noic->GetYaxis()->SetTitleFont(nicefont);
	h_S_ave_abund_dgf_noic->GetXaxis()->SetLabelFont(nicefont);
	h_S_ave_abund_dgf_noic->GetXaxis()->SetTitleFont(nicefont);	
	
	char t_string_S_ave_abund_gf_no_ic[100]="CALET 5 Year Abundances for Average Solar Modulation";
	
	TH1F *h_S_ave_abund_gf_no_ic = new TH1F("h_S_ave_abund_gf_no_ic",t_string_S_ave_abund_gf_no_ic,99,0.5,99.5);
	char l_string_S_ave_abund_gf_no_ic[100]="Before Interactions";

		// h_S_ave_abund_gf_no_ic
	h_S_ave_abund_gf_no_ic->GetXaxis()->SetTitle("Z");
	h_S_ave_abund_gf_no_ic->GetXaxis()->CenterTitle();
	h_S_ave_abund_gf_no_ic->GetYaxis()->SetTitle("counts");
	h_S_ave_abund_gf_no_ic->GetYaxis()->CenterTitle();
		//	h_S_ave_abund_gf_no_ic->GetYaxis()->SetTitleOffset(1.4);
	h_S_ave_abund_gf_no_ic->GetYaxis()->SetLabelFont(nicefont);
	h_S_ave_abund_gf_no_ic->GetYaxis()->SetTitleFont(nicefont);
	h_S_ave_abund_gf_no_ic->GetXaxis()->SetLabelFont(nicefont);
	h_S_ave_abund_gf_no_ic->GetXaxis()->SetTitleFont(nicefont);	
	 
	double CALET_abund[100];
	for (i=0; i<100; i++) {
		CALET_abund[i]=0.0;
	}
		
	for(i=0;i<el_index;i++)
		{
		if(i>3 && i<92) 
			{
			if(i<26 || i == 27)
				{
				el_dat[i].abund_S_max_tot=el_dat[i].abund_S_max+S_max[i+1].I_E_max*geo_factor_cm*time_in_orbit*el_dat[i].I_factor;
				el_dat[i].abund_S_min_tot=el_dat[i].abund_S_min+S_min[i+1].I_E_max*geo_factor_cm*time_in_orbit*el_dat[i].I_factor;
				}
			else
				{
				el_dat[i].abund_S_max_tot=el_dat[i].abund_S_max+S_max[26].I_E_max*geo_factor_cm*time_in_orbit*el_dat[i].TOA_abund*el_dat[i].I_factor;
				el_dat[i].abund_S_min_tot=el_dat[i].abund_S_min+S_min[26].I_E_max*geo_factor_cm*time_in_orbit*el_dat[i].TOA_abund*el_dat[i].I_factor;
				}
				h_S_min_abund->Fill(el_dat[i].Z,el_dat[i].abund_S_min*el_dat[i].I_factor);
				h_S_max_abund->Fill(el_dat[i].Z,el_dat[i].abund_S_max*el_dat[i].I_factor);
				h_S_ave_abund->Fill(el_dat[i].Z,el_dat[i].abund_S_min*0.5*el_dat[i].I_factor);
				h_S_ave_abund->Fill(el_dat[i].Z,el_dat[i].abund_S_max*0.5*el_dat[i].I_factor);
				h_S_min_abund_gf->Fill(el_dat[i].Z,el_dat[i].abund_S_min_gf*el_dat[i].I_factor);
				h_S_max_abund_gf->Fill(el_dat[i].Z,el_dat[i].abund_S_max_gf*el_dat[i].I_factor);
				h_S_ave_abund_gf->Fill(el_dat[i].Z,el_dat[i].abund_S_min_gf*0.5*el_dat[i].I_factor);
				h_S_ave_abund_gf->Fill(el_dat[i].Z,el_dat[i].abund_S_max_gf*0.5*el_dat[i].I_factor);
				h_S_ave_abund_dgf->Fill(el_dat[i].Z,el_dat[i].abund_S_min_dgf*0.5);
				h_S_ave_abund_dgf->Fill(el_dat[i].Z,el_dat[i].abund_S_max_dgf*0.5);
				h_S_min_abund_dgf->Fill(el_dat[i].Z,el_dat[i].abund_S_min_dgf);
				h_S_max_abund_dgf->Fill(el_dat[i].Z,el_dat[i].abund_S_max_dgf);
				h_S_ave_abund_dgf_noic->Fill(el_dat[i].Z,el_dat[i].abund_S_min_dgf_noic*0.5);
				h_S_ave_abund_dgf_noic->Fill(el_dat[i].Z,el_dat[i].abund_S_max_dgf_noic*0.5);
				CALET_abund[i]+=el_dat[i].abund_S_min_dgf*0.5;
				CALET_abund[i]+=el_dat[i].abund_S_max_dgf*0.5;				
				h_S_ave_abund_gf_no_ic->Fill(el_dat[i].Z,el_dat[i].abund_S_min_gf*0.5);
				h_S_ave_abund_gf_no_ic->Fill(el_dat[i].Z,el_dat[i].abund_S_max_gf*0.5);
            
		fprintf(stderr,"%2d el_dat[%2d].abund_S_max: %.4e %.4e %.4e %.4f %.4e %.4f el_dat[%2d].abund_S_min: %.4e %.4e %.4e %.4f %.4e %.4f\n",
									el_dat[i].Z,i,el_dat[i].abund_S_max,el_dat[i].abund_S_max_tot,el_dat[i].abund_S_max_tot-el_dat[i].abund_S_max,
									(el_dat[i].abund_S_max_tot-el_dat[i].abund_S_max)/el_dat[i].abund_S_max,el_dat[i].abund_S_max_no_mag,
									el_dat[i].abund_S_max_tot/el_dat[i].abund_S_max_no_mag,
									i,el_dat[i].abund_S_min,el_dat[i].abund_S_min_tot,el_dat[i].abund_S_min_tot-el_dat[i].abund_S_min,
									(el_dat[i].abund_S_min_tot-el_dat[i].abund_S_min)/el_dat[i].abund_S_min,el_dat[i].abund_S_min_no_mag,
									el_dat[i].abund_S_min_tot/el_dat[i].abund_S_min_no_mag);
			}
		}

	
		// Interlude for CHD limiting factors plots
	
	double CHD_Z[100];
	double CHD_S[100];			// PEs (photo electrons)
	double CHD_S_30[100];		// PEs (photo electrons) for 30 degrees incidence
	double CHD_S_45[100];		// PEs (photo electrons) for 45 degrees incidence
	double CHD_PMT[100];		// fC (10^(-15)C)
	double CHD_PMT_30[100];		// fC (10^(-15)C) for 30 degrees incidence
	double CHD_PMT_45[100];		// fC (10^(-15)C) for 45 degrees incidence
	double CHD_CSA[100];		// mV
	double CHD_CSA_30[100];		// mV for 30 degrees incidence
	double CHD_CSA_45[100];		// mV for 45 degrees incidence
	double CHD_ADC[100];		// counts
	double CHD_ADC_30[100];		// counts for 30 degrees incidence
	double CHD_ADC_45[100];		// counts for 45 degrees incidence

	double p_SI=1e-12;
	double f_SI=1e-15;
	double m_SI=1e-3;
	
	double charge_e=1.60217646e-19;	// Charge of an electron
	
	double PMT_max_ideal=375000;	// (fC) from Calet-CHD range.pptx
	double PMT_max_reduced=320000;	// (fC) from Calet-CHD range.pptx for linearity < 2%
	double CSA_max=3000;			// (mV) from Calet-CHD range.pptx
	double ADC_max=32768;			// (counts) from Calet-CHD range.pptx
	
	double	PMT_gain=5000;	// From Calet-CHD range.pptx
	double CSA_cap=100;		// (pF) From Calet-CHD range.pptx (10^(-12)C)
	
	CHD_Z[0]=1;			// for the MIP (minimum ionizing particle)
	CHD_Z[1]=26;
	CHD_S[0]=300;		// From Calet-CHD range.pptx
	CHD_S[1]=81000;		// From Calet-CHD range.pptx
	
	CHD_S_30[0]=CHD_S[0]/cos(30*TMath::Pi()/180);
	CHD_S_30[1]=CHD_S[1]/cos(30*TMath::Pi()/180);
	
	CHD_S_45[0]=CHD_S[0]/cos(45*TMath::Pi()/180);
	CHD_S_45[1]=CHD_S[1]/cos(45*TMath::Pi()/180);
	
	
	CHD_PMT[0]=CHD_S[0]*charge_e*PMT_gain/f_SI;
	CHD_PMT[1]=CHD_S[1]*charge_e*PMT_gain/f_SI;
	
	CHD_PMT_30[0]=(CHD_S[0]*charge_e*PMT_gain/f_SI)/cos(30*TMath::Pi()/180);
	CHD_PMT_30[1]=(CHD_S[1]*charge_e*PMT_gain/f_SI)/cos(30*TMath::Pi()/180);

	CHD_PMT_45[0]=(CHD_S[0]*charge_e*PMT_gain/f_SI)/cos(45*TMath::Pi()/180);
	CHD_PMT_45[1]=(CHD_S[1]*charge_e*PMT_gain/f_SI)/cos(45*TMath::Pi()/180);
	
	fprintf(stderr,"CHD_PMT[%d]: %e CHD_PMT[%d]: %e\n",0,CHD_PMT[0],1,CHD_PMT[1]);
	fprintf(stderr,"CHD_PMT_30[%d]: %e CHD_PMT_30[%d]: %e\n",0,CHD_PMT_30[0],1,CHD_PMT_45[1]);
	fprintf(stderr,"CHD_PMT_45[%d]: %e CHD_PMT_45[%d]: %e\n",0,CHD_PMT_30[0],1,CHD_PMT_45[1]);
	
	//  Simply derive the charge collected by the charge sensitive amplifier (CSA) using the capacitor value. C=Q/V, V=Q/C
	
	CHD_CSA[0]=(CHD_PMT[0]/(CSA_cap*p_SI/f_SI))/m_SI;			
	CHD_CSA[1]=(CHD_PMT[1]/(CSA_cap*p_SI/f_SI))/m_SI;
	
	CHD_CSA_30[0]=(CHD_PMT_30[0]/(CSA_cap*p_SI/f_SI))/m_SI;			
	CHD_CSA_30[1]=(CHD_PMT_30[1]/(CSA_cap*p_SI/f_SI))/m_SI;
	
	
	CHD_CSA_45[0]=(CHD_PMT_45[0]/(CSA_cap*p_SI/f_SI))/m_SI;			
	CHD_CSA_45[1]=(CHD_PMT_45[1]/(CSA_cap*p_SI/f_SI))/m_SI;
	
	
	fprintf(stderr,"CHD_CSA[%d]: %e CHD_CSA[%d]: %e\n",0,CHD_CSA[0],1,CHD_CSA[1]);
	fprintf(stderr,"CHD_CSA_30[%d]: %e CHD_CSA_30[%d]: %e\n",0,CHD_CSA_30[0],1,CHD_CSA_30[1]);
	fprintf(stderr,"CHD_CSA_45[%d]: %e CHD_CSA_45[%d]: %e\n",0,CHD_CSA_45[0],1,CHD_CSA_45[1]);

	CHD_ADC[0]=32;		// From Calet-CHD range.pptx
	CHD_ADC[1]=8600;	// From Calet-CHD range.pptx

	CHD_ADC_30[0]=CHD_ADC[0]/cos(30*TMath::Pi()/180);
	CHD_ADC_45[0]=CHD_ADC[0]/cos(45*TMath::Pi()/180);
	
	CHD_ADC_30[1]=CHD_ADC[1]/cos(30*TMath::Pi()/180);
	CHD_ADC_45[1]=CHD_ADC[1]/cos(45*TMath::Pi()/180);
	
	double CHD_ADC_factor=0.0; // Derived quantity

	CHD_ADC_factor=CHD_ADC[0]/CHD_CSA[0];
	
	fprintf(stderr,"Z: %2.0f CHD_ADC_factor: %f\n",CHD_Z[0],CHD_ADC_factor);
	
	CHD_ADC_factor=CHD_ADC[1]/CHD_CSA[1];
	
	fprintf(stderr,"Z: %2.0f CHD_ADC_factor: %f\n",CHD_Z[1],CHD_ADC_factor);
	
	int CHD_count=0;
	
	for (CHD_count=2; CHD_count<68; CHD_count++) {
		CHD_Z[CHD_count]=CHD_count+25;
		CHD_S[CHD_count]=pow((CHD_Z[CHD_count]/CHD_Z[1]),1.71)*CHD_S[1];
		CHD_PMT[CHD_count]=CHD_S[CHD_count]*charge_e*PMT_gain/f_SI;
		CHD_CSA[CHD_count]=(CHD_PMT[CHD_count]/(CSA_cap*p_SI/f_SI))/m_SI;
		CHD_ADC[CHD_count]=CHD_CSA[CHD_count]*CHD_ADC_factor;

		CHD_S_30[CHD_count]=pow((CHD_Z[CHD_count]/CHD_Z[1]),1.71)*CHD_S_30[1];
		CHD_PMT_30[CHD_count]=CHD_S_30[CHD_count]*charge_e*PMT_gain/f_SI;
		CHD_CSA_30[CHD_count]=(CHD_PMT_30[CHD_count]/(CSA_cap*p_SI/f_SI))/m_SI;
		CHD_ADC_30[CHD_count]=CHD_CSA_30[CHD_count]*CHD_ADC_factor;

		CHD_S_45[CHD_count]=pow((CHD_Z[CHD_count]/CHD_Z[1]),1.71)*CHD_S_45[1];
		CHD_PMT_45[CHD_count]=CHD_S_45[CHD_count]*charge_e*PMT_gain/f_SI;
		CHD_CSA_45[CHD_count]=(CHD_PMT_45[CHD_count]/(CSA_cap*p_SI/f_SI))/m_SI;
		CHD_ADC_45[CHD_count]=CHD_CSA_45[CHD_count]*CHD_ADC_factor;
	}
	
	TGraph *gr_CHD_S = new TGraph(CHD_count,CHD_Z,CHD_S);
	TGraph *gr_CHD_PMT = new TGraph(CHD_count,CHD_Z,CHD_PMT);
	TGraph *gr_CHD_CSA = new TGraph(CHD_count,CHD_Z,CHD_CSA);
	TGraph *gr_CHD_ADC = new TGraph(CHD_count,CHD_Z,CHD_ADC);
	
	
	TGraph *gr_CHD_S_30 = new TGraph(CHD_count,CHD_Z,CHD_S_30);
	TGraph *gr_CHD_PMT_30 = new TGraph(CHD_count,CHD_Z,CHD_PMT_30);
	TGraph *gr_CHD_CSA_30 = new TGraph(CHD_count,CHD_Z,CHD_CSA_30);
	TGraph *gr_CHD_ADC_30 = new TGraph(CHD_count,CHD_Z,CHD_ADC_30);
	
	
	TGraph *gr_CHD_S_45 = new TGraph(CHD_count,CHD_Z,CHD_S_45);
	TGraph *gr_CHD_PMT_45 = new TGraph(CHD_count,CHD_Z,CHD_PMT_45);
	TGraph *gr_CHD_CSA_45 = new TGraph(CHD_count,CHD_Z,CHD_CSA_45);
	TGraph *gr_CHD_ADC_45 = new TGraph(CHD_count,CHD_Z,CHD_ADC_45);
	
	TGraph *gr_CHD_S_45_eval = new TGraph(CHD_count,CHD_S_45,CHD_Z);
	TGraph *gr_CHD_PMT_45_eval = new TGraph(CHD_count,CHD_PMT_45,CHD_Z);
	TGraph *gr_CHD_CSA_45_eval = new TGraph(CHD_count,CHD_CSA_45,CHD_Z);
	TGraph *gr_CHD_ADC_45_eval = new TGraph(CHD_count,CHD_ADC_45,CHD_Z);


	
	double CHD_PMT_45_Z_max_ideal=gr_CHD_S_45_eval->Eval(PMT_max_ideal);
	double CHD_PMT_45_Z_max_reduced=gr_CHD_S_45_eval->Eval(PMT_max_reduced);
	double CHD_CSA_45_Z_max=gr_CHD_CSA_45_eval->Eval(CSA_max);
	double CHD_ADC_45_Z_max=gr_CHD_ADC_45_eval->Eval(ADC_max);

	fprintf(stderr,"CHD_PMT_45_Z_max_ideal: %.2f\n",CHD_PMT_45_Z_max_ideal);
	fprintf(stderr,"CHD_PMT_45_Z_max_reduced: %.2f\n",CHD_PMT_45_Z_max_reduced);
	fprintf(stderr,"CHD_CSA_45_Z_max: %.2f\n",CHD_CSA_45_Z_max);
	fprintf(stderr,"CHD_ADC_45_Z_max: %.2f\n",CHD_ADC_45_Z_max);
	
	
	
	
	TGraph *gr_CHD_S_m = new TGraph(2,CHD_Z,CHD_S);
	TGraph *gr_CHD_PMT_m = new TGraph(2,CHD_Z,CHD_PMT);
	TGraph *gr_CHD_CSA_m = new TGraph(2,CHD_Z,CHD_CSA);
	TGraph *gr_CHD_ADC_m = new TGraph(2,CHD_Z,CHD_ADC);
	
	
	TGraph *gr_CHD_S_30_m = new TGraph(2,CHD_Z,CHD_S_30);
	TGraph *gr_CHD_PMT_30_m = new TGraph(2,CHD_Z,CHD_PMT_30);
	TGraph *gr_CHD_CSA_30_m = new TGraph(2,CHD_Z,CHD_CSA_30);
	TGraph *gr_CHD_ADC_30_m = new TGraph(2,CHD_Z,CHD_ADC_30);
	
	
	TGraph *gr_CHD_S_45_m = new TGraph(2,CHD_Z,CHD_S_45);
	TGraph *gr_CHD_PMT_45_m = new TGraph(2,CHD_Z,CHD_PMT_45);
	TGraph *gr_CHD_CSA_45_m = new TGraph(2,CHD_Z,CHD_CSA_45);
	TGraph *gr_CHD_ADC_45_m = new TGraph(2,CHD_Z,CHD_ADC_45);
	
	double CHD_Z_range[2];
	CHD_Z_range[0]=0;
	CHD_Z_range[1]=101;
	
	double CHD_PMT_max_reduced[2];
	CHD_PMT_max_reduced[0]=PMT_max_reduced;
	CHD_PMT_max_reduced[1]=PMT_max_reduced;

	double CHD_PMT_max_ideal[2];
	CHD_PMT_max_ideal[0]=PMT_max_ideal;
	CHD_PMT_max_ideal[1]=PMT_max_ideal;
	
	double CHD_CSA_max[2];
	CHD_CSA_max[0]=CSA_max;
	CHD_CSA_max[1]=CSA_max;

	double CHD_ADC_max[2];
	CHD_ADC_max[0]=ADC_max;
	CHD_ADC_max[1]=ADC_max;
	
	TGraph *gr_CHD_PMT_max_reduced = new TGraph(2,CHD_Z_range,CHD_PMT_max_reduced);
	TGraph *gr_CHD_PMT_max_ideal = new TGraph(2,CHD_Z_range,CHD_PMT_max_ideal);
	TGraph *gr_CHD_CSA_max = new TGraph(2,CHD_Z_range,CHD_CSA_max);
	TGraph *gr_CHD_ADC_max = new TGraph(2,CHD_Z_range,CHD_ADC_max);

	
	
	double label_x_CHD[2];
	double label_y_CHD[2];
	label_x_CHD[0] = 0;
	label_x_CHD[1] = 100;
	label_y_CHD[0] = 0.1;
	label_y_CHD[1] = 1e6;
	
	TGraph *gr_label_S_CHD = new TGraph(2,label_x_CHD,label_y_CHD);
	
	gr_label_S_CHD->SetMarkerStyle(1);
	gr_label_S_CHD->SetMarkerSize(1);
		//	gr_label_S_CHD->SetTitle("Limiting Factors in CHD Measurable Range");
	gr_label_S_CHD->SetTitle("");
	gr_label_S_CHD->GetXaxis()->SetTitle("Z");
	gr_label_S_CHD->GetXaxis()->CenterTitle();
	gr_label_S_CHD->GetYaxis()->SetTitle("Signal");
	gr_label_S_CHD->GetYaxis()->CenterTitle();
	gr_label_S_CHD->GetXaxis()->SetTitleOffset(1.2);
	gr_label_S_CHD->GetYaxis()->SetTitleOffset(1.2);
	gr_label_S_CHD->GetYaxis()->SetLabelFont(nicefont);
	gr_label_S_CHD->GetYaxis()->SetTitleFont(nicefont);
	gr_label_S_CHD->GetXaxis()->SetLabelFont(nicefont);
	gr_label_S_CHD->GetXaxis()->SetTitleFont(nicefont);
	gr_label_S_CHD->GetXaxis()->SetRangeUser(26,80);
	gr_label_S_CHD->GetYaxis()->SetRangeUser(1e2,1e6);
	
	
	
	TCanvas *C_CHD = new TCanvas("C_CHD","CALET CHD Measureable Range",0,0,800,800);
	gr_label_S_CHD->Draw("AP");
	
	gr_CHD_PMT_max_reduced->SetLineColor(2);
	gr_CHD_PMT_max_reduced->SetLineStyle(2);
	gr_CHD_PMT_max_reduced->SetLineWidth(2);
	gr_CHD_PMT_max_reduced->Draw("L");
	gr_CHD_PMT_max_ideal->SetLineColor(2);
	gr_CHD_PMT_max_ideal->SetLineStyle(2);
	gr_CHD_PMT_max_ideal->SetLineWidth(2);
	gr_CHD_PMT_max_ideal->Draw("L");
	gr_CHD_CSA_max->SetLineColor(3);
	gr_CHD_CSA_max->SetLineStyle(2);
	gr_CHD_CSA_max->SetLineWidth(2);
	gr_CHD_CSA_max->Draw("L");
	gr_CHD_ADC_max->SetLineColor(4);
	gr_CHD_ADC_max->SetLineStyle(2);
	gr_CHD_ADC_max->SetLineWidth(2);
	gr_CHD_ADC_max->Draw("L");
	
	gr_CHD_S->SetLineColor(1);
	gr_CHD_S->Draw("L");
	gr_CHD_PMT->SetLineColor(2);
	gr_CHD_PMT->Draw("L");
	gr_CHD_CSA->SetLineColor(3);
	gr_CHD_CSA->Draw("L");
	gr_CHD_ADC->SetLineColor(4);
	gr_CHD_ADC->Draw("L");
	
		//	gr_CHD_S_m->SetMarkerStyle(23);
		//	gr_CHD_S_m->Draw("P");
/*
	gr_CHD_PMT_m->SetLineColor(2);
	gr_CHD_PMT_m->SetMarkerColor(2);
	gr_CHD_PMT_m->SetMarkerStyle(23);
	gr_CHD_PMT_m->Draw("P");
	gr_CHD_CSA_m->SetLineColor(3);
	gr_CHD_CSA_m->SetMarkerColor(3);
	gr_CHD_CSA_m->SetMarkerStyle(23);
	gr_CHD_CSA_m->Draw("P");
	gr_CHD_ADC_m->SetLineColor(4);
	gr_CHD_ADC_m->SetMarkerColor(4);
	gr_CHD_ADC_m->SetMarkerStyle(23);
	gr_CHD_ADC_m->Draw("P");
	*/
	
	gr_CHD_S_30->SetLineColor(1);
	gr_CHD_S_30->Draw("L");
	gr_CHD_PMT_30->SetLineColor(2);
	gr_CHD_PMT_30->Draw("L");
	gr_CHD_CSA_30->SetLineColor(3);
	gr_CHD_CSA_30->Draw("L");
	gr_CHD_ADC_30->SetLineColor(4);
	gr_CHD_ADC_30->Draw("L");

	
	gr_CHD_S_30_m->SetMarkerStyle(21);
		//	gr_CHD_S_30_m->Draw("P");
	gr_CHD_PMT_30_m->SetMarkerStyle(21);
	gr_CHD_PMT_30_m->SetMarkerColor(2);
	gr_CHD_PMT_30_m->SetLineColor(2);
		//	gr_CHD_PMT_30_m->Draw("P");
	gr_CHD_CSA_30_m->SetMarkerStyle(21);
	gr_CHD_CSA_30_m->SetMarkerColor(3);
	gr_CHD_CSA_30_m->SetLineColor(3);
		//	gr_CHD_CSA_30_m->Draw("P");
	gr_CHD_ADC_30_m->SetMarkerStyle(21);
	gr_CHD_ADC_30_m->SetMarkerColor(4);
	gr_CHD_ADC_30_m->SetLineColor(4);
		//	gr_CHD_ADC_30_m->Draw("P");
	
	
	gr_CHD_S_45->SetLineColor(1);
	gr_CHD_S_45->Draw("L");
	gr_CHD_PMT_45->SetLineColor(2);
	gr_CHD_PMT_45->Draw("L");
	gr_CHD_CSA_45->SetLineColor(3);
	gr_CHD_CSA_45->Draw("L");
	gr_CHD_ADC_45->SetLineColor(4);
	gr_CHD_ADC_45->Draw("L");
	
		//	gr_CHD_S_45_m->SetMarkerStyle(22);
		//	gr_CHD_S_45_m->Draw("P");
/*
	gr_CHD_PMT_45_m->SetMarkerStyle(22);
	gr_CHD_PMT_45_m->SetMarkerColor(2);
	gr_CHD_PMT_45_m->SetLineColor(2);
	gr_CHD_PMT_45_m->Draw("P");
	gr_CHD_CSA_45_m->SetLineColor(3);
	gr_CHD_CSA_45_m->SetMarkerColor(3);
	gr_CHD_CSA_45_m->SetMarkerStyle(22);
	gr_CHD_CSA_45_m->Draw("P");
	gr_CHD_ADC_45_m->SetLineColor(4);
	gr_CHD_ADC_45_m->SetMarkerColor(4);
	gr_CHD_ADC_45_m->SetMarkerStyle(22);
	gr_CHD_ADC_45_m->Draw("P");
	*/

	
	
	TPaveText *pPMT_max = new TPaveText(0.118,0.828,0.354,0.858,"NDC");
	pPMT_max->SetBorderSize(0);
	pPMT_max->SetFillColor(kWhite);
	TText *tPMT_max = pPMT_max->AddText("PMT Maximum Signal");
	tPMT_max->SetTextFont(nicefont);
	tPMT_max->SetTextFont(nicefont);
	tPMT_max->SetTextColor(2);
	
	pPMT_max->Draw();
	

	TPaveText *pCSA_max = new TPaveText(0.118,0.408,0.354,0.438,"NDC");
	pCSA_max->SetBorderSize(0);
	pCSA_max->SetFillColor(kWhite);
	TText *tCSA_max = pCSA_max->AddText("CSA Maximum Signal");
	tCSA_max->SetTextFont(nicefont);
	tCSA_max->SetTextFont(nicefont);
	tCSA_max->SetTextColor(3);
	
	pCSA_max->Draw();

	
	TPaveText *pADC_max = new TPaveText(0.118,0.613,0.354,0.642,"NDC");
	pADC_max->SetBorderSize(0);
	pADC_max->SetFillColor(kWhite);
	TText *tADC_max = pADC_max->AddText("ADC Maximum Signal");
	tADC_max->SetTextFont(nicefont);
	tADC_max->SetTextFont(nicefont);
	tADC_max->SetTextColor(4);
	
	pADC_max->Draw();
	
		// ADC angle ranges
	
	TPaveText *p_ADC_0 = new TPaveText(0.903,0.627,0.936,0.655,"NDC");
	p_ADC_0->SetBorderSize(0);
	p_ADC_0->SetFillColor(kWhite);
	TText *t_ADC_0 = p_ADC_0->AddText(" 0^{o}");
	t_ADC_0->SetTextFont(nicefont);
	t_ADC_0->SetTextFont(nicefont);
	t_ADC_0->SetTextColor(4);
	
	p_ADC_0->Draw();
	
	
	TPaveText *p_ADC_30 = new TPaveText(0.903,0.655,0.936,0.683,"NDC");
	p_ADC_30->SetBorderSize(0);
	p_ADC_30->SetFillColor(kWhite);
	TText *t_ADC_30 = p_ADC_30->AddText("30^{o}");
	t_ADC_30->SetTextFont(nicefont);
	t_ADC_30->SetTextFont(nicefont);
	t_ADC_30->SetTextColor(4);
	
	p_ADC_30->Draw();
	
	
	
	TPaveText *p_ADC_45 = new TPaveText(0.903,0.683,0.936,0.711,"NDC");
	p_ADC_45->SetBorderSize(0);
	p_ADC_45->SetFillColor(kWhite);
	TText *t_ADC_45 = p_ADC_45->AddText("45^{o}");
	t_ADC_45->SetTextFont(nicefont);
	t_ADC_45->SetTextFont(nicefont);
	t_ADC_45->SetTextColor(4);
	
	p_ADC_45->Draw();
	
	
		// PMT angle ranges
	
	
	TPaveText *p_PMT_0 = new TPaveText(0.903,0.807,0.936,0.835,"NDC");
	p_PMT_0->SetBorderSize(0);
	p_PMT_0->SetFillColor(kWhite);
	TText *t_PMT_0 = p_PMT_0->AddText(" 0^{o}");
	t_PMT_0->SetTextFont(nicefont);
	t_PMT_0->SetTextFont(nicefont);
	t_PMT_0->SetTextColor(2);
	
	p_PMT_0->Draw();
	
	
	TPaveText *p_PMT_30 = new TPaveText(0.903,0.835,0.936,0.863,"NDC");
	p_PMT_30->SetBorderSize(0);
	p_PMT_30->SetFillColor(kWhite);
	TText *t_PMT_30 = p_PMT_30->AddText("30^{o}");
	t_PMT_30->SetTextFont(nicefont);
	t_PMT_30->SetTextFont(nicefont);
	t_PMT_30->SetTextColor(2);
	
	p_PMT_30->Draw();
	
	
	
	TPaveText *p_PMT_45 = new TPaveText(0.903,0.863,0.936,0.891,"NDC");
	p_PMT_45->SetBorderSize(0);
	p_PMT_45->SetFillColor(kWhite);
	TText *t_PMT_45 = p_PMT_45->AddText("45^{o}");
	t_PMT_45->SetTextFont(nicefont);
	t_PMT_45->SetTextFont(nicefont);
	t_PMT_45->SetTextColor(2);
	
	p_PMT_45->Draw();
	
		// CSA angle ranges
	
	TPaveText *p_CSA_0 = new TPaveText(0.903,0.400,0.936,0.428,"NDC");
	p_CSA_0->SetBorderSize(0);
	p_CSA_0->SetFillColor(kWhite);
	TText *t_CSA_0 = p_CSA_0->AddText(" 0^{o}");
	t_CSA_0->SetTextFont(nicefont);
	t_CSA_0->SetTextFont(nicefont);
	t_CSA_0->SetTextColor(3);
	
	p_CSA_0->Draw();
	
	
	TPaveText *p_CSA_30 = new TPaveText(0.903,0.428,0.936,0.456,"NDC");
	p_CSA_30->SetBorderSize(0);
	p_CSA_30->SetFillColor(kWhite);
	TText *t_CSA_30 = p_CSA_30->AddText("30^{o}");
	t_CSA_30->SetTextFont(nicefont);
	t_CSA_30->SetTextFont(nicefont);
	t_CSA_30->SetTextColor(3);
	
	p_CSA_30->Draw();
	
	
	
	TPaveText *p_CSA_45 = new TPaveText(0.903,0.456,0.936,0.484,"NDC");
	p_CSA_45->SetBorderSize(0);
	p_CSA_45->SetFillColor(kWhite);
	TText *t_CSA_45 = p_CSA_45->AddText("45^{o}");
	t_CSA_45->SetTextFont(nicefont);
	t_CSA_45->SetTextFont(nicefont);
	t_CSA_45->SetTextColor(3);
	
	p_CSA_45->Draw();
	
	

	gPad->SetLogy(1);
	C_CHD->SetGrid();


	
	TLegend *leg_CHD = new TLegend(0.50,0.14,0.85,0.35);
	leg_CHD->AddEntry(gr_CHD_S,"Scintillation Signal (PE)","l");
	leg_CHD->AddEntry(gr_CHD_PMT,"PMT Signal (fC)","l");
	leg_CHD->AddEntry(gr_CHD_CSA,"CSA Signal (mV)","l");
	leg_CHD->AddEntry(gr_CHD_ADC,"ADC Signal (counts)","l");

	
	/*
	leg_CHD->AddEntry(gr_CHD_PMT_m,"PMT Signal (fC)","lp");
	leg_CHD->AddEntry(gr_CHD_PMT_30_m,"PMT Signal (fC) 30^{o}","l");
	leg_CHD->AddEntry(gr_CHD_PMT_45_m,"PMT Signal (fC) 45^{o}","lp");
	leg_CHD->AddEntry(gr_CHD_CSA_m,"CSA Signal (mV)","lp");
	leg_CHD->AddEntry(gr_CHD_CSA_30_m,"CSA Signal (mV) 30^{o}","l");
	leg_CHD->AddEntry(gr_CHD_CSA_45_m,"CSA Signal (mV) 45^{o}","lp");
	leg_CHD->AddEntry(gr_CHD_ADC_m,"ADC Signal (counts)","lp");
	leg_CHD->AddEntry(gr_CHD_ADC_30_m,"ADC Signal (counts) 30^{o}","l");
	leg_CHD->AddEntry(gr_CHD_ADC_45_m,"ADC Signal (counts) 45^{o}","lp");
	*/
	
	leg_CHD->SetFillColor(10);
	leg_CHD->SetTextFont(nicefont);
	leg_CHD->SetBorderSize(0); // remove ugly border from stats 
	
	leg_CHD->Draw();
	
	
	
	
	
	
	
	
	
  	TH1F *h_S_min_abund_gf_2 = (TH1F*)h_S_min_abund_gf->Clone();
	h_S_min_abund_gf_2->SetName("h_S_min_abund_gf_2");
	h_S_min_abund_gf_2->SetTitle(t_string_abund);
	
  	TH1F *h_S_ave_abund_gf_2 = (TH1F*)h_S_ave_abund_gf->Clone();
	h_S_ave_abund_gf_2->SetName("h_S_ave_abund_gf_2");
	h_S_ave_abund_gf_2->SetTitle(t_string_abund);
  	
	TH1F *h_S_ave_abund_gf_3 = (TH1F*)h_S_ave_abund_gf->Clone();
	h_S_ave_abund_gf_3->SetName("h_S_ave_abund_gf_3");
	h_S_ave_abund_gf_3->SetTitle(t_string_abund);

	
	
  	TH1F *h_S_ave_abund_dgf_2 = (TH1F*)h_S_ave_abund_dgf->Clone();
	h_S_ave_abund_dgf_2->SetName("h_S_ave_abund_dgf_2");
	h_S_ave_abund_dgf_2->SetTitle(t_string_abund);
	
	
  	TH1F *h_S_ave_abund_dgf_3 = (TH1F*)h_S_ave_abund_dgf->Clone();
	h_S_ave_abund_dgf_3->SetName("h_S_ave_abund_dgf_3");
	h_S_ave_abund_dgf_3->SetTitle(t_string_abund);
	
  	TH1F *h_S_ave_abund_dgf_4 = (TH1F*)h_S_ave_abund_dgf->Clone();
	h_S_ave_abund_dgf_4->SetName("h_S_ave_abund_dgf_4");
	h_S_ave_abund_dgf_4->SetTitle(t_string_abund);
	
		//	TH1F *hnew = (TH1F*)h->Clone(); // renaming is recommended, because otherwise you
		//	hnew->SetName("hnew");
	
	TCanvas *C_h_abund = new TCanvas("C_h_abund","Impact of Solar Modulation on CALET 5 Year Abundances",0,0,1200,600);
	h_S_min_abund_dgf->GetXaxis()->SetRangeUser(25,60);
	h_S_min_abund_dgf->SetLineColor(4);
	h_S_min_abund_dgf->Draw("");
	h_S_max_abund_dgf->SetLineColor(2);
	h_S_max_abund_dgf->Draw("same");
	h_S_ave_abund_dgf_4->SetLineColor(1);
	h_S_ave_abund_dgf_4->Draw("same");
	gPad->SetLogy(1);
	C_h_abund->SetGrid();
	
	
	TLegend *leg_abund = new TLegend(0.157,0.724,0.407,0.874);
	leg_abund->AddEntry(h_S_max_abund_dgf,"Solar Maximum","l");
	leg_abund->AddEntry(h_S_min_abund_dgf,"Solar Minimum","l");
	leg_abund->AddEntry(h_S_ave_abund_dgf_4,"Average Modulation","l");
	
	
	leg_abund->SetFillColor(10);
	leg_abund->SetTextFont(nicefont);
	leg_abund->SetBorderSize(0); // remove ugly border from stats 
	leg_abund->SetTextSize(0.04);

	leg_abund->Draw();

	
	TCanvas *C_h_abund_dgf_comp = new TCanvas("C_h_abund_dgf_comp","Impact of Differential Geometry Factor CALET 5 Year Abundances",0,0,600,600);

	h_S_ave_abund_dgf->SetLineColor(2);
	h_S_ave_abund_dgf->Draw("");
	h_S_ave_abund_gf_3->Draw("same");
	gPad->SetLogy(1);
	C_h_abund_dgf_comp->SetGrid();
	

	
	TCanvas *C_h_abund_dgf = new TCanvas("C_h_abund_dgf","Impact of Interaction Correction for Differential Geometry Factor CALET 5 Year Abundances",0,0,600,600);
	
	h_S_ave_abund_dgf_2->SetLineColor(2);
	h_S_ave_abund_dgf_2->Draw("");
	h_S_ave_abund_dgf_noic->Draw("same");
	gPad->SetLogy(1);
	C_h_abund_dgf->SetGrid();
	
		// Comparison of CALET Average 5 Year Abundances before and after interaction corrections
	
		// Graphs of limiting maximum charges
	double x_CHD_PMT_45_Z_max_ideal[2];
	double y_CHD_PMT_45_Z_max_ideal[2];
	x_CHD_PMT_45_Z_max_ideal[0]=CHD_PMT_45_Z_max_ideal;
	x_CHD_PMT_45_Z_max_ideal[1]=CHD_PMT_45_Z_max_ideal;
	y_CHD_PMT_45_Z_max_ideal[0]=0;
	y_CHD_PMT_45_Z_max_ideal[1]=1e7;
	
	TGraph *gr_CHD_PMT_45_Z_max_ideal = new TGraph(2,x_CHD_PMT_45_Z_max_ideal,y_CHD_PMT_45_Z_max_ideal);
	
	
	double x_CHD_PMT_45_Z_max_reduced[2];
	double y_CHD_PMT_45_Z_max_reduced[2];
	x_CHD_PMT_45_Z_max_reduced[0]=CHD_PMT_45_Z_max_reduced;
	x_CHD_PMT_45_Z_max_reduced[1]=CHD_PMT_45_Z_max_reduced;
	y_CHD_PMT_45_Z_max_reduced[0]=0;
	y_CHD_PMT_45_Z_max_reduced[1]=1e7;
	
	TGraph *gr_CHD_PMT_45_Z_max_reduced = new TGraph(2,x_CHD_PMT_45_Z_max_reduced,y_CHD_PMT_45_Z_max_reduced);
	
	
	double x_CHD_CSA_45_Z_max[2];
	double y_CHD_CSA_45_Z_max[2];
	x_CHD_CSA_45_Z_max[0]=CHD_CSA_45_Z_max;
	x_CHD_CSA_45_Z_max[1]=CHD_CSA_45_Z_max;
	y_CHD_CSA_45_Z_max[0]=0;
	y_CHD_CSA_45_Z_max[1]=1e7;
	
	TGraph *gr_CHD_CSA_45_Z_max = new TGraph(2,x_CHD_CSA_45_Z_max,y_CHD_CSA_45_Z_max);
	
		
	double x_CHD_ADC_45_Z_max[2];
	double y_CHD_ADC_45_Z_max[2];
	x_CHD_ADC_45_Z_max[0]=CHD_ADC_45_Z_max;
	x_CHD_ADC_45_Z_max[1]=CHD_ADC_45_Z_max;
	y_CHD_ADC_45_Z_max[0]=0;
	y_CHD_ADC_45_Z_max[1]=1e7;
	
	TGraph *gr_CHD_ADC_45_Z_max = new TGraph(2,x_CHD_ADC_45_Z_max,y_CHD_ADC_45_Z_max);
	
	

	
	
	
	
	TCanvas *C_h_abund_ave = new TCanvas("C_h_abund_ave","CALET Average 5 Year Abundances Before and After Interaction Corrections",0,0,1200,600);
		//	h_S_ave_abund->Draw("");
	h_S_ave_abund_dgf_noic->GetXaxis()->SetRangeUser(25,60);
	h_S_ave_abund_dgf_noic->Draw("");
	h_S_ave_abund_dgf_3->SetLineColor(2);
	h_S_ave_abund_dgf_3->Draw("same");
	
	
	gr_CHD_PMT_45_Z_max_reduced->SetLineColor(2);
	gr_CHD_PMT_45_Z_max_reduced->SetLineStyle(2);
	gr_CHD_PMT_45_Z_max_reduced->SetLineWidth(2);
	gr_CHD_PMT_45_Z_max_reduced->Draw("L");
	gr_CHD_PMT_45_Z_max_ideal->SetLineColor(2);
	gr_CHD_PMT_45_Z_max_ideal->SetLineStyle(2);
	gr_CHD_PMT_45_Z_max_ideal->SetLineWidth(2);
	gr_CHD_PMT_45_Z_max_ideal->Draw("L");
	gr_CHD_CSA_45_Z_max->SetLineColor(3);
	gr_CHD_CSA_45_Z_max->SetLineStyle(2);
	gr_CHD_CSA_45_Z_max->SetLineWidth(2);
	gr_CHD_CSA_45_Z_max->Draw("L");
	gr_CHD_ADC_45_Z_max->SetLineColor(4);
	gr_CHD_ADC_45_Z_max->SetLineStyle(2);
	gr_CHD_ADC_45_Z_max->SetLineWidth(2);
	gr_CHD_ADC_45_Z_max->Draw("L");
	
	
	gPad->SetLogy(1);
	C_h_abund_ave->SetGrid();
	
	TLegend *leg_ave_abund = new TLegend(0.157,0.724,0.407,0.874);
	leg_ave_abund->AddEntry(h_S_ave_abund_dgf_noic,l_string_S_ave_abund_dgf_noic,"l");
	leg_ave_abund->AddEntry(h_S_ave_abund_dgf_3,l_string_S_ave_abund_dgf,"l");
	leg_ave_abund->AddEntry(gr_CHD_PMT_max_reduced,"PMT Limit at 45^{o}","l");
	leg_ave_abund->AddEntry(gr_CHD_CSA_max,"CSA Limit at 45^{o}","l");
	leg_ave_abund->AddEntry(gr_CHD_ADC_max,"ADC Limit at 45^{o}","l");

	
	leg_ave_abund->SetFillColor(10);
	leg_ave_abund->SetTextFont(nicefont);
	leg_ave_abund->SetBorderSize(0); // remove ugly border from stats 
	leg_ave_abund->SetTextSize(0.04);
	
	leg_ave_abund->Draw();
	
		// Plot comparing anticipated CALET 5 year abundances with TIGER results at instrument

	double TI_Z[30];
	double TI_abund[30];
	double TI_abund_r[30];
	double x_error_low[30];
	double x_error_high[30];
	double TI_err_lo_r[30];
	double TI_err_lo[30];
	double TI_err_hi_r[30];
	double TI_err_hi[30];
	double dummy1[30];
	double dummy2[30];

	double TI_Fe_count=564627.12;
	int TI_index=0;
	
	char tiger_instrument[200] ="/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/st_tigeriss_2022.txt";
	
	in_1.open(tiger_instrument);
	while(1) {
		in_1 >> TI_Z[TI_index] >> TI_abund_r[TI_index] >> TI_err_lo_r[TI_index] >> TI_err_hi_r[TI_index] >> dummy1[TI_index] >> dummy2[TI_index];
		x_error_low[TI_index]=0;
		x_error_high[TI_index]=0;
		TI_abund[TI_index]=TI_abund_r[TI_index]*TI_Fe_count;
		TI_err_lo[TI_index]=TI_err_lo_r[TI_index]*TI_Fe_count;
		TI_err_hi[TI_index]=TI_err_hi_r[TI_index]*TI_Fe_count;
		if(in_1.eof() || !in_1.good()) break;
		TI_index++;
	}
	in_1.close();
	in_1.clear();
	
	for (i=0; i<100; i++) {
		for (j=0; j<TI_index; j++) {
			if (TI_Z[j]==el_dat[i].Z) {
				fprintf(stderr,"%d %e %e %.2f\n",el_dat[i].Z,CALET_abund[i],TI_abund[j],CALET_abund[i]/TI_abund[j]);
			}
		}
	}
	
	/*
	TGraphAsymmErrors *gr_TIGER_instrument = 
    new TGraphAsymmErrors(TI_index,TI_Z,TI_abund,
						  x_error_low,x_error_high,
						  TI_err_lo,TI_err_hi);
	 */
	
	TGraph *gr_TIGER_instrument = new TGraph(TI_index,TI_Z,TI_abund);
	
	
	
	TCanvas *C_h_abund_ave_comp = new TCanvas("C_h_abund_ave_comp","CALET Average 5 Year Abundances Compared to TIGER",0,0,800,400);
		//	h_S_ave_abund->Draw("");
	h_S_ave_abund_dgf->GetXaxis()->SetRangeUser(25,40);
	h_S_ave_abund_dgf->SetLineColor(1);
	h_S_ave_abund_dgf->SetTitle("");
	
	h_S_ave_abund_dgf->Draw("");
	
	gr_TIGER_instrument->SetMarkerStyle(21);
	gr_TIGER_instrument->Draw("P");
	
	gPad->SetLogy(1);
	C_h_abund_ave_comp->SetGrid();

	
	TLegend *leg_abund_ave_comp = new TLegend(0.384,0.715,0.869,0.861);
	leg_abund_ave_comp->AddEntry(h_S_ave_abund_gf_2,"CALET Average 5 Years","l");
	leg_abund_ave_comp->AddEntry(gr_TIGER_instrument,"TIGER Instrument","p");
	
	
	leg_abund_ave_comp->SetFillColor(10);
	leg_abund_ave_comp->SetTextFont(nicefont);
	leg_abund_ave_comp->SetBorderSize(0); // remove ugly border from stats 
	leg_abund_ave_comp->SetTextSize(0.04);
	
	leg_abund_ave_comp->Draw();
	
    //in_2.close();
    //in_2.clear();

//	in_2.open(abund_smax);
    
    in_2.open("test.txt");

    // Try using C I/O here!
    
	ACE_Z_index_smax=0;

    fprintf(stderr,"ACE-CRIS SMAX Abundances\n");

    while(1) {
        in_2 >> A_max[ACE_Z_index_smax].Z >> A_max[ACE_Z_index_smax].I >> A_max[ACE_Z_index_smax].pu >> A_max[ACE_Z_index_smax].s_stat >> A_max[ACE_Z_index_smax].s_sys;
        
		if(A_max[ACE_Z_index_smax].Z == 26) ref_index=ACE_Z_index_smax;
        
		if(in_2.eof() || !in_2.good()) break;
        fprintf(stderr,"%d\t%e\t%f\t%f\t%f\n",A_max[ACE_Z_index_smax].Z,A_max[ACE_Z_index_smax].I,A_max[ACE_Z_index_smax].pu,A_max[ACE_Z_index_smax].s_stat,A_max[ACE_Z_index_smax].s_sys);
		ACE_Z_index_smax++;
	}
    
    fprintf(stderr,"ACE_Z_index_smax: %d\n",ACE_Z_index_smax);
    
	in_2.close();
	in_2.clear();
    

    
	in_2.open(abund_smin);
    
	ACE_Z_index_smin=0;
    
    fprintf(stderr,"ACE-CRIS SMIN Abundances\n");
    
	while(1) {
		in_2 >> A_min[ACE_Z_index_smin].Z >> A_min[ACE_Z_index_smin].I >> A_min[ACE_Z_index_smin].pu >> A_min[ACE_Z_index_smin].s_stat >> A_min[ACE_Z_index_smin].s_sys;
        
        
		if(A_min[ACE_Z_index_smin].Z == 26) ref_index=ACE_Z_index_smin;
        
		if(in_2.eof() || !in_2.good()) break;
        fprintf(stderr,"%d\t%e\t%f\t%f\t%f\n",A_min[ACE_Z_index_smin].Z,A_min[ACE_Z_index_smin].I,A_min[ACE_Z_index_smin].pu,A_min[ACE_Z_index_smin].s_stat,A_min[ACE_Z_index_smin].s_sys);
		ACE_Z_index_smin++;
	}
    
	in_2.close();
	in_2.clear();

	for(i=0;i<ACE_Z_index_smin;i++)
		{
			A_min[i].rel_a=A_min[i].I/A_min[ref_index].I;
			A_max[i].rel_a=A_max[i].I/A_max[ref_index].I;
            fprintf(stderr,"%d %f %f\n",A_min[i].Z,A_min[i].rel_a,A_max[i].rel_a);
		}
	fprintf(stderr,"ACE_Z_index_smin: %d\n",ACE_Z_index_smin);


  for(Z_index=5;Z_index<33;Z_index++)
		{
		S_max[Z_index].rel_ab_I=S_max[Z_index].I_spec/S_max[26].I_spec;
		S_min[Z_index].rel_ab_I=S_min[Z_index].I_spec/S_min[26].I_spec;
		}


	fprintf(stderr," Z Smin_rel_E Smax_rel_E Smin_rel_I Smax_rel_I  c_min   i_min  c_max   i_max\n");
	for(i=0;i<28;i++)
		{
		if(i<ACE_Z_index_smin)
			{
			fprintf(stderr,"%2d %2d %.4e %.4e %.4e %.4e %.2e %.2f %.2e %.2f\n",i+5,A_min[i].Z,A_min[i].rel_a,A_max[i].rel_a,S_min[i+5].rel_ab_I,S_max[i+5].rel_ab_I,S_min[i+5].c_pl,S_min[i+5].i_pl,S_max[i+5].c_pl,S_max[i+5].i_pl);
			}
		else
			{
			fprintf(stderr,"%2d %.4e %.4e %.2e %.2f %.2e %.2f\n",i+5,S_min[i+5].rel_ab_I,S_max[i+5].rel_ab_I,S_min[i+5].c_pl,
							S_min[i+5].i_pl,S_max[i+5].c_pl,S_max[i+5].i_pl);
			}
		}


   // Solar Max Plot Label for relative abundances taken at 160 MeV/nuc
  
  double label_x_max_abund[2];
  double label_y_max_abund[2];
  label_x_max_abund[0] = 0;
  label_x_max_abund[1] = 30;
  label_y_max_abund[0] = 1e-3;
  label_y_max_abund[1] = 2e1;

  TGraph *gr_label_S_max_abund = new TGraph(2,label_x_max_abund,label_y_max_abund);

  gr_label_S_max_abund->SetMarkerStyle(1);
  gr_label_S_max_abund->SetMarkerSize(1);
  gr_label_S_max_abund->SetTitle("ACE-CRIS Solar Maximum Abundances");
  gr_label_S_max_abund->GetXaxis()->SetTitle("Z");
  gr_label_S_max_abund->GetXaxis()->CenterTitle();
  gr_label_S_max_abund->GetYaxis()->SetTitle("Relative Abundance (Fe=1)");
  gr_label_S_max_abund->GetYaxis()->CenterTitle();
  gr_label_S_max_abund->GetXaxis()->SetTitleOffset(1.2);
  gr_label_S_max_abund->GetYaxis()->SetTitleOffset(1.4);
  gr_label_S_max_abund->GetYaxis()->SetLabelFont(nicefont);
  gr_label_S_max_abund->GetYaxis()->SetTitleFont(nicefont);
  gr_label_S_max_abund->GetXaxis()->SetLabelFont(nicefont);
  gr_label_S_max_abund->GetXaxis()->SetTitleFont(nicefont);
  gr_label_S_max_abund->GetXaxis()->SetRangeUser(5,30);
  gr_label_S_max_abund->GetYaxis()->SetRangeUser(1e-3,2e1);

	plot_count=0;
	for(i=0;i<ACE_Z_index_smax;i++) {
 	  x_plot[i]=A_max[i].Z;
 	  y_plot[i]=A_max[i].rel_a;
		plot_count++;
		}
  TGraph *gr_S_max_abund = new TGraph(plot_count,x_plot,y_plot);
  TGraph *gr_S_max_abund_2 = new TGraph(plot_count,x_plot,y_plot);
  
  TCanvas *C_S_max_abund = new TCanvas("C_S_max_abund","ACE-CRIS S_Max Abundances",0,0,800,800);
  
  C_S_max_abund->SetLeftMargin(0.124365);
  C_S_max_abund->SetRightMargin(0.0748731);
  C_S_max_abund->SetBottomMargin(0.117801);
  C_S_max_abund->SetTopMargin(0.0824607);

  C_S_max_abund->SetGrid();
  
  // gPad->SetLogx(1);
  gPad->SetLogy(1);
  
  gr_label_S_max_abund->Draw("AP");

	gr_S_max_abund->Draw("PL");

   // Solar Min Plot Label for relative abundances taken at 160 MeV/nuc
  
  double label_x_min_abund[2];
  double label_y_min_abund[2];
  label_x_min_abund[0] = 0;
  label_x_min_abund[1] = 30;
  label_y_min_abund[0] = 1e-3;
  label_y_min_abund[1] = 2e1;

  TGraph *gr_label_S_min_abund = new TGraph(2,label_x_min_abund,label_y_min_abund);

  gr_label_S_min_abund->SetMarkerStyle(1);
  gr_label_S_min_abund->SetMarkerSize(1);
  gr_label_S_min_abund->SetTitle("ACE-CRIS Solar Minimum Abundances");
  gr_label_S_min_abund->GetXaxis()->SetTitle("Z");
  gr_label_S_min_abund->GetXaxis()->CenterTitle();
  gr_label_S_min_abund->GetYaxis()->SetTitle("Relative Abundance (Fe=1)");
  gr_label_S_min_abund->GetYaxis()->CenterTitle();
  gr_label_S_min_abund->GetXaxis()->SetTitleOffset(1.2);
  gr_label_S_min_abund->GetYaxis()->SetTitleOffset(1.4);
  gr_label_S_min_abund->GetYaxis()->SetLabelFont(nicefont);
  gr_label_S_min_abund->GetYaxis()->SetTitleFont(nicefont);
  gr_label_S_min_abund->GetXaxis()->SetLabelFont(nicefont);
  gr_label_S_min_abund->GetXaxis()->SetTitleFont(nicefont);
  gr_label_S_min_abund->GetXaxis()->SetRangeUser(5,30);
  gr_label_S_min_abund->GetYaxis()->SetRangeUser(1e-3,2e1);

    fprintf(stderr,"ACE Relative Abundance!!!\n");
	plot_count=0;
	for(i=0;i<ACE_Z_index_smin;i++) {
 	  x_plot[i]=A_min[i].Z;
 	  y_plot[i]=A_min[i].rel_a;
        fprintf(stderr,"%f %f\n",x_plot[i],y_plot[i]);
		plot_count++;
		}
  TGraph *gr_S_min_abund = new TGraph(plot_count,x_plot,y_plot);
  TGraph *gr_S_min_abund_2 = new TGraph(plot_count,x_plot,y_plot);
  
  TCanvas *C_S_min_abund = new TCanvas("C_S_min_abund","ACE-CRIS S_Min Abundances",0,0,800,800);
  
  C_S_min_abund->SetLeftMargin(0.124365);
  C_S_min_abund->SetRightMargin(0.0748731);
  C_S_min_abund->SetBottomMargin(0.117801);
  C_S_min_abund->SetTopMargin(0.0824607);

  C_S_min_abund->SetGrid();
  
  // gPad->SetLogx(1);
  gPad->SetLogy(1);
  
  gr_label_S_min_abund->Draw("AP");

	gr_S_min_abund->Draw("PL");


   // Solar Comp Plot Label for relative abundances taken at 160 MeV/nuc
  
  double label_x_comp_abund[2];
  double label_y_comp_abund[2];
  label_x_comp_abund[0] = 0;
  label_x_comp_abund[1] = 30;
  label_y_comp_abund[0] = 1e-3;
  label_y_comp_abund[1] = 2e1;

  TGraph *gr_label_comp_abund = new TGraph(2,label_x_comp_abund,label_y_comp_abund);

  gr_label_comp_abund->SetMarkerStyle(1);
  gr_label_comp_abund->SetMarkerSize(1);
  gr_label_comp_abund->SetTitle("Comparison of ACE-CRIS Abundances for Solar Minimum and Maximum");
  gr_label_comp_abund->GetXaxis()->SetTitle("Z");
  gr_label_comp_abund->GetXaxis()->CenterTitle();
  gr_label_comp_abund->GetYaxis()->SetTitle("Relative Abundance (Fe=1)");
  gr_label_comp_abund->GetYaxis()->CenterTitle();
  gr_label_comp_abund->GetXaxis()->SetTitleOffset(1.2);
  gr_label_comp_abund->GetYaxis()->SetTitleOffset(1.4);
  gr_label_comp_abund->GetYaxis()->SetLabelFont(nicefont);
  gr_label_comp_abund->GetYaxis()->SetTitleFont(nicefont);
  gr_label_comp_abund->GetXaxis()->SetLabelFont(nicefont);
  gr_label_comp_abund->GetXaxis()->SetTitleFont(nicefont);
  gr_label_comp_abund->GetXaxis()->SetRangeUser(5,30);
  gr_label_comp_abund->GetYaxis()->SetRangeUser(1e-3,2e1);

  TCanvas *C_comp_abund = new TCanvas("C_comp_abund","ACE-CRIS Comparison Abundances",0,0,800,800);
  
  C_comp_abund->SetLeftMargin(0.124365);
  C_comp_abund->SetRightMargin(0.0748731);
  C_comp_abund->SetBottomMargin(0.117801);
  C_comp_abund->SetTopMargin(0.0824607);

  C_comp_abund->SetGrid();
  
  // gPad->SetLogx(1);
  gPad->SetLogy(1);
  
  gr_label_comp_abund->Draw("AP");
	gr_S_min_abund_2->SetLineColor(2);
	gr_S_max_abund_2->SetLineColor(4);
	gr_S_min_abund_2->Draw("PL");
	gr_S_max_abund_2->Draw("PL");

  TLegend *leg_comp_abund = new TLegend(0.60,0.73,0.87,0.82);
  leg_comp_abund->AddEntry(gr_S_min_abund_2,"Solar Minimum","l");
  leg_comp_abund->AddEntry(gr_S_max_abund_2,"Solar Maximum","l");

  leg_comp_abund->SetFillColor(10);
  leg_comp_abund->SetTextFont(nicefont);
  leg_comp_abund->SetBorderSize(0); // remove ugly border from stats 

	leg_comp_abund->Draw();

	plot_count=0;
	for(i=5;i<33;i++) {
 	  x_plot[i-5]=i;
 	  y_plot[i-5]=S_max[i].rel_ab_I;
		plot_count++;
		}
  TGraph *gr_S_max_abund_I = new TGraph(plot_count,x_plot,y_plot);

	plot_count=0;
	for(i=5;i<33;i++) {
 	  x_plot[i-5]=i;
 	  y_plot[i-5]=S_min[i].rel_ab_I;
		plot_count++;
		}
  TGraph *gr_S_min_abund_I = new TGraph(plot_count,x_plot,y_plot);


   // Solar Comp Plot Label for relative abundances taken at 160 MeV/nuc and integrated from spectra
  
  double label_x_comp_abund_2[2];
  double label_y_comp_abund_2[2];
  label_x_comp_abund_2[0] = 0;
  label_x_comp_abund_2[1] = 35;
  label_y_comp_abund_2[0] = 1e-5;
  label_y_comp_abund_2[1] = 2e1;

  TGraph *gr_label_comp_abund_2 = new TGraph(2,label_x_comp_abund_2,label_y_comp_abund_2);

  gr_label_comp_abund_2->SetMarkerStyle(1);
  gr_label_comp_abund_2->SetMarkerSize(1);
  gr_label_comp_abund_2->SetTitle("Comparison of ACE-CRIS Abundances for Solar Minimum and Maximum");
  gr_label_comp_abund_2->GetXaxis()->SetTitle("Z");
  gr_label_comp_abund_2->GetXaxis()->CenterTitle();
  gr_label_comp_abund_2->GetYaxis()->SetTitle("Relative Abundance (Fe=1)");
  gr_label_comp_abund_2->GetYaxis()->CenterTitle();
  gr_label_comp_abund_2->GetXaxis()->SetTitleOffset(1.2);
  gr_label_comp_abund_2->GetYaxis()->SetTitleOffset(1.4);
  gr_label_comp_abund_2->GetYaxis()->SetLabelFont(nicefont);
  gr_label_comp_abund_2->GetYaxis()->SetTitleFont(nicefont);
  gr_label_comp_abund_2->GetXaxis()->SetLabelFont(nicefont);
  gr_label_comp_abund_2->GetXaxis()->SetTitleFont(nicefont);
  gr_label_comp_abund_2->GetXaxis()->SetRangeUser(5,35);
  gr_label_comp_abund_2->GetYaxis()->SetRangeUser(1e-5,2e1);

  TCanvas *C_comp_abund_2 = new TCanvas("C_comp_abund_2","ACE-CRIS Comparison Abundances",0,0,800,800);
  
  C_comp_abund_2->SetLeftMargin(0.124365);
  C_comp_abund_2->SetRightMargin(0.0748731);
  C_comp_abund_2->SetBottomMargin(0.117801);
  C_comp_abund_2->SetTopMargin(0.0824607);

  C_comp_abund_2->SetGrid();
  
  // gPad->SetLogx(1);
  gPad->SetLogy(1);
  
  gr_label_comp_abund_2->Draw("AP");
	gr_S_min_abund_2->SetLineColor(2);
	gr_S_max_abund_2->SetLineColor(3);
	gr_S_min_abund_I->SetLineColor(4);
	gr_S_max_abund_I->SetLineColor(6);
	gr_S_min_abund_2->Draw("PL");
	gr_S_max_abund_2->Draw("PL");
	gr_S_min_abund_I->Draw("PL");
	gr_S_max_abund_I->Draw("PL");


  TLegend *leg_comp_abund_2 = new TLegend(0.60,0.78,0.87,0.87);
  leg_comp_abund_2->AddEntry(gr_S_min_abund_2,"Solar Minimum 160 MeV/nuc","l");
  leg_comp_abund_2->AddEntry(gr_S_max_abund_2,"Solar Maximum 160 MeV/nuc","l");
  leg_comp_abund_2->AddEntry(gr_S_min_abund_I,"Solar Minimum Integral","l");
  leg_comp_abund_2->AddEntry(gr_S_max_abund_I,"Solar Maximum Integral","l");

  leg_comp_abund_2->SetFillColor(10);
  leg_comp_abund_2->SetTextFont(nicefont);
  leg_comp_abund_2->SetBorderSize(0); // remove ugly border from stats 

	leg_comp_abund_2->Draw();

		// ISS Orbit
	int flip_count=0;
	double rate=10;
	int ISS_max=40500; // 20100
	int ISS_count=0;
	double omega_E=1.0/15.74;
	double long_plot[ISS_max];
	double lat_plot[ISS_max];
	double long_calc;
	double long_calc_last=0.0;
	double lat_calc;
	double lat_calc_2;
	double lat_calc_last=0.0;
	double omega=0.0;
	double phase=TMath::Pi();
	int orbit_count=0;
	for (ISS_count=0; ISS_count<ISS_max; ISS_count++) {
		omega=(TMath::Pi()/180.0)*ISS_count/rate;
		long_calc=(180.0/TMath::Pi())*atan(sin(omega-phase)/(cos(omega-phase)*cos(ISS_angle_rad))) + 
		omega_E*ISS_count/rate;
		lat_calc_2=(180.0/TMath::Pi())*atan(cos(omega-phase)*sin(ISS_angle_rad))/
		sqrt(pow(cos(omega-phase)*cos(ISS_angle_rad),2.0) + pow(sin(omega-phase),2.0));
		lat_calc=(180.0/TMath::Pi())*atan2(cos(omega-phase)*sin(ISS_angle_rad),
										   sqrt(pow(cos(omega-phase)*cos(ISS_angle_rad),2.0) 
																				+ pow(sin(omega-phase),2.0)));
		if (lat_calc < 0 && lat_calc_last >= 0) {
			flip_count++;
		}
		if (lat_calc > 0 && lat_calc_last <= 0) {
			flip_count++;
		}
		long_calc+=180.0*flip_count;
		
		if (long_calc - 360*orbit_count>=180.0) {
			orbit_count++;
		}

			//fprintf(stderr,"lat: %.2f %.2f long: %.2f\n",lat_calc,lat_calc_2,long_calc);
		long_plot[ISS_count]=long_calc-360*orbit_count;
		lat_plot[ISS_count]=lat_calc;
		long_calc_last=long_calc;
		lat_calc_last=lat_calc;
	}
	
	TGraph *gr_ISS_orbit = new TGraph(ISS_max,long_plot,lat_plot);
	
	gr_ISS_orbit->SetMarkerStyle(1);
	gr_ISS_orbit->SetMarkerSize(1);
	gr_ISS_orbit->SetTitle("ISS Orbit");
	gr_ISS_orbit->GetXaxis()->SetTitle("Longitude (degrees)");
	gr_ISS_orbit->GetXaxis()->CenterTitle();
	gr_ISS_orbit->GetYaxis()->SetTitle("Relative Abundance (Fe=1)");
	gr_ISS_orbit->GetYaxis()->CenterTitle();
	gr_ISS_orbit->GetXaxis()->SetTitleOffset(1.2);
	gr_ISS_orbit->GetYaxis()->SetTitleOffset(1.4);
	gr_ISS_orbit->GetYaxis()->SetLabelFont(nicefont);
	gr_ISS_orbit->GetYaxis()->SetTitleFont(nicefont);
	gr_ISS_orbit->GetXaxis()->SetLabelFont(nicefont);
	gr_ISS_orbit->GetXaxis()->SetTitleFont(nicefont);
		//	gr_ISS_orbit->GetXaxis()->SetRangeUser(5,35);
		//	gr_ISS_orbit->GetYaxis()->SetRangeUser(1e-5,2e1);

		// Contour Map	
	
	TGraph2D *g_rigidity_cutoff = 
    new TGraph2D(c_count,x_cont,y_cont,z_cont);
	
		//	g_rigidity_cutoff->SetTitle("Geomagnetic Vertical Cutoff Rigidities");
	g_rigidity_cutoff->SetTitle("");
	g_rigidity_cutoff->GetXaxis()->SetTitle("Longitude (degrees)");
	g_rigidity_cutoff->GetXaxis()->CenterTitle();
	g_rigidity_cutoff->GetXaxis()->SetTitleOffset(1.2);
	g_rigidity_cutoff->GetYaxis()->SetTitle("Latitude (degrees)");
	g_rigidity_cutoff->GetYaxis()->CenterTitle();
	g_rigidity_cutoff->GetYaxis()->SetTitleOffset(1.2);
	g_rigidity_cutoff->GetYaxis()->SetLabelFont(nicefont);
	g_rigidity_cutoff->GetYaxis()->SetTitleFont(nicefont);
	g_rigidity_cutoff->GetXaxis()->SetLabelFont(nicefont);
	g_rigidity_cutoff->GetXaxis()->SetTitleFont(nicefont);
	g_rigidity_cutoff->GetZaxis()->SetTitle("Rigidity (GV)");
	g_rigidity_cutoff->GetZaxis()->SetTitleFont(nicefont);
	g_rigidity_cutoff->GetZaxis()->SetTitleOffset(0.7);
	g_rigidity_cutoff->GetZaxis()->CenterTitle();
		//g_rigidity_cutoff->GetZaxis()->SetTickSize(0.5);

	g_rigidity_cutoff->GetZaxis()->SetRangeUser(-0.1,16);
	
	
	
	
	TH2F *hm = new TH2F("hm","Mercator",  180, -180, 180, 161, -80.5, 80.5);
	
	
	double x;
	double y;
	double x_m[22000];
	double y_m[22000];
	
	int xy_count=0;
	
	char earth_file[100] = "/usr/local/root-v5-34-00-patches/tutorials/graphics/earth.dat";
	
	in_1.open(earth_file);
	
	while (1) {
		in_1 >> x >> y;
		if (!in_1.good()) break;
		hm->Fill(x,y, 1);
		x_m[xy_count]=x;
		y_m[xy_count]=y;
		xy_count++;	
	}
	in_1.close();
	in_1.clear();
	
	
	
	TGraph *gr_mercator = new TGraph(xy_count,x_m,y_m);
	
	gr_mercator->SetMarkerStyle(1);
	gr_mercator->SetMarkerSize(1);
	gr_mercator->SetTitle("Magnetic Cutoff Rigidities");
	gr_mercator->GetXaxis()->SetTitle("Longitude (degrees)");
	gr_mercator->GetXaxis()->CenterTitle();
	gr_mercator->GetYaxis()->SetTitle("Latitude (degrees)");
	gr_mercator->GetYaxis()->CenterTitle();
	gr_mercator->GetXaxis()->SetTitleOffset(1.2);
	gr_mercator->GetYaxis()->SetTitleOffset(1.2);
	gr_mercator->GetYaxis()->SetLabelFont(nicefont);
	gr_mercator->GetYaxis()->SetTitleFont(nicefont);
	gr_mercator->GetXaxis()->SetLabelFont(nicefont);
	gr_mercator->GetXaxis()->SetTitleFont(nicefont);
	gr_mercator->GetXaxis()->SetRangeUser(-180,180);
	gr_mercator->GetYaxis()->SetRangeUser(-85,85);
	
	
	
	
	for (i=0; i<long_bins_calc; i++) {
		for (j=0; j<lat_bins; j++) {
			
			if (alt_450[j+1].cut_off[i] > el_dat[plot_Z-1].R_cut && alt_450[j].cut_off[i] <= el_dat[plot_Z-1].R_cut && alt_450[j].latitude > 10) {
				northern_hemisphere_cut[i]=alt_450[j+1].latitude + (alt_450[j].latitude-alt_450[j+1].latitude)*
				(alt_450[j].cut_off[i]-el_dat[plot_Z-1].R_cut)/(alt_450[j].cut_off[i]-alt_450[j+1].cut_off[i]); 
			}
			if (alt_450[j+1].cut_off[i] < el_dat[plot_Z-1].R_cut && alt_450[j].cut_off[i] >= el_dat[plot_Z-1].R_cut && alt_450[j].latitude < -10) {
				southern_hemisphere_cut[i]=alt_450[j+1].latitude + (alt_450[j].latitude-alt_450[j+1].latitude)*
				(alt_450[j].cut_off[i]-el_dat[plot_Z-1].R_cut)/(alt_450[j].cut_off[i]-alt_450[j+1].cut_off[i]);
			}
		}
	}
	double northern_hemisphere_cut_x[long_bins_calc];
	double northern_hemisphere_cut_y[long_bins_calc];
	double southern_hemisphere_cut_x[long_bins_calc];
	double southern_hemisphere_cut_y[long_bins_calc];
	for (i=0; i<long_bins_calc; i++) {
		if (longitude[i]< 180) {
		northern_hemisphere_cut_x[i+long_bins_calc/2]=longitude[i];
		northern_hemisphere_cut_y[i+long_bins_calc/2]=northern_hemisphere_cut[i];
		}
		else if (longitude[i]>= 180) {
			northern_hemisphere_cut_x[i-long_bins_calc/2]=longitude[i]-360;
			northern_hemisphere_cut_y[i-long_bins_calc/2]=northern_hemisphere_cut[i];
		}
		if (longitude[i]< 180) {
			southern_hemisphere_cut_x[i+long_bins_calc/2]=longitude[i];
			southern_hemisphere_cut_y[i+long_bins_calc/2]=southern_hemisphere_cut[i];
		}
		else if (longitude[i]>= 180) {
			southern_hemisphere_cut_x[i-long_bins_calc/2]=longitude[i]-360;
			southern_hemisphere_cut_y[i-long_bins_calc/2]=southern_hemisphere_cut[i];
		}
		fprintf(stderr,"northern_hemisphere_cut_x[%2d]: %7.2f northern_hemisphere_cut_y[%2d]: %7.2f southern_hemisphere_cut_x[%2d]: %7.2f southern_hemisphere_cut_y[%2d]: %7.2f\n"
				,i,northern_hemisphere_cut_x[i],i,northern_hemisphere_cut_y[i]
				,i,southern_hemisphere_cut_x[i],i,southern_hemisphere_cut_y[i]);
	}
	
	TGraph *gr_northern_hemisphere_cut = new TGraph(long_bins_calc,northern_hemisphere_cut_x,northern_hemisphere_cut_y);
	
	gr_northern_hemisphere_cut->SetMarkerStyle(1);
	gr_northern_hemisphere_cut->SetMarkerSize(1);
	gr_northern_hemisphere_cut->SetTitle("Northern Rigidity Threshold");
	gr_northern_hemisphere_cut->GetXaxis()->SetTitle("Longitude (degrees)");
	gr_northern_hemisphere_cut->GetXaxis()->CenterTitle();
	gr_northern_hemisphere_cut->GetYaxis()->SetTitle("Latitude (degrees)");
	gr_northern_hemisphere_cut->GetYaxis()->CenterTitle();
	gr_northern_hemisphere_cut->GetXaxis()->SetTitleOffset(1.2);
	gr_northern_hemisphere_cut->GetYaxis()->SetTitleOffset(1.2);
	gr_northern_hemisphere_cut->GetYaxis()->SetLabelFont(nicefont);
	gr_northern_hemisphere_cut->GetYaxis()->SetTitleFont(nicefont);
	gr_northern_hemisphere_cut->GetXaxis()->SetLabelFont(nicefont);
	gr_northern_hemisphere_cut->GetXaxis()->SetTitleFont(nicefont);
	gr_northern_hemisphere_cut->GetXaxis()->SetRangeUser(-180,180);
	gr_northern_hemisphere_cut->GetYaxis()->SetRangeUser(-85,85);
	
	
	TGraph *gr_southern_hemisphere_cut = new TGraph(long_bins_calc,southern_hemisphere_cut_x,southern_hemisphere_cut_y);
	
	gr_southern_hemisphere_cut->SetMarkerStyle(1);
	gr_southern_hemisphere_cut->SetMarkerSize(1);
	gr_southern_hemisphere_cut->SetTitle("Southern Rigidity Threshold");
	gr_southern_hemisphere_cut->GetXaxis()->SetTitle("Longitude (degrees)");
	gr_southern_hemisphere_cut->GetXaxis()->CenterTitle();
	gr_southern_hemisphere_cut->GetYaxis()->SetTitle("Latitude (degrees)");
	gr_southern_hemisphere_cut->GetYaxis()->CenterTitle();
	gr_southern_hemisphere_cut->GetXaxis()->SetTitleOffset(1.2);
	gr_southern_hemisphere_cut->GetYaxis()->SetTitleOffset(1.2);
	gr_southern_hemisphere_cut->GetYaxis()->SetLabelFont(nicefont);
	gr_southern_hemisphere_cut->GetYaxis()->SetTitleFont(nicefont);
	gr_southern_hemisphere_cut->GetXaxis()->SetLabelFont(nicefont);
	gr_southern_hemisphere_cut->GetXaxis()->SetTitleFont(nicefont);
	gr_southern_hemisphere_cut->GetXaxis()->SetRangeUser(-180,180);
	gr_southern_hemisphere_cut->GetYaxis()->SetRangeUser(-85,85);
	
	
	TPaletteAxis *p = (TPaletteAxis*)g_rigidity_cutoff->FindObject("palette");
	
		//	c1->Update();
	
	/*
	p->SetX1NDC(0.90);
	p->SetX2NDC(0.94);
	p->SetY1NDC(0.10);
	p->SetY2NDC(0.90);
	*/
		//c1->Modified();
	
	TCanvas *c1 = new TCanvas("c1","earth_projections",1200,800);

		//	hm->Draw("contz,same");
	g_rigidity_cutoff->Draw("COLZ");
	gStyle->SetPalette(1);
	gr_mercator->Draw("P");
	gr_mercator->SetMarkerStyle(1);
	gr_ISS_orbit->SetMarkerColor(1);
	gr_ISS_orbit->SetMarkerStyle(6);
	gr_ISS_orbit->SetLineWidth(2);
	gr_ISS_orbit->Draw("P");
	gr_northern_hemisphere_cut->SetLineStyle(2);
	gr_northern_hemisphere_cut->SetLineWidth(4);
	gr_northern_hemisphere_cut->Draw("C");
	gr_southern_hemisphere_cut->SetLineStyle(2);
	gr_southern_hemisphere_cut->SetLineWidth(4);
	gr_southern_hemisphere_cut->Draw("C");

	/*
	C_Evsbeta->SetLeftMargin(0.148477);
	C_Evsbeta->SetRightMargin(0.0507614);
	C_Evsbeta->SetBottomMargin(0.125654);
	C_Evsbeta->SetTopMargin(0.0746073);
	*/
	char tstring_gr_threh[100];
	sprintf(tstring_gr_threh,"Rigidity Threshold for _{%d}%s",plot_Z,chem_s);
	
	
	
	TLegend *leg_northern_shouthern_hemisphere_cut = new TLegend(0.685619,0.0129534,0.955686,0.0518135);
	leg_northern_shouthern_hemisphere_cut->AddEntry(gr_northern_hemisphere_cut,tstring_gr_threh,"l");	
	leg_northern_shouthern_hemisphere_cut->SetFillColor(10);
	leg_northern_shouthern_hemisphere_cut->SetTextFont(nicefont);
	leg_northern_shouthern_hemisphere_cut->SetBorderSize(0); // remove ugly border from stats 
	leg_northern_shouthern_hemisphere_cut->SetTextSize(0.03);

	
	
	leg_northern_shouthern_hemisphere_cut->Draw();
    

    
    // East-West Dependent Analysis
    
    // Read in geomagnetic latitude residence times
	char geo_lat_file[200] = "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/ISS_Position/geomagnetic_latitude_fractions_original.txt";
	int g_m_l_count=0;
	double angle_g_m_l[90];
	double orbit_fraction_g_m_l[90];
	
	char tstring_h_g_m_l[100];
	sprintf(tstring_h_g_m_l,"ISS Geomagnetic Latitude Residence");
	
    //	TH1F *h_g_m_l = new TH1F("h_g_m_l_p",tstring_h_g_m_l_p,100,0,20);
	TH1F *h_g_m_l= new TH1F("h_g_m_l",tstring_h_g_m_l,90,0,90);
	
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
	h_g_m_l->Draw("");

     
    
    // Interpolated Vertical Cutoff Plot
    
    
	TH2F *h_geo_lat_bl = new TH2F("h_geo_lat_bl","",long_bins_1,-180,180,lat_bins_1,-90,90);
    //	h_geo_lat_bl->SetTitle("Geomagnetic Vertical Cutoff Rigidities");
	h_geo_lat_bl->SetTitle("");
	h_geo_lat_bl->GetXaxis()->SetTitle("Longitude (degrees)");
	h_geo_lat_bl->GetXaxis()->CenterTitle();
	h_geo_lat_bl->GetXaxis()->SetTitleOffset(1.2);
	h_geo_lat_bl->GetYaxis()->SetTitle("Latitude (degrees)");
	h_geo_lat_bl->GetYaxis()->CenterTitle();
	h_geo_lat_bl->GetYaxis()->SetTitleOffset(1.2);
	h_geo_lat_bl->GetYaxis()->SetLabelFont(nicefont);
	h_geo_lat_bl->GetYaxis()->SetTitleFont(nicefont);
	h_geo_lat_bl->GetXaxis()->SetLabelFont(nicefont);
	h_geo_lat_bl->GetXaxis()->SetTitleFont(nicefont);
	h_geo_lat_bl->GetZaxis()->SetTitle("Vertical Cutoff Rigidity (GV)");
	h_geo_lat_bl->GetZaxis()->SetTitleFont(nicefont);
	h_geo_lat_bl->GetZaxis()->SetTitleOffset(0.7);
	h_geo_lat_bl->GetZaxis()->CenterTitle();
    //h_geo_lat_bl->GetZaxis()->SetTickSize(0.5);
	
	h_geo_lat_bl->GetZaxis()->SetRangeUser(-0.1,16);
    /*
    double cut_off[long_bins_1][lat_bins_1];
	double cut_off_bl[long_bins_1][lat_bins_1];
	double cut_off_eval[long_bins_1][lat_bins_1];
    int index_long;
	int index_lat;
    
    
    for (index_long=0; index_long<long_bins_1; index_long++) {
		for (index_lat=0; index_lat<lat_bins_1; index_lat++) {
			cut_off[index_long][index_lat]=0.0;
			cut_off_bl[index_long][index_lat]=0.0;
			cut_off_eval[index_long][index_lat]=0.0;
		}
	}
	
	
	for (long_index=0; long_index<long_bins_calc; long_index++) {
		for (lat_index=0; lat_index<lat_bins; lat_index++) {
			cut_off[15*(long_index)][5*(lat_bins-lat_index-1)]=alt_450[lat_bins-lat_index-1].cut_off[long_index];
            //fprintf(stderr,"long: %3.0f lat: %3.0f cut_off: %.4f\n",longitude[long_index],alt_450[lat_bins-lat_index-1].latitude,cut_off[15*(long_index)][5*(lat_bins-lat_index-1)]);
		}
	}
    
    int lat_index_lo;
	int lat_index_hi;
	int long_index_lo;
	int long_index_hi;
	double lat_weight;
	double long_weight;
	double lat_plot_1[lat_bins_1];
	double long_plot_1[long_bins_1];
	lat_plot_1[lat_bins_1]=90.0;
	
	
	for (index_long=0; index_long<long_bins_1; index_long++) {
		long_index_lo=(int)(index_long/15);
		long_index_hi=(int)(index_long/15 + 1);

		long_weight=(index_long-longitude[long_index_lo])/(longitude[long_index_hi]-longitude[long_index_lo]);
		for (index_lat=0; index_lat<lat_bins_1; index_lat++) {
			lat_index_lo=lat_bins-index_lat/5 -1;
			lat_index_hi=lat_bins-index_lat/5 -2;
			if (lat_index_hi == -1 ) {
				lat_index_hi = 0;
			}
			if (lat_index_hi != lat_index_lo) {
				lat_weight=((index_lat-90)-alt_450[lat_index_lo].latitude)/(alt_450[lat_index_hi].latitude-alt_450[lat_index_lo].latitude);
			}
			else {
				lat_weight=0;
			}
            
			
			cut_off[index_long][index_lat]=0.5*
            (
             (0.5*(alt_450[lat_index_lo].cut_off[long_index_lo]+alt_450[lat_index_lo].cut_off[long_index_hi]) -
              lat_weight*0.5*((alt_450[lat_index_lo].cut_off[long_index_lo]+alt_450[lat_index_lo].cut_off[long_index_hi]) -
                              (alt_450[lat_index_hi].cut_off[long_index_lo]+alt_450[lat_index_hi].cut_off[long_index_hi]))) +
             (0.5*(alt_450[lat_index_lo].cut_off[long_index_lo]+alt_450[lat_index_hi].cut_off[long_index_lo]) -
              long_weight*0.5*((alt_450[lat_index_lo].cut_off[long_index_lo]+alt_450[lat_index_hi].cut_off[long_index_lo]) -
                               (alt_450[lat_index_lo].cut_off[long_index_hi]+alt_450[lat_index_hi].cut_off[long_index_hi])))
             );
			
			
			cut_off_bl[index_long][index_lat] = (1.0/((longitude[long_index_hi]-longitude[long_index_lo])*(alt_450[lat_index_hi].latitude-alt_450[lat_index_lo].latitude)))*
            (
             alt_450[lat_index_lo].cut_off[long_index_lo]*(longitude[long_index_hi]-index_long)*(alt_450[lat_index_hi].latitude-(index_lat-90)) +
             alt_450[lat_index_lo].cut_off[long_index_hi]*(index_long-longitude[long_index_lo])*(alt_450[lat_index_hi].latitude-(index_lat-90)) +
             alt_450[lat_index_hi].cut_off[long_index_lo]*(longitude[long_index_hi]-index_long)*((index_lat-90)-alt_450[lat_index_lo].latitude) +
             alt_450[lat_index_hi].cut_off[long_index_hi]*(index_long-longitude[long_index_lo])*((index_lat-90)-alt_450[lat_index_lo].latitude)
             );
			
			
			lat_plot_1[index_lat]=index_lat-90;
			long_plot_1[index_long]=index_long;
			if (long_plot_1[index_long]>= 180) {
				long_plot_1[index_long]=index_long-360;
			}
			h_geo_lat_bl->Fill(long_plot_1[index_long],lat_plot_1[index_lat],cut_off_bl[index_long][index_lat]);
			
			fprintf(stderr,"index_long: %3d long_index_lo: %2d long_index_hi: %2d index_lat: %3d lat_index_lo: %2d lat_index_hi: %2d long_weight: %.4f lat_weight: %.4f cut_off: %.4f\n",
					index_long,long_index_lo,long_index_hi,index_lat,lat_index_lo,lat_index_hi,long_weight,lat_weight,cut_off_bl[index_long][index_lat]); // output the biliner interpolation of geomagnetic vertical cut_off map
		}
	}

   
    
    

    
    
    TCanvas *c_i_bl = new TCanvas("c_i_bl","bilinear interpolated cut_off",1200,800);
	h_geo_lat_bl->GetZaxis()->SetRangeUser(-0.1,16);
	h_geo_lat_bl->Draw("COLZ");
	gr_mercator->Draw("P");
	gr_ISS_orbit->Draw("P");
    
    gr_northern_hemisphere_cut->SetLineStyle(2);
	gr_northern_hemisphere_cut->SetLineWidth(4);
	gr_northern_hemisphere_cut->Draw("C");
	gr_southern_hemisphere_cut->SetLineStyle(2);
	gr_southern_hemisphere_cut->SetLineWidth(4);
	gr_southern_hemisphere_cut->Draw("C");
    
	leg_northern_shouthern_hemisphere_cut->Draw();
    */
    
    
    // Write integral spectra files
    char S_max_out[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/S_max_integral_spectra_2026-06-06.txt";
    
    char S_min_out[200] =
    "/Users/william/Desktop/TIGERISS_Sim/TIGERISS_Abundances/tables/S_min_integral_spectra_2026-06-06.txt";
    
    int S_Fe_index=0;
    
    // S_max file
    if((outfile = fopen(S_max_out,"w")) == NULL)
    {    fprintf(stderr,"Error opening %s for output\n",S_max_out);
        exit(0);    }
    
    
    fprintf(outfile," Z  A X  ");
    for (i=0; i<E_index; i++) {
        fprintf(outfile,"%.1f ",S_max[0].EI[i]);
    }
    fprintf(outfile,"\n");
    for (i=0; i<92; i++) {
        fprintf(outfile,"%d %f %s ",el_dat[i].Z,el_dat[i].A,el_dat[i].chem_s);
        if (el_dat[i].Z==26) {
            S_Fe_index=i+1;
        }
        for (j=0; j<E_index; j++) {
            if (el_dat[i].Z > 30) {
                if (j==E_index-1) {
                    fprintf(outfile,"%e ",S_max[S_Fe_index].I_E_hi*el_dat[i].TOA_abund/el_dat[el_dat_Fe_index].TOA_abund);
                }
                else {
                    fprintf(outfile,"%e ",S_max[S_Fe_index].EI_int[j]*el_dat[i].TOA_abund/el_dat[el_dat_Fe_index].TOA_abund);
                }
            }
            else {
                if (j==E_index-1) {
                    fprintf(outfile,"%e ",S_max[i+1].I_E_hi);
                }
                else {
                    fprintf(outfile,"%e ",S_max[i+1].EI_int[j]);
                }
            }
        }
        fprintf(outfile,"\n");
    }
    
    fclose(outfile);
    
    // S_min file
    if((outfile = fopen(S_min_out,"w")) == NULL)
    {    fprintf(stderr,"Error opening %s for output\n",S_min_out);
        exit(0);    }
    
    
    fprintf(outfile," Z  A X  ");
    for (i=0; i<E_index; i++) {
        fprintf(outfile,"%.1f ",S_min[0].EI[i]);
    }
    fprintf(outfile,"\n");
    for (i=0; i<92; i++) {
        fprintf(outfile,"%d %f %s ",el_dat[i].Z,el_dat[i].A,el_dat[i].chem_s);
        if (el_dat[i].Z==26) {
            S_Fe_index=i+1;
        }
        for (j=0; j<E_index; j++) {
            if (el_dat[i].Z > 30) {
                if (j==E_index-1) {
                    fprintf(outfile,"%e ",S_min[S_Fe_index].I_E_hi*el_dat[i].TOA_abund/el_dat[el_dat_Fe_index].TOA_abund);
                }
                else {
                    fprintf(outfile,"%e ",S_min[S_Fe_index].EI_int[j]*el_dat[i].TOA_abund/el_dat[el_dat_Fe_index].TOA_abund);
                }
            }
            else {
                if (j==E_index-1) {
                    fprintf(outfile,"%e ",S_min[i+1].I_E_hi);
                }
                else {
                    fprintf(outfile,"%e ",S_min[i+1].EI_int[j]);
                }
            }
        }
        fprintf(outfile,"\n");
    }
    
    fclose(outfile);
    
    
    fprintf(stderr,"E_index: %d\n",E_index);
    
	
	theApp.Run();
  
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



