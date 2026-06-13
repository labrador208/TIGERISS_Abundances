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
//#include <TPaletteAxis.h>
#include <TGraphAsymmErrors.h>

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
	
	
	
	
	
//	TApplication theApp("App", &argc, argv);
	
	ifstream in_1;
	ofstream out_1;
	
	TH2F *h_dgf_b1 = new TH2F("h_dgf_b1","",90,0,90,360,0,360);
	TH2F *h_dgf_g_b1 = new TH2F("h_dgf_g_b1","",180,0,180,180,0,180);
	TH2F *h_dgf_g_theta_b1 = new TH2F("h_dgf_g_theta_b1","",36,0,180,18,0,90);

    
    char extension[20];
    
	if(argc!=2) {
		fprintf(stderr,"usage: %s file extension\n",argv[0]);
		return 1;
	}
	
	sprintf(extension,"%s",argv[1]);
    
	// Detect if extension contains "_rotated" - if so, remove it for file lookup
	char file_ext[50];
	sprintf(file_ext, "%s", extension);
	char* rotated_pos = strstr(file_ext, "_rotated");
	bool is_rotated = (rotated_pos != NULL);
	if (rotated_pos != NULL) {
		*rotated_pos = '\0';  // Truncate at "_rotated"
	}
	
	// Extract angle number for title (e.g., "_deg000" -> "000")
	char angle_str[20]="";
	const char* deg_pos = strstr(file_ext, "deg");
	if (deg_pos != NULL) {
		sprintf(angle_str, "%s", deg_pos + 3);  // Skip "deg"
	} else {
		sprintf(angle_str, "N/A");
	}
	
	char dir[200]="./";
    
	char filename[200]="";
	char title_str[200]="";
	
	// Differential Geometry Factor interlude:

//char geom_file[200]="/data/calet/abundance_estimate/rigidity/code/geofac_dir/dgf_2d.txt";

    char geom_file[200];
	sprintf(geom_file,"%sdgf_2d%s.txt",dir,extension);
    
	int dgf_index=0;
	int dgf_count=0;

	double theta_dgf[1300];
	double phi_dgf[1300];
	double geo_factor_dgf[1300];
	double geo_factor_dgf_plot[1300];

    double geo_fac_tot=0.0;

  	double angle_dgf_1D[100];
  	double angle_dgf_b1_1D[100];
	double geo_factor_dgf_1D[100];  
	double geo_factor_dgf_b1_1D[100];  
	int dgf_index_b1=0;
	int dgf_count_b1=0;
	int dgf_count_1D=0;
	int dgf_count_b1_1D=0;
	angle_dgf_b1_1D[0]=0.0;
	geo_factor_dgf_b1_1D[0]=0.0;


	in_1.open(geom_file);
	while(1) {
		in_1 >> theta_dgf[dgf_count] >> phi_dgf[dgf_count] >> geo_factor_dgf[dgf_count];
		if(in_1.eof() || !in_1.good()) break;
		
		geo_factor_dgf_plot[dgf_count] = geo_factor_dgf[dgf_count]/25.0;

		geo_fac_tot+=geo_factor_dgf[dgf_count];
		if(phi_dgf[dgf_count] > 2.5) {
		geo_factor_dgf_1D[dgf_count_1D]+=geo_factor_dgf[dgf_count]/5.0;
		}     
		if(phi_dgf[dgf_count] == 2.5) {
		dgf_count_1D++;
		geo_factor_dgf_1D[dgf_count_1D]=geo_factor_dgf[dgf_count]/5.0;
		angle_dgf_1D[dgf_count_1D]=theta_dgf[dgf_count];
		}
		dgf_count++;
	}
	in_1.close();
	in_1.clear();

    fprintf(stderr,"dgf_count: %d\n",dgf_count);
    fprintf(stderr,"%s geo_fac_tot= %.4f\n",geom_file,geo_fac_tot);
	
	
		// Contour Map	
	
	TGraph2D *g_dgf = 
    new TGraph2D(dgf_count,theta_dgf,phi_dgf,geo_factor_dgf_plot);
	g_dgf->SetTitle("");

		//	g_dgf->SetTitle("Geomagnetic Vertical Cutoff Rigidities");

	//g_dgf->GetZaxis()->CenterTitle();
		//g_dgf->GetZaxis()->SetTickSize(0.5);
	
		//g_dgf->GetZaxis()->SetRangeUser(-0.1,16);
	



	TCanvas *C_dgf = new TCanvas("C_dgf","5 degree gamma-theta DGF",0,0,800,800);
	C_dgf->SetLeftMargin(0.139594);
	C_dgf->SetRightMargin(0.0596447);
	C_dgf->SetBottomMargin(0.11911);
	C_dgf->SetTopMargin(0.15);
	gStyle->SetPalette(1);

	
	g_dgf->Draw("SURF1");

	gPad->Update();
	sprintf(title_str,"Angle %s ",angle_str);
	TPaveText *pt_5d = new TPaveText(0.2, 0.92, 0.8, 1.0, "NDC");
	pt_5d->AddText(title_str);
	pt_5d->SetFillStyle(0);
	pt_5d->SetBorderSize(0);
	pt_5d->SetTextSize(0.03);
	pt_5d->SetTextFont(42);
	pt_5d->SetTextAlign(22);
	pt_5d->Draw();

	g_dgf->GetXaxis()->SetTitle("#theta (degrees)");
	g_dgf->GetXaxis()->CenterTitle();
	g_dgf->GetXaxis()->SetTitleOffset(1.7);
	g_dgf->GetYaxis()->SetTitle("#phi (degrees)");
	g_dgf->GetYaxis()->CenterTitle();
	g_dgf->GetYaxis()->SetTitleOffset(1.7);
	g_dgf->GetYaxis()->SetLabelFont(nicefont);
	g_dgf->GetYaxis()->SetTitleFont(nicefont);
	g_dgf->GetXaxis()->SetLabelFont(nicefont);
	g_dgf->GetXaxis()->SetTitleFont(nicefont);
	g_dgf->GetZaxis()->SetTitle("DGF (cm^{2}sr/(1^{o} bins))");
	g_dgf->GetZaxis()->SetTitleFont(nicefont);
	g_dgf->GetZaxis()->SetTitleOffset(1.30);


	C_dgf->SetGrid();

	// writing plot files
	sprintf(filename,"plots/dgf_2D/dgf_theta_phi_b5_2D%s.png",extension);
	C_dgf->Print(filename);
	sprintf(filename,"plots/dgf_2D/dgf_theta_phi_b5_2D%s.eps",extension);
	C_dgf->Print(filename);
	
	
		// Differential Geometry Factor interlude:
	
	char geom_file_b1[200]="/data/calet/abundance_estimate/rigidity/code/geofac_dir/dgf_2d_bin1.txt";
	
    

	
	double theta_dgf_b1;
	double phi_dgf_b1;
	double geo_factor_dgf_b1;
	
    sprintf(geom_file,"%sdgf_2d_bin1%s.txt",dir,extension);

    geo_fac_tot=0.0;
    
	in_1.open(geom_file);
	while(1) {
		in_1 >> theta_dgf_b1 >> phi_dgf_b1 >> geo_factor_dgf_b1;
		if(in_1.eof() || !in_1.good()) break;
		h_dgf_b1->Fill(theta_dgf_b1,phi_dgf_b1,geo_factor_dgf_b1);
        geo_fac_tot+=geo_factor_dgf_b1;
    if(phi_dgf_b1 > 0.5) {
      geo_factor_dgf_b1_1D[dgf_count_b1_1D]+=geo_factor_dgf_b1;
    }     
    if(phi_dgf_b1 == 0.5) {
      dgf_count_b1_1D++;
      geo_factor_dgf_b1_1D[dgf_count_b1_1D]=geo_factor_dgf_b1;
      angle_dgf_b1_1D[dgf_count_b1_1D]=theta_dgf_b1;
    }
		
		dgf_count_b1++;
	}
	in_1.close();
	in_1.clear();
	
    fprintf(stderr,"dgf_count_b1: %d\n",dgf_count_b1);
    fprintf(stderr,"%s geo_fac_tot= %.4f\n",geom_file,geo_fac_tot);
	
	
		// Contour Map	

	
		//	h_dgf_b1->SetTitle("Geomagnetic Vertical Cutoff Rigidities");
	h_dgf_b1->GetXaxis()->SetTitle("#theta (TESTEST)");
	h_dgf_b1->GetXaxis()->CenterTitle();
	h_dgf_b1->GetXaxis()->SetTitleOffset(1.7);
	h_dgf_b1->GetYaxis()->SetTitle("#phi (degrees)");
	h_dgf_b1->GetYaxis()->CenterTitle();
	h_dgf_b1->GetYaxis()->SetTitleOffset(1.7);
	h_dgf_b1->GetYaxis()->SetLabelFont(nicefont);
	h_dgf_b1->GetYaxis()->SetTitleFont(nicefont);
	h_dgf_b1->GetXaxis()->SetLabelFont(nicefont);
	h_dgf_b1->GetXaxis()->SetTitleFont(nicefont);
	h_dgf_b1->GetZaxis()->SetTitle("DGF (cm^{2}sr/(1^{o} bins))");
	h_dgf_b1->GetZaxis()->SetTitleFont(nicefont);
	h_dgf_b1->GetZaxis()->SetTitleOffset(1.7);
	h_dgf_b1->GetZaxis()->CenterTitle();
		//h_dgf_b1->GetZaxis()->SetTickSize(0.5);
	
		//h_dgf_b1->GetZaxis()->SetRangeUser(-0.1,16);
	
	
	
	
	TCanvas *C_dgf_b1 = new TCanvas("C_dgf_b1","1 degree theta-phi DGF",0,0,800,800);
	C_dgf_b1->SetLeftMargin(0.139594);
	C_dgf_b1->SetRightMargin(0.0596447);
	C_dgf_b1->SetBottomMargin(0.11911);
	C_dgf_b1->SetTopMargin(0.15);
	
	h_dgf_b1->Draw("SURF1");
	gPad->Update();
	sprintf(title_str,"Angle offset %s ",angle_str);
	TPaveText *pt_b1_2d = new TPaveText(0.2, 0.92, 0.8, 1.0, "NDC");
	pt_b1_2d->AddText(title_str);
	pt_b1_2d->SetFillStyle(0);
	pt_b1_2d->SetBorderSize(0);
	pt_b1_2d->SetTextSize(0.03);
	pt_b1_2d->SetTextFont(42);
	pt_b1_2d->SetTextAlign(22);
	pt_b1_2d->Draw();
	
	C_dgf_b1->SetGrid();
	
	// writing plot files
	sprintf(filename,"plots/dgf_2D/dgf_theta_phi_b1_2D%s.png",extension);
	C_dgf_b1->Print(filename);
	sprintf(filename,"plots/dgf_2D/dgf_theta_phi_b1_2D%s.eps",extension);
	C_dgf_b1->Print(filename);
	

	
		// Differential Geometry Factor interlude:
	
	char geom_file_g[200]="/data/calet/abundance_estimate/rigidity/code/geofac_dir/dgf_2d_g.txt";
	
	int g_dgf_index=0;
	int g_dgf_count=0;
	
	double gamma_g_dgf[1300];
	double delta_g_dgf[1300];
	double geo_factor_g_dgf[1300];
	double geo_factor_g_dgf_plot[1300];

    sprintf(geom_file,"%sdgf_2d_g%s.txt",dir,extension);

    geo_fac_tot=0.0;
	
	in_1.open(geom_file);
	while(1) {
		in_1 >> gamma_g_dgf[g_dgf_count] >> delta_g_dgf[g_dgf_count] >> geo_factor_g_dgf[g_dgf_count];
		if(in_1.eof() || !in_1.good()) break;
			//fprintf(stderr,"%f %f %f\n",gamma_g_dgf[g_dgf_count],delta_g_dgf[g_dgf_count],geo_factor_g_dgf[g_dgf_count]);
        geo_fac_tot+=geo_factor_g_dgf[g_dgf_count];
        geo_factor_g_dgf_plot[g_dgf_count] = geo_factor_g_dgf[g_dgf_count]/25;
		g_dgf_count++;
	}
	in_1.close();
	in_1.clear();
    
    fprintf(stderr,"g_dgf_count: %d\n",g_dgf_count);
    fprintf(stderr,"%s geo_fac_tot= %.4f\n",geom_file,geo_fac_tot);
	

		// Contour Map	
	
	TGraph2D *g_g_dgf = 
    new TGraph2D(dgf_count,gamma_g_dgf,delta_g_dgf,geo_factor_g_dgf_plot);
	g_g_dgf->SetTitle("");

		//	g_g_dgf->SetTitle("Geomagnetic Vertical Cutoff Rigidities");

		//g_g_dgf->GetZaxis()->SetTickSize(0.5);
	
		//g_g_dgf->GetZaxis()->SetRangeUser(-0.1,16);
	
	
	TCanvas *C_g_dgf = new TCanvas("C_g_dgf","5 degree gamma-delta DGF",0,0,800,800);
	C_g_dgf->SetLeftMargin(0.139594);
	C_g_dgf->SetRightMargin(0.0596447);
	C_g_dgf->SetBottomMargin(0.11911);
	C_g_dgf->SetTopMargin(0.15);

	
	g_g_dgf->Draw("SURF1");

	gPad->Update();
	sprintf(title_str,"Angle offset %s ",angle_str);
	TPaveText *pt_g2d = new TPaveText(0.2, 0.92, 0.8, 1.0, "NDC");
	pt_g2d->AddText(title_str);
	pt_g2d->SetFillStyle(0);
	pt_g2d->SetBorderSize(0);
	pt_g2d->SetTextSize(0.03);
	pt_g2d->SetTextFont(42);
	pt_g2d->SetTextAlign(22);
	pt_g2d->Draw();


	g_g_dgf->GetXaxis()->SetTitle("#gamma (Testtestestestest)");
	g_g_dgf->GetXaxis()->CenterTitle();
	g_g_dgf->GetXaxis()->SetTitleOffset(1.7);
	g_g_dgf->GetYaxis()->SetTitle("#delta (degrees)");
	g_g_dgf->GetYaxis()->CenterTitle();
	g_g_dgf->GetYaxis()->SetTitleOffset(1.7);
	g_g_dgf->GetYaxis()->SetLabelFont(nicefont);
	g_g_dgf->GetYaxis()->SetTitleFont(nicefont);
	g_g_dgf->GetXaxis()->SetLabelFont(nicefont);
	g_g_dgf->GetXaxis()->SetTitleFont(nicefont);
	g_g_dgf->GetZaxis()->SetTitle("DGF (cm^{2}sr/(1^{o} bins))");
	g_g_dgf->GetZaxis()->SetTitleFont(nicefont);
	g_g_dgf->GetZaxis()->SetTitleOffset(1.7);
	g_g_dgf->GetZaxis()->CenterTitle();

	
	C_g_dgf->SetGrid();

    // writing plot files
    sprintf(filename,"plots/dgf_2D/dgf_gamma_theta_b5_2D%s.png",extension);
    C_g_dgf->Print(filename);
    sprintf(filename,"plots/dgf_2D/dgf_gamma_theta_b5_2D%s.eps",extension);
    C_g_dgf->Print(filename);


		// Differential Geometry Factor interlude:
	
	char geom_file_g_b1[200]="/data/calet/abundance_estimate/rigidity/code/geofac_dir/dgf_2d_g_bin1.txt";
	
	int g_dgf_count_b1=0;
	
	double gamma_g_dgf_b1;
	double delta_g_dgf_b1;
	double geo_factor_g_dgf_b1;
	
    sprintf(geom_file,"%sdgf_2d_g_bin1%s.txt",dir,extension);

    geo_fac_tot=0.0;
    
	in_1.open(geom_file);
	while(1) {
		in_1 >> gamma_g_dgf_b1 >> delta_g_dgf_b1 >> geo_factor_g_dgf_b1;
		if(in_1.eof() || !in_1.good()) break;
		h_dgf_g_b1->Fill(gamma_g_dgf_b1,delta_g_dgf_b1,geo_factor_g_dgf_b1);
        geo_fac_tot+=geo_factor_g_dgf_b1;
		g_dgf_count_b1++;
	}
	in_1.close();
	in_1.clear();
    
    fprintf(stderr,"g_dgf_count_b1: %d\n",g_dgf_count_b1);
    fprintf(stderr,"%s geo_fac_tot= %.4f\n",geom_file,geo_fac_tot);
	
	
		// Contour Map	
	

		//	h_dgf_g_b1->SetTitle("Geomagnetic Vertical Cutoff Rigidities");
	h_dgf_g_b1->GetXaxis()->SetTitle("#gamma (degrees)");
	h_dgf_g_b1->GetXaxis()->CenterTitle();
	h_dgf_g_b1->GetXaxis()->SetTitleOffset(1.7);
	h_dgf_g_b1->GetYaxis()->SetTitle("#delta (degrees)");
	h_dgf_g_b1->GetYaxis()->CenterTitle();
	h_dgf_g_b1->GetYaxis()->SetTitleOffset(1.7);
	h_dgf_g_b1->GetYaxis()->SetLabelFont(nicefont);
	h_dgf_g_b1->GetYaxis()->SetTitleFont(nicefont);
	h_dgf_g_b1->GetXaxis()->SetLabelFont(nicefont);
	h_dgf_g_b1->GetXaxis()->SetTitleFont(nicefont);
	h_dgf_g_b1->GetZaxis()->SetTitle("DGF (cm^{2}sr/(1^{o} bins))");
	h_dgf_g_b1->GetZaxis()->SetTitleFont(nicefont);
	h_dgf_g_b1->GetZaxis()->SetTitleOffset(1.7);
	h_dgf_g_b1->GetZaxis()->CenterTitle();
		//h_dgf_g_b1->GetZaxis()->SetTickSize(0.5);
	
		//h_dgf_g_b1->GetZaxis()->SetRangeUser(-0.1,16);
	
	
	TCanvas *C_g_dgf_b1 = new TCanvas("C_g_dgf_b1","1 degree gamma-delta DGF",0,0,800,800);
	C_g_dgf_b1->SetLeftMargin(0.139594);
	C_g_dgf_b1->SetRightMargin(0.0596447);
	C_g_dgf_b1->SetBottomMargin(0.11911);
	C_g_dgf_b1->SetTopMargin(0.15);
	
			h_dgf_g_b1->Draw("SURF1");
			gPad->Update();
	sprintf(title_str,"Angle offset %s ",angle_str);
	TPaveText *pt_gb1 = new TPaveText(0.2, 0.92, 0.8, 1.0, "NDC");
	pt_gb1->AddText(title_str);
	pt_gb1->SetFillStyle(0);
	pt_gb1->SetBorderSize(0);
	pt_gb1->SetTextSize(0.03);
	pt_gb1->SetTextFont(42);
	pt_gb1->SetTextAlign(22);
	pt_gb1->Draw();
	
	C_g_dgf_b1->SetGrid();

    // writing plot files
    sprintf(filename,"plots/dgf_2D/dgf_gamma_delta_b1_2D%s.png",extension);
    C_g_dgf_b1->Print(filename);
    sprintf(filename,"plots/dgf_2D/dgf_gamma_delta_b1_2D%s.eps",extension);
    C_g_dgf_b1->Print(filename);

	
    // Differential Geometry Factor for gamma versus theta interlude:

    // Contour Map
	
    
    //	h_dgf_g_theta_b1->SetTitle("Geomagnetic Vertical Cutoff Rigidities");
	h_dgf_g_theta_b1->GetXaxis()->SetTitle("#gamma (degrees)");
	h_dgf_g_theta_b1->GetXaxis()->CenterTitle();
	h_dgf_g_theta_b1->GetXaxis()->SetTitleOffset(1.7);
	h_dgf_g_theta_b1->GetYaxis()->SetTitle("#theta (degrees)");
	h_dgf_g_theta_b1->GetYaxis()->CenterTitle();
	h_dgf_g_theta_b1->GetYaxis()->SetTitleOffset(1.7);
	h_dgf_g_theta_b1->GetYaxis()->SetLabelFont(nicefont);
	h_dgf_g_theta_b1->GetYaxis()->SetTitleFont(nicefont);
	h_dgf_g_theta_b1->GetXaxis()->SetLabelFont(nicefont);
	h_dgf_g_theta_b1->GetXaxis()->SetTitleFont(nicefont);
	h_dgf_g_theta_b1->GetZaxis()->SetTitle("DGF (cm^{2}sr/(1^{o} bins))");
	h_dgf_g_theta_b1->GetZaxis()->SetTitleFont(nicefont);
	h_dgf_g_theta_b1->GetZaxis()->SetTitleOffset(1.7);
	h_dgf_g_theta_b1->GetZaxis()->CenterTitle();
    //h_dgf_g_theta_b1->GetZaxis()->SetTickSize(0.5);
	
    //h_dgf_g_theta_b1->GetZaxis()->SetRangeUser(-0.1,16);


    
	char geom_file_g_theta_b1[200]="/data/calet/abundance_estimate/rigidity/code/geofac_dir/dgf_2d_g_theta_bin1.txt";
	
	int g_dgf_count_theta_b1=0;
	
	double gamma_g_dgf_theta_b1;
	double theta_g_dgf_theta_b1;
	double geo_factor_g_dgf_theta_b1;
	
    geo_fac_tot=0.0;
    
    sprintf(geom_file,"%sdgf_2d_g_theta_bin1%s.txt",dir,extension);
    
	
	in_1.open(geom_file);
	while(1) {
		in_1 >> gamma_g_dgf_theta_b1 >> theta_g_dgf_theta_b1 >> geo_factor_g_dgf_theta_b1;
		if(in_1.eof() || !in_1.good()) break;
		h_dgf_g_theta_b1->Fill(gamma_g_dgf_theta_b1,theta_g_dgf_theta_b1,geo_factor_g_dgf_theta_b1/25);
		g_dgf_count_theta_b1++;
        geo_fac_tot+=geo_factor_g_dgf_theta_b1;
	}
	in_1.close();
	in_1.clear();
	
    fprintf(stderr,"g_dgf_count_theta_b1: %d\n",g_dgf_count_theta_b1);
    fprintf(stderr,"%s geo_fac_tot= %.4f\n",geom_file,geo_fac_tot);

	
	
	TCanvas *C_g_dgf_theta_b1 = new TCanvas("C_g_dgf_theta_b1","1 degree gamma-theta DGF",0,0,800,800);
	C_g_dgf_theta_b1->SetLeftMargin(0.139594);
	C_g_dgf_theta_b1->SetRightMargin(0.0596447);
	C_g_dgf_theta_b1->SetBottomMargin(0.11911);
	C_g_dgf_theta_b1->SetTopMargin(0.0811518);
	C_g_dgf_theta_b1->SetTheta(30.);
	C_g_dgf_theta_b1->SetPhi(125.);
	
    h_dgf_g_theta_b1->Draw("SURF1");
    gPad->Update();
    sprintf(title_str,"Angle offset %s ",angle_str);
    TPaveText *pt_gtheta = new TPaveText(0.2, 0.92, 0.8, 1.0, "NDC");
    pt_gtheta->AddText(title_str);
    pt_gtheta->SetFillStyle(0);
    pt_gtheta->SetBorderSize(0);
    pt_gtheta->SetTextSize(0.03);
    pt_gtheta->SetTextFont(42);
    pt_gtheta->SetTextAlign(22);
    pt_gtheta->Draw();

	
	
	C_g_dgf_theta_b1->SetGrid();

    // writing plot files
    sprintf(filename,"plots/dgf_2D/dgf_gamma_theta_b1_2D%s.png",extension);
    C_g_dgf_theta_b1->Print(filename);
    sprintf(filename,"plots/dgf_2D/dgf_gamma_theta_b1_2D%s.eps",extension);
    C_g_dgf_theta_b1->Print(filename);


  // Do 1-D differential geometry factors

		// 1 degree bins
		// Differential Geometry Factor interlude:
	std::cout<< extension<<std::endl;
    sprintf(geom_file,"%sdgf_1d_bin1%s.txt",dir,extension);
	


    geo_fac_tot=0;
	

	
/*
	
	in_1.open(geom_file);
	while(1) {
		in_1 >> angle_dgf_b1[dgf_count_b1] >> geo_factor_dgf_b1_1D[dgf_count_b1];
		if(in_1.eof() || !in_1.good()) break;
		
		geo_fac_tot+=geo_factor_dgf_b1_1D[dgf_count_b1];
		dgf_count_b1++;
	}
	in_1.close();
	in_1.clear();

    fprintf(stderr,"g_dgf_count_b1: %d\n",g_dgf_count_b1);
    fprintf(stderr,"%s geo_fac_tot= %.4f\n",geom_file,geo_fac_tot);
*/	
	
	TGraph *gr_dgf_1D = new TGraph(dgf_count_1D,angle_dgf_1D,geo_factor_dgf_1D);
	gr_dgf_1D->SetTitle("");
	
	gr_dgf_1D->SetMarkerStyle(1);
	gr_dgf_1D->SetMarkerSize(1);
		//	gr_dgf_1D->SetTitle("CALET Differential Geometry Factor");
	gr_dgf_1D->GetXaxis()->SetTitle("Incidence Angle (degrees)");
	gr_dgf_1D->GetXaxis()->CenterTitle();
	gr_dgf_1D->GetYaxis()->SetTitle("Differential Geometry Factor (cm^{2}sr/(1^{o} bins))");
	gr_dgf_1D->GetYaxis()->CenterTitle();
	gr_dgf_1D->GetXaxis()->SetTitleOffset(1.2);
	gr_dgf_1D->GetYaxis()->SetTitleOffset(1.3);
	gr_dgf_1D->GetYaxis()->SetLabelFont(nicefont);
	gr_dgf_1D->GetYaxis()->SetTitleFont(nicefont);
	gr_dgf_1D->GetXaxis()->SetLabelFont(nicefont);
	gr_dgf_1D->GetXaxis()->SetTitleFont(nicefont);
	gr_dgf_1D->GetXaxis()->SetRangeUser(0,90);
	//gr_dgf_1D->GetYaxis()->SetRangeUser(0,40);
	
	
	
	TCanvas *C_dgf_1D = new TCanvas("C_dgf_1D","5 degree 1D DGF",0,0,800,800);
	C_dgf_1D->SetLeftMargin(0.139594);
	C_dgf_1D->SetRightMargin(0.0596447);
	C_dgf_1D->SetBottomMargin(0.11911);
	C_dgf_1D->SetTopMargin(0.15);
	
	gr_dgf_1D->SetLineWidth(2);
	gr_dgf_1D->Draw("AL");
	gPad->Update();
	
	// Add title as text
	sprintf(title_str,"Angle offset %s ",angle_str);
	TPaveText *pt = new TPaveText(0.2, 0.92, 0.8, 1.0, "NDC");
	pt->AddText(title_str);
	pt->SetFillStyle(0);
	pt->SetBorderSize(0);
	pt->SetTextSize(0.03);
	pt->SetTextFont(42);
	pt->SetTextAlign(22);
	pt->Draw();
	
	C_dgf_1D->SetGrid();

    // writing plot files
    sprintf(filename,"plots/dgf_2D/dgf_b5_1D%s.png",extension);
    C_dgf_1D->Print(filename);
    sprintf(filename,"plots/dgf_2D/dgf_b5_1D%s.eps",extension);
    C_dgf_1D->Print(filename);

	  
	TGraph *gr_dgf_b1_1D = new TGraph(dgf_count_b1_1D,angle_dgf_b1_1D,geo_factor_dgf_b1_1D);
	gr_dgf_b1_1D->SetTitle("");
	
	gr_dgf_b1_1D->SetMarkerStyle(1);
	gr_dgf_b1_1D->SetMarkerSize(1);
		//	gr_dgf_b1_1D->SetTitle("CALET Differential Geometry Factor");
	gr_dgf_b1_1D->GetXaxis()->SetTitle("Incidence Angle (degrees)");
	gr_dgf_b1_1D->GetXaxis()->CenterTitle();
	gr_dgf_b1_1D->GetYaxis()->SetTitle("Differential Geometry Factor (cm^{2}sr/(1^{o} bins))");
	gr_dgf_b1_1D->GetYaxis()->CenterTitle();
	gr_dgf_b1_1D->GetXaxis()->SetTitleOffset(1.2);
	gr_dgf_b1_1D->GetYaxis()->SetTitleOffset(1.3);
	gr_dgf_b1_1D->GetYaxis()->SetLabelFont(nicefont);
	gr_dgf_b1_1D->GetYaxis()->SetTitleFont(nicefont);
	gr_dgf_b1_1D->GetXaxis()->SetLabelFont(nicefont);
	gr_dgf_b1_1D->GetXaxis()->SetTitleFont(nicefont);
	gr_dgf_b1_1D->GetXaxis()->SetRangeUser(0,90);
	//gr_dgf_b1_1D->GetYaxis()->SetRangeUser(0,40);
	
	
	
	TCanvas *C_dgf_b1_1D = new TCanvas("C_dgf_b1_1D","1 degree 1D DGF",0,0,800,800);
	C_dgf_b1_1D->SetLeftMargin(0.139594);
	C_dgf_b1_1D->SetRightMargin(0.0596447);
	C_dgf_b1_1D->SetBottomMargin(0.11911);
	C_dgf_b1_1D->SetTopMargin(0.15);
	
	gr_dgf_b1_1D->SetLineWidth(2);
	gr_dgf_b1_1D->Draw("AL");
	gPad->Update();
	
	// Add title as text
	sprintf(title_str,"Angle %s ",angle_str);
	TPaveText *pt_b1 = new TPaveText(0.2, 0.92, 0.8, 1.0, "NDC");
	pt_b1->AddText(title_str);
	pt_b1->SetFillStyle(0);
	pt_b1->SetBorderSize(0);
	pt_b1->SetTextSize(0.03);
	pt_b1->SetTextFont(42);
	pt_b1->SetTextAlign(22);
	pt_b1->Draw();

	C_dgf_b1_1D->SetGrid();

    // writing plot files
    sprintf(filename,"plots/dgf_2D/dgf_b1_1D%s.png",extension);
    C_dgf_b1_1D->Print(filename);
    sprintf(filename,"plots/dgf_2D/dgf_b1_1D%s.eps",extension);
    C_dgf_b1_1D->Print(filename);

	
//	theApp.Run();
	
}
