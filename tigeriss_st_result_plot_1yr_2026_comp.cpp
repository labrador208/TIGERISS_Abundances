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

FILE *outfile;

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
    
    
    
    
    
    char input_file[200];
    if(argc < 2) {
        fprintf(stderr,"Usage: %s <input_file>\n", argv[0]);
        fprintf(stderr,"Example: %s output/2026-06-02_av_5_test.txt\n", argv[0]);
        return 1;
    }
    strcpy(input_file, argv[1]);
    fprintf(stderr,"Reading: %s\n",input_file);
    
    TApplication theApp("App", &argc, argv);
    
    ifstream in_1, in_2;
    ofstream out_1;

     double Z_1[100];
     double S_av_events_int_1[100];
     double S_max_events_int_1[100];
     double S_min_events_int_1[100];
     double S_av_events_int_top_1[100];
     double S_max_events_int_top_1[100];
     double S_min_events_int_top_1[100];
     double S_av_events_int_Si2_1[100];
     double S_max_events_int_Si2_1[100];
     double S_min_events_int_Si2_1[100];
     double S_av_events_int_Si2_top_1[100];
     double S_max_events_int_Si2_top_1[100];
     double S_min_events_int_Si2_top_1[100];     
     double new_av_1y_1[100];     
     double new_av_3y_1[100];     
     double new_max_1y_1[100];  

     double Z_temp_1;    
     double S_av_events_int_temp_1;
     double S_max_events_int_temp_1;
     double S_min_events_int_temp_1;
     double S_av_events_int_top_temp_1;
     double S_max_events_int_top_temp_1;
     double S_min_events_int_top_temp_1;
     double S_av_events_int_Si2_temp_1;
     double S_max_events_int_Si2_temp_1;
     double S_min_events_int_Si2_temp_1;
     double S_av_events_int_Si2_top_temp_1;
     double S_max_events_int_Si2_top_temp_1;
     double S_min_events_int_Si2_top_temp_1; 
    
    int index_1;
    
    for (index_1 = 0; index_1 < 100; index_1++) {
        Z_1[index_1] = 0.0;
        S_av_events_int_1[index_1] = 0.0;
        S_max_events_int_1[index_1] = 0.0;
        S_min_events_int_1[index_1] = 0.0;
        S_av_events_int_top_1[index_1] = 0.0;
        S_max_events_int_top_1[index_1] = 0.0;
        S_min_events_int_top_1[index_1] = 0.0;
        S_av_events_int_Si2_1[index_1] = 0.0;
        S_max_events_int_Si2_1[index_1] = 0.0;
        S_min_events_int_Si2_1[index_1] = 0.0;
        S_av_events_int_Si2_top_1[index_1] = 0.0;
        S_max_events_int_Si2_top_1[index_1] = 0.0;
        S_min_events_int_Si2_top_1[index_1] = 0.0;   
        new_av_1y_1[index_1] = 0.0;     
        new_av_3y_1[index_1] = 0.0;     
        new_max_1y_1[index_1] = 0.0;     

    }   


        
    in_1.open(input_file);  // Use command-line argument instead of hardcoded filename
    fprintf(stderr,"Reading dataset 1: %s\n",input_file);
 
    
    index_1 = 0;
    
    while (1) {
        in_1 >> Z_temp_1 >> S_av_events_int_temp_1 >> S_max_events_int_temp_1 >> S_min_events_int_temp_1 >> S_av_events_int_top_temp_1 >> S_max_events_int_top_temp_1 >> S_min_events_int_top_temp_1 >> S_av_events_int_Si2_temp_1 >> S_max_events_int_Si2_temp_1 >> S_min_events_int_Si2_temp_1 >> S_av_events_int_Si2_top_temp_1 >> S_max_events_int_Si2_top_temp_1 >> S_min_events_int_Si2_top_temp_1; 



        fprintf(stdout,"%2f %e %e %e ",Z_temp_1,S_av_events_int_temp_1,S_max_events_int_temp_1,S_min_events_int_temp_1);
        fprintf(stdout,"%e %e %e ",S_av_events_int_top_temp_1,S_max_events_int_top_temp_1,S_min_events_int_top_temp_1);
        fprintf(stdout,"%e %e %e ",S_av_events_int_Si2_temp_1,S_max_events_int_Si2_temp_1,S_min_events_int_Si2_temp_1);
        fprintf(stdout,"%e %e %e\n",S_av_events_int_Si2_top_temp_1,S_max_events_int_Si2_top_temp_1,S_min_events_int_Si2_top_temp_1);  

        if(in_1.eof() || !in_1.good()) break;
            Z_1[index_1] = Z_temp_1;
            S_av_events_int_1[index_1] = S_av_events_int_temp_1;
            S_max_events_int_1[index_1] = S_max_events_int_temp_1;
            S_min_events_int_1[index_1] = S_min_events_int_temp_1;
            S_av_events_int_top_1[index_1] = S_av_events_int_top_temp_1;
            S_max_events_int_top_1[index_1] = S_max_events_int_top_temp_1;
            S_min_events_int_top_1[index_1] = S_min_events_int_top_temp_1;
            S_av_events_int_Si2_1[index_1] = S_av_events_int_Si2_temp_1;
            S_max_events_int_Si2_1[index_1] = S_max_events_int_Si2_temp_1;
            S_min_events_int_Si2_1[index_1] = S_min_events_int_Si2_temp_1;
            S_av_events_int_Si2_top_1[index_1] = S_av_events_int_Si2_top_temp_1;
            S_max_events_int_Si2_top_1[index_1] = S_max_events_int_Si2_top_temp_1;
            S_min_events_int_Si2_top_1[index_1] = S_min_events_int_Si2_top_temp_1; 
        
                
        new_av_1y_1[index_1] = S_av_events_int_top_1[index_1]/5.0;     
        new_av_3y_1[index_1] = S_av_events_int_top_1[index_1]*3.0/5.0;     
        new_max_1y_1[index_1] = S_max_events_int_top_1[index_1]/5.0;   

        fprintf(stdout,"%2f %e %e %e\n",Z_temp_1,new_av_1y_1[index_1],new_av_3y_1[index_1],new_max_1y_1[index_1]);        
  
        index_1++;

    }
    


    
    
    in_1.close();
    in_1.clear(); 

    TGraph *gr_tigeriss_sa_1y_1 = new TGraph(index_1,Z_1,new_av_1y_1);
    gr_tigeriss_sa_1y_1->SetTitle("It Works! (Kinda)");
   // gr_tigeriss_sa_1y_1->AddText(0.5,0.5,"It works (kinda)");
	gr_tigeriss_sa_1y_1->GetXaxis()->SetTitle("Z");
	gr_tigeriss_sa_1y_1->GetXaxis()->CenterTitle();
	gr_tigeriss_sa_1y_1->GetXaxis()->SetLabelSize(0.05);
	gr_tigeriss_sa_1y_1->GetYaxis()->SetLabelSize(0.05);
	gr_tigeriss_sa_1y_1->GetXaxis()->SetTitleSize(0.05);
	gr_tigeriss_sa_1y_1->GetYaxis()->SetTitleSize(0.05);
	gr_tigeriss_sa_1y_1->GetYaxis()->SetTitle("Counts");
	gr_tigeriss_sa_1y_1->GetYaxis()->CenterTitle();
	gr_tigeriss_sa_1y_1->GetXaxis()->SetTitleOffset(1.0);
	gr_tigeriss_sa_1y_1->GetYaxis()->SetTitleOffset(1.0);
	gr_tigeriss_sa_1y_1->GetYaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sa_1y_1->GetYaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sa_1y_1->GetXaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sa_1y_1->GetXaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sa_1y_1->GetXaxis()->SetRangeUser(5,85);
	gr_tigeriss_sa_1y_1->GetYaxis()->SetRangeUser(0.5,2E8); 
    TPaveText *pt_b1_2d = new TPaveText(0.2, 0.92, 0.8, 1.0, "NDC");
	pt_b1_2d->AddText("It Works! (Kinda)");
	pt_b1_2d->SetFillStyle(0);
	pt_b1_2d->SetBorderSize(0);
	pt_b1_2d->SetTextSize(0.03);
	pt_b1_2d->SetTextFont(42);
	pt_b1_2d->SetTextAlign(22);
	pt_b1_2d->Draw();
    TGraph *gr_tigeriss_sa_3y_1 = new TGraph(index_1,Z_1,new_av_3y_1);
    gr_tigeriss_sa_3y_1->SetTitle("");
	gr_tigeriss_sa_3y_1->GetXaxis()->SetTitle("Z");
	gr_tigeriss_sa_3y_1->GetXaxis()->CenterTitle();
	gr_tigeriss_sa_3y_1->GetXaxis()->SetLabelSize(0.05);
	gr_tigeriss_sa_3y_1->GetYaxis()->SetLabelSize(0.05);
	gr_tigeriss_sa_3y_1->GetXaxis()->SetTitleSize(0.05);
	gr_tigeriss_sa_3y_1->GetYaxis()->SetTitleSize(0.05);
	gr_tigeriss_sa_3y_1->GetYaxis()->SetTitle("Counts");
	gr_tigeriss_sa_3y_1->GetYaxis()->CenterTitle();
	gr_tigeriss_sa_3y_1->GetXaxis()->SetTitleOffset(1.0);
	gr_tigeriss_sa_3y_1->GetYaxis()->SetTitleOffset(1.0);
	gr_tigeriss_sa_3y_1->GetYaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sa_3y_1->GetYaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sa_3y_1->GetXaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sa_3y_1->GetXaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sa_3y_1->GetXaxis()->SetRangeUser(5,85);
	gr_tigeriss_sa_3y_1->GetYaxis()->SetRangeUser(0.5,2E8); 

     double Z_2[100];
     double S_av_events_int_2[100];
     double S_max_events_int_2[100];
     double S_min_events_int_2[100];
     double S_av_events_int_top_2[100];
     double S_max_events_int_top_2[100];
     double S_min_events_int_top_2[100];
     double S_av_events_int_Si2_2[100];
     double S_max_events_int_Si2_2[100];
     double S_min_events_int_Si2_2[100];
     double S_av_events_int_Si2_top_2[100];
     double S_max_events_int_Si2_top_2[100];
     double S_min_events_int_Si2_top_2[100];     
     double new_av_1y_2[100];     
     double new_av_3y_2[100];     
     double new_max_1y_2[100];  

     double Z_temp_2;    
     double S_av_events_int_temp_2;
     double S_max_events_int_temp_2;
     double S_min_events_int_temp_2;
     double S_av_events_int_top_temp_2;
     double S_max_events_int_top_temp_2;
     double S_min_events_int_top_temp_2;
     double S_av_events_int_Si2_temp_2;
     double S_max_events_int_Si2_temp_2;
     double S_min_events_int_Si2_temp_2;
     double S_av_events_int_Si2_top_temp_2;
     double S_max_events_int_Si2_top_temp_2;
     double S_min_events_int_Si2_top_temp_2; 
    
    int index_2;
    
    for (index_2 = 0; index_2 < 100; index_2++) {
        Z_2[index_2] = 0.0;
        S_av_events_int_2[index_2] = 0.0;
        S_max_events_int_2[index_2] = 0.0;
        S_min_events_int_2[index_2] = 0.0;
        S_av_events_int_top_2[index_2] = 0.0;
        S_max_events_int_top_2[index_2] = 0.0;
        S_min_events_int_top_2[index_2] = 0.0;
        S_av_events_int_Si2_2[index_2] = 0.0;
        S_max_events_int_Si2_2[index_2] = 0.0;
        S_min_events_int_Si2_2[index_2] = 0.0;
        S_av_events_int_Si2_top_2[index_2] = 0.0;
        S_max_events_int_Si2_top_2[index_2] = 0.0;
        S_min_events_int_Si2_top_2[index_2] = 0.0;   
        new_av_1y_2[index_2] = 0.0;     
        new_av_3y_2[index_2] = 0.0;     
        new_max_1y_2[index_2] = 0.0;     

    }   


        
    // Disabled: second dataset (requires additional files)
    // in_1.open("output/2023-07-15_av_4_150x60x42.txt");
 
    
    index_2 = 0;
    
    while (in_1.good() && in_1.is_open()) {
        in_1 >> Z_temp_2 >> S_av_events_int_temp_2 >> S_max_events_int_temp_2 >> S_min_events_int_temp_2 >> S_av_events_int_top_temp_2 >> S_max_events_int_top_temp_2 >> S_min_events_int_top_temp_2 >> S_av_events_int_Si2_temp_2 >> S_max_events_int_Si2_temp_2 >> S_min_events_int_Si2_temp_2 >> S_av_events_int_Si2_top_temp_2 >> S_max_events_int_Si2_top_temp_2 >> S_min_events_int_Si2_top_temp_2; 



        fprintf(stdout,"%2f %e %e %e ",Z_temp_2,S_av_events_int_temp_2,S_max_events_int_temp_2,S_min_events_int_temp_2);
        fprintf(stdout,"%e %e %e ",S_av_events_int_top_temp_2,S_max_events_int_top_temp_2,S_min_events_int_top_temp_2);
        fprintf(stdout,"%e %e %e ",S_av_events_int_Si2_temp_2,S_max_events_int_Si2_temp_2,S_min_events_int_Si2_temp_2);
        fprintf(stdout,"%e %e %e\n",S_av_events_int_Si2_top_temp_2,S_max_events_int_Si2_top_temp_2,S_min_events_int_Si2_top_temp_2);  

        if(in_1.eof() || !in_1.good()) break;
            Z_2[index_2] = Z_temp_2;
            S_av_events_int_2[index_2] = S_av_events_int_temp_2;
            S_max_events_int_2[index_2] = S_max_events_int_temp_2;
            S_min_events_int_2[index_2] = S_min_events_int_temp_2;
            S_av_events_int_top_2[index_2] = S_av_events_int_top_temp_2;
            S_max_events_int_top_2[index_2] = S_max_events_int_top_temp_2;
            S_min_events_int_top_2[index_2] = S_min_events_int_top_temp_2;
            S_av_events_int_Si2_2[index_2] = S_av_events_int_Si2_temp_2;
            S_max_events_int_Si2_2[index_2] = S_max_events_int_Si2_temp_2;
            S_min_events_int_Si2_2[index_2] = S_min_events_int_Si2_temp_2;
            S_av_events_int_Si2_top_2[index_2] = S_av_events_int_Si2_top_temp_2;
            S_max_events_int_Si2_top_2[index_2] = S_max_events_int_Si2_top_temp_2;
            S_min_events_int_Si2_top_2[index_2] = S_min_events_int_Si2_top_temp_2; 
        
                
        new_av_1y_2[index_2] = S_av_events_int_top_2[index_2]/5.0;     
        new_av_3y_2[index_2] = S_av_events_int_top_2[index_2]*3.0/5.0;     
        new_max_1y_2[index_2] = S_max_events_int_top_2[index_2]/5.0;   

        fprintf(stdout,"%2f %e %e %e\n",Z_temp_2,new_av_1y_2[index_2],new_av_3y_2[index_2],new_max_1y_2[index_2]);        
  
        index_2++;

    }
    


    
    
    in_1.close();
    in_1.clear(); 

    TGraph *gr_tigeriss_sa_1y_2 = new TGraph(index_2,Z_2,new_av_1y_2);    
    TGraph *gr_tigeriss_sa_3y_2 = new TGraph(index_2,Z_2,new_av_3y_2);

     double Z_3[100];
     double S_av_events_int_3[100];
     double S_max_events_int_3[100];
     double S_min_events_int_3[100];
     double S_av_events_int_top_3[100];
     double S_max_events_int_top_3[100];
     double S_min_events_int_top_3[100];
     double S_av_events_int_Si2_3[100];
     double S_max_events_int_Si2_3[100];
     double S_min_events_int_Si2_3[100];
     double S_av_events_int_Si2_top_3[100];
     double S_max_events_int_Si2_top_3[100];
     double S_min_events_int_Si2_top_3[100];     
     double new_av_1y_3[100];     
     double new_av_3y_3[100];     
     double new_max_1y_3[100];  

     double Z_temp_3;    
     double S_av_events_int_temp_3;
     double S_max_events_int_temp_3;
     double S_min_events_int_temp_3;
     double S_av_events_int_top_temp_3;
     double S_max_events_int_top_temp_3;
     double S_min_events_int_top_temp_3;
     double S_av_events_int_Si2_temp_3;
     double S_max_events_int_Si2_temp_3;
     double S_min_events_int_Si2_temp_3;
     double S_av_events_int_Si2_top_temp_3;
     double S_max_events_int_Si2_top_temp_3;
     double S_min_events_int_Si2_top_temp_3; 
    
    int index_3;
    
    for (index_3 = 0; index_3 < 100; index_3++) {
        Z_3[index_3] = 0.0;
        S_av_events_int_3[index_3] = 0.0;
        S_max_events_int_3[index_3] = 0.0;
        S_min_events_int_3[index_3] = 0.0;
        S_av_events_int_top_3[index_3] = 0.0;
        S_max_events_int_top_3[index_3] = 0.0;
        S_min_events_int_top_3[index_3] = 0.0;
        S_av_events_int_Si2_3[index_3] = 0.0;
        S_max_events_int_Si2_3[index_3] = 0.0;
        S_min_events_int_Si2_3[index_3] = 0.0;
        S_av_events_int_Si2_top_3[index_3] = 0.0;
        S_max_events_int_Si2_top_3[index_3] = 0.0;
        S_min_events_int_Si2_top_3[index_3] = 0.0;   
        new_av_1y_3[index_3] = 0.0;     
        new_av_3y_3[index_3] = 0.0;     
        new_max_1y_3[index_3] = 0.0;     

    }   


        
    // Disabled: third dataset (requires additional files)
    // in_1.open("output/2026-01-25_av_4_9x9det.txt");
    // fprintf(stderr,"Reading dataset 3: %s\n",input_file);
    
    index_3 = 0;  // Keep as 0 if dataset disabled
    
    while (in_1.good() && in_1.is_open()) {
        in_1 >> Z_temp_3 >> S_av_events_int_temp_3 >> S_max_events_int_temp_3 >> S_min_events_int_temp_3 >> S_av_events_int_top_temp_3 >> S_max_events_int_top_temp_3 >> S_min_events_int_top_temp_3 >> S_av_events_int_Si2_temp_3 >> S_max_events_int_Si2_temp_3 >> S_min_events_int_Si2_temp_3 >> S_av_events_int_Si2_top_temp_3 >> S_max_events_int_Si2_top_temp_3 >> S_min_events_int_Si2_top_temp_3; 



        fprintf(stdout,"%2f %e %e %e ",Z_temp_3,S_av_events_int_temp_3,S_max_events_int_temp_3,S_min_events_int_temp_3);
        fprintf(stdout,"%e %e %e ",S_av_events_int_top_temp_3,S_max_events_int_top_temp_3,S_min_events_int_top_temp_3);
        fprintf(stdout,"%e %e %e ",S_av_events_int_Si2_temp_3,S_max_events_int_Si2_temp_3,S_min_events_int_Si2_temp_3);
        fprintf(stdout,"%e %e %e\n",S_av_events_int_Si2_top_temp_3,S_max_events_int_Si2_top_temp_3,S_min_events_int_Si2_top_temp_3);  

        if(in_1.eof() || !in_1.good()) break;
            Z_3[index_3] = Z_temp_3;
            S_av_events_int_3[index_3] = S_av_events_int_temp_3;
            S_max_events_int_3[index_3] = S_max_events_int_temp_3;
            S_min_events_int_3[index_3] = S_min_events_int_temp_3;
            S_av_events_int_top_3[index_3] = S_av_events_int_top_temp_3;
            S_max_events_int_top_3[index_3] = S_max_events_int_top_temp_3;
            S_min_events_int_top_3[index_3] = S_min_events_int_top_temp_3;
            S_av_events_int_Si2_3[index_3] = S_av_events_int_Si2_temp_3;
            S_max_events_int_Si2_3[index_3] = S_max_events_int_Si2_temp_3;
            S_min_events_int_Si2_3[index_3] = S_min_events_int_Si2_temp_3;
            S_av_events_int_Si2_top_3[index_3] = S_av_events_int_Si2_top_temp_3;
            S_max_events_int_Si2_top_3[index_3] = S_max_events_int_Si2_top_temp_3;
            S_min_events_int_Si2_top_3[index_3] = S_min_events_int_Si2_top_temp_3; 
        
                
        new_av_1y_3[index_3] = S_av_events_int_top_3[index_3]/5.0;     
        new_av_3y_3[index_3] = S_av_events_int_top_3[index_3]*3.0/5.0;     
        new_max_1y_3[index_3] = S_max_events_int_top_3[index_3]/5.0;   

        fprintf(stdout,"%2f %e %e %e\n",Z_temp_3,new_av_1y_3[index_3],new_av_3y_3[index_3],new_max_1y_3[index_3]);        
  
        index_3++;

    }
    


    
    
    in_1.close();
    in_1.clear(); 

    TGraph *gr_tigeriss_sa_1y_3 = new TGraph(index_3,Z_3,new_av_1y_3);   
    gr_tigeriss_sa_1y_3->SetTitle("");
	gr_tigeriss_sa_1y_3->GetXaxis()->SetTitle("Z");
	gr_tigeriss_sa_1y_3->GetXaxis()->CenterTitle();
	gr_tigeriss_sa_1y_3->GetXaxis()->SetLabelSize(0.05);
	gr_tigeriss_sa_1y_3->GetYaxis()->SetLabelSize(0.05);
	gr_tigeriss_sa_1y_3->GetXaxis()->SetTitleSize(0.05);
	gr_tigeriss_sa_1y_3->GetYaxis()->SetTitleSize(0.05);
	gr_tigeriss_sa_1y_3->GetYaxis()->SetTitle("Counts");
	gr_tigeriss_sa_1y_3->GetYaxis()->CenterTitle();
	gr_tigeriss_sa_1y_3->GetXaxis()->SetTitleOffset(1.0);
	gr_tigeriss_sa_1y_3->GetYaxis()->SetTitleOffset(1.0);
	gr_tigeriss_sa_1y_3->GetYaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sa_1y_3->GetYaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sa_1y_3->GetXaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sa_1y_3->GetXaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sa_1y_3->GetXaxis()->SetRangeUser(5,85);
	gr_tigeriss_sa_1y_3->GetYaxis()->SetRangeUser(0.5,2E8); 
    TGraph *gr_tigeriss_sa_3y_3 = new TGraph(index_3,Z_3,new_av_3y_3);
    
    
    char st_tigeriss_in[200] = "tables/st_tigeriss_AP_2022.txt";  // Workspace-relative path

    
     int counts_UH_tiger=0.0;
     int counts_UH_CALET_full=0.0;
     int counts_UH_CALET_UH=0.0;
    
     double Z[100];
     double Z_st_plot[100];
     double Z_st_plot_err_hi[100];
     double Z_st_plot_err_lo[100];
     double st_nathan[100];
     double st_nathan_err_hi[100];
     double st_nathan_err_lo[100];
     double st_nathan_plot[100];
     double st_nathan_plot_err_hi[100];
     double st_nathan_plot_err_lo[100];
     double tigeriss_sm[100];
     double tigeriss_sa[100];
     double tigeriss_sa_1m[100];
     double tigeriss_sa_3y[100];
        
    int index=0;
    int st_plot=0;
    
    in_1.open(st_tigeriss_in);
    fprintf(stderr,"%s file\n",st_tigeriss_in);
    
    while (1) {
        in_1 >> Z[index] >> st_nathan[index] >> st_nathan_err_hi[index] >> st_nathan_err_lo[index] >> tigeriss_sa[index] >> tigeriss_sm[index];

        if(in_1.eof() || !in_1.good()) break;

        fprintf(stderr,"Z: %.0f ST Nathan: %.2f +%.2f -%.2f T_ISS SA: %.2f T_ISS SM: %.2f\n",
                Z[index],st_nathan[index],st_nathan_err_hi[index],st_nathan_err_lo[index],tigeriss_sa[index],tigeriss_sm[index]);
        tigeriss_sa_1m[index] = tigeriss_sa[index]/12.0;
        tigeriss_sa_3y[index] = tigeriss_sa[index]*3.0;

        if(Z[index] > 16) {
           st_nathan_plot[st_plot] = st_nathan[index];
           st_nathan_plot_err_hi[st_plot] = st_nathan_err_hi[index];
           st_nathan_plot_err_lo[st_plot] = st_nathan_err_lo[index];
           Z_st_plot[st_plot] = Z[index];
           Z_st_plot_err_hi[st_plot] = 0.0;
           Z_st_plot_err_lo[st_plot] = 0.0;
           st_plot++;
        }    
        index++;

    }
    


    
    
    in_1.close();
    in_1.clear();
    



    
    TGraphAsymmErrors *gr_st_nathan = new TGraphAsymmErrors(st_plot,Z_st_plot,st_nathan_plot,Z_st_plot_err_lo,Z_st_plot_err_hi,st_nathan_plot_err_lo,st_nathan_plot_err_hi);

    TGraph *gr_tigeriss_sa = new TGraph(index,Z,tigeriss_sa);
    
    TGraph *gr_tigeriss_sm = new TGraph(index,Z,tigeriss_sm);

	gr_tigeriss_sm->SetMarkerStyle(1);
	gr_tigeriss_sm->SetMarkerSize(1);
    //	gr_tigeriss_sm->SetTitle("TIGERISS Solar Maximum Modulation - 1 Year");
	gr_tigeriss_sm->SetTitle("It works (kinda)");
	gr_tigeriss_sm->GetXaxis()->SetTitle("Z");
	gr_tigeriss_sm->GetXaxis()->CenterTitle();
	gr_tigeriss_sm->GetYaxis()->SetTitle("Counts");
	gr_tigeriss_sm->GetYaxis()->CenterTitle();
	gr_tigeriss_sm->GetXaxis()->SetTitleOffset(1.2);
	gr_tigeriss_sm->GetYaxis()->SetTitleOffset(1.3);
	gr_tigeriss_sm->GetYaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sm->GetYaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sm->GetXaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sm->GetXaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sm->GetXaxis()->SetRangeUser(5,85);
	gr_tigeriss_sm->GetYaxis()->SetRangeUser(0.5,1E8);
	
	TGraph *gr_tigeriss_sa_1m = new TGraph(index,Z,tigeriss_sa_1m);
	
    TGraph *gr_tigeriss_sa_3y = new TGraph(index,Z,tigeriss_sa_3y);

	gr_tigeriss_sa_3y->SetMarkerStyle(1);
	gr_tigeriss_sa_3y->SetMarkerSize(1);
    //	gr_tigeriss_sa_3y->SetTitle("TIGERISS Solar Average Modulation - 3 Years");
	gr_tigeriss_sa_3y->SetTitle("");
	gr_tigeriss_sa_3y->GetXaxis()->SetTitle("Z");
	gr_tigeriss_sa_3y->GetXaxis()->CenterTitle();
	gr_tigeriss_sa_3y->GetXaxis()->SetLabelSize(0.05);
	gr_tigeriss_sa_3y->GetYaxis()->SetLabelSize(0.05);
	gr_tigeriss_sa_3y->GetXaxis()->SetTitleSize(0.05);
	gr_tigeriss_sa_3y->GetYaxis()->SetTitleSize(0.05);
	gr_tigeriss_sa_3y->GetYaxis()->SetTitle("Counts");
	gr_tigeriss_sa_3y->GetYaxis()->CenterTitle();
	gr_tigeriss_sa_3y->GetXaxis()->SetTitleOffset(1.0);
	gr_tigeriss_sa_3y->GetYaxis()->SetTitleOffset(1.0);
	gr_tigeriss_sa_3y->GetYaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sa_3y->GetYaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sa_3y->GetXaxis()->SetLabelFont(nicefont);
	gr_tigeriss_sa_3y->GetXaxis()->SetTitleFont(nicefont);
	gr_tigeriss_sa_3y->GetXaxis()->SetRangeUser(5,85);
	gr_tigeriss_sa_3y->GetYaxis()->SetRangeUser(0.5,2E8);
    
  	
	TCanvas *C_st_tigeriss = new TCanvas("C_st_tigeriss","Saturation Curves",0,0,1600,800);
	C_st_tigeriss->SetLeftMargin(0.139594);
	C_st_tigeriss->SetRightMargin(0.0596447);
	C_st_tigeriss->SetBottomMargin(0.11911);
	C_st_tigeriss->SetTopMargin(0.0811518);
        C_st_tigeriss->SetGrid();
        C_st_tigeriss->SetLogy();
 
    gr_tigeriss_sa_3y_1->Draw("AL");
    gr_tigeriss_sa_3y_1->SetLineColor(4);
    gr_tigeriss_sa_3y_1->SetLineWidth(2);  
    
    //gr_tigeriss_sa_3y->Draw("AL");
    gr_tigeriss_sa_3y->SetLineColor(6);
    gr_tigeriss_sa_3y->SetLineWidth(2);
    
    //gr_tigeriss_sa_1y_3->Draw("AL");
    gr_tigeriss_sa_1y_3->SetLineColor(4);
    gr_tigeriss_sa_1y_3->SetLineWidth(2);  
    
    

    //gr_tigeriss_sa_1m->Draw("L");
    gr_tigeriss_sa_1m->SetLineColor(6);
    gr_tigeriss_sa_1m->SetLineWidth(2);


    //gr_tigeriss_sa_3y_1->Draw("L");
    //gr_tigeriss_sa_3y_1->SetLineColor(1);
    //gr_tigeriss_sa_3y_1->SetLineWidth(2);
    
    //gr_tigeriss_sa_3y_2->Draw("L");
    gr_tigeriss_sa_3y_2->SetLineColor(3);
    gr_tigeriss_sa_3y_2->SetLineWidth(2);
    
         gr_tigeriss_sa->Draw("L");
    gr_tigeriss_sa->SetLineColor(6);
    gr_tigeriss_sa->SetLineWidth(2);

    
    gr_tigeriss_sa_1y_1->Draw("L");
    gr_tigeriss_sa_1y_1->SetLineColor(1);
    gr_tigeriss_sa_1y_1->SetLineWidth(2);

    //gr_tigeriss_sa_1y_2->Draw("L");
    gr_tigeriss_sa_1y_2->SetLineColor(3);
    gr_tigeriss_sa_1y_2->SetLineWidth(2); 
    
    gr_st_nathan->SetMarkerColor(2);
    gr_st_nathan->SetMarkerStyle(21);
    gr_st_nathan->SetMarkerSize(1);
    gr_st_nathan->SetLineColor(2);
    gr_st_nathan->SetLineWidth(2);
    gr_st_nathan->Draw("P"); 
    

    


	TLegend *leg_st_tigeriss = new TLegend(0.392596,0.695604,0.81074,0.881319);
	leg_st_tigeriss->AddEntry(gr_st_nathan,"SuperTIGER Walsh","p");
	//leg_st_tigeriss->AddEntry(gr_tigeriss_sa_1m,"TIGERISS Solar Average 1 Month","l");
	leg_st_tigeriss->AddEntry(gr_tigeriss_sa,"TIGERISS Proposal Solar Average 1 Year","l");
	leg_st_tigeriss->AddEntry(gr_tigeriss_sa_1y_1,"TIGERISS SOX Solar Average 1 Year","l");
	leg_st_tigeriss->AddEntry(gr_tigeriss_sa_3y_1,"TIGERISS SOX Solar Average 3 Years","l");
	//leg_st_tigeriss->AddEntry(gr_tigeriss_sa_1y_2,"TIGERISS EFU Standard Solar Average 1 Year","l");
	//leg_st_tigeriss->AddEntry(gr_tigeriss_sa_1y_3,"TIGERISS EFU Wide Solar Average 1 Year","l");	
    
    leg_st_tigeriss->SetTextSize(0.05);
    leg_st_tigeriss->SetFillColor(10);
	leg_st_tigeriss->SetTextFont(nicefont);
	leg_st_tigeriss->SetBorderSize(0); // remove ugly border from stats
    leg_st_tigeriss->Draw();
    
	theApp.Run();
    
}

    


