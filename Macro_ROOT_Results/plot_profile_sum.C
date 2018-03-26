#include <TCanvas.h>

#include <TROOT.h>
#include "TH1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include <math.h>


void plot_profile_sum(){

  double true_gamma_1[63] = {0.07, 0.07, 0.05, 0.03, 0.07, 0.06, 0.08, 0.05, 0.11, 0.09, 0.05, 0.12, 0.12, 0.08, 0.17, 0.13, 0.14, 0.18, 0.23, 0.23, 0.13, 0.2, 0.25, 0.22, 0.27, 0.33, 0.36, 0.44, 0.44, 0.35, 0.34, 0.3, 0.38, 0.43, 0.33, 0.41, 0.45, 0.435, 0.395, 0.376, 0.36, 0.48, 0.376, 0.46, 0.34, 0.35, 0.32, 0.35, 0.24, 0.28, 0.25, 0.28, 0.3, 0.28, 0.25, 0.22, 0.15, 0.16, 0.17, 0.22, 0.14, 0.15, 0.2};
//vector <double> true_gamma (true_gamma_1, true_gamma_1+sizeof(true_gamma_1)/sizeof(double));
  double background[63] = {0.06, 0.05, 0.1, 0.1, 0.13, 0.07, 0.08, 0.1, 0.11, 0.11, 0.1, 0.09, 0.07, 0.17, 0.24, 0.12, 0.12, 0.15, 0.11, 0.11, 0.1, 0.1, 0.12, 0.2, 0.19, 0.15, 0.17, 0.1, 0.16, 0.19, 0.2, 0.15, 0.17, 0.25, 0.27, 0.22, 0.24, 0.17, 0.17, 0.15, 0.25, 0.21, 0.11, 0.18, 0.21, 0.15, 0.15, 0.15, 0.19, 0.13, 0.17, 0.15, 0.17, 0.25, 0.22, 0.16, 0.19, 0.23, 0.11,0.11, 0.16, 0.2};
//for(Int_t l = 0; l<true_gamma_1.size(); l++){
//  true_gamma.push_back(true_gamma_1.at(l));
//}

//for(Int_t ll = 0; ll<true_gamma_2.size(); ll++){
  //true_gamma.push_back(true_gamma_2.at(ll));
//}

//cout<<"true gamma size "<< true_gamma.size()<<endl;

  TH1F *profile = new TH1F("prof", "prof", 63, -180., 130.);
  TH1F *profile_back = new TH1F("prof_back", "prof_back", 63, -180., 130.);
TH1F *profile_sum = new TH1F("prof_sum", "prof_sum", 63, -180., 130.);

  for(Int_t x = 0; x <63; x++){
      profile->Fill(profile->GetBinCenter(x), true_gamma_1[x]*1e-6);
      profile_back->Fill(profile->GetBinCenter(x), background[x]*1e-6);
      profile_sum->Fill(profile->GetBinCenter(x), (true_gamma_1[x]+background[x])*1e-6);


  }

  TCanvas *c = new TCanvas(" ", " ", 650, 500);
  c->SetLeftMargin(0.15);
  TH1 *frame = c->DrawFrame(-180, 0, 120, 1.e-6);   //200,-180,130);
    //frame->SetMinimum(0);
    //frame->SetMaximum(300);//
    frame->SetDirectory(0);
    frame->SetStats(0);
    frame->SetTitle("");
    //	frame->	GetYaxis()->SetLabelSize(20);
    //frame->GetYaxis()->SetLabelSize(0.001);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitle("Beam axis [mm]");
    frame->GetXaxis()->SetTickLength(0.02);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitle("#splitline{Number of events}{per incident proton}");
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetTitleOffset(1.5);
    //frame->Draw("");
    //profile->SetLineColor(kBlue);
    //profile->Draw("same hist");
    //profile_back->SetLineColor(kBlack);
    //profile_back->Draw("same hist");
    profile_sum->SetLineColor(kBlue);
    profile_sum->SetLineWidth(2);

    profile_sum->Draw("same hist");


        c->Modified();
        c->Update();

}
