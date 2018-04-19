#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <Riostream.h>
#include <TH3F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TGraph.h>


void deviation_histos(){

gROOT->Reset();
//gROOT->SetStyle("Plain");
//gROOT->ForceStyle();

//GLOBAL HISTO OPTIONS

gStyle -> SetStatW(0.28);
gStyle -> SetStatH(0.13);
gStyle -> SetStatColor(0);
gStyle -> SetStatX(0.87);
gStyle -> SetStatY(0.85);
gStyle -> SetStatFont(0);
gStyle -> SetOptStat(111);
gStyle -> SetPalette(1);

TFile *output_histos = new TFile("out_histos.root", "RECREATE");

TH1F *line_cone = new TH1F("LC", "LC", 40, -20, 20);
TH1F *MLEM = new TH1F("MLEM", "MLEM", 40, -20, 20);

Double_t positions[40];
Double_t contents_LC[40] = {0., 0., 0., 0., 2., 6., 5., 7., 9., 7., 19., 17., 17., 39., 43., 50., 57., 70., 64., 72., 94., 99., 59., 54., 60., 39., 27., 24., 15., 14., 7., 4., 6., 0., 0., 3., 2., 2., 0., 0.};
Double_t contents_MLEM[40] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 6., 13., 49., 85., 140., 186., 179., 148., 112., 45., 22., 7., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

for(Int_t l = 0; l<40; l++){
  positions[l] = -20. + l;
  line_cone->Fill(positions[l], contents_LC[l]);
  MLEM->Fill(positions[l], contents_MLEM[l]);

}

TCanvas *c = new TCanvas(" line_cone", " line_cone", 650, 500);
c->SetTitle(" ");
TCanvas *c1 = new TCanvas(" MLEM", " MLEM", 650, 500);
c1->SetTitle(" ");

line_cone->SetLineColor(kBlue);
MLEM->SetLineColor(kBlue);
line_cone->SetLineWidth(2);
MLEM->SetLineWidth(2);
line_cone->GetXaxis()->SetTitleSize(0.05);
MLEM->GetXaxis()->SetTitleSize(0.05);
line_cone->GetXaxis()->SetTickLength(0.02);
line_cone->GetXaxis()->SetLabelSize(0.05);
line_cone->GetYaxis()->SetTitleSize(0.05);
MLEM->GetYaxis()->SetTitleSize(0.05);
line_cone->GetYaxis()->SetLabelSize(0.05);
MLEM->GetYaxis()->SetLabelSize(0.05);
line_cone->GetXaxis()->SetTitle("Falloff position deviation [mm]");
line_cone->GetYaxis()->SetTitle("Entries");
MLEM->GetXaxis()->SetTitle("Falloff position deviation [mm]");
MLEM->GetYaxis()->SetTitle("Entries");
MLEM->SetTitle(" ");
line_cone->SetTitle(" ");

line_cone->SetStats(kFALSE);
MLEM->SetStats(kFALSE);


c->cd();
line_cone->Draw("hist");
c1->cd();
MLEM->Draw("hist");

c->Modified();
c->Update();
c1->Modified();
c1->Update();

c->Write("LC_deviations");
c1->Write("MLEM_deviations");

output_histos->Close();


}
