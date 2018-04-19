
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TFrame.h"
#include "TFile.h"
#include <iostream>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "Riostream.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>





//----------------------------------------------------------------------------------------------------------------------------
// Creation d'un histogramme pour tracer la FWHM de reconstruction en fonction de la resolution en energie de l'absorbeur NaI
//----------------------------------------------------------------------------------------------------------------------------

void graphique_2D_Precision_Chi2()
{

TCanvas *c1 = new TCanvas("c1","c1",700,800); // Création du canvas

//Double_t number_datas = 19;




	Double_t number_protons[11]	= {1.e8, 2.e8,3.e8,4.e8, 5.e8,7.e8,1.e9,2.e9,3.e9,4.e9,5.e9};

    //bin 1 mm shift Nurbs 0.1 mm correction Moyenne Up and Down
    Double_t datas_MLEM_0_1mm_Moyenne[11]		={2.05293, 1.40637, 1.16388, 1.0135, 0.898718, 0.767922, 0.643808, 0.432282, 0.364576, 0.333195, 0.285363};

   //
    //bin 1 mm shift Nurbs 0.1 mm correction Moyenne Up and Down LINE CONE
    Double_t datas_LineCone_0_1mm_Moyenne[11]		={5.45785, 3.66188, 2.88942, 2.45234, 2.16918, 1.84977, 1.55683, 1.08817, 0.884143, 0.768433, 0.693309};


    //--------------------------------------------
    //-------------------ERROR -------------------------
    //--------------------------------------------


    //bin 1 mm shift Nurbs 0.1 mm correction Moyenne Up and Down
    Double_t datas_error_MLEM_0_1mm_Moyenne[11]		={0.0459049, 0.0314474, 0.0260252, 0.0226625, 0.0200959, 0.0171713, 0.014396, 0.00966612, 0.00815217, 0.00745047, 0.0063809};

    //bin 1 mm shift Nurbs 0.1 mm correction Moyenne Up and Down LINE CONE
    Double_t datas_error_LineCone_0_1mm_Moyenne[11]		={0.122041, 0.081882, 0.0646093, 0.054836, 0.0485044, 0.0413621, 0.0348118, 0.0243323, 0.01977, 0.0171827, 0.0155029};



	char name[80];	sprintf(name, ""); // Nom du graphique


	TH1 *frame = new TH1F("frame","",1000,9e7,1e10);
	frame->SetMinimum(0.1);
	frame->SetMaximum(10);
	frame->SetDirectory(0);
	frame->SetStats(0);
	frame->SetTitle(name);
	frame->GetXaxis()->SetTitle("Number of incident protons");
	frame->GetXaxis()->SetTickLength(0.02);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.05);
	frame->GetYaxis()->SetTitle("FRP [mm]");

	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetTitleSize(0.05);

	frame->Draw(" ");

  //gStyle->SetOptFit(0100);


	 // Échelle logarithmique sur les axes y
	c1->SetGrid();
	c1->SetLogx();
	c1->SetLogy();


	int number_datas = 11;



    TGraph *gr7 = new TGraphErrors(number_datas,number_protons,datas_MLEM_0_1mm_Moyenne,0,datas_error_MLEM_0_1mm_Moyenne); // Création du graphique

    gr7->SetMarkerColor(4); // Options de mise en forme  du graphique
    gr7->SetMarkerStyle(21);//25
    gr7->SetMarkerSize(1.3);

    gr7->Draw("P");
    TF1 * f1 = new TF1("f1","[0]*TMath::Power(x,[1])",1e8,1e10);
    f1->SetLineWidth(1);
    f1->SetLineColor(4);
    f1->SetLineStyle(2);

    f1->SetParameter(0,1);
    f1->SetParameter(1,0.5);

    gr7->Fit(f1);

    TGraph *gr8 = new TGraphErrors(number_datas,number_protons,datas_LineCone_0_1mm_Moyenne,0,datas_error_LineCone_0_1mm_Moyenne); // Création du graphique

    gr8->SetMarkerColor(2); // Options de mise en forme  du graphique
    gr8->SetMarkerStyle(22);
    gr8->SetMarkerSize(1.6);


    gr8->Draw("P");

    TF1 * f2 = new TF1("f2","[0]*TMath::Power(x,[1])",1e8,1e10);
    f2->SetLineColor(2);
    f2->SetLineWidth(1);
    f2->SetLineStyle(3);
    f2->SetParameter(0,1);
    f2->SetParameter(1,0.5);

    gr8->Fit(f2);






    leg = new TLegend(0.65,0.75,0.9,0.9);
    leg->SetFillColor(10);

    // leg->AddEntry(gr8,"Line cone - bin 1mm - shift 0.1 mm","p");
   // leg->AddEntry(gr7,"MLEM - bin 1mm - shift 0.1 mm","p");

    leg->AddEntry(gr8,"Line cone","p");
    leg->AddEntry(gr7,"MLEM","p");





	leg->Draw();
	c1->Update();


}
