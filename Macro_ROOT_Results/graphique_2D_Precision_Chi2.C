
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
#include "TGraph.h"
#include "TGraphErrors.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>





//----------------------------------------------------------------------------------------------------------------------------
// Creation d'un histogramme pour tracer la FWHM de reconstruction en fonction de la resolution en energie de l'absorbeur NaI
//----------------------------------------------------------------------------------------------------------------------------

void graphique_2D()
{

TCanvas *c1 = new TCanvas("c1","c1",700,800); // Création du canvas

//Double_t number_datas = 19;




	Double_t number_protons[11]	= {1.e8, 2.e8,3.e8,4.e8, 5.e8,7.e8,1.e9,2.e9,3.e9,4.e9,5.e9};

    //bin 1 mm shift Nurbs 0.1 mm correction Moyenne Up and Down
    Double_t datas_MLEM_0_1mm_Moyenne[11]		={ 0.665094, 0.42195, 0.37683, 0.297494, 0.26716, 0.220854, 0.184756, 0.112511, 0.0893288, 0.0736367, 0.0609205};

   //
    //bin 1 mm shift Nurbs 0.1 mm correction Moyenne Up and Down LINE CONE
    Double_t datas_LineCone_0_1mm_Moyenne[11]		={2.27691, 1.34725, 1.0893, 0.873027, 0.803503, 0.627224, 0.480087, 0.348813, 0.302566, 0.23498, 0.210634};


    //--------------------------------------------
    //-------------------ERROR -------------------------
    //--------------------------------------------


    //bin 1 mm shift Nurbs 0.1 mm correction Moyenne Up and Down
    Double_t datas_error_MLEM_0_1mm_Moyenne[11]		={0.0148719, 0.00943509, 0.00842618, 0.00665216, 0.00597387, 0.00493845, 0.00413128, 0.00251583, 0.00199745, 0.00164657, 0.00136222};

    //bin 1 mm shift Nurbs 0.1 mm correction Moyenne Up and Down LINE CONE
    Double_t datas_error_LineCone_0_1mm_Moyenne[11]		={0.0509133, 0.0301253, 0.0243575, 0.0195215, 0.0179669, 0.0140251, 0.0107351, 0.0077997, 0.00676557, 0.00525431, 0.00470992};



	char name[80];	sprintf(name, " "); // Nom du graphique



	TH1 *frame = new TH1F("frame","",1000,9e7,1e10);
	frame->SetMinimum(0.1);
	frame->SetMaximum(5);
	frame->SetDirectory(0);
	frame->SetStats(0);
	frame->SetTitle(name);
	frame->GetXaxis()->SetTitle("Number of incident protons");
	frame->GetXaxis()->SetTickLength(0.02);
	frame->GetXaxis()->SetLabelSize(0.03);
	frame->GetYaxis()->SetTitle("Camera precision [mm]");
	frame->GetYaxis()->SetRangeUser(1e-1, 1e1);

	frame->GetYaxis()->SetLabelSize(0.03);
	frame->Draw(" ");


	 // Échelle logarithmique sur les axes y
	c1->SetGrid();
	c1->SetLogx();
	c1->SetLogy();


	int number_datas = 11;

    TGraph *gr7 = new TGraphErrors(number_datas,number_protons,datas_MLEM_0_1mm_Moyenne,0,datas_error_MLEM_0_1mm_Moyenne); // Création du graphique

    gr7->SetMarkerColor(6); // Options de mise en forme  du graphique
    gr7->SetMarkerStyle(21);//25
    gr7->SetMarkerSize(1);

    gr7->Draw("P");

    TGraph *gr8 = new TGraphErrors(number_datas,number_protons,datas_LineCone_0_1mm_Moyenne,0,datas_error_LineCone_0_1mm_Moyenne); // Création du graphique

    gr8->SetMarkerColor(1); // Options de mise en forme  du graphique
    gr8->SetMarkerStyle(26);
    gr8->SetMarkerSize(1);

    gr8->Draw("P");



    TLegend *leg = new TLegend(0.3,0.75,0.85,0.89);
    leg->SetFillColor(10);


    leg->AddEntry(gr8,"Line cone - bin 1mm - shift 0.1 mm","p");
    leg->AddEntry(gr7,"MLEM - bin 1mm - shift 0.1 mm","p");





	leg->Draw();
	c1->Update();


}
