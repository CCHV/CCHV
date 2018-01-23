
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

//#include <math>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>





//----------------------------------------------------------------------------------------------------------------------------
// Creation d'un histogramme pour tracer la FWHM de reconstruction en fonction de la resolution en energie de l'absorbeur NaI
//----------------------------------------------------------------------------------------------------------------------------

void graphique_2D_Efficacite_CC()
{

TCanvas *c1 = new TCanvas("c1","c1",700,800); // Création du canvas

//Double_t number_datas = 19;




	Double_t Energy_abscisse [19]	= {-300,-200,-150,-100,-80,-60,-40,-20,-10,0,10,20,40,60,80,100,150,200,300};

	//Avec Cut energy Si 50 keV et BGO 100 keV
	Double_t datas_300keV [20]	= {8394, 9471, 8346, 6356, 5552, 5136, 4702, 4585, 4607, 4529, 4493, 4451, 4839, 5089, 5647, 6390, 8364, 9487, 8413} ; //300 keV
	Double_t datas_500kev [20]		={13333, 20368, 23288, 23482, 22281, 20490, 18807, 17929, 17416, 17313, 17683, 17619, 18924, 20526, 22196, 23520, 23108, 20372, 13431};// 500 kev
	Double_t datas_1Mev[20]	= {9445, 18301, 24931, 31245, 33388, 35296, 37258, 38002, 38008, 38384, 38232, 38021, 37205, 35570, 33471, 31061, 24808, 18343, 9598}   ; // 1 MeV
	Double_t datas_2Mev [20]	=  {6144, 12638, 18884, 27231, 30670, 33715, 36735, 38710, 39135, 39267, 39583, 38691, 37380, 34377, 30562, 27051, 18893, 12580, 6128}; // 2 MeV
    Double_t datas_4Mev [20]	=  {4097, 8120, 13049, 20169, 23438, 27023, 29471, 31723, 31670, 32395, 31868, 31404, 29947, 26794, 23748, 20208, 12929, 8082, 4072}; // 4 MeV

    Double_t datas_6Mev [20]	=  {3427, 6890, 11149, 17901, 20680, 23544, 26003, 27605, 28136, 28107, 28268, 27654, 25613, 23431, 20591, 17654, 10920, 6803, 3474}; // 6 MeV

    //Sans aucun Cut énergie : Si 0 keV et BGO 0 keV
  /*  Double_t datas_300keV [19]	= { 15613, 29301, 39674, 52957, 58117, 62280, 66329, 68280, 69533, 69090, 68703, 68353, 65598, 62375, 57524, 52683, 40086, 29266, 15623} ; // 300 keV
    Double_t datas_500kev [19]		={13777, 26258, 37279, 50398, 55703, 60407, 63904, 66698, 67294, 67664, 67691, 66622, 64324, 60400, 55458, 50101, 37185, 26296, 13839}; // 500 kev
    Double_t datas_1Mev[19]	= {9870, 19827, 28999, 40547, 45297, 49288, 53738, 55867, 56374, 56727, 56481, 55632, 53474, 49697, 45368, 40480, 28756, 19876, 10061}   ; // 1 Mev
    Double_t datas_2Mev [19]	=  {6551, 13305, 20377, 30405, 34491, 38189, 41996, 44327, 44630, 45210, 45163, 44350, 42573, 39127, 34372, 30274, 20471, 13334, 6543}; // 2 MeV
        Double_t datas_4Mev [19]	=  {4341, 8456, 13460, 21169, 24686, 28286, 30743, 33250, 33215, 33868, 33612, 32961, 31398, 28147, 24991, 21265, 13377, 8388, 4314 }; // 4 MeV
    Double_t datas_6Mev [19]	=  {3659, 7165, 11594, 18370, 21365, 24191, 26696, 28389, 28988, 28932, 29103, 28478, 26399, 24184, 21298, 18262, 11309, 7146, 3744}; // 6 MeV
*/


    //Différence de cut  à 300 keV
   /* Double_t datas_LYSO [19]	= { 15613, 29301, 39674, 52957, 58117, 62280, 66329, 68280, 69533, 69090, 68703, 68353, 65598, 62375, 57524, 52683, 40086, 29266, 15623} ; // sans cut
    Double_t datas_BGO [19]		={ 8744, 9943, 8927, 7086, 6351, 5894, 5521, 5447, 5462, 5371, 5338, 5276, 5676, 5886, 6407, 7113, 8954, 10001, 8717}; // cut Si 50 keV BGO 0 KeV
    Double_t datas_LaBr3[19]	= {15345, 28885, 39205, 52289, 57435, 61605, 65566, 67470, 68725, 68313, 67964, 67590, 64865, 61651, 56858, 52067, 39602, 28856, 15353}   ; //cut Si 0 keV BGO 100 KeV
    Double_t datas_LYSO1 [19]	=  {6551, 13305, 20377, 30405, 34491, 38189, 41996, 44327, 44630, 45210, 45163, 44350, 42573, 39127, 34372, 30274, 20471, 13334, 6543}; // rien en rapport*/


    //just cut BGO 100 keV at 300 keV
   // 15345, 28885, 39205, 52289, 57435, 61605, 65566, 67470, 68725, 68313, 67964, 67590, 64865, 61651, 56858, 52067, 39602, 28856, 15353
    //just cut Si 50 keV at 300 keV
  //  8744, 9943, 8927, 7086, 6351, 5894, 5521, 5447, 5462, 5371, 5338, 5276, 5676, 5886, 6407, 7113, 8954, 10001, 8717




	// Calcul des incertitudes
	Double_t incertitude_300kev [20];
    Double_t incertitude_500kev [20];
    Double_t incertitude_1Mev [20];
    Double_t incertitude_2Mev [20];
    Double_t incertitude_4Mev [20];
    Double_t incertitude_6Mev [20];


    for (Int_t i=0;i<19;i++)
    {
        datas_300keV[i] = datas_300keV[i]*(1/(9*pow(10,7)*2)); // simulation sur 2 pi et 9 simu de 1e7
        datas_500kev[i] = datas_500kev[i]*(1/(9*pow(10,7)*2));
        datas_1Mev[i] = datas_1Mev[i]*(1/(9*pow(10,7)*2));
        datas_2Mev[i] = datas_2Mev[i]*(1/(9*pow(10,7)*2));
        datas_4Mev[i] = datas_4Mev[i]*(1/(9*pow(10,7)*2));
        datas_6Mev[i] = datas_6Mev[i]*(1/(9*pow(10,7)*2));
    }



	char name[80];	sprintf(name, "Detection efficiency vs source position"); // Nom du graphique


	TH1 *frame = new TH1F("frame","",1000,-300,300);
	frame->SetMinimum(1e-5);
	frame->SetMaximum(5e-4);
	frame->SetDirectory(0);
	frame->SetStats(0);
	frame->SetTitle(name);
	frame->GetXaxis()->SetTitle("position [mm]");
	frame->GetXaxis()->SetTickLength(0.02);
	frame->GetXaxis()->SetLabelSize(0.03);
	frame->GetYaxis()->SetTitle("Efficiency");

	frame->GetYaxis()->SetLabelSize(0.03);
	frame->Draw(" ");


	 // Échelle logarithmique sur les axes y
	c1->SetGrid();
	c1->SetLogy();

	int number_datas = 19;

	TGraph *gr = new TGraphErrors(number_datas,Energy_abscisse,datas_300keV,0,0); // Création du graphique

	gr->SetMarkerColor(4); // Options de mise en forme  du graphique
	gr->SetMarkerStyle(26);
	gr->SetMarkerSize(0.8);

	gr->Draw("P");


	TGraph *gr1 = new TGraphErrors(number_datas,Energy_abscisse,datas_500kev,0,0); // Création du graphique

	gr1->SetMarkerColor(3); // Options de mise en forme  du graphique
	gr1->SetMarkerStyle(25);
	gr1->SetMarkerSize(0.8);
	gr1->Draw("P"); // Dessine le graph


	TGraph *gr2 = new TGraphErrors(number_datas,Energy_abscisse,datas_1Mev,0,0); // Création du graphique

	gr2->SetMarkerColor(2); // Options de mise en forme  du graphique
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerSize(0.8);
	gr2->Draw("P"); // Dessine le graph


	TGraph *gr3 = new TGraphErrors(number_datas,Energy_abscisse,datas_2Mev,0,0); // Création du graphique

	gr3->SetMarkerColor(6); // Options de mise en forme  du graphique
	gr3->SetMarkerStyle(22);
	gr3->SetMarkerSize(0.8);

	gr3->Draw("P");


	TGraph *gr4 = new TGraphErrors(number_datas,Energy_abscisse,datas_4Mev,0,0); // Création du graphique

	gr4->SetMarkerColor(1); // Options de mise en forme  du graphique
	gr4->SetMarkerStyle(26);
	gr4->SetMarkerSize(1);
	gr4->Draw("P"); // Dessine le graph


	TGraph *gr5 = new TGraphErrors(number_datas,Energy_abscisse,datas_6Mev,0,0); // Création du graphique

	gr5->SetMarkerColor(9); // Options de mise en forme  du graphique
	gr5->SetMarkerStyle(27);
	gr5->SetMarkerSize(1);
	gr5->Draw("P"); // Dessine le graph




	leg = new TLegend(0.3,0.75,0.85,0.89);
	leg->SetFillColor(10);


  /*  leg->AddEntry(gr,"300 keV - No cuts","p");
    leg->AddEntry(gr1,"300 keV - Cut Si 50 keV ","p");
    leg->AddEntry(gr2,"300 keV - Cut BGO 100 keV ","p");*/


    //No cuts energie Si 0keV et BGO 0 keV
	leg->AddEntry(gr,"300 keV","p");
	leg->AddEntry(gr1,"500 keV","p");
    leg->AddEntry(gr2,"1 MeV","p");
	leg->AddEntry(gr3,"2 MeV","p");
    leg->AddEntry(gr4,"4 MeV","p");
    leg->AddEntry(gr5,"6 MeV","p");
	//

	leg->Draw();
	c1->Update();


}
