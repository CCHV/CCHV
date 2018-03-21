
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
	Double_t datas_300keV_cutSingle [19]	= {8394, 9471, 8346, 6356, 5552, 5136, 4702, 4585, 4607, 4529, 4493, 4451, 4839, 5089, 5647, 6390, 8364, 9487, 8413} ; //300 keV
	Double_t datas_500kev_cutSingle [19]		={13333, 20368, 23288, 23482, 22281, 20490, 18807, 17929, 17416, 17313, 17683, 17619, 18924, 20526, 22196, 23520, 23108, 20372, 13431};// 500 kev
	Double_t datas_1Mev_cutSingle[19]	= {9445, 18301, 24931, 31245, 33388, 35296, 37258, 38002, 38008, 38384, 38232, 38021, 37205, 35570, 33471, 31061, 24808, 18343, 9598}   ; // 1 MeV
	Double_t datas_2Mev_cutSingle [19]	=  {6144, 12638, 18884, 27231, 30670, 33715, 36735, 38710, 39135, 39267, 39583, 38691, 37380, 34377, 30562, 27051, 18893, 12580, 6128}; // 2 MeV
	Double_t datas_4Mev_cutSingle [19]	=  {4097, 8120, 13049, 20169, 23438, 27023, 29471, 31723, 31670, 32395, 31868, 31404, 29947, 26794, 23748, 20208, 12929, 8082, 4072}; // 4 MeV
	Double_t datas_6Mev_cutSingle [19]	=  {3427, 6890, 11149, 17901, 20680, 23544, 26003, 27605, 28136, 28107, 28268, 27654, 25613, 23431, 20591, 17654, 10920, 6803, 3474}; // 6 MeV

    //Sans aucun Cut énergie : Si 0 keV et BGO 0 keV
  Double_t datas_300keV [19]	= { 15613, 29301, 39674, 52957, 58117, 62280, 66329, 68280, 69533, 69090, 68703, 68353, 65598, 62375, 57524, 52683, 40086, 29266, 15623} ; // 300 keV
  Double_t datas_500kev [19]		={13777, 26258, 37279, 50398, 55703, 60407, 63904, 66698, 67294, 67664, 67691, 66622, 64324, 60400, 55458, 50101, 37185, 26296, 13839}; // 500 kev
  Double_t datas_1Mev[19]	= {9870, 19827, 28999, 40547, 45297, 49288, 53738, 55867, 56374, 56727, 56481, 55632, 53474, 49697, 45368, 40480, 28756, 19876, 10061}   ; // 1 Mev
  Double_t datas_2Mev[19]	=  {6551, 13305, 20377, 30405, 34491, 38189, 41996, 44327, 44630, 45210, 45163, 44350, 42573, 39127, 34372, 30274, 20471, 13334, 6543}; // 2 MeV
  Double_t datas_4Mev[19]	=  {4341, 8456, 13460, 21169, 24686, 28286, 30743, 33250, 33215, 33868, 33612, 32961, 31398, 28147, 24991, 21265, 13377, 8388, 4314 }; // 4 MeV
  Double_t datas_6Mev[19]	=  {3659, 7165, 11594, 18370, 21365, 24191, 26696, 28389, 28988, 28932, 29103, 28478, 26399, 24184, 21298, 18262, 11309, 7146, 3744}; // 6 MeV


//Avec Cut sum energy Si + BGO 90% of tot
Double_t datas_300keV_cutSum [19]	= {12650, 23534, 31529, 41798, 45871, 49158, 52102, 53644, 54523, 54154, 53794, 53567, 51496, 49187, 45458, 41659, 31976, 23480, 12631} ; //300 keV
Double_t datas_500kev_cutSum [19]		={11813, 22112, 30732, 41455, 45486, 49254, 51971, 54185, 54566, 54797, 54846, 53963, 52211, 49056, 45229, 41061, 30972, 22152, 11968};// 500 kev
Double_t datas_1Mev_cutSum[19]	= {7859, 15302, 22130, 30136, 33334, 35925, 38809, 40412, 40555, 40924, 40882, 39994, 38861, 36433, 33227, 30146, 21914, 15469, 7908}   ; // 1 MeV
Double_t datas_2Mev_cutSum [19]	=  {3479, 7216, 11041, 16430, 18505, 20643, 22553, 23845, 24302, 24336, 24302, 23841, 22945, 21089, 18672, 16450, 11103, 7290, 3503}; // 2 MeV
Double_t datas_4Mev_cutSum [19]	=  {1106, 2482, 4525, 7823, 9410, 11176, 12375, 13524, 13579, 13866, 13577, 13326, 12564, 11204, 9621, 7753, 4497, 2502, 1061}; // 4 MeV
Double_t datas_6Mev_cutSum [19]	=  {958, 1861, 3418, 6378, 7540, 9086, 10262, 11034, 11274, 11287, 11403, 10989, 10236, 8968, 7634, 6398, 3266, 1799, 935}; // 6 MeV


//Avec Cut sum energy Si + BGO 90% of tot
/*Double_t datas_300keV_cutSumSingle [19]	= {} ; //300 keV
Double_t datas_500kev_cutSumSingle [19]		={};// 500 kev
Double_t datas_1Mev_cutSumSingle[19]	= {}   ; // 1 MeV
Double_t datas_2Mev_cutSumSingle [19]	=  {}; // 2 MeV
Double_t datas_4Mev_cutSumSingle [19]	=  {}; // 4 MeV
Double_t datas_6Mev_cutSumSingle [19]	=  {}; // 6 MeV
*/
Double_t ratio_noCut_cutSum[6][19];
Double_t ratio_cutSingle_cutSingleSum[6][19];


	for(Int_t rat = 0; rat <19; rat ++){

		ratio_noCut_cutSum[0][rat]=datas_300keV_cutSum[rat]/datas_300keV[rat];
		ratio_noCut_cutSum[1][rat]=datas_500kev_cutSum[rat]/datas_500kev[rat];
		ratio_noCut_cutSum[2][rat]=datas_1Mev_cutSum[rat]/datas_1Mev[rat];
		ratio_noCut_cutSum[3][rat]=datas_2Mev_cutSum[rat]/datas_2Mev[rat];
		ratio_noCut_cutSum[4][rat]=datas_4Mev_cutSum[rat]/datas_4Mev[rat];
		ratio_noCut_cutSum[5][rat]=datas_6Mev_cutSum[rat]/datas_6Mev[rat];

		/*ratio_cutSingle_cutSingleSum[0][rat]=datas_300keV_cutSumSingle[rat]/datas_300keV_cutSingle[rat];
		ratio_cutSingle_cutSingleSum[1][rat]=datas_500kev_cutSumSingle[rat]/datas_500kev_cutSingle[rat];
		ratio_cutSingle_cutSingleSum[2][rat]=datas_1Mev_cutSumSingle[rat]/datas_1Mev_cutSingle[rat];
		ratio_cutSingle_cutSingleSum[3][rat]=datas_2Mev_cutSumSingle[rat]/datas_2Mev_cutSingle[rat];
		ratio_cutSingle_cutSingleSum[4][rat]=datas_4Mev_cutSumSingle[rat]/datas_4Mev_cutSingle[rat];
		ratio_cutSingle_cutSingleSum[5][rat]=datas_6Mev_cutSumSingle[rat]/datas_6Mev_cutSingle[rat];
*/
	}

	std::ofstream out_txt("output_efficiency_ratio_sumCut.txt", std::ios::out);
	std::ofstream out_txt_1("output_efficiency_ratio_sumSingleCut.txt", std::ios::out);
	Double_t energies[6] = {0.3, 0.5,1., 2., 4., 6.};
	out_txt<<"\t"<<energies[0]<<"\t"<<energies[1]<<"\t"<<energies[2]<<"\t"<<energies[3]<<"\t"<<energies[4]<<"\t"<<energies[5]<<endl;
	out_txt_1<<"\t"<<energies[0]<<"\t"<<energies[1]<<"\t"<<energies[2]<<"\t"<<energies[3]<<"\t"<<energies[4]<<"\t"<<energies[5]<<endl;

		for(Int_t pos = 0; pos<19; pos++){
			for(Int_t en = 0; en<6; en ++){
				if(en==0){out_txt<<Energy_abscisse[pos]; out_txt_1<<Energy_abscisse[pos];}
				out_txt<<"\t"<<std::setprecision(3)<<ratio_noCut_cutSum[en][pos];
				//out_txt_1<<"\t"<<ratio_cutSingle_cutSingleSum[en][pos];
				if(en==5){out_txt<<std::endl; out_txt_1<<std::endl;}
			}
	}
out_txt.close();
out_txt_1.close();

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
				datas_300keV_cutSum[i] = datas_300keV_cutSum[i]*(1/(9*pow(10,7)*2)); // simulation sur 2 pi et 9 simu de 1e7
        datas_500kev_cutSum[i] = datas_500kev_cutSum[i]*(1/(9*pow(10,7)*2));
        datas_1Mev_cutSum[i] = datas_1Mev_cutSum[i]*(1/(9*pow(10,7)*2));
        datas_2Mev_cutSum[i] = datas_2Mev_cutSum[i]*(1/(9*pow(10,7)*2));
        datas_4Mev_cutSum[i] = datas_4Mev_cutSum[i]*(1/(9*pow(10,7)*2));
        datas_6Mev_cutSum[i] = datas_6Mev_cutSum[i]*(1/(9*pow(10,7)*2));

				datas_300keV_cutSingle[i] = datas_300keV_cutSingle[i]*(1/(9*pow(10,7)*2)); // simulation sur 2 pi et 9 simu de 1e7
        datas_500kev_cutSingle[i] = datas_500kev_cutSingle[i]*(1/(9*pow(10,7)*2));
        datas_1Mev_cutSingle[i] = datas_1Mev_cutSingle[i]*(1/(9*pow(10,7)*2));
        datas_2Mev_cutSingle[i] = datas_2Mev_cutSingle[i]*(1/(9*pow(10,7)*2));
        datas_4Mev_cutSingle[i] = datas_4Mev_cutSingle[i]*(1/(9*pow(10,7)*2));
        datas_6Mev_cutSingle[i] = datas_6Mev_cutSingle[i]*(1/(9*pow(10,7)*2));
    }



	char name[80];	sprintf(name, "Detection efficiency vs source position"); // Nom du graphique
	TLegend *leg;
	TLegend *leg_cut;

	TH1 *frame = new TH1F("frame","",1000,-300,300);
	frame->SetMinimum(1e-6);
	frame->SetMaximum(5e-4);
	frame->SetDirectory(0);
	frame->SetStats(0);
	frame->SetTitle(name);
	frame->GetXaxis()->SetTitle("Source position [mm]");
	frame->GetXaxis()->SetTickLength(0.02);
	frame->GetXaxis()->SetLabelSize(0.03);
	frame->GetYaxis()->SetTitle("Efficiency");

	frame->GetYaxis()->SetLabelSize(0.03);
	frame->Draw(" ");


	 // Échelle logarithmique sur les axes y
	c1->SetGrid();
	c1->SetLogy();

	int number_datas = 19;


	//No cut
	TGraph *gr = new TGraphErrors(number_datas,Energy_abscisse,datas_300keV,0,0); // Création du graphique
	gr->SetMarkerColor(4); // Options de mise en forme  du graphique
	gr->SetMarkerStyle(22);
	gr->SetMarkerSize(1);
	//gr->Draw("P");

	TGraph *gr1 = new TGraphErrors(number_datas,Energy_abscisse,datas_500kev,0,0); // Création du graphique
	gr1->SetMarkerColor(3); // Options de mise en forme  du graphique
	gr1->SetMarkerStyle(22);
	gr1->SetMarkerSize(1);
	//gr1->Draw("P"); // Dessine le graph

	TGraph *gr2 = new TGraphErrors(number_datas,Energy_abscisse,datas_1Mev,0,0); // Création du graphique
	gr2->SetMarkerColor(2); // Options de mise en forme  du graphique
	gr2->SetMarkerStyle(22);
	gr2->SetMarkerSize(1);
	//gr2->Draw("P"); // Dessine le graph

	TGraph *gr3 = new TGraphErrors(number_datas,Energy_abscisse,datas_2Mev,0,0); // Création du graphique
	gr3->SetMarkerColor(6); // Options de mise en forme  du graphique
	gr3->SetMarkerStyle(22);
	gr3->SetMarkerSize(1);
	//gr3->Draw("P");

	TGraph *gr4 = new TGraphErrors(number_datas,Energy_abscisse,datas_4Mev,0,0); // Création du graphique
	gr4->SetMarkerColor(1); // Options de mise en forme  du graphique
	gr4->SetMarkerStyle(22);
	gr4->SetMarkerSize(1);
	//gr4->Draw("P"); // Dessine le graph

	TGraph *gr5 = new TGraphErrors(number_datas,Energy_abscisse,datas_6Mev,0,0); // Création du graphique
	gr5->SetMarkerColor(kOrange+1); // Options de mise en forme  du graphique
	gr5->SetMarkerStyle(22);
	gr5->SetMarkerSize(1);
	//gr5->Draw("P"); // Dessine le graph

	//Cut on single layers
	TGraph *gr_cutSingle = new TGraphErrors(number_datas,Energy_abscisse,datas_300keV_cutSingle,0,0); // Création du graphique
	gr_cutSingle->SetMarkerColor(4); // Options de mise en forme  du graphique
	gr_cutSingle->SetMarkerStyle(22);
	gr_cutSingle->SetMarkerSize(1);
	//gr_cutSingle->Draw("P");

	TGraph *gr1_cutSingle = new TGraphErrors(number_datas,Energy_abscisse,datas_500kev_cutSingle,0,0); // Création du graphique
	gr1_cutSingle->SetMarkerColor(3); // Options de mise en forme  du graphique
	gr1_cutSingle->SetMarkerStyle(22);
	gr1_cutSingle->SetMarkerSize(1);
	//gr1_cutSingle->Draw("P"); // Dessine le graph

	TGraph *gr2_cutSingle = new TGraphErrors(number_datas,Energy_abscisse,datas_1Mev_cutSingle,0,0); // Création du graphique
	gr2_cutSingle->SetMarkerColor(2); // Options de mise en forme  du graphique
	gr2_cutSingle->SetMarkerStyle(22);
	gr2_cutSingle->SetMarkerSize(1);
	//gr2_cutSingle->Draw("P"); // Dessine le graph

	TGraph *gr3_cutSingle = new TGraphErrors(number_datas,Energy_abscisse,datas_2Mev_cutSingle,0,0); // Création du graphique
	gr3_cutSingle->SetMarkerColor(6); // Options de mise en forme  du graphique
	gr3_cutSingle->SetMarkerStyle(22);
	gr3_cutSingle->SetMarkerSize(1);
	//gr3_cutSingle->Draw("P");

	TGraph *gr4_cutSingle = new TGraphErrors(number_datas,Energy_abscisse,datas_4Mev_cutSingle,0,0); // Création du graphique
	gr4_cutSingle->SetMarkerColor(1); // Options de mise en forme  du graphique
	gr4_cutSingle->SetMarkerStyle(22);
	gr4_cutSingle->SetMarkerSize(1);
	//gr4_cutSingle->Draw("P"); // Dessine le graph

	TGraph *gr5_cutSingle = new TGraphErrors(number_datas,Energy_abscisse,datas_6Mev_cutSingle,0,0); // Création du graphique
	gr5_cutSingle->SetMarkerColor(kOrange+1); // Options de mise en forme  du graphique
	gr5_cutSingle->SetMarkerStyle(22);
	gr5_cutSingle->SetMarkerSize(1);
	//gr5_cutSingle->Draw("P"); // Dessine le graph

	//Cut BGO+Si energy deposit
	TGraph *gr_cutSum = new TGraphErrors(number_datas,Energy_abscisse,datas_300keV_cutSum,0,0); // Création du graphique
	gr_cutSum->SetMarkerColor(4); // Options de mise en forme  du graphique
	gr_cutSum->SetMarkerStyle(28);
	gr_cutSum->SetMarkerSize(1);
	//gr_cutSum->Draw("P");

	TGraph *gr1_cutSum = new TGraphErrors(number_datas,Energy_abscisse,datas_500kev_cutSum,0,0); // Création du graphique
	gr1_cutSum->SetMarkerColor(3); // Options de mise en forme  du graphique
	gr1_cutSum->SetMarkerStyle(28);
	gr1_cutSum->SetMarkerSize(1);
	//gr1_cutSum->Draw("P"); // Dessine le graph

	TGraph *gr2_cutSum = new TGraphErrors(number_datas,Energy_abscisse,datas_1Mev_cutSum,0,0); // Création du graphique
	gr2_cutSum->SetMarkerColor(2); // Options de mise en forme  du graphique
	gr2_cutSum->SetMarkerStyle(28);
	gr2_cutSum->SetMarkerSize(1);
	//gr2_cutSum->Draw("P"); // Dessine le graph

	TGraph *gr3_cutSum = new TGraphErrors(number_datas,Energy_abscisse,datas_2Mev_cutSum,0,0); // Création du graphique
	gr3_cutSum->SetMarkerColor(6); // Options de mise en forme  du graphique
	gr3_cutSum->SetMarkerStyle(28);
	gr3_cutSum->SetMarkerSize(1);
	//gr3_cutSum->Draw("P");

	TGraph *gr4_cutSum = new TGraphErrors(number_datas,Energy_abscisse,datas_4Mev_cutSum,0,0); // Création du graphique
	gr4_cutSum->SetMarkerColor(1); // Options de mise en forme  du graphique
	gr4_cutSum->SetMarkerStyle(28);
	gr4_cutSum->SetMarkerSize(1);
	//gr4_cutSum->Draw("P"); // Dessine le graph

	TGraph *gr5_cutSum = new TGraphErrors(number_datas,Energy_abscisse,datas_6Mev_cutSum,0,0); // Création du graphique
	gr5_cutSum->SetMarkerColor(kOrange+1); // Options de mise en forme  du graphique
	gr5_cutSum->SetMarkerStyle(28);
	gr5_cutSum->SetMarkerSize(1);
	//gr5_cutSum->Draw("P"); // Dessine le graph

	//Cut single layer and sum BGO+Si energy deposit
	/*TGraph *gr_cutSingleSum = new TGraphErrors(number_datas,Energy_abscisse,datas_300keV_cutSingleSum,0,0); // Création du graphique
	gr_cutSingleSum->SetMarkerColor(4); // Options de mise en forme  du graphique
	gr_cutSingleSum->SetMarkerStyle(26);
	gr_cutSingleSum->SetMarkerSize(1);
	//gr_cutSingleSum->Draw("P");

	TGraph *gr1_cutSingleSum = new TGraphErrors(number_datas,Energy_abscisse,datas_500kev_cutSingleSum,0,0); // Création du graphique
	gr1_cutSingleSum->SetMarkerColor(3); // Options de mise en forme  du graphique
	gr1_cutSingleSum->SetMarkerStyle(26);
	gr1_cutSingleSum->SetMarkerSize(1);
	//gr1_cutSingleSum->Draw("P"); // Dessine le graph

	TGraph *gr2_cutSingleSum = new TGraphErrors(number_datas,Energy_abscisse,datas_1Mev_cutSingleSum,0,0); // Création du graphique
	gr2_cutSingleSum->SetMarkerColor(2); // Options de mise en forme  du graphique
	gr2_cutSingleSum->SetMarkerStyle(26);
	gr2_cutSingleSum->SetMarkerSize(1);
	//gr2_cutSingleSum->Draw("P"); // Dessine le graph

	TGraph *gr3_cutSingleSum = new TGraphErrors(number_datas,Energy_abscisse,datas_2Mev_cutSingleSum,0,0); // Création du graphique
	gr3_cutSingleSum->SetMarkerColor(6); // Options de mise en forme  du graphique
	gr3_cutSingleSum->SetMarkerStyle(26);
	gr3_cutSingleSum->SetMarkerSize(1);
	//gr3_cutSingleSum->Draw("P");

	TGraph *gr4_cutSingleSum = new TGraphErrors(number_datas,Energy_abscisse,datas_4Mev_cutSingleSum,0,0); // Création du graphique
	gr4_cutSingleSum->SetMarkerColor(1); // Options de mise en forme  du graphique
	gr4_cutSingleSum->SetMarkerStyle(26);
	gr4_cutSingleSum->SetMarkerSize(1);
	//gr4_cutSingleSum->Draw("P"); // Dessine le graph

	TGraph *gr5_cutSingleSum = new TGraphErrors(number_datas,Energy_abscisse,datas_6Mev_cutSingleSum,0,0); // Création du graphique
	gr5_cutSingleSum->SetMarkerColor(kOrange+1); // Options de mise en forme  du graphique
	gr5_cutSingleSum->SetMarkerStyle(26);
	gr5_cutSingleSum->SetMarkerSize(1);
	//gr5_cutSingleSum->Draw("P"); // Dessine le graph
*/

	leg = new TLegend(0.3,0.12,0.5,0.45);
	//leg->SetHeader("#splitline{50 keV threshold on scatterer}{100 keV threshold on absorber}"); // option "C" allows to center the header
	leg->SetHeader("No energy cut");
	leg->SetFillColor(10);
	leg->SetTextSize(0.03);
	leg_cut = new TLegend(0.5,0.12,0.7,0.45);
	leg_cut->SetHeader("> 90 \% E dep"); // option "C" allows to center the header
	leg_cut->SetFillColor(10);
	leg_cut->SetTextSize(0.03);

  /*  leg->AddEntry(gr,"300 keV - No cuts","p");
    leg->AddEntry(gr1,"300 keV - Cut Si 50 keV ","p");
    leg->AddEntry(gr2,"300 keV - Cut BGO 100 keV ","p");*/

		//Choose here what to plot
		//no cut
		c1->cd();
		gr->Draw("P");
		gr1->Draw("P");
		gr2->Draw("P");
		gr3->Draw("P");
		gr4->Draw("P");
		gr5->Draw("P");
		//single layer cuts
		/*gr_cutSingle->Draw("P");
		gr1_cutSingle->Draw("P");
		gr2_cutSingle->Draw("P");
		gr3_cutSingle->Draw("P");
		gr4_cutSingle->Draw("P");
		gr5_cutSingle->Draw("P");*/
		//sum of energy deposit cut
		c1->cd();
		gr_cutSum->Draw("P");
		gr1_cutSum->Draw("P");
		gr2_cutSum->Draw("P");
		gr3_cutSum->Draw("P");
		gr4_cutSum->Draw("P");
		gr5_cutSum->Draw("P");
		//single layer and sum cut
		/*gr_cutSingleSum->Draw("P");
		gr1_cutSingleSum->Draw("P");
		gr2_cutSingleSum->Draw("P");
		gr3_cutSingleSum->Draw("P");
		gr4_cutSingleSum->Draw("P");
		gr5_cutSingleSum->Draw("P");*/


		//No cuts energie Si 0keV et BGO 0 keV
	leg->AddEntry(gr,"300 keV","p");
	leg->AddEntry(gr1,"500 keV","p");
  leg->AddEntry(gr2,"1 MeV","p");
	leg->AddEntry(gr3,"2 MeV","p");
  leg->AddEntry(gr4,"4 MeV","p");
  leg->AddEntry(gr5,"6 MeV","p");
	leg_cut->AddEntry(gr_cutSum,"300 keV","p");
	leg_cut->AddEntry(gr1_cutSum,"500 keV","p");
	leg_cut->AddEntry(gr2_cutSum,"1 MeV","p");
	leg_cut->AddEntry(gr3_cutSum,"2 MeV","p");
  leg_cut->AddEntry(gr4_cutSum,"4 MeV","p");
  leg_cut->AddEntry(gr5_cutSum,"6 MeV","p");
	//

	leg->Draw();
	leg_cut->Draw();
	c1->Update();


}
