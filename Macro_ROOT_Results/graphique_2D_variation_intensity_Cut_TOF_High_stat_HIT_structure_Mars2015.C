
// Mise à jour le 18 aout 2014 JLL

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



//-------------------------------------------------------------------------------------------
//Datas avec changement de la structure du faisceau HIT : Mesure faite lors du week end de Paques 4 avril 2015
//-------------------------------------------------------------------------------------------


void graphique_2D_variation_intensity_Cut_TOF_High_stat_HIT_structure_Mars2015()
{

TCanvas *c1 = new TCanvas("c1","c1",700,800); // Création du canvas
c1->SetLeftMargin(0.12);
    const int b =23;//23 protons , 23 carbon // number of values

	const int b_sans =23; // sans cut TOF

    Int_t choix = 4; // choix = 1 : %;  choix = 2 : absolu ; choix = 3 : % PTCOG; choix = 4 : absolu PTCOG

	Int_t choix_p_c = 1; // choix_p_c = 1 : protons, choix_p_c = 2 : carbone

	Int_t choix_histo = 1; //choisit de mettre l'histo sur la courbe ou pas: choix_histo =0: pas d'histo, choix_histo =1: pas d'histo,

	Int_t comparaison = 1; // si on veut comparer avec ou sans cut TOF


//-----------------	data 05/01/2015 - Prise en compte fenêtre coïnc [-20;20] ns et CUT TOF [0;6] ns + distribution Poisson ------------

	//Protons

//Double_t sigma_ns [b]	 = {0.001,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,80,100,150,200,217};//b=23

	//Carbon
	Double_t sigma_ns [b]	 = {0.0001,0.001,0.005,0.01,0.05,0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,60,70};//b =23


	//-----------------	data data 05/01/2015 - Prise en compte fenêtre coïnc [-20;20] ns et NO cut TOF + distribution Poisson ------------

	//Protons
	//Double_t sigma_ns_sans_cut [b_sans]	 = {0.001,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,80,100,150,200,217};//b=23

	//Carbon
	Double_t sigma_ns_sans_cut [b_sans]	 = {0.0001,0.001,0.005,0.01,0.05,0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,60,70};//b =23





	//Définition variables signal to background ratio SANS cut TOF
	Double_t Ratio_Background_Signal_proton_absolu_sans_cut [b_sans];
	Double_t Ratio_Background_Signal_proton_relatif_sans_cut[b_sans];
	Double_t Ratio_Background_Signal_proton_histo_sans_cut[b_sans];

	Double_t Ratio_Background_Signal_carbon_absolu_sans_cut  [b_sans];
	Double_t Ratio_Background_Signal_carbon_relatif_sans_cut [b_sans];
	Double_t Ratio_Background_Signal_carbon_histo_sans_cut [b_sans];


	//Définition variables signal to background ratio AVEC cut TOF
	Double_t Ratio_Background_Signal_proton_absolu[b];
	Double_t Ratio_Background_Signal_proton_relatif[b];
	Double_t Ratio_Background_Signal_proton_histo[b];

	Double_t Ratio_Background_Signal_carbon_absolu  [b];
	Double_t Ratio_Background_Signal_carbon_relatif [b];
	Double_t Ratio_Background_Signal_carbon_histo [b];

	//Définition variables signal to background ratio comparaison AVEC et SANS cut TOF

	Double_t Ratio_Background_Signal_proton_absolu_comparaison[b];
	Double_t Ratio_Background_Signal_proton_relatif_comparaison[b];
	Double_t Ratio_Background_Signal_proton_histo_comparaison[b];

	Double_t Ratio_Background_Signal_carbon_absolu_comparaison [b];
	Double_t Ratio_Background_Signal_carbon_relatif_comparaison [b];
	Double_t Ratio_Background_Signal_carbon_histo_comparaison [b];


	Double_t Ratio_Background_Signal_proton_absolu_comparaison_gamma[b];
	Double_t Ratio_Background_Signal_proton_relatif_comparaison_gamma[b];
	Double_t Ratio_Background_Signal_proton_histo_comparaison_gamma[b];

	Double_t Ratio_Background_Signal_proton_absolu_comparaison_bad[b];
	Double_t Ratio_Background_Signal_proton_relatif_comparaison_bad[b];
	Double_t Ratio_Background_Signal_proton_histo_comparaison_bad[b];

	Double_t Ratio_Background_Signal_carbon_absolu_comparaison_gamma  [b];
	Double_t Ratio_Background_Signal_carbon_relatif_comparaison_gamma [b];
	Double_t Ratio_Background_Signal_carbon_histo_comparaison_gamma [b];

	Double_t Ratio_Background_Signal_carbon_absolu_comparaison_bad  [b];
	Double_t Ratio_Background_Signal_carbon_relatif_comparaison_bad [b];
	Double_t Ratio_Background_Signal_carbon_histo_comparaison_bad [b];



	for(int i =0;i<b_sans;i++)
	{
		Ratio_Background_Signal_proton_absolu_sans_cut [i]=0;
		Ratio_Background_Signal_proton_relatif_sans_cut[i]=0;
		Ratio_Background_Signal_proton_histo_sans_cut[i]=0;

		Ratio_Background_Signal_carbon_absolu_sans_cut  [i]=0;
		Ratio_Background_Signal_carbon_relatif_sans_cut [i]=0;
		Ratio_Background_Signal_carbon_histo_sans_cut [i]=0;
	}

	for(int i =0;i<b;i++)
	{

		Ratio_Background_Signal_proton_absolu [i]=0;
		Ratio_Background_Signal_proton_relatif[i]=0;
		Ratio_Background_Signal_proton_histo[i]=0;

		Ratio_Background_Signal_carbon_absolu  [i]=0;
		Ratio_Background_Signal_carbon_relatif [i]=0;
		Ratio_Background_Signal_carbon_histo [i]=0;

		//protons
		Ratio_Background_Signal_proton_absolu_comparaison [i]=0; // Rapport des ratio S/B avec et sans cut
		Ratio_Background_Signal_proton_relatif_comparaison[i]=0;
		Ratio_Background_Signal_proton_histo_comparaison[i]=0;

		Ratio_Background_Signal_proton_absolu_comparaison_gamma [i]=0; // Rapport des taux avec et sans cut pour les bons events gamma
		Ratio_Background_Signal_proton_relatif_comparaison_gamma[i]=0;
		Ratio_Background_Signal_proton_histo_comparaison_gamma[i]=0;

		Ratio_Background_Signal_proton_absolu_comparaison_bad  [i]=0; // Rapport des taux avec et sans cut pour les bad events - background
		Ratio_Background_Signal_proton_relatif_comparaison_bad [i]=0;
		Ratio_Background_Signal_proton_histo_comparaison_bad [i]=0;

		//carbon ions
		Ratio_Background_Signal_carbon_absolu_comparaison [i]=0; // Rapport des ratio S/B avec et sans cut
		Ratio_Background_Signal_carbon_relatif_comparaison[i]=0;
		Ratio_Background_Signal_carbon_histo_comparaison[i]=0;

		Ratio_Background_Signal_carbon_absolu_comparaison_gamma [i]=0; // Rapport des taux avec et sans cut pour les bons events gamma
		Ratio_Background_Signal_carbon_relatif_comparaison_gamma[i]=0;
		Ratio_Background_Signal_carbon_histo_comparaison_gamma[i]=0;

		Ratio_Background_Signal_carbon_absolu_comparaison_bad  [i]=0; // Rapport des taux avec et sans cut pour les bad events - background
		Ratio_Background_Signal_carbon_relatif_comparaison_bad [i]=0;
		Ratio_Background_Signal_carbon_histo_comparaison_bad [i]=0;


	}

	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- PROTONS -------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------




	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- Avec CUT TOF --------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------



	//if (choix_p_c ==1)
	//{


	//-------------------
	//-------------------
	//Coinc Pourcentage
	//-------------------
	//-------------------



	//Variation intensity
	Double_t datas_BGO_pourcentage_Fortuitous[b]={0.0458295,0,0.275862,1.14574,2.48258,4.36331,6.22585,8.7926,10.3577,12.3427,12.8862,14.7882,16.8607,17.537,24.8183,29.5368,38.248,50.6064,61.7862,66.7969,74.8592,79.6517,80.992};
	Double_t datas_BGO_pourcentage_True_gamma[b]	={99.9542,100,99.7241,98.8543,97.5174,95.6367,93.7742,91.2074,89.6423,87.6573,87.1138,85.2118,83.1393,82.463,75.1817,70.4632,61.752,49.3936,38.2138,33.2031,25.1408,20.3483,19.008};
	Double_t datas_BGO_pourcentage_True_Other_gamma[b]= {33.9597,33.8622,33.2414,33.3639,32.7962,31.7008,31.6445,30.3324,29.7005,29.8417,30.1626,28.0576,28.0511,27.2798,25.5106,24.4289,21.2778,17.503,12.5366,11.1057,7.86061,6.73322,6.01019};
	Double_t datas_BGO_pourcentage_Same_gamma[b]	={56.4161,56.4823,56.2299,55.7745,55.0087,53.6955,53.3705,51.83,50.8319,48.8023,48.4553,48.721,46.9771,47.1941,42.3676,39.2449,35.0101,27.1819,21.9436,18.922,14.7124,11.6929,11.1573};
	Double_t datas_BGO_pourcentage_Same_Particle[b]	={9.57837,9.65549,10.2529,9.71586,9.71254,10.2404,8.75912,9.04501,9.10982,9.0134,8.49593,8.43325,8.11115,7.98909,7.30357,6.78934,5.46401,4.70868,3.73353,3.17538,2.5678,1.92226,1.84049};

		Double_t Sum_bad_events_pourcentage[b];

		for(int u=0;u<b;u++)
		{
			Sum_bad_events_pourcentage[u]= datas_BGO_pourcentage_Fortuitous[u]+datas_BGO_pourcentage_True_gamma[u]-datas_BGO_pourcentage_Same_gamma[u];
		}

	//--------------------------------
	//Reconstructed events (line-cone)
	//--------------------------------


	Double_t datas_BGO_pourcentage_total_histo[b]={1438,1446,1438,1418,1497,1451,1489,1510,1537,1571,1585,1586,1639,1586,1739,1879,2048,2362,3012,3266,3977,4798,4934};
	Double_t datas_BGO_pourcentage_Fortuitous_histo[b]={0,0,3,12,27,39,75,98,116,150,156,177,215,211,338,429,639,999,1642,1921,2736,3573,3748} ;
	Double_t datas_BGO_pourcentage_True_gamma_histo[b]	={1438,1446,1435,1406,1470,1412,1414,1412,1421,1421,1429,1409,1424,1375,1401,1450,1409,1363,1370,1345,1241,1225,1186};
	Double_t datas_BGO_pourcentage_True_Other_gamma_histo[b]= {466,471,461,463,476,450,455,445,454,479,480,453,473,437,450,492,464,467,437,436,374,379,360};
	Double_t datas_BGO_pourcentage_Same_gamma_histo[b]	={843,855,847,814,870,825,847,841,845,815,822,832,834,820,825,845,836,785,818,796,754,745,732};
	Double_t datas_BGO_pourcentage_Same_Particle_histo[b]	={129,120,127,129,124,137,112,126,122,127,127,124,117,118,126,113,109,111,115,113,113,101,94};


	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_pourcentage_Fortuitous_histo[ff] = datas_BGO_pourcentage_Fortuitous_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_True_gamma_histo[ff]= datas_BGO_pourcentage_True_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_True_Other_gamma_histo[ff]= datas_BGO_pourcentage_True_Other_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_Same_gamma_histo[ff] = datas_BGO_pourcentage_Same_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_Same_Particle_histo[ff] = datas_BGO_pourcentage_Same_Particle_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;

	}



		Double_t Sum_bad_events_pourcentage_histo[b];

		for(int u=0;u<b;u++)
		{
			Sum_bad_events_pourcentage_histo[u]= datas_BGO_pourcentage_Fortuitous_histo[u]+datas_BGO_pourcentage_True_gamma_histo[u]-datas_BGO_pourcentage_Same_gamma_histo[u];
		}


	//-------------------
	//-------------------
	//Coinc Absolu
	//-------------------
	//-------------------

	//Variation intensity
	Double_t datas_BGO_absolu_Fortuitous[b]={1e-08,0,6e-08,2.5e-07,5.7e-07,9.8e-07,1.45e-06,2.09e-06,2.49e-06,3.04e-06,3.17e-06,3.7e-06,4.49e-06,4.5e-06,7.17e-06,9.31e-06,1.323e-05,2.128e-05,3.376e-05,4.102e-05,5.714e-05,7.5e-05,7.789e-05};
	Double_t datas_BGO_absolu_True_gamma[b]	={2.181e-05,2.206e-05,2.169e-05,2.157e-05,2.239e-05,2.148e-05,2.184e-05,2.168e-05,2.155e-05,2.159e-05,2.143e-05,2.132e-05,2.214e-05,2.116e-05,2.172e-05,2.221e-05,2.136e-05,2.077e-05,2.088e-05,2.039e-05,1.919e-05,1.916e-05,1.828e-05};
	Double_t datas_BGO_absolu_True_Other_gamma[b]= {7.41e-06,7.47e-06,7.23e-06,7.28e-06,7.53e-06,7.12e-06,7.37e-06,7.21e-06,7.14e-06,7.35e-06,7.42e-06,7.02e-06,7.47e-06,7e-06,7.37e-06,7.7e-06,7.36e-06,7.36e-06,6.85e-06,6.82e-06,6e-06,6.34e-06,5.78e-06};
	Double_t datas_BGO_absolu_Same_gamma[b]	={1.231e-05,1.246e-05,1.223e-05,1.217e-05,1.263e-05,1.206e-05,1.243e-05,1.232e-05,1.222e-05,1.202e-05,1.192e-05,1.219e-05,1.251e-05,1.211e-05,1.224e-05,1.237e-05,1.211e-05,1.143e-05,1.199e-05,1.162e-05,1.123e-05,1.101e-05,1.073e-05};
	Double_t datas_BGO_absolu_Same_Particle[b]	={2.09e-06,2.13e-06,2.23e-06,2.12e-06,2.23e-06,2.3e-06,2.04e-06,2.15e-06,2.19e-06,2.22e-06,2.09e-06,2.11e-06,2.16e-06,2.05e-06,2.11e-06,2.14e-06,1.89e-06,1.98e-06,2.04e-06,1.95e-06,1.96e-06,1.81e-06,1.77e-06};
file://localhost/Users/JeanLuc/Desktop/The%CC%80se/Article/Article/Clinic_Applicability_Compton_Camera/Analysis/Simulations/Timing_resolution_and_window_time/graphique_2D_variation_intensity_Cut_TOF_Distribution_poisson_Seminar2ndYear.C
		Double_t Sum_bad_events_absolu[b];

		for(int u=0;u<b;u++)
		{
			Sum_bad_events_absolu[u]= datas_BGO_absolu_Fortuitous[u]+datas_BGO_absolu_True_gamma[u]-datas_BGO_absolu_Same_gamma[u];
		}



	//--------------------------------
	//Reconstructed events (line-cone)
	//--------------------------------


	Double_t datas_BGO_absolu_total_histo[b]= {1438,1446,1438,1418,1497,1451,1489,1510,1537,1571,1585,1586,1639,1586,1739,1879,2048,2362,3012,3266,3977,4798,4934};
	Double_t datas_BGO_absolu_Fortuitous_histo[b]= {0,0,3,12,27,39,75,98,116,150,156,177,215,211,338,429,639,999,1642,1921,2736,3573,3748} ;
	Double_t datas_BGO_absolu_True_gamma_histo[b]= {1438,1446,1435,1406,1470,1412,1414,1412,1421,1421,1429,1409,1424,1375,1401,1450,1409,1363,1370,1345,1241,1225,1186};
	Double_t datas_BGO_absolu_True_Other_gamma_histo[b]= {466,471,461,463,476,450,455,445,454,479,480,453,473,437,450,492,464,467,437,436,374,379,360};
	Double_t datas_BGO_absolu_Same_gamma_histo[b]= {843,855,847,814,870,825,847,841,845,815,822,832,834,820,825,845,836,785,818,796,754,745,732};
	Double_t datas_BGO_absolu_Same_Particle_histo[b]= {129,120,127,129,124,137,112,126,122,127,127,124,117,118,126,113,109,111,115,113,113,101,94};


/*
	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_absolu_Fortuitous[ff] = datas_BGO_absolu_Fortuitous[ff]/ (100);
		datas_BGO_absolu_True_gamma[ff]= datas_BGO_absolu_True_gamma[ff]/ (100);
		datas_BGO_absolu_True_Other_gamma[ff]= datas_BGO_absolu_True_Other_gamma[ff]/ (100);
		datas_BGO_absolu_Same_gamma[ff] = datas_BGO_absolu_Same_gamma[ff]/ (100);
		datas_BGO_absolu_Same_Particle[ff] = datas_BGO_absolu_Same_Particle[ff]/ (100);

	}
	*/



	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_absolu_Fortuitous_histo[ff] = datas_BGO_absolu_Fortuitous_histo[ff]/ (1e8);
		datas_BGO_absolu_True_gamma_histo[ff]= datas_BGO_absolu_True_gamma_histo[ff]/ (1e8);
		datas_BGO_absolu_True_Other_gamma_histo[ff]= datas_BGO_absolu_True_Other_gamma_histo[ff]/ (1e8);
		datas_BGO_absolu_Same_gamma_histo[ff] = datas_BGO_absolu_Same_gamma_histo[ff]/ (1e8);
		datas_BGO_absolu_Same_Particle_histo[ff] = datas_BGO_absolu_Same_Particle_histo[ff]/ (1e8);

	}


		Double_t Sum_bad_events_absolu_histo[b];

		for(int u=0;u<b;u++)
		{
			Sum_bad_events_absolu_histo[u]= datas_BGO_absolu_Fortuitous_histo[u]+datas_BGO_absolu_True_gamma_histo[u]-datas_BGO_absolu_Same_gamma_histo[u];
		}



	//}






	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- SANS CUT TOF --------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------




	//if (choix_p_c ==1)
	//{


		//-------------------
		//-------------------
		//Coinc Pourcentage
		//-------------------
		//-------------------



		//Variation intensity
		Double_t datas_BGO_pourcentage_Fortuitous_sans_cut[b_sans]={0.0679117,0,0.27416,1.25042,2.49918,4.43575,5.62419,7.83307,9.64467,10.7352,12.112,13.0137,15.5247,16.2547,23.0169,27.9586,37.1479,48.8351,60.4767,65.681,74.1818,78.8033,80.0333};
		Double_t datas_BGO_pourcentage_True_gamma_sans_cut[b_sans]	={99.9321,100,99.7258,98.7496,97.5008,95.5643,94.3758,92.1669,90.3553,89.2648,87.888,86.9863,84.4753,83.7453,76.9831,72.0414,62.8521,51.1649,39.5233,34.319,25.8182,21.1967,19.9667};
		Double_t datas_BGO_pourcentage_True_Other_gamma_sans_cut[b_sans]= {29.1341,28.4381,29.4037,28.4218,28.1158,27.1363,27.8934,27.1589,25.698,25.3153,25.7152,25.4616,24.4001,24.3094,23.0976,21.8017,18.6738,14.3611,11.2332,9.88223,7.32597,6.11069,5.59767};
		Double_t datas_BGO_pourcentage_Same_gamma_sans_cut[b_sans]	={46.3158,46.8556,46.3331,46.4008,45.3798,44.7489,43.4655,42.1188,42.2906,41.7718,41.053,40.1132,39.4334,38.5868,35.4934,33.1567,29.3857,24.3356,18.4984,16.0522,12.2707,10.0347,9.53769};
		Double_t datas_BGO_pourcentage_Same_Particle_sans_cut[b_sans]	={24.4822,24.7063,23.989,23.927,24.0053,23.6791,23.0169,22.8892,22.3668,22.1778,21.1199,21.4116,20.6418,20.8491,18.392,17.083,14.7926,12.4681,9.79161,8.38454,6.2215,5.05128,4.83132};

		Double_t Sum_bad_events_pourcentage_sans_cut[b_sans];

		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_pourcentage_sans_cut[u]= datas_BGO_pourcentage_Fortuitous_sans_cut[u]+datas_BGO_pourcentage_True_gamma_sans_cut[u]-datas_BGO_pourcentage_Same_gamma_sans_cut[u];
		}

		//--------------------------------
		//Reconstructed events (line-cone)
		//--------------------------------


		Double_t datas_BGO_pourcentage_total_histo_sans_cut[b_sans]={1856,1823,1828,1870,1934,1915,1917,1926,1945,1998,2026,2041,2093,2060,2216,2282,2576,2970,3682,4050,5042,5857,6008};
		Double_t datas_BGO_pourcentage_Fortuitous_histo_sans_cut[b_sans]={0,0,4,18,33,60,87,114,137,152,187,202,240,252,413,487,787,1206,1935,2367,3414,4298,4481} ;
		Double_t datas_BGO_pourcentage_True_gamma_histo_sans_cut[b_sans]	={1856,1823,1824,1852,1901,1855,1830,1812,1808,1846,1839,1839,1853,1808,1803,1795,1789,1764,1747,1683,1628,1559,1527};
		Double_t datas_BGO_pourcentage_True_Other_gamma_histo_sans_cut[b_sans]= {519,507,545,528,542,517,533,535,511,529,534,534,532,513,532,539,521,497,494,492,445,436,428};
		Double_t datas_BGO_pourcentage_Same_gamma_histo_sans_cut[b_sans]	={935,928,909,937,945,933,909,887,902,920,916,918,934,897,888,892,888,899,879,830,836,788,774};
		Double_t datas_BGO_pourcentage_Same_Particle_histo_sans_cut[b_sans]	={402,388,370,387,414,405,388,390,395,397,389,387,387,398,383,364,380,368,374,361,347,335,325};


		for (int ff= 0; ff<b_sans; ff++)
		{
			datas_BGO_pourcentage_Fortuitous_histo_sans_cut[ff] = datas_BGO_pourcentage_Fortuitous_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;
			datas_BGO_pourcentage_True_gamma_histo_sans_cut[ff]= datas_BGO_pourcentage_True_gamma_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;
			datas_BGO_pourcentage_True_Other_gamma_histo_sans_cut[ff]= datas_BGO_pourcentage_True_Other_gamma_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;
			datas_BGO_pourcentage_Same_gamma_histo_sans_cut[ff] = datas_BGO_pourcentage_Same_gamma_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;
			datas_BGO_pourcentage_Same_Particle_histo_sans_cut[ff] = datas_BGO_pourcentage_Same_Particle_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;

		}



		Double_t Sum_bad_events_pourcentage_histo_sans_cut[b_sans];

		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_pourcentage_histo_sans_cut[u]= datas_BGO_pourcentage_Fortuitous_histo_sans_cut[u]+datas_BGO_pourcentage_True_gamma_histo_sans_cut[u]-datas_BGO_pourcentage_Same_gamma_histo_sans_cut[u];
		}


		//-------------------
		//-------------------
		//Coinc Absolu
		//-------------------
		//-------------------

		//Variation intensity
		Double_t datas_BGO_absolu_Fortuitous_sans_cut[b_sans]={2e-08,0,8e-08,3.7e-07,7.6e-07,1.36e-06,1.73e-06,2.44e-06,3.04e-06,3.49e-06,3.98e-06,4.37e-06,5.37e-06,5.59e-06,8.56e-06,1.108e-05,1.675e-05,2.683e-05,4.237e-05,5.131e-05,7.321e-05,9.298e-05,9.608e-05};
		Double_t datas_BGO_absolu_True_gamma_sans_cut[b_sans]	={2.943e-05,2.894e-05,2.91e-05,2.922e-05,2.965e-05,2.93e-05,2.903e-05,2.871e-05,2.848e-05,2.902e-05,2.888e-05,2.921e-05,2.922e-05,2.88e-05,2.863e-05,2.855e-05,2.834e-05,2.811e-05,2.769e-05,2.681e-05,2.548e-05,2.501e-05,2.397e-05};
		Double_t datas_BGO_absolu_True_Other_gamma_sans_cut[b_sans]= {8.58e-06,8.23e-06,8.58e-06,8.41e-06,8.55e-06,8.32e-06,8.58e-06,8.46e-06,8.1e-06,8.23e-06,8.45e-06,8.55e-06,8.44e-06,8.36e-06,8.59e-06,8.64e-06,8.42e-06,7.89e-06,7.87e-06,7.72e-06,7.23e-06,7.21e-06,6.72e-06};
		Double_t datas_BGO_absolu_Same_gamma_sans_cut[b_sans]	={1.364e-05,1.356e-05,1.352e-05,1.373e-05,1.38e-05,1.372e-05,1.337e-05,1.312e-05,1.333e-05,1.358e-05,1.349e-05,1.347e-05,1.364e-05,1.327e-05,1.32e-05,1.314e-05,1.325e-05,1.337e-05,1.296e-05,1.254e-05,1.211e-05,1.184e-05,1.145e-05};
		Double_t datas_BGO_absolu_Same_Particle_sans_cut[b_sans]	={7.21e-06,7.15e-06,7e-06,7.08e-06,7.3e-06,7.26e-06,7.08e-06,7.13e-06,7.05e-06,7.21e-06,6.94e-06,7.19e-06,7.14e-06,7.17e-06,6.84e-06,6.77e-06,6.67e-06,6.85e-06,6.86e-06,6.55e-06,6.14e-06,5.96e-06,5.8e-06};

		Double_t Sum_bad_events_absolu_sans_cut[b_sans];

		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_absolu_sans_cut[u]= datas_BGO_absolu_Fortuitous_sans_cut[u]+datas_BGO_absolu_True_gamma_sans_cut[u]-datas_BGO_absolu_Same_gamma_sans_cut[u];
		}



		//--------------------------------
		//Reconstructed events (line-cone)
		//--------------------------------


		Double_t datas_BGO_absolu_total_histo_sans_cut[b_sans]={1856,1823,1828,1870,1934,1915,1917,1926,1945,1998,2026,2041,2093,2060,2216,2282,2576,2970,3682,4050,5042,5857,6008};
		Double_t datas_BGO_absolu_Fortuitous_histo_sans_cut[b_sans]={0,0,4,18,33,60,87,114,137,152,187,202,240,252,413,487,787,1206,1935,2367,3414,4298,4481} ;
		Double_t datas_BGO_absolu_True_gamma_histo_sans_cut[b_sans]	={1856,1823,1824,1852,1901,1855,1830,1812,1808,1846,1839,1839,1853,1808,1803,1795,1789,1764,1747,1683,1628,1559,1527};
		Double_t datas_BGO_absolu_True_Other_gamma_histo_sans_cut[b_sans]= {519,507,545,528,542,517,533,535,511,529,534,534,532,513,532,539,521,497,494,492,445,436,428};
		Double_t datas_BGO_absolu_Same_gamma_histo_sans_cut[b_sans]	={935,928,909,937,945,933,909,887,902,920,916,918,934,897,888,892,888,899,879,830,836,788,774};
		Double_t datas_BGO_absolu_Same_Particle_histo_sans_cut[b_sans]	={402,388,370,387,414,405,388,390,395,397,389,387,387,398,383,364,380,368,374,361,347,335,325};


		/*
		 for (int ff= 0; ff<b; ff++)
		 {
		 datas_BGO_absolu_Fortuitous[ff] = datas_BGO_absolu_Fortuitous[ff]/ (100);
		 datas_BGO_absolu_True_gamma[ff]= datas_BGO_absolu_True_gamma[ff]/ (100);
		 datas_BGO_absolu_True_Other_gamma[ff]= datas_BGO_absolu_True_Other_gamma[ff]/ (100);
		 datas_BGO_absolu_Same_gamma[ff] = datas_BGO_absolu_Same_gamma[ff]/ (100);
		 datas_BGO_absolu_Same_Particle[ff] = datas_BGO_absolu_Same_Particle[ff]/ (100);

		 }
		 */



		for (int ff= 0; ff<b_sans; ff++)
		{
			datas_BGO_absolu_Fortuitous_histo_sans_cut[ff] = datas_BGO_absolu_Fortuitous_histo_sans_cut[ff]/ (1e8);
			datas_BGO_absolu_True_gamma_histo_sans_cut[ff]= datas_BGO_absolu_True_gamma_histo_sans_cut[ff]/ (1e8);
			datas_BGO_absolu_True_Other_gamma_histo_sans_cut[ff]= datas_BGO_absolu_True_Other_gamma_histo_sans_cut[ff]/ (1e8);
			datas_BGO_absolu_Same_gamma_histo_sans_cut[ff] = datas_BGO_absolu_Same_gamma_histo_sans_cut[ff]/ (1e8);
			datas_BGO_absolu_Same_Particle_histo_sans_cut[ff] = datas_BGO_absolu_Same_Particle_histo_sans_cut[ff]/ (1e8);

		}


		Double_t Sum_bad_events_absolu_histo_sans_cut[b_sans];

		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_absolu_histo_sans_cut[u]= datas_BGO_absolu_Fortuitous_histo_sans_cut[u]+datas_BGO_absolu_True_gamma_histo_sans_cut[u]-datas_BGO_absolu_Same_gamma_histo_sans_cut[u];
		}



	//}




	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- IONS CARBONE --------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------



	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- Avec CUT TOF --------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------




	//if (choix_p_c ==2)
	//{


	//-------------------
	//-------------------
	//Coinc Pourcentage
	//-------------------
	//-------------------


/*

	 //Variation intensity
	 Double_t datas_BGO_pourcentage_Fortuitous[b]={0.0916514,0.0674252,0.10073,0.218542,0.936904,1.59061,8.03453,15.1658,26.749,34.9497,41.4776,46.7463,51.2818,55.2236,58.4531,61.3696,63.0908,71.9593,77.176,83.3408,88.7818,90.5962,91.8277,};
	 Double_t datas_BGO_pourcentage_True_gamma[b]	={99.9083,99.9326,99.8993,99.7815,99.0631,98.4094,91.9655,84.8342,73.251,65.0503,58.5224,53.2537,48.7182,44.7764,41.5469,38.6304,36.9092,28.0407,22.824,16.6592,11.2182,9.40381,8.17227};
	 Double_t datas_BGO_pourcentage_True_Other_gamma[b]= {50.6249,50.4762,50.6841,50.1975,49.8549,49.7502,46.6905,42.9577,36.9213,32.686,29.452,26.7247,24.3705,22.2983,20.6409,19.2467,18.3551,13.9141,11.1329,8.01088,5.24641,4.39092,3.76348};
	 Double_t datas_BGO_pourcentage_Same_gamma[b]	={22.3796,22.2672,22.018,22.5183,22.3696,21.8105,20.2691,18.6851,16.2873,14.4127,12.916,11.7869,10.7252,9.89942,9.184,8.4611,8.08081,5.9827,4.90767,3.4794,2.30871,1.95526,1.59813};
	 Double_t datas_BGO_pourcentage_Same_Particle[b]	={26.9038,27.1892,27.1972,27.0656,26.8386,26.8488,25.0058,23.1914,20.0424,17.9516,16.1543,14.7421,13.6224,12.5786,11.722,10.9226,10.4733,8.14388,6.78343,5.16896,3.66304,3.05763,2.81066};

		Double_t Sum_bad_events_pourcentage[b];

		for(int u=0;u<b;u++)
		{
			Sum_bad_events_pourcentage[u]= datas_BGO_pourcentage_Fortuitous[u]+datas_BGO_pourcentage_True_gamma[u]-datas_BGO_pourcentage_Same_gamma[u];
		}

	 //--------------------------------
	 //Reconstructed events (line-cone)
	 //--------------------------------

	 Double_t datas_BGO_pourcentage_total_histo[b]={5294,5279,5280,5292,5313,5263,5568,5876,6558,6987,7612,8141,8600,9025,9618,9947,10303,12231,13799,16069,19363,20706,21881};
	 Double_t datas_BGO_pourcentage_Fortuitous_histo[b]={4,1,4,14,32,63,349,703,1464,2114,2759,3316,3886,4423,5042,5499,5874,8174,10048,12843,16654,18288,19651} ;
	 Double_t datas_BGO_pourcentage_True_gamma_histo[b]	={5290,5278,5276,5278,5281,5200,5219,5173,5094,4873,4853,4825,4714,4602,4576,4448,4429,4057,3751,3226,2709,2418,2230};
	 Double_t datas_BGO_pourcentage_True_Other_gamma_histo[b]= {2435,2426,2451,2410,2415,2393,2413,2388,2348,2230,2217,2207,2141,2068,2082,2009,1998,1819,1655,1423,1161,1026,949};
	 Double_t datas_BGO_pourcentage_Same_gamma_histo[b]	={1492,1474,1453,1507,1497,1431,1441,1431,1409,1357,1339,1334,1311,1273,1262,1220,1223,1084,1025,823,679,597,503};
	 Double_t datas_BGO_pourcentage_Same_Particle_histo[b]	={1363,1378,1372,1361,1369,1376,1365,1354,1337,1286,1297,1284,1262,1261,1232,1219,1208,1154,1071,980,869,795,778};



	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_pourcentage_Fortuitous_histo[ff] = datas_BGO_pourcentage_Fortuitous_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_True_gamma_histo[ff]= datas_BGO_pourcentage_True_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_True_Other_gamma_histo[ff]= datas_BGO_pourcentage_True_Other_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_Same_gamma_histo[ff] = datas_BGO_pourcentage_Same_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_Same_Particle_histo[ff] = datas_BGO_pourcentage_Same_Particle_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;

	}


		Double_t Sum_bad_events_pourcentage_histo[b];

		for(int u=0;u<b;u++)
		{
			Sum_bad_events_pourcentage_histo[u]= datas_BGO_pourcentage_Fortuitous_histo[u]+datas_BGO_pourcentage_True_gamma_histo[u]-datas_BGO_pourcentage_Same_gamma_histo[u];
		//	cout<<"test "<<Sum_bad_events_pourcentage_histo[u]<<endl;
		}

	//-------------------
	//-------------------
	//Coinc Absolu
	//-------------------
	//-------------------


	 //Variation intensity
	 Double_t datas_BGO_absolu_Fortuitous[b]={5.5e-07,4e-07,6e-07,1.3e-06,5.65e-06,9.55e-06,5.165e-05,0.0001045,0.000208,0.00029875,0.0003913,0.0004777,0.0005611,0.00064515,0.0007278,0.00080655,0.0008557,0.00118535,0.00145855,0.0018695,0.0024419,0.0026874,0.0028816};
	 Double_t datas_BGO_absolu_True_gamma[b]	={0.00059955,0.00059285,0.00059505,0.00059355,0.0005974,0.00059085,0.0005912,0.00058455,0.0005696,0.00055605,0.0005521,0.0005442,0.00053305,0.0005231,0.0005173,0.0005077,0.0005006,0.0004619,0.00043135,0.0003737,0.00030855,0.00027895,0.00025645};
	 Double_t datas_BGO_absolu_True_Other_gamma[b]= {0.0003038,0.00029945,0.0003019,0.0002986,0.00030065,0.0002987,0.00030015,0.000296,0.0002871,0.0002794,0.00027785,0.0002731,0.00026665,0.0002605,0.000257,0.00025295,0.00024895,0.0002292,0.0002104,0.0001797,0.0001443,0.00013025,0.0001181};
	 Double_t datas_BGO_absolu_Same_gamma[b]	={0.0001343,0.0001321,0.00013115,0.00013395,0.0001349,0.00013095,0.0001303,0.00012875,0.00012665,0.0001232,0.00012185,0.00012045,0.00011735,0.00011565,0.00011435,0.0001112,0.0001096,9.855e-05,9.275e-05,7.805e-05,6.35e-05,5.8e-05,5.015e-05};
	 Double_t datas_BGO_absolu_Same_Particle[b]	={0.00016145,0.0001613,0.000162,0.000161,0.00016185,0.0001612,0.00016075,0.0001598,0.00015585,0.00015345,0.0001524,0.00015065,0.00014905,0.00014695,0.00014595,0.00014355,0.00014205,0.00013415,0.0001282,0.00011595,0.00010075,9.07e-05,8.82e-05};

		Double_t Sum_bad_events_absolu[b];

		for(int u=0;u<b;u++)
		{
			Sum_bad_events_absolu[u]= datas_BGO_absolu_Fortuitous[u]+datas_BGO_absolu_True_gamma[u]-datas_BGO_absolu_Same_gamma[u];
		}


	//--------------------------------
	//Reconstructed events (line-cone)
	//--------------------------------


	 Double_t datas_BGO_absolu_total_histo[b]={5294,5279,5280,5292,5313,5263,5568,5876,6558,6987,7612,8141,8600,9025,9618,9947,10303,12231,13799,16069,19363,20706,21881};
        Double_t datas_BGO_absolu_Fortuitous_histo[b]={4,1,4,14,32,63,349,703,1464,2114,2759,3316,3886,4423,5042,5499,5874,8174,10048,12843,16654,18288,19651} ;
	 Double_t datas_BGO_absolu_True_gamma_histo[b]	={5290,5278,5276,5278,5281,5200,5219,5173,5094,4873,4853,4825,4714,4602,4576,4448,4429,4057,3751,3226,2709,2418,2230};
        Double_t datas_BGO_absolu_True_Other_gamma_histo[b]={2435,2426,2451,2410,2415,2393,2413,2388,2348,2230,2217,2207,2141,2068,2082,2009,1998,1819,1655,1423,1161,1026,949};
	 Double_t datas_BGO_absolu_Same_gamma_histo[b]	={1492,1474,1453,1507,1497,1431,1441,1431,1409,1357,1339,1334,1311,1273,1262,1220,1223,1084,1025,823,679,597,503};
	 Double_t datas_BGO_absolu_Same_Particle_histo[b]	={1363,1378,1372,1361,1369,1376,1365,1354,1337,1286,1297,1284,1262,1261,1232,1219,1208,1154,1071,980,869,795,778};





	/*
	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_absolu_Fortuitous[ff] = datas_BGO_absolu_Fortuitous[ff]/ (100);
		datas_BGO_absolu_True_gamma[ff]= datas_BGO_absolu_True_gamma[ff]/ (100);
		datas_BGO_absolu_True_Other_gamma[ff]= datas_BGO_absolu_True_Other_gamma[ff]/ (100);
		datas_BGO_absolu_Same_gamma[ff] = datas_BGO_absolu_Same_gamma[ff]/ (100);
		datas_BGO_absolu_Same_Particle[ff] = datas_BGO_absolu_Same_Particle[ff]/ (100);

	}
	*/

/*

	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_absolu_Fortuitous_histo[ff] = datas_BGO_absolu_Fortuitous_histo[ff]/ (2e7);
		datas_BGO_absolu_True_gamma_histo[ff]= datas_BGO_absolu_True_gamma_histo[ff]/ (2e7);
		datas_BGO_absolu_True_Other_gamma_histo[ff]= datas_BGO_absolu_True_Other_gamma_histo[ff]/ (2e7);
		datas_BGO_absolu_Same_gamma_histo[ff] = datas_BGO_absolu_Same_gamma_histo[ff]/ (2e7);
		datas_BGO_absolu_Same_Particle_histo[ff] = datas_BGO_absolu_Same_Particle_histo[ff]/ (2e7);

	}


		Double_t Sum_bad_events_absolu_histo[b];
		for(int u=0;u<b;u++)
		{
			Sum_bad_events_absolu_histo[u]= datas_BGO_absolu_Fortuitous_histo[u]+datas_BGO_absolu_True_gamma_histo[u]-datas_BGO_absolu_Same_gamma_histo[u];
		}



	//}




	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- SANS CUT TOF --------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------


	//if (choix_p_c ==2)
	//{


		//-------------------
		//-------------------
		//Coinc Pourcentage
		//-------------------
		//-------------------




		//Variation intensity
		Double_t datas_BGO_pourcentage_Fortuitous_sans_cut[b_sans]={0.0913743,0.0600018,0.0597372,0.161693,0.666451,1.30275,6.12885,11.6696,20.8771,28.1004,34.058,39.1236,43.4917,47.2119,50.6452,53.3086,55.5481,64.7907,70.822,77.6289,84.7799,86.8133,88.2589 };
		Double_t datas_BGO_pourcentage_True_gamma_sans_cut[b_sans]	={99.9086,99.94,99.9403,99.8383,99.3335,98.6972,93.8712,88.3304,79.1229,71.8996,65.942,60.8764,56.5083,52.7881,49.3548,46.6914,44.4519,35.2093,29.178,22.3711,15.2201,13.1867,11.7411 };
        Double_t datas_BGO_pourcentage_True_Other_gamma_sans_cut[b_sans]= {37.4497,37.4227,37.3909,37.3556,36.9834,36.867,35.0089,33.0702,29.4157,26.6104,24.4377,22.4645,20.7062,19.2652,17.9339,16.9947,16.0539,12.4243,10.0482,7.58789,4.93234,4.32076,3.75553};
		Double_t datas_BGO_pourcentage_Same_gamma_sans_cut[b_sans]	={12.1482,12.2127,12.3978,12.455,12.2229,12.2018,11.3429,10.7313,9.63127,8.68169,7.94894,7.27113,6.68903,6.24303,5.87086,5.53039,5.2595,4.07162,3.34458,2.51626,1.71654,1.3948,1.23414 };
		Double_t datas_BGO_pourcentage_Same_Particle_sans_cut[b_sans]	={ 50.3107,50.3046,50.1516,50.0277,50.1273,49.6284,47.5193,44.529,40.0758,36.6076,33.5553,31.1408,29.113,27.2799,25.5501,24.1664,23.1385,18.7134,15.7852,12.267,8.57125,7.47113,6.75141};

		Double_t Sum_bad_events_pourcentage_sans_cut[b_sans];

		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_pourcentage_sans_cut[u]= datas_BGO_pourcentage_Fortuitous_sans_cut[u]+datas_BGO_pourcentage_True_gamma_sans_cut[u]-datas_BGO_pourcentage_Same_gamma_sans_cut[u];
		}

		//--------------------------------
		//Reconstructed events (line-cone)
		//--------------------------------

		Double_t datas_BGO_pourcentage_total_histo_sans_cut[b_sans]={ 10040,9986,10073,9988,9959,10036,10410,10784,11690,12250,12888,13695,14058,14908,15573,15948,16448,18732,20640,23402,28050,29252,30498};
		Double_t datas_BGO_pourcentage_Fortuitous_histo_sans_cut[b_sans]={8,4,6,20,43,100,478,944,1924,2747,3581,4418,5047,5984,6698,7277,7899,10825,13253,16765,22454,24189,25734} ;
		Double_t datas_BGO_pourcentage_True_gamma_histo_sans_cut[b_sans]	={ 10032,9982,10067,9968,9916,9936,9932,9840,9766,9503,9307,9277,9011,8924,8875,8671,8549,7907,7387,6637,5596,5063,4764};
		Double_t datas_BGO_pourcentage_True_Other_gamma_histo_sans_cut[b_sans]= {3194,3222,3246,3177,3172,3162,3172,3158,3114,2993,2956,2936,2803,2763,2775,2676,2642,2364,2161,1881,1496,1376,1286 };
		Double_t datas_BGO_pourcentage_Same_gamma_histo_sans_cut[b_sans]	={1455,1458,1518,1487,1471,1480,1450,1432,1420,1367,1341,1334,1286,1258,1261,1243,1200,1083,987,854,723,608,537};
		Double_t datas_BGO_pourcentage_Same_Particle_histo_sans_cut[b_sans]	={ 5383,5302,5303,5304,5273,5294,5310,5250,5232,5143,5010,5007,4922,4903,4839,4752,4707,4460,4239,3902,3377,3079,2941};



		for (int ff= 0; ff<b_sans; ff++)
		{
			datas_BGO_pourcentage_Fortuitous_histo_sans_cut[ff] = datas_BGO_pourcentage_Fortuitous_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;
			datas_BGO_pourcentage_True_gamma_histo_sans_cut[ff]= datas_BGO_pourcentage_True_gamma_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;
			datas_BGO_pourcentage_True_Other_gamma_histo_sans_cut[ff]= datas_BGO_pourcentage_True_Other_gamma_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;
			datas_BGO_pourcentage_Same_gamma_histo_sans_cut[ff] = datas_BGO_pourcentage_Same_gamma_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;
			datas_BGO_pourcentage_Same_Particle_histo_sans_cut[ff] = datas_BGO_pourcentage_Same_Particle_histo_sans_cut[ff]/ datas_BGO_pourcentage_total_histo_sans_cut[ff]*100;

		}


		Double_t Sum_bad_events_pourcentage_histo_sans_cut[b_sans];

		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_pourcentage_histo_sans_cut[u]= datas_BGO_pourcentage_Fortuitous_histo_sans_cut[u]+datas_BGO_pourcentage_True_gamma_histo_sans_cut[u]-datas_BGO_pourcentage_Same_gamma_histo_sans_cut[u];
			//	cout<<"test "<<Sum_bad_events_pourcentage_histo[u]<<endl;
		}

		//-------------------
		//-------------------
		//Coinc Absolu
		//-------------------
		//-------------------


		//Variation intensity
		Double_t datas_BGO_absolu_Fortuitous_sans_cut[b_sans]={1e-06,6.5e-07,6.5e-07,1.75e-06,7.2e-06,1.42e-05,7.035e-05,0.00014115,0.000278,0.00040265,0.00052165,0.00064595,0.00075195,0.0008674,0.00098515,0.0010767,0.00115965,0.00157775,0.00195235,0.00248195,0.0033363,0.0035975,0.00382245 };
		Double_t datas_BGO_absolu_True_gamma_sans_cut[b_sans]	={0.0010934,0.00108265,0.00108745,0.00108055,0.00107315,0.0010758,0.0010775,0.0010684,0.0010536,0.00103025,0.00101,0.0010051,0.000977,0.00096985,0.00096005,0.00094305,0.000928,0.0008574,0.00080435,0.00071525,0.00059895,0.00054645,0.0005085 };
		Double_t datas_BGO_absolu_True_Other_gamma_sans_cut[b_sans]= {0.00040985,0.0004054,0.00040685,0.0004043,0.00039955,0.00040185,0.00040185,0.0004,0.0003917,0.0003813,0.0003743,0.0003709,0.000358,0.00035395,0.00034885,0.00034325,0.00033515,0.00030255,0.000277,0.0002426,0.0001941,0.00017905,0.00016265 };
		Double_t datas_BGO_absolu_Same_gamma_sans_cut[b_sans]	={0.00013295,0.0001323,0.0001349,0.0001348,0.00013205,0.000133,0.0001302,0.0001298,0.00012825,0.0001244,0.00012175,0.00012005,0.00011565,0.0001147,0.0001142,0.0001117,0.0001098,9.915e-05,9.22e-05,8.045e-05,6.755e-05,5.78e-05,5.345e-05 };
		Double_t datas_BGO_absolu_Same_Particle_sans_cut[b_sans]	={0.0005506,0.00054495,0.0005457,0.00054145,0.00054155,0.00054095,0.00054545,0.0005386,0.00053365,0.00052455,0.00051395,0.00051415,0.00050335,0.0005012,0.000497,0.0004881,0.00048305,0.0004557,0.00043515,0.0003922,0.0003373,0.0003096,0.0002924 };

		Double_t Sum_bad_events_absolu_sans_cut[b_sans];

		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_absolu_sans_cut[u]= datas_BGO_absolu_Fortuitous_sans_cut[u]+datas_BGO_absolu_True_gamma_sans_cut[u]-datas_BGO_absolu_Same_gamma_sans_cut[u];
		}


		//--------------------------------
		//Reconstructed events (line-cone)
		//--------------------------------


		Double_t datas_BGO_absolu_total_histo_sans_cut[b_sans]={ 10040,9986,10073,9988,9959,10036,10410,10784,11690,12250,12888,13695,14058,14908,15573,15948,16448,18732,20640,23402,28050,29252,30498};
		Double_t datas_BGO_absolu_Fortuitous_histo_sans_cut[b_sans]={8,4,6,20,43,100,478,944,1924,2747,3581,4418,5047,5984,6698,7277,7899,10825,13253,16765,22454,24189,25734} ;
		Double_t datas_BGO_absolu_True_gamma_histo_sans_cut[b_sans]	={ 10032,9982,10067,9968,9916,9936,9932,9840,9766,9503,9307,9277,9011,8924,8875,8671,8549,7907,7387,6637,5596,5063,4764};
		Double_t datas_BGO_absolu_True_Other_gamma_histo_sans_cut[b_sans]= {3194,3222,3246,3177,3172,3162,3172,3158,3114,2993,2956,2936,2803,2763,2775,2676,2642,2364,2161,1881,1496,1376,1286 };
		Double_t datas_BGO_absolu_Same_gamma_histo_sans_cut[b_sans]	={1455,1458,1518,1487,1471,1480,1450,1432,1420,1367,1341,1334,1286,1258,1261,1243,1200,1083,987,854,723,608,537};
		Double_t datas_BGO_absolu_Same_Particle_histo_sans_cut[b_sans]	={ 5383,5302,5303,5304,5273,5294,5310,5250,5232,5143,5010,5007,4922,4903,4839,4752,4707,4460,4239,3902,3377,3079,2941};





		/*
		 for (int ff= 0; ff<b; ff++)
		 {
		 datas_BGO_absolu_Fortuitous[ff] = datas_BGO_absolu_Fortuitous[ff]/ (100);
		 datas_BGO_absolu_True_gamma[ff]= datas_BGO_absolu_True_gamma[ff]/ (100);
		 datas_BGO_absolu_True_Other_gamma[ff]= datas_BGO_absolu_True_Other_gamma[ff]/ (100);
		 datas_BGO_absolu_Same_gamma[ff] = datas_BGO_absolu_Same_gamma[ff]/ (100);
		 datas_BGO_absolu_Same_Particle[ff] = datas_BGO_absolu_Same_Particle[ff]/ (100);

		 }
		 */

/*

		for (int ff= 0; ff<b_sans; ff++)
		{
			datas_BGO_absolu_Fortuitous_histo_sans_cut[ff] = datas_BGO_absolu_Fortuitous_histo_sans_cut[ff]/ (2e7);
			datas_BGO_absolu_True_gamma_histo_sans_cut[ff]= datas_BGO_absolu_True_gamma_histo_sans_cut[ff]/ (2e7);
			//datas_BGO_absolu_True_Other_gamma_histo[ff]= datas_BGO_absolu_True_Other_gamma_histo[ff]/ (1e8);
			datas_BGO_absolu_Same_gamma_histo_sans_cut[ff] = datas_BGO_absolu_Same_gamma_histo_sans_cut[ff]/ (2e7);
			datas_BGO_absolu_Same_Particle_histo_sans_cut[ff] = datas_BGO_absolu_Same_Particle_histo_sans_cut[ff]/ (2e7);

		}


		Double_t Sum_bad_events_absolu_histo_sans_cut[b_sans];
		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_absolu_histo_sans_cut[u]= datas_BGO_absolu_Fortuitous_histo_sans_cut[u]+datas_BGO_absolu_True_gamma_histo_sans_cut[u]-datas_BGO_absolu_Same_gamma_histo_sans_cut[u];
		}



	//}

*/
	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- Traitement datas ----------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------


	// Calcul du ratio signal gamma true vs le background assimilé à tout le reste - pour proton et carbon

	if (choix_p_c ==1)
	{

		for(int l=0;l<b;l++)
		{
			Ratio_Background_Signal_proton_absolu [l]=datas_BGO_absolu_Same_gamma[l]/Sum_bad_events_absolu[l];
			Ratio_Background_Signal_proton_relatif[l]=datas_BGO_pourcentage_Same_gamma[l]/Sum_bad_events_pourcentage[l];
			Ratio_Background_Signal_proton_histo[l]=datas_BGO_pourcentage_Same_gamma_histo[l]/Sum_bad_events_pourcentage_histo[l];

			for(int h=0;h<b_sans;h++)
			{

					Ratio_Background_Signal_proton_absolu_comparaison_gamma [l]=datas_BGO_absolu_Same_gamma_sans_cut[l]/datas_BGO_absolu_Same_gamma[l];
					Ratio_Background_Signal_proton_relatif_comparaison_gamma[l]=datas_BGO_pourcentage_Same_gamma_sans_cut[l]/datas_BGO_pourcentage_Same_gamma[l];
					Ratio_Background_Signal_proton_histo_comparaison_gamma[l]=datas_BGO_pourcentage_Same_gamma_histo_sans_cut[l]/datas_BGO_pourcentage_Same_gamma_histo[l];

					Ratio_Background_Signal_proton_absolu_comparaison_bad [l]=Sum_bad_events_absolu_sans_cut[l]/Sum_bad_events_absolu[l];
					Ratio_Background_Signal_proton_relatif_comparaison_bad[l]=Sum_bad_events_pourcentage_sans_cut[l]/Sum_bad_events_pourcentage[l];
					Ratio_Background_Signal_proton_histo_comparaison_bad[l]=Sum_bad_events_pourcentage_histo_sans_cut[l]/Sum_bad_events_pourcentage_histo[l];

			}
		}
		for(int p=0;p<b_sans;p++)
		{
			Ratio_Background_Signal_proton_absolu_sans_cut [p]=datas_BGO_absolu_Same_gamma_sans_cut[p]/Sum_bad_events_absolu_sans_cut[p];
			Ratio_Background_Signal_proton_relatif_sans_cut[p]=datas_BGO_pourcentage_Same_gamma_sans_cut[p]/Sum_bad_events_pourcentage_sans_cut[p];
			Ratio_Background_Signal_proton_histo_sans_cut[p]=datas_BGO_pourcentage_Same_gamma_histo_sans_cut[p]/Sum_bad_events_pourcentage_histo_sans_cut[p];
		}

		for(int l=0;l<b;l++)
		{
			for(int h=0;h<b_sans;h++)
			{

					Ratio_Background_Signal_proton_absolu_comparaison[l]=Ratio_Background_Signal_proton_absolu[l]/Ratio_Background_Signal_proton_absolu_sans_cut[l];
					Ratio_Background_Signal_proton_relatif_comparaison[l]=Ratio_Background_Signal_proton_relatif[l]/Ratio_Background_Signal_proton_relatif_sans_cut[l];
					Ratio_Background_Signal_proton_histo_comparaison[l]= Ratio_Background_Signal_proton_histo[l]/Ratio_Background_Signal_proton_histo_sans_cut[l];

			}
		}



	}
	if (choix_p_c ==2)
	{
		for(int l=0;l<b;l++)
		{
			Ratio_Background_Signal_carbon_absolu[l]=datas_BGO_absolu_Same_gamma[l]/Sum_bad_events_absolu[l];
			Ratio_Background_Signal_carbon_relatif[l]=datas_BGO_pourcentage_Same_gamma[l]/Sum_bad_events_pourcentage[l];
			Ratio_Background_Signal_carbon_histo[l]=datas_BGO_pourcentage_Same_gamma_histo[l]/Sum_bad_events_pourcentage_histo[l];

			for(int h=0;h<b_sans;h++)
			{

					Ratio_Background_Signal_carbon_absolu_comparaison_gamma [l]=datas_BGO_absolu_Same_gamma_sans_cut[l]/datas_BGO_absolu_Same_gamma[l];
					Ratio_Background_Signal_carbon_relatif_comparaison_gamma[l]=datas_BGO_pourcentage_Same_gamma_sans_cut[l]/datas_BGO_pourcentage_Same_gamma[l];
					Ratio_Background_Signal_carbon_histo_comparaison_gamma[l]=datas_BGO_pourcentage_Same_gamma_histo_sans_cut[l]/datas_BGO_pourcentage_Same_gamma_histo[l];

					Ratio_Background_Signal_carbon_absolu_comparaison_bad [l]=Sum_bad_events_absolu_sans_cut[l]/Sum_bad_events_absolu[l];
					Ratio_Background_Signal_carbon_relatif_comparaison_bad[l]=Sum_bad_events_pourcentage_sans_cut[l]/Sum_bad_events_pourcentage[l];
					Ratio_Background_Signal_carbon_histo_comparaison_bad[l]=Sum_bad_events_pourcentage_histo_sans_cut[l]/Sum_bad_events_pourcentage_histo[l];

			}
		}
		for(int p=0;p<b_sans;p++)
		{
			Ratio_Background_Signal_carbon_absolu_sans_cut  [p]=datas_BGO_absolu_Same_gamma_sans_cut[p]/Sum_bad_events_absolu_sans_cut[p];
			Ratio_Background_Signal_carbon_relatif_sans_cut [p]=datas_BGO_pourcentage_Same_gamma_sans_cut[p]/Sum_bad_events_pourcentage_sans_cut[p];
			Ratio_Background_Signal_carbon_histo_sans_cut [p]=datas_BGO_pourcentage_Same_gamma_histo_sans_cut[p]/Sum_bad_events_pourcentage_histo_sans_cut[p];
		}

		for(int l=0;l<b;l++)
		{
			for(int h=0;h<b_sans;h++)
			{

					Ratio_Background_Signal_carbon_absolu_comparaison[l]=Ratio_Background_Signal_carbon_absolu[l]/Ratio_Background_Signal_carbon_absolu_sans_cut[l];
					Ratio_Background_Signal_carbon_relatif_comparaison[l]=Ratio_Background_Signal_carbon_relatif[l]/Ratio_Background_Signal_carbon_relatif_sans_cut[l];
					Ratio_Background_Signal_carbon_histo_comparaison[l]= Ratio_Background_Signal_carbon_histo[l]/Ratio_Background_Signal_carbon_histo_sans_cut[l];

			}
		}


	}






	//Affichage des datas calculées
	if (choix_p_c ==1)
	{
		//Avec Cut
		cout <<" Ratio_Background_Signal_proton_absolu "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_absolu [g]<<endl;}
		cout <<" Ratio_Background_Signal_proton_relatif "<<endl;
		cout<<endl<<endl;
		for(int g =0;g<b;g++) {cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_relatif [g]<<endl;}
		cout <<" Ratio_Background_Signal_proton_histo "<<endl;
		cout<<endl<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_histo [g]<<endl;}
		cout<<endl<<endl;

		//Sans Cut
		cout <<" Ratio_Background_Signal_proton_absolu_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++) {cout<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_proton_absolu_sans_cut [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_proton_relatif_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++) {cout<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_proton_relatif_sans_cut [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_proton_histo_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++){ cout<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_proton_histo_sans_cut [g]<<endl;}
		cout<<endl<<endl;

		//Comparaison
		//Bons events -gamma
		cout <<" Ratio_Background_Signal_proton_absolu_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_absolu_comparaison_gamma [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_proton_relatif_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++) {cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_relatif_comparaison_gamma [g]<<endl;}
	    cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_proton_histo_comparaison_gamma "<<endl;
		cout<<endl<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_histo_comparaison_gamma [g]<<endl;}
		cout<<endl<<endl;
		//Bad events - background
		cout <<" Ratio_Background_Signal_proton_absolu_comparaison_bad "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_absolu_comparaison_bad [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_proton_relatif_comparaison_bad "<<endl;
		for(int g =0;g<b;g++) {cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_relatif_comparaison_bad [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_proton_histo_comparaison_bad "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_histo_comparaison_bad [g]<<endl;}
		cout<<endl<<endl;
		//Rapport ratios S/B avec et sans cut TOF
		cout <<" Ratio_Background_Signal_proton_absolu_comparaison "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_absolu_comparaison [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_proton_relatif_comparaison "<<endl;
		for(int g =0;g<b;g++) {cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_relatif_comparaison [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_proton_histo_comparaison "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_histo_comparaison [g]<<endl;}
		cout<<endl<<endl;

	}

	if (choix_p_c ==2)
	{
		//Avec cut TOF
		cout <<" Ratio_Background_Signal_carbon_absolu "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_absolu [g]<<endl;}
		cout <<" Ratio_Background_Signal_carbon_relatif "<<endl;
		cout<<endl<<endl;
		for(int g =0;g<b;g++) {cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_relatif [g]<<endl;}
		cout <<" Ratio_Background_Signal_carbon_histo "<<endl;
		cout<<endl<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_histo [g]<<endl;}
		cout<<endl<<endl;

		//Sans Cut TOF
		cout <<" Ratio_Background_Signal_carbon_absolu_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++) {cout<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_carbon_absolu_sans_cut [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_carbon_relatif_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++) {cout<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_carbon_relatif_sans_cut [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_carbon_histo_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++){ cout<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_carbon_histo_sans_cut [g]<<endl;}
		cout<<endl<<endl;

		//Comparaison
		//Bons events -gamma
		cout <<" Ratio_Background_Signal_carbon_absolu_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_absolu_comparaison_gamma [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_carbon_relatif_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++) {cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_relatif_comparaison_gamma [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_carbon_histo_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_histo_comparaison_gamma [g]<<endl;}
		cout<<endl<<endl;
		//Bad events - background
		cout <<" Ratio_Background_Signal_carbon_absolu_comparaison_bad "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_absolu_comparaison_bad [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_carbon_relatif_comparaison_bad "<<endl;
		for(int g =0;g<b;g++) {cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_relatif_comparaison_bad [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_carbon_histo_comparaison_bad "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_histo_comparaison_bad [g]<<endl;}
		cout<<endl<<endl;
		//Rapport ratios S/B avec et sans cut TOF
		cout <<" Ratio_Background_Signal_carbon_absolu_comparaison "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_absolu_comparaison [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_carbon_relatif_comparaison "<<endl;
		for(int g =0;g<b;g++) {cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_relatif_comparaison [g]<<endl;}
		cout<<endl<<endl;
		cout <<" Ratio_Background_Signal_carbon_histo_comparaison "<<endl;
		for(int g =0;g<b;g++){ cout<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_histo_comparaison [g]<<endl;}
		cout<<endl<<endl;
	}

	//----------------------Sortie des datas en fichier texte----------------------------

	/*if (choix_p_c ==1)
	{
	ofstream
	fichier("2014_09_16_Signal_to_Background_Ratio_Coincidences_Hadrontherapy_Avec_Sans_CUT_TOF_Protons.txt",ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt

		//Avec Cut
		fichier <<"Ratio_Background_Signal_proton_absolu "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_absolu [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_relatif "<<endl;
		for(int g =0;g<b;g++) {fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_relatif [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_histo "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_histo [g]<<endl;}
		fichier<<endl<<endl;

		//Sans Cut
		fichier <<"Ratio_Background_Signal_proton_absolu_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++) {fichier<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_proton_absolu_sans_cut [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_relatif_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++) {fichier<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_proton_relatif_sans_cut [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_histo_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++){ fichier<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_proton_histo_sans_cut [g]<<endl;}
		fichier<<endl<<endl;

		//Comparaison
		//Bons events -gamma
		fichier <<"Ratio_Background_Signal_proton_absolu_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_absolu_comparaison_gamma [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_relatif_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++) {fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_relatif_comparaison_gamma [g]<<endl;}
	    fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_histo_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_histo_comparaison_gamma [g]<<endl;}
		fichier<<endl<<endl;
		//Bad events - background
		fichier <<"Ratio_Background_Signal_proton_absolu_comparaison_bad "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_absolu_comparaison_bad [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_relatif_comparaison_bad "<<endl;
		for(int g =0;g<b;g++) {fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_relatif_comparaison_bad [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_histo_comparaison_bad "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_histo_comparaison_bad [g]<<endl;}
		fichier<<endl<<endl;
		//Rapport ratios S/B avec et sans cut TOF
		fichier <<"Ratio_Background_Signal_proton_absolu_comparaison "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_absolu_comparaison [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_relatif_comparaison "<<endl;
		for(int g =0;g<b;g++) {fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_relatif_comparaison [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_proton_histo_comparaison "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_proton_histo_comparaison [g]<<endl;}
		fichier<<endl<<endl;
		fichier.close();
	}
*/
	if (choix_p_c ==2)
	{
		ofstream
		fichier("2014_09_16_Signal_to_Background_Ratio_Coincidences_Hadrontherapy_Avec_Sans_CUT_TOF_Carbon_Ions.txt",ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt

		//Avec cut TOF
		fichier <<"Ratio_Background_Signal_carbon_absolu "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_absolu [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_carbon_relatif "<<endl;
		for(int g =0;g<b;g++) {fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_relatif [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_carbon_histo "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_histo [g]<<endl;}
		fichier<<endl<<endl;

		//Sans Cut TOF
		fichier <<"Ratio_Background_Signal_carbon_absolu_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++) {fichier<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_carbon_absolu_sans_cut [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_carbon_relatif_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++) {fichier<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_carbon_relatif_sans_cut [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_carbon_histo_sans_cut "<<endl;
		for(int g =0;g<b_sans;g++){ fichier<< "num "<<sigma_ns_sans_cut[g]<<" value "<<Ratio_Background_Signal_carbon_histo_sans_cut [g]<<endl;}
		fichier<<endl<<endl;

		//Comparaison
		//Bons events -gamma
		fichier <<"Ratio_Background_Signal_carbon_absolu_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_absolu_comparaison_gamma [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_carbon_relatif_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++) {fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_relatif_comparaison_gamma [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_carbon_histo_comparaison_gamma "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_histo_comparaison_gamma [g]<<endl;}
		fichier<<endl<<endl;
		//Bad events - background
		fichier <<"Ratio_Background_Signal_carbon_absolu_comparaison_bad "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_absolu_comparaison_bad [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_carbon_relatif_comparaison_bad "<<endl;
		for(int g =0;g<b;g++) {fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_relatif_comparaison_bad [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<" Ratio_Background_Signal_carbon_histo_comparaison_bad "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_histo_comparaison_bad [g]<<endl;}
		fichier<<endl<<endl;
		//Rapport ratios S/B avec et sans cut TOF
		fichier <<"Ratio_Background_Signal_carbon_absolu_comparaison "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_absolu_comparaison [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_carbon_relatif_comparaison "<<endl;
		for(int g =0;g<b;g++) {fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_relatif_comparaison [g]<<endl;}
		fichier<<endl<<endl;
		fichier <<"Ratio_Background_Signal_carbon_histo_comparaison "<<endl;
		for(int g =0;g<b;g++){ fichier<< "num "<<sigma_ns[g]<<" value "<<Ratio_Background_Signal_carbon_histo_comparaison [g]<<endl;}
		fichier<<endl<<endl;
		fichier.close();
	}




	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- GRAPHIQUES ----------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------



/*
	if(choix==1)
	{

		char name[150];
	      sprintf(name, " Variation of the beam intensity carbon ion");

	      TH1 *frame = new TH1F("frame","",0,0,50);
	      frame->SetMinimum(0);
	      frame->SetMaximum(100);
	      frame->SetDirectory(0);
	      frame->SetStats(0);
	      frame->SetTitle(name);
	      frame->GetXaxis()->SetTitle("Number of carbon ions per bunch");
	      frame->GetXaxis()->SetTickLength(0.02);
	      frame->GetXaxis()->SetLabelSize(0.03);
	      frame->GetYaxis()->SetTitle("Coincidences %");
	      frame->GetYaxis()->SetLabelSize(0.03);
	      frame->Draw(" ");

	      // Échelle logarithmique sur les axes y
	      c1->SetGrid();
		 // c1->SetLogy();


	      Double_t RMS1=0;


	      TGraph *gr = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_Fortuitous,0,0); // Création du graphique

	      gr->SetLineColor(1); // Options de mise en forme  du graphique
	      gr->SetLineStyle(21);
	      gr->SetLineWidth(2);
	      gr->Draw("l");

		TGraph *gr1 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_True_gamma,0,0); // Création du graphique

		gr1->SetLineColor(4); // Options de mise en forme  du graphique
		gr1->SetLineStyle(21);
		gr1->SetLineWidth(2);
		gr1->Draw("l");


	      TGraph *gr2 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_Same_gamma,0,0);//error_FWHM_BGO); // Création du graphique

	      gr2->SetLineColor(2); // Options de mise en forme  du graphique
	      gr2->SetLineStyle(21);
	      gr2->SetLineWidth(2);
	      gr2->Draw("l"); // Dessine le graph


	      TGraph *gr3 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_Same_Particle,0,0); // Création du graphique

	      gr3->SetLineColor(6); // Options de mise en forme  du graphique
	      gr3->SetLineStyle(21);
	      gr3->SetLineWidth(2);
	      gr3->Draw("l"); // Dessine le graph



		/*TGraph *gr4 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_True_Other_gamma,0,0); // Création du graphique

		gr4->SetLineColor(8); // Options de mise en forme  du graphique
		gr4->SetLineStyle(21);
		gr4->SetLineWidth(2);
		gr4->Draw("l"); // Dessine le graph*/


	      //RMS1 = gr3->GetRMS(2);
	      //cout<<"RMS1 "<<RMS1<<endl;


/*
		TGraph *gr5 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_Fortuitous_histo,0,0); // Création du graphique

		gr5->SetLineColor(1); // Options de mise en forme  du graphique
		gr5->SetLineStyle(7);
		gr5->SetLineWidth(2);
		gr5->Draw("l");

		TGraph *gr6 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_True_gamma_histo,0,0); // Création du graphique

		gr6->SetLineColor(4); // Options de mise en forme  du graphique
		gr6->SetLineStyle(7);
		gr6->SetLineWidth(2);
		gr6->Draw("l");


		TGraph *gr7 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_Same_gamma_histo,0,0);//error_FWHM_BGO); // Création du graphique

		gr7->SetLineColor(2); // Options de mise en forme  du graphique
		gr7->SetLineStyle(7);
		gr7->SetLineWidth(2);
		gr7->Draw("l"); // Dessine le graph


		TGraph *gr8 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_Same_Particle_histo,0,0); // Création du graphique

		gr8->SetLineColor(6); // Options de mise en forme  du graphique
		gr8->SetLineStyle(7);
		gr8->SetLineWidth(2);
		gr8->Draw("l"); // Dessine le graph



	}

	if(choix==2) // Absolu rates
	{
	      char name[150];
	      sprintf(name,  " Variation of the beam intensity carbon ion");

	      TH1 *frame = new TH1F("frame","",0,0,50);
	      frame->SetMinimum(0.000000001);
	      frame->SetMaximum(0.001);
	      frame->SetDirectory(0);
	      frame->SetStats(0);
	      frame->SetTitle(name);
		frame->GetXaxis()->SetTitle("Number of carbon ions per bunch");
		frame->GetXaxis()->SetTickLength(0.02);
		frame->GetXaxis()->SetLabelSize(0.03);
		frame->GetYaxis()->SetTitle("Coincidences per incident carbon ion");

	      frame->GetYaxis()->SetLabelSize(0.03);
	      frame->Draw(" ");

	      // Échelle logarithmique sur les axes y
	      c1->SetGrid();
		  c1->SetLogy();

	      Double_t RMS1=0;

		TGraph *gr = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Fortuitous,0,0); // Création du graphique


		gr->SetLineColor(1); // Options de mise en forme  du graphique
		gr->SetLineStyle(21);
		gr->SetLineWidth(2);
		gr->Draw("l");

		TGraph *gr1 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_True_gamma,0,0);//error_FWHM_BGO); // Création du graphique

		gr1->SetLineColor(4); // Options de mise en forme  du graphique
		gr1->SetLineStyle(21);
		gr1->SetLineWidth(2);
		gr1->Draw("l");


		TGraph *gr2 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Same_gamma,0,0); // Création du graphique

		gr2->SetLineColor(2); // Options de mise en forme  du graphique
		gr2->SetLineStyle(21);
		gr2->SetLineWidth(2);
		gr2->Draw("l"); // Dessine le graph


		TGraph *gr3 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Same_Particle,0,0); // Création du graphique

		gr3->SetLineColor(6); // Options de mise en forme  du graphique
		gr3->SetLineStyle(21);
		gr3->SetLineWidth(2);
		gr3->Draw("l"); // Dessine le graph

		/*TGraph *gr4 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_True_Other_gamma,0,0); // Création du graphique

		gr4->SetLineColor(8); // Options de mise en forme  du graphique
		gr4->SetLineStyle(21);
		gr4->SetLineWidth(2);
		gr4->Draw("l"); // Dessine le graph*/




/*

		TGraph *gr5 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Fortuitous_histo,0,0); // Création du graphique

		gr5->SetLineColor(1); // Options de mise en forme  du graphique
		gr5->SetLineStyle(7);
		gr5->SetLineWidth(2);
		gr5->Draw("l");

		TGraph *gr6 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_True_gamma_histo,0,0); // Création du graphique

		gr6->SetLineColor(4); // Options de mise en forme  du graphique
		gr6->SetLineStyle(7);
		gr6->SetLineWidth(2);
		gr6->Draw("l");


		TGraph *gr7 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Same_gamma_histo,0,0);//error_FWHM_BGO); // Création du graphique

		gr7->SetLineColor(2); // Options de mise en forme  du graphique
		gr7->SetLineStyle(7);
		gr7->SetLineWidth(2);
		gr7->Draw("l"); // Dessine le graph


		TGraph *gr8 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Same_Particle_histo,0,0); // Création du graphique

		gr8->SetLineColor(6); // Options de mise en forme  du graphique
		gr8->SetLineStyle(7);
		gr8->SetLineWidth(2);
		gr8->Draw("l"); // Dessine le graph



	}

	/*leg = new TLegend(0.3,0.75,0.85,0.89);
	leg->SetFillColor(10);
	leg->AddEntry(gr, "Random","l");
	leg->AddEntry(gr1,"True","l");
	leg->AddEntry(gr2,"True - gamma","l");
	leg->AddEntry(gr3,"True - same particles","l");
	leg->AddEntry(gr5, "Random (reconstructed)","l");
	leg->AddEntry(gr6,"True  (reconstructed)","l");
	leg->AddEntry(gr7,"True - gamma (reconstructed)","l");
	leg->AddEntry(gr8,"True - same particles (reconstructed)","l");
	//leg->AddEntry(gr4,"True - other particles","l");
	leg->Draw();
	c1->Update();
	*/




	//----------------------------------------------------------------------------
	//-----------------------Presentation PTCOG relatif---------------------------
	//----------------------------------------------------------------------------
/*
	if(choix==3)
	{
		if(choix_histo==0)
		{
		char name[150];
		sprintf(name, " Variation of the beam intensity proton");


		TH1 *frame = new TH1F("frame","",0,0,217);
		frame->SetMinimum(0);
		frame->SetMaximum(100);
		frame->SetDirectory(0);
		frame->SetStats(0);
		frame->SetTitle(name);
		frame->GetXaxis()->SetTitle("Number of protons per bunch");
		frame->GetXaxis()->SetTickLength(0.02);
		frame->GetXaxis()->SetLabelSize(0.03);
		frame->GetYaxis()->SetTitle("Coincidences %");
		frame->GetYaxis()->SetLabelSize(0.03);
		frame->Draw(" ");

		// Échelle logarithmique sur les axes y
		c1->SetGrid();
		// c1->SetLogy();


		Double_t RMS1=0;



			if (comparaison == 1)
			{
				TGraph *gr2 = new TGraphErrors(b_sans,sigma_ns_sans_cut,Sum_bad_events_pourcentage_sans_cut,0,0); // Création du graphique

				gr2->SetLineColor(6); // Options de mise en forme  du graphique
				gr2->SetLineStyle(21);
				gr2->SetLineWidth(2);
				gr2->Draw("l");

				TGraph *gr3 = new TGraphErrors(b_sans,sigma_ns_sans_cut,datas_BGO_pourcentage_Same_gamma_sans_cut,0,0); // Création du graphique

				gr3->SetLineColor(2); // Options de mise en forme  du graphique
				gr3->SetLineStyle(21);
				gr3->SetLineWidth(2);
				gr3->Draw("l");
			}



		TGraph *gr = new TGraphErrors(b,sigma_ns,Sum_bad_events_pourcentage,0,0); // Création du graphique

		gr->SetLineColor(1); // Options de mise en forme  du graphique
		gr->SetLineStyle(21);
		gr->SetLineWidth(2);
		gr->Draw("l");

		TGraph *gr1 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_Same_gamma,0,0); // Création du graphique

		gr1->SetLineColor(4); // Options de mise en forme  du graphique
		gr1->SetLineStyle(21);
		gr1->SetLineWidth(2);
		gr1->Draw("l");

		if (comparaison == 1)
		{
				leg = new TLegend(0.3,0.75,0.85,0.89);
				leg->SetFillColor(10);
				leg->AddEntry(gr, "Background - Cut TOF","l");
				leg->AddEntry(gr1,"True gamma - Cut TOF","l");
				leg->AddEntry(gr2, "Background ","l");
				leg->AddEntry(gr3,"True gamma","l");
			    leg->Draw();
			    c1->Update();
		}

		else
		{
			leg = new TLegend(0.3,0.75,0.85,0.89);
			leg->SetFillColor(10);
			leg->AddEntry(gr, "Background - Cut TOF","l");
			leg->AddEntry(gr1,"True gamma - Cut TOF","l");

		//leg->AddEntry(gr4,"True - other particles","l");
		leg->Draw();
		c1->Update();
		}
		}

		if(choix_histo==1)
		{

			char name[150];
			sprintf(name, " Variation of the beam intensity proton");

			TH1 *frame = new TH1F("frame","",0,0,217);
			frame->SetMinimum(0);
			frame->SetMaximum(100);
			frame->SetDirectory(0);
			frame->SetStats(0);
			frame->SetTitle(name);
			frame->GetXaxis()->SetTitle("Number of protons per bunch");
			frame->GetXaxis()->SetTickLength(0.02);
			frame->GetXaxis()->SetLabelSize(0.03);
			frame->GetYaxis()->SetTitle("Coincidences %");
			frame->GetYaxis()->SetLabelSize(0.03);
			frame->Draw(" ");

			// Échelle logarithmique sur les axes y
			c1->SetGrid();
			// c1->SetLogy();


			Double_t RMS1=0;

			if (comparaison == 1)
			{
				TGraph *gr4 = new TGraphErrors(b_sans,sigma_ns_sans_cut,Sum_bad_events_pourcentage_sans_cut,0,0); // Création du graphique

				gr4->SetLineColor(6); // Options de mise en forme  du graphique
				gr4->SetLineStyle(21);
				gr4->SetLineWidth(2);
				gr4->Draw("l");

				TGraph *gr5 = new TGraphErrors(b_sans,sigma_ns_sans_cut,datas_BGO_pourcentage_Same_gamma_sans_cut,0,0); // Création du graphique

				gr5->SetLineColor(2); // Options de mise en forme  du graphique
				gr5->SetLineStyle(21);
				gr5->SetLineWidth(2);
				gr5->Draw("l");


				TGraph *gr6 = new TGraphErrors(b_sans,sigma_ns_sans_cut,Sum_bad_events_pourcentage_histo_sans_cut,0,0); // Création du graphique

				gr6->SetLineColor(6); // Options de mise en forme  du graphique
				gr6->SetLineStyle(7);
				gr6->SetLineWidth(2);
				gr6->Draw("l");

				TGraph *gr7 = new TGraphErrors(b_sans,sigma_ns_sans_cut,datas_BGO_pourcentage_Same_gamma_histo_sans_cut,0,0); // Création du graphique

				gr7->SetLineColor(2); // Options de mise en forme  du graphique
				gr7->SetLineStyle(7);
				gr7->SetLineWidth(2);
				gr7->Draw("l");
			}



				TGraph *gr = new TGraphErrors(b,sigma_ns,Sum_bad_events_pourcentage,0,0); // Création du graphique

				gr->SetLineColor(1); // Options de mise en forme  du graphique
				gr->SetLineStyle(21);
				gr->SetLineWidth(2);
				gr->Draw("l");

				TGraph *gr1 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_Same_gamma,0,0); // Création du graphique

				gr1->SetLineColor(4); // Options de mise en forme  du graphique
				gr1->SetLineStyle(21);
				gr1->SetLineWidth(2);
				gr1->Draw("l");


				TGraph *gr2 = new TGraphErrors(b,sigma_ns,Sum_bad_events_pourcentage_histo,0,0); // Création du graphique

				gr2->SetLineColor(1); // Options de mise en forme  du graphique
				gr2->SetLineStyle(7);
				gr2->SetLineWidth(2);
				gr2->Draw("l");

				TGraph *gr3 = new TGraphErrors(b,sigma_ns,datas_BGO_pourcentage_Same_gamma_histo,0,0); // Création du graphique

				gr3->SetLineColor(4); // Options de mise en forme  du graphique
				gr3->SetLineStyle(7);
				gr3->SetLineWidth(2);
				gr3->Draw("l");


			if (comparaison == 1)
			{
				leg = new TLegend(0.3,0.75,0.85,0.89);
				leg->SetFillColor(10);
				leg->AddEntry(gr, "Background - Cut TOF","l");
				leg->AddEntry(gr1,"True gamma - Cut TOF","l");
				leg->AddEntry(gr2, "Background (reconstructed) - Cut TOF","l");
				leg->AddEntry(gr3,"True gamma  (reconstructed) - Cut TOF","l");
				leg->AddEntry(gr4, "Background ","l");
				leg->AddEntry(gr5,"True gamma","l");
				leg->AddEntry(gr6, "Background (reconstructed)","l");
				leg->AddEntry(gr7,"True gamma  (reconstructed)","l");

			}
			else
			{
			leg = new TLegend(0.3,0.75,0.85,0.89);
			leg->SetFillColor(10);
			leg->AddEntry(gr, "Background - Cut TOF","l");
			leg->AddEntry(gr1,"True gamma - Cut TOF","l");
			leg->AddEntry(gr2, "Background (reconstructed) - Cut TOF","l");
			leg->AddEntry(gr3,"True gamma  (reconstructed) - Cut TOF","l");
			}
			//leg->AddEntry(gr4,"True - other particles","l");
			leg->Draw();
			c1->Update();
		}


	}
*/
	//----------------------------------------------------------------------------
	//-----------------------Presentation PTCOG Absolu---------------------------
	//----------------------------------------------------------------------------

	//if(choix==4)
	//{

		/*if(choix_histo==0)
		{
			char name[150];
			sprintf(name,  "Variation of the beam intensity carbon ion");//Variation of the beam intensity proton

			TH1 *frame = new TH1F("frame","",0,0,70);
			frame->SetMinimum(0.00001);
			frame->SetMaximum(0.01);
			frame->SetDirectory(0);
			frame->SetStats(0);
			frame->SetTitle(name);
			//	frame->	GetXaxis()->SetLabelSize(20);
			//	frame->	GetYaxis()->SetLabelSize(20);
			//frame->GetYaxis()->SetLabelSize(0.001);
			frame->GetXaxis()->SetTitleSize(0.04);
			frame->GetXaxis()->SetTitle("Number of carbon ions per bunch");//carbon ions per bunch
			frame->GetXaxis()->SetTickLength(0.02);
			frame->GetXaxis()->SetLabelSize(0.04);
			frame->GetYaxis()->SetTitleSize(0.04);
			frame->GetYaxis()->SetTitle("Coincidences per incident carbon ion");//incident carbon ion
			frame->GetYaxis()->SetLabelSize(0.04);
			frame->GetYaxis()->SetTitleOffset(1.6);
			frame->Draw(" ");

			// Échelle logarithmique sur les axes y
			c1->SetGrid();
			c1->SetLogy();

			Double_t RMS1=0;


			if(comparaison==1)
			{
				TGraph *gr = new TGraphErrors(b,sigma_ns,Sum_bad_events_absolu,0,0); // Background Cut TOF

				gr->SetLineColor(6); // Options de mise en forme  du graphique
				gr->SetLineStyle(21);
				gr->SetLineWidth(4);
				gr->Draw("l");

				TGraph *gr1 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Same_gamma,0,0); // Gamma  Cut TOF
				gr1->SetLineColor(2); // Options de mise en forme  du graphique
				gr1->SetLineStyle(21);
				gr1->SetLineWidth(4);
				gr1->Draw("l");
			}

			TGraph *gr2 = new TGraphErrors(b_sans,sigma_ns_sans_cut,Sum_bad_events_absolu_sans_cut,0,0); // Background NO TOF

			gr2->SetLineColor(1); // Options de mise en forme  du graphique
			gr2->SetLineStyle(21);
			gr2->SetLineWidth(4);
			gr2->Draw("l");

			TGraph *gr3 = new TGraphErrors(b_sans,sigma_ns_sans_cut,datas_BGO_absolu_Same_gamma_sans_cut,0,0); //  Gamma  NO TOF

			gr3->SetLineColor(4); // Options de mise en forme  du graphique
			gr3->SetLineStyle(21);
			gr3->SetLineWidth(4);
			gr3->Draw("l");





			if (comparaison == 1)
			{
				leg = new TLegend(0.3,0.75,0.85,0.89);
				leg->SetFillColor(10);
				leg->AddEntry(gr, "Background - Cut TOF","l");
				leg->AddEntry(gr1,"True gamma - Cut TOF","l");
				leg->AddEntry(gr2,"Background","l");
				leg->AddEntry(gr3,"True gamma","l");
			    leg->Draw();
			    c1->Update();
			}

			else
			{
				leg = new TLegend(0.3,0.75,0.85,0.89);
				leg->SetFillColor(10);
				leg->AddEntry(gr2, "Background","l");
				leg->AddEntry(gr3,"True gamma","l");

				//leg->AddEntry(gr4,"True - other particles","l");
				leg->Draw();
				c1->Update();
			}
		}

*/
		//if( choix_histo==1)
		//{
			char name[150];
			sprintf(name,  "");//Variation of the beam intensity proton");//carbon ion

			TH1 *frame = new TH1F("frame","",0,0.1,70);//217
			frame->SetMinimum(0.00001);//   0.000002
			frame->SetMaximum(0.01);// 0.01 0.0005
			frame->SetDirectory(0);
			frame->SetStats(0);
			frame->SetTitle(name);
		//	frame->	GetXaxis()->SetLabelSize(20);
		//	frame->	GetYaxis()->SetLabelSize(20);
			//frame->GetYaxis()->SetLabelSize(0.001);
			frame->GetXaxis()->SetTitleSize(0.04);
			frame->GetXaxis()->SetTitle("Number of protons per bunch");//carbon ions
			frame->GetXaxis()->SetTickLength(0.02);
			frame->GetXaxis()->SetLabelSize(0.04);
			frame->GetYaxis()->SetTitleSize(0.04);
			frame->GetYaxis()->SetTitle("Coincidences per incident proton");//carbon ion
			frame->GetYaxis()->SetLabelSize(0.04);
			frame->GetYaxis()->SetTitleOffset(1.6);

			frame->Draw(" ");

			// Échelle logarithmique sur les axes y
			c1->SetGrid();
			c1->SetLogy();
            c1->SetLogx();

			Double_t RMS1=0;
			cout<<"1 !!!"<<endl;

			//if(comparaison ==1)
			//{

				cout<<"3 !!!"<<endl;
				TGraph *gr = new TGraphErrors(b,sigma_ns,Sum_bad_events_absolu,0,0); // Background cut TOF

				gr->SetLineColor(6); // Options de mise en forme  du graphique
				gr->SetLineStyle(8);
				gr->SetLineWidth(2);
                gr->SetMarkerSize(2);
                gr->SetMarkerStyle(24);
                gr->SetMarkerColor(2);
				gr->Draw("p");

				TGraph *gr1 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Same_gamma,0,0); // Gamma  cut TOF

				//gr1->SetLineColor(2); // Options de mise en forme  du graphique
				//gr1->SetLineStyle(9);
				//gr1->SetLineWidth(2);
                gr1->SetMarkerSize(2);
                gr1->SetMarkerStyle(25);
                gr1->SetMarkerColor(4);
				gr1->Draw("p");

				TGraph *gr2 = new TGraphErrors(b,sigma_ns,Sum_bad_events_absolu_histo,0,0); // Background cut TOF - reconstruit

				 /*gr2->SetLineColor(8); // Options de mise en forme  du graphique
				 gr2->SetLineStyle(7);
				 gr2->SetLineWidth(2);*/
                 gr2->SetMarkerSize(2);
                 gr2->SetMarkerStyle(26);
                 gr2->SetMarkerColor(2);
				 gr2->Draw("P");

				 TGraph *gr3 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Same_gamma_histo,0,0); // Gamma cut TOF - reconstruit

				 /*gr3->SetLineColor(2); // Options de mise en forme  du graphique
				 gr3->SetLineStyle(9);
				 gr3->SetLineWidth(2);*/
                 gr3->SetMarkerSize(2);
                 gr3->SetMarkerStyle(28);
                 gr3->SetMarkerColor(4);
				 gr3->Draw("P");

			//}



			cout<<"2 !!!"<<endl;
			TGraph *gr4 = new TGraphErrors(b_sans,sigma_ns_sans_cut,Sum_bad_events_absolu_sans_cut,0,0); // Background NO cut TOF

			/*gr4->SetLineColor(1); // Options de mise en forme  du graphique
			gr4->SetLineStyle(3);
			gr4->SetLineWidth(2);*/
            gr4->SetMarkerSize(1);
            gr4->SetMarkerStyle(20);
            gr4->SetMarkerColor(2);
			gr4->Draw("P");

			TGraph *gr5 = new TGraphErrors(b_sans,sigma_ns_sans_cut,datas_BGO_absolu_Same_gamma_sans_cut,0,0); // Gamma NO cut TOF

			gr5->SetLineColor(4); // Options de mise en forme  du graphique
			gr5->SetLineStyle(5);
			gr5->SetLineWidth(2);
            gr5->SetMarkerSize(1);
            gr5->SetMarkerStyle(21);//21
            gr5->SetMarkerColor(4);
			gr5->Draw("P");

			TGraph *gr6 = new TGraphErrors(b_sans,sigma_ns_sans_cut,Sum_bad_events_absolu_histo_sans_cut,0,0); //Background NO cut TOF - reconstruit

			 gr6->SetLineColor(1); // Options de mise en forme  du graphique
			 gr6->SetLineStyle(10);
			 gr6->SetLineWidth(2);
             gr6->SetMarkerSize(1);
             gr6->SetMarkerStyle(22);
             gr6->SetMarkerColor(2);
			 gr6->Draw("p");

			 TGraph *gr7 = new TGraphErrors(b_sans,sigma_ns_sans_cut,datas_BGO_absolu_Same_gamma_histo_sans_cut,0,0); // Gamma NO cut TOF - reconstruit

			 gr7->SetLineColor(4); // Options de mise en forme  du graphique
			 gr7->SetLineStyle(4);
             gr7->SetLineWidth(2);
             gr7->SetMarkerSize(1);
             gr7->SetMarkerStyle(34);
             gr7->SetMarkerColor(4);
			 gr7->Draw("p");



			//if (comparaison == 1)
			//{
				leg = new TLegend(0.12,0.76,0.7,0.9);
				leg->SetFillColor(10);
				leg->AddEntry(gr4, "Background ","p");
				leg->AddEntry(gr5,"True gamma","p");
				leg->AddEntry(gr6, "Background (reconstructed)","p");
				leg->AddEntry(gr7,"True gamma  (reconstructed)","p");

                /*leg->AddEntry(gr, "Background - Cut TOF","p");
				leg->AddEntry(gr1,"True gamma - Cut TOF","p");
				leg->AddEntry(gr2, "Background (reconstructed) - Cut TOF","P");
				leg->AddEntry(gr3,"True gamma  (reconstructed) - Cut TOF","P");
				leg->AddEntry(gr4, "Background ","p");
				leg->AddEntry(gr5,"True gamma","p");
				leg->AddEntry(gr6, "Background (reconstructed)","p");
				leg->AddEntry(gr7,"True gamma  (reconstructed)","p");*/
				cout<<"4 !!!"<<endl;
			//}
			/*else
			{
				leg = new TLegend(0.3,0.75,0.85,0.89);
				leg->SetFillColor(10);
				leg->AddEntry(gr4, "Background ","p");
				leg->AddEntry(gr5,"True gamma ","p");
				leg->AddEntry(gr6, "Background (reconstructed) ","p");
				leg->AddEntry(gr7,"True gamma  (reconstructed) ","p");
				cout<<"5 !!!"<<endl;
			}*/
			//leg->AddEntry(gr4,"True - other particles","l");
			leg->Draw();
			c1->Update();
		//}


//	}







}


//-------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------- 30/05/2014 - PTCOG 53 ----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------



/*
//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------- PROTONS -------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------


if (choix_p_c ==1)
{


	//-------------------
	//-------------------
	//Coinc Pourcentage
	//-------------------
	//-------------------



	//Variation intensity
	Double_t datas_BGO_pourcentage_Fortuitous[b]={0,0.0684697,0.849473,3.83598,5.05689,10.3056,16.1336,21.2097,24.9559,27.5146,31.6932,34.9129,36.7713,39.108,49.2698,56.6647,65.0378,71.3597,75.23,80.1723,84.6409,88.3753,90.229,90.7323};
	Double_t datas_BGO_pourcentage_True_gamma[b]	={100,99.9315,99.1505,96.164,94.9431,89.6944,83.8664,78.7903,75.0441,72.4854,68.3068,65.0871,63.2287,60.892,50.7302,43.3353,34.9622,28.6403,24.77,19.8277,15.3591,11.6247,9.77101,9.2677};
	Double_t datas_BGO_pourcentage_True_Other_gamma[b]= {1.0274,1.02705,1.01937,0.992063,4.83565,4.52367,4.23681,3.96197,3.88203,3.75977,3.4419,3.35023,3.30984,3.22114,2.62517,2.18637,1.75063,1.42398,1.34661,0.990179,0.873984,0.634939,0.47799,0.423007};
	Double_t datas_BGO_pourcentage_Same_gamma[b]	={73.6301,73.5707,72.9867,70.8333,51.359,48.4721,45.2301,42.3138,40.1311,38.7207,36.3825,34.8468,33.6109,32.4592,26.8255,22.7232,18.0605,14.6788,12.3566,9.86154,7.58808,5.69665,4.72988,4.62561};
	Double_t datas_BGO_pourcentage_Same_Particle[b]	={24.726,24.7176,24.5328,23.7434,30.5942,28.9694,27.1324,25.6999,24.5021,23.877,22.638,21.2916,20.884,19.9257,17.0202,14.7165,12.1788,10.0214,8.92366,7.16471,5.50813,4.19535,3.59604,3.32363};

	Double_t Sum_bad_events_pourcentage[b];

	for(int u=0;u<b;u++)
	{
		Sum_bad_events_pourcentage[u]= datas_BGO_pourcentage_Fortuitous[u]+datas_BGO_pourcentage_True_gamma[u]-datas_BGO_pourcentage_Same_gamma[u];
	}

	//--------------------------------
	//Reconstructed events (line-cone)
	//--------------------------------


	Double_t datas_BGO_pourcentage_total_histo[b]={1782,1852,1943,2010,2073,2128,2236,2285,2337,2402,2782,3092,3555,4034,4526,5224,6117,6895,7226,7358};
	Double_t datas_BGO_pourcentage_Fortuitous_histo[b]={67,145,240,310,385,443,556,603,657,731,1131,1468,1992,2548,3077,3885,4865,5819,6259,6446} ;
	Double_t datas_BGO_pourcentage_True_gamma_histo[b]	={1715,1707,1703,1700,1688,1685,1680,1682,1680,1671,1651,1624,1563,1486,1449,1339,1252,1076,967,912};
	Double_t datas_BGO_pourcentage_True_Other_gamma_histo[b]= {64,64,64,64,63,64,62,62,64,61,63,60,56,54,58,53,57,51,42,34};
	Double_t datas_BGO_pourcentage_Same_gamma_histo[b]	={1095,1090,1086,1082,1075,1072,1065,1071,1063,1064,1050,1030,979,924,884,821,759,653,569,571};
	Double_t datas_BGO_pourcentage_Same_Particle_histo[b]	={445,442,443,444,440,441,444,440,444,437,430,426,420,401,408,366,344,292,278,236};

	//pour ptcog, non prise en compte de la fenetre en temps -20 /+20 ns pour les protons : perte de certains évènements

	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_pourcentage_Fortuitous_histo[ff] = datas_BGO_pourcentage_Fortuitous_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_True_gamma_histo[ff]= datas_BGO_pourcentage_True_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_True_Other_gamma_histo[ff]= datas_BGO_pourcentage_True_Other_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_Same_gamma_histo[ff] = datas_BGO_pourcentage_Same_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_Same_Particle_histo[ff] = datas_BGO_pourcentage_Same_Particle_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;

	}



	Double_t Sum_bad_events_pourcentage_histo[b];

	for(int u=0;u<b;u++)
	{
		Sum_bad_events_pourcentage_histo[u]= datas_BGO_pourcentage_Fortuitous_histo[u]+datas_BGO_pourcentage_True_gamma_histo[u]-datas_BGO_pourcentage_Same_gamma_histo[u];
	}


	//-------------------
	//-------------------
	//Coinc Absolu
	//-------------------
	//-------------------

	//Variation intensity
	Double_t datas_BGO_absolu_Fortuitous[b]={1.6e-06,3.44e-06,5.75e-06,8.03e-06,9.9e-06,1.127e-05,1.372e-05,1.584e-05,1.722e-05,1.894e-05,2.834e-05,3.758e-05,5.164e-05,6.665e-05,7.933e-05,9.959e-05,0.00012493,0.00014893,0.00016234,0.00016516};
	Double_t datas_BGO_absolu_True_gamma[b]	={3.004e-05,2.994e-05,2.989e-05,2.983e-05,2.977e-05,2.969e-05,2.957e-05,2.953e-05,2.961e-05,2.949e-05,2.918e-05,2.874e-05,2.776e-05,2.675e-05,2.612e-05,2.463e-05,2.267e-05,1.959e-05,1.758e-05,1.687e-05};
	Double_t datas_BGO_absolu_True_Other_gamma[b]= {1.53e-06,1.51e-06,1.51e-06,1.5e-06,1.54e-06,1.54e-06,1.49e-06,1.52e-06,1.55e-06,1.56e-06,1.51e-06,1.45e-06,1.39e-06,1.33e-06,1.42e-06,1.23e-06,1.29e-06,1.07e-06,8.6e-07,7.7e-07};
	Double_t datas_BGO_absolu_Same_gamma[b]	={1.625e-05,1.618e-05,1.612e-05,1.602e-05,1.592e-05,1.586e-05,1.575e-05,1.581e-05,1.574e-05,1.572e-05,1.543e-05,1.507e-05,1.434e-05,1.371e-05,1.303e-05,1.225e-05,1.12e-05,9.6e-06,8.51e-06,8.42e-06};
	Double_t datas_BGO_absolu_Same_Particle[b]	={9.68e-06,9.67e-06,9.67e-06,9.73e-06,9.72e-06,9.78e-06,9.8e-06,9.66e-06,9.78e-06,9.65e-06,9.79e-06,9.76e-06,9.67e-06,9.36e-06,9.41e-06,8.9e-06,8.13e-06,7.07e-06,6.47e-06,6.05e-06};

	Double_t Sum_bad_events_absolu[b];

	for(int u=0;u<b;u++)
	{
		Sum_bad_events_absolu[u]= datas_BGO_absolu_Fortuitous[u]+datas_BGO_absolu_True_gamma[u]-datas_BGO_absolu_Same_gamma[u];
	}



	//--------------------------------
	//Reconstructed events (line-cone)
	//--------------------------------


	Double_t datas_BGO_absolu_total_histo[b]={1782,1852,1943,2010,2073,2128,2236,2285,2337,2402,2782,3092,3555,4034,4526,5224,6117,6895,7226,7358};
	Double_t datas_BGO_absolu_Fortuitous_histo[b]={67,145,240,310,385,443,556,603,657,731,1131,1468,1992,2548,3077,3885,4865,5819,6259,6446} ;
	Double_t datas_BGO_absolu_True_gamma_histo[b]	={1715,1707,1703,1700,1688,1685,1680,1682,1680,1671,1651,1624,1563,1486,1449,1339,1252,1076,967,912};
	Double_t datas_BGO_absolu_True_Other_gamma_histo[b]= {64,64,64,64,63,64,62,62,64,61,63,60,56,54,58,53,57,51,42,34};
	Double_t datas_BGO_absolu_Same_gamma_histo[b]	={1095,1090,1086,1082,1075,1072,1065,1071,1063,1064,1050,1030,979,924,884,821,759,653,569,571};
	Double_t datas_BGO_absolu_Same_Particle_histo[b]	={445,442,443,444,440,441,444,440,444,437,430,426,420,401,408,366,344,292,278,236};


	/*
	 for (int ff= 0; ff<b; ff++)
	 {
	 datas_BGO_absolu_Fortuitous[ff] = datas_BGO_absolu_Fortuitous[ff]/ (100);
	 datas_BGO_absolu_True_gamma[ff]= datas_BGO_absolu_True_gamma[ff]/ (100);
	 datas_BGO_absolu_True_Other_gamma[ff]= datas_BGO_absolu_True_Other_gamma[ff]/ (100);
	 datas_BGO_absolu_Same_gamma[ff] = datas_BGO_absolu_Same_gamma[ff]/ (100);
	 datas_BGO_absolu_Same_Particle[ff] = datas_BGO_absolu_Same_Particle[ff]/ (100);

	 }
	 */

/*

	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_absolu_Fortuitous_histo[ff] = datas_BGO_absolu_Fortuitous_histo[ff]/ (1e8);
		datas_BGO_absolu_True_gamma_histo[ff]= datas_BGO_absolu_True_gamma_histo[ff]/ (1e8);
		datas_BGO_absolu_True_Other_gamma_histo[ff]= datas_BGO_absolu_True_Other_gamma_histo[ff]/ (1e8);
		datas_BGO_absolu_Same_gamma_histo[ff] = datas_BGO_absolu_Same_gamma_histo[ff]/ (1e8);
		datas_BGO_absolu_Same_Particle_histo[ff] = datas_BGO_absolu_Same_Particle_histo[ff]/ (1e8);

	}


	Double_t Sum_bad_events_absolu_histo[b];

	for(int u=0;u<b;u++)
	{
		Sum_bad_events_absolu_histo[u]= datas_BGO_absolu_Fortuitous_histo[u]+datas_BGO_absolu_True_gamma_histo[u]-datas_BGO_absolu_Same_gamma_histo[u];
	}



}

//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------- IONS CARBONE -------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------


if (choix_p_c ==2)
{


	//-------------------
	//-------------------
	//Coinc Pourcentage
	//-------------------
	//-------------------




	//Variation intensity
	Double_t datas_BGO_pourcentage_Fortuitous[b]={0.221047,6.03455,17.1383,24.2244,31.2202,36.0823,41.2525,45.0256,48.7508,51.6842,62.6233,69.2741,76.938,81.4909};
	Double_t datas_BGO_pourcentage_True_gamma[b]	={99.779,93.9655,82.8617,75.7756,68.7798,63.9177,58.7475,54.9744,51.2492,48.3158,37.3767,30.7259,23.062,18.5091};
	Double_t datas_BGO_pourcentage_True_Other_gamma[b]= {25.4164,23.9886,21.1069,19.3387,17.5979,16.3316,15.0152,14.1242,13.066,12.3565,9.52952,7.93349,5.96452,4.825};
	Double_t datas_BGO_pourcentage_Same_gamma[b]	={22.8626,21.4873,18.8735,17.191,15.3764,14.3406,13.1057,12.2029,11.3431,10.5757,8.07621,6.56742,4.87345,3.8659};
	Double_t datas_BGO_pourcentage_Same_Particle[b]	={38.2845,36.0689,31.8433,29.264,26.6698,24.7185,22.8364,21.4127,20.0344,18.9489,14.8651,12.2026,9.22022,7.39531};

	Double_t Sum_bad_events_pourcentage[b];

	for(int u=0;u<b;u++)
	{
		Sum_bad_events_pourcentage[u]= datas_BGO_pourcentage_Fortuitous[u]+datas_BGO_pourcentage_True_gamma[u]-datas_BGO_pourcentage_Same_gamma[u];
	}

	//--------------------------------
	//Reconstructed events (line-cone)
	//--------------------------------

	Double_t datas_BGO_pourcentage_total_histo[b]={12547,12973,13701,14464,15153,15790,16431,17082,17748,18474,21204,23600,27344,29454};
	Double_t datas_BGO_pourcentage_Fortuitous_histo[b]={26,568,1675,2638,3640,4455,5408,6204,7122,8023,11591,14693,19504,22647} ;
	Double_t datas_BGO_pourcentage_True_gamma_histo[b]	={12521,12405,12026,11826,11513,11335,11023,10878,10626,10451,9613,8907,7840,6807};
	Double_t datas_BGO_pourcentage_True_Other_gamma_histo[b]= {2375,2360,2296,2263,2197,2175,2139,2121,2073,2035,1891,1780,1580,1420};
	Double_t datas_BGO_pourcentage_Same_gamma_histo[b]	={3733,3703,3586,3530,3410,3362,3255,3211,3134,3071,2820,2608,2297,1994};
	Double_t datas_BGO_pourcentage_Same_Particle_histo[b]	={5230,5164,5000,4915,4811,4709,4573,4524,4394,4352,3988,3663,3210,2693};



	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_pourcentage_Fortuitous_histo[ff] = datas_BGO_pourcentage_Fortuitous_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_True_gamma_histo[ff]= datas_BGO_pourcentage_True_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_True_Other_gamma_histo[ff]= datas_BGO_pourcentage_True_Other_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_Same_gamma_histo[ff] = datas_BGO_pourcentage_Same_gamma_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;
		datas_BGO_pourcentage_Same_Particle_histo[ff] = datas_BGO_pourcentage_Same_Particle_histo[ff]/ datas_BGO_pourcentage_total_histo[ff]*100;

	}


	Double_t Sum_bad_events_pourcentage_histo[b];

	for(int u=0;u<b;u++)
	{
		Sum_bad_events_pourcentage_histo[u]= datas_BGO_pourcentage_Fortuitous_histo[u]+datas_BGO_pourcentage_True_gamma_histo[u]-datas_BGO_pourcentage_Same_gamma_histo[u];
		//	cout<<"test "<<Sum_bad_events_pourcentage_histo[u]<<endl;
	}

	//-------------------
	//-------------------
	//Coinc Absolu
	//-------------------
	//-------------------


	//Variation intensity
	Double_t datas_BGO_absolu_Fortuitous[b]={5.6e-07,1.614e-05,5.057e-05,7.715e-05,0.00010751,0.00013139,0.00015987,0.00018397,0.00020995,0.00023277,0.00033955,0.00042952,0.00057118,0.00066358};
	Double_t datas_BGO_absolu_True_gamma[b]	={0.00025278,0.00025132,0.0002445,0.00024133,0.00023685,0.00023275,0.00022767,0.00022462,0.00022071,0.0002176,0.00020266,0.00019051,0.00017121,0.00015072};
	Double_t datas_BGO_absolu_True_Other_gamma[b]= {6.439e-05,6.416e-05,6.228e-05,6.159e-05,6.06e-05,5.947e-05,5.819e-05,5.771e-05,5.627e-05,5.565e-05,5.167e-05,4.919e-05,4.428e-05,3.929e-05};
	Double_t datas_BGO_absolu_Same_gamma[b]	={5.792e-05,5.747e-05,5.569e-05,5.475e-05,5.295e-05,5.222e-05,5.079e-05,4.986e-05,4.885e-05,4.763e-05,4.379e-05,4.072e-05,3.618e-05,3.148e-05};
	Double_t datas_BGO_absolu_Same_Particle[b]	={9.699e-05,9.647e-05,9.396e-05,9.32e-05,9.184e-05,9.001e-05,8.85e-05,8.749e-05,8.628e-05,8.534e-05,8.06e-05,7.566e-05,6.845e-05,6.022e-05};

	Double_t Sum_bad_events_absolu[b];

	for(int u=0;u<b;u++)
	{
		Sum_bad_events_absolu[u]= datas_BGO_absolu_Fortuitous[u]+datas_BGO_absolu_True_gamma[u]-datas_BGO_absolu_Same_gamma[u];
	}


	//--------------------------------
	//Reconstructed events (line-cone)
	//--------------------------------


	Double_t datas_BGO_absolu_total_histo[b]={12547,12973,13701,14464,15153,15790,16431,17082,17748,18474,21204,23600,27344,29454};
	Double_t datas_BGO_absolu_Fortuitous_histo[b]={26,568,1675,2638,3640,4455,5408,6204,7122,8023,11591,14693,19504,22647} ;
	Double_t datas_BGO_absolu_True_gamma_histo[b]	={12521,12405,12026,11826,11513,11335,11023,10878,10626,10451,9613,8907,7840,6807};
	Double_t datas_BGO_absolu_True_Other_gamma_histo[b]= {2375,2360,2296,2263,2197,2175,2139,2121,2073,2035,1891,1780,1580,1420};
	Double_t datas_BGO_absolu_Same_gamma_histo[b]	={3733,3703,3586,3530,3410,3362,3255,3211,3134,3071,2820,2608,2297,1994};
	Double_t datas_BGO_absolu_Same_Particle_histo[b]	={5230,5164,5000,4915,4811,4709,4573,4524,4394,4352,3988,3663,3210,2693};





	/*
	 for (int ff= 0; ff<b; ff++)
	 {
	 datas_BGO_absolu_Fortuitous[ff] = datas_BGO_absolu_Fortuitous[ff]/ (100);
	 datas_BGO_absolu_True_gamma[ff]= datas_BGO_absolu_True_gamma[ff]/ (100);
	 datas_BGO_absolu_True_Other_gamma[ff]= datas_BGO_absolu_True_Other_gamma[ff]/ (100);
	 datas_BGO_absolu_Same_gamma[ff] = datas_BGO_absolu_Same_gamma[ff]/ (100);
	 datas_BGO_absolu_Same_Particle[ff] = datas_BGO_absolu_Same_Particle[ff]/ (100);

	 }
	 */

/*

	for (int ff= 0; ff<b; ff++)
	{
		datas_BGO_absolu_Fortuitous_histo[ff] = datas_BGO_absolu_Fortuitous_histo[ff]/ (1e8);
		datas_BGO_absolu_True_gamma_histo[ff]= datas_BGO_absolu_True_gamma_histo[ff]/ (1e8);
		datas_BGO_absolu_True_Other_gamma_histo[ff]= datas_BGO_absolu_True_Other_gamma_histo[ff]/ (1e8);
		datas_BGO_absolu_Same_gamma_histo[ff] = datas_BGO_absolu_Same_gamma_histo[ff]/ (1e8);
		datas_BGO_absolu_Same_Particle_histo[ff] = datas_BGO_absolu_Same_Particle_histo[ff]/ (1e8);

	}


	Double_t Sum_bad_events_absolu_histo[b];
	for(int u=0;u<b;u++)
	{
		Sum_bad_events_absolu_histo[u]= datas_BGO_absolu_Fortuitous_histo[u]+datas_BGO_absolu_True_gamma_histo[u]-datas_BGO_absolu_Same_gamma_histo[u];
	}



}
*/


//--------------------------------------------------------------------
//--------------------------------------------------------------------
//--------------------------------------------------------------------



//Datas de décembre 2013

/*
//Coinc Pourcentage

//Timing resolutions : Si = 15 ns, BGO = 3 ns
Double_t datas_BGO_pourcentage_Fortuitous[b]={95.94,95.69,95.17,94.96,94.93,94.92,94.92,94.92,94.92};
Double_t datas_BGO_pourcentage_True_gamma[b]	={0.94,1.93,2.47,2.67,2.69,2.69,2.69,2.69,2.69};
Double_t datas_BGO_pourcentage_Same_gamma[b]	={0.52,0.62,0.64,0.64,0.63,0.63,0.63,0.63,0.63};
Double_t datas_BGO_pourcentage_Same_Particle[b]	={2.59,1.74,1.7,1.74,1.75,1.75,1.75,1.76,1.76};

//Timing resolutions : Si = 1 ns, BGO = 1 ns
/*Double_t datas_Si_pourcentage_Fortuitous [b]	={96.28,90.33,80.12,78.23,77.5,77.3,78.79,78.47,78.4,78.43,78.43} ;
 Double_t datas_Si_pourcentage_True_gamma[b]	={2.08,6.69,15.69,17.22,17.6,,17.65,16.0,16.24,16.32,16.33,16.33};
 Double_t datas_Si_pourcentage_Same_gamma[b]	={1.18,1.99,2.75,2.77,2.81,2.82,2.53,2.49,2.49,2.49,2.49};
 Double_t datas_Si_pourcentage_Same_Particle[b]	={0.45,2.17,1.42,1.77,2.05,2.19,2.65,2.78,2.77,2.75,2.74};*/

/*

//Coinc absolu

//Timing resolutions : Si = 15 ns, BGO = 3 ns
Double_t datas_BGO_absolu_Fortuitous[b]={573,830,861,845,835,833, 832,831,831};
Double_t datas_BGO_absolu_True_gamma[b]	={5.6,16.77,22.4,23.77,23.67,23.67,23.67,23.67,23.67};
Double_t datas_BGO_absolu_Same_gamma[b]	={3.14,5.38,5.8,5.59,5.52,5.52,5.49,5.5,5.5};
Double_t datas_BGO_absolu_Same_Particle[b]	={15.5,15.16,15.4,15.49,15.44,15.39,15.38,15.38,15.38};

//Timing resolutions : Si = 1 ns, BGO = 1 ns
/*Double_t datas_Si_absolu_Fortuitous [b]	={181,203,191,181,175,174,190,183,180,180,180} ;
 Double_t datas_Si_absolu_True_gamma[b]	={3.9,14.96,37.57,39.9,39.75,39.68,38.69,37.84,37.59,37.52,37.49};
 Double_t datas_Si_absolu_Same_gamma[b]	={2.21,4.95,6.59,6.43,6.35,6.35,6.1,5.81,5.74,5.71,5.71};
 Double_t datas_Si_absolu_Same_Particle[b]	={0.86,1.76,3.42,4.1,4.64,4.93,6.4,6.48,6.39,6.31,6.3};*/


//Double_t sigma_ns [b]	 = {0.5,1,2,3,4,5,10,20,30,40,50};//Timing resolutions : Si = 1 ns, BGO = 1 ns

//Double_t sigma_ns [b]	 = {1,5,10,20,30,40,50,70,100};//Timing resolutions : Si = 15 ns, BGO = 3 ns, b=9
