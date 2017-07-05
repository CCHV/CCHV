
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

#include <math>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>





//----------------------------------------------------------------------------------------------------------------------------
// Creation d'un histogramme pour tracer la FWHM de reconstruction en fonction de la resolution en energie de l'absorbeur NaI
//----------------------------------------------------------------------------------------------------------------------------

void graphique_2D_variation_intensity_Cut_TOF_High_stat()
{

TCanvas *c1 = new TCanvas("c1","c1",700,800); // Création du canvas

    const int b =23;//23 protons , 23 carbon // number of values
	
	const int b_sans =23; // sans cut TOF

    Int_t choix = 4; // choix = 1 : %;  choix = 2 : absolu ; choix = 3 : % PTCOG; choix = 4 : absolu PTCOG
	
	Int_t choix_p_c = 1; // choix_p_c = 1 : protons, choix_p_c = 2 : carbone
	
	Int_t choix_histo = 1; //choisit de mettre l'histo sur la courbe ou pas: choix_histo =0: pas d'histo, choix_histo =1: pas d'histo, 
	
	Int_t comparaison = 1; // si on veut comparer avec ou sans cut TOF


//-----------------	data 05/01/2015 - Prise en compte fenêtre coïnc [-20;20] ns et CUT TOF [0;6] ns + distribution Poisson ------------
	
	//Protons
	
Double_t sigma_ns [b]	 = {0.001,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,80,100,150,200,217};//b=23
	
	//Carbon
	//Double_t sigma_ns [b]	 = {0.0001,0.001,0.005,0.01,0.05,0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,60,70};//b =23
	
	
	//-----------------	data data 05/01/2015 - Prise en compte fenêtre coïnc [-20;20] ns et NO cut TOF + distribution Poisson ------------
	
	//Protons
	Double_t sigma_ns_sans_cut [b_sans]	 = {0.001,0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,80,100,150,200,217};//b=23
	
	//Carbon
	//Double_t sigma_ns_sans_cut [b_sans]	 = {0.0001,0.001,0.005,0.01,0.05,0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,60,70};//b =23
	

	
	
	
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

	
	
	if (choix_p_c ==1)
	{
	
	
	//-------------------
	//-------------------
	//Coinc Pourcentage
	//-------------------
	//-------------------	
	

	
	//Variation intensity
	Double_t datas_BGO_pourcentage_Fortuitous[b]={0.530786,0.634473,1.1036,5.32423,9.78155,16.9532,24.0439,29.3809,33.959,38.5106,41.8584,45.1866,47.7008,50.6222,60.3779,66.457,74.8971,83.2051,88.0096,90.0316,92.9967,94.565,94.9203};
	Double_t datas_BGO_pourcentage_True_gamma[b]	={99.4692,99.3655,98.8964,94.6758,90.2185,83.0468,75.9561,70.6191,66.041,61.4894,58.1416,54.8134,52.2992,49.3778,39.6221,33.543,25.1029,16.7949,11.9904,9.9684,7.00332,5.43505,5.07969};
	Double_t datas_BGO_pourcentage_True_Other_gamma[b]= {34.2534,34.4025,34.1047,32.3891,30.8445,28.6878,25.9561,23.9769,22.697,21.1797,19.8334,18.9445,18.0186,17.1192,13.6499,11.1856,8.30088,5.28064,3.62414,3.12586,2.02113,1.55114,1.4722};
	Double_t datas_BGO_pourcentage_Same_gamma[b]	={56.2987,55.9746,55.6426,53.2082,51.0923,46.6968,43.1384,39.9528,37.3179,34.6618,32.9389,30.9008,29.4849,27.9789,22.4895,19.1987,14.432,9.90558,7.16791,5.92196,4.22042,3.33253,3.0904};
	Double_t datas_BGO_pourcentage_Same_Particle[b]	={8.9172,8.98837,9.14916,9.0785,8.28171,7.66214,6.86164,6.6894,6.02618,5.64792,5.36927,4.96805,4.79574,4.27979,3.48274,3.15876,2.37009,1.60867,1.1983,0.920583,0.761764,0.551381,0.517095};
	
		Double_t Sum_bad_events_pourcentage[b];
		
		for(int u=0;u<b;u++)
		{
			Sum_bad_events_pourcentage[u]= datas_BGO_pourcentage_Fortuitous[u]+datas_BGO_pourcentage_True_gamma[u]-datas_BGO_pourcentage_Same_gamma[u];
		}
	
	//--------------------------------
	//Reconstructed events (line-cone)
	//--------------------------------
	

	Double_t datas_BGO_pourcentage_total_histo[b]={1604,1616,1596,1636,1702,1790,1863,1968,2059,2190,2252,2350,2443,2514,2961,3323,3888,4821,5615,5924,6481,6537,6529};
	Double_t datas_BGO_pourcentage_Fortuitous_histo[b]={7,8,14,64,135,226,331,438,541,660,749,844,949,1038,1521,1932,2627,3759,4707,5132,5873,6061,6087} ;
	Double_t datas_BGO_pourcentage_True_gamma_histo[b]	={1597,1608,1582,1572,1567,1564,1532,1530,1518,1530,1503,1506,1494,1476,1440,1391,1261,1062,908,792,608,476,442};
	Double_t datas_BGO_pourcentage_True_Other_gamma_histo[b]= {433,442,437,424,427,435,420,413,420,421,408,421,417,408,396,376,337,273,231,200,143,111,104};
	Double_t datas_BGO_pourcentage_Same_gamma_histo[b]	={1037,1043,1017,1016,1017,1009,996,991,983,986,979,970,959,953,931,904,822,708,609,530,417,331,312};
	Double_t datas_BGO_pourcentage_Same_Particle_histo[b]	={127,123,128,132,123,120,116,126,115,123,116,115,118,115,113,111,102,81,68,62,48,34,26};
	
		
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
	Double_t datas_BGO_absolu_Fortuitous[b]={1.5e-07,1.8e-07,3.1e-07,1.56e-06,3e-06,5.62e-06,8.55e-06,1.12e-05,1.375e-05,1.691e-05,1.91e-05,2.192e-05,2.417e-05,2.685e-05,3.866e-05,4.86e-05,6.731e-05,9.517e-05,0.00012045,0.00013105,0.00015138,0.00015607,0.00015603};
	Double_t datas_BGO_absolu_True_gamma[b]	={2.811e-05,2.819e-05,2.778e-05,2.774e-05,2.767e-05,2.753e-05,2.701e-05,2.692e-05,2.674e-05,2.7e-05,2.653e-05,2.659e-05,2.65e-05,2.619e-05,2.537e-05,2.453e-05,2.256e-05,1.921e-05,1.641e-05,1.451e-05,1.14e-05,8.97e-06,8.35e-06};
	Double_t datas_BGO_absolu_True_Other_gamma[b]= {9.68e-06,9.76e-06,9.58e-06,9.49e-06,9.46e-06,9.51e-06,9.23e-06,9.14e-06,9.19e-06,9.3e-06,9.05e-06,9.19e-06,9.13e-06,9.08e-06,8.74e-06,8.18e-06,7.46e-06,6.04e-06,4.96e-06,4.55e-06,3.29e-06,2.56e-06,2.42e-06};
	Double_t datas_BGO_absolu_Same_gamma[b]	={1.591e-05,1.588e-05,1.563e-05,1.559e-05,1.567e-05,1.548e-05,1.534e-05,1.523e-05,1.511e-05,1.522e-05,1.503e-05,1.499e-05,1.494e-05,1.484e-05,1.44e-05,1.404e-05,1.297e-05,1.133e-05,9.81e-06,8.62e-06,6.87e-06,5.5e-06,5.08e-06};
	Double_t datas_BGO_absolu_Same_Particle[b]	={2.52e-06,2.55e-06,2.57e-06,2.66e-06,2.54e-06,2.54e-06,2.44e-06,2.55e-06,2.44e-06,2.48e-06,2.45e-06,2.41e-06,2.43e-06,2.27e-06,2.23e-06,2.31e-06,2.13e-06,1.84e-06,1.64e-06,1.34e-06,1.24e-06,9.1e-07,8.5e-07};
file://localhost/Users/JeanLuc/Desktop/The%CC%80se/Article/Article/Clinic_Applicability_Compton_Camera/Analysis/Simulations/Timing_resolution_and_window_time/graphique_2D_variation_intensity_Cut_TOF_Distribution_poisson_Seminar2ndYear.C
		Double_t Sum_bad_events_absolu[b];
		
		for(int u=0;u<b;u++)
		{
			Sum_bad_events_absolu[u]= datas_BGO_absolu_Fortuitous[u]+datas_BGO_absolu_True_gamma[u]-datas_BGO_absolu_Same_gamma[u];
		}
		
		

	//--------------------------------
	//Reconstructed events (line-cone)
	//--------------------------------
	
	
	Double_t datas_BGO_absolu_total_histo[b]= {1604,1616,1596,1636,1702,1790,1863,1968,2059,2190,2252,2350,2443,2514,2961,3323,3888,4821,5615,5924,6481,6537,6529};
	Double_t datas_BGO_absolu_Fortuitous_histo[b]= {7,8,14,64,135,226,331,438,541,660,749,844,949,1038,1521,1932,2627,3759,4707,5132,5873,6061,6087} ;
	Double_t datas_BGO_absolu_True_gamma_histo[b]= {1597,1608,1582,1572,1567,1564,1532,1530,1518,1530,1503,1506,1494,1476,1440,1391,1261,1062,908,792,608,476,442};
	Double_t datas_BGO_absolu_True_Other_gamma_histo[b]= {433,442,437,424,427,435,420,413,420,421,408,421,417,408,396,376,337,273,231,200,143,111,104};
	Double_t datas_BGO_absolu_Same_gamma_histo[b]= {1037,1043,1017,1016,1017,1009,996,991,983,986,979,970,959,953,931,904,822,708,609,530,417,331,312};
	Double_t datas_BGO_absolu_Same_Particle_histo[b]= {127,123,128,132,123,120,116,126,115,123,116,115,118,115,113,111,102,81,68,62,48,34,26};
	
	
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
		
		
		
	}
	
	
	
	
	
	
	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- SANS CUT TOF --------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	
	
	
	
	if (choix_p_c ==1)
	{
		
		
		//-------------------
		//-------------------
		//Coinc Pourcentage
		//-------------------
		//-------------------	
		
		
		
		//Variation intensity
		Double_t datas_BGO_pourcentage_Fortuitous_sans_cut[b_sans]={0.493182,0.576868,0.987224,4.97624,9.64209,16.141,23.4682,28.5961,32.9984,37.691,41.1957,44.4785,46.8489,49.2368,59.4365,65.5635,74.1489,82.5292,87.5589,89.7777,92.7326,94.3128,94.7026};
		Double_t datas_BGO_pourcentage_True_gamma_sans_cut[b_sans]	={99.5068,99.4231,99.0128,95.0238,90.3579,83.859,76.5318,71.4039,67.0016,62.309,58.8043,55.5215,53.1511,50.7632,40.5635,34.4365,25.8511,17.4708,12.4411,10.2223,7.26738,5.68716,5.29744};
		Double_t datas_BGO_pourcentage_True_Other_gamma_sans_cut[b_sans]= {29.4749,29.3914,28.9199,27.9843,26.3622,24.8324,22.5896,20.7462,19.7703,18.4493,17.1558,16.377,15.6163,14.9646,12.0256,9.91698,7.27647,4.64181,3.23884,2.73127,1.78461,1.41413,1.31291};
		Double_t datas_BGO_pourcentage_Same_gamma_sans_cut[b_sans]	={46.4462,45.9763,45.5865,43.7517,41.6667,38.5895,35.6069,33.0817,31.05,28.7493,27.3913,25.6476,24.603,23.5405,18.86,16.0014,12.1586,8.3845,6.00012,4.96543,3.62595,2.84358,2.65126};
		Double_t datas_BGO_pourcentage_Same_Particle_sans_cut[b_sans]	={23.5857,24.0554,24.5064,23.2877,22.3291,20.437,18.3353,17.576,16.1813,15.1104,14.2572,13.4969,12.9317,12.2581,9.67784,8.51814,6.41601,4.44444,3.20211,2.52557,1.85682,1.42945,1.33327};
		
		Double_t Sum_bad_events_pourcentage_sans_cut[b_sans];
		
		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_pourcentage_sans_cut[u]= datas_BGO_pourcentage_Fortuitous_sans_cut[u]+datas_BGO_pourcentage_True_gamma_sans_cut[u]-datas_BGO_pourcentage_Same_gamma_sans_cut[u];
		}
		
		//--------------------------------
		//Reconstructed events (line-cone)
		//--------------------------------
		
		
		Double_t datas_BGO_pourcentage_total_histo_sans_cut[b_sans]={1902,1915,1900,1946,2017,2117,2209,2332,2428,2587,2675,2800,2887,2965,3491,3902,4577,5663,6573,6989,7597,7655,7639};
		Double_t datas_BGO_pourcentage_Fortuitous_histo_sans_cut[b_sans]={7,8,14,70,155,255,387,504,616,763,878,997,1103,1188,1765,2226,3054,4363,5472,6015,6842,7052,7072} ;
		Double_t datas_BGO_pourcentage_True_gamma_histo_sans_cut[b_sans]	={1895,1907,1886,1876,1862,1862,1822,1828,1812,1824,1797,1803,1784,1777,1726,1676,1523,1300,1101,974,755,603,567};
		Double_t datas_BGO_pourcentage_True_Other_gamma_histo_sans_cut[b_sans]= {451,462,453,447,443,455,440,435,439,447,430,443,435,428,419,394,350,287,245,213,153,127,113};
		Double_t datas_BGO_pourcentage_Same_gamma_histo_sans_cut[b_sans]	={1041,1045,1022,1020,1014,1012,1001,999,988,986,986,973,963,963,934,910,824,718,608,539,426,336,319};
		Double_t datas_BGO_pourcentage_Same_Particle_histo_sans_cut[b_sans]	={403,400,411,409,405,395,381,394,385,391,381,387,386,386,373,372,349,295,248,222,176,140,135};
		
		
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
		Double_t datas_BGO_absolu_Fortuitous_sans_cut[b_sans]={1.7e-07,2e-07,3.4e-07,1.78e-06,3.61e-06,6.5e-06,1.015e-05,1.326e-05,1.609e-05,1.998e-05,2.274e-05,2.61e-05,2.862e-05,3.129e-05,4.557e-05,5.765e-05,7.928e-05,0.0001129,0.00014301,0.00015712,0.00017979,0.00018474,0.0001861};
		Double_t datas_BGO_absolu_True_gamma_sans_cut[b_sans]	={3.43e-05,3.447e-05,3.41e-05,3.399e-05,3.383e-05,3.377e-05,3.31e-05,3.311e-05,3.267e-05,3.303e-05,3.246e-05,3.258e-05,3.247e-05,3.226e-05,3.11e-05,3.028e-05,2.764e-05,2.39e-05,2.032e-05,1.789e-05,1.409e-05,1.114e-05,1.041e-05};
		Double_t datas_BGO_absolu_True_Other_gamma_sans_cut[b_sans]= {1.016e-05,1.019e-05,9.96e-06,1.001e-05,9.87e-06,1e-05,9.77e-06,9.62e-06,9.64e-06,9.78e-06,9.47e-06,9.61e-06,9.54e-06,9.51e-06,9.22e-06,8.72e-06,7.78e-06,6.35e-06,5.29e-06,4.78e-06,3.46e-06,2.77e-06,2.58e-06};
		Double_t datas_BGO_absolu_Same_gamma_sans_cut[b_sans]	={1.601e-05,1.594e-05,1.57e-05,1.565e-05,1.56e-05,1.554e-05,1.54e-05,1.534e-05,1.514e-05,1.524e-05,1.512e-05,1.505e-05,1.503e-05,1.496e-05,1.446e-05,1.407e-05,1.3e-05,1.147e-05,9.8e-06,8.69e-06,7.03e-06,5.57e-06,5.21e-06};
		Double_t datas_BGO_absolu_Same_Particle_sans_cut[b_sans]	={8.13e-06,8.34e-06,8.44e-06,8.33e-06,8.36e-06,8.23e-06,7.93e-06,8.15e-06,7.89e-06,8.01e-06,7.87e-06,7.92e-06,7.9e-06,7.79e-06,7.42e-06,7.49e-06,6.86e-06,6.08e-06,5.23e-06,4.42e-06,3.6e-06,2.8e-06,2.62e-06};
		
		Double_t Sum_bad_events_absolu_sans_cut[b_sans];
		
		for(int u=0;u<b_sans;u++)
		{
			Sum_bad_events_absolu_sans_cut[u]= datas_BGO_absolu_Fortuitous_sans_cut[u]+datas_BGO_absolu_True_gamma_sans_cut[u]-datas_BGO_absolu_Same_gamma_sans_cut[u];
		}
		
		
		
		//--------------------------------
		//Reconstructed events (line-cone)
		//--------------------------------
		
		
		Double_t datas_BGO_absolu_total_histo_sans_cut[b_sans]={1902,1915,1900,1946,2017,2117,2209,2332,2428,2587,2675,2800,2887,2965,3491,3902,4577,5663,6573,6989,7597,7655,7639};
		Double_t datas_BGO_absolu_Fortuitous_histo_sans_cut[b_sans]={7,8,14,70,155,255,387,504,616,763,878,997,1103,1188,1765,2226,3054,4363,5472,6015,6842,7052,7072} ;
		Double_t datas_BGO_absolu_True_gamma_histo_sans_cut[b_sans]	={1895,1907,1886,1876,1862,1862,1822,1828,1812,1824,1797,1803,1784,1777,1726,1676,1523,1300,1101,974,755,603,567};
		Double_t datas_BGO_absolu_True_Other_gamma_histo_sans_cut[b_sans]= {451,462,453,447,443,455,440,435,439,447,430,443,435,428,419,394,350,287,245,213,153,127,113};
		Double_t datas_BGO_absolu_Same_gamma_histo_sans_cut[b_sans]	={1041,1045,1022,1020,1014,1012,1001,999,988,986,986,973,963,963,934,910,824,718,608,539,426,336,319};
		Double_t datas_BGO_absolu_Same_Particle_histo_sans_cut[b_sans]	={403,400,411,409,405,395,381,394,385,391,381,387,386,386,373,372,349,295,248,222,176,140,135};
		
		
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
		
		
		
	}
	
	
	
	
	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- IONS CARBONE --------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------	
	
	
	
	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- Avec CUT TOF --------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	
	
	
	
	if (choix_p_c ==2)
	{
	
	
	//-------------------
	//-------------------
	//Coinc Pourcentage
	//-------------------
	//-------------------	
	

	
	
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

	
	
	}
	
	
	
	
	//--------------------------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------- SANS CUT TOF --------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------------------
	
	
	if (choix_p_c ==2)
	{
		
		
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
		
		
		
	}
	
	
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
	
	if (choix_p_c ==1)
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
	
	//----------------------------------------------------------------------------
	//-----------------------Presentation PTCOG Absolu---------------------------
	//----------------------------------------------------------------------------
	
	if(choix==4)
	{

		if(choix_histo==0)
		{	
			char name[150];
			sprintf(name,  "Variation of the beam intensity proton");//Variation of the beam intensity proton

			TH1 *frame = new TH1F("frame","",0,0,70);
			frame->SetMinimum(0.000001);
			frame->SetMaximum(0.01);
			frame->SetDirectory(0);
			frame->SetStats(0);
			frame->SetTitle(name);
			//	frame->	GetXaxis()->SetLabelSize(20);
			//	frame->	GetYaxis()->SetLabelSize(20);
			//frame->GetYaxis()->SetLabelSize(0.001);
			frame->GetXaxis()->SetTitleSize(0.04);
			frame->GetXaxis()->SetTitle("Number of protons per bunch");//carbon ions per bunch
			frame->GetXaxis()->SetTickLength(0.02);
			frame->GetXaxis()->SetLabelSize(0.04);
			frame->GetYaxis()->SetTitleSize(0.04);
			frame->GetYaxis()->SetTitle("Coincidences per incident proton");//incident carbon ion
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
		
		
		if( choix_histo==1)
		{
			char name[150];
			sprintf(name,  "");//Variation of the beam intensity proton");//carbon ion
			
			TH1 *frame = new TH1F("frame","",0,0.1,220);
			frame->SetMinimum(0.000001);//  0.00001
			frame->SetMaximum(0.0005);// 0.01
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
			
			if(comparaison ==1)
			{
				
				cout<<"3 !!!"<<endl;
				TGraph *gr = new TGraphErrors(b,sigma_ns,Sum_bad_events_absolu,0,0); // Background cut TOF
				
				gr->SetLineColor(6); // Options de mise en forme  du graphique
				gr->SetLineStyle(8);
				gr->SetLineWidth(2);
                gr->SetMarkerSize(2);
                gr->SetMarkerStyle(24);
                gr->SetMarkerColor(2);//7
				gr->Draw("p");
				
				TGraph *gr1 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Same_gamma,0,0); // Gamma  cut TOF
				
				//gr1->SetLineColor(2); // Options de mise en forme  du graphique
				//gr1->SetLineStyle(9);
				//gr1->SetLineWidth(2);
                gr1->SetMarkerSize(2);
                gr1->SetMarkerStyle(25);
                gr1->SetMarkerColor(1);//9
				gr1->Draw("p");
				
				TGraph *gr2 = new TGraphErrors(b,sigma_ns,Sum_bad_events_absolu_histo,0,0); // Background cut TOF - reconstruit
				 
				 /*gr2->SetLineColor(8); // Options de mise en forme  du graphique
				 gr2->SetLineStyle(7);
				 gr2->SetLineWidth(2);*/
                 gr2->SetMarkerSize(2);
                 gr2->SetMarkerStyle(26);
                 gr2->SetMarkerColor(4);//38
				 gr2->Draw("P");
				 
				 TGraph *gr3 = new TGraphErrors(b,sigma_ns,datas_BGO_absolu_Same_gamma_histo,0,0); // Gamma cut TOF - reconstruit
				 
				 /*gr3->SetLineColor(2); // Options de mise en forme  du graphique
				 gr3->SetLineStyle(9);
				 gr3->SetLineWidth(2);*/
                 gr3->SetMarkerSize(2);
                 gr3->SetMarkerStyle(28);
                 gr3->SetMarkerColor(6);//15
				 gr3->Draw("P");
				
			}

			

			cout<<"2 !!!"<<endl;
			TGraph *gr4 = new TGraphErrors(b_sans,sigma_ns_sans_cut,Sum_bad_events_absolu_sans_cut,0,0); // Background NO cut TOF
			
			/*gr4->SetLineColor(1); // Options de mise en forme  du graphique
			gr4->SetLineStyle(3);
			gr4->SetLineWidth(2);*/
            gr4->SetMarkerSize(1);
            gr4->SetMarkerStyle(20);
            gr4->SetMarkerColor(2);//2
			gr4->Draw("P");
			
			TGraph *gr5 = new TGraphErrors(b_sans,sigma_ns_sans_cut,datas_BGO_absolu_Same_gamma_sans_cut,0,0); // Gamma NO cut TOF
			
			gr5->SetLineColor(4); // Options de mise en forme  du graphique
			gr5->SetLineStyle(5);
			gr5->SetLineWidth(2);
            gr5->SetMarkerSize(1);
            gr5->SetMarkerStyle(21);//21
            gr5->SetMarkerColor(1);
			gr5->Draw("P");
			
			TGraph *gr6 = new TGraphErrors(b_sans,sigma_ns_sans_cut,Sum_bad_events_absolu_histo_sans_cut,0,0); //Background NO cut TOF - reconstruit
			 
			 gr6->SetLineColor(1); // Options de mise en forme  du graphique
			 gr6->SetLineStyle(10);
			 gr6->SetLineWidth(2);
             gr6->SetMarkerSize(1);
             gr6->SetMarkerStyle(22);
             gr6->SetMarkerColor(4);
			 gr6->Draw("p");
			 
			 TGraph *gr7 = new TGraphErrors(b_sans,sigma_ns_sans_cut,datas_BGO_absolu_Same_gamma_histo_sans_cut,0,0); // Gamma NO cut TOF - reconstruit
			 
			 gr7->SetLineColor(4); // Options de mise en forme  du graphique
			 gr7->SetLineStyle(4);
             gr7->SetLineWidth(2);
             gr7->SetMarkerSize(1);
             gr7->SetMarkerStyle(34);
             gr7->SetMarkerColor(6);
			 gr7->Draw("p");
			
			
			
			if (comparaison == 1)
			{
				leg = new TLegend(0.3,0.75,0.85,0.89);
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
			}
			else
			{
				leg = new TLegend(0.3,0.75,0.85,0.89);
				leg->SetFillColor(10);
				leg->AddEntry(gr4, "Background ","p");
				leg->AddEntry(gr5,"True gamma ","p");
				leg->AddEntry(gr6, "Background (reconstructed) ","p");
				leg->AddEntry(gr7,"True gamma  (reconstructed) ","p");
				cout<<"5 !!!"<<endl;
			}
			//leg->AddEntry(gr4,"True - other particles","l");
			leg->Draw();
			c1->Update();
		}
		
		
	}
	
	
	
	
	
	
	
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
