#include <iostream>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include <TMath.h>
#include "../header/commonUtility.h"
#include "../header/cutsAndBin.h"
#include "../header/HiEvtPlaneList.h"
#include "../header/Style.h"
#include "../header/tdrstyle.C"
#include "../header/CMS_lumi_v2mass.C"
#include "../header/rootFitHeaders.h"
using namespace std;
using namespace RooFit;

using namespace hi;

double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);
void GetHistSqrt(TH1D* h1 =0, TH1D* h2=0);
double GetHistBkg(TH1D* h1 =0, TH1D* h2=0);

void PR_pT45_65_y16_24_c20_120(
    bool isPtW = true, bool isPtWAcc = true, bool isTnP = true, // added
		double ptLow = 4.5, double ptHigh = 6.5,
		double yLow = 1.6, double yHigh = 2.4,
		int cLow = 20, int cHigh = 120,
		int PR=0, // 0: PR, 1: NP, 2: Inc.
		float ctauCut=0.0395,
		float massLow = 3.4, float massHigh = 4.0, 
		bool dimusign=true, bool fAccW = true, bool fEffW = true, bool isMC = false, 
		int dtype = 1, 
		int weight_PR = 0 // PR : 0, NP : 1
		) //PR 0: PR, 1: NP, 2: Inc.
{
	//Basic Setting
	gStyle->SetOptStat(0);
	TString DATE="210913";//210702
	//TString DATE="Corr";
	//TString DATE="210503";
	gStyle->SetEndErrorSize(0);
	gSystem->mkdir(Form("figs"), kTRUE);
	gSystem->mkdir(Form("roots"), kTRUE);
	gSystem->mkdir(Form("figs/v2mass_hist"), kTRUE);
	gSystem->mkdir(Form("figs/decayL"), kTRUE);
	gSystem->mkdir(Form("figs/mass_dist"), kTRUE);
	gSystem->mkdir(Form("figs/q_vector"), kTRUE);
	gSystem->mkdir(Form("roots/%s",DATE.Data()), kTRUE);
	gStyle->SetOptStat(000000000);
	gROOT->ForceStyle();
	setTDRStyle();
	writeExtraText= true;
	int iPeriod = 2;
	int iPos = 33;
	TH1::SetDefaultSumw2();
	//TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh, 1) ;
	TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0, cLow, cHigh) ;
	TString dimusignString;
	if(dimusign) dimusignString = "OS";
	else if(!dimusign) dimusignString = "SS";
	TString sample;
	TString wName;
	TString cutName;
	TString prName;

	if(PR==0) prName="Prompt";
	else if(PR==1) prName="NonPrompt";
	else if(PR==2) prName="Inclusive";


	if(weight_PR==0) wName="prompt";
	else if(weight_PR==1) wName="nprompt";

	cout<<"Data type: "<<dtype<<endl;

	TChain *tree = new TChain("mmepevt");
	if(!isMC){
		//TString f1 = "roots/OniaFlowSkim_JpsiTrig_DBAllPD_isMC0_HFNom_201127.root";
		//###TString f1 = "../../OniaFlowSkim_JpsiTrig_DBAllPD_isMC0_HFNom_201127.root";
		TString f1 = "../../../OniaFlowSkim_JpsiTrig_DBAllPD_isMC0_HFUp_210729.root";
		//TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_JpsiTrig_DBPD_isMC0_HFNom_AddEP_200217.root";
		//TString f2 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_JpsiTrig_DBPeriPD_isMC0_HFNom_AddEP_Peri_200217.root";
		tree->Add(f1.Data());
		//tree->Add(f2.Data());
		sample="RD";}
	else if(dtype==1){
		TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/";
		tree->Add(f1.Data());
		sample="MC_PR";}
	else if(dtype==2){
		TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/";
		tree->Add(f1.Data());
		sample="MC_NP";}

	//READ Input Skimmed File
	//TString f1 = "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/skimmedFiles/OniaFlowSkim_JpsiTrig_DBPD_isMC0_HFNom_190716.root";
	//TString f2 = "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/skimmedFiles/OniaFlowSkim_JpsiTrig_DBPsiPeriPD_isMC0_HFNom_190716.root";
	//TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_JpsiTrig_AllPD_HFNom.root";
	//TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_MC_Jpsi_HFNom.root";
	//TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_BJpsiTrig_DBPD_isMC1_HFNom_200217.root";

	//Get Correction histograms
	TFile *fEff1 = new TFile(Form("../../make_RooDataSet/mc_eff_vs_pt_cent_0_to_20_rap_%s_pbpb_psi2S_PtW%d_tnp%d.root",wName.Data(),isPtW,isTnP),"read");
	TFile *fEff2 = new TFile(Form("../../make_RooDataSet/mc_eff_vs_pt_cent_20_to_120_rap_%s_pbpb_psi2S_PtW%d_tnp%d.root",wName.Data(),isPtW,isTnP),"read");
	TH1D* hEffPt1[5];
	hEffPt1[0] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy0_1p2",isTnP,isPtW));
	hEffPt1[1] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy1p2_1p6",isTnP,isPtW));
	hEffPt1[2] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy1p6_1p8",isTnP,isPtW));
	hEffPt1[3] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy1p8_2p4",isTnP,isPtW));
  hEffPt1[4] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy0_0p6",isTnP,isPtW)); // added

	TH1D* hEffPt2[5];
	hEffPt2[0] = (TH1D*) fEff2 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy0_1p2",isTnP,isPtW));
	hEffPt2[1] = (TH1D*) fEff2 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p2_1p6",isTnP,isPtW));
	hEffPt2[2] = (TH1D*) fEff2 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p6_1p8",isTnP,isPtW));
	hEffPt2[3] = (TH1D*) fEff2 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p8_2p4",isTnP,isPtW));
  hEffPt2[4] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy0_0p6",isTnP,isPtW)); // added

  TFile *fAcc = new TFile(Form("../../make_RooDataSet/acceptance_Prompt_GenOnly_wgt%d_2021Psi2Sv2_20210601.root",isPtWAcc),"read"); // added
	TH1D* hAccPt[3];
	hAccPt[0] = (TH1D*) fAcc -> Get("hAccPt_2021_ally");
	hAccPt[1] = (TH1D*) fAcc -> Get("hAccPt_2021_midy");
	hAccPt[2] = (TH1D*) fAcc -> Get("hAccPt_2021_Fory");


	//SetBranchAddress
	const int nMaxDimu = 1000;
	float mass[nMaxDimu];
	float qxa[nMaxDimu];
	float qya[nMaxDimu];
	float qxb[nMaxDimu]; 
	float qyb[nMaxDimu];
	float qxc[nMaxDimu];
	float qyc[nMaxDimu];
	float qxdimu[nMaxDimu];
	float qydimu[nMaxDimu];
	//  float P[nMaxDimu];
	//  float Px[nMaxDimu];
	//  float Py[nMaxDimu];
	//  float Pz[nMaxDimu];
	float pt[nMaxDimu];
	float y[nMaxDimu]; 
	float pt1[nMaxDimu];
	float pt2[nMaxDimu];
	float eta[nMaxDimu];
	float eta1[nMaxDimu];
	float eta2[nMaxDimu];
	float ctau3D[nMaxDimu];
	float ctau3DErr[nMaxDimu];
	Int_t cBin;
	Int_t event; 
	Int_t nDimu; 
	float vz;
	int recoQQsign[nMaxDimu];
	double weight;

	TBranch *b_event;
	TBranch *b_cBin;
	TBranch *b_nDimu;
	TBranch *b_vz;
	TBranch *b_mass;
	TBranch *b_recoQQsign;
	//  TBranch *b_P;
	//  TBranch *b_Px;
	//  TBranch *b_Py;
	//  TBranch *b_Pz;
	TBranch *b_pt;
	TBranch *b_y;
	TBranch *b_eta;
	TBranch *b_eta1;
	TBranch *b_eta2;
	TBranch *b_ctau3D;
	TBranch *b_ctau3DErr;
	TBranch *b_pt1;
	TBranch *b_pt2;
	TBranch *b_qxa;
	TBranch *b_qxb;
	TBranch *b_qxc;
	TBranch *b_qxdimu;
	TBranch *b_qya;
	TBranch *b_qyb;
	TBranch *b_qyc;
	TBranch *b_qydimu;
	TBranch *b_weight;

	tree -> SetBranchAddress("event", &event, &b_event);
	tree -> SetBranchAddress("cBin", &cBin, &b_cBin);
	tree -> SetBranchAddress("nDimu", &nDimu, &b_nDimu);
	tree -> SetBranchAddress("vz", &vz, &b_vz);
	tree -> SetBranchAddress("recoQQsign", recoQQsign, &b_recoQQsign);
	tree -> SetBranchAddress("mass", mass, &b_mass);
	tree -> SetBranchAddress("y", y, &b_y);
	//  tree -> SetBranchAddress("P", P, &b_P);
	//  tree -> SetBranchAddress("Px", Px, &b_Px);
	//  tree -> SetBranchAddress("Py", Py, &b_Py);
	//  tree -> SetBranchAddress("Pz", Pz, &b_Pz);
	tree -> SetBranchAddress("pt", pt, &b_pt);
	tree -> SetBranchAddress("pt1", pt1, &b_pt1);
	tree -> SetBranchAddress("pt2", pt2, &b_pt2);
	tree -> SetBranchAddress("eta", eta, &b_eta);
	tree -> SetBranchAddress("eta1", eta1, &b_eta1);
	tree -> SetBranchAddress("eta2", eta2, &b_eta2);
	tree -> SetBranchAddress("ctau3D", ctau3D, &b_ctau3D);
	tree -> SetBranchAddress("ctau3DErr", ctau3DErr, &b_ctau3DErr);
	tree -> SetBranchAddress("qxa", qxa, &b_qxa);
	tree -> SetBranchAddress("qxb", qxb, &b_qxb);
	tree -> SetBranchAddress("qxc", qxc, &b_qxc);
	tree -> SetBranchAddress("qxdimu", qxdimu, &b_qxdimu);
	tree -> SetBranchAddress("qya", qya, &b_qya);
	tree -> SetBranchAddress("qyb", qyb, &b_qyb);
	tree -> SetBranchAddress("qyc", qyc, &b_qyc);
	tree -> SetBranchAddress("qydimu", qydimu, &b_qydimu);
	tree -> SetBranchAddress("weight", &weight, &b_weight);

	//crosscheck+tnp+genweight
	//  else if(PR==1){
	//  else if(ptLow==3&&ptHigh==4.5) ctauCut=0.1135;
	//  else if(ptLow==4.5&&ptHigh==6.5) ctauCut=0.0935;
	//  else if(ptLow==6.5&&ptHigh==9) ctauCut=0.0755;
	//  else if(ptLow==9&&ptHigh==12)  ctauCut=0.0635;
	//  else if(ptLow==12&&ptHigh==30) ctauCut=0.0495;}
	//  cout<<"pt["<<ptLow<<" - "<<ptHigh<<" GeV/c], "<<"ctau cut: "<<ctauCut<<endl;
	
	/* ==== pol1, manual bkg params ==== */
	// float massBinDiff[nMassBin+1]={3.4,3.47,3.58,3.62,3.66,3.72,3.80,3.90,4.0}; // pT 4.5 ~ 6.5
	// vector<float> massBinDiff={3.4,3.45,3.50,3.55,3.60,3.65,3.70,3.75,3.80,3.85,3.90,3.95,4.0}; //pT6.5 ~ 10
	// vector<float> massBinDiff={3.4,3.46,3.52,3.58,3.64,3.7,3.76,3.82,3.88,3.94,4}; //pT6.5 ~ 50, c0~ 20
	// vector<float> massBinDiff={3.4,3.465,3.53,3.595,3.66,3.725,3.79,3.855,3.92,4}; // pT4.0 ~ 6.5
	

	// vector<float> massBinDiff;
	// for (auto i = 3.4; i <= 4.0; i += 0.05) {
		// massBinDiff.push_back(i);
		// // cout << i << endl;
	// }
	// if (massBinDiff.back() != 4.0) {
		// massBinDiff.push_back(4.0);
	// }
//
	// for (auto ele : massBinDiff)
		// cout << ele << ",";
	// cout << endl;


	// here
	vector<float> massBinDiff={3.4,3.5,3.58,3.65,3.7,3.75,3.8,3.9,4};
	// vector<float> massBinDiff={3.4,3.48,3.58,3.64,3.70,3.75,3.80,3.90,4.0};

	const int nMassBin = massBinDiff.size() - 1;
	float massBin_[nMassBin+1];

	kineLabel = kineLabel + Form("_m%.1f-%.1f",massLow,massHigh) + "_" + dimusignString;
	cout<<kineLabel<<endl;

	for(int i=0; i<=nMassBin; i++){
		massBin_[i] = massBinDiff[i];
		//massBin[i] = 2.6 + massBinDiff[i]*0.02;
	}

	float* massBin;
	int nMassBin_;
	massBin = massBin_; nMassBin_ = nMassBin;

	int bfevt =-1;
	int afevt =-1;

	double mass_low_SB1 = 2.6; 
	double mass_high_SB1 = 2.9; 
	double mass_low_SB = 2.9;//3.0 
	double mass_high_SB = 3.3; 
	double mass_low_SB2 = 3.3; //3.3
	double mass_high_SB2 = 3.5; 

	int nQBin = 200;
	const int nMass_div = 3;

	TString fSB[nMass_div] = {"SB1 (2.6<m<2.9)","SB2 (3.3<m<3.5)","S (2.9<m<3.3)"};

	//Define drawing histogram
	TH1D* h_v2_1[nMass_div];
	TH1D* h_v2_2[nMass_div];
	TH1D* h_v2_3[nMass_div];
	TH1D* h_v2_4[nMass_div];
	TH1D* h_ljpsi[nMass_div];

	double Q_avg_low = -6500;
	double Q_avg_high = 27000;
	double Q_avg_low_dimu = -150;
	double Q_avg_high_dimu = 150;
	for(int imass=0; imass<nMass_div;imass++){
		h_v2_1[imass] = new TH1D(Form("h_v2_1_%d",imass),";#LTQ_{2}Q_{2A}^{*}#GT;counts",nQBin,Q_avg_low_dimu,Q_avg_high_dimu);
		h_v2_2[imass] = new TH1D(Form("h_v2_2_%d",imass),";#LTQ_{2A}Q_{2B}^{*}#GT;counts",nQBin,Q_avg_low,Q_avg_high);
		h_v2_3[imass] = new TH1D(Form("h_v2_3_%d",imass),";#LTQ_{2A}Q_{2C}^{*}#GT;counts",nQBin,Q_avg_low,Q_avg_high);
		h_v2_4[imass] = new TH1D(Form("h_v2_4_%d",imass),";#LTQ_{2B}Q_{2C}^{*}#GT;counts",nQBin,Q_avg_low,Q_avg_high);
		h_ljpsi[imass] = new TH1D(Form("h_ljBkg_%d",imass),";l_{J/#psi(mm)};Counts/(0.03 mm)",1050,-1.5,2);
	}

	TH1D* h_v2_num_q1 = new TH1D("h_v2_num_q1",";m_{#mu^{+}#mu^{-}};#LTQ_{2}Q_{2A}^{*}#GT" ,nMassBin,massBin);
	TH1D* h_v2_den_q2 = new TH1D("h_v2_num_q2",";m_{#mu^{+}#mu^{-}};#LTQ_{2A}Q_{2B}^{*}#GT",nMassBin,massBin);
	TH1D* h_v2_den_q3 = new TH1D("h_v2_num_q3",";m_{#mu^{+}#mu^{-}};#LTQ_{2A}Q_{2C}^{*}#GT",nMassBin,massBin);
	TH1D* h_v2_den_q4 = new TH1D("h_v2_num_q4",";m_{#mu^{+}#mu^{-}};#LTQ_{2B}Q_{2C}^{*}#GT",nMassBin,massBin);

	TH1D* h_mass = new TH1D("h_mass",";m_{#mu^{+}#mu^{-}};Counts/(0.025 GeV)",36,massLow,massHigh);
	TH1D* h_decayL = new TH1D("h_decayL",";l_{J/#psi(mm)};Counts/(0.1 mm)",100,-4,6);
	//TH1D* h_ljpsi_bkg = new TH1D("h_ljpsi_bkg",";l_{J/#psi(mm)};Counts/(0.03 mm)",1050,-1.5,2);

	RooRealVar* massVar = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
	RooRealVar* weightVar = new RooRealVar("weight", "weight var", 0, 10000,"");

	RooArgSet* argSet = new RooArgSet(*massVar, *weightVar);
	RooDataSet* dataSet = new RooDataSet("dataset","a dataset",*argSet); 

	const static int countMax = 1000000;
	vector<vector<double>> v2_1(nMassBin,vector<double> (countMax,0));
	vector<vector<double>> v2_2(nMassBin,vector<double> (countMax,0));
	vector<vector<double>> v2_3(nMassBin,vector<double> (countMax,0));
	vector<vector<double>> v2_4(nMassBin,vector<double> (countMax,0));
	vector<vector<double>> v2_1_raw(nMassBin,vector<double> (countMax,0));
	vector<vector<double>> v2_2_raw(nMassBin,vector<double> (countMax,0));
	vector<vector<double>> v2_3_raw(nMassBin,vector<double> (countMax,0));
	vector<vector<double>> v2_4_raw(nMassBin,vector<double> (countMax,0));
	vector<vector<double>> weight_dimu(nMassBin,vector<double> (countMax,0));

	vector<double> v2_1_avg(nMassBin,0);
	vector<double> v2_2_avg(nMassBin,0);
	vector<double> v2_3_avg(nMassBin,0);
	vector<double> v2_4_avg(nMassBin,0);

	int dbcount=0;
	vector<unsigned int> count(nMassBin,0);
	vector<unsigned int> count_ss(nMassBin,0);

	vector<double> weight_s(nMassBin,0);

	double weight_acc = 1;
	double weight_eff = 1;

	Int_t nEvt = tree->GetEntries();
	cout << "nEvt : " << nEvt << endl;

	int nDimu_all=0;
	int nDimuPass=0;
	int nDimu_one=0;
	int nDimu_more=0;
	//Begin Loop
	for(int i=0; i<nEvt; i++){
		tree->GetEntry(i);
		nDimuPass=0;
		if(cBin>=cLow&&cBin<cHigh){ 
			if(fabs(vz)<15){
				for(int j=0; j<nDimu; j++){
					//if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && (double)pt1[j]>SiMuPtCut&&(double)pt2[j]>SiMuPtCut&&abs((double)eta1[j])<2.4&&abs((double)eta2[j])<2.4)) continue;
					if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && IsAcceptanceQQ(pt1[j],eta1[j]) && IsAcceptanceQQ(pt2[j],eta2[j])) ) continue;
					nDimuPass++;
				}
				nDimu_all++;
				if(nDimuPass>1) {nDimu_more++; continue;}
				if(nDimuPass==1) nDimu_one++;
				//		if(PR==0 && ctau3D[i]>ctauCut) continue;
				//		else if(PR==1 && ctau3D[i]<ctauCut) continue;
				// Fill Dimuon Loop

				for(int j=0; j<nDimu; j++){
					if(PR==0 && ctau3D[j]>ctauCut) continue;
					else if(PR==1 && ctau3D[j]<ctauCut) continue;

					if(pt[j]>ptLow&&pt[j]<ptHigh&&recoQQsign[j]==0&&mass[j]>massLow&&mass[j]<massHigh&&abs(y[j])>yLow&&abs(y[j])<yHigh && IsAcceptanceQQ(pt1[j],eta1[j]) && IsAcceptanceQQ(pt2[j],eta2[j])
							//( ((abs(eta1[j]) <= 1.2) && (pt1[j] >=3.5)) || ((abs(eta2[j]) <= 1.2) && (pt2[j] >=3.5)) ||
							//  ((abs(eta1[j]) > 1.2) && (abs(eta1[j]) <= 2.1) && (pt1[j] >= 5.47-1.89*(abs(eta1[j])))) || 
							//  ((abs(eta2[j]) > 1.2)  && (abs(eta2[j]) <= 2.1) && (pt2[j] >= 5.47-1.89*(abs(eta2[j])))) ||
							//  ((abs(eta1[j]) > 2.1) && (abs(eta1[j]) <= 2.4) && (pt1[j] >= 1.5)) || 
							//  ((abs(eta2[j]) > 2.1)  && (abs(eta2[j]) <= 2.4) && (pt2[j] >= 1.5)) )
					  ){
						weight_acc=1;
						weight_eff=1;
						if(fAccW){
							//if ( abs((double)y[j])<2.4 ) { weight_acc = getAccWeight(hAccPt[1], pt[j]); }
							if ( abs((double)y[j])<1.6 ) { weight_acc = getAccWeight(hAccPt[1], pt[j]); }
							else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<2.4 ) { weight_acc = getAccWeight(hAccPt[2], pt[j]); }
							//cout<<"weight acc : "<<weight_acc<<endl;
							//if (mass[j]>3.0&&mass[j]<3.06) {cout << "Acc Weight : " << (double)1.0/weight_acc << " y : " << y[j] << " pT : " << pt[j] << endl;}
						}
						if(fEffW){
							//if ( abs((double)y[j])<2.4) { weight_eff = getEffWeight(hEffPt[0], pt[j]); }
							//if ( abs((double)y[j])<0.6){ 
							//  if(getEffWeight(hEffPt[0], pt[j])<=0){
							//    //weight_eff = getEffWeight(hEffPt[0], pt[j]); 
							//    weight_eff=hEffPt[0]->GetBinContent(0);
							//  }
							//  else  weight_eff = getEffWeight(hEffPt[0], pt[j]);}
							if(cBin >= 0 && cBin < 20) {
								if ( abs((double)y[j])>=0.0 && abs((double)y[j])<1.2 ) { 
                  weight_eff = getEffWeight(hEffPt1[0], pt[j]);
                }
								else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) { weight_eff = getEffWeight(hEffPt1[1], pt[j]); }
								else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<1.8 ) { weight_eff = getEffWeight(hEffPt1[2], pt[j]); }
								else if ( abs((double)y[j])>=1.8 && abs((double)y[j])<2.4 ) { weight_eff = getEffWeight(hEffPt1[3], pt[j]); }
							}else if(cBin >= 20 && cBin < 120) {
								if ( abs((double)y[j])>=0.0 && abs((double)y[j])<1.2 ) { 
                  weight_eff = getEffWeight(hEffPt2[0], pt[j]);
                }
								else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) { weight_eff = getEffWeight(hEffPt2[1], pt[j]); }
								else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<1.8 ) { weight_eff = getEffWeight(hEffPt2[2], pt[j]); }
								else if ( abs((double)y[j])>=1.8 && abs((double)y[j])<2.4 ) { weight_eff = getEffWeight(hEffPt2[3], pt[j]); }
							} 

						}
						double weight_ = weight * weight_eff * weight_acc;
						//weight_=1.;
						for(int imbin=0; imbin<nMassBin; imbin++){
							if(mass[j]>=massBin[imbin] && mass[j]<massBin[imbin+1]){
								v2_1[imbin][count[imbin]] = (qxa[j]*qxdimu[j] + qya[j]*qydimu[j])*weight_;
								v2_2[imbin][count[imbin]] = (qxa[j]*qxb[j] + qya[j]*qyb[j])*weight_;
								v2_3[imbin][count[imbin]] = (qxa[j]*qxc[j] + qya[j]*qyc[j])*weight_;
								v2_4[imbin][count[imbin]] = (qxb[j]*qxc[j] + qyb[j]*qyc[j])*weight_;
								v2_1_raw[imbin][count[imbin]] = (qxa[j]*qxdimu[j] + qya[j]*qydimu[j]);
								v2_2_raw[imbin][count[imbin]] = (qxa[j]*qxb[j] + qya[j]*qyb[j]);
								v2_3_raw[imbin][count[imbin]] = (qxa[j]*qxc[j] + qya[j]*qyc[j]);
								v2_4_raw[imbin][count[imbin]] = (qxb[j]*qxc[j] + qyb[j]*qyc[j]);

								v2_1_avg[imbin] += v2_1[imbin][count[imbin]];
								v2_2_avg[imbin] += v2_2[imbin][count[imbin]];
								v2_3_avg[imbin] += v2_3[imbin][count[imbin]];
								v2_4_avg[imbin] += v2_4[imbin][count[imbin]];

								weight_dimu[imbin][count[imbin]] = weight_;

								count[imbin]++;
								weight_s[imbin] += weight_;
							}
						}
						if(mass[j]>=mass_low_SB1 && mass[j]<mass_high_SB1){
							h_v2_1[0]->Fill(qxa[j]*qxdimu[j] + qya[j]*qydimu[j], weight_);
							h_v2_2[0]->Fill(qxa[j]*qxb[j] + qya[j]*qyb[j], weight_);
							h_v2_3[0]->Fill(qxa[j]*qxc[j] + qya[j]*qyc[j], weight_);
							h_v2_4[0]->Fill(qxb[j]*qxc[j] + qyb[j]*qyc[j], weight_);
							h_ljpsi[0]->Fill(ctau3D[j],weight_);
						}
						else if(mass[j]>=mass_low_SB2 && mass[j]<mass_high_SB2){
							h_v2_1[1]->Fill(qxa[j]*qxdimu[j] + qya[j]*qydimu[j], weight_);
							h_v2_2[1]->Fill(qxa[j]*qxb[j] + qya[j]*qyb[j], weight_);
							h_v2_3[1]->Fill(qxa[j]*qxc[j] + qya[j]*qyc[j], weight_);
							h_v2_4[1]->Fill(qxb[j]*qxc[j] + qyb[j]*qyc[j], weight_);
							h_ljpsi[1]->Fill(ctau3D[j],weight_);
						}
						else if(mass[j]>=mass_low_SB && mass[j]<mass_high_SB){
							h_v2_1[2]->Fill(qxa[j]*qxdimu[j] + qya[j]*qydimu[j], weight_);
							h_v2_2[2]->Fill(qxa[j]*qxb[j] + qya[j]*qyb[j], weight_);
							h_v2_3[2]->Fill(qxa[j]*qxc[j] + qya[j]*qyc[j], weight_);
							h_v2_4[2]->Fill(qxb[j]*qxc[j] + qyb[j]*qyc[j], weight_);
							h_ljpsi[2]->Fill(ctau3D[j],weight_);
						}
						h_decayL->Fill(ctau3D[j]);
						massVar->setVal((double)mass[j]);
						weightVar->setVal((double)weight_);
						dataSet->add(*argSet);
						h_mass->Fill(mass[j],weight_);
					}
				}
			}
		}
	}
	cout<<"How many Jpsi??: "<<float(h_mass->GetEntries())<<endl;

	TLine *lcut;
	//for(int bins=0; bins<h_decayL->GetNbinsX(); bins++){
	//  h_decayL->Integral(0,bins);
	//  if(h_decayL->Integral(0,bins)>=0.79*h_decayL->GetEntries()
	//      &&h_decayL->Integral(0,bins)<=0.81*h_decayL->GetEntries()){ 
	//    cout<<"THE BIN IS: "<<bins<<endl;
	//    lcut = new TLine(h_decayL->GetBinCenter(bins), 0, h_decayL->GetBinCenter(bins), h_decayL->GetMaximum());
	//  }
	//}
	//for(int bins=0; bins<h_decayL->GetNbinsX(); bins++){
	//  h_decayL->Integral(0,bins);
	//  if(h_decayL->Integral(0,bins)<=0.89*h_decayL->GetEntries() || h_decayL->Integral(0,bins)>0.90*h_decayL->GetEntries()) continue;
	//    cout<<"ANSWER: "<<bins<<"th bin"<<", ljpsi: "<<h_decayL->GetBinCenter(bins)<<endl;
	//    lcut = new TLine(h_decayL->GetBinCenter(bins), 0, h_decayL->GetBinCenter(bins), 1);
	//}

	int nMassFrameBin = 36;
	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataSet);
	ws->data("dataset")->Print();
	RooDataSet *reducedDS = new RooDataSet("reducedDS","A sample",*dataSet->get(),Import(*dataSet),WeightVar(*ws->var("weight")));
	reducedDS->SetName("reducedDS");
	ws->import(*reducedDS);
	ws->var("mass")->setRange(massLow, massHigh);
	ws->var("mass")->Print();

	RooPlot* myPlot = ws->var("mass")->frame(nMassFrameBin);
	ws->data("reducedDS")->plotOn(myPlot,Name("massDataHist"));

	RooHist* hist = (RooHist*) myPlot->findObject("massDataHist");
	int nHistP = hist->GetN();
	cout << "nHistP : " << nHistP << endl;
	if(nHistP != nMassFrameBin) cout << "ERROR::: INCONSISTENT NUMBER OF BINS" << endl;

	TGraphAsymmErrors *g_mass = new TGraphAsymmErrors();
	g_mass->SetName("g_mass");
	Double_t xp, yp, ypl, yph;
	for(int ip = 0; ip < nHistP; ip++){
		xp=0; yp=0; ypl=0; yph=0;
		hist->GetPoint(ip,xp,yp);
		ypl = hist->GetErrorYlow(ip);
		yph = hist->GetErrorYhigh(ip);
		g_mass->SetPoint(ip,xp,yp);
		g_mass->SetPointError(ip,0,0,ypl,yph);
	}

	TF1* fitmc1;

	vector<double> v2_1_err(nMassBin,0);
	vector<double> v2_2_err(nMassBin,0);
	vector<double> v2_3_err(nMassBin,0);
	vector<double> v2_4_err(nMassBin,0);

	for(int ibin=0; ibin<nMassBin; ibin++){
		//v2_1_avg[ibin] = v2_1_avg[ibin]/count[ibin];
		//v2_2_avg[ibin] = v2_2_avg[ibin]/count[ibin];
		//v2_3_avg[ibin] = v2_3_avg[ibin]/count[ibin];
		//v2_4_avg[ibin] = v2_4_avg[ibin]/count[ibin];

		v2_1_avg[ibin] = v2_1_avg[ibin]/weight_s[ibin];
		v2_2_avg[ibin] = v2_2_avg[ibin]/weight_s[ibin];
		v2_3_avg[ibin] = v2_3_avg[ibin]/weight_s[ibin];
		v2_4_avg[ibin] = v2_4_avg[ibin]/weight_s[ibin];

		for(int icount=0; icount<count[ibin]; icount++){
			v2_1_err[ibin] += (v2_1_raw[ibin][icount]-v2_1_avg[ibin])*(v2_1_raw[ibin][icount]-v2_1_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
			v2_2_err[ibin] += (v2_2_raw[ibin][icount]-v2_2_avg[ibin])*(v2_2_raw[ibin][icount]-v2_2_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
			v2_3_err[ibin] += (v2_3_raw[ibin][icount]-v2_3_avg[ibin])*(v2_3_raw[ibin][icount]-v2_3_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
			v2_4_err[ibin] += (v2_4_raw[ibin][icount]-v2_4_avg[ibin])*(v2_4_raw[ibin][icount]-v2_4_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
		}

		//v2_1_err[ibin] = TMath::Sqrt(v2_1_err[ibin])/count[ibin];
		//v2_2_err[ibin] = TMath::Sqrt(v2_2_err[ibin])/count[ibin];
		//v2_3_err[ibin] = TMath::Sqrt(v2_3_err[ibin])/count[ibin];
		//v2_4_err[ibin] = TMath::Sqrt(v2_4_err[ibin])/count[ibin];

		v2_1_err[ibin] = TMath::Sqrt(v2_1_err[ibin])/weight_s[ibin];
		v2_2_err[ibin] = TMath::Sqrt(v2_2_err[ibin])/weight_s[ibin];
		v2_3_err[ibin] = TMath::Sqrt(v2_3_err[ibin])/weight_s[ibin];
		v2_4_err[ibin] = TMath::Sqrt(v2_4_err[ibin])/weight_s[ibin];

		h_v2_num_q1->SetBinContent(ibin+1,v2_1_avg[ibin]);
		h_v2_num_q1->SetBinError(ibin+1,v2_1_err[ibin]);
		h_v2_den_q2->SetBinContent(ibin+1,v2_2_avg[ibin]);
		h_v2_den_q2->SetBinError(ibin+1,v2_2_err[ibin]);
		h_v2_den_q3->SetBinContent(ibin+1,v2_3_avg[ibin]);
		h_v2_den_q3->SetBinError(ibin+1,v2_3_err[ibin]);
		h_v2_den_q4->SetBinContent(ibin+1,v2_4_avg[ibin]);
		h_v2_den_q4->SetBinError(ibin+1,v2_4_err[ibin]);

		cout << ibin << "th Bin : " << count[ibin] << ",  weight sum  : " << weight_s[ibin] << endl;
		cout << "v2_1_avg " << ibin << " : " << v2_1_avg[ibin] << endl;
		cout << "v2_2_avg " << ibin << " : " << v2_2_avg[ibin] << endl;
		cout << "v2_3_avg " << ibin << " : " << v2_3_avg[ibin] << endl;
		cout << "v2_4_avg " << ibin << " : " << v2_4_avg[ibin] << endl;

		cout << "h_v2_num_q1 " << ibin << ", val : " << h_v2_num_q1->GetBinContent(ibin+1) << " err : " << h_v2_num_q1->GetBinError(ibin+1) << endl;
		cout << "h_v2_den_q2 " << ibin << ", val : " << h_v2_den_q2->GetBinContent(ibin+1) << " err : " << h_v2_den_q2->GetBinError(ibin+1) << endl;
		cout << "h_v2_den_q3 " << ibin << ", val : " << h_v2_den_q3->GetBinContent(ibin+1) << " err : " << h_v2_den_q3->GetBinError(ibin+1) << endl;
		cout << "h_v2_den_q4 " << ibin << ", val : " << h_v2_den_q4->GetBinContent(ibin+1) << " err : " << h_v2_den_q4->GetBinError(ibin+1) << endl;
	}

	TH1D* h_v2_den_ = (TH1D*)h_v2_den_q2->Clone("h_v2_den_");
	h_v2_den_->Multiply(h_v2_den_q3);
	h_v2_den_->Divide(h_v2_den_q4);

	TH1D* h_v2_den = (TH1D*) h_v2_den_->Clone("h_v2_den"); //h_v2_den->Reset();
	GetHistSqrt(h_v2_den_,h_v2_den);

	TH1D* h_v2_final = (TH1D*) h_v2_num_q1 -> Clone("h_v2_SplusB");
	h_v2_final->Divide(h_v2_den);

	SetHistStyle(h_v2_final,0,0);
	SetGraphStyle2(g_mass,0,0);
	//g_mass->GetYaxis()->SetLimits(0,10000);
	//g_mass->GetYaxis()->SetRangeUser(0,);
	g_mass->GetYaxis()->SetTitle("Events/(0.025 GeV/c^{2})");
	g_mass->SetMarkerSize(1);
	g_mass->GetYaxis()->SetLabelSize(0.055);
	g_mass->GetXaxis()->SetLabelSize(0.04);
	g_mass->GetYaxis()->SetTitleSize(0.055);
	g_mass->GetYaxis()->SetTitleOffset(1.7);
	g_mass->GetXaxis()->SetTitleSize(0.065);

	h_v2_final->GetYaxis()->SetRangeUser(-0.05,0.26);
	h_v2_final->GetYaxis()->SetTitle("v_{2}^{S+B}");
	h_v2_final->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
	h_v2_final->GetYaxis()->SetLabelSize(0.055);
	h_v2_final->GetXaxis()->SetLabelSize(0.055);
	h_v2_final->GetXaxis()->SetTitleSize(0.07);
	h_v2_final->GetYaxis()->SetTitleSize(0.07);
	h_v2_final->GetYaxis()->SetTitleOffset(1.2);
	h_v2_final->GetXaxis()->SetTitleOffset(1.0);
	h_v2_final->GetXaxis()->SetLabelOffset(0.011);
	//TLine *l=new TLine(2.6, 0, 3.5, 0);
	//l->SetLineColor(kBlack);

	float pos_x = 0.43;
	float pos_x_mass = 0.73;
	float pos_y = 0.7;
	float pos_y_diff = 0.071;
	int text_color = 1;
	float text_size = 16;
	TString perc = "%";

	//TGaxis::SetMaxDigits(2);
	TCanvas* c_mass_v2 = new TCanvas("c_mass_v2","",590,750);
	TPad* pad1 = new TPad("pad1","pad1",0,0.5,1.0,1.0);
	TPad* pad2 = new TPad("pad2","pad2",0,0.0,1.0,0.5);
	c_mass_v2->cd();
	pad1->SetTicks(1,1);
	pad1->SetBottomMargin(0);
	pad1->SetLeftMargin(0.19);
	pad1->SetTopMargin(0.08);
	pad1->cd();
	// pad1->SetLogy();
	g_mass->GetXaxis()->SetTitleSize(0);
	g_mass->GetXaxis()->SetLabelSize(0);
	g_mass->Draw("AP");

	Double_t YMax = h_mass->GetBinContent(h_mass->GetMaximumBin());
	// Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
	Double_t YMin = 1e99;
	for (int i=1; i<=h_mass->GetNbinsX(); i++) if (h_mass->GetBinContent(i)>0) YMin = min(YMin, h_mass->GetBinContent(i));
	Double_t Yup,Ydown;
	Ydown = YMin/(TMath::Power((YMax/YMin), (0.1/(1.0-0.1-0.4))));
	Yup = YMax*TMath::Power((YMax/YMin), (0.4/(1.0-0.1-0.4)));
	g_mass->GetYaxis()->SetRangeUser(Ydown,Yup);
	cout<<"Range: "<<YMin<<"-"<<YMax<<endl;
	cout<<"Range: "<<Ydown<<"-"<<Yup<<endl;
	g_mass->GetXaxis()->SetLimits(massLow,massHigh);
	g_mass->GetXaxis()->SetRangeUser(massLow,massHigh);

	if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.1f GeV/c", ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	else if(ptLow ==6.5) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	else if(ptHigh==6.5) drawText(Form("%1.f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	else if(ptLow!=0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
	else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
	//drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
	drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
	drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
	//    drawText(Form("#font[12]{l}_{J/#psi} < %.4f", ctauCut),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.17);
	pad2->SetLeftMargin(0.19);
	pad2->SetTicks();
	pad2->cd();
	////////////pad2//////////
	// jumSun(massLow,0,massHigh,0,1,1);
	h_v2_final->Draw("P");
	// l->Draw("same");

	//CMS_lumi_v2mass(pad1,iPeriod,iPos,1);  
	CMS_lumi_v2mass(pad1,iPeriod,iPos);  
	pad1->Update();
	pad2->Update();
	c_mass_v2->cd();
	pad1->Draw();
	pad2->Draw();
	//gSystem->mkdir(Form("figs/q_vector",sample.Data()),1);
	c_mass_v2->SaveAs(Form("figs/v2mass_hist/Psi2S_%s_v2Mass_%s.pdf", prName.Data(),kineLabel.Data()));



	TCanvas* c_decayL = new TCanvas("c_decayL","",600,600);
	c_decayL->SetLogy();
	if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	else if(ptLow!=0) drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
	else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
	//drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
	drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
	drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
	h_decayL->Draw("P");
	//lcut->Draw("same");
	c_decayL->Update();
	c_decayL->SaveAs(Form("figs/decayL/Psi2S_%s_decayL_%s.pdf", prName.Data(),kineLabel.Data()));

	TCanvas* c_mass = new TCanvas("c_mass","",600,600);
	c_mass->cd();
	h_mass->Draw("P");
	if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.1f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	else if(ptLow!=0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
	else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
	//drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
	drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
	drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
	//   drawText(Form("#font[12]{l}_{J/#psi} < %.4f", ctauCut),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);
	//  drawText(Form("#font[12]{l}_{J/#psi} > %.4f", ctauCut),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);
	//CMS_lumi_v2mass(c_mass,iPeriod,iPos,0);
	CMS_lumi_v2mass(c_mass,iPeriod,iPos);
	c_mass->Update();
	c_mass->SaveAs(Form("figs/mass_dist/Psi2S_%s_MassDist_%s.pdf", prName.Data(),kineLabel.Data()));

	TLegend *leg_v2_1 = new TLegend(0.38,0.64,0.77,0.9);
	TLegend *leg_v2_2 = new TLegend(0.38,0.64,0.77,0.9);
	TLegend *leg_v2_3 = new TLegend(0.38,0.64,0.77,0.9);
	TLegend *leg_v2_4 = new TLegend(0.38,0.64,0.77,0.9);
	SetLegendStyle(leg_v2_1);
	SetLegendStyle(leg_v2_2);
	SetLegendStyle(leg_v2_3);
	SetLegendStyle(leg_v2_4);

	double mean_v2_1[nMass_div];
	double mean_v2_2[nMass_div];
	double mean_v2_3[nMass_div];
	double mean_v2_4[nMass_div];

	double mean_err_v2_1[nMass_div];
	double mean_err_v2_2[nMass_div];
	double mean_err_v2_3[nMass_div];
	double mean_err_v2_4[nMass_div];

	for(int i=0;i<nMass_div;i++){
		SetHistStyle(h_v2_1[i],i,i);
		SetHistStyle(h_v2_2[i],i,i);
		SetHistStyle(h_v2_3[i],i,i);
		SetHistStyle(h_v2_4[i],i,i);

		scaleInt(h_v2_1[i]);
		scaleInt(h_v2_2[i]);
		scaleInt(h_v2_3[i]);
		scaleInt(h_v2_4[i]);

		mean_v2_1[i] = h_v2_1[i]->GetMean();
		mean_v2_2[i] = h_v2_2[i]->GetMean();
		mean_v2_3[i] = h_v2_3[i]->GetMean();
		mean_v2_4[i] = h_v2_4[i]->GetMean();

		mean_err_v2_1[i] = h_v2_1[i]->GetMeanError();
		mean_err_v2_2[i] = h_v2_2[i]->GetMeanError();
		mean_err_v2_3[i] = h_v2_3[i]->GetMeanError();
		mean_err_v2_4[i] = h_v2_4[i]->GetMeanError();

		leg_v2_1->AddEntry(h_v2_1[i],(fSB[i]+Form(" mean:%.2f",mean_v2_1[i])+Form(" err:%.2f",mean_err_v2_1[i])).Data(),"l");
		leg_v2_2->AddEntry(h_v2_2[i],(fSB[i]+Form(" mean:%.2f",mean_v2_2[i])+Form(" err:%.2f",mean_err_v2_2[i])).Data(),"l");
		leg_v2_3->AddEntry(h_v2_3[i],(fSB[i]+Form(" mean:%.2f",mean_v2_3[i])+Form(" err:%.2f",mean_err_v2_3[i])).Data(),"l");
		leg_v2_4->AddEntry(h_v2_4[i],(fSB[i]+Form(" mean:%.2f",mean_v2_4[i])+Form(" err:%.2f",mean_err_v2_4[i])).Data(),"l");
	}
	//cout<<"what?? : "<<h_v2_1[0]->GetEntries()<<endl;

	//TCanvas *c_ljBkg = new TCanvas("c_ljBkg","",600,600);
	//c_ljBkg->cd();
	//h_ljpsi[0]->Draw("hist");
	//h_ljpsi[0]->SetLineColor(kBlue+2);
	//h_ljpsi[1]->Draw("hist same");
	//h_ljpsi[1]->SetLineColor(kRed+2);
	//h_ljpsi[2]->Draw("hist same");
	//h_ljpsi[2]->SetLineColor(kBlack);

	//TH1D* h_ljpsi_bkg = (TH1D*) h_ljpsi[0]->Clone("h_ljpsi_bkg");
	//h_ljpsi_bkg->Scale(h_ljpsi[0]->GetEntries()/(h_ljpsi[0]->GetEntries()+h_ljpsi[1]->GetEntries()));
	//h_ljpsi_bkg->Draw("EP same");
	//h_ljpsi_bkg->SetMarkerColor(kRed+2);
	//h_ljpsi_bkg->SetMarkerStyle(4);

	TCanvas *c_qq_1 = new TCanvas("c_qqa","",600,600);
	c_qq_1->cd();
	h_v2_1[2]->Draw("hist");
	h_v2_1[1]->Draw("hist same");
	h_v2_1[0]->Draw("hist same");
	leg_v2_1->Draw("same");
	leg_v2_1->SetFillStyle(0);
	c_qq_1->SaveAs(Form("figs/q_vector/Psi2S_%s_c_qqa_%s_Eff%d_Acc%d_PtW%d_TnP%d.pdf",prName.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

	TCanvas *c_qq_2 = new TCanvas("c_qaqb","",600,600);
	c_qq_2->cd();
	h_v2_2[2]->Draw("hist");
	h_v2_2[1]->Draw("hist same");
	h_v2_2[0]->Draw("hist same");
	leg_v2_2->Draw("same");
	c_qq_2->SaveAs(Form("figs/q_vector/Psi2S_%s_c_qaqb_%s_Eff%d_Acc%d_PtW%d_TnP%d.pdf",prName.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

	TCanvas *c_qq_3 = new TCanvas("c_qaqc","",600,600);
	c_qq_3->cd();
	h_v2_3[2]->Draw("hist");
	h_v2_3[1]->Draw("hist same");
	h_v2_3[0]->Draw("hist same");
	leg_v2_3->Draw("same");
	c_qq_3->SaveAs(Form("figs/q_vector/Psi2S_%s_c_qaqc_%s_Eff%d_Acc%d_PtW%d_TnP%d.pdf",prName.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

	TCanvas *c_qq_4 = new TCanvas("c_qbqc","",600,600);
	c_qq_4->cd();
	h_v2_4[2]->Draw("hist");
	h_v2_4[1]->Draw("hist same");
	h_v2_4[0]->Draw("hist same");
	leg_v2_4->Draw("same");
	c_qq_4->SaveAs(Form("figs/q_vector/Psi2S_%s_c_qbqc_%s_Eff%d_Acc%d_PtW%d_TnP%d.pdf",prName.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

	TFile *wf; 
	wf = new TFile(Form("roots/%s/Psi2S_%s_%s_Eff%d_Acc%d_PtW%d_PtWAcc%d_TnP%d.root",DATE.Data(),prName.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isPtWAcc,isTnP),"recreate");
	wf->cd();
	h_v2_final->Write();
	h_decayL->Write();
	g_mass->Write();
	h_mass->Write();

	// Get # of entries for every bins.
	// To compare events number point by point
	// Due to there was events number issue
	// ofstream fout(Form("entry_check_v2mass_hist_Psi2S_%s_%s.txt", prName.Data(),kineLabel.Data()));
	// for (int i = 1; i <= h_mass->GetNbinsX(); i++)
	// fout << h_mass->GetBinContent(i) << endl;

}

void GetHistSqrt(TH1D* h1, TH1D* h2){
	if(h1->GetNbinsX() != h2->GetNbinsX()){ cout << "Inconsistent # of bins b/w histograms !! " << endl;}
	double content;
	double err;
	for(int i=1; i<=h1->GetNbinsX(); i++){
		content=0;err=0;
		content = h1->GetBinContent(i);
		err = h1->GetBinError(i);
		err = 0.5*err*TMath::Power(content,-0.5);
		h2->SetBinContent(i,TMath::Sqrt(content));
		h2->SetBinError(i,err);
	}
} 

double GetHistBkg(TH1D* h1, TH1D* h2){
	double Nh1;
	double Nh2;
	double Nh12;
	Nh1 = h1->GetEntries();
	Nh2 = h2->GetEntries();
	Nh12 = Nh1+Nh2;
	return Nh1/(Nh1+Nh2);

}

double getAccWeight(TH1D* h, double pt){
	double binN = h->FindBin(pt);
	double weight_ = 1./(h->GetBinContent(binN));
	return weight_;
}


double getEffWeight(TH1D *h, double pt){
	double binN = h->FindBin(pt);
	TF1 *eff1 = (TF1*)h->FindObject("f1");
	double eff = eff1->Eval(pt);
	double weight_ = 1./h->GetBinContent(binN);;
	//double weight_ = 1./eff;
	return weight_;
}

