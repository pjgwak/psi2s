#include <iostream>
#include "../header/commonUtility.h"
#include "../header/cutsAndBin.h"
#include "../header/HiEvtPlaneList.h"
#include "../header/Style.h"
#include "../header/tdrstyle.C"
#include "../header/CMS_lumi_v2mass.C"
#include "../header/rootFitHeaders.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"

using namespace std;
using namespace RooFit;

using namespace hi;

double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);
void GetHistSqrt(TH1D* h1 =0, TH1D* h2=0);

void makeRooDataSet_Psi_2S_0604_v3(
		bool isMC = false, 
		bool fAccW = true, bool fEffW = true, 
		int state=1) //state 0: inclusive, state 1: Prompt, state 2: NonPrompt
{
  //Basic Setting
  gStyle->SetOptStat(0);
  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  TString dimusignString;
  TString bCont;
  TString outName;
  if (state==1) { bCont = "prompt"; outName="PR"; }
  else if (state==2) {bCont = "nprompt"; outName="NP"; }
  else if (state==0) {outName="Inclusive"; }
 
  //READ Input Skimmed File
  TFile *rf;
  if(isMC){
    if(state==1) rf = new TFile("../primary_input/OniaFlowSkim_JpsiTrig_isMC1_Prompt_HFNom_200407.root","read");
    else if(state==2) rf = new TFile("../primary_input/OniaFlowSkim_JpsiTrig_isMC1_NonPrompt_HFNom_200407.root","read");
  }
  //else if(!isMC) rf = new TFile("./skimmedFiles/OniaFlowSkim_JpsiTrig_AllPD_isMC0_HFNom_200407.root","read");
  else if(!isMC) rf = new TFile("../primary_input/OniaFlowSkim_JpsiTrig_DBAllPD_isMC0_HFNom_201127.root","read");
//  else if(!isMC) rf = new TFile("OniaFlowSkim_JpsiTrig_JPsi_isMC0_HFNom_200619_RooDataSet_test100k.root","read");
  TTree *tree = (TTree*) rf -> Get("mmepevt");

  
  //Get Correction histograms
  bool isTnP = true;
  bool isPtW = true;
//  TFile *fEff = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Efficiency/mc_eff_vs_pt_TnP%d_PtW1_OfficialMC_Y%dS_muPtCut3.5.root",isTnP,state),"read");
  TFile *fEff1 = new TFile(Form("../primary_input/mc_eff_vs_pt_cent_0_to_20_rap_%s_pbpb_psi2S_PtW%d_tnp%d.root",bCont.Data(),isPtW,isTnP),"read");
  TFile *fEff2 = new TFile(Form("../primary_input/mc_eff_vs_pt_cent_20_to_120_rap_%s_pbpb_psi2S_PtW%d_tnp%d.root",bCont.Data(),isPtW,isTnP),"read");
  // cout << "File Open : " << Form("../mc_eff_vs_pt_cent_0_to_20_rap_%s_pbpb_psi2S_PtW%d_tnp%d.root",bCont.Data(),isPtW,isTnP) << endl;
  TH1D* hEffPt1[15];
  TH1D* hEffPt2[15];
  hEffPt1[0] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy0_1p2",isTnP,isPtW));
  hEffPt1[1] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy1p2_1p6",isTnP,isPtW));
  hEffPt1[2] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy1p6_1p8",isTnP,isPtW));
  hEffPt1[3] = (TH1D*) fEff1 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy1p8_2p4",isTnP,isPtW));
  hEffPt2[0] = (TH1D*) fEff2 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy0_1p2",isTnP,isPtW));
  hEffPt2[1] = (TH1D*) fEff2 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p2_1p6",isTnP,isPtW));
  hEffPt2[2] = (TH1D*) fEff2 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p6_1p8",isTnP,isPtW));
  hEffPt2[3] = (TH1D*) fEff2 -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p8_2p4",isTnP,isPtW));

  // TFile *fEff = new TFile(Form("../primary_input/mc_eff_vs_pt_cent_20_to_120_rap_%s_pbpb_psi2S_PtW%d_tnp%d.root",bCont.Data(),isPtW,isTnP),"read");
  // cout << "File Open : " << Form("../primary_input/mc_eff_vs_pt_cent_20_to_120_rap_%s_pbpb_psi2S_PtW%d_tnp%d.root",bCont.Data(),isPtW,isTnP) << endl;
  // TH1D* hEffPt[15];
  // hEffPt[0] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy0_1p2",isTnP,isPtW));
  // hEffPt[1] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p2_1p6",isTnP,isPtW));
  // hEffPt[2] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p6_1p8",isTnP,isPtW));
  // hEffPt[3] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p8_2p4",isTnP,isPtW));

  // TFile *fAcc = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Acceptance/acceptance_wgt_%dS_pt0_50_20190813_dNdptWeighted.root",state),"read");
  TFile *fAcc = new TFile("../primary_input/acceptance_Prompt_GenOnly_wgt1_2021Psi2Sv2_20210603.root","read");
  TH1D* hAccPt[6];
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
  double SumEff_4;
  double SumEff_5;
  double SumEff_6;

  TBranch *b_event;
  TBranch *b_cBin;
  TBranch *b_nDimu;
  TBranch *b_vz;
  TBranch *b_mass;
  TBranch *b_recoQQsign;
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
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",1.0,6.0,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* evtWeight = new RooRealVar("weight","corr weight", 0, 1000000,"");
  RooRealVar* recoQQ = new RooRealVar("recoQQsign","qq sign",-1,3,"");
  RooRealVar* ctau3DVar = new RooRealVar("ctau3D","ctau3D dimuon",-100000.0,100000.0,"mm");
  RooRealVar* ctau3DErrVar = new RooRealVar("ctau3DErr","ctau3DErr dimuon",-100000.0,100000.0,"mm");
  RooRealVar* ctau3DResVar = new RooRealVar("ctau3DRes","ctau3D Resolution dimuon",-100000.0,100000.0,"mm");
  RooRealVar* ctau3DSignVar = new RooRealVar("ctau3DSign","ctau3D Significance dimuon",-100000.0,100000.0,"mm");
  RooRealVar* NumDimu = new RooRealVar("NumDimu","number of dimuon",0,100,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var,*evtWeight);
  argSet->add(*cBinVar); argSet->add(*recoQQ); argSet->add(*NumDimu); argSet->add(*ctau3DVar); argSet->add(*ctau3DErrVar); argSet->add(*ctau3DResVar);
  
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);


  int nDimuPass=0;
  int nDimu_one=0;
  int nDimu_more=0;
  int nDimu_all=0;

  double weight_acc = 1;
  double weight_eff = 1;
  

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;
  
  double SiMuPtCut = 0;
  SumEff_4=0;
  SumEff_5=0;
  SumEff_6=0;
  int N_4=0;
  int N_5=0;
  int N_6=0;
  //Begin Loop
  for(int i=0; i<nEvt; i++){
	  tree->GetEntry(i);
	  nDimuPass=0;

	  if(fabs(vz)<15) {
		  //Remove double candidate
		  for(int j=0; j<nDimu; j++){
			  //if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && (double)pt1[j]>SiMuPtCut&&(double)pt2[j]>SiMuPtCut&&abs((double)eta1[j])<2.4&&abs((double)eta2[j])<2.4)) continue;
			  if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && IsAcceptanceQQ(pt1[j],eta1[j]) && IsAcceptanceQQ(pt2[j],eta2[j])) ) continue;
			  nDimuPass++;
		  }

		  nDimu_all++;
		  if(nDimuPass>1) {nDimu_more++; continue;}
		  if(nDimuPass==1) nDimu_one++;

		  // Fill Dimuon Loop
		  for(int j=0; j<nDimu; j++){
			  //if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && (double)pt1[j]>SiMuPtCut&&(double)pt2[j]>SiMuPtCut&&abs((double)eta1[j])<2.4&&abs((double)eta2[j])<2.4)) continue;
			  if( ((double)pt[j]<50 && abs((double)y[j])<2.4 && IsAcceptanceQQ(pt1[j],eta1[j])&&IsAcceptanceQQ(pt2[j],eta2[j])) ) { 
			  //if( ((double)pt[j]<50 && abs((double)y[j])<2.4) ) { 
				  weight_acc = 1;
				  weight_eff = 1;
				  if(fAccW){ 
					  if ( abs((double)y[j])<1.6) { weight_acc = getAccWeight(hAccPt[1], pt[j]); }
					  else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<2.4 ) { weight_acc = getAccWeight(hAccPt[2], pt[j]); }
				  }
				  if(fEffW){
					  if(cBin >= 0 && cBin < 20) {
						  if ( abs((double)y[j])<1.2) { weight_eff = getEffWeight(hEffPt1[0], pt[j]); }
						  else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) { weight_eff = getEffWeight(hEffPt1[1], pt[j]); }
						  else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<1.8 ) { weight_eff = getEffWeight(hEffPt1[2], pt[j]); }
						  else if ( abs((double)y[j])>=1.8 && abs((double)y[j])<2.4 ) { weight_eff = getEffWeight(hEffPt1[3], pt[j]); }
					  }else if (cBin >= 20 && cBin < 120) {
						  if ( abs((double)y[j])<1.2) { weight_eff = getEffWeight(hEffPt2[0], pt[j]); }
						  else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) { weight_eff = getEffWeight(hEffPt2[1], pt[j]); }
						  else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<1.8 ) { weight_eff = getEffWeight(hEffPt2[2], pt[j]); }
						  else if ( abs((double)y[j])>=1.8 && abs((double)y[j])<2.4 ) { weight_eff = getEffWeight(hEffPt2[3], pt[j]); }
					  }
				  }
				  //SumEff=+1./weight_eff;

				  double weight_ = weight * weight_eff * weight_acc;
				  if ( pt[j]>6.5 && pt[j]<7 && abs((double)y[j])<0.6 && cBin<180 ) { 
					  SumEff_4+=1./weight_eff; N_4++; 
					  //			if ( 1./weight_eff<=0 ){ cout << "pt :" <<  pt[j] << "		Eff : " << 1./weight_eff << "			Acc : " << 1./weight_acc << "	weight : " << weight_ << endl; }
				  }
				  //		if ( pt[j]>7 && pt[j]<8 && abs((double)y[j])<0.6 && cBin<180 ) { 	SumEff_5+=1./weight_eff; N_5++; }
				  //		if ( pt[j]>8 && pt[j]<8.5 && abs((double)y[j])<0.6 && cBin<180 ) { 	SumEff_6+=1./weight_eff; N_6++; }
				  recoQQ->setVal((int)recoQQsign[j]);     
				  massVar->setVal( (double)mass[j] ) ;
				  ptVar->setVal(   (double)pt[j]   ) ;
				  yVar->setVal(    (double)y[j]    ) ;
				  pt1Var->setVal(  (double)pt1[j]  ) ;
				  eta1Var->setVal( (double)eta1[j] ) ;
				  pt2Var->setVal(  (double)pt2[j]  ) ;
				  eta2Var->setVal( (double)eta2[j] ) ;
				  cBinVar->setVal( (double)cBin ) ;
				  ctau3DVar->setVal( (double)ctau3D[j] ) ;
				  ctau3DErrVar->setVal( (double)ctau3DErr[j] ) ;
				  ctau3DResVar->setVal( (double)ctau3D[j]/ctau3DErr[j] ) ;
				  ctau3DSignVar->setVal( (double)ctau3DErr[j]/ctau3D[j] ) ;
				  evtWeight->setVal( (double)weight_ ) ;
				  NumDimu->setVal((int)nDimu);
				  dataSet->add( *argSet);
			  }
		  }
	  }
  }

 // cout << "Sum of Eff bin 4 : " << SumEff_4  << " N 4 : " << N_4 << endl;
 // cout << "Sum of Eff bin 5 : " << SumEff_5  << " N 5 : " << N_5 << endl;
 // cout << "Sum of Eff bin 6 : " << SumEff_6  << " N 6 : " << N_6 << endl;
 // double avgEff_4=SumEff_4/(double)N_4;
 // double avgEff_5=SumEff_5/(double)N_5;
 // double avgEff_6=SumEff_6/(double)N_6;
 // cout << "Average of Eff bin 4: " << avgEff_4 << endl;
 // cout << "Average of Eff bin 5: " << avgEff_5 << endl;
 // cout << "Average of Eff bin 6: " << avgEff_6 << endl;
 // cout << "Average Eff : " << (avgEff_4+avgEff_5+avgEff_6)/3 << endl;
  cout << "All : " << nDimu_all << endl;
  cout << "more than one dimuon : " << nDimu_more << endl;
  cout << "one dimuon : " << nDimu_one << endl;
  
  
  if (isMC && state==1) {TFile *wf = new TFile(Form("OniaRooDataSet_isMC%d_PR_Psi_2S_20201123.root",isMC),"recreate");  wf->cd();}
  else if (isMC && state==2) {TFile *wf = new TFile(Form("OniaRooDataSet_isMC%d_BtoPsi_2S_20201123.root",isMC),"recreate");  wf->cd();}
  else if (!isMC) {TFile *wf = new TFile(Form("OniaRooDataSet_isMC%d_Psi_2S_%sw_Effw%d_Accw%d_PtW%d_TnP%d_20210604.root",isMC,outName.Data(),fEffW,fAccW,isPtW,isTnP),"recreate");  wf->cd();}
 dataSet->Write();
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

double getAccWeight(TH1D* h, double pt){
  double binN = h->FindBin(pt);
  double weight_ = 1./(h->GetBinContent(binN));
  return weight_;
} 

/*
double getEffWeight(TH1D *h, double pt){
  double binN = h->FindBin(pt);
  TF1 *eff1 = (TF1*)h->FindObject("f1");
  double eff = eff1->Eval(pt);
  double weight_ = 1./eff;
  return weight_;
} 
*/
double getEffWeight(TH1D *h, double pt){
  double binN = h->FindBin(pt);
  double weight_ = 1./(h->GetBinContent(binN));
  return weight_;
} 
