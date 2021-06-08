#include <iostream>
#include "../header/rootFitHeaders.h"
#include "../header/commonUtility.h"
#include "../header/JpsiUtility.h"
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../header/cutsAndBin.h"
#include "../header/CMS_lumi_v2mass.C"
#include "../header/tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"
using namespace std;
using namespace RooFit;

void MassFit_weight_pt_4p0_6p5_y_1p6_2p4_Cent_20_120_PR_Psi_2S(
		float ptLow=4.0, float ptHigh=6.5,
		float yLow=1.6, float yHigh=2.4,
		int cLow=20, int cHigh=120,
		double sl1_mean = 0.512, double sl2_mean = 0.369, double sl3_mean = 0.32,
//###		double sl1_mean = 0.512, double sl2_mean = 0.369, double sl3_mean = 0.32,
		double N_Jpsi_high =60500, double N_Bkg_high = 1900000,
//###		double N_Jpsi_high =42000, 45000, 46500, 48000, 48500, 51500, 54500, 60500 double N_Bkg_high = 1900000,
		int PR=0, //0=PR, 1=NP, 2=Inc.
		int PRw=1, bool fEffW = true, bool fAccW = true, bool isPtW = true, bool isTnP = true
		)
{
	TString DATE = "210604";
	gStyle->SetEndErrorSize(0);
	gSystem->mkdir(Form("roots/2DFit_%s/Mass",DATE.Data()),kTRUE);
	gSystem->mkdir(Form("figs/2DFit_%s/Mass",DATE.Data()),kTRUE);

	TString bCont;
	if(PR==0) bCont="Prompt";
	else if(PR==1) bCont="NonPrompt";
	else if(PR==2) bCont="Inclusive";

	TString fname;
	if (PRw==1) fname="PR";
	else if (PRw==2) fname="NP";
	
	RooMsgService::instance().getStream(0).removeTopic(Caching);
	RooMsgService::instance().getStream(1).removeTopic(Caching);
	RooMsgService::instance().getStream(0).removeTopic(Plotting);
	RooMsgService::instance().getStream(1).removeTopic(Plotting);
	RooMsgService::instance().getStream(0).removeTopic(Integration);
	RooMsgService::instance().getStream(1).removeTopic(Integration);
	RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
	RooMsgService::instance().getStream(1).removeTopic(Fitting);
	RooMsgService::instance().getStream(1).removeTopic(Minimization);
	RooMsgService::instance().getStream(1).removeTopic(InputArguments);
	RooMsgService::instance().getStream(1).removeTopic(Eval);
	RooMsgService::instance().getStream(1).removeTopic(DataHandling);
	// RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
	RooMsgService::instance().setGlobalKillBelow(ERROR);
	RooMsgService::instance().setSilentMode(true);

	TFile* f1 = new TFile("../make_RooDataSet/roots/OniaRooDataSet_isMC0_pt_4.0_6.5_y_1.6_2.4_Cent_20_120_CtauCtu_0.0395_PR_Psi_2S_20210604.root","read");
	//###TFile* f1 = new TFile(Form("../make_RooDataSet/roots/OniaRooDataSet_isMC0_Psi2S_PRw_Effw1_Accw1_PtW1_TnP1_20210604.root"));
	// TFile* f1 = new TFile(Form("../data/OniaRooDataSet_isMC0_JPsi_%sw_Effw%d_Accw%d_PtW%d_TnP%d_20210111.root",fname.Data(),fEffW,fAccW,isPtW,isTnP));


	// cout << "Input file: "
	// << Form("../data/OniaRooDataSet_isMC0_JPsi_%sw_Effw%d_Accw%d_PtW%d_TnP%d_20210111.root",
	// fname.Data(),fEffW,fAccW,isPtW,isTnP) << endl;
	// exit(0);

	TString kineCut;
	kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>3.4 && mass<4.0 && cBin>=%d && cBin<%d",ptLow, ptHigh, yLow, yHigh, cLow, cHigh);

	TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

	TString OS="recoQQsign==0 &&";

	kineCut = OS+accCut+kineCut;

	RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	ws->data("dataset")->Print();
	cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<", Cent: "<<cLow<<"-"<<cHigh<<"%"<<endl;
	cout << "####################################" << endl;
	RooDataSet *datasetW = new RooDataSet("datasetW","A sample",*dataset->get(),Import(*dataset),WeightVar(*ws->var("weight")));
	RooDataSet *dsAB = (RooDataSet*)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
	cout << "******** New Combined Dataset ***********" << endl;
	dsAB->SetName("dsAB");
	ws->import(*dsAB);
	ws->var("mass")->setRange(massLow, massHigh);
	ws->var("mass")->Print();
	//***********************************************************************
	//****************************** MASS FIT *******************************
	//***********************************************************************

	//         The order is {sigma_1,  x, alpha_1, n_1,   f, m_lambda}
	//pt3-4.5
	//double paramsupper[8] = {0.4,    1.0,     4.9, 3.9, 1.0,     25.0};
	//double paramslower[8] = {0.01,   0.0,     1.1, 1.1, 0.0,    -25.0};//pt3-4.5 m_lambda==-25.0
	//double paramsupper[8] = {0.4,    1.0,     4.9, 2.9, 1.0,     25.0};
	//double paramslower[8] = {0.01,   0.0,     1., 1., 0.0,      0.0};//pt3-4.5 m_lambda==-25.0
	//Cent.10-20
	double paramsupper[8] = {0.4,    1.0,     4.9, 3.9, 1.0,     25.0};
	double paramslower[8] = {0.01,   0.0,     1.1, 1.1, 0.0,      0.0};//pt3-4.5 m_lambda==-25.0
	//SIGNAL: initial params
	double sigma_1_init = 0.04;
	double x_init = 0.35;
	double alpha_1_init = 2.1;
	double n_1_init = 1.4;
	double f_init = 0.4;
	double m_lambda_init = 5;

	double psi_2S_mass = 3.686;

	//SIGNAL
	RooRealVar    mean("m_{J/#Psi}","mean of the signal gaussian mass PDF",psi_2S_mass, psi_2S_mass -0.1, psi_2S_mass + 0.1 ) ;
	RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init, paramslower[1], paramsupper[1]);
	RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init, paramslower[0], paramsupper[0]);
	RooFormulaVar sigma_2_A("sigma_2_A","@0*@1",RooArgList(sigma_1_A, *x_A) );
	RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init, paramslower[2], paramsupper[2]);
	RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0",RooArgList(alpha_1_A) );
	RooRealVar    n_1_A("n_1_A","power order", n_1_init , paramslower[3], paramsupper[3]);
	RooFormulaVar n_2_A("n_2_A","1.0*@0",RooArgList(n_1_A) );
	RooRealVar   *f = new RooRealVar("f","cb fraction", f_init, paramslower[4], paramsupper[4]);

	//Set up crystal ball shapes
	RooCBShape* cb_1_A = new RooCBShape("cball_1_A", "cystal Ball", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
	RooAddPdf*  pdfMASS_Jpsi;

	//DOUBLE CRYSTAL BALL
	RooCBShape* cb_2_A = new RooCBShape("cball_2_A", "cystal Ball", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
	pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
	pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );

	//BACKGROUND
	RooRealVar m_lambda_A("#lambda_A","m_lambda",  m_lambda_init, paramslower[5], paramsupper[5]);
	RooRealVar *sl1 = new RooRealVar("sl1","sl1", sl1_mean, -5., 5.);
	RooRealVar *sl2 = new RooRealVar("sl2","sl2", sl2_mean, -5., 5.);
	RooRealVar *sl3 = new RooRealVar("sl3","sl3", sl3_mean, -5., 5.);

	//THIS IS THE BACKGROUND FUNCTION
	RooChebychev *pdfMASS_bkg;
	pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3));

	//Build the model
	RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0, N_Jpsi_high);
	RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0, N_Bkg_high);
	RooAddPdf* pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *pdfMASS_bkg),RooArgList(*N_Jpsi,*N_Bkg));
	//pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *bkg_1order),RooArgList(*N_Jpsi,*N_Bkg));
	//pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","PR Jpsi + NP Jpsi + Bkg",RooArgList(*cb_1_A, *cb_2_A, *bkg),RooArgList(*N_JpsiPR,*N_JpsiNP,*N_Bkg));
	ws->import(*pdfMASS_Tot);

	TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,4,550,520);
	c_A->cd();
	TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.16, 0.98, 1.0);
	pad_A_1->SetTicks(1,1);
	pad_A_1->Draw(); pad_A_1->cd();
	nMassBin=30;
	RooPlot* myPlot_A = ws->var("mass")->frame(nMassBin); // bins
	myPlot_A->SetTitle("");
	ws->data("dsAB")->plotOn(myPlot_A,Name("dataHist_A"));

	pad_A_1->cd();
	// gPad->SetLogy();
	RooPlot* myPlot2_A = (RooPlot*)myPlot_A->Clone();
	dsAB->plotOn(myPlot2_A,Name("dataOS"),MarkerSize(.8));
	bool isWeighted = ws->data("dsAB")->isWeighted();
	cout << endl << "********* Starting Mass Dist. Fit **************" << endl << endl;
	RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU));
	cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;
	ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_tot"), LineColor(kBlack));
	//ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("Sig_A"),Components(RooArgSet(*pdfMASS_Jpsi)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
	ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_Bkg"),Components(RooArgSet(*pdfMASS_bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

	//make a pretty plot
	myPlot2_A->SetFillStyle(4000);
	myPlot2_A->GetYaxis()->SetTitleOffset(1.43);
	//myPlot2_A->GetYaxis()->CenterTitle();
	//myPlot2_A->GetYaxis()->SetTitleSize(0.058);
	//myPlot2_A->GetYaxis()->SetLabelSize(0.054);
	//myPlot2_A->GetYaxis()->SetRangeUser(ws->var("N_Jpsi")->getVal()/100, ws->var("N_Jpsi")->getVal());
	TH1* h = ws->data("dsAB")->createHistogram("hist", *ws->var("mass"), Binning(myPlot_A->GetNbinsX(),myPlot_A->GetXaxis()->GetXmin(),myPlot_A->GetXaxis()->GetXmax()));
	Double_t YMax = h->GetBinContent(h->GetMaximumBin());
	// Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
	Double_t YMin = 1e99;
	for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
	double Ydown;
	double Yup;
	Ydown = YMin/(TMath::Power((YMax/YMin), (0.1/(1.0-0.1-0.4))));
	Yup = YMax*TMath::Power((YMax/YMin), (0.4/(1.0-0.1-0.4)));
	myPlot2_A->GetYaxis()->SetRangeUser(Ydown,Yup);

	//myPlot2_A->SetMinimum(2*10);
	myPlot2_A->GetXaxis()->SetLabelSize(0);
	myPlot2_A->GetXaxis()->SetTitleSize(0);
	myPlot2_A->GetXaxis()->CenterTitle();
	myPlot2_A->GetXaxis()->SetRangeUser(massLow,massHigh);
	myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
	myPlot2_A->Draw();

	drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
	if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
	else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
	drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
	drawText(Form("n_{J/#psi} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
	drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);

	TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.006, 0.98, 0.227);
	c_A->cd();
	pad_A_2->Draw();
	pad_A_2->cd();
	pad_A_2->SetTopMargin(0); // Upper and lower plot are joined
	pad_A_2->SetBottomMargin(0.67);
	pad_A_2->SetBottomMargin(0.4);
	pad_A_2->SetFillStyle(4000);
	pad_A_2->SetFrameFillStyle(4000);
	pad_A_2->SetTicks(1,1);

	RooPlot* frameTMP = (RooPlot*)myPlot2_A->Clone("TMP");
	RooHist* hpull_A = frameTMP->pullHist("dataOS","pdfMASS_tot", true);
	hpull_A->SetMarkerSize(0.8);
	RooPlot* pullFrame_A = ws->var("mass")->frame(Title("Pull Distribution")) ;
	pullFrame_A->addPlotable(hpull_A,"P") ;
	pullFrame_A->SetTitle("");
	pullFrame_A->SetTitleSize(0);
	pullFrame_A->GetYaxis()->SetTitleOffset(0.3) ;
	pullFrame_A->GetYaxis()->SetTitle("Pull") ;
	pullFrame_A->GetYaxis()->SetTitleSize(0.15) ;
	pullFrame_A->GetYaxis()->SetLabelSize(0.15) ;
	pullFrame_A->GetYaxis()->SetRangeUser(-3.8,3.8);
	pullFrame_A->GetYaxis()->CenterTitle();

	pullFrame_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
	pullFrame_A->GetXaxis()->SetTitleOffset(1.05) ;
	pullFrame_A->GetXaxis()->SetLabelOffset(0.04) ;
	pullFrame_A->GetXaxis()->SetLabelSize(0.15) ;
	pullFrame_A->GetXaxis()->SetTitleSize(0.15) ;
	pullFrame_A->GetXaxis()->CenterTitle();

	pullFrame_A->GetYaxis()->SetTickSize(0.04);
	pullFrame_A->GetYaxis()->SetNdivisions(404);
	pullFrame_A->GetXaxis()->SetTickSize(0.03);
	pullFrame_A->Draw() ;

	TLine *l1 = new TLine(massLow,0,massHigh,0);
	l1->SetLineStyle(1);
	l1->Draw("same");

	printChi2(ws, pad_A_2, frameTMP, fitMass, "mass", "dataOS", "pdfMASS_tot", nMassBin, false);

	TH1D* outh = new TH1D("fitResults","fit result",20,0,20);
	outh->GetXaxis()->SetBinLabel(1,"Jpsi");

	float temp1 = ws->var("N_Jpsi")->getVal();
	float temp1err=ws->var("N_Jpsi")->getError();

	outh->SetBinContent(1,temp1);
	outh->SetBinError(1,temp1err);

	fitMass->Print();
	Double_t theNLL = fitMass->minNll();
	cout << " *** NLL : " << theNLL << endl;
	//RooRealVar *Par1 = ws->var("sl1");
	//RooRealVar *Par2 = ws->var("sl2");
	//cout << "Chebychev Par1 : " << Par1->getVal() << " +/- " << Par1->getError() << endl;
	//cout << "Chebychev Par2 : " << Par2->getVal() << " +/- " << Par2->getError() << endl;

	//ws->Print();
	RooArgSet* fitargs = new RooArgSet();
	fitargs->add(fitMass->floatParsFinal());
	//RooDataSet *datasetMass = new RooDataSet("datasetMass","dataset with Mass Fit result", RooArgList(*ws->var("N_Jpsi"),*ws->var("N_Bkg")) );
	RooDataSet *datasetMass = new RooDataSet("datasetMass","dataset with Mass Fit result", *fitargs );
	datasetMass->add(*fitargs);
	RooWorkspace *wsmass = new RooWorkspace("workspaceMass");

	c_A->Update();
	TString kineLabel = getKineLabel (ptLow, ptHigh,yLow, yHigh, 0.0, cLow, cHigh);

	TFile* outFile;
	//if (PR==2) {
		outFile = new TFile(Form("roots/2DFit_%s/Mass/MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP),"recreate");
		c_A->SaveAs(Form("figs/2DFit_%s/Mass/Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
	//}


	pdfMASS_Tot->Write();
	datasetMass->Write();
	outh->Write();
	// ws->Write();

	outFile->Close();
	fitMass->Print("v");
	// Get # of entries for every bins.
	// To compare events number point by point
	// Due to there was events number issue
	// ofstream fout(Form("entry_check_mass_fit_of_2D_fit_Psi2S_Inclusive_%s.txt", kineLabel.Data()));
	// for (int i = 1; i <= h->GetNbinsX(); i++)
	// fout << h->GetBinContent(i) << endl;
}
