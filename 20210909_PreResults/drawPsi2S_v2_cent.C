void drawPsi2S_v2_cent(){
	gROOT->Macro("~/rootlogon.C");
	gStyle->SetEndErrorSize(3);
	gStyle->SetPadTopMargin(0.07);

	double xbin[] = {5, 35.0};
	double ybin1[] = {0.083, 0.117}; // inclusive
  double yerr1[] = {0.037, 0.043}; // inclusive

  double xerr[] = {0.0, 0.0};
  double xerr2[] = {2, 2};
  //double xerr2[] = {5, 25};
  double ybin2[] = {0.108, 0.196}; // prompt 
  double yerr2[] = {0.067, 0.111}; // prompt
  double yerrzero[] = {0.0, 0.0};
  double prpsys[] = {0.030232433, 0.052839379}; // prompt w/o bcont
  //double prpsys[] = {0.01570, 0.0495}; // prompt

  TGraphErrors *v2cent = new TGraphErrors(2, xbin, ybin1, xerr, yerr1);
  TGraphErrors *v2centPrp = new TGraphErrors(2, xbin, ybin2, xerr, yerr2);
  TGraphErrors *v2centPrpsys = new TGraphErrors(2, xbin, ybin2, xerr2, prpsys);
  TGraphErrors *v2cent1 = new TGraphErrors(2, xbin, ybin1, xerr, yerrzero);
  TGraphErrors *v2centPrp1 = new TGraphErrors(2, xbin, ybin2, xerr, yerrzero);
  TH1F *hPad = new TH1F("hPad",";Centrality (%);v_{2}", 10, 0.0, 100.0);
	hPad->GetXaxis()->CenterTitle();
	hPad->GetYaxis()->CenterTitle();
	hPad->SetMaximum(0.6);
	hPad->SetMinimum(-0.3);


	TCanvas *c1 = new TCanvas("c1","",880,800);
	c1->cd();
	hPad->Draw();

	TLatex *lt1 = new TLatex();lt1->SetNDC();
	TLatex *CMS = new TLatex();CMS->SetNDC();
	TLatex *Pre = new TLatex();Pre->SetNDC();
	//CMS->SetTextSize(0.075*1.55);
	CMS->SetTextSize(0.05);
	CMS->SetTextFont(61);
	Pre->SetTextSize(0.04);
	Pre->SetTextFont(52);
	CMS->DrawLatex(0.60,0.87,"CMS");
	Pre->DrawLatex(0.70,0.87,"Preliminary");
	//CMS->DrawLatex(0.2,0.85,"CMS");
	//Pre->DrawLatex(0.3,0.85,"Preliminary");
	//lt1->DrawLatex(0.2,0.9,"CMS Preliminary");
	lt1->DrawLatex(0.2,0.87,"Scalar Product");
	lt1->DrawLatex(0.62,0.94,"PbPb 1.7 nb^{-1} (5.02 TeV)");
	//lt1->DrawLatex(0.62,0.95,"PbPb #sqrt{s_{NN}} = 5.02 TeV");
	lt1->DrawLatex(0.2,0.80,"|y| < 2.4");
	lt1->DrawLatex(0.2,0.74,"6.5 < p_{T} < 50.0 (GeV/c)");

  v2centPrpsys-> SetFillColorAlpha(kBlue-10,0.3);
  v2centPrpsys-> SetLineColor(kViolet-6);
  v2centPrpsys-> Draw("5");

	v2cent->SetMarkerColor(kAzure+1);
	v2cent->SetLineColor(kAzure+1);
	v2cent->SetMarkerSize(1.6);
	v2cent->Draw("p same");
	v2centPrp->SetMarkerStyle(33);
	v2centPrp->SetMarkerColor(kViolet+1);
	v2centPrp->SetLineColor(kViolet+1);
	v2centPrp->SetMarkerSize(2.6);
	v2centPrp->Draw("p same");


	v2cent1->SetMarkerStyle(24);
	v2cent1->SetMarkerSize(1.7);
	v2cent1->Draw("p same");
	v2centPrp1->SetMarkerStyle(27);
	v2centPrp1->SetMarkerSize(2.6);
	v2centPrp1->Draw("p same");


	TLegend *leg1 = new TLegend(0.60,0.65,0.80,0.77);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(0.035);
	//leg1->AddEntry(v2cent,"Fixed parameter #psi(2S)");
	leg1->AddEntry(v2cent,"Inclusive #psi(2S)");
	//leg1->AddEntry(v2centPrp,"Free parameter #psi(2S)");
	leg1->AddEntry(v2centPrp,"Prompt #psi(2S)");
	leg1->Draw("same");

	//c1->SaveAs("plot_jpsi_v2_cent_fixed.png");
	//c1->SaveAs("plot_psi2s_v2_cent_weighted.png");

  c1->cd();
	hPad->Draw();

	CMS->DrawLatex(0.60,0.87,"CMS");
	Pre->DrawLatex(0.70,0.87,"Preliminary");
	lt1->DrawLatex(0.2,0.87,"Scalar Product");
	lt1->DrawLatex(0.62,0.94,"PbPb 1.7 nb^{-1} (5.02 TeV)");
	lt1->DrawLatex(0.2,0.80,"|y| < 2.4");
	lt1->DrawLatex(0.2,0.74,"6.5 < p_{T} < 50.0 (GeV/c)");

  v2centPrpsys-> Draw("5");
	v2centPrp->Draw("p same");
	v2centPrp1->Draw("p same");

	TLegend *leg2 = new TLegend(0.60,0.75,0.80,0.77);
	leg2->SetBorderSize(0);
	leg2->SetFillColor(0);
	leg2->SetFillStyle(0);
	leg2->SetTextSize(0.035);
	leg2->AddEntry(v2centPrp,"Prompt #psi(2S)");
	leg2->Draw("same");

	c1->SaveAs("plot_psi2s_v2_prompt_only_cent_weighted.png");

}
