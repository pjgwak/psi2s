void simpleDrawCent(){
	gROOT->Macro("~/rootlogon.C");

	double xbin[] = {5, 35.0};
	double ybin1[] = {0.000, 0.064}; // inclusive
	double ybin2[] = {0.156, 0.165}; // prompt 
	double xerr[] = {0.0, 0.0};
	double yerr1[] = {0.193, 0.025}; // inclusive
	double yerr2[] = {0.054, 0.077}; // prompt

	TGraphErrors *v2cent = new TGraphErrors(2, xbin, ybin1, xerr, yerr1);
	TGraphErrors *v2centPrp = new TGraphErrors(2, xbin, ybin2, xerr, yerr2);
	TH1F *hPad = new TH1F("hPad",";Centrality (%);v_{2}", 10, 0.0, 100.0);
	hPad->GetXaxis()->CenterTitle();
	hPad->GetYaxis()->CenterTitle();
	hPad->SetMaximum(0.35);
	hPad->SetMinimum(-0.08);

	TLatex *lt1 = new TLatex();
	lt1->SetNDC();

	TCanvas *c1 = new TCanvas("c1","",880,800);
	c1->cd();
	hPad->Draw();
	lt1->DrawLatex(0.2,0.9,"CMS Preliminary");
	lt1->DrawLatex(0.62,0.9,"PbPb #sqrt{s_{NN}} = 5.02 TeV");
	lt1->DrawLatex(0.2,0.82,"|y| < 2.4");
	lt1->DrawLatex(0.2,0.74,"6.5 < p_{T} < 50.0 (GeV/c)");
	v2cent->SetMarkerColor(kBlue+1);
	v2cent->SetLineColor(kBlue+1);
	v2cent->SetMarkerSize(1.6);
	v2cent->Draw("pz same");
	v2centPrp->SetMarkerStyle(33);
	v2centPrp->SetMarkerColor(kYellow+2);
	v2centPrp->SetLineColor(kYellow+2);
	v2centPrp->SetMarkerSize(2.2);
	v2centPrp->Draw("pz same");

	TLegend *leg1 = new TLegend(0.61,0.70,0.91,0.81);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(0.03);
	leg1->AddEntry(v2cent,"Inclusive #psi(2S)");
	leg1->AddEntry(v2centPrp,"Prompt #psi(2S)");
	leg1->Draw("same");

	c1->SaveAs("plot_jpsi_v2_cent.png");
}
