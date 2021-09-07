void simpleDrawPt(){
	gROOT->Macro("rootlogon.C");

	double xbin[] = {4.75, 8.25, 30.0};
	double xerr[] = {0.0, 0.0, 0.0};
	double ybin[] = {0, 0, 0}; // inclusive
	double yerr[] = {0, 0, 0}; // inclusive
	double ybin1[] = {0.154, 0.230, 0.139}; // prompt
	double yerr2[] = {0.150, 0.118, 0.072}; // prompt

	TGraphErrors *v2pt = new TGraphErrors(3, xbin, ybin, xerr, yerr);
	TGraphErrors *v2ptPrp = new TGraphErrors(3, xbin, ybin1, xerr, yerr2);
	TH1F *hPad = new TH1F("hPad",";p_{T} (GeV/c);v_{2}", 10, 0.0, 50.0);
	hPad->SetMaximum(0.5);
	hPad->SetMinimum(-0.15);
	hPad->GetXaxis()->CenterTitle();
	hPad->GetYaxis()->CenterTitle();

	TCanvas *c1 = new TCanvas("c1","",880,800);
	c1->cd();
	hPad->Draw();
	v2pt->SetMarkerColor(kBlue+1);
	v2pt->SetLineColor(kBlue+1);
	v2pt->SetMarkerSize(1.6);
	v2pt->Draw("pz same");

	v2ptPrp->SetMarkerStyle(33);
	v2ptPrp->SetMarkerColor(kYellow+2);
	v2ptPrp->SetLineColor(kYellow+2);
	v2ptPrp->SetMarkerSize(2.2);
	v2ptPrp->Draw("pz same");

	TLatex *lt1 = new TLatex();
	lt1->SetNDC();
	lt1->DrawLatex(0.2,0.9,"CMS Preliminary");
	lt1->DrawLatex(0.62,0.9,"PbPb #sqrt{s_{NN}} = 5.02 TeV");
	lt1->DrawLatex(0.2,0.82,"|y| < 2.4");
	lt1->DrawLatex(0.2,0.74,"Cent. 10-60 \%");

	TLegend *leg1 = new TLegend(0.61,0.70,0.91,0.81);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(0.03);
	leg1->AddEntry(v2pt,"Inclusive #psi(2S)");
	leg1->AddEntry(v2ptPrp,"Prompt #psi(2S)");
	leg1->Draw("same");

	c1->SaveAs("plot_psi2s_v2_pt.png");
}
