void AccPt(){
	gROOT->Macro("rootlogon.C");

	double xbin[] = {4.75, 8.25, 30.0};
	double xerr[] = {0.0, 0.0, 0.0};
	double ybin[] = {0.175, 0.23, 0.139}; // Nominal
    double ybin1[] = {0.136, 0.228, 0.139}; // Acc
	double yerr[] = {0.063, 0.123, 0.072}; // Nominal
    double yerr2[] = {0.063, 0.123, 0.072}; // Acc

	TGraphErrors *v2pt = new TGraphErrors(3, xbin, ybin, xerr, yerr);
	TGraphErrors *v2ptPrp = new TGraphErrors(3, xbin, ybin1, xerr, yerr2);
	TH1F *hPad = new TH1F("hPad",";p_{T} (GeV/c);v_{2}", 10, 0.0, 50.0);
	hPad->SetMaximum(0.5);
	hPad->SetMinimum(-0.15);
	hPad->GetXaxis()->CenterTitle();
	hPad->GetYaxis()->CenterTitle();
    
    TLatex *lt1 = new TLatex();
    lt1->SetNDC();
    
    TCanvas *c1 = new TCanvas("c1","",880,800);
    c1->cd();
    hPad->Draw();
    lt1->DrawLatex(0.2,0.9,"CMS Preliminary");
    lt1->DrawLatex(0.62,0.9,"PbPb #sqrt{s_{NN}} = 5.02 TeV");
    lt1->DrawLatex(0.2,0.82,"|y| < 2.4");
    lt1->DrawLatex(0.2,0.74,"6.5 < p_{T} < 50.0 (GeV/c)");
    v2pt->SetMarkerColor(kBlack);
    v2pt->SetLineColor(kBlack);
    v2pt->SetMarkerSize(1.8);
    v2pt->Draw("pz same");
    
    v2ptPrp->SetMarkerStyle(24);
    v2ptPrp->SetMarkerColor(kRed+2);
    v2ptPrp->SetLineColor(kRed+2);
    v2ptPrp->SetMarkerSize(1.8);
    v2ptPrp->Draw("pz same");

    TLegend *leg1 = new TLegend(0.61,0.70,0.91,0.81);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(v2pt,"Nominal");
    leg1->AddEntry(v2ptPrp,"#splitline{Acceptance}{(w/o p_{T} weight)}");
    leg1->Draw("same");

    c1->SaveAs("acc_pT.png");
}
