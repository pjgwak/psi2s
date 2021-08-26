void simpleDrawCent_x(){
	gROOT->Macro("./rootlogon.C");

	double xbin[] = {5, 35};
	double ybin1[] = {0.5717, 0.5430};
	double xerr[] = {5.0, 25.0};
	double yerr1[] = {0.2070, 0.0291};
    
    double xbin2[] = {30};
    double xerr2[] = {30};
    double ybin2[] = {0.5626};
    double yerr2[] = {0.0470};

    double average = 0;
    average = (ybin1[0] + ybin1[1]) / 2;
    
	TGraphErrors *v2cent = new TGraphErrors(2, xbin, ybin1, xerr, yerr1);
    TGraphErrors *v2cent2 = new TGraphErrors(2, xbin2, ybin2, xerr2, yerr2);
	TH1F *hPad = new TH1F("hPad",";Centrality (%);#sigma_{1}/#sigma_{2}", 10, 0.0, 100.0);
	hPad->GetXaxis()->CenterTitle();
	hPad->GetYaxis()->CenterTitle();
	hPad->SetMaximum(1.1);
	hPad->SetMinimum(0.2);

	TLatex *lt1 = new TLatex();
	lt1->SetNDC();

	TCanvas *c1 = new TCanvas("c1","",880,800);
	c1->cd();
	hPad->Draw();
    
    auto avg_func = new TF1("avg_func", Form("%f", average), -0.01, 100.01);
    avg_func->SetLineStyle(10);
    avg_func->SetLineColor(kGreen+2);
    avg_func->Draw("same");
    
    lt1->DrawLatex(0.2,0.9,"CMS Preliminary");
    lt1->DrawLatex(0.62,0.9,"PbPb #sqrt{s_{NN}} = 5.02 TeV");
    lt1->DrawLatex(0.2,0.82,"|y| < 2.4");
    lt1->DrawLatex(0.2,0.76,"Cent. 10-60 \%");
    lt1->DrawLatex(0.2,0.70, Form("avg. %.4f", average));
    lt1->DrawLatex(0.2,0.64, Form("int. %.4f", ybin2[0]));
    
	v2cent->SetMarkerColor(kBlue+1);
	v2cent->SetLineColor(kBlue+1);
	v2cent->SetMarkerSize(1.6);
	v2cent->Draw("pz same");
    
    v2cent2->SetMarkerColor(kRed+1);
    v2cent2->SetLineColor(kRed+1);
    v2cent2->SetMarkerSize(1.6);
    v2cent2->SetMarkerStyle(24);
    v2cent2->Draw("pz same");

	TLegend *leg1 = new TLegend(0.61,0.70,0.91,0.81);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(0.03);
	leg1->AddEntry(v2cent,"MC #psi(2S)");
    leg1->AddEntry(v2cent2, "Integrated");
	leg1->Draw("same");

	c1->SaveAs("plot_jpsi_mc_cent_x.pdf");
}
