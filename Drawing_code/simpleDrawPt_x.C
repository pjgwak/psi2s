void simpleDrawPt_x(){
    gROOT->Macro("./rootlogon.C");

	double xbin[] = {8.25, 30.0};
	double xerr[] = {1.75, 20.0};
	double ybin[] = {0.5493, 0.5017};
	double yerr[] = {0.0379, 0.0058};
    
    // mid rapidity, integrated
    double xbin2[] = {28.25};
    double xerr2[] = {21.75};
    double ybin2[] = {0.5430};
    double yerr2[] = {0.0291};
    
    // forward
    double xbin3[] = {5.25};
    double xerr3[] = {1.25};
    double ybin3[] = {1.6793};
    double yerr3[] = {1.0551};
    
    TH1F *hPad = new TH1F("hPad",";p_{T} (GeV/c);#sigma_{1}/#sigma_{2}", 10, 0.0, 50.0);
    hPad->SetMaximum(5);
    hPad->SetMinimum(0.2);
    hPad->GetXaxis()->CenterTitle();
    hPad->GetYaxis()->CenterTitle();
    
    double average = 0; // average for 0 < |y| < 2.4
    average = (ybin[0] + ybin[1]) / 2;
    
    TGraphErrors *v2pt = new TGraphErrors(2, xbin, ybin, xerr, yerr);
    TGraphErrors *v2pt2 = new TGraphErrors(1, xbin2, ybin2, xerr2, yerr2);
    TGraphErrors *v2pt3 = new TGraphErrors(1, xbin3, ybin3, xerr3, yerr3);

    TCanvas *c1 = new TCanvas("c1","",880,800);
    c1->cd();
    hPad->Draw();

    auto avg_func = new TF1("avg_func", Form("%f", average), 6.5, 50.01);
    avg_func->SetLineStyle(10);
    avg_func->SetLineColor(kGreen+2);
    avg_func->Draw("same");
    
    auto avg_func2 = new TF1("avg_func2", Form("%f", ybin3[0]), 0, 6.5);
    avg_func2->SetLineStyle(2);
    avg_func2->SetLineColor(kGray+2);
    avg_func2->Draw("same");
    
    v2pt->SetMarkerColor(kBlue+1);
    v2pt->SetLineColor(kBlue+1);
    v2pt->SetMarkerSize(1.6);
    v2pt->Draw("pz same");
   
    v2pt2->SetMarkerColor(kRed+1);
    v2pt2->SetLineColor(kRed+1);
    v2pt2->SetMarkerSize(1.6);
    v2pt2->SetMarkerStyle(24);
    v2pt2->Draw("pz same");
    
    v2pt3->SetMarkerColor(kBlack);
    v2pt3->SetLineColor(kBlack);
    v2pt3->SetMarkerSize(1.6);
    v2pt3->Draw("pz same");
    
    TLatex *lt1 = new TLatex();
    lt1->SetNDC();
    lt1->DrawLatex(0.2,0.9,"CMS Preliminary");
    lt1->DrawLatex(0.62,0.9,"PbPb #sqrt{s_{NN}} = 5.02 TeV");
    lt1->DrawLatex(0.2,0.82,"|y| < 2.4");
    lt1->DrawLatex(0.2,0.76,"Cent. 10-60 \%");
    lt1->DrawLatex(0.2,0.70, Form("avg. %.4f", average));
    lt1->DrawLatex(0.2,0.64, Form("int. %.4f", ybin2[0]));

    TLegend *leg1 = new TLegend(0.61,0.70,0.91,0.81);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(v2pt,"MC #psi(2S)");
    leg1->AddEntry(v2pt2, "Integrated");
    leg1->AddEntry(v2pt3, "MC #psi(2S); forward");
    leg1->Draw("same");

    c1->SaveAs("plot_psi2s_mc_pt_x.pdf");
}
