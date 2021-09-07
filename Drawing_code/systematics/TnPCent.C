void TnPCent(){
    gROOT->Macro("./rootlogon.C");

    double xbin[] = {5, 35.0};
    double ybin1[] = {0.108, 0.196}; // Nominal
    double ybin2[] = {0.093, 0.208}; // TnP
    double xerr[] = {0.0, 0.0};
    double yerr1[] = {0.067, 0.111}; // Nominal staticstical
    double yerr2[] = {0.067, 0.111}; // TnP staticstical

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
    v2cent->SetMarkerColor(kBlack);
    v2cent->SetLineColor(kBlack);
    v2cent->SetMarkerSize(1.8);
    v2cent->Draw("pz same");
    
    v2centPrp->SetMarkerStyle(24);
    v2centPrp->SetMarkerColor(kRed+2);
    v2centPrp->SetLineColor(kRed+2);
    v2centPrp->SetMarkerSize(1.8);
    v2centPrp->Draw("pz same");

    TLegend *leg1 = new TLegend(0.61,0.70,0.91,0.81);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(v2cent,"Nominal");
    leg1->AddEntry(v2centPrp,"w/o Tag and Probe");
    leg1->Draw("same");

    c1->SaveAs("TnP_cent.png");
}
