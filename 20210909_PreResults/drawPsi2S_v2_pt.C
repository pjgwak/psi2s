void drawPsi2S_v2_pt(){
  gROOT->Macro("~/rootlogon.C");
  gStyle->SetEndErrorSize(3);
  gStyle->SetPadTopMargin(0.07);

  double xlowinc[] = {4.75};
  double xbin[] = {8.25, 30.0};
  double xerr[] = {0.0, 0.0};
  double xlowerr[] = {0.0};
  double xerr2[] = {1.0, 1.0};
  //double xerr2[] = {3.5/2.0, 40.0/2.0};
  //double xerr2[] = {3.5/4.0, 40.0/4.0};
  double xlowerr2[] = {1.0};
  //double xlowerr2[] = {2.5/2.0};

  double yinclow[] = {0.175}; // inclusive
  //double yinclow[] = {0.017}; // inclusive
  double yinc[] = {0.152, 0.111}; // inclusive
  double yinclowerr[] = {0.063}; // inclusive
  //double yinclowerr[] = {0.010}; // inclusive
  double yincerr[] = {0.024, 0.042}; // inclusive

  double xlowprp[] = {5.25};
  double yprplow[] = {0.175}; // prompt
  double yprplowerr[] = {0.063}; // prompt 
  double yprp[] = {0.230, 0.139}; // prompt
  double yprperr[] = {0.123, 0.072}; // prompt 

  double syslowprp[] = {0.097467943}; //0.08697126}; 
  double syshighprp[] = {0.051526692,0.018466185};
  //double syshighprp[] = {0.015937377,0.004123106}; // wo bcont
  //double syshighprp[] = {0.011532563,0.001};

  TGraphErrors *v2ptinclow = new TGraphErrors(1, xlowinc, yinclow, xlowerr, yinclowerr);
  TGraphErrors *v2ptinc = new TGraphErrors(2, xbin, yinc, xerr, yincerr);
  TGraphErrors *v2ptPrplow = new TGraphErrors(1, xlowprp, yprplow, xlowerr, yprplowerr);
  TGraphErrors *v2ptPrp = new TGraphErrors(2, xbin, yprp, xerr, yprperr);
  TGraphErrors *v2ptPrplowsys = new TGraphErrors(1, xlowprp, yprplow, xlowerr2, syslowprp);
  TGraphErrors *v2ptPrpsys = new TGraphErrors(2, xbin, yprp, xerr2, syshighprp);
  TH1F *hPad = new TH1F("hPad",";p_{T} (GeV/c);v_{2}", 10, 0.0, 50.0);
  //hPad->SetMaximum(0.8);
  //hPad->SetMaximum(0.25);
  hPad->SetMaximum(0.6);
  //hPad->SetMinimum(-0.1);
  hPad->SetMinimum(-0.15);
  hPad->GetXaxis()->CenterTitle();
  hPad->GetYaxis()->CenterTitle();

  TCanvas *c1 = new TCanvas("c1","",880,800);
  c1->cd();
  hPad->Draw();
  v2ptinc->SetMarkerColor(kAzure+1);
  v2ptinc->SetLineColor(kAzure+1);
  v2ptinc->SetMarkerSize(1.8);
  v2ptinc->Draw("p same");

  v2ptinclow->SetMarkerColor(kAzure+1);
  v2ptinclow->SetLineColor(kAzure+1);
  v2ptinclow->SetMarkerStyle(24);
  v2ptinclow->SetMarkerSize(1.8);
  v2ptinclow->Draw("p same");

  v2ptPrplow->SetMarkerStyle(27);
  v2ptPrplow->SetMarkerColor(kViolet+2);
  v2ptPrplow->SetLineColor(kViolet+2);
  v2ptPrplow->SetMarkerSize(2.6);
  v2ptPrp->SetMarkerStyle(33);
  v2ptPrp->SetMarkerColor(kViolet+2);
  v2ptPrp->SetLineColor(kViolet+2);
  v2ptPrp->SetMarkerSize(2.6);
  v2ptPrp->Draw("p same");
  v2ptPrplow->Draw("p same");

  TLatex *lt1 = new TLatex(); lt1->SetNDC();
  TLatex *CMS = new TLatex(); CMS->SetNDC();
  TLatex *Pre = new TLatex(); Pre->SetNDC();
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
  lt1->DrawLatex(0.2,0.74,"Cent. 10-60 \%");

  TLegend *leg1 = new TLegend(0.60,0.66,0.74,0.85);
  //TLegend *leg1 = new TLegend(0.55,0.66,0.68,0.85);
  //TLegend *leg1 = new TLegend(0.52,0.61,0.65,0.80);
  //TLegend *leg1 = new TLegend(0.52,0.55,0.65,0.74);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.03);
  //leg1->AddEntry(v2pt,"Inclusive #psi(2S)");
  //leg1->AddEntry(v2pt,"fixed parameter #psi(2S)");
  //leg1->AddEntry(v2ptPrp,"free parameter 0-10 %");
  leg1->AddEntry(v2ptinc,"Inclusive #psi(2S) |y|<2.4");
  leg1->AddEntry(v2ptinclow,"Inclusive #psi(2S) 1.6<|y|<2.4");
  //leg1->AddEntry(v2ptPrp,"Inclusive #psi(2S) free parameter");
  leg1->AddEntry(v2ptPrp,"Prompt #psi(2S) |y|<2.4");
  leg1->AddEntry(v2ptPrplow,"Prompt #psi(2S) 1.6<|y|<2.4");
  leg1->Draw("same");

  //c1->SaveAs("plot_psi2s_v2_pt_fixed.png");
  //c1->SaveAs("plot_psi2s_v2_pt_weighted.png");

  v2ptPrplowsys-> SetFillColorAlpha(kBlue-10,0.3);
  v2ptPrplowsys-> SetLineColor(kBlue-3);
  //v2ptPrplowsys->SetFillStyle(3018);
  //v2ptPrpsys->SetFillStyle(3018);
  v2ptPrpsys-> SetFillColorAlpha(kMagenta-10,0.3);
  v2ptPrpsys-> SetLineColor(kViolet-6);

  // Only prompt
  c1->cd();
  hPad->Draw();
  v2ptPrplowsys->Draw("5"); // A5 or 5
  v2ptPrpsys->Draw("5");
  v2ptPrp->Draw("p same");
  v2ptPrplow->Draw("p same");

  CMS->DrawLatex(0.60,0.87,"CMS");
  Pre->DrawLatex(0.70,0.87,"Preliminary");
  lt1->DrawLatex(0.2,0.87,"Scalar Product");
  lt1->DrawLatex(0.62,0.94,"PbPb 1.7 nb^{-1} (5.02 TeV)");
  lt1->DrawLatex(0.2,0.80,"|y| < 2.4");
  lt1->DrawLatex(0.2,0.74,"Cent. 10-60 \%");

  TLegend *leg2 = new TLegend(0.60,0.66,0.74,0.78);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(v2ptPrp,"Prompt #psi(2S) |y|<2.4");
  leg2->AddEntry(v2ptPrplow,"Prompt #psi(2S) 1.6<|y|<2.4");
  leg2->Draw("same");

  c1->SaveAs("plot_psi2s_v2_prompt_only_pt_weighted.png");
}
