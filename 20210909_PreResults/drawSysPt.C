void drawSysPt(){
  gROOT->Macro("~/rootlogon.C");
  gStyle->SetEndErrorSize(3);
  gStyle->SetPadTopMargin(0.07);

  double systot[] = {0.097467943, 0.051526692, 0.018466185};
  double sysacc[] = {0.039, 0.002, 0.0};
  double syseff[] = {0.037, 0.002, 0.001};
  double systnp[] = {0.052, 0.011, 0.004};
  double syssgn[] = {0.041, 0.0, 0.0};
  double sysbkg[] = {0.0, 0.005, 0.0};
  double sysv2bkg[] = {0.017, 0.01, 0.0};
  double sysbcont[] = {0.044, 0.049, 0.018};

  double sysbin[] = {4.75, 8.25, 30.0};
  double hsysbin[] = {4.0, 6.5, 10.0, 50.0};

  TH1F *hSysTot = new TH1F("hSysTot",";p_{T} (GeV/c);Variations",3, hsysbin);
  TH1F *hSysAcc = new TH1F("hSysAcc",";p_{T} (GeV/c);Systematics",3, hsysbin);
  TH1F *hSysEff = new TH1F("hSysEff",";p_{T} (GeV/c);Systematics",3, hsysbin);
  TH1F *hSysTnp = new TH1F("hSysTnp",";p_{T} (GeV/c);Systematics",3, hsysbin);
  TH1F *hSysSgn = new TH1F("hSysSgn",";p_{T} (GeV/c);Systematics",3, hsysbin);
  TH1F *hSysBkg = new TH1F("hSysBkg",";p_{T} (GeV/c);Systematics",3, hsysbin);
  TH1F *hSysv2Bkg = new TH1F("hSysv2Bkg",";p_{T} (GeV/c);Systematics",3, hsysbin);
  TH1F *hSysbCont = new TH1F("hSysbCont",";p_{T} (GeV/c);Systematics",3, hsysbin);
  for(int i = 0; i < 3; i++){
    hSysTot->SetBinContent(i+1,systot[i]);
    hSysAcc->SetBinContent(i+1,sysacc[i]);
    hSysEff->SetBinContent(i+1,syseff[i]);
    hSysTnp->SetBinContent(i+1,systnp[i]);
    hSysSgn->SetBinContent(i+1,syssgn[i]);
    hSysBkg->SetBinContent(i+1,sysbkg[i]);
    hSysv2Bkg->SetBinContent(i+1,sysv2bkg[i]);
    hSysbCont->SetBinContent(i+1,sysbcont[i]);
  }
  hSysTot->SetMaximum(0.13);
  hSysTot->SetLineColor(kBlack);
  hSysAcc->SetLineColor(kRed);
  hSysEff->SetLineColor(kViolet);
  hSysTnp->SetLineColor(kBlue);
  hSysSgn->SetLineColor(kPink);
  hSysBkg->SetLineColor(kGreen+1);
  hSysv2Bkg->SetLineColor(kMagenta);
  hSysbCont->SetLineColor(kOrange);

  hSysTot->SetLineWidth(6);
  hSysAcc->SetLineWidth(4);
  hSysEff->SetLineWidth(4);
  hSysTnp->SetLineWidth(4);
  hSysSgn->SetLineWidth(4);
  hSysBkg->SetLineWidth(4);
  hSysv2Bkg->SetLineWidth(4);
  hSysbCont->SetLineWidth(4);


  TH1F *hPad = new TH1F("hPad",";p_{T} (GeV/c);v_{2}", 10, 0.0, 50.0);
  hPad->SetMaximum(0.6);
  hPad->SetMinimum(-0.15);
  hPad->GetXaxis()->CenterTitle();
  hPad->GetYaxis()->CenterTitle();

  TCanvas *c1 = new TCanvas("c1","",880,800);
  c1->cd();
  hPad->Draw();

  TLegend *leg1 = new TLegend(0.60,0.50,0.74,0.85);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(hSysTot,"Total Systematics","L");
  leg1->AddEntry(hSysAcc,"Acceptance","L");
  leg1->AddEntry(hSysEff,"Efficiency","L");
  leg1->AddEntry(hSysTnp,"Tag & Probe","L");
  leg1->AddEntry(hSysSgn,"Signal PDF","L");
  leg1->AddEntry(hSysBkg,"Backgournd PDF","L");
  leg1->AddEntry(hSysv2Bkg,"v_{2} backgournd PDF","L");
  leg1->AddEntry(hSysbCont,"b Contamination","L");

  // Only prompt
  c1->cd();
  hSysTot->GetYaxis()->CenterTitle();
  hSysTot->GetXaxis()->CenterTitle();
  hSysTot->Draw();
  hSysAcc->Draw("same");
  hSysEff->Draw("same");
  hSysTnp->Draw("same");
  hSysSgn->Draw("same");
  hSysBkg->Draw("same");
  hSysv2Bkg->Draw("same");
  hSysbCont->Draw("same");
  
  leg1->Draw("same");

  c1->SaveAs("plot_psi2s_v2_sys_pt.png");
}
