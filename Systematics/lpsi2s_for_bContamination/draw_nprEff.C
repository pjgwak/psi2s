void draw_nprEff(){
  gROOT->Macro("~/rootlogon.C");
  double prpeff[] = {0.907455, 0.907455, 0.905827};
  double npreff[] = {0.0768952, 0.0785394, 0.0538812};
  double npreff2[] = {0.153331, 0.158044, 0.106735};
  double xptbin[] = {4.75, 8.25, 30.0};

  double prpeffct[] = {0.905339, 0.901534};
  double npreffct[] = {0.0616692, 0.0565486}; 
  double npreffct2[] = {0.12333, 0.11310};

  double xctbin[] = {5, 35};

  TGraph *gPrpEff = new TGraph(3, xptbin, prpeff);
  TGraph *gNPrpEff = new TGraph(3, xptbin, npreff);
  TGraph *gNPrpEff2 = new TGraph(3, xptbin, npreff2);

  gPrpEff->SetMarkerSize(1.6);
  gNPrpEff->SetMarkerSize(1.6);
  gNPrpEff2->SetMarkerSize(1.6);
  gNPrpEff->SetMarkerColor(kBlue);
  gNPrpEff2->SetMarkerColor(kRed);

  TH1F *hPad = new TH1F("hPad",";p_{T} (Gev/C);Efficiency",100,0.0,50.0);

  TLine *line = new TLine(0.0,0.9,50.0,0.9);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  hPad->GetXaxis()->CenterTitle();
  hPad->GetYaxis()->CenterTitle();

  TCanvas *c1 = new TCanvas("c1","",880,800);
  c1->cd();
  hPad->Draw();
  gPrpEff->Draw("p same");
  gNPrpEff->Draw("p same");
  gNPrpEff2->Draw("p same");

	TLegend *leg1 = new TLegend(0.47,0.51,0.67,0.71);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(0.035);
	leg1->AddEntry(gPrpEff,"Prompt #psi(2S)");
	leg1->AddEntry(gNPrpEff,"Nonprompt #psi(2S) default");
	leg1->AddEntry(gNPrpEff2,"Nonprompt #psi(2S) loose cut");
	leg1->Draw("same");


  line->Draw("same");

  c1->SaveAs("plot_nprp_eff.png");


  TGraph *gPrpEffct = new TGraph(2, xctbin, prpeffct);
  TGraph *gNPrpEffct = new TGraph(2, xctbin, npreffct);
  TGraph *gNPrpEffct2 = new TGraph(2, xctbin, npreffct2);

  gPrpEffct->SetMarkerSize(1.6);
  gNPrpEffct->SetMarkerSize(1.6);
  gNPrpEffct2->SetMarkerSize(1.6);
  gNPrpEffct->SetMarkerColor(kBlue);
  gNPrpEffct2->SetMarkerColor(kRed);

  TH1F *hPad1 = new TH1F("hPad1",";Centrality (%);Efficiency",100,0.0,100.0);
  hPad1->GetXaxis()->CenterTitle();
  hPad1->GetYaxis()->CenterTitle();
  c1->cd();
  hPad1->Draw();
  gPrpEffct->Draw("p same");
  gNPrpEffct->Draw("p same");
  gNPrpEffct2->Draw("p same");

	leg1->Draw("same");

  TLine *line1 = new TLine(0.0,0.9,100.0,0.9);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);

  line1->Draw("same");

  c1->SaveAs("plot_nprp_cent_eff.png");

}
