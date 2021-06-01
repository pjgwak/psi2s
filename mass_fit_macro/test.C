void test()
{
	TFile f("MassFitPlot_20210524_pt3.0-6.5_y1.6-2.4_muPt0.0_centrality0-20w_Effw0_Accw0_PtW1_TnP1.root") ;
	RooWorkspace* w = (RooWorkspace*)f.Get("wsMassPlot");
 w->Draw("myPlot_A") ;
}
