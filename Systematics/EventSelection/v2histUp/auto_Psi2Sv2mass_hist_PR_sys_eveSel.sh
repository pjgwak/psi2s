#!/bin/bash
#								  pt        /|y|       /cent.      /ctau cut


# root -l -b -q 'Psi2S_v2mass_hist_weight_PR.C(3.0, 6.5, 1.6, 2.4, 20, 120, 0, 0.0405)'
# root -l -b -q 'Psi2S_v2mass_hist_weight_PR.C(4.0, 6.5, 1.6, 2.4, 20, 120, 0, 0.0395)'
# root -l -b -q 'Psi2S_v2mass_hist_weight_PR.C(4.5, 6.5, 1.6, 2.4, 20, 120, 0, 0.0395)'
# root -l -b -q 'Psi2S_v2mass_hist_weight_PR.C(6.5, 10., 0, 2.4, 20, 120, 0, 0.0285)'
# root -l -b -q 'Psi2S_v2mass_hist_weight_PR.C(10., 50., 0, 2.4, 20, 120, 0, 0.0185)'
# root -l -b -q 'Psi2S_v2mass_hist_weight_PR.C(6.5, 50., 0, 2.4, 0, 20, 0, 0.0205)'
# root -l -b -q 'Psi2S_v2mass_hist_weight_PR.C(6.5, 50., 0, 2.4, 20, 120, 0, 0.0195)'
# root -l -b -q 'Psi2S_v2mass_hist_weight_PR.C(6.5, 30., 0, 2.4, 20, 120, 0, 0.0215)'


root -l -b -q 'PR_pT4_65_y16_24_c20_120.C(true,true,true)'
root -l -b -q 'PR_pT10_50_y0_24_c20_120.C(true,true,true)'
root -l -b -q 'PR_pT65_10_y0_24_c20_120.C(true,true,true)'
root -l -b -q 'PR_pT65_30_y0_24_c20_120.C(true,true,true)'
root -l -b -q 'PR_pT65_50_y0_24_c0_20.C(true,true,true)'
root -l -b -q 'PR_pT65_50_y0_24_c20_120.C(true,true,true)'


