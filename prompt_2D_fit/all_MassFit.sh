#!/bin/bash
#						 pt	        |y|     cent.  /BKG func paras   /Max # of Jpsi, Bkg 


# root -l -b -q MassFit_weight_pt_10p0_50p0_y_0p0_2p4_Cent_20_120.C
# root -l -b -q MassFit_weight_pt_4p0_6p5_y_1p6_2p4_Cent_20_120.C
# root -l -b -q MassFit_weight_pt_6p5_10p0_y_0p0_2p4_Cent_20_120.C
# root -l -b -q MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_0_20.C
# root -l -b -q MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_20_120.C

# Weight
# root -l -b -q 'MassFit_weight_pt_4p0_6p5_y_1p6_2p4_Cent_20_120_PR_Psi_2S.C(4.0, 6.5, 1.6, 2.4, 20, 120, 0.312, 0.369, 0.32, 24000, 2000000)'
# root -l -b -q 'MassFit_weight_pt_4p5_6p5_y_1p6_2p4_Cent_20_120_PR_Psi_2S.C(4.5, 6.5, 1.6, 2.4, 20, 120, 0.312, 0.369, 0.32, 18000, 1000000)'
# root -l -b -q 'MassFit_weight_pt_6p5_10p0_y_0p0_2p4_Cent_20_120_PR_Psi_2S.C(6.5, 10, 0, 2.4, 20, 120, 0.212, 0.49, 0.32, 8500, 800000)'
# root -l -b -q 'MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_0_20_PR_Psi_2S.C(6.5, 50.0, 0, 2.4, 0, 20, 0.212, 0.2, 0.09, 9000, 500000)'
# root -l -b -q 'MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_20_120_PR_Psi_2S.C(6.5, 50.0, 0, 2.4, 20, 120, 0.412, 0.2, 0.09, 8500, 800000)'
# root -l -b -q 'MassFit_weight_pt_10p0_50p0_y_0p0_2p4_Cent_20_120_PR_Psi_2S.C(10, 50, 0, 2.4, 20, 120, 0.042, 0.079, 0.4, 1400, 22000)'

# Signal systematic study
root -l -b -q signal_systematic_MassFit_weight_pt_10p0_50p0_y_0p0_2p4_Cent_20_120.C
root -l -b -q signal_systematic_MassFit_weight_pt_4p0_6p5_y_1p6_2p4_Cent_20_120.C
root -l -b -q signal_systematic_MassFit_weight_pt_6p5_10p0_y_0p0_2p4_Cent_20_120.C
root -l -b -q signal_systematic_MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_0_20.C
root -l -b -q signal_systematic_MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_20_120.C


# Background systematic study
# root -l -b -q bkg_systematic_MassFit_weight_pt_10p0_50p0_y_0p0_2p4_Cent_20_120.C
# root -l -b -q bkg_systematic_MassFit_weight_pt_4p0_6p5_y_1p6_2p4_Cent_20_120.C
# root -l -b -q bkg_systematic_MassFit_weight_pt_6p5_10p0_y_0p0_2p4_Cent_20_120.C
# root -l -b -q bkg_systematic_MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_0_20.C
# root -l -b -q bkg_systematic_MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_20_120.C

# MC prompt
# root -l -b -q mc_MassFit_weight_pt_4p0_6p5_y_1p6_2p4_Cent_20_120.C
# root -l -b -q mc_MassFit_weight_pt_6p5_10p0_y_0p0_2p4_Cent_20_120.C
# root -l -b -q mc_MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_0_20.C
# root -l -b -q mc_MassFit_weight_pt_6p5_50p0_y_0p0_2p4_Cent_20_120.C
# root -l -b -q mc_MassFit_weight_pt_10p0_50p0_y_0p0_2p4_Cent_20_120.C
