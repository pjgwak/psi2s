#!/bin/bash

#root -l -b -q 'Psi2Sv2mass_hist.C('0,20,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0475,1')'
#root -l -b -q 'Psi2Sv2mass_hist.C('20,120,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0495,1')'
#root -l -b -q 'Psi2Sv2mass_hist.C('120,200,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0525,1')'
#root -l -b -q 'Psi2Sv2mass_hist.C('20,120,3.0,6.5,1.6,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0685,1')'
#root -l -b -q 'Psi2Sv2mass_hist.C('20,120,6.5,10.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0515,1')'
#root -l -b -q 'Psi2Sv2mass_hist.C('20,120,10.0,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0395,1')'

# root -l -b -q 'Psi2Sv2mass_hist.C('0,20,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0205,0')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0195,0')'
# root -l -b -q 'Psi2Sv2mass_hist.C('120,200,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0195,0')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,3.0,6.5,1.6,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0405,0')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,6.5,10.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0285,0')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,10.0,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0185,0')'
#
#
# root -l -b -q 'Psi2Sv2mass_hist.C('0,20,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0205,1')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0195,1')'
# root -l -b -q 'Psi2Sv2mass_hist.C('120,200,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0195,1')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,3.0,6.5,1.6,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0405,1')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,6.5,10.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0285,1')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,10.0,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0185,1')'
#
# root -l -b -q 'Psi2Sv2mass_hist.C('0,20,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0205,2')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0195,2')'
# root -l -b -q 'Psi2Sv2mass_hist.C('120,200,6.5,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0195,2')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,3.0,6.5,1.6,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0405,2')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,6.5,10.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0285,2')'
# root -l -b -q 'Psi2Sv2mass_hist.C('20,120,10.0,50.0,0.0,2.4,0.0,3.4,4.0,true,false,false,false,2,0.0185,2')'



#									pt      / |y|      / cent.
root -l -b -q 'Psi2S_v2mass_hist_weight_inclusive.C(3.0, 6.5, 1.6, 2.4, 20, 120)'
root -l -b -q 'Psi2S_v2mass_hist_weight_inclusive.C(6.5, 10., 0, 2.4, 20, 120)'
root -l -b -q 'Psi2S_v2mass_hist_weight_inclusive.C(10., 50., 0, 2.4, 20, 120)'
root -l -b -q 'Psi2S_v2mass_hist_weight_inclusive.C(6.5, 50., 0, 2.4, 0, 20)'
root -l -b -q 'Psi2S_v2mass_hist_weight_inclusive.C(6.5, 50., 0, 2.4, 20, 120)'
