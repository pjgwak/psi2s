#!/bin/bash

#									                                   / pt      / |y|     / cent.  / psi cut, isPtW(Eff), isPtWAcc(Acc), isTnP(TnP)
# root -l -q -b 'makeRooDataSet_Psi_2S_PR_20210901.C(4.0, 6.5, 1.6, 2.4, 20, 120, 0.0395,true,false,true)'

root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901.C(4.0, 6.5, 1.6, 2.4, 20, 120, 0.0395,true,false,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901.C(4.5, 6.5, 1.6, 2.4, 20, 120, 0.0395,true,false,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901.C(6.5, 10., 0, 2.4, 20, 120, 0.0285,true,false,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901.C(10., 50., 0, 2.4, 20, 120, 0.0185,true,false,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901.C(6.5, 50., 0, 2.4, 0, 20, 0.0205,true,false,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901.C(6.5, 50., 0, 2.4, 20, 120, 0.0195,true,false,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901.C(6.5, 30., 0, 2.4, 20, 120, 0.0215,true,false,true)'

