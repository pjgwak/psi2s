#!/bin/bash

#									pt      / |y|      / cent.
# root -l -b -q 'makeRooDataSet_Psi_2S_PR.C(3.0, 6.5, 1.6, 2.4, 20, 120, 0.0405)' bool isTnP = true, bool isPtW = true, bool isPtWAcc = true

root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210831.C(4.0, 6.5, 1.6, 2.4, 20, 120, 0.0395,false,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210831.C(4.5, 6.5, 1.6, 2.4, 20, 120, 0.0395,false,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210831.C(6.5, 10., 0, 2.4, 20, 120, 0.0285,false,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210831.C(10., 50., 0, 2.4, 20, 120, 0.0185,false,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210831.C(6.5, 50., 0, 2.4, 0, 20, 0.0205,false,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210831.C(6.5, 50., 0, 2.4, 20, 120, 0.0195,false,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210831.C(6.5, 30., 0, 2.4, 20, 120, 0.0215,false,true,true)'
