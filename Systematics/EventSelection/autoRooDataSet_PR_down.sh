#!/bin/bash

#									pt      / |y|      / cent.
# root -l -b -q 'makeRooDataSet_Psi_2S_PR.C(3.0, 6.5, 1.6, 2.4, 20, 120, 0.0405)' bool isTnP = true, bool isPtW = true, bool isPtWAcc = true
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901_down.C(4.0, 6.5, 1.6, 2.4, 20, 120, 0.0395,true,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901_down.C(6.5, 10., 0, 2.4, 20, 120, 0.0285,true,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901_down.C(10., 50., 0, 2.4, 20, 120, 0.0185,true,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901_down.C(6.5, 50., 0, 2.4, 0, 20, 0.0205,true,true,true)'
root -l -b -q 'makeRooDataSet_Psi_2S_PR_20210901_down.C(6.5, 50., 0, 2.4, 20, 120, 0.0195,true,true,true)'
