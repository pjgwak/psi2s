#!/bin/bash
#						 pt	        |y|     cent.  /BKG func paras   /Max # of Jpsi, Bkg 
# root -l -b -q 'MassFit.C(3.0, 6.5, 1.6, 2.4, 20, 120, 0.512, 0.695, 0.32, 5000, 200000)'
# root -l -b -q 'MassFit.C(6.5, 10, 0, 2.4, 20, 120, 0.512, 0.369, 0.32, 2290, 50000)'
# root -l -b -q 'MassFit.C(6.5, 50.0, 0, 2.4, 0, 20, 0.512, 0.8, 0.39, 4000, 50000)'
# root -l -b -q 'MassFit.C(6.5, 50.0, 0, 2.4, 20, 120, 0.512, 0.8, 0.39, 5000, 100000)'
# root -l -b -q 'MassFit.C(10, 50, 0, 2.4, 20, 120, 0.2, 0.3, 0.1, 5000, 200000)'


# Weight
root -l -b -q 'MassFit_weight.C(3.0, 6.5, 1.6, 2.4, 20, 120, 0.3, 0.2, 0.2, 58000, 6000000)'
root -l -b -q 'MassFit_weight.C(6.5, 10, 0, 2.4, 20, 120, 0.512, 0.369, 0.32, 20000, 700000)'
root -l -b -q 'MassFit_weight.C(6.5, 50.0, 0, 2.4, 0, 20, 0.0012, 0.004, 0.0039, 27000, 700000)'
root -l -b -q 'MassFit_weight.C(6.5, 50.0, 0, 2.4, 20, 120, 0.002, 0.03042, 0.0351, 35000, 650000)'
root -l -b -q 'MassFit_weight.C(10, 50, 0, 2.4, 20, 120, 0.3, 0.3, 0.3, 6500, 50000)'


# Preveious bins
# root -l -b -q 'MassFit.C(6.5, 50.0, 0, 2.4, 20, 180, 0.2, 0.2, 0.2, 6500, 100000)'
#
# root -l -b -q 'MassFit.C(6.5, 50.0, 1.6, 2.4, 20, 180, 0.2, 0.2, 0.2, 3000, 100000)'
# root -l -b -q 'MassFit.C(4.5, 50.0, 1.6, 2.4, 20, 180, 0.512, 0.695, 0.32, 5000, 1000000)'
#
# root -l -b -q 'MassFit.C(3.0, 50.0, 1.6, 2.4, 0, 20 ,-0.3, -0.1, -0.1, 1500, 100000)'
# root -l -b -q 'MassFit.C(3.0, 50.0, 1.6, 2.4, 20, 120, 0.312, 0.195, 0.32, 4900, 200000)'
# root -l -b -q 'MassFit.C(3.0, 50.0, 1.6, 2.4, 20, 180, 0.512, 0.695, 0.32, 5000, 200000)'
#
# root -l -b -q 'MassFit.C(3.0, 6.5, 1.6, 2.4, 0, 20, 0.3, 0.2, 0.3, 400, 150000)'






