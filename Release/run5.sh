#!/bin/bash
eval `scramv1 runtime -sh`


jpsi_cut='abs(sv_pid)==443'
d0_cut='abs(sv_pid)==421'
#d0_cut='abs(sv_pid)==421'
#dstar_cut='abs(sv_pid)==413&&sv_dau_pt[0]>4&&sv_dau_pt[1]>4&&sv_dau_pt[2]>4'
dstar_cut='abs(sv_pid)==413'

./massPlot.py  -c ${jpsi_cut}  -b '[50,2.5,3.5]' -p sv_mass -x 'Mass of J/psi Cands [GeV/c^2]' &
./massPlot.py  -c ${d0_cut}    -b '[40,1.3,2.5]' -p sv_mass -x 'Mass of D0 Cands [GeV/c^2] ' &
#./massPlot.py  -c ${dstar_cut} -b [60,1.9,2.1] -p sv_mass -x 'Mass of Dstar Cands [GeV/c^2] ' & 
./massPlot.py  -c ${dstar_cut} -b '[35,0.135,0.17]' -p sv_diffmass -x 'M_{K #pi #pi}-M_{K #pi} [GeV/c^2] ' & 


