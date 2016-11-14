#!/bin/bash
eval `scramv1 runtime -sh`

cut_jpsi='abs(lsv_pid)==443&&(sv_mass-3.09)<0.1'
cut_d0='abs(lsv_pid)==421&&(sv_mass-1.8648)<0.15'
cut_dstar='abs(lsv_pid)==413 && abs(sv_diffmass-0.145)<0.05'

./massPlot.py -c "${cut_jpsi}"  -b [15,5,245] -p 'llsv_mass[0]' -x 'Invariant mass of l+J/#psi[GeV/c^2]' & 
./massPlot.py -c "${cut_d0}"    -b [60,5,245] -p 'llsv_mass[0]' -x 'Invariant mass of l+D^{0}[GeV/c^2]'  &
./massPlot.py -c "${cut_dstar}" -b [20,5,245] -p 'llsv_mass[0]' -x 'Invariant mass of l+D^{*}[GeV/c^2]'   &

