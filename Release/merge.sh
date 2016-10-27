#!/bin/bash

eval `scramv1 runtime -sh`
TTLLDIR=`ls -td  TTLL*/ | head -n 1`
DY1=`ls -td DY10to50*/ | head -n 1 | cut -d'/' -f1`  
DY2=`ls -td DY50*/ | head -n 1 | cut -d'/' -f1`  
ST1=`ls -td ST_tW_top*/ | head -n 1 | cut -d'/' -f1`  
ST2=`ls -td ST_tW_antitop*/ | head -n 1 | cut -d'/' -f1`  
TMASS1=`ls -td TT_mass_1665*/ | head -n 1 | cut -d'/' -f1`  
TMASS2=`ls -td TT_mass_1695*/ | head -n 1 | cut -d'/' -f1`  
TMASS3=`ls -td TT_mass_1755*/ | head -n 1 | cut -d'/' -f1`  
TMASS4=`ls -td TT_mass_1785*/ | head -n 1 | cut -d'/' -f1`  


echo ${TTLLDIR}
hadd hisOut_TTLL_powheg.root ${TTLLDIR}/ANA*/ANA*.root
echo ${DY1}
hadd histOut_DYJets_10to50.root ${DY1}/ANA*/ANA*.root
echo ${DY2}
hadd histOut_DYJets.root ${DY2}/ANA*/ANA*.root
echo ${ST1}
hadd histOut_SingleTop_tW.root ${ST1}/ANA*/ANA*.root
echo ${ST2}
hadd histOut_SingleTbar_tW.root ${ST2}/ANA*/ANA*.root
echo ${TMASS1}
hadd histOut_TT_powheg_mtop1665.root ${TMASS1}/ANA*/ANA*.root
echo ${TMASS2}
hadd histOut_TT_powheg_mtop1695.root ${TMASS2}/ANA*/ANA*.root
echo ${TMASS3}
hadd histOut_TT_powheg_mtop1755.root ${TMASS3}/ANA*/ANA*.root
echo ${TMASS4}
hadd histOut_TT_powheg_mtop1785.root ${TMASS4}/ANA*/ANA*.root
