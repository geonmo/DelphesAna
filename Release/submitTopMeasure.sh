#!/bin/bash

eval `scramv1 runtime -sh`

echo "Now waiting next job for TTLL"

~/bin/submitCondor -o TTLL -f filelist_TTLL.txt $1
TTLLDIR=`ls -td TTLL*/ | head -n 1 | cut -d'/' -f1`  
cd ${TTLLDIR}
condor_submit run.jdl
cd ..

echo "Now waiting next job for DY10to50"
~/bin/submitCondor -o DY10to50 -f filelist_DY10to50.txt $1
DY1=`ls -td DY10to50*/ | head -n 1 | cut -d'/' -f1`  
cd ${DY1}
condor_submit run.jdl
cd ..

echo "Now waiting next job for DY50"
~/bin/submitCondor -o DY50 -f filelist_DY50.txt $1
DY2=`ls -td DY50*/ | head -n 1 | cut -d'/' -f1`  
cd ${DY2}
condor_submit run.jdl
cd ..

echo "Now waiting next job for SingleTop_tW_top"
~/bin/submitCondor -o ST_tW_top -f filelist_ST_tW_top.txt $1
ST1=`ls -td ST_tW_top*/ | head -n 1 | cut -d'/' -f1`  
cd ${ST1}
condor_submit run.jdl
cd ..

echo "Now waiting next job for SingleTop_tW_antitop"
~/bin/submitCondor -o ST_tW_antitop -f filelist_ST_tW_atop.txt $1
ST2=`ls -td ST_tW_antitop*/ | head -n 1 | cut -d'/' -f1`  
cd ${ST2}
condor_submit run.jdl
cd ..

echo "Now waiting next job for TT mass 166.5"
~/bin/submitCondor -o TT_mass_1665 -f filelist_TT_mass_1665.txt $1
TMASS1=`ls -td TT_mass_1665*/ | head -n 1 | cut -d'/' -f1`  
cd ${TMASS1}
condor_submit run.jdl
cd ..

echo "Now waiting next job for TT mass 169.5"
~/bin/submitCondor -o TT_mass_1695 -f filelist_TT_mass_1695.txt $1
TMASS2=`ls -td TT_mass_1695*/ | head -n 1 | cut -d'/' -f1`  
cd ${TMASS2}
condor_submit run.jdl
cd ..

echo "Now waiting next job for TT mass 175.5"
~/bin/submitCondor -o TT_mass_1755 -f filelist_TT_mass_1755.txt $1
TMASS3=`ls -td TT_mass_1755*/ | head -n 1 | cut -d'/' -f1`  
cd ${TMASS3}
condor_submit run.jdl
cd ..

echo "Now waiting next job for TT mass 178.5"
~/bin/submitCondor -o TT_mass_1785 -f filelist_TT_mass_1785.txt $1
TMASS4=`ls -td TT_mass_1785*/ | head -n 1 | cut -d'/' -f1`  
cd ${TMASS4}
condor_submit run.jdl
cd ..
