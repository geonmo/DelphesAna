#!/bin/bash

eval `scramv1 runtime -sh`

echo "Now waiting next job for TTLL"
TTLLDIR=`~/bin/submitCondor -o TTLL -f filelist_TTv6.txt measureTop|grep TTLL`  
cd ${TTLLDIR}
condor_submit run.jdl
cd ..

echo "Now waiting next job for DY10to50"
DY1=`~/bin/submitCondor -o DY10to50 -f filelist_DY10to50v1.txt measureTop | grep DY`
cd ${DY1}
condor_submit run.jdl
cd ..

echo "Now waiting next job for DY50"
DY2=`~/bin/submitCondor -o DY50 -f filelist_DY50v1.txt measureTop|grep DY `
cd ${DY2}
condor_submit run.jdl
cd ..
