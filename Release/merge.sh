#!/bin/bash

eval `scramv1 runtime -sh`
TTLLDIR=`ls -l | grep TTLL | egrep -v filelist | tail -n 1 | awk '{print $9}'`
DY50DIR=`ls -l | grep DY50 | egrep -v filelist | tail -n 1 | awk '{print $9}'`
DY10to50DIR=`ls -l | grep DY10to50 | egrep -v filelist | tail -n 1 | awk '{print $9}'`
echo $TTLLDIR
hadd TTLL.root ${TTLLDIR}/ANA*/ANA*.root
echo $DY50DIR
hadd DY50.root ${DY50DIR}/ANA*/ANA*.root
echo $DY10to50DIR
hadd DY10to50.root ${DY10to50DIR}/ANA*/ANA*.root

