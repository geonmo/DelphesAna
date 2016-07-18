#!/bin/bash
eval `scramv1 runtime -sh`
pwd
echo "This is job number $1"

tar -zxvf delphesProd.tar.gz

#cp DelphesCMSFWLite.cpp  
chmod +x delphesCode/DelphesCMSFWLite 
cp delphesCode/MinBias.pileup .

echo "========= Acquire input filename ========="
#inputfile=`python -c "import PSet; print PSet.process.source.fileNames[0]"`
inputfile=`python -c "import pset_delphes; print pset_delphes.process.source.fileNames[0]"`
echo $inputfile

pfnInput=`edmFileUtil -d $inputfile`
echo $pfnInput
delphesCode/DelphesCMSFWLite delphes_card_CMSPU_mod.tcl Delphes.root $pfnInput

#cmsRun -j FrameworkJobReport.xml -p PSet.py
cmsRun -j FrameworkJobReport.xml -p pset_delphes.py
