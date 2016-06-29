#!/bin/bash
eval `scramv1 runtime -sh`

echo "This is job number $1"

tar -zxvf delphesProd.tar.gz 
chmod +x Delphes-3.3.2/DelphesCMSFWLite 


echo "========= Acquire input filename ========="
inputfile=`python -c "import PSet; print PSet.process.source.fileNames[0]"`
echo $inputfile

pfnInput=`edmFileUtil -d $inputfile`
echo $pfnInput
Delphes-3.3.2/DelphesCMSFWLite Delphes-3.3.2/cards/delphes_card_CMS_PileUp.tcl Delphes.root $pfnInput

cmsRun -j FrameworkJobReport.xml -p PSet.py
