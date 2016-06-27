#!/bin/bash
eval `scramv1 runtime -sh`

echo "This is job number $1"

tar -zxvf delphesProd.tar.gz 
chmod +x Delphes3.3.2/DelphesCMSFWLite 


echo "========= Acquire input filename ========="
inputfile=`python -c "import PSet; print PSet.Process.source.fileNames[0]"`
echo $inputfile

pfnInput=`edmFileUtil -d $inputfile`

Delphes3.3.2/DelphesCMSFWLite Delphes-3.3.2/cards/delphes_card_CMS_PileUp.tcl Delphes.root $inputfile 

cmsRun -j FrameworkJobReport.xml -p PSet.py
