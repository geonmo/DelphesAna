#!/bin/bash
eval `scramv1 runtime -sh`
pwd
echo "This is job number $1"

tar -zxf delphesProd.tar.gz

#cp DelphesCMSFWLite.cpp  
chmod +x delphesCode/DelphesCMSFWLite 
cp delphesCode/MinBias.pileup .

echo "========= Acquire input filename ========="
#inputfile=`python -c "import PSet; print PSet.process.source.fileNames[0]"`
inputfile=`python -c "import pset_delphes; print pset_delphes.process.source.fileNames[0]"`
echo $inputfile

pfnInput=`edmFileUtil -d $inputfile`
echo $pfnInput
delphesCode/DelphesCMSFWLite delphes_card_CMSPU_mod_noTauTagging.tcl Delphes.root $pfnInput

if [ $? != 0 ] 
then
exitCode=$?
exitMessage="My arbitrary exit message"
errorType="My arbitrary error type"
  if [ -e FrameworkJobReport.xml ]
  then
cat << EOF > FrameworkJobReport.xml.tmp
<FrameworkJobReport>
<FrameworkError ExitStatus="$exitCode" Type="$errorType" >
$exitMessage
</FrameworkError>
EOF
tail -n+2 FrameworkJobReport.xml >> FrameworkJobReport.xml.tmp
mv FrameworkJobReport.xml.tmp FrameworkJobReport.xml
  else
cat << EOF > FrameworkJobReport.xml
<FrameworkJobReport>
<FrameworkError ExitStatus="$exitCode" Type="$errorType" >
$exitMessage
</FrameworkError>
</FrameworkJobReport>
EOF
fi
else 
echo "Successfully prod."
fi


#cmsRun -j FrameworkJobReport.xml -p PSet.py
cmsRun -j FrameworkJobReport.xml -p pset_delphes.py
