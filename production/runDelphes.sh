#!/bin/bash
eval `scramv1 runtime -sh`
pwd
echo "This is job number $1"

tar -zxf delphesProd.tar.gz

#cp DelphesCMSFWLite.cpp  
chmod +x delphesCode/DelphesCMSFWLite 
cp delphesCode/MinBias.pileup .

echo "========= Acquire input filename ========="

if [ -e PSet.py ] 
then
pythonFile="PSet"
#inputfile=`python -c "import PSet; print PSet.process.source.fileNames[0]"`
else
pythonFile="pset_delphes"
fi
inputfile=`python -c "import ${pythonFile}; print ${pythonFile}.process.source.fileNames[0]"`

echo $inputfile

pfnInput=`edmFileUtil -d $inputfile`
echo $pfnInput
delphesCode/DelphesCMSFWLite delphes_card_CMSPU_mod_noTauTagging.tcl Delphes.root $pfnInput
exitCode=$?
echo "Program is running. (${exitCode})"
cmsRun -j FrameworkJobReport.xml -p ${pythonFile}.py

if [ ${exitCode} != 0 ] 
then
echo "Program is failed. (${exitCode})"
exitMessage="Delphes critical failed."
errorType="Fatal Exception"
  if [ -e FrameworkJobReport.xml ]
  then
cat << EOF > FrameworkJobReport.xml.tmp
<FrameworkJobReport>
<FrameworkError ExitStatus="8020" Type="$errorType" >
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


