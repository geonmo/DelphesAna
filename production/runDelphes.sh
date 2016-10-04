#!/bin/bash
eval `scramv1 runtime -sh`
pwd
echo "This is job number $1"

echo "========= Acquire input filename ========="
inputfile=`python -c "import PSet; print PSet.process.source.fileNames[0]"`
#inputfile=`python -c "import pset_delphes; print pset_delphes.process.source.fileNames[0]"`
echo $inputfile

cp delphes_card_CMSPU_mod.tcl $1

mkdir $1
cd $1

tar -zxf delphesProd.tar.gz

#cp DelphesCMSFWLite.cpp  
chmod +x delphesCode/DelphesCMSFWLite 
cp delphesCode/MinBias.pileup .


pfnInput=`edmFileUtil -d $inputfile`
echo $pfnInput
delphesCode/DelphesCMSFWLite delphes_card_CMSPU_mod.tcl Delphes.root $pfnInput
jobFin=$?

if [ "$jobFin" -ne 0 ]
then
  echo "delphes job failed"
  exitCode=8028
  exitMessage="Something is wrong during delphes simulation"
  errorType="File open error or logical error"
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
fi


cd ..
cmsRun -j FrameworkJobReport.xml -p PSet.py
#cmsRun -j FrameworkJobReport.xml -p pset_delphes.py
cp $1/Delphes.root .
rm -rf $1
