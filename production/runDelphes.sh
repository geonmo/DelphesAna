#!/bin/bash
eval `scramv1 runtime -sh`
pwd
echo "This is job number $1"

echo "========= Acquire input filename ========="

if [ -e PSet.py ]; then
PYCFG='PSet'
else
PYCFG='pset_delphes'
fi
inputfile=`python -c "import ${PYCFG}; print ${PYCFG}.process.source.fileNames[0]"`
echo $inputfile

tar -zxf delphesProd.tar.gz

#cp DelphesCMSFWLite.cpp  
chmod +x delphesCode/DelphesCMSFWLite 
cp delphesCode/MinBias.pileup .



DELPHES_CARD="delphes_card_CMSnoPU_noTauTagging.tcl"

pfnInput=`edmFileUtil -d $inputfile`
echo $pfnInput
delphesCode/DelphesCMSFWLite ${DELPHES_CARD} Delphes.root $pfnInput
jobFin=$?
if [ ${jobFin} -ne 0 ]; then
  rm Delphes.root
  delphesCode/DelphesCMSFWLite ${DELPHES_CARD} Delphes.root root://cmsxrootd.fnal.gov/${inputfile}
  jobFin=$?
fi

cmsRun -j FrameworkJobReport.xml -p ${PYCFG}.py
if [ "$jobFin" -ne 0 ]
then
  echo "delphes job failed"
  exitCode=8028
  exitMessage="delphes simulation failed"
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


