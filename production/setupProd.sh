#!/bin/bash
eval `scramv1 runtime -sh`
wget https://github.com/delphes/delphes/archive/master.zip
unzip master.zip
#tar -zxf Delphes-3.3.2.tar.gz
mv delphes-master delphesCode

cd delphesCode
./configure
cp ../FastJetFinder.cc modules/
make -j 20

cd ..
rm delphesProd.tar.gz
tar -czvf delphesProd.tar.gz delphesCode delphes_card_CMSPU_mod_noTauTagging.tcl
