#!/bin/bash
eval `scramv1 runtime -sh`
if [ ! -e master.zip ]; then
  wget https://github.com/delphes/delphes/archive/master.zip
fi
unzip master.zip
#tar -zxf Delphes-3.3.2.tar.gz
mv delphes-master delphesCode

cd delphesCode
./configure
cp ../FastJetFinder.* modules/
cp ../DelphesCMSFWLite.cpp readers/
make -j 20

cd ..
rm delphesProd.tar.gz
tar -czvf delphesProd.tar.gz delphesCode delphes_card_CMS*.tcl
