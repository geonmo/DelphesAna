#!/bin/bash
eval `scramv1 runtime -sh`
wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.3.2.tar.gz
tar -zxf Delphes-3.3.2.tar.gz

cp DelphesCMSFWLite.cpp Delphes-3.3.2/reader 

cd Delphes-3.3.2
make -j 20

cd -
tar -czvf delphesProd.tar.gz Delphes-3.3.2/