#!/bin/bash
eval `scramv1 runtime -sh`

./delphesCode/DelphesCMSFWLite delphes_card_CMSnoPU_noTauTagging.tcl Delphes_noPU.root root://cms-xrdr.sdfarm.kr//xrd//store/user/geonmo/13TeV_miniaod/ttbar_test/miniaod.root 
