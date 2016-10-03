import FWCore.ParameterSet.Config as cms
process = cms.Process("Delphes")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr//xrd//store/user/geonmo/13TeV_miniaod/ttbar_test/miniaod1.root']
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
