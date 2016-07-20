import FWCore.ParameterSet.Config as cms
process = cms.Process("Delphes")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['file:/pnfs/user/geonmo/CMSSW_8_1_0_pre8/src/DelphesAna/production/00491874-02FA-E511-BCE8-008CFA197454.root']
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
