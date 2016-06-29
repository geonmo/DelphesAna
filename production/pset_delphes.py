import FWCore.ParameterSet.Config as cms
process = cms.Process("Delphes")
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.source.fileNames = ['/store/user/geonmo/TTLL13TeV_AODSIM/AODSIM.root']
