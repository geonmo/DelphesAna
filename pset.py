import FWCore.ParameterSet.Config as cms
process = cms.Process("Delphes")
#process.source = cms.Source("EmptySource")
print process.source.fileNames
