from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True

import os
index = 1
baseName = 'DelphesProduction'
# if it finds crab_UserScriptTest1/, it will use crab_UserScriptTest2/ automatically!
while os.path.isdir("crab_%s%s" % (baseName, index)):
    index += 1
config.General.requestName = baseName + str(index)

print config.General.requestName

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'pset.py'
config.JobType.inputFiles = ['runDelphes.sh', 'delphesProd.tar.gz']
config.JobType.outputFiles = ['Delphes.root']
config.JobType.scriptExe = 'runDelphes.sh'

config.section_('Data')
config.Data.inputDataset = '/TTTo2L2Nu_13TeV-powheg/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/MINIAODSIM'
config.Data.publication = False
config.Data.totalUnits = 1
config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.ignoreLocality = True

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_KR_KISTI'
