from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()


config.General.workArea = 'crab_projects'
import os
index = 1
baseName = 'DelphesProduction'
# if it finds crab_UserScriptTest1/, it will use crab_UserScriptTest2/ automatically!
while os.path.isdir("%s/crab_%s%s" % (config.General.workArea, baseName, index)):
    index += 1
config.General.requestName = baseName + str(index)

config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'pset_delphes.py'
config.JobType.inputFiles = ['runDelphes.sh', 'delphesProd.tar.gz']
config.JobType.outputFiles = ['Delphes.root']
config.JobType.scriptExe = 'runDelphes.sh'

### TT to Dilepton
config.Data.inputDataset = '/TTTo2L2Nu_13TeV-powheg/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/AODSIM'

### TT to Lepton+Jet
#config.Data.inputDataset = '/TTToSemiLeptonic_13TeV-powheg/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/AODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 1
config.Data.outLFNDirBase = '/store/user/%s/DelphesAna_201607/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'DelphesProduction_%d'%(index)

config.Site.storageSite = 'T3_KR_UOS'
