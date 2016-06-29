from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()




config.General.requestName = 'tutorial_May2015_MC_analysis'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'pset_delphes.py'
config.JobType.inputFiles = ['runDelphes.sh', 'delphesProd.tar.gz']
config.JobType.outputFiles = ['Delphes.root']
config.JobType.scriptExe = 'runDelphes.sh'

config.Data.inputDataset = '/TTTo2L2Nu_13TeV-powheg/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'DelphesProduction_from_ttDilepton'

config.Site.storageSite = 'T3_KR_KISTI'
