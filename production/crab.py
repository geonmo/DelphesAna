from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True

import os
index = 1
baseName = 'UserScriptTest'
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
config.Data.publication = False
config.Data.totalUnits = 1
config.Data.unitsPerJob = 1
config.Data.splitting = 'EventBased'
config.Data.ignoreLocality = True

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_US_UCSD'
config.Site.whitelist = ['T2_US_Caltech','T2_US_Florida', 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt', 'T2_US_Wisconsin']
