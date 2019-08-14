from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.requestName = # Fill

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runLongLivedAnalysis_cfg.py'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.inputFiles = ['PUreweighting/2016/2016DataPileupHistogram.root',
                             'PUreweighting/2016/2016MCPileupHistogram.root']
config.JobType.outputFiles = ['output.root']

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.inputDataset = # Fill
config.Data.publication = False
config.Data.outLFNDirBase = # Fill 

config.section_('Site')
config.Site.storageSite = 'T2_ES_IFCA'
