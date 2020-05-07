from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.requestName = 'TTJets_DiLept_NTuple'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runLongLivedAnalysis_cfg.py'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.inputFiles = ['PUreweighting/2016DataPileupHistogram.root',
                             'PUreweighting/2016MCPileupHistogram.root']
config.JobType.outputFiles = ['output.root']
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 1
#config.Data.totalUnits = 40
config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/fernance-TTJets_RunIISummer16MiniAODv3_modified-bd3e7bcff6c9bcad356ea4ed7e4f08b4/USER'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/fernance/' 

config.section_('Site')
config.Site.storageSite = 'T2_ES_IFCA'
