from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.requestName = 'DoubleMuon-Run2016B-07Aug17_ver2-v1_NTuple_rsb3'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runLongLivedMuonData_cfg.py'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.inputFiles = ['PUreweighting/2016DataPileupHistogram.root',
                             'PUreweighting/2016MCPileupHistogram.root']
config.JobType.outputFiles = ['output.root']
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/DoubleMuon/fernance-DoubleMuon-Run2016B-07Aug17_ver2-v1_modified-42be9799ae70470560a96f9e21287209/USER'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/fernance/' 

config.section_('Site')
config.Site.storageSite = 'T2_ES_IFCA'
