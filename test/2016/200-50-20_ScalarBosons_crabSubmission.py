from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.requestName = '200-50-20_ScalarBosons_NTuple'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runLongLivedAnalysis_cfg.py'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['output.root']
config.section_('Data')
config.JobType.inputFiles = ['PUreweighting/2016/2016DataPileupHistogram.root',
                             'PUreweighting/2016/2016MCPileupHistogram.root']


config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.inputDataset = '/H2ToLLPXToLeptons_MH_200_MX_50_ctau_20mm_TuneCP2_13TeV_pythia8_80X_26062019-1529/fernance-ScalarBosons_RunIISummer16MiniAODv3_26062019-1529-bd3e7bcff6c9bcad356ea4ed7e4f08b4/USER'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/fernance/' 

config.section_('Site')
config.Site.storageSite = 'T2_ES_IFCA'
