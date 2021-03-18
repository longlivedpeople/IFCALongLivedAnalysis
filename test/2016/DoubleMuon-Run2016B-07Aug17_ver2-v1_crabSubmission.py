import CRABClient
from CRABClient.UserUtilities import config

config = config()

# All output/log files go in directory workArea/requestName/
#config.General.workArea = 'crab_projects'

config.General.requestName = 'DoubleMuon-Run2016B-07Aug17_ver2-v1_NTuples_rsb'

config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# Set pluginName = Analysis if you are reading a dataset, or to PrivateMC if not (so you are generating events)
config.JobType.pluginName = 'Analysis'
# CMSSW cfg file you wish to run
config.JobType.psetName = 'runLongLivedMuon2016Data_cfg.py'
# Increase virtual memory limit (sum needed by all threads) from default of 2000 MB.
config.JobType.maxMemoryMB = 2500
# Number of threads to use.
#config.JobType.numCores = 2
# To allow use of SL7 CMSSW versions that were SL6 at time of original DATA/MC production.
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['PUreweighting/2016DataPileupHistogram.root',
                             'PUreweighting/2016MCPileupHistogram.root']
config.JobType.outputFiles = ['output.root']


# Input dataset (small)
config.Data.inputDataset = '/DoubleMuon/fernance-DoubleMuon-Run2016B-07Aug17_ver2-v1_modified-42be9799ae70470560a96f9e21287209/USER'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'
config.Data.inputDBS = 'phys03'

# Units of "totalUnits" and "unitsPerJob" (e.g. files, events, lumi sections)
config.Data.splitting = 'FileBased'
# Total number of these units to be processed.
#config.Data.totalUnits = 999
# Requested number in each subjob.
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/fernance/'
config.Data.publication = False

# Output dataset stored in my dcache area in AAA/tomalin-outputDatasetTag-encodedDataAndTime/USER 
# where "AAA" is name between first pair of slashes in input dataset.
config.Data.outputDatasetTag = 'NTuples-Galapago'

config.Site.storageSite = 'T2_ES_IFCA'
