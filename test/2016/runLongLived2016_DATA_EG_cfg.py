import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '106X_dataRun2_v32'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #['/store/data/Run2016E/DoubleEG/AOD/21Feb2020_UL2016_HIPM-v1/20000/02A8F48F-E2C7-9345-88AE-D918A9F7BE5C.root']
       ['/store/data/Run2016C/DoubleEG/MINIAOD/21Feb2020_UL2016_HIPM-v1/10000/099E723B-E1C1-1049-B0E8-F255C718A937.root']
    ),
    skipEvents = cms.untracked.uint32(0)
)


process.load("L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff")
process.es_prefer_l1GtTriggerMaskAlgoTrig = cms.ESPrefer("L1GtTriggerMaskAlgoTrigTrivialProducer","l1GtTriggerMaskAlgoTrig")

process.load('HLTrigger.HLTfilters.hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(False)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('NOT L1_DoubleEG_23_10')
process.hltLevel1GTSeed.L1GtObjectMapTag = cms.InputTag("gtStage2Digis")


process.longlivedanalyzer.isData = True
process.longlivedanalyzer.Era = 2016
process.longlivedanalyzer.DSAMode = False
process.longlivedanalyzer.FilterByElectron = True

process.p = cms.Path(process.hltLevel1GTSeed * process.longlivedanalyzer)
#process.p = cms.Path(process.longlivedanalyzer)


