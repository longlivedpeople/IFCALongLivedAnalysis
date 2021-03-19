import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '106X_dataRun2_v32'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       [
'file:/eos/user/f/fernance/LLP_Analysis/miniAOD_extended/DoubleMuon_test/problematic/EXO-RunIISummer16MiniAODv3-08121_854.root'
       ]
    ),
    skipEvents = cms.untracked.uint32(0)
)

process.longlivedanalyzer.isData = True
process.longlivedanalyzer.DSAMode = True
process.longlivedanalyzer.Era = 2016

process.p = cms.Path(process.longlivedanalyzer)


