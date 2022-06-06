import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '106X_mc2017_realistic_v9'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       [
          '/store/user/fernance/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-10_TuneCP5_13TeV-powheg-pythia8/private-RunIISummer20UL17MiniAODv2ext/210715_154417/0000/EXO-RunIISummer20UL17MiniAODv2_1.root'
       ]
    )
)

process.longlivedanalyzer.isData = False
process.longlivedanalyzer.DSAMode = False
process.longlivedanalyzer.FilterByLL = True
process.longlivedanalyzer.Era = 2017

process.p = cms.Path(process.longlivedanalyzer)


