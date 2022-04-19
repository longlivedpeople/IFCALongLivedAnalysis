import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_preVFP_v9'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       [
          '/store/user/fernance/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/fernance-RunIISummer20UL16MiniAODAPVext-skim/220407_143620/0000/output_skim_1.root'
       ]
    )
)

process.longlivedanalyzer.isData = False
process.longlivedanalyzer.DSAMode = True
process.longlivedanalyzer.Era = 2016

process.p = cms.Path(process.longlivedanalyzer)


