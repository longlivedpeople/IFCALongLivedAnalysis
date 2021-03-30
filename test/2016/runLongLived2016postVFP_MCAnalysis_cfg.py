import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_v15'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       [
          'das://H2ToLLPXToLeptons_MH_400_MX_50_ctau_400mm_TuneCP2_13TeV_pythia8_80X_13082019-1313/fernance-400-50-400_HXX_RunIISummer16MiniAODv3_230220-1650-bd3e7bcff6c9bcad356ea4ed7e4f08b4/USER'
       ]
    )
)

process.longlivedanalyzer.isData = False
process.longlivedanalyzer.DSAMode = True
process.longlivedanalyzer.Era = 2016

process.p = cms.Path(process.longlivedanalyzer)


