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
       '/store/user/fernance/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVext/210330_095618/0000/EXO-RunIISummer20UL16MiniAODAPV_10.root'
       #'/store/mc/RunIISummer20UL16MiniAODAPV/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v8-v1/50000/001FEF0E-58A2-894A-8280-2F4AA900A5A6.root'    
       ]
    )
)

process.longlivedanalyzer.isData = False
process.longlivedanalyzer.DSAMode = True
process.longlivedanalyzer.FilterByLL = True
process.longlivedanalyzer.Era = 2016

process.p = cms.Path(process.longlivedanalyzer)


