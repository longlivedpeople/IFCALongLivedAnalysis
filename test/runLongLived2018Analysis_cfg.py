import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v21'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       [
          #'/store/mc/RunIIAutumn18MiniAOD/ggH_HToSSTo4l_MH-110_MS-10_ctauS-1_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/280000/23989087-91AE-344D-A1E8-64171B95AF1C.root',
          '/store/user/fernance/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAODext/201114_192331/0000/RunIIAutumn18MiniAOD_1.root',
#          'file:/eos/user/f/fernance/LLP_Analysis/testfiles/2018/signalMINIAODSIM/23989087-91AE-344D-A1E8-64171B95AF1C.root',
       ]
    )
)

process.longlivedanalyzer.isData = False
process.longlivedanalyzer.DSAMode = True
process.longlivedanalyzer.BSMode = False
process.longlivedanalyzer.Era = 2018

process.p = cms.Path(process.longlivedanalyzer)


