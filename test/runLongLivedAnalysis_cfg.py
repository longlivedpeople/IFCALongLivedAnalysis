import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '94X_mc2017_realistic_v12'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root"
        #'file:/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/Samples/DisplacedSUSY/miniAOD/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP5_14TeV_pythia8_miniAOD_N10000.root'
        #'file:test/DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_994.root'
        'dbs:/DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1/fernance-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_v1_MiniAODv3-bd3e7bcff6c9bcad356ea4ed7e4f08b4/USER'
        #'file:/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_9/src/2016/MiniAODv3/HIG-RunIISummer15GS-01425_4_MiniAODv3.root'
        #'file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root'
    )
)

process.p = cms.Path(process.longlivedanalyzer)


