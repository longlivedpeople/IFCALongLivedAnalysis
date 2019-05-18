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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root"
       #'file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root'
       # 'file:test/merged.root'
       [
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_1.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_10.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_11.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_12.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_13.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_14.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_15.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_16.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_17.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_18.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_19.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_2.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_20.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_21.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_22.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_23.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_24.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_25.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_26.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_27.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_28.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_29.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_3.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_30.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_31.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_32.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_33.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_34.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_35.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_36.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_37.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_38.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_39.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_4.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_40.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_41.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_42.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_43.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_44.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_45.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_46.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_47.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_48.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_49.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_5.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_50.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_51.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_52.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_53.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_54.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_55.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_56.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_57.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_58.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_59.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_6.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_60.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_61.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_62.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_63.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_64.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_65.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_66.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_67.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_68.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_69.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_7.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_70.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_71.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_72.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_73.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_74.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_75.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_76.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_77.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_78.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_79.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_8.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_80.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_81.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_82.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_83.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_84.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_85.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_86.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_87.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_88.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_89.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_9.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_90.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_91.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_92.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_93.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_94.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_95.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_96.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_97.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_98.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/mini/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_99.root'
       ]
    )
)

process.p = cms.Path(process.longlivedanalyzer)


