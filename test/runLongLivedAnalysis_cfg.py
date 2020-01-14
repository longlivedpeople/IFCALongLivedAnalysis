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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root"
       #'file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root'
       # 'file:test/merged.root'
       [
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_100.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_101.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_102.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_103.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_104.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_105.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_106.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_107.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_108.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_109.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_10.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_110.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_111.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_112.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_113.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_114.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_115.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_116.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_117.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_118.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_119.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_11.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_120.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_121.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_122.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_123.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_124.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_125.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_126.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_127.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_128.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_129.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_12.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_130.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_131.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_132.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_133.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_134.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_135.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_136.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_137.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_138.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_139.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_13.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_140.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_141.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_142.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_143.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_144.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_145.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_146.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_147.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_148.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_149.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_14.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_150.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_151.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_152.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_153.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_154.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_155.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_156.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_157.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_158.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_159.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_15.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_160.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_161.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_162.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_163.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_164.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_165.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_166.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_167.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_168.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_169.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_16.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_170.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_171.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_172.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_173.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_174.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_175.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_176.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_177.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_178.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_179.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_17.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_180.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_181.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_182.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_183.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_184.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_185.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_186.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_187.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_188.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_189.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_18.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_190.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_191.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_192.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_193.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_194.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_195.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_196.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_197.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_198.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_199.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_19.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_1.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_200.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_201.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_202.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_203.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_204.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_205.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_206.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_207.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_208.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_209.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_20.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_210.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_211.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_212.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_213.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_214.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_215.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_216.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_217.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_218.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_219.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_21.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_220.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_221.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_222.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_223.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_224.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_225.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_226.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_227.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_228.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_229.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_22.root',
'file:/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_pythia8_80X_v1_MiniAODv3_230.root'
       ]
    )
)

process.p = cms.Path(process.longlivedanalyzer)


