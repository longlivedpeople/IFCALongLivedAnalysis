import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '94X_mcRun2_asymptotic_v3'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root"
       #'file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root'
       # 'file:test/merged.root'
       [
          #'file:/eos/user/f/fernance/LLPNTuples/2016/test/HXX_400_50_400_test.root',
          #'file:/eos/user/f/fernance/LLP_Analysis/miniAOD_extended/DY_test/EXO-RunIISummer16MiniAODv3-08121_324.root',
          #'file:/eos/user/f/fernance/LLP_Analysis/miniAOD_extended/DY_test/EXO-RunIISummer16MiniAODv3-08121_325.root',
          #'/store/user/fernance/H2ToLLPXToLeptons_MH_400_MX_50_ctau_400mm_TuneCP2_13TeV_pythia8_80X_13082019-1313/400-50-400_HXX_RunIISummer16MiniAODv3_230220-1650/200226_084221/0000/EXO-RunIISummer16MiniAODv3-08121_1.root'
          'das://H2ToLLPXToLeptons_MH_400_MX_50_ctau_400mm_TuneCP2_13TeV_pythia8_80X_13082019-1313/fernance-400-50-400_HXX_RunIISummer16MiniAODv3_230220-1650-bd3e7bcff6c9bcad356ea4ed7e4f08b4/USER'
       ]
    )
)

process.longlivedanalyzer.isData = False
process.longlivedanalyzer.DSAMode = True
process.longlivedanalyzer.BSMode = False
process.longlivedanalyzer.Era = 2016

process.p = cms.Path(process.longlivedanalyzer)


