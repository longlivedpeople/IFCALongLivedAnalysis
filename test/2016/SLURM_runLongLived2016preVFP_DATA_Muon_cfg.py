import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '106X_dataRun2_v32'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       [
       INPUT
       ]
    )
)

# Add lumilist
process.source.lumisToProcess = LumiList.LumiList(filename = '/gpfs/users/fernanc/UL-LL-analysis/CMSSW_10_6_20/src/MyAnalysis/IFCALongLivedAnalysis/test/2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt').getVLuminosityBlockRange()

process.longlivedanalyzer.nameOfOutput = OUTPUT
process.longlivedanalyzer.isData = True
process.longlivedanalyzer.DSAMode = True
process.longlivedanalyzer.Era = 2016

process.p = cms.Path(process.longlivedanalyzer)


