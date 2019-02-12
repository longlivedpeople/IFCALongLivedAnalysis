import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/user/pablom/CRAB_PrivateMC/stop-stop-205-22p5_171p5/180429_085641/0000/SUS-RunIISummer15GS-00210_459.root'
    )
)

process.longlivedanalyzer = cms.EDAnalyzer('IFCALongLivedAnalysis',
    nameOfFile = cms.string('output.root')
)


process.p = cms.Path(process.longlivedanalyzer)
