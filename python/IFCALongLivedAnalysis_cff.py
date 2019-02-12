import FWCore.ParameterSet.Config as cms


longlivedanalyzer = cms.EDAnalyzer('LongLivedAnalysis',
    nameOfFile = cms.string('output.root'),
    MuonCollection = cms.InputTag("slimmedMuons")
)


