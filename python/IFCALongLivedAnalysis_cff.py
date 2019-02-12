import FWCore.ParameterSet.Config as cms


longlivedanalyzer = cms.EDAnalyzer('IFCALongLivedAnalysis',
    nameOfFile = cms.string('output.root'),
    MuonCollection = cms.InputTag("slimmedMuons")
)


