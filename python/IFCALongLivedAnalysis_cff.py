import FWCore.ParameterSet.Config as cms


longlivedanalyzer = cms.EDAnalyzer('LongLivedAnalysis',
    nameOfOutput = cms.string('output.root'),
    MuonCollection = cms.InputTag("slimmedMuons"),
    PhotonCollection = cms.InputTag("slimmedPhotons"),
    IsoTrackCollection = cms.InputTag("isolatedTracks")
)


