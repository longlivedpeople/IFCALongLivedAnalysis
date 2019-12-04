import FWCore.ParameterSet.Config as cms


longlivedanalyzer = cms.EDAnalyzer('LongLivedAnalysis',
    isData = cms.bool(False),
    BSMode = cms.bool(True),
    nameOfOutput = cms.string('output.root'),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ElectronCollection = cms.InputTag("slimmedElectrons"),
    MuonCollection = cms.InputTag("slimmedMuons"),
    PackedPFCandidateCollection = cms.InputTag("packedPFCandidates"),
    LostTracksCollection = cms.InputTag("lostTracks"),
    genParticleCollection = cms.InputTag("prunedGenParticles"),
    PhotonCollection = cms.InputTag("slimmedPhotons"),
    IsoTrackCollection = cms.InputTag("isolatedTracks"),
    METCollection = cms.InputTag("slimmedMETs"),
    PrimaryVertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    #bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    objects = cms.InputTag("slimmedPatTrigger"),
    theGenEventInfoProduct = cms.InputTag("generator"),
    thePileUpSummary = cms.InputTag("slimmedAddPileupInfo")
)


