import FWCore.ParameterSet.Config as cms


longlivedanalyzer = cms.EDAnalyzer('LongLivedAnalysis',
    nameOfOutput = cms.string('output.root'),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ElectronCollection = cms.InputTag("slimmedElectrons"),
    genParticleCollection = cms.InputTag("prunedGenParticles"),
    PhotonCollection = cms.InputTag("slimmedPhotons"),
    IsoTrackCollection = cms.InputTag("isolatedTracks"),
    PrimaryVertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("slimmedPatTrigger"),
)


