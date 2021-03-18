import FWCore.ParameterSet.Config as cms


longlivedanalyzer = cms.EDAnalyzer('LongLivedAnalysis',
    Era = cms.double(2018),
    isData = cms.bool(False),
    BSMode = cms.bool(False),
    DSAMode = cms.bool(False),
    nameOfOutput = cms.string('output.root'),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ElectronCollection = cms.InputTag("slimmedElectrons"),
    MuonCollection = cms.InputTag("slimmedMuons"),
    StandAloneCollection = cms.InputTag("standAloneMuons", "", "RECO"),
    DisplacedStandAloneCollection = cms.InputTag("displacedStandAloneMuons"),
    RefittedStandAloneCollection = cms.InputTag("refittedStandAloneMuons"),
    GlobalMuonCollection = cms.InputTag("globalMuons"),
    DisplacedGlobalMuonCollection = cms.InputTag("displacedGlobalMuons"),
    PackedPFCandidateCollection = cms.InputTag("packedPFCandidates"),
    LostTracksCollection = cms.InputTag("lostTracks"),
    EleLostTracksCollection = cms.InputTag("lostTracks", "eleTracks", "PAT"),
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


