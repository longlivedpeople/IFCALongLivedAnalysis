import FWCore.ParameterSet.Config as cms


longlivedanalyzer = cms.EDAnalyzer('LongLivedAnalysis',


    # Main details
    Era = cms.double(2018),
    nameOfOutput = cms.string('output.root'),


    # Running modes

    isData         = cms.bool(False),
    DSAMode        = cms.bool(False),
    doCMSElectrons = cms.bool(False), 
    doCMSMuons     = cms.bool(False), 


    # Collections

    EventInfo                     = cms.InputTag("generator"),
    RunInfo                       = cms.InputTag("generator"),
    BeamSpot                      = cms.InputTag("offlineBeamSpot"),
    theGenEventInfoProduct        = cms.InputTag("generator"),
    thePileUpSummary              = cms.InputTag("slimmedAddPileupInfo"),

    ElectronCollection            = cms.InputTag("slimmedElectrons"),
    MuonCollection                = cms.InputTag("slimmedMuons"),
    DisplacedGlobalMuonCollection = cms.InputTag("displacedGlobalMuons"),
    PackedPFCandidateCollection   = cms.InputTag("packedPFCandidates"),
    LostTracksCollection          = cms.InputTag("lostTracks"),
    EleLostTracksCollection       = cms.InputTag("lostTracks", "eleTracks", "PAT"),
    genParticleCollection         = cms.InputTag("prunedGenParticles"),
    PhotonCollection              = cms.InputTag("slimmedPhotons"),
    IsoTrackCollection            = cms.InputTag("isolatedTracks"),
    METCollection                 = cms.InputTag("slimmedMETs"),
    PrimaryVertexCollection       = cms.InputTag("offlineSlimmedPrimaryVertices"),

    # Trigger

    prescales  = cms.InputTag("patTrigger"),
    bits       = cms.InputTag("TriggerResults","","HLT"),
    objects    = cms.InputTag("slimmedPatTrigger")


)


