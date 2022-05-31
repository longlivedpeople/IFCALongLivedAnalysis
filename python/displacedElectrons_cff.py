import FWCore.ParameterSet.Config as cms


# Cloned from isolatedTracks setup
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock

tkAssocParamBlock = TrackAssociatorParameterBlock.clone()
tkAssocParamBlock.TrackAssociatorParameters.useMuon = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.useCalo = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.useHO = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.usePreshower = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("reducedEcalRecHitsEE")
tkAssocParamBlock.TrackAssociatorParameters.EBRecHitCollectionLabel = cms.InputTag("reducedEcalRecHitsEB")
tkAssocParamBlock.TrackAssociatorParameters.HBHERecHitCollectionLabel = cms.InputTag("reducedHcalRecHits","hbhereco")
tkAssocParamBlock.TrackAssociatorParameters.HORecHitCollectionLabel = cms.InputTag("reducedHcalRecHits","horeco")

displacedElectrons = cms.EDAnalyzer('displacedElectronAnalyzer',

    tkAssocParamBlock,
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
    PackedPFCandidateCollection   = cms.InputTag("packedPFCandidates"),
    LostTracksCollection          = cms.InputTag("lostTracks"),
    EleLostTracksCollection       = cms.InputTag("lostTracks", "eleTracks", "PAT"),
    genParticleCollection         = cms.InputTag("prunedGenParticles"),
    PhotonCollection              = cms.InputTag("slimmedPhotons"),
    IsoTrackCollection            = cms.InputTag("isolatedTracks"),
    PrimaryVertexCollection       = cms.InputTag("offlineSlimmedPrimaryVertices"),


)


