//=======================================================================================================================================================================================================================// 
//                                                                                                                                                                                                                       // 
//$$$$$$\ $$$$$$$$\  $$$$$$\   $$$$$$\                      $$\                                    $$\       $$\                            $$\  $$$$$$\                      $$\                     $$                 //
//\_$$  _|$$  _____|$$  __$$\ $$  __$$\                     $$ |                                   $$ |      \__|                           $$ |$$  __$$\                     $$ |                    \__|               //
//  $$ |  $$ |      $$ /  \__|$$ /  $$ |                    $$ |      $$$$$$\  $$$$$$$\   $$$$$$\  $$ |      $$\ $$\    $$\  $$$$$$\   $$$$$$$ |$$ /  $$ |$$$$$$$\   $$$$$$\  $$ |$$\   $$\  $$$$$$$\ $$\  $$$$$$$       //
//  $$ |  $$$$$\    $$ |      $$$$$$$$ |      $$$$$$\       $$ |     $$  __$$\ $$  __$$\ $$  __$$\ $$ |      $$ |\$$\  $$  |$$  __$$\ $$  __$$ |$$$$$$$$ |$$  __$$\  \____$$\ $$ |$$ |  $$ |$$  _____|$$ |$$  _____|     //
//  $$ |  $$  __|   $$ |      $$  __$$ |      \______|      $$ |     $$ /  $$ |$$ |  $$ |$$ /  $$ |$$ |      $$ | \$$\$$  / $$$$$$$$ |$$ /  $$ |$$  __$$ |$$ |  $$ | $$$$$$$ |$$ |$$ |  $$ |\$$$$$$\  $$ |\$$$$$$        // 
//  $$ |  $$ |      $$ |  $$\ $$ |  $$ |                    $$ |     $$ |  $$ |$$ |  $$ |$$ |  $$ |$$ |      $$ |  \$$$  /  $$   ____|$$ |  $$ |$$ |  $$ |$$ |  $$ |$$  __$$ |$$ |$$ |  $$ | \____$$\ $$ | \____$$       //
//$$$$$$\ $$ |      \$$$$$$  |$$ |  $$ |                    $$$$$$$$\\$$$$$$  |$$ |  $$ |\$$$$$$$ |$$$$$$$$\ $$ |   \$  /   \$$$$$$$\ \$$$$$$$ |$$ |  $$ |$$ |  $$ |\$$$$$$$ |$$ |\$$$$$$$ |$$$$$$$  |$$ |$$$$$$$  |     //
//\______|\__|       \______/ \__|  \__|                    \________|\______/ \__|  \__| \____$$ |\________|\__|    \_/     \_______| \_______|\__|  \__|\__|  \__| \_______|\__| \____$$ |\_______/ \__|\_______/      //
//                                                                                       $$\   $$ |                                                                               $$\   $$ |                             // 
//                                                                                       \$$$$$$  |                                                                               \$$$$$$  |                             //
//=======================================================================================================================================================================================================================//
//                                                                                                                                                                                                                       //
// Authors of the code: Celia Fernandez Madrazo                                                                                                                                                                          //
//                      Pablo Martinez Ruiz Del Arbol                                                                                                                                                                    //
//                      Jesus Vizan Garcia                                                                                                                                                                               //
//=======================================================================================================================================================================================================================//
//                                                                                                                                                                                                                       //
// Description: Main analyzer                                                                                                                                                                                            //
//                                                                                                                                                                                                                       //
//=======================================================================================================================================================================================================================//



#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrackAssociator/interface/DetIdInfo.h"
#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"

#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Candidate/interface/Candidate.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "MyAnalysis/IFCALongLivedAnalysis/interface/llCandidateDataFormat.h"
#include "MyAnalysis/IFCALongLivedAnalysis/interface/trackPairDataFormat.h"


//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"



//#include "DataFormats/MuonReco/interface/MuonFwd.h" 

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

//=======================================================================================================================================================================================================================//

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// DATA DEFINITION //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////// BRANCHES /////////////////////////////////////



/////////////////////////////////////// OUTPUT //////////////////////////////////////



//=======================================================================================================================================================================================================================//
//

class displacedElectronAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit displacedElectronAnalyzer(const edm::ParameterSet&);
      ~displacedElectronAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // Parameters
      std::string output_filename;
      bool   _isData;
      bool   _DSAMode;
      bool   _doCMSElectrons;
      bool   _doCMSMuons;

      double _Era;

      edm::ParameterSet parameters;

      // Tokens
      edm::EDGetTokenT<edm::View<pat::Electron> > theElectronCollection;   
      edm::EDGetTokenT<edm::View<pat::Photon> > thePhotonCollection;
      edm::EDGetTokenT<edm::View<pat::IsolatedTrack> >  theIsoTrackCollection;
      edm::EDGetTokenT<edm::View<reco::Vertex> > thePrimaryVertexCollection;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > thePackedPFCandidateCollection;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > theLostTracksCollection;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > theEleLostTracksCollection;
      edm::EDGetTokenT<edm::View<pat::MET> > theMETCollection;

      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone> > triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales>  triggerPrescales_;

      edm::EDGetTokenT<reco::BeamSpot> theBeamSpot;

      // Gen collection
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  theGenParticleCollection;      
      edm::EDGetTokenT<GenEventInfoProduct>  theGenEventInfoProduct;

      // PU reweighting
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  thePileUpSummary;
      edm::LumiReWeighting lumi_weights;

      // "Global" variables
      std::vector<int> iT; // track indexes
      edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;

      // Magnetic field
      //const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
      edm::ESHandle<MagneticField> bField;

      // Tools to extrapolate the tracks
      TrackDetectorAssociator trackAssociator_;
      TrackAssociatorParameters trackAssocParameters_;

      // Class functions
      std::string getPathVersion(const edm:: TriggerNames &names, const std::string &rawPath);
      bool passIsotrackSelection(const pat::IsolatedTrack &track);
      float getDeltaR(float phi1, float eta1, float phi2, float eta2);
      float dxy_value(const pat::IsolatedTrack & track, const reco::Vertex pv);
      float computeDxy(const pat::IsolatedTrack & track, const reco::Vertex pv);
      float computeDxy(const reco::Track & track, const reco::Vertex pv);
      float computeDxy(const reco::GenParticle &p, const reco::Vertex &pv);
      float computeDxyError(const pat::IsolatedTrack & track, const reco::Vertex pv);
      float computeDxyError(const reco::Track & track, const reco::Vertex pv);
      float computeRelIso(const reco::Track & track,  edm::Handle<edm::View<pat::PackedCandidate> > pfs, edm::Handle<edm::View<pat::PackedCandidate> > lts, bool isPF, bool includeLostTracks);
      TVector3 extrapolateToECAL(reco::GenParticle gen);

      //-> EVENT INFO
      Int_t Event_event = 0;
      Int_t Event_luminosityBlock = 0;
      Int_t Event_run = 0;
      Int_t nPU = 0;
      Int_t nPUTrue = 0;
      Float_t wPU = 0.;
      Float_t genWeight = 0.;

      //-> PRIMARY VERTEX SELECTION
      Int_t nPV = 0;
      Int_t nTruePV = 0;
      Int_t PV_passAcceptance = 0;
      Float_t PV_vx = 0.;
      Float_t PV_vy = 0.;
      Float_t PV_vz = 0.;

      //-> BEAM SPOT
      Float_t BeamSpot_x0 = 0.;
      Float_t BeamSpot_y0 = 0.;
      Float_t BeamSpot_z0 = 0.;
      Float_t BeamSpot_BeamWidthX = 0.;
      Float_t BeamSpot_BeamWidthY = 0.;


      //-> ISOTRACK SELECTION
      Int_t nIsoTrack = 0;
      // Primitive:
      Float_t IsoTrackSel_pt[500] = {0.};
      Float_t IsoTrackSel_eta[500] = {0.};
      Float_t IsoTrackSel_etaExtra[500] = {0.};
      Float_t IsoTrackSel_phiExtra[500] = {0.};
      Float_t IsoTrackSel_phi[500] = {0.};
      Int_t IsoTrackSel_charge[500] = {0};
      Float_t IsoTrackSel_dxy[500] = {0.};
      Float_t IsoTrackSel_dxyError[500] = {0.};
      Float_t IsoTrackSel_dxy_PV[500] = {0.};
      Float_t IsoTrackSel_dxyError_PV[500] = {0.};
      Float_t IsoTrackSel_dxy_BS[500] = {0.};
      Float_t IsoTrackSel_dxyError_BS[500] = {0.};
      Float_t IsoTrackSel_dz[500] = {0.};
      Float_t IsoTrackSel_dzError[500] = {0.};
      Float_t IsoTrackSel_vx[500] = {0.};
      Float_t IsoTrackSel_vy[500] = {0.};
      Float_t IsoTrackSel_vz[500] = {0.};
      Float_t IsoTrackSel_pfIsolationDR03[500] = {0.};
      Float_t IsoTrackSel_miniPFIsolation[500] = {0.};
      Float_t IsoTrackSel_relPfIsolationDR03[500] = {0.};
      Float_t IsoTrackSel_relMiniPFIsolation[500] = {0.};
      Int_t IsoTrackSel_isHighPurityTrack[500] = {0};
      Int_t IsoTrackSel_numberOfValidTrackerHits[500] = {0};
      Int_t IsoTrackSel_numberOfValidPixelHits[500] = {0};
      Int_t IsoTrackSel_numberOfValidPixelBarrelHits[500] = {0};
      Int_t IsoTrackSel_numberOfValidPixelEndcapHits[500] = {0};
      Int_t IsoTrackSel_numberOfValidStripHits[500] = {0};
      Int_t IsoTrackSel_numberOfValidStripTIBHits[500] = {0};
      Int_t IsoTrackSel_numberOfValidStripTIDHits[500] = {0};
      Int_t IsoTrackSel_numberOfValidStripTOBHits[500] = {0};
      Int_t IsoTrackSel_numberOfValidStripTECHits[500] = {0};
      Int_t IsoTrackSel_fromPV[500] = {0};
      Float_t IsoTrackSel_PVx[500] = {0.};
      Float_t IsoTrackSel_PVy[500] = {0.};
      Float_t IsoTrackSel_PVz[500] = {0.};

      //-> PHOTON SELECTION
      Int_t nPhoton = 0;
      Float_t PhotonSel_et[100] = {0.};
      Float_t PhotonSel_eta[100] = {0.};
      Float_t PhotonSel_phi[100] = {0.};
      Float_t PhotonSel_hadronicOverEm[100] = {0.};
      Float_t PhotonSel_full5x5_sigmaIetaIeta[100] = {0.};
      Int_t PhotonSel_isEB[100] = {0};
      Int_t PhotonSel_isEE[100] = {0};
      Float_t PhotonSel_r9[100] = {0.};

      //-> ELECTRON SELECTION
      Int_t nElectron = 0;
      Float_t ElectronSel_pt[100] = {0.};
      Float_t ElectronSel_et[100] = {0.};
      Float_t ElectronSel_eta[100] = {0.};
      Float_t ElectronSel_phi[100] = {0.};
      Float_t ElectronSel_dxy[100] = {0.};
      Float_t ElectronSel_dxyError[100] = {0.};
      Float_t ElectronSel_dxySignificance[100] = {0.};
      Float_t ElectronSel_dB[100] = {0.};
      Float_t ElectronSel_edB[100] = {0.};
      Float_t ElectronSel_isLoose[100] = {0.};
      Float_t ElectronSel_isMedium[100] = {0.};
      Float_t ElectronSel_isTight[100] = {0.};


      //-> GENLEPTON SELECTION
      Int_t nGenLepton = 0;
      Int_t nGenLepton_PFS = 0;
      Int_t nGenLepton_HPFS = 0;
      Int_t nGenLepton_PTDP = 0;
      Int_t nGenLepton_HDP = 0;
      Float_t GenLeptonSel_pt[50] = {0.};
      Float_t GenLeptonSel_E[50] = {0.};
      Float_t GenLeptonSel_et[50] = {0.};
      Float_t GenLeptonSel_eta[50] = {0.};
      Float_t GenLeptonSel_phi[50] = {0.};
      Float_t GenLeptonSel_etaAtECAL[100] = {0.};
      Float_t GenLeptonSel_phiAtECAL[100] = {0.};
      Int_t GenLeptonSel_pdgId[50] = {0};
      Float_t GenLeptonSel_dxy[50] = {0.};
      Float_t GenLeptonSel_vx[50] = {0.};
      Float_t GenLeptonSel_vy[50] = {0.};
      Float_t GenLeptonSel_vz[50] = {0.};
      Int_t GenLeptonSel_motherPdgId[50] = {0};
      Int_t GenLeptonSel_fromHardProcessFinalState[50] = {0};
      Int_t GenLeptonSel_isPromptFinalState[50] = {0};
      Int_t GenLeptonSel_isDirectPromptTauDecayProductFinalState[50] = {0};
      Int_t GenLeptonSel_isDirectHadronDecayProduct[50] = {0};
      Int_t GenLeptonSel_hasRadiated[50] = {0};

      Int_t nHardProcessParticle = 0;
      Float_t HardProcessParticle_pt[30] = {0.};
      Float_t HardProcessParticle_E[30] = {0.};
      Float_t HardProcessParticle_eta[30] = {0.};
      Float_t HardProcessParticle_phi[30] = {0.};
      Float_t HardProcessParticle_vx[30] = {0.};
      Float_t HardProcessParticle_vy[30] = {0.};
      Float_t HardProcessParticle_vz[30] = {0.};
      Int_t HardProcessParticle_pdgId[30] = {0};


      //-> ELECTRON CANDIDATE SELECTION
      Int_t nElectronCandidate = 0;
      Float_t ElectronCandidate_pt[50] = {0.};
      Float_t ElectronCandidate_et[50] = {0.};
      Float_t ElectronCandidate_eta[50] = {0.};
      Float_t ElectronCandidate_phi[50] = {0.};
      Float_t ElectronCandidate_dxy[50] = {0.};
      Float_t ElectronCandidate_dxyError[50] = {0.};
      Float_t ElectronCandidate_dxy_PV[50] = {0.};
      Float_t ElectronCandidate_dxyError_PV[50] = {0.};
      Float_t ElectronCandidate_dxy_BS[50] = {0.};
      Float_t ElectronCandidate_dxyError_BS[50] = {0.};
      Float_t ElectronCandidate_relPFiso[50] = {0.};
      Float_t ElectronCandidate_relTrkiso[50] = {0.};
      Float_t ElectronCandidate_relTrkiso_noLT[50] = {0.};
      Int_t ElectronCandidate_photonIdx[50] = {0};
      Int_t ElectronCandidate_isotrackIdx[50] = {0};
      Int_t ElectronCandidate_pvAssociationQuality[50] = {0};
      Float_t ElectronCandidate_ptDiff[50] = {0.};
      Float_t ElectronCandidate_dR[50] = {0.};


      // -> All EE candidates
      Int_t nEE = 0;
      Int_t EE_idxA[20] = {0};
      Int_t EE_idxB[20] = {0};
      Float_t EE_Lxy_PV[20] = {0.};
      Float_t EE_Ixy_PV[20] = {0.};
      Float_t EE_Lxy_0[20] = {0.};
      Float_t EE_Ixy_0[20] = {0.};
      Float_t EE_Lxy_BS[20] = {0.};
      Float_t EE_Ixy_BS[20] = {0.};
      Float_t EE_trackDxy[20] = {0.};
      Float_t EE_trackIxy[20] = {0.};
      Float_t EE_trackDxy_PV[20] = {0.};
      Float_t EE_trackIxy_PV[20] = {0.};
      Float_t EE_trackDxy_BS[20] = {0.};
      Float_t EE_trackIxy_BS[20] = {0.};
      Float_t EE_vx[20] = {0.};
      Float_t EE_vy[20] = {0.};
      Float_t EE_mass[20] = {0.};
      Float_t EE_normalizedChi2[20] = {0.};
      Float_t EE_leadingPt[20] = {0.};
      Float_t EE_subleadingPt[20] = {0.};
      Float_t EE_leadingEt[20] = {0.};
      Float_t EE_subleadingEt[20] = {0.};
      Float_t EE_cosAlpha[20] = {0.};
      Float_t EE_dR[20] = {0.};
      Float_t EE_dPhi[20] = {0.};
      Float_t EE_lldPhi[20] = {0.};
      Float_t EE_relisoA[20] = {0.};
      Float_t EE_relisoB[20] = {0.};





      // Output
      TH1F *counts, *sum2Weights;
      TFile *file_out;
      TTree *tree_out;
};
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
displacedElectronAnalyzer::displacedElectronAnalyzer(const edm::ParameterSet& iConfig)
{
   usesResource("TFileService");
  
   parameters = iConfig;

   counts = new TH1F("counts", "", 1, 0, 1);
   sum2Weights = new TH1F("sum2Weights", "", 1, 0, 1);

   theElectronCollection = consumes<edm::View<pat::Electron> >  (parameters.getParameter<edm::InputTag>("ElectronCollection"));
   thePhotonCollection = consumes<edm::View<pat::Photon> > (parameters.getParameter<edm::InputTag>("PhotonCollection"));
   theIsoTrackCollection = consumes<edm::View<pat::IsolatedTrack> >  (parameters.getParameter<edm::InputTag>("IsoTrackCollection"));
   thePrimaryVertexCollection = consumes<edm::View<reco::Vertex> >  (parameters.getParameter<edm::InputTag>("PrimaryVertexCollection"));
   thePackedPFCandidateCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("PackedPFCandidateCollection"));
   theLostTracksCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("LostTracksCollection"));
   theEleLostTracksCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("EleLostTracksCollection"));

   theBeamSpot = consumes<reco::BeamSpot>  (parameters.getParameter<edm::InputTag>("BeamSpot"));

   theGenParticleCollection = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genParticleCollection"));

   theGenEventInfoProduct = consumes<GenEventInfoProduct> (parameters.getParameter<edm::InputTag>("theGenEventInfoProduct"));

   thePileUpSummary = consumes<std::vector<PileupSummaryInfo> > (parameters.getParameter<edm::InputTag>("thePileUpSummary"));

   edm::ParameterSet assocParameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
   edm::ConsumesCollector iC = consumesCollector();
   trackAssocParameters_.loadParameters(assocParameters, iC);
   trackAssociator_.useDefaultPropagator();

}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
displacedElectronAnalyzer::~displacedElectronAnalyzer()
{

}
//=======================================================================================================================================================================================================================//



//=======================================================================================================================================================================================================================//
void displacedElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////// MAIN CODE /////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////



   //////////////////////////////// GET THE COLLECTIONS ////////////////////////////////
   
   // Handle:
   edm::Handle<edm::View<pat::Electron> > electrons;
   edm::Handle<edm::View<pat::Photon> > photons;
   edm::Handle<edm::View<pat::IsolatedTrack> > isotracks;
   edm::Handle<edm::View<reco::Vertex> > primaryvertices;
   edm::Handle<edm::View<pat::PackedCandidate> > packedPFCandidates;
   edm::Handle<edm::View<pat::PackedCandidate> > lostTracks;
   edm::Handle<edm::View<pat::PackedCandidate> > eleLostTracks;
   edm::Handle<reco::BeamSpot> beamSpot;
   edm::Handle<edm::View<reco::GenParticle> > genParticles;
   edm::Handle<GenEventInfoProduct> genEvtInfo;
   edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
   
   // Get tokens:
   iEvent.getByToken(theElectronCollection, electrons);
   iEvent.getByToken(thePhotonCollection, photons);
   iEvent.getByToken(theIsoTrackCollection, isotracks);
   iEvent.getByToken(thePrimaryVertexCollection, primaryvertices);
   iEvent.getByToken(thePackedPFCandidateCollection, packedPFCandidates);
   iEvent.getByToken(theLostTracksCollection, lostTracks);
   iEvent.getByToken(theEleLostTracksCollection, eleLostTracks);
   iEvent.getByToken(theBeamSpot, beamSpot);

   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);

   iSetup.get<IdealMagneticFieldRecord>().get(bField);

   if (!_isData){
      
       iEvent.getByToken(theGenParticleCollection, genParticles);    
       iEvent.getByToken(theGenEventInfoProduct, genEvtInfo);
       iEvent.getByToken(thePileUpSummary, puInfoH);
   }
   
   //auto const& bField = iSetup.getData(bFieldToken_);

   //// -----------------
   //// --
   //// ---- Event info 
   //// --
   //// -----------------

   // General event information:
   Event_event = iEvent.id().event();
   Event_run = iEvent.id().run();
   Event_luminosityBlock = iEvent.id().luminosityBlock();
   

   counts->Fill(0.5);

   // Monte Carlo (only) information:
   if (!_isData){
       genWeight = (float) genEvtInfo->weight();
       sum2Weights->Fill(0.5, genWeight*genWeight);
   }


   //// ------------------------
   //// --
   //// ---- PU info (only MC)
   //// --
   //// ------------------------

   nPUTrue = -1; // True number of primary vertices (default)
   wPU = 1; // Pile-up reweighting factor

   // Compute the lumi weight for each event:
   if (!_isData){

       for(size_t i=0;i<puInfoH->size();++i) {
           if( puInfoH->at(i).getBunchCrossing() == 0) {
                nPU = puInfoH->at(i).getPU_NumInteractions();
                nPUTrue = puInfoH->at(i).getTrueNumInteractions();
                continue;                                          
           }
       }
       wPU = lumi_weights.weight(nPUTrue); // PU weight

   }

   //// ------------------------
   //// --
   //// ---- Beam spot object 
   //// --
   //// ------------------------

   reco::BeamSpot beamSpotObject = *beamSpot;
   BeamSpot_x0 = beamSpotObject.x0();
   BeamSpot_y0 = beamSpotObject.y0();
   BeamSpot_z0 = beamSpotObject.z0();
   BeamSpot_BeamWidthX = beamSpotObject.BeamWidthX();
   BeamSpot_BeamWidthY = beamSpotObject.BeamWidthY();
   GlobalPoint _BSpoint(beamSpotObject.x0(), beamSpotObject.y0(), beamSpotObject.z0());
   GlobalPoint _0point(0.0, 0.0, 0.0);


   //// ----------------------------
   //// --
   //// ---- Slimmed PV Collection
   //// --
   //// ----------------------------


   // General PV's information:
   nTruePV = 0;
   nPV = primaryvertices->size();


   for (size_t i = 0; i < primaryvertices->size(); i ++){

       const reco::Vertex &current_vertex = (*primaryvertices)[i];
       if(current_vertex.isValid()){ nTruePV++; }

   }

   // The PV information:
   const reco::Vertex &thePrimaryVertex = (*primaryvertices)[0];

   PV_vx = thePrimaryVertex.x();
   PV_vy = thePrimaryVertex.y();
   PV_vz = thePrimaryVertex.z();
   PV_passAcceptance = false;

   if (!thePrimaryVertex.isFake() && thePrimaryVertex.ndof() > 4 && fabs(thePrimaryVertex.z()) < 25 && thePrimaryVertex.position().rho() <= 2) {
       PV_passAcceptance = true;
   }

   GlobalPoint _PVpoint(thePrimaryVertex.x(), thePrimaryVertex.y(), thePrimaryVertex.z());



   //// ---------------------------
   //// --
   //// ---- Isotracks Collection
   //// --
   //// ---------------------------

   // Isotrack init:
   iT.clear();
   
   for (size_t i = 0; i < isotracks->size(); i++){

       const pat::IsolatedTrack & isotrack = (*isotracks)[i];

       if (!passIsotrackSelection(isotrack)){ continue; } // some quality cuts

       iT.push_back(i);

   }

   nIsoTrack = iT.size(); // number of isotracks

   // Sort the isotracks by pt:
   std::sort( std::begin(iT), std::end(iT), [&](int i1, int i2){ return isotracks->at(i1).pt() > isotracks->at(i2).pt(); });


   // Loop over the isotracks:
   for (size_t i = 0; i < iT.size(); ++i){

       
       const pat::IsolatedTrack & isotrack = (*isotracks)[iT.at(i)];

       // Basic features:
       IsoTrackSel_pt[i] = isotrack.pt();
       IsoTrackSel_eta[i] = isotrack.eta();
       IsoTrackSel_phi[i] = isotrack.phi();
       IsoTrackSel_charge[i] = isotrack.charge();      

       // Track extrapolations to ECAL surface (used to do the electron matching)
       // To be noticed: deltaPhi and deltaEta should be added to the raw values (checked: 16/01/2020)
       IsoTrackSel_etaExtra[i] = isotrack.eta() + isotrack.deltaEta();
       IsoTrackSel_phiExtra[i] = isotrack.phi() + isotrack.deltaPhi();
 
       // Isolation info:
       const pat::PFIsolation &pfiso = isotrack.pfIsolationDR03();
       const pat::PFIsolation &minipfiso = isotrack.miniPFIsolation();

       double neutralIso = fmax(0.0, pfiso.photonIso() + pfiso.neutralHadronIso() - 0.5*pfiso.puChargedHadronIso());
       double chargedIso = pfiso.chargedHadronIso();
       IsoTrackSel_pfIsolationDR03[i] = neutralIso + chargedIso;
       IsoTrackSel_relPfIsolationDR03[i] = IsoTrackSel_pfIsolationDR03[i]/isotrack.pt();

       double miniNeutralIso = fmax(0.0, minipfiso.photonIso() + minipfiso.neutralHadronIso() - 0.5*minipfiso.puChargedHadronIso());
       double miniChargedIso = minipfiso.chargedHadronIso();
       IsoTrackSel_miniPFIsolation[i] = miniNeutralIso + miniChargedIso;
       IsoTrackSel_relMiniPFIsolation[i] = IsoTrackSel_miniPFIsolation[i]/isotrack.pt();


       // Quality info:
       IsoTrackSel_isHighPurityTrack[i] = isotrack.isHighPurityTrack();


       // Hit info:
       const reco::HitPattern &hits = isotrack.hitPattern();

       IsoTrackSel_numberOfValidTrackerHits[i] = hits.numberOfValidTrackerHits();
       IsoTrackSel_numberOfValidPixelHits[i] = hits.numberOfValidPixelHits();
       IsoTrackSel_numberOfValidPixelBarrelHits[i] = hits.numberOfValidPixelBarrelHits();
       IsoTrackSel_numberOfValidPixelEndcapHits[i] = hits.numberOfValidPixelEndcapHits();
       IsoTrackSel_numberOfValidStripHits[i] = hits.numberOfValidStripHits();
       IsoTrackSel_numberOfValidStripTIBHits[i] = hits.numberOfValidStripTIBHits();
       IsoTrackSel_numberOfValidStripTIDHits[i] = hits.numberOfValidStripTIDHits();
       IsoTrackSel_numberOfValidStripTOBHits[i] = hits.numberOfValidStripTOBHits();
       IsoTrackSel_numberOfValidStripTECHits[i] = hits.numberOfValidStripTECHits();

       // Info extracted form the packedCandidate of the isotrack
       //fromPV: 0 NoPV, 1 PVLoose, 2 PVTight, 3 used in the PV fit
       IsoTrackSel_fromPV[i] = isotrack.fromPV(); 
       
       const pat::PackedCandidateRef &pckCand = isotrack.packedCandRef(); // access the packed candidate
        
       //This catches all the isotracks (selected at this point)
       if (isotrack.fromPV() > -1){ // check it has a PV

           IsoTrackSel_vx[i] = (*pckCand).vx();
           IsoTrackSel_vy[i] = (*pckCand).vy();
           IsoTrackSel_vz[i] = (*pckCand).vz();

           const reco::VertexRef &PV = (*pckCand).vertexRef(); // access the PV of the candidate
           IsoTrackSel_PVx[i] = (*PV).x();
           IsoTrackSel_PVy[i] = (*PV).y();
           IsoTrackSel_PVz[i] = (*PV).z();
           
           // Access to the track and computation of impact parameters:
           //const reco::Track *trref = (*pckCand).pseudoTrack();
           //const reco::Track &ctr = *trref;
           const reco::Track ctr = (*pckCand).pseudoTrack();
           reco::TransientTrack _isotk = theTransientTrackBuilder->build(ctr);
           TrajectoryStateClosestToPoint _trajPV = _isotk.trajectoryStateClosestToPoint( _PVpoint );
           TrajectoryStateClosestToPoint _trajBS = _isotk.trajectoryStateClosestToPoint( _BSpoint );
                   
           // Impact parameter info:
           IsoTrackSel_dxy[i] = (*pckCand).dxy();
           IsoTrackSel_dxyError[i] = (*pckCand).dxyError();
           IsoTrackSel_dxy_PV[i] = (*pckCand).dxy(thePrimaryVertex.position());
           IsoTrackSel_dxyError_PV[i] = _trajPV.perigeeError().transverseImpactParameterError();
           IsoTrackSel_dxy_BS[i] = ctr.dxy(beamSpotObject);
           IsoTrackSel_dxyError_BS[i] = _trajBS.perigeeError().transverseImpactParameterError();
           IsoTrackSel_dz[i] = (*pckCand).dz(thePrimaryVertex.position());
           IsoTrackSel_dzError[i] = (*pckCand).dzError(); 



       }else{

           IsoTrackSel_vx[i] = -99;
           IsoTrackSel_vy[i] = -99;
           IsoTrackSel_vz[i] = -99;

           IsoTrackSel_PVx[i] = -99;
           IsoTrackSel_PVy[i] = -99;
           IsoTrackSel_PVz[i] = -99;

       }


   }

   std::cout << "Hey" << std::endl;


   //// ----------------------------
   //// --
   //// ---- CMS Photon Collection
   //// --
   //// ----------------------------
   
   std::vector<int> iP; // photon indexes


   // Select good photons
   for (size_t i = 0; i < photons->size(); i++){

       const pat::Photon & photon = (*photons)[i];
       iP.push_back(i);

   }

   // Sort good lepton indexes by pt
   std::sort( std::begin(iP), std::end(iP), [&](int i1, int i2){ return photons->at(i1).et() > photons->at(i2).et(); });


   // Loop over the photons
   nPhoton = iP.size();

   for (size_t i = 0; i < iP.size(); i++){

       const pat::Photon & photon = (*photons)[iP.at(i)];

       PhotonSel_et[i] = photon.et();
       PhotonSel_eta[i] = photon.eta();
       PhotonSel_phi[i] = photon.phi();
       PhotonSel_hadronicOverEm[i] = photon.hadronicOverEm();
       PhotonSel_full5x5_sigmaIetaIeta[i] = photon.full5x5_sigmaIetaIeta();
       PhotonSel_isEB[i] = photon.isEB();
       PhotonSel_isEE[i] = photon.isEE();
       PhotonSel_r9[i] = photon.r9();

   }



   //// ------------------------------
   //// --
   //// ---- CMS Electron Collection
   //// --
   //// ------------------------------
    

   if (true) {
  
      std::vector<int> iE; // electron indexes


      // Select good electrons
      for (size_t i = 0; i < electrons->size(); i++){

         //const pat::Electron & electron = (*electrons)[i];

         iE.push_back(i);

      }

      // Sort good lepton indexes by pt
      std::sort( std::begin(iE), std::end(iE), [&](int i1, int i2){ return electrons->at(i1).pt() > electrons->at(i2).pt(); });


      nElectron = iE.size();

       for (size_t i = 0; i < iE.size(); i++){

          const pat::Electron & electron = (* electrons)[iE.at(i)];

          ElectronSel_pt[i] = electron.pt();
          ElectronSel_et[i] = electron.et();
          ElectronSel_eta[i] = electron.eta();
          ElectronSel_phi[i] = electron.phi();
          ElectronSel_dxyError[i] = electron.dxyError();
          ElectronSel_dxy[i] = electron.gsfTrack()->dxy();
          ElectronSel_dxySignificance[i] = fabs(ElectronSel_dxy[i])/electron.dxyError();

          ElectronSel_dB[i] = electron.dB();
          ElectronSel_edB[i] = electron.edB();

          ElectronSel_isLoose[i] = electron.electronID("cutBasedElectronID-Summer16-80X-V1-loose");
          ElectronSel_isMedium[i] = electron.electronID("cutBasedElectronID-Summer16-80X-V1-medium");
          ElectronSel_isTight[i] = electron.electronID("cutBasedElectronID-Summer16-80X-V1-tight");

      }
   }
   

   //// ------------------------------
   //// --
   //// ---- GenParticle Collections
   //// --
   //// ------------------------------

   std::vector<int> iGL; // Generated lepton indexes

   reco::GenParticleRef mref;
   reco::GenParticle m;

   if (!_isData)
   {

       //// *** Get the leptons that can be idenitified in the detector ***
       // pdgId: 11, -11, 13 and -13
       for(size_t i = 0; i < genParticles->size(); i++) {

	   const reco::GenParticle &genparticle = (*genParticles)[i];

           if ( ( abs(genparticle.pdgId()) == 11 ) && genparticle.status() == 1) {
               iGL.push_back(i);
           }
       }

       nGenLepton = iGL.size();
       std::sort( std::begin(iGL), std::end(iGL), [&](int i1, int i2){ return genParticles->at(i1).pt() > genParticles->at(i2).pt(); });

       for(size_t i = 0; i < iGL.size(); i++) {

          const reco::GenParticle &genparticle = (*genParticles)[iGL.at(i)];
    
          GenLeptonSel_pt[i] = genparticle.pt();
          GenLeptonSel_E[i] = genparticle.energy();
          GenLeptonSel_et[i] = genparticle.et();
          GenLeptonSel_eta[i] = genparticle.eta();
          GenLeptonSel_phi[i] = genparticle.phi();
          GenLeptonSel_pdgId[i] = genparticle.pdgId();

          // Bottom-up to get the real decaying particle:
          if (genparticle.mother()->pdgId() == genparticle.pdgId()) {

              GenLeptonSel_hasRadiated[i] = 1;
              mref = genparticle.motherRef();
              m = *mref;
              while (m.pdgId() == m.mother()->pdgId()) {
                  mref = m.motherRef();
                  m = *mref;
              }

              GenLeptonSel_vx[i] = m.vx();
              GenLeptonSel_vy[i] = m.vy();
              GenLeptonSel_vz[i] = m.vz();
	      GenLeptonSel_dxy[i] = computeDxy(m, thePrimaryVertex); // should be computed here or before?

              if(m.numberOfMothers() != 0){
                  GenLeptonSel_motherPdgId[i] = m.motherRef()->pdgId();
              } else {
                  GenLeptonSel_motherPdgId[i] = 0; 
              }
          }else{

              GenLeptonSel_hasRadiated[i] = 0;
              GenLeptonSel_vx[i] = genparticle.vx();
              GenLeptonSel_vy[i] = genparticle.vy();
              GenLeptonSel_vz[i] = genparticle.vz();
	      GenLeptonSel_dxy[i] = computeDxy(genparticle, thePrimaryVertex); // should be computed here or before?

              GenLeptonSel_motherPdgId[i] = genparticle.motherRef()->pdgId();
          }

   
          // Get interaction with ECAL       
          TVector3 intECAL;
          intECAL = extrapolateToECAL(genparticle);
          if (intECAL.Px() < 0.0000001 && intECAL.Py() < 0.0000001 && intECAL.Pz() < 0.0000001 ) {
            GenLeptonSel_etaAtECAL[i] = -99;
            GenLeptonSel_phiAtECAL[i] = -99;
          } else {
            GenLeptonSel_etaAtECAL[i] = intECAL.Eta();
            GenLeptonSel_phiAtECAL[i] = intECAL.Phi();
          }

          // Flags:
          GenLeptonSel_isPromptFinalState[i] = genparticle.isPromptFinalState();
          GenLeptonSel_fromHardProcessFinalState[i] = genparticle.fromHardProcessFinalState(); // has to be done with the last one
          GenLeptonSel_isDirectPromptTauDecayProductFinalState[i] = genparticle.isDirectPromptTauDecayProductFinalState(); 
          GenLeptonSel_isDirectHadronDecayProduct[i] = genparticle.statusFlags().isDirectHadronDecayProduct(); 


       }


       // Counters initialization
       nGenLepton_PFS = 0; 
       nGenLepton_HPFS = 0; 
       nGenLepton_PTDP = 0; 
       nGenLepton_HDP = 0; 

       for(size_t i = 0; i < iGL.size(); i++) {

           if (GenLeptonSel_isPromptFinalState[i]) { nGenLepton_PFS++; }
           if (GenLeptonSel_fromHardProcessFinalState[i]) { nGenLepton_HPFS++; }
           if (GenLeptonSel_isDirectPromptTauDecayProductFinalState[i]) { nGenLepton_PTDP++; }
           if (GenLeptonSel_isDirectHadronDecayProduct[i]) { nGenLepton_HDP++; }

       }


       // -> Hard Process Collection
       nHardProcessParticle = 0;
       for(size_t i = 0; i < genParticles->size(); i++) {

          const reco::GenParticle &genparticle = (*genParticles)[i];
       
           if (genparticle.isHardProcess()){

               HardProcessParticle_pt[nHardProcessParticle] = genparticle.pt();
               HardProcessParticle_E[nHardProcessParticle] = genparticle.energy();
               HardProcessParticle_eta[nHardProcessParticle] = genparticle.eta();
               HardProcessParticle_phi[nHardProcessParticle] = genparticle.phi();
               HardProcessParticle_vx[nHardProcessParticle] = genparticle.vx();
               HardProcessParticle_vy[nHardProcessParticle] = genparticle.vy();
               HardProcessParticle_vz[nHardProcessParticle] = genparticle.vz();
               HardProcessParticle_pdgId[nHardProcessParticle] = genparticle.pdgId();

               nHardProcessParticle++;
           }       
       }



   }


   //// -----------------------------------
   //// --
   //// ---- Displaced Electrons reconstruction
   //// --
   //// -----------------------------------


   // Variable initiallization:

   std::vector<int> matched_tracks, matched_SC, matched_triggerObjects; // std vectors with matched objects to avoid overlapping

   float dRMin = 99999; // dR to minimize as high as possible in the beginning
   float dRThreshold = 0.3; // Maximum dR to do the lepton matching
   float dR; // Computation of dR
   int tmin, scmin, li; // minimum track, minimum SC, minimum trigger object, reconstructed lepton index


   // while loop that looks for the minimum matchings and stops when the dRThreshold is reached:

   while (1){

       dRMin = 99999; // redefinicion

       // Loop over the tracks
       for (size_t t = 0; t < iT.size(); t++){

           const pat::IsolatedTrack & isotrack = (*isotracks)[iT.at(t)];

           // pass if the track is associated already:
           if(std::find(matched_tracks.begin(), matched_tracks.end(), t) != matched_tracks.end()){ continue; }
           // pass if the track does not fulfil the prerequisites:
           if(!passIsotrackSelection(isotrack)){ continue; }
           // Reject tracks falling in Barrel-Endcap transition (only to do the matching)
           if (fabs(isotrack.eta()) > 1.4442 and fabs(isotrack.eta()) < 1.566) { continue; }


           // Loop over the superclusters:
           for (size_t sc = 0; sc < iP.size(); sc++){

               const pat::Photon & photon = (*photons)[iP.at(sc)];

               // pass if the SC is associated to other track:
               if(std::find(matched_SC.begin(), matched_SC.end(), sc) != matched_SC.end()){ continue; }

               // pass if the SC does not fulfil the prerequisites:
               //if (!goodPhoton(photon)){ continue; }
              
               // ------------ SC matching -------------
               dR = getDeltaR(isotrack.phi() + isotrack.deltaPhi(), isotrack.eta() + isotrack.deltaEta(), photon.phi(), photon.eta());
                      
               if (dR < dRMin){

                   dRMin = dR;
                   scmin = sc;
                   tmin = t;

               }
           }
       }

       

       if (dRMin > dRThreshold){ break; } // Here we go out the while loop


       li = matched_SC.size();
       ElectronCandidate_pt[li] = (*isotracks)[iT.at(tmin)].pt();
       ElectronCandidate_eta[li] = (*isotracks)[iT.at(tmin)].eta();
       ElectronCandidate_phi[li] = (*isotracks)[iT.at(tmin)].phi();
       ElectronCandidate_et[li] = (*photons)[iP.at(scmin)].et();
       ElectronCandidate_photonIdx[li] = scmin;
       ElectronCandidate_isotrackIdx[li] = tmin;
       ElectronCandidate_dxy[li] = IsoTrackSel_dxy[tmin];
       ElectronCandidate_dxyError[li] = IsoTrackSel_dxyError[tmin];
       ElectronCandidate_dxy_PV[li] = IsoTrackSel_dxy_PV[tmin];
       ElectronCandidate_dxyError_PV[li] = IsoTrackSel_dxyError_PV[tmin];
       ElectronCandidate_dxy_BS[li] = IsoTrackSel_dxy_BS[tmin];
       ElectronCandidate_dxyError_BS[li] = IsoTrackSel_dxyError_BS[tmin];
       ElectronCandidate_dR[li] = dRMin;

       // re-compute isolation
       const pat::PackedCandidateRef &e_pck = (*isotracks)[iT.at(tmin)].packedCandRef();
       ElectronCandidate_relPFiso[li] = computeRelIso(*(*e_pck).bestTrack(), packedPFCandidates, lostTracks, true, false);
       ElectronCandidate_relTrkiso[li] = computeRelIso(*(*e_pck).bestTrack(), packedPFCandidates, lostTracks, false, true);
       ElectronCandidate_relTrkiso_noLT[li] = computeRelIso(*(*e_pck).bestTrack(), packedPFCandidates, lostTracks, false, false);
 
       ElectronCandidate_pvAssociationQuality[li] = (*isotracks)[iT.at(tmin)].packedCandRef()->pvAssociationQuality();
       ElectronCandidate_ptDiff[li] = (*isotracks)[iT.at(tmin)].pt() - (*isotracks)[iT.at(tmin)].packedCandRef()->pseudoTrack().pt();

       matched_SC.push_back(scmin); matched_tracks.push_back(tmin);

   }

   nElectronCandidate = matched_SC.size();


   ////////////////////////////////////////////////////////////////////////////////////////////
   //// ---------------------------------------------------------------------------------- ////
   //// -------------------------- LL CANDIDATES RECONSTRUCTION -------------------------- ////                                   
   //// ---------------------------------------------------------------------------------- ////
   ////////////////////////////////////////////////////////////////////////////////////////////

   double minChi2 = 10000;
   int min_i = 99;
   int min_j = 99;

   std::vector<pat::IsolatedTrack> leptonTracks;


   /////////////////////////////////////////////////////////
   // --------------------------------------------------- //
   // ----------- eeCandidates reconstruction ----------- //
   // --------------------------------------------------- //
   /////////////////////////////////////////////////////////
   nEE = 0;
   std::vector<double> pairedE; // electrons that are already paired


   while(2*nEE < nElectronCandidate - 1 ){

     // Init control variables:
     minChi2 = 10000;
     min_i = 99;
     min_j = 99;

     for (int i = 0; i < nElectronCandidate; i++) {
       for (int j = i+1; j < nElectronCandidate; j++) {
 
         if (i == j) { continue; }
         if ( std::find(pairedE.begin(), pairedE.end(), i) != pairedE.end() ) {continue;}
         if ( std::find(pairedE.begin(), pairedE.end(), j) != pairedE.end() ) {continue;}

         const pat::IsolatedTrack & it_i = (*isotracks)[iT.at(ElectronCandidate_isotrackIdx[i])];
         const pat::IsolatedTrack & it_j = (*isotracks)[iT.at(ElectronCandidate_isotrackIdx[j])];
         const pat::PackedCandidateRef &pckCand_i = it_i.packedCandRef(); 
         const pat::PackedCandidateRef &pckCand_j = it_j.packedCandRef(); 
         const reco::Track *itrref_i = (*pckCand_i).bestTrack();
         const reco::Track *itrref_j = (*pckCand_j).bestTrack();
         const reco::Track &itr_i = *itrref_i;
         const reco::Track &itr_j = *itrref_j;

         //llCandidate testcandidate(thePrimaryVertex, theTransientTrackBuilder, it_i, it_j, true);
         trackPair testcandidate(thePrimaryVertex, beamSpotObject, theTransientTrackBuilder, itr_i, itr_j, true);

         if (!testcandidate.hasValidVertex) { continue ;} 

         // Check if the Chi2 is lower:
         if (testcandidate.normalizedChi2 < minChi2) {
           minChi2 = testcandidate.normalizedChi2;
           min_i = i;
           min_j = j;
         }

       } // end j electron loop
     } // end i electron loop

     if (min_i == 99 || min_j == 99) { break; }
     pairedE.push_back(min_i);
     pairedE.push_back(min_j);

     // -> Get LLP Candidate variables:
     const pat::IsolatedTrack &it_A = (*isotracks)[iT.at(ElectronCandidate_isotrackIdx[min_i])];
     const pat::IsolatedTrack &it_B = (*isotracks)[iT.at(ElectronCandidate_isotrackIdx[min_j])]; 
     const pat::PackedCandidateRef &pckCand_A = it_A.packedCandRef(); 
     const pat::PackedCandidateRef &pckCand_B = it_B.packedCandRef(); 
     const reco::Track *itrref_A = (*pckCand_A).bestTrack();
     const reco::Track *itrref_B = (*pckCand_B).bestTrack();
     const reco::Track &itr_A = *itrref_A;
     const reco::Track &itr_B = *itrref_B;


     trackPair eeCandidate(thePrimaryVertex, beamSpotObject, theTransientTrackBuilder, itr_A, itr_B, true);

     // Additionally, for electrons we have to redefine:
     eeCandidate.leadingEt = (ElectronCandidate_et[min_i] > ElectronCandidate_et[min_j])? ElectronCandidate_et[min_i]: ElectronCandidate_et[min_j];
     eeCandidate.subleadingEt = (ElectronCandidate_et[min_i] < ElectronCandidate_et[min_j])? ElectronCandidate_et[min_i]: ElectronCandidate_et[min_j];

     eeCandidate.relisoA = ElectronCandidate_relTrkiso[min_i];
     eeCandidate.relisoB = ElectronCandidate_relTrkiso[min_j];

     eeCandidate.trackDxy = (fabs(ElectronCandidate_dxy[min_i])/ElectronCandidate_dxyError[min_i] < fabs(ElectronCandidate_dxy[min_j])/ElectronCandidate_dxyError[min_j]) ? ElectronCandidate_dxy[min_i] : ElectronCandidate_dxy[min_j];
     eeCandidate.trackIxy = (fabs(ElectronCandidate_dxy[min_i])/ElectronCandidate_dxyError[min_i] < fabs(ElectronCandidate_dxy[min_j])/ElectronCandidate_dxyError[min_j]) ? fabs(ElectronCandidate_dxy[min_i])/ElectronCandidate_dxyError[min_i] : fabs(ElectronCandidate_dxy[min_j])/ElectronCandidate_dxyError[min_j];

     eeCandidate.trackDxy_PV = (fabs(ElectronCandidate_dxy_PV[min_i])/ElectronCandidate_dxyError_PV[min_i] < fabs(ElectronCandidate_dxy_PV[min_j])/ElectronCandidate_dxyError_PV[min_j]) ? ElectronCandidate_dxy_PV[min_i] : ElectronCandidate_dxy_PV[min_j];
     eeCandidate.trackIxy_PV = (fabs(ElectronCandidate_dxy_PV[min_i])/ElectronCandidate_dxyError_PV[min_i] < fabs(ElectronCandidate_dxy_PV[min_j])/ElectronCandidate_dxyError_PV[min_j]) ? fabs(ElectronCandidate_dxy_PV[min_i])/ElectronCandidate_dxyError_PV[min_i] : fabs(ElectronCandidate_dxy_PV[min_j])/ElectronCandidate_dxyError_PV[min_j];
     
     eeCandidate.trackDxy_BS = (fabs(ElectronCandidate_dxy_BS[min_i])/ElectronCandidate_dxyError_BS[min_i] < fabs(ElectronCandidate_dxy_BS[min_j])/ElectronCandidate_dxyError_BS[min_j]) ? ElectronCandidate_dxy_BS[min_i] : ElectronCandidate_dxy_BS[min_j];
     eeCandidate.trackIxy_BS = (fabs(ElectronCandidate_dxy_BS[min_i])/ElectronCandidate_dxyError_BS[min_i] < fabs(ElectronCandidate_dxy_BS[min_j])/ElectronCandidate_dxyError_BS[min_j]) ? fabs(ElectronCandidate_dxy_BS[min_i])/ElectronCandidate_dxyError_BS[min_i] : fabs(ElectronCandidate_dxy_BS[min_j])/ElectronCandidate_dxyError_BS[min_j];

     TLorentzVector l1 = TLorentzVector(); 
     TLorentzVector l2 = TLorentzVector();
     l1.SetPtEtaPhiM(ElectronCandidate_pt[min_i], ElectronCandidate_eta[min_i], ElectronCandidate_phi[min_i], 0.501/1000.0);
     l2.SetPtEtaPhiM(ElectronCandidate_pt[min_j], ElectronCandidate_eta[min_j], ElectronCandidate_phi[min_j], 0.501/1000.0);
     TVector3 vl1 = l1.Vect(); 
     TVector3 vl2 = l2.Vect(); 
     TVector3 vl1l2 = vl1 + vl2;
     TVector3 vec = TVector3(eeCandidate.vx - PV_vx, eeCandidate.vy - PV_vy, 0.0);

     eeCandidate.mass = (l1 + l2).M();
     eeCandidate.dPhi = fabs(vec.DeltaPhi(vl1l2));
     eeCandidate.lldPhi = fabs(l1.DeltaPhi(l2));
     eeCandidate.dR = fabs(l1.DeltaR(l2));


        EE_idxA[nEE] = min_i;
        EE_idxB[nEE] = min_j;
        EE_Lxy_PV[nEE] = eeCandidate.Lxy_PV;
        EE_Ixy_PV[nEE] = eeCandidate.Ixy_PV;
        EE_Lxy_0[nEE] = eeCandidate.Lxy_0;
        EE_Ixy_0[nEE] = eeCandidate.Ixy_0;
        EE_Lxy_BS[nEE] = eeCandidate.Lxy_BS;
        EE_Ixy_BS[nEE] = eeCandidate.Ixy_BS;
        EE_trackDxy[nEE] = eeCandidate.trackDxy;
        EE_trackIxy[nEE] = eeCandidate.trackIxy;
        EE_trackDxy_PV[nEE] = eeCandidate.trackDxy_PV;
        EE_trackIxy_PV[nEE] = eeCandidate.trackIxy_PV;
        EE_trackDxy_BS[nEE] = eeCandidate.trackDxy_BS;
        EE_trackIxy_BS[nEE] = eeCandidate.trackIxy_BS;
        EE_vx[nEE] = eeCandidate.vx;
        EE_vy[nEE] = eeCandidate.vy;
        EE_normalizedChi2[nEE] = eeCandidate.normalizedChi2;
        EE_mass[nEE] = eeCandidate.mass;
        EE_leadingPt[nEE] = eeCandidate.leadingPt;
        EE_subleadingPt[nEE] = eeCandidate.subleadingPt;
        EE_cosAlpha[nEE] = eeCandidate.cosAlpha;
        EE_dPhi[nEE] = eeCandidate.dPhi;
        EE_lldPhi[nEE] = eeCandidate.lldPhi;
        EE_dR[nEE] = eeCandidate.dR;
        EE_relisoA[nEE] = eeCandidate.relisoA;
        EE_relisoB[nEE] = eeCandidate.relisoB;
        EE_leadingEt[nEE] = eeCandidate.leadingEt;
        EE_subleadingEt[nEE] = eeCandidate.subleadingEt;

     nEE++;

   } // end while 



   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////// FILL THE TREE ///////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////
   tree_out->Fill();

}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void displacedElectronAnalyzer::beginJob()
{

    // Output file definition
    output_filename = parameters.getParameter<std::string>("nameOfOutput");
    file_out = new TFile(output_filename.c_str(), "RECREATE");
    //file_out->cd();

    
    // Output Tree definition
    tree_out = new TTree("Events", "Events");

    // Analyzer parameters
    _isData         = parameters.getParameter<bool>("isData");
    _DSAMode        = parameters.getParameter<bool>("DSAMode");
    _doCMSElectrons = parameters.getParameter<bool>("doCMSElectrons");
    _doCMSMuons     = parameters.getParameter<bool>("doCMSMuons");

    _Era     = parameters.getParameter<double>("Era");

    // PU reweighting
    if (_Era == 2016) {
      lumi_weights = edm::LumiReWeighting("2016MCPileupHistogram.root", "2016DataPileupHistogram.root", "pileup", "pileup");
    } else if (_Era == 2017) {
      lumi_weights = edm::LumiReWeighting("2017MCPileupHistogram.root", "2017DataPileupHistogram.root", "pileup", "pileup");
    } else if (_Era == 2018) {
      lumi_weights = edm::LumiReWeighting("2018MCPileupHistogram.root", "2018DataPileupHistogram.root", "pileup", "pileup");
    }


    ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Event_event", &Event_event, "Event_event/I");
    tree_out->Branch("Event_run", &Event_run, "Event_run/I");
    tree_out->Branch("Event_luminosityBlock", &Event_luminosityBlock, "Event_luminosityBlock/I");
    
    tree_out->Branch("nPU", &nPU, "nPU/I");
    tree_out->Branch("nPUTrue", &nPUTrue, "nPUTrue/I");
    tree_out->Branch("wPU", &wPU, "wPU/F");
    tree_out->Branch("genWeight", &genWeight, "genWeight/F");


    ///////////////////////////////// BEAM SPOT BRANCHES ////////////////////////////////

    tree_out->Branch("BeamSpot_x0", &BeamSpot_x0, "BeamSpot_x0/F");
    tree_out->Branch("BeamSpot_y0", &BeamSpot_y0, "BeamSpot_y0/F");
    tree_out->Branch("BeamSpot_z0", &BeamSpot_z0, "BeamSpot_z0/F");


    ////////////////////////////// PRIMARY VERTEX BRANCHES //////////////////////////////

    tree_out->Branch("nPV", &nPV, "nPV/I");
    tree_out->Branch("nTruePV", &nTruePV, "nTruePV/I");
    tree_out->Branch("PV_passAcceptance", &PV_passAcceptance, "PV_passAcceptance/I");
    tree_out->Branch("PV_vx", &PV_vx, "PV_vx/F");
    tree_out->Branch("PV_vy", &PV_vy, "PV_vy/F");
    tree_out->Branch("PV_vz", &PV_vz, "PV_vz/F");


    ///////////////////////////////// ISOTRACK BRANCHES /////////////////////////////////
    
    tree_out->Branch("nIsoTrack", &nIsoTrack, "nIsoTrack/I");
    tree_out->Branch("IsoTrackSel_pt", IsoTrackSel_pt, "IsoTrackSel_pt[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_eta", IsoTrackSel_eta, "IsoTrackSel_eta[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_etaExtra", IsoTrackSel_etaExtra, "IsoTrackSel_etaExtra[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_phi", IsoTrackSel_phi, "IsoTrackSel_phi[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_phiExtra", IsoTrackSel_phiExtra, "IsoTrackSel_phiExtra[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_charge", IsoTrackSel_charge, "IsoTrackSel_charge[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_dxy", IsoTrackSel_dxy, "IsoTrackSel_dxy[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxyError", IsoTrackSel_dxyError, "IsoTrackSel_dxyError[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxy_PV", IsoTrackSel_dxy_PV, "IsoTrackSel_dxy_PV[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxyError_PV", IsoTrackSel_dxyError_PV, "IsoTrackSel_dxyError_PV[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxy_BS", IsoTrackSel_dxy_BS, "IsoTrackSel_dxy_BS[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxyError_BS", IsoTrackSel_dxyError_BS, "IsoTrackSel_dxyError_BS[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dz", IsoTrackSel_dz, "IsoTrackSel_dz[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dzError", IsoTrackSel_dzError, "IsoTrackSel_dzError[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_vx", IsoTrackSel_vx, "IsoTrackSel_vx[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_vy", IsoTrackSel_vy, "IsoTrackSel_vy[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_vz", IsoTrackSel_vz, "IsoTrackSel_vz[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_pfIsolationDR03", IsoTrackSel_pfIsolationDR03, "IsoTrackSel_pfIsolationDR03[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_miniPFIsolation", IsoTrackSel_miniPFIsolation, "IsoTrackSel_miniPFIsolation[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_relPfIsolationDR03", IsoTrackSel_relPfIsolationDR03, "IsoTrackSel_relPfIsolationDR03[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_relMiniPFIsolation", IsoTrackSel_relMiniPFIsolation, "IsoTrackSel_relMiniPFIsolation[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_isHighPurityTrack", IsoTrackSel_isHighPurityTrack, "IsoTrackSel_isHighPurityTrack[nIsoTrack]/I");
    
    tree_out->Branch("IsoTrackSel_numberOfValidTrackerHits", IsoTrackSel_numberOfValidTrackerHits, "IsoTrackSel_numberOfValidTrackerHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidPixelHits", IsoTrackSel_numberOfValidPixelHits, "IsoTrackSel_numberOfValidPixelHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidPixelBarrelHits", IsoTrackSel_numberOfValidPixelBarrelHits, "IsoTrackSel_numberOfValidPixelBarrelHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidPixelEndcapHits", IsoTrackSel_numberOfValidPixelEndcapHits, "IsoTrackSel_numberOfValidPixelEndcapHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripHits", IsoTrackSel_numberOfValidStripHits, "IsoTrackSel_numberOfValidStripHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripTIBHits", IsoTrackSel_numberOfValidStripTIBHits, "IsoTrackSel_numberOfValidStripTIBHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripTIDHits", IsoTrackSel_numberOfValidStripTIDHits, "IsoTrackSel_numberOfValidStripTIDHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripTOBHits", IsoTrackSel_numberOfValidStripTOBHits, "IsoTrackSel_numberOfValidStripTOBHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripTECHits", IsoTrackSel_numberOfValidStripTECHits, "IsoTrackSel_numberOfValidStripTECHits[nIsoTrack]/I");
    

    
    ////////////////////////////////// PHOTON BRANCHES //////////////////////////////////
    
    tree_out->Branch("nPhoton", &nPhoton, "nPhoton/I");
    tree_out->Branch("PhotonSel_et", PhotonSel_et, "PhotonSel_et[nPhoton]/F");
    tree_out->Branch("PhotonSel_eta", PhotonSel_eta, "PhotonSel_eta[nPhoton]/F");
    tree_out->Branch("PhotonSel_phi", PhotonSel_phi, "PhotonSel_phi[nPhoton]/F");
    
    tree_out->Branch("PhotonSel_hadronicOverEm", PhotonSel_hadronicOverEm, "PhotonSel_hadronicOverEm[nPhoton]/F");
    tree_out->Branch("PhotonSel_full5x5_sigmaIetaIeta", PhotonSel_full5x5_sigmaIetaIeta, "PhotonSel_full5x5_sigmaIetaIeta[nPhoton]/F");
    tree_out->Branch("PhotonSel_r9", PhotonSel_r9, "PhotonSel_r9[nPhoton]/F");

    ///////////////////////////////// ELECTRON BRANCHES /////////////////////////////////
    
    if (true) {
       tree_out->Branch("nElectron", &nElectron, "nElectron/I");
       tree_out->Branch("ElectronSel_pt", ElectronSel_pt, "ElectronSel_pt[nElectron]/F");
       tree_out->Branch("ElectronSel_et", ElectronSel_et, "ElectronSel_et[nElectron]/F");
       tree_out->Branch("ElectronSel_eta", ElectronSel_eta, "ElectronSel_eta[nElectron]/F");
       tree_out->Branch("ElectronSel_phi", ElectronSel_phi, "ElectronSel_phi[nElectron]/F");
       tree_out->Branch("ElectronSel_dB", ElectronSel_dB, "ElectronSel_dB[nElectron]/F");
       tree_out->Branch("ElectronSel_edB", ElectronSel_edB, "ElectronSel_edB[nElectron]/F");
       tree_out->Branch("ElectronSel_isLoose", ElectronSel_isLoose, "ElectronSel_isLoose[nElectron]/F");
       tree_out->Branch("ElectronSel_isMedium", ElectronSel_isMedium, "ElectronSel_isMedium[nElectron]/F");
    }
  
    
    //////////////////////////////// GENPARTICLE BRANCHES ///////////////////////////////

    tree_out->Branch("nGenLepton", &nGenLepton, "nGenLepton/I");
    tree_out->Branch("GenLeptonSel_pt", GenLeptonSel_pt, "GenLeptonSel_pt[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_E", GenLeptonSel_E, "GenLeptonSel_E[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_et", GenLeptonSel_et, "GenLeptonSel_et[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_eta", GenLeptonSel_eta, "GenLeptonSel_eta[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_etaAtECAL", GenLeptonSel_etaAtECAL, "GenLeptonSel_etaAtECAL[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_phiAtECAL", GenLeptonSel_phiAtECAL, "GenLeptonSel_phiAtECAL[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_phi", GenLeptonSel_phi, "GenLeptonSel_phi[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_dxy", GenLeptonSel_dxy, "GenLeptonSel_dxy[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_vx", GenLeptonSel_vx, "GenLeptonSel_vx[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_vy", GenLeptonSel_vy, "GenLeptonSel_vy[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_vz", GenLeptonSel_vz, "GenLeptonSel_vz[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_pdgId", GenLeptonSel_pdgId, "GenLeptonSel_pdgId[nGenLepton]/I");
    tree_out->Branch("GenLeptonSel_motherPdgId", GenLeptonSel_motherPdgId, "GenLeptonSel_motherPdgId[nGenLepton]/I");
    tree_out->Branch("GenLeptonSel_isPromptFinalState", GenLeptonSel_isPromptFinalState, "GenLeptonSel_isPromptFinalState[nGenLepton]/I");
    tree_out->Branch("GenLeptonSel_fromHardProcessFinalState", GenLeptonSel_fromHardProcessFinalState, "GenLeptonSel_fromHardProcessFinalState[nGenLepton]/I");
    tree_out->Branch("GenLeptonSel_isDirectPromptTauDecayProductFinalState", GenLeptonSel_isDirectPromptTauDecayProductFinalState, "GenLeptonSel_isDirectPromptTauDecayProductFinalState[nGenLepton]/I");
    tree_out->Branch("GenLeptonSel_isDirectHadronDecayProduct", GenLeptonSel_isDirectHadronDecayProduct, "GenLeptonSel_isDirectHadronDecayProduct[nGenLepton]/I");
    tree_out->Branch("GenLeptonSel_hasRadiated", GenLeptonSel_hasRadiated, "GenLeptonSel_hasRadiated[nGenLepton]/I");

    
    tree_out->Branch("nHardProcessParticle", &nHardProcessParticle, "nHardProcessParticle/I");
    tree_out->Branch("HardProcessParticle_E", HardProcessParticle_E, "HardProcessParticle_E[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_pt", HardProcessParticle_pt, "HardProcessParticle_pt[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_eta", HardProcessParticle_eta, "HardProcessParticle_eta[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_phi", HardProcessParticle_phi, "HardProcessParticle_phi[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_vx", HardProcessParticle_vx, "HardProcessParticle_vx[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_vy", HardProcessParticle_vy, "HardProcessParticle_vy[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_vz", HardProcessParticle_vz, "HardProcessParticle_vz[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_pdgId", HardProcessParticle_pdgId, "HardProcessParticle_pdgId[nHardProcessParticle]/I");

    //////////////////////////// ELECTRON CANDIDATE BRANCHES ////////////////////////////
    
    tree_out->Branch("nElectronCandidate", &nElectronCandidate, "nElectronCandidate/I");
    tree_out->Branch("ElectronCandidate_pt", ElectronCandidate_pt, "ElectronCandidate_pt[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_et", ElectronCandidate_et, "ElectronCandidate_et[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_eta", ElectronCandidate_eta, "ElectronCandidate_eta[nElectronCandidate]/F");    
    tree_out->Branch("ElectronCandidate_phi", ElectronCandidate_phi, "ElectronCandidate_phi[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_photonIdx", ElectronCandidate_photonIdx, "ElectronCandidate_photonIdx[nElectronCandidate]/I");
    tree_out->Branch("ElectronCandidate_isotrackIdx", ElectronCandidate_isotrackIdx, "ElectronCandidate_isotrackIdx[nElectronCandidate]/I");
    tree_out->Branch("ElectronCandidate_dxy", ElectronCandidate_dxy, "ElectronCandidate_dxy[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_dxyError", ElectronCandidate_dxyError, "ElectronCandidate_dxyError[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_dxy_PV", ElectronCandidate_dxy_PV, "ElectronCandidate_dxy_PV[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_dxyError_PV", ElectronCandidate_dxyError_PV, "ElectronCandidate_dxyError_PV[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_dxy_BS", ElectronCandidate_dxy_BS, "ElectronCandidate_dxy_BS[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_dxyError_BS", ElectronCandidate_dxyError_BS, "ElectronCandidate_dxyError_BS[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_relPFiso", ElectronCandidate_relPFiso, "ElectronCandidate_relPFiso[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_relTrkiso", ElectronCandidate_relTrkiso, "ElectronCandidate_relTrkiso[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_relTrkiso_noLT", ElectronCandidate_relTrkiso_noLT, "ElectronCandidate_relTrkiso_noLT[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_dR", ElectronCandidate_dR, "ElectronCandidate_dR[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_pvAssociationQuality", ElectronCandidate_pvAssociationQuality, "ElectronCandidate_pvAssociationQuality[nElectronCandidate]/I");
    

    ////////////////////////////// LL BRANCHES /////////////////////////////

    tree_out->Branch("nEE", &nEE, "nEE/I");

       tree_out->Branch("EE_idxA", EE_idxA, "EE_idxA[nEE]/I");
       tree_out->Branch("EE_idxB", EE_idxB, "EE_idxB[nEE]/I");
       tree_out->Branch("EE_Lxy_PV", EE_Lxy_PV, "EE_Lxy_PV[nEE]/F");
       tree_out->Branch("EE_Ixy_PV", EE_Ixy_PV, "EE_Ixy_PV[nEE]/F");
       tree_out->Branch("EE_Lxy_0", EE_Lxy_0, "EE_Lxy_0[nEE]/F");
       tree_out->Branch("EE_Ixy_0", EE_Ixy_0, "EE_Ixy_0[nEE]/F");
       tree_out->Branch("EE_Lxy_BS", EE_Lxy_BS, "EE_Lxy_BS[nEE]/F");
       tree_out->Branch("EE_Ixy_BS", EE_Ixy_BS, "EE_Ixy_BS[nEE]/F");
       tree_out->Branch("EE_trackDxy", EE_trackDxy, "EE_trackDxy[nEE]/F");
       tree_out->Branch("EE_trackIxy", EE_trackIxy, "EE_trackIxy[nEE]/F");
       tree_out->Branch("EE_trackDxy_PV", EE_trackDxy_PV, "EE_trackDxy_PV[nEE]/F");
       tree_out->Branch("EE_trackIxy_PV", EE_trackIxy_PV, "EE_trackIxy_PV[nEE]/F");
       tree_out->Branch("EE_trackDxy_BS", EE_trackDxy_BS, "EE_trackDxy_BS[nEE]/F");
       tree_out->Branch("EE_trackIxy_BS", EE_trackIxy_BS, "EE_trackIxy_BS[nEE]/F");
       tree_out->Branch("EE_vx", EE_vx, "EE_vx[nEE]/F");
       tree_out->Branch("EE_vy", EE_vy, "EE_vy[nEE]/F");
       tree_out->Branch("EE_mass", EE_mass, "EE_mass[nEE]/F");
       tree_out->Branch("EE_normalizedChi2", EE_normalizedChi2, "EE_normalizedChi2[nEE]/F");
       tree_out->Branch("EE_leadingPt", EE_leadingPt, "EE_leadingPt[nEE]/F");
       tree_out->Branch("EE_subleadingPt", EE_subleadingPt, "EE_subleadingPt[nEE]/F");
       tree_out->Branch("EE_leadingEt", EE_leadingEt, "EE_leadingEt[nEE]/F");
       tree_out->Branch("EE_subleadingEt", EE_subleadingEt, "EE_subleadingEt[nEE]/F");
       tree_out->Branch("EE_cosAlpha", EE_cosAlpha, "EE_cosAlpha[nEE]/F");
       tree_out->Branch("EE_dPhi", EE_dPhi, "EE_dPhi[nEE]/F");
       tree_out->Branch("EE_lldPhi", EE_lldPhi, "EE_lldPhi[nEE]/F");
       tree_out->Branch("EE_dR", EE_dR, "EE_dR[nEE]/F");
       tree_out->Branch("EE_relisoA", EE_relisoA, "EE_relisoA[nEE]/F");
       tree_out->Branch("EE_relisoB", EE_relisoB, "EE_relisoB[nEE]/F");


}
//=======================================================================================================================================================================================================================//

//=======================================================================================================================================================================================================================//
void displacedElectronAnalyzer::endJob() 
{

    file_out->cd();
    tree_out->Write();
    counts->Write();
    sum2Weights->Write();
    file_out->Close();

}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void displacedElectronAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//=======================================================================================================================================================================================================================//

float displacedElectronAnalyzer::getDeltaR(float phi1, float eta1, float phi2, float eta2)
{

    float dPhi = fabs(phi1 - phi2);
    if (dPhi > 3.14) {dPhi = 2*3.14 - dPhi;}
    float dEta = eta1 - eta2;

    float dR = sqrt(dPhi*dPhi + dEta*dEta);

    return dR;

}


bool displacedElectronAnalyzer::passIsotrackSelection( const pat::IsolatedTrack &track) {

   // Quality cuts:
   const reco::HitPattern &hits = track.hitPattern();
   //if (hits.numberOfValidTrackerHits() < 6) { return false; }
   //if (!track.isHighPurityTrack()) { return false;}

   // Isotrack must have packed candidate:
   const pat::PackedCandidateRef &pckCand = track.packedCandRef();
   if (!pckCand.isNonnull()) { return false; }
   if (!pckCand->hasTrackDetails()) { return false; }

   // Preselection cuts:
   //if (track.pt() < 15) { return false; }
   if (fabs(track.eta()) > 2.4) { return false; }

   // To be noticed: Isolation cuts are applied later with the LLCandidate selection

   return true;
}




//=======================================================================================================================================================================================================================//


float displacedElectronAnalyzer::computeDxy(const pat::IsolatedTrack & track, const reco::Vertex pv) {

   double vx = track.vx();
   double vy = track.vy();
   double phi = track.phi();
   double PVx = pv.x();
   double PVy = pv.y();

   double dxy = -(vx - PVx)*sin(phi) + (vy - PVy)*cos(phi);
   return dxy;
}



float displacedElectronAnalyzer::computeDxy(const reco::Track & track, const reco::Vertex pv) {

   double vx = track.vx();
   double vy = track.vy();
   double phi = track.phi();
   double PVx = pv.x();
   double PVy = pv.y();

   double dxy = -(vx - PVx)*sin(phi) + (vy - PVy)*cos(phi);
   return dxy;
}

float displacedElectronAnalyzer::computeDxy(const reco::GenParticle &p, const reco::Vertex &pv)
{

    float vx = p.vx();
    float vy = p.vy();
    float phi = p.phi();

    float pv_x = pv.x();
    float pv_y = pv.y();

    float dxy = -(vx-pv_x)*sin(phi) + (vy-pv_y)*cos(phi);
    return dxy;

}

//=======================================================================================================================================================================================================================//


float displacedElectronAnalyzer::computeDxyError(const pat::IsolatedTrack & track, const reco::Vertex pv) {

   // Trajectory computation [following steps in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTransientTracks#Examples_including_calculation_o]
   // Used to compute uncertainty in the transverse impact parameter with respect to primary vertex.
   
   // Get packedCandidate and transient track:
   const pat::PackedCandidateRef &pCand = track.packedCandRef();
   reco::TransientTrack isotk = theTransientTrackBuilder->build((*pCand).pseudoTrack());
   
   // Define the new point of reference and the trajectory:
   GlobalPoint vert(pv.x(), pv.y(), pv.z());
   TrajectoryStateClosestToPoint  traj = isotk.trajectoryStateClosestToPoint(vert);

   float sigmaXY = traj.perigeeError().transverseImpactParameterError();

   return sigmaXY;

}

//=======================================================================================================================================================================================================================//


float displacedElectronAnalyzer::computeDxyError(const reco::Track & track, const reco::Vertex pv) {

   // Trajectory computation [following steps in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTransientTracks#Examples_including_calculation_o]
   // Used to compute uncertainty in the transverse impact parameter with respect to primary vertex.
   
   // Get packedCandidate and transient track:
   reco::TransientTrack isotk = theTransientTrackBuilder->build(track);
   
   // Define the new point of reference and the trajectory:
   GlobalPoint vert(pv.x(), pv.y(), pv.z());
   TrajectoryStateClosestToPoint  traj = isotk.trajectoryStateClosestToPoint(vert);

   float sigmaXY = traj.perigeeError().transverseImpactParameterError();

   return sigmaXY;

}


//=======================================================================================================================================================================================================================//


float displacedElectronAnalyzer::computeRelIso(const reco::Track & track, edm::Handle<edm::View<pat::PackedCandidate> > pfs, edm::Handle<edm::View<pat::PackedCandidate> > lts, bool isPF, bool includeLostTracks) {


   // Contributions to isolation:
   double charged = 0, neutral = 0, pileup  = 0, trackiso = 0;

   // Loop over packed Candidates
   for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {

      const pat::PackedCandidate &pf = (*pfs)[i];
      
      // Reject pf candidate if it is the same track:
      if (fabs(pf.pt() - track.pt()) < 0.01) { continue; }

      // Only count tracks within a 0.3 cone
      double _dR = getDeltaR(track.phi(), track.eta(), pf.phi(), pf.eta());
      if (_dR > 0.3 || _dR < 0.03) { continue; }

      // PF
      if (pf.charge() == 0) {
         if (pf.pt() > 0.5) neutral += pf.pt();
      } else if (pf.fromPV() >= 2) {
         charged += pf.pt();
      } else {
         if (pf.pt() > 0.5) pileup += pf.pt();
      }

      // track
      if (pf.charge() != 0 and pf.fromPV() >= 2) {trackiso += pf.pt(); }

   }

   // Loop over lost Candidates
   if (includeLostTracks) {

      for (unsigned int i = 0, n = lts->size(); i < n; ++i) {

         const pat::PackedCandidate &lt = (*lts)[i];
      
         // Reject pf candidate if it is the same track:
         if (fabs(lt.pt() - track.pt()) < 0.01) { continue; }

         // Only count tracks within a 0.3 cone
         double _dR = getDeltaR(track.phi(), track.eta(), lt.phi(), lt.eta());
         if (_dR > 0.3 || _dR < 0.03) { continue; }

         // PF
         if (lt.charge() == 0) {
            if (lt.pt() > 0.5) neutral += lt.pt();
         } else if (lt.fromPV() >= 2) {
            charged += lt.pt();
         } else {
            if (lt.pt() > 0.5) pileup += lt.pt();
         }

         // track
         if (lt.charge() != 0 and lt.fromPV() >= 2) {trackiso += lt.pt(); }

      }
   }


   // do deltaBeta:
   double iso = charged + std::max(0.0, neutral-0.5*pileup);
   
   if (isPF){
     return iso/track.pt();
   } else {
     return trackiso/track.pt();
   }
}


TVector3 displacedElectronAnalyzer::extrapolateToECAL(reco::GenParticle gen) {

   Double_t px, py, pz, pt, e, q;
   Double_t x, y, z, t, r, phi;
   Double_t x_c, y_c, r_c, phi_c, phi_0;
   Double_t x_t, y_t, z_t;
   Double_t t1, t2, t3, t4, t5, t6;
   Double_t t_z, t_r, t_ra, t_rb;
   Double_t delta, gammam, omega, asinrho;
   Double_t thefRadius, thefHalfLength;

   const Double_t c_light = 2.99792458E8;
   const Double_t fBz = 3.8;
   

   // set ECAL radius:
   thefRadius = 1.290;
   thefHalfLength = 3.170;

   x = gen.vx()*1.0E-2;
   y = gen.vy()*1.0E-2;
   z = gen.vz()*1.0E-2;

   q = gen.pdgId() > 0 ? -1 : 1;

   // check that particle position is inside the cylinder
   if(TMath::Hypot(x, y) > thefRadius || TMath::Abs(z) > thefHalfLength)
   {
     TVector3 P(0, 0, 0);
     return P;
   }

   px = gen.px();
   py = gen.py();
   pz = gen.pz();
   pt = gen.pt();
   e = gen.energy();

   if (pt*pt < 1.0E-9)
   {
     TVector3 P(0, 0, 0);
     return P;
   }

   // 1. Init particle trajectory
   gammam = e*1.0E9 / (c_light*c_light);      // gammam in [eV/c^2]
   omega = q * fBz / (gammam);                // omega is here in [89875518/s]
   r = pt / (q * fBz) * 1.0E9/c_light;        // in [m]

   phi_0 = TMath::ATan2(py, px); // [rad] in [-pi, pi]

   // 2. helix axis coordinates
   x_c = x + r*TMath::Sin(phi_0);
   y_c = y - r*TMath::Cos(phi_0);
   r_c = TMath::Hypot(x_c, y_c);
   phi_c = TMath::ATan2(y_c, x_c);
   phi = phi_c;
   if(x_c < 0.0) phi += TMath::Pi();

   // 3. Time evaluation
   t_r = 0.0; // in [ns]
   int sign_pz = (pz > 0.0) ? 1 : -1;
   if(pz == 0.0) t_z = 1.0E99;
   else t_z = gammam / (pz*1.0E9/c_light) * (-z + thefHalfLength*sign_pz);

   if(r_c + TMath::Abs(r)  < thefRadius)
   {
     // helix does not cross the cylinder sides
     t = t_z;
   }
   else
   {
     asinrho = TMath::ASin((thefRadius*thefRadius - r_c*r_c - r*r) / (2*TMath::Abs(r)*r_c));
     delta = phi_0 - phi;
     if(delta <-TMath::Pi()) delta += 2*TMath::Pi();
     if(delta > TMath::Pi()) delta -= 2*TMath::Pi();
     t1 = (delta + asinrho) / omega;
     t2 = (delta + TMath::Pi() - asinrho) / omega;
     t3 = (delta + TMath::Pi() + asinrho) / omega;
     t4 = (delta - asinrho) / omega;
     t5 = (delta - TMath::Pi() - asinrho) / omega;
     t6 = (delta - TMath::Pi() + asinrho) / omega;

     if(t1 < 0.0) t1 = 1.0E99;
     if(t2 < 0.0) t2 = 1.0E99;
     if(t3 < 0.0) t3 = 1.0E99;
     if(t4 < 0.0) t4 = 1.0E99;
     if(t5 < 0.0) t5 = 1.0E99;
     if(t6 < 0.0) t6 = 1.0E99;

     t_ra = TMath::Min(t1, TMath::Min(t2, t3));
     t_rb = TMath::Min(t4, TMath::Min(t5, t6));
     t_r = TMath::Min(t_ra, t_rb);
     t = TMath::Min(t_r, t_z);
   }

   // 4. position in terms of x(t), y(t), z(t)
   x_t = x_c + r * TMath::Sin(omega * t - phi_0);
   y_t = y_c + r * TMath::Cos(omega * t - phi_0);
   z_t = z + pz*1.0E9 / c_light / gammam * t;

   TVector3 P(x_t, y_t, z_t);

   return P;


}


DEFINE_FWK_MODULE(displacedElectronAnalyzer);
