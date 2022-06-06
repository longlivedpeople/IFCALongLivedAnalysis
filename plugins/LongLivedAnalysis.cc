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


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// FUNCTIONS ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

float dxy_value(const reco::GenParticle &p, const reco::Vertex &pv)
{

    float vx = p.vx();
    float vy = p.vy();
    float phi = p.phi();

    float pv_x = pv.x();
    float pv_y = pv.y();
  
    float dxy = -(vx-pv_x)*sin(phi) + (vy-pv_y)*cos(phi);
    return dxy;

}


bool isLongLivedLepton(const reco::GenParticle &p)
{

    bool LLP = false;

    if (!( abs(p.pdgId()) == 11 || abs(p.pdgId()) == 13)){ return false; }
    if (abs(p.mother()->pdgId()) == 1000022 or abs(p.mother()->pdgId()) == 54){ LLP = true; }

    return LLP;

}



float getDeltaR(float phi1, float eta1, float phi2, float eta2)
{

    float dPhi = fabs(phi1 - phi2);
    if (dPhi > 3.14) {dPhi = 2*3.14 - dPhi;}
    float dEta = eta1 - eta2;

    float dR = sqrt(dPhi*dPhi + dEta*dEta);

    return dR;

}



/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// DATA DEFINITION //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////// BRANCHES /////////////////////////////////////

//-> EVENT INFO
Int_t Event_event;
Int_t Event_luminosityBlock;
Int_t Event_run;
Int_t nPU;
Int_t nPUTrue;
Float_t wPU;
Float_t genWeight;

//-> TRIGGER TAGS
std::vector<std::string> HLTPaths_;
bool triggerPass[200] = {false};

//-> PRIMARY VERTEX SELECTION
Int_t nPV;
Int_t nTruePV;
Int_t PV_passAcceptance;
Float_t PV_vx;
Float_t PV_vy;
Float_t PV_vz;


//-> BEAM SPOT
Float_t BeamSpot_x0;
Float_t BeamSpot_y0;
Float_t BeamSpot_z0;
Float_t BeamSpot_BeamWidthX;
Float_t BeamSpot_BeamWidthY;



//-> ISOTRACK SELECTION
const Int_t nIsoTrackMax = 500;
Int_t nIsoTrack;
// Primitive:
Float_t IsoTrackSel_pt[nIsoTrackMax];
Float_t IsoTrackSel_eta[nIsoTrackMax];
Float_t IsoTrackSel_etaExtra[nIsoTrackMax];
Float_t IsoTrackSel_phiExtra[nIsoTrackMax];
Float_t IsoTrackSel_phi[nIsoTrackMax];
Int_t IsoTrackSel_charge[nIsoTrackMax];
Float_t IsoTrackSel_dxy[nIsoTrackMax];
Float_t IsoTrackSel_dxyError[nIsoTrackMax];
Float_t IsoTrackSel_dxy_PV[nIsoTrackMax];
Float_t IsoTrackSel_dxyError_PV[nIsoTrackMax];
Float_t IsoTrackSel_dxy_BS[nIsoTrackMax];
Float_t IsoTrackSel_dxyError_BS[nIsoTrackMax];
Float_t IsoTrackSel_dz[nIsoTrackMax];
Float_t IsoTrackSel_dzError[nIsoTrackMax];
Float_t IsoTrackSel_vx[nIsoTrackMax];
Float_t IsoTrackSel_vy[nIsoTrackMax];
Float_t IsoTrackSel_vz[nIsoTrackMax];
Float_t IsoTrackSel_pfIsolationDR03[nIsoTrackMax];
Float_t IsoTrackSel_miniPFIsolation[nIsoTrackMax];
Float_t IsoTrackSel_relPfIsolationDR03[nIsoTrackMax];
Float_t IsoTrackSel_relMiniPFIsolation[nIsoTrackMax];
Int_t IsoTrackSel_isHighPurityTrack[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidTrackerHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidPixelHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidPixelBarrelHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidPixelEndcapHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidStripHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidStripTIBHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidStripTIDHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidStripTOBHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidStripTECHits[nIsoTrackMax];
Int_t IsoTrackSel_fromPV[nIsoTrackMax];
Float_t IsoTrackSel_PVx[nIsoTrackMax];
Float_t IsoTrackSel_PVy[nIsoTrackMax];
Float_t IsoTrackSel_PVz[nIsoTrackMax];


//-> PHOTON SELECTION
const Int_t nPhotonMax = 100;
Int_t nPhoton;
Float_t PhotonSel_et[nPhotonMax];
Float_t PhotonSel_eta[nPhotonMax];
Float_t PhotonSel_phi[nPhotonMax];
Float_t PhotonSel_hadronicOverEm[nPhotonMax];
Float_t PhotonSel_full5x5_sigmaIetaIeta[nPhotonMax];
Int_t PhotonSel_isEB[nPhotonMax];
Int_t PhotonSel_isEE[nPhotonMax];
Float_t PhotonSel_r9[nPhotonMax];
Float_t PhotonSel_ecalIso[nPhotonMax];
Float_t PhotonSel_hcalIso[nPhotonMax];
Float_t PhotonSel_caloIso[nPhotonMax];
Float_t PhotonSel_relIso[nPhotonMax];

//-> ELECTRON SELECTION
const Int_t nElectronMax = 100;
Int_t nElectron;
Float_t ElectronSel_pt[nElectronMax];
Float_t ElectronSel_et[nElectronMax];
Float_t ElectronSel_eta[nElectronMax];
Float_t ElectronSel_phi[nElectronMax];
Float_t ElectronSel_dxy[nElectronMax];
Float_t ElectronSel_dxyError[nElectronMax];
Float_t ElectronSel_dxySignificance[nElectronMax];
Float_t ElectronSel_dB[nElectronMax];
Float_t ElectronSel_edB[nElectronMax];
Float_t ElectronSel_isLoose[nElectronMax];
Float_t ElectronSel_isMedium[nElectronMax];
Float_t ElectronSel_isTight[nElectronMax];


// -> MUON SELECTION
const Int_t nMuonMax = 100;
Int_t nMuon;
Float_t MuonSel_pt[nMuonMax];
Float_t MuonSel_eta[nMuonMax];
Float_t MuonSel_phi[nMuonMax];
Float_t MuonSel_relIso[nMuonMax];
Float_t MuonSel_dB[nMuonMax];
Float_t MuonSel_edB[nMuonMax];
Float_t MuonSel_dBSignificance[nMuonMax];
Int_t MuonSel_isMuon[nMuonMax];
Int_t MuonSel_isGlobalMuon[nMuonMax];
Int_t MuonSel_isTrackerMuon[nMuonMax];
Int_t MuonSel_isStandAloneMuon[nMuonMax];
Int_t MuonSel_isLooseMuon[nMuonMax];
Int_t MuonSel_isMediumMuon[nMuonMax];
Int_t MuonSel_isGoodMediumMuon[nMuonMax];
Int_t MuonSel_isPFMuon[nMuonMax];

// -> AOD MUON COLLECTIONS
// DisplacedStandAlone:
Int_t nDSA;
Float_t DSA_pt[200];
Float_t DSA_eta[200];
Float_t DSA_phi[200];
Float_t DSA_dxy[200];
Int_t DSA_q[200];
// DisplacedGlobalMuons 
Int_t nDGM;
Int_t DGM_idx[200];
Float_t DGM_pt[200];
Float_t DGM_ptError[200];
Float_t DGM_eta[200];
Float_t DGM_phi[200];
Float_t DGM_dxy[200];
Float_t DGM_dz[200];
Float_t DGM_dxyError[200];
Float_t DGM_dxy_PV[200];
Float_t DGM_dxyError_PV[200];
Float_t DGM_dxy_BS[200];
Float_t DGM_dxyError_BS[200];
Float_t DGM_relPFiso[200];
Float_t DGM_relPFiso_noLT[200];
Float_t DGM_normChi2[200];
Int_t DGM_numberOfValidHits[200];
Int_t DGM_numberOfLostHits[200];
Int_t DGM_muonHits[200];
Int_t DGM_outerTrackerHits[200];
Int_t DGM_trackerHits[200];
Int_t DGM_CSCHits[200];
Int_t DGM_DTHits[200];
Int_t DGM_charge[200];



//-> GENLEPTON SELECTION
const Int_t nGenLeptonMax = 30;
Int_t nGenLepton;
Int_t nGenLepton_PFS;
Int_t nGenLepton_HPFS;
Int_t nGenLepton_PTDP;
Int_t nGenLepton_HDP;
Float_t GenLeptonSel_pt[nGenLeptonMax];
Float_t GenLeptonSel_E[nGenLeptonMax];
Float_t GenLeptonSel_et[nGenLeptonMax];
Float_t GenLeptonSel_eta[nGenLeptonMax];
Float_t GenLeptonSel_phi[nGenLeptonMax];
Int_t GenLeptonSel_pdgId[nGenLeptonMax];
Float_t GenLeptonSel_dxy[nGenLeptonMax];
Float_t GenLeptonSel_vx[nGenLeptonMax];
Float_t GenLeptonSel_vy[nGenLeptonMax];
Float_t GenLeptonSel_vz[nGenLeptonMax];
Int_t GenLeptonSel_motherPdgId[nGenLeptonMax];
Int_t GenLeptonSel_fromHardProcessFinalState[nGenLeptonMax];
Int_t GenLeptonSel_isPromptFinalState[nGenLeptonMax];
Int_t GenLeptonSel_isDirectPromptTauDecayProductFinalState[nGenLeptonMax];
Int_t GenLeptonSel_isDirectHadronDecayProduct[nGenLeptonMax];

Int_t nHardProcessParticle;
Float_t HardProcessParticle_pt[30];
Float_t HardProcessParticle_E[30];
Float_t HardProcessParticle_eta[30];
Float_t HardProcessParticle_phi[30];
Float_t HardProcessParticle_vx[30];
Float_t HardProcessParticle_vy[30];
Float_t HardProcessParticle_vz[30];
Int_t HardProcessParticle_pdgId[30];


// -> GENERATION ACCEPTANCE CRITERIA
Bool_t passAcceptanceCriteria;


// -> MET 

Int_t MET_isPFMET;
Float_t MET_pt; // default Type 1 corrected MET
Float_t MET_phi;

//-> ELECTRON CANDIDATE SELECTION
const Int_t nElectronCandidateMax = 1000;
Int_t nElectronCandidate;
Float_t ElectronCandidate_pt[nElectronCandidateMax];
Float_t ElectronCandidate_et[nElectronCandidateMax];
Float_t ElectronCandidate_eta[nElectronCandidateMax];
Float_t ElectronCandidate_phi[nElectronCandidateMax];
Float_t ElectronCandidate_dxy[nElectronCandidateMax];
Float_t ElectronCandidate_dxyError[nElectronCandidateMax];
Float_t ElectronCandidate_dxy_PV[nElectronCandidateMax];
Float_t ElectronCandidate_dxyError_PV[nElectronCandidateMax];
Float_t ElectronCandidate_dxy_BS[nElectronCandidateMax];
Float_t ElectronCandidate_dxyError_BS[nElectronCandidateMax];
Float_t ElectronCandidate_relPFiso[nElectronCandidateMax];
Float_t ElectronCandidate_relTrkiso[nElectronCandidateMax];
Float_t ElectronCandidate_relTrkiso_noLT[nElectronCandidateMax];
Int_t ElectronCandidate_photonIdx[nElectronCandidateMax];
Int_t ElectronCandidate_isotrackIdx[nElectronCandidateMax];
Int_t ElectronCandidate_pvAssociationQuality[nElectronCandidateMax];
Float_t ElectronCandidate_ptDiff[nElectronCandidateMax];

//-> MUON CANDIDATE SELECTION
const Int_t nMuonCandidateMax = 100;
Int_t nMuonCandidate;
Float_t MuonCandidate_pt[nMuonCandidateMax];
Float_t MuonCandidate_eta[nMuonCandidateMax];
Float_t MuonCandidate_phi[nMuonCandidateMax];
Float_t MuonCandidate_dxy[nMuonCandidateMax];
Float_t MuonCandidate_dxyError[nMuonCandidateMax];
Float_t MuonCandidate_triggerPt[nMuonCandidateMax];
Int_t MuonCandidate_muonTriggerObjectIdx[nMuonCandidateMax];
Int_t MuonCandidate_isotrackIdx[nMuonCandidateMax];
Int_t MuonCandidate_pvAssociationQuality[nMuonCandidateMax];
Float_t MuonCandidate_ptDiff[nMuonCandidateMax];


// -> All EE candidates
Int_t nEE;
Int_t EE_idxA[20];
Int_t EE_idxB[20];
Float_t EE_Lxy_PV[20];
Float_t EE_Ixy_PV[20];
Float_t EE_Lxy_0[20];
Float_t EE_Ixy_0[20];
Float_t EE_Lxy_BS[20];
Float_t EE_Ixy_BS[20];
Float_t EE_trackDxy[20];
Float_t EE_trackIxy[20];
Float_t EE_trackDxy_PV[20];
Float_t EE_trackIxy_PV[20];
Float_t EE_trackDxy_BS[20];
Float_t EE_trackIxy_BS[20];
Float_t EE_vx[20];
Float_t EE_vy[20];
Float_t EE_mass[20];
Float_t EE_normalizedChi2[20];
Float_t EE_leadingPt[20];
Float_t EE_subleadingPt[20];
Float_t EE_leadingEt[20];
Float_t EE_subleadingEt[20];
Float_t EE_cosAlpha[20];
Float_t EE_dR[20];
Float_t EE_dPhi[20];
Float_t EE_lldPhi[20];
Float_t EE_relisoA[20];
Float_t EE_relisoB[20];


// -> All DGM pairs
Int_t nDMDM;
Int_t DMDM_idxA[20];
Int_t DMDM_idxB[20];
Float_t DMDM_Lxy_PV[20];
Float_t DMDM_Ixy_PV[20];
Float_t DMDM_Lxy_0[20];
Float_t DMDM_Ixy_0[20];
Float_t DMDM_Lxy_BS[20];
Float_t DMDM_Ixy_BS[20];
Float_t DMDM_trackDxy[20];
Float_t DMDM_trackIxy[20];
Float_t DMDM_trackDxy_PV[20];
Float_t DMDM_trackIxy_PV[20];
Float_t DMDM_trackDxy_BS[20];
Float_t DMDM_trackIxy_BS[20];
Float_t DMDM_vx[20];
Float_t DMDM_vy[20];
Float_t DMDM_normalizedChi2[20];
Float_t DMDM_mass[20];
Float_t DMDM_leadingPt[20];
Float_t DMDM_subleadingPt[20];
Float_t DMDM_cosAlpha[20];
Float_t DMDM_dPhi[20];
Float_t DMDM_lldPhi[20];
Float_t DMDM_dR[20];
Float_t DMDM_relisoA[20];
Float_t DMDM_relisoB[20];



/////////////////////////////////////// OUTPUT //////////////////////////////////////

TFile *file_out;
TTree *tree_out;


//=======================================================================================================================================================================================================================//
class LongLivedAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit LongLivedAnalysis(const edm::ParameterSet&);
      ~LongLivedAnalysis();

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

      // Filter variables
      bool _filterByEE;
      bool _filterByMM;
      bool _filterByLL;
      bool _filterByElectron;
      bool _filterByMuon;

      edm::ParameterSet parameters;

      // Tokens
      edm::EDGetTokenT<edm::View<pat::Electron> > theElectronCollection;   
      edm::EDGetTokenT<edm::View<pat::Muon> > theMuonCollection;   
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

      // Displaced objects (AOD) tokens:
      edm::EDGetTokenT<std::vector<reco::Track> > theDSACollection; // DisplacedStandAlone
      edm::EDGetTokenT<std::vector<reco::Track> > theDGMCollection; // DisplacedGlobalMuons


      // Gen collection
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  theGenParticleCollection;      
      edm::EDGetTokenT<GenEventInfoProduct>  theGenEventInfoProduct;

      // PU reweighting
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  thePileUpSummary;
      edm::LumiReWeighting lumi_weights;

      // "Global" variables
      std::vector<int> iT; // track indexes
      edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;


      // Class functions
      std::string getPathVersion(const edm:: TriggerNames &names, const std::string &rawPath);
      bool passIsotrackSelection(const pat::IsolatedTrack &track);
      bool passPhotonSelection(const pat::Photon &photon);
      bool passMuonSelection(const pat::Muon &muon);
      bool passDGMSelection(const reco::Track &muon);
      float computeDxy(const pat::IsolatedTrack & track, const reco::Vertex pv);
      float computeDxy(const reco::Track & track, const reco::Vertex pv);
      reco::Vertex getSVCandidate(const pat::PackedCandidateRef &pckCandA, const pat::PackedCandidateRef &pckCandB);
      float computeDxyError(const pat::IsolatedTrack & track, const reco::Vertex pv);
      float computeDxyError(const reco::Track & track, const reco::Vertex pv);
      float computeRelIso(const reco::Track & track,  edm::Handle<edm::View<pat::PackedCandidate> > pfs, edm::Handle<edm::View<pat::PackedCandidate> > lts, bool isPF, bool includeLostTracks);

      // Histograms
      TH1F *counts, *sum2Weights;
};
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
LongLivedAnalysis::LongLivedAnalysis(const edm::ParameterSet& iConfig)
{
   usesResource("TFileService");
  
   parameters = iConfig;

   counts = new TH1F("counts", "", 1, 0, 1);
   sum2Weights = new TH1F("sum2Weights", "", 1, 0, 1);

   theElectronCollection = consumes<edm::View<pat::Electron> >  (parameters.getParameter<edm::InputTag>("ElectronCollection"));
   theMuonCollection = consumes<edm::View<pat::Muon> >  (parameters.getParameter<edm::InputTag>("MuonCollection"));
   thePhotonCollection = consumes<edm::View<pat::Photon> > (parameters.getParameter<edm::InputTag>("PhotonCollection"));
   theIsoTrackCollection = consumes<edm::View<pat::IsolatedTrack> >  (parameters.getParameter<edm::InputTag>("IsoTrackCollection"));
   thePrimaryVertexCollection = consumes<edm::View<reco::Vertex> >  (parameters.getParameter<edm::InputTag>("PrimaryVertexCollection"));
   thePackedPFCandidateCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("PackedPFCandidateCollection"));
   theLostTracksCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("LostTracksCollection"));
   theEleLostTracksCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("EleLostTracksCollection"));
   theMETCollection = consumes<edm::View<pat::MET> >  (parameters.getParameter<edm::InputTag>("METCollection"));


   triggerBits_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("bits"));
   triggerObjects_ = consumes<edm::View<pat::TriggerObjectStandAlone> > (parameters.getParameter<edm::InputTag>("objects"));
   triggerPrescales_ = consumes<pat::PackedTriggerPrescales > (parameters.getParameter<edm::InputTag>("prescales"));

   theBeamSpot = consumes<reco::BeamSpot>  (parameters.getParameter<edm::InputTag>("BeamSpot"));


   theGenParticleCollection = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genParticleCollection"));

   theGenEventInfoProduct = consumes<GenEventInfoProduct> (parameters.getParameter<edm::InputTag>("theGenEventInfoProduct"));

   thePileUpSummary = consumes<std::vector<PileupSummaryInfo> > (parameters.getParameter<edm::InputTag>("thePileUpSummary"));

   
   theDGMCollection = consumes<std::vector<reco::Track> >  (parameters.getParameter<edm::InputTag>("DisplacedGlobalMuonCollection"));




}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
LongLivedAnalysis::~LongLivedAnalysis()
{

}
//=======================================================================================================================================================================================================================//



//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////// MAIN CODE /////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////



   //////////////////////////////// GET THE COLLECTIONS ////////////////////////////////
   
   // Handle:
   edm::Handle<edm::View<pat::Electron> > electrons;
   edm::Handle<edm::View<pat::Muon> > muons;
   edm::Handle<edm::View<pat::Photon> > photons;
   edm::Handle<edm::View<pat::IsolatedTrack> > isotracks;
   edm::Handle<edm::View<reco::Vertex> > primaryvertices;
   edm::Handle<edm::View<pat::PackedCandidate> > packedPFCandidates;
   edm::Handle<edm::View<pat::PackedCandidate> > lostTracks;
   edm::Handle<edm::View<pat::PackedCandidate> > eleLostTracks;
   edm::Handle<edm::View<pat::MET> > METs;
   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<edm::View<pat::TriggerObjectStandAlone>  >triggerObjects;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   edm::Handle<reco::BeamSpot> beamSpot;
   edm::Handle<edm::View<reco::GenParticle> > genParticles;
   edm::Handle<GenEventInfoProduct> genEvtInfo;
   edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
   
   edm::Handle<std::vector<reco::Track> > DGMs;

   // Get tokens:
   iEvent.getByToken(theElectronCollection, electrons);
   iEvent.getByToken(theMuonCollection, muons);
   iEvent.getByToken(thePhotonCollection, photons);
   iEvent.getByToken(theIsoTrackCollection, isotracks);
   iEvent.getByToken(thePrimaryVertexCollection, primaryvertices);
   iEvent.getByToken(thePackedPFCandidateCollection, packedPFCandidates);
   iEvent.getByToken(theLostTracksCollection, lostTracks);
   iEvent.getByToken(theEleLostTracksCollection, eleLostTracks);
   iEvent.getByToken(theMETCollection, METs);
   iEvent.getByToken(triggerBits_, triggerBits);
   iEvent.getByToken(triggerObjects_, triggerObjects);
   iEvent.getByToken(theBeamSpot, beamSpot);

   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);

   
   if (_DSAMode){
      iEvent.getByToken(theDGMCollection, DGMs);
   }
   



   if (!_isData){
      
       iEvent.getByToken(theGenParticleCollection, genParticles);    
       iEvent.getByToken(theGenEventInfoProduct, genEvtInfo);
       iEvent.getByToken(thePileUpSummary, puInfoH);


   }
   

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


   
   //// -------------------------
   //// --
   //// ---- Trigger acceptance 
   //// --
   //// -------------------------



   // Check if trigger fired:
   //
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   unsigned int ipath = 0;
   for (auto path : HLTPaths_) {
     std::string path_v = path + "_v";
     bool fired = false;
     for (unsigned int itrg = 0; itrg < triggerBits->size(); ++itrg) {
       TString TrigPath = names.triggerName(itrg);
       if (!triggerBits->accept(itrg))
         continue;
       if (!TrigPath.Contains(path_v)){
         continue;
       }
       //std::cout << path << "\t" << TrigPath << std::endl;
       fired = true;
     }
     triggerPass[ipath] = fired;
     ipath++;
   }

    
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



   //// ----------------------------
   //// --
   //// ---- CMS Photon Collection
   //// --
   //// ----------------------------
   
   std::vector<int> iP; // photon indexes


   // Select good photons
   for (size_t i = 0; i < photons->size(); i++){

       const pat::Photon & photon = (*photons)[i];
       
       if (!passPhotonSelection(photon)) { continue;}

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
    

   if (_doCMSElectrons) {
  
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
   


   //// --------------------------
   //// --
   //// ---- CMS Muon Collection
   //// --
   //// --------------------------

   if (_doCMSMuons) {

     std::vector<int> iM; // muon indexes


     // Select good muons
     for (size_t i = 0; i < muons->size(); i++){

       const pat::Muon & muon = (*muons)[i];
       
       if (fabs(muon.eta()) < 2.0) { iM.push_back(i);}

     }

     // Sort good lepton indexes by pt
     std::sort( std::begin(iM), std::end(iM), [&](int i1, int i2){ return muons->at(i1).pt() > muons->at(i2).pt(); });

     nMuon = iM.size();

     for (size_t i = 0; i < iM.size(); i++){

       const pat::Muon & muon = (* muons)[iM.at(i)];

       // Kinematics:
       MuonSel_pt[i] = muon.pt();
       MuonSel_eta[i] = muon.eta();
       MuonSel_phi[i] = muon.phi();

       // Isolation
       // Defined as PFIsolation in a cone R < 0.3
       const reco::MuonPFIsolation & muonpfiso = muon.pfIsolationR03();
       float muonNeutralIso = fmax(0.0, muonpfiso.sumNeutralHadronEt + muonpfiso.sumPhotonEt - 0.5*muonpfiso.sumPUPt);
       float muonChargedIso = muonpfiso.sumChargedHadronPt;
       MuonSel_relIso[i] = (muonNeutralIso + muonChargedIso)/muon.pt();

       MuonSel_isMuon[i] = muon.isMuon();
       MuonSel_isGlobalMuon[i] = muon.isGlobalMuon();
       MuonSel_isTrackerMuon[i] = muon.isTrackerMuon();
       MuonSel_isStandAloneMuon[i] = muon.isStandAloneMuon();
       MuonSel_isLooseMuon[i] = muon.isLooseMuon();
       MuonSel_isMediumMuon[i] = muon.isMediumMuon();
       MuonSel_isPFMuon[i] = muon.isPFMuon();

       MuonSel_dB[i] = muon.dB();
       MuonSel_edB[i] = muon.edB();
       MuonSel_dBSignificance[i] = muon.dB()/muon.edB();

     }

   }

   //// ---------------------------
   //// --
   //// ---- AOD Muon Collections
   //// --
   //// ---------------------------

   std::vector<int> iDGM; // needs to be declared outside the conditional
   nDGM = 0; // initialize for avoiding crushing

   if (_DSAMode){


     //// -----------------------------
     //// ---- Displaced Global Muons
     //// -----------------------------

     for (size_t i = 0; i < DGMs->size(); i++){

       const reco::Track &muon = (*DGMs)[i];
       if (passDGMSelection(muon)) {iDGM.push_back(i);}

     }
     nDGM = iDGM.size();

     for (size_t i = 0; i < iDGM.size(); i++){

       const reco::Track &muon = (*DGMs)[iDGM.at(i)];

       // Main variables
       DGM_pt[i] = muon.pt();
       DGM_ptError[i] = muon.ptError();
       DGM_eta[i] = muon.eta();
       DGM_phi[i] = muon.phi();
       DGM_relPFiso[i] = computeRelIso(muon, packedPFCandidates, lostTracks, true, true);
       DGM_relPFiso_noLT[i] = computeRelIso(muon, packedPFCandidates, lostTracks, true, false);
       DGM_numberOfValidHits[i] = muon.numberOfValidHits();
       DGM_numberOfLostHits[i] = muon.numberOfLostHits();
       DGM_muonHits[i] = muon.hitPattern().numberOfValidMuonHits();
       DGM_outerTrackerHits[i] = muon.hitPattern().numberOfValidStripHits();
       DGM_trackerHits[i] = muon.hitPattern().numberOfValidTrackerHits();
       DGM_CSCHits[i] = muon.hitPattern().numberOfValidMuonCSCHits();
       DGM_DTHits[i] = muon.hitPattern().numberOfValidMuonDTHits();
       DGM_charge[i] = muon.charge();
       DGM_normChi2[i] = muon.normalizedChi2();

       // Transverse impact parameters
       reco::TransientTrack _dgmtk = theTransientTrackBuilder->build(muon);
       TrajectoryStateClosestToPoint _trajPV = _dgmtk.trajectoryStateClosestToPoint( _PVpoint );
       TrajectoryStateClosestToPoint _trajBS = _dgmtk.trajectoryStateClosestToPoint( _BSpoint );
       DGM_dxy[i] = muon.dxy();
       DGM_dz[i] = muon.dz();
       DGM_dxyError[i] = muon.dxyError();
       DGM_dxy_PV[i] = -_trajPV.perigeeParameters().transverseImpactParameter();
       DGM_dxyError_PV[i] = _trajPV.perigeeError().transverseImpactParameterError();
       DGM_dxy_BS[i] = -_trajBS.perigeeParameters().transverseImpactParameter();
       DGM_dxyError_BS[i] = _trajBS.perigeeError().transverseImpactParameterError();


       DGM_idx[i] = iDGM.at(i); // only to use in analyzer not in Galapago

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

           if ( ( abs(genparticle.pdgId()) == 11 || abs(genparticle.pdgId()) == 13 ) && genparticle.status() == 1) {
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

              mref = genparticle.motherRef();
              m = *mref;
              while (m.pdgId() == m.mother()->pdgId()) {
                  mref = m.motherRef();
                  m = *mref;
              }

              GenLeptonSel_vx[i] = m.vx();
              GenLeptonSel_vy[i] = m.vy();
              GenLeptonSel_vz[i] = m.vz();
	      GenLeptonSel_dxy[i] = dxy_value(m, thePrimaryVertex); // should be computed here or before?

              if(m.numberOfMothers() != 0){
                  GenLeptonSel_motherPdgId[i] = m.motherRef()->pdgId();
              } else {
                  GenLeptonSel_motherPdgId[i] = 0; 
              }
          }else{

              GenLeptonSel_vx[i] = genparticle.vx();
              GenLeptonSel_vy[i] = genparticle.vy();
              GenLeptonSel_vz[i] = genparticle.vz();
	      GenLeptonSel_dxy[i] = dxy_value(genparticle, thePrimaryVertex); // should be computed here or before?

              GenLeptonSel_motherPdgId[i] = genparticle.motherRef()->pdgId();
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


   //// ----------------------
   //// --
   //// ---- MET Collections
   //// --
   //// ----------------------

   //CHECK: DO WE WANT TO USE PUPPIMET
   const pat::MET &met = (*METs)[0]; // access the MET object

   MET_pt = met.pt();
   MET_phi = met.phi();


   //// -----------------------------------
   //// --
   //// ---- Displaced Electrons reconstruction
   //// --
   //// -----------------------------------


   // Variable initiallization:

   std::vector<int> matched_tracks, matched_SC, matched_triggerObjects; // std vectors with matched objects to avoid overlapping

   float dRMin = 99999; // dR to minimize as high as possible in the beginning
   float dRThreshold = 0.1; // Maximum dR to do the lepton matching
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
   //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);

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




   /////////////////////////////////////////////////////////////
   // ------------------------------------------------------- //
   // ----------- dGMdGMCandidates reconstruction ----------- //
   // ------------------------------------------------------- //
   /////////////////////////////////////////////////////////////
   nDMDM = 0;
   std::vector<double> pairedDM; // displaced global muons that are already paired


   while(2*nDMDM < nDGM - 1 ){

     // Init control variables:
     minChi2 = 10000;
     min_i = 99;
     min_j = 99;

     for (int i = 0; i < nDGM; i++) {
       for (int j = i+1; j < nDGM; j++) {
 
         if (i == j) { continue; }
         if ( std::find(pairedDM.begin(), pairedDM.end(), i) != pairedDM.end() ) {continue;}
         if ( std::find(pairedDM.begin(), pairedDM.end(), j) != pairedDM.end() ) {continue;}

         const reco::Track & tr_i = (*DGMs)[DGM_idx[i]];
         const reco::Track & tr_j = (*DGMs)[DGM_idx[j]];

         if ( deltaR(tr_i.eta(), tr_i.phi(), tr_j.eta(), tr_j.phi()) < 0.1 ) // Ghost suppresion (to be optimized)
           continue; 

         trackPair testcandidate(thePrimaryVertex, beamSpotObject, theTransientTrackBuilder, tr_i, tr_j, false);

         if (!testcandidate.hasValidVertex) { continue ;} 

         // Check if the Chi2 is lower:
         if (testcandidate.normalizedChi2 < minChi2) {
           minChi2 = testcandidate.normalizedChi2;
           min_i = i;
           min_j = j;
         }

       } // end j muon loop
     } // end i muon loop

     if (min_i == 99 || min_j == 99) { break; }
     pairedDM.push_back(min_i);
     pairedDM.push_back(min_j);

     // -> Get LLP Candidate variables:
     const reco::Track & tr_i = (*DGMs)[DGM_idx[min_i]];
     const reco::Track & tr_j = (*DGMs)[DGM_idx[min_j]];

     trackPair dmdmCandidate(thePrimaryVertex, beamSpotObject, theTransientTrackBuilder, tr_i, tr_j, false);

     dmdmCandidate.relisoA = DGM_relPFiso[min_i];
     dmdmCandidate.relisoB = DGM_relPFiso[min_j];

     dmdmCandidate.trackDxy = (fabs(DGM_dxy[min_i])/DGM_dxyError[min_i] < fabs(DGM_dxy[min_j])/DGM_dxyError[min_j]) ? DGM_dxy[min_i] : DGM_dxy[min_j];
     dmdmCandidate.trackIxy = (fabs(DGM_dxy[min_i])/DGM_dxyError[min_i] < fabs(DGM_dxy[min_j])/DGM_dxyError[min_j]) ? fabs(DGM_dxy[min_i])/DGM_dxyError[min_i] : fabs(DGM_dxy[min_j])/DGM_dxyError[min_j];

     dmdmCandidate.trackDxy_PV = (fabs(DGM_dxy_PV[min_i])/DGM_dxyError_PV[min_i] < fabs(DGM_dxy_PV[min_j])/DGM_dxyError_PV[min_j]) ? DGM_dxy_PV[min_i] : DGM_dxy_PV[min_j];
     dmdmCandidate.trackIxy_PV = (fabs(DGM_dxy_PV[min_i])/DGM_dxyError_PV[min_i] < fabs(DGM_dxy_PV[min_j])/DGM_dxyError_PV[min_j]) ? fabs(DGM_dxy_PV[min_i])/DGM_dxyError_PV[min_i] : fabs(DGM_dxy_PV[min_j])/DGM_dxyError_PV[min_j];

     dmdmCandidate.trackDxy_BS = (fabs(DGM_dxy_BS[min_i])/DGM_dxyError_BS[min_i] < fabs(DGM_dxy_BS[min_j])/DGM_dxyError_BS[min_j]) ? DGM_dxy_BS[min_i] : DGM_dxy_BS[min_j];
     dmdmCandidate.trackIxy_BS = (fabs(DGM_dxy_BS[min_i])/DGM_dxyError_BS[min_i] < fabs(DGM_dxy_BS[min_j])/DGM_dxyError_BS[min_j]) ? fabs(DGM_dxy_BS[min_i])/DGM_dxyError_BS[min_i] : fabs(DGM_dxy_BS[min_j])/DGM_dxyError_BS[min_j];


        DMDM_idxA[nDMDM] = min_i;
        DMDM_idxB[nDMDM] = min_j;
        DMDM_Lxy_PV[nDMDM] = dmdmCandidate.Lxy_PV;
        DMDM_Ixy_PV[nDMDM] = dmdmCandidate.Ixy_PV;
        DMDM_Lxy_0[nDMDM] = dmdmCandidate.Lxy_0;
        DMDM_Ixy_0[nDMDM] = dmdmCandidate.Ixy_0;
        DMDM_Lxy_BS[nDMDM] = dmdmCandidate.Lxy_BS;
        DMDM_Ixy_BS[nDMDM] = dmdmCandidate.Ixy_BS;
        DMDM_trackDxy[nDMDM] = dmdmCandidate.trackDxy;
        DMDM_trackIxy[nDMDM] = dmdmCandidate.trackIxy;
        DMDM_trackDxy_PV[nDMDM] = dmdmCandidate.trackDxy_PV;
        DMDM_trackIxy_PV[nDMDM] = dmdmCandidate.trackIxy_PV;
        DMDM_trackDxy_BS[nDMDM] = dmdmCandidate.trackDxy_BS;
        DMDM_trackIxy_BS[nDMDM] = dmdmCandidate.trackIxy_BS;
        DMDM_vx[nDMDM] = dmdmCandidate.vx;
        DMDM_vy[nDMDM] = dmdmCandidate.vy;
        DMDM_normalizedChi2[nDMDM] = dmdmCandidate.normalizedChi2;
        DMDM_mass[nDMDM] = dmdmCandidate.mass;
        DMDM_leadingPt[nDMDM] = dmdmCandidate.leadingPt;
        DMDM_subleadingPt[nDMDM] = dmdmCandidate.subleadingPt;
        DMDM_cosAlpha[nDMDM] = dmdmCandidate.cosAlpha;
        DMDM_dPhi[nDMDM] = dmdmCandidate.dPhi;
        DMDM_lldPhi[nDMDM] = dmdmCandidate.lldPhi;
        DMDM_dR[nDMDM] = dmdmCandidate.dR;
        DMDM_relisoA[nDMDM] = dmdmCandidate.relisoA;
        DMDM_relisoB[nDMDM] = dmdmCandidate.relisoB;

     nDMDM++;

   } // end while 


   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////// FILL THE TREE ///////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////

   bool filter = true;

   if (_filterByEE)
       filter = (nEE > 0) ? true : false;

   if (_filterByMM)
       filter = (nDMDM > 0) ? true : false;

   if (_filterByLL)
       filter = (nDMDM > 0 || nEE > 0) ? true : false;

   if (_filterByElectron)
       filter = (nElectronCandidate > 1) ? true : false;

   if (_filterByMuon)
       filter = (nDGM > 1) ? true : false;

   // Fill the tree
   if (filter)
     tree_out->Fill();

}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::beginJob()
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

    // Filter parameters
   _filterByEE       = parameters.getParameter<bool>("FilterByEE");
   _filterByMM       = parameters.getParameter<bool>("FilterByMM");
   _filterByLL       = parameters.getParameter<bool>("FilterByLL");
   _filterByElectron = parameters.getParameter<bool>("FilterByElectron");
   _filterByMuon     = parameters.getParameter<bool>("FilterByMuon");


    // PU reweighting
    if (_Era == 2016) {
      lumi_weights = edm::LumiReWeighting("2016MCPileupHistogram.root", "2016DataPileupHistogram.root", "pileup", "pileup");
    } else if (_Era == 2017) {
      lumi_weights = edm::LumiReWeighting("2017MCPileupHistogram.root", "2017DataPileupHistogram.root", "pileup", "pileup");
    } else if (_Era == 2018) {
      lumi_weights = edm::LumiReWeighting("2018MCPileupHistogram.root", "2018DataPileupHistogram.root", "pileup", "pileup");
    }


   // Load HLT paths for every year:
   //
   if (_Era == 2016) {

     HLTPaths_.push_back("HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10"); // signal
     HLTPaths_.push_back("HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10"); // backup
     HLTPaths_.push_back("HLT_DoubleMu28NoFiltersNoVtxDisplaced"); // backup
     HLTPaths_.push_back("HLT_DoubleMu23NoFiltersNoVtxDisplaced"); // backup
     HLTPaths_.push_back("HLT_DoubleMu33NoFiltersNoVtx"); // backup
     HLTPaths_.push_back("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55"); // backup
     HLTPaths_.push_back("HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55"); // backup
     HLTPaths_.push_back("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70"); // backup
     HLTPaths_.push_back("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"); // signal
     HLTPaths_.push_back("HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55"); // backup
     HLTPaths_.push_back("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55"); // signal
     HLTPaths_.push_back("HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15"); // signal
     HLTPaths_.push_back("HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60"); // backup (prescaled)
     HLTPaths_.push_back("HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15"); // backup (prescaled)
     HLTPaths_.push_back("HLT_DoublePhoton85"); // backup
     HLTPaths_.push_back("HLT_DoublePhoton60"); // backup
     HLTPaths_.push_back("HLT_PFMET120_PFMHT90_IDTight"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET120_PFMHT100_IDTight"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET120_PFMHT110_IDTight"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET120_PFMHT120_IDTight"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_MET200"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET170_HBHECleaned"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET300"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"); // MET path for trigger efficiency

   } else if (_Era == 2017) {
  
     HLTPaths_.push_back("HLT_DoubleMu43NoFiltersNoVtx"); // signal (?)
     HLTPaths_.push_back("HLT_DoubleMu48NoFiltersNoVtx"); // backup (?)
     HLTPaths_.push_back("HLT_DoublePhoton33_CaloIdL"); // prescaled
     HLTPaths_.push_back("HLT_DoublePhoton70"); // signal
     HLTPaths_.push_back("HLT_DoublePhoton85"); // backup
     HLTPaths_.push_back("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"); // signal
     HLTPaths_.push_back("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95"); // backup
     HLTPaths_.push_back("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55"); // signal
     HLTPaths_.push_back("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55"); // prescaled
     HLTPaths_.push_back("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55"); // prescaled
     HLTPaths_.push_back("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55"); // prescaled
     HLTPaths_.push_back("HLT_Photon60_R9Id90_CaloIdL_IsoL"); // prescaled
     HLTPaths_.push_back("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL"); // prescaled
     HLTPaths_.push_back("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15"); // no sensitive
     HLTPaths_.push_back("HLT_HISinglePhoton20_Eta3p1ForPPRef"); // prescaled
     HLTPaths_.push_back("HLT_HISinglePhoton30_Eta3p1ForPPRef"); // prescaled
     HLTPaths_.push_back("HLT_HISinglePhoton40_Eta3p1ForPPRef"); // prescaled
     HLTPaths_.push_back("HLT_HISinglePhoton50_Eta3p1ForPPRef"); // prescaled
     HLTPaths_.push_back("HLT_HISinglePhoton60_Eta3p1ForPPRef"); // prescaled
     HLTPaths_.push_back("HLT_Photon30_HoverELoose"); // prescaled
     HLTPaths_.push_back("HLT_Photon40_HoverELoose"); // prescaled
     HLTPaths_.push_back("HLT_Photon50_HoverELoose"); // prescaled
     HLTPaths_.push_back("HLT_Photon60_HoverELoose"); // prescaled
     HLTPaths_.push_back("HLT_PFMET120_PFMHT120_IDTight"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET120_PFMHT120_IDTight_PFHT60"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_CaloMET350_HBHECleaned"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET250_HBHECleaned"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"); // MET path for trigger efficiency

   } else if (_Era == 2018) {

     HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX3BX"); // backup
     HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX"); // backup
     HLTPaths_.push_back("HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4"); // backup
     HLTPaths_.push_back("HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4"); // backup
     HLTPaths_.push_back("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed"); // signal
     HLTPaths_.push_back("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed"); // backup
     HLTPaths_.push_back("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4"); // backup
     HLTPaths_.push_back("HLT_DoubleL2Mu23NoVtx_2Cha"); // signal
     HLTPaths_.push_back("HLT_DoubleL2Mu25NoVtx_2Cha"); // backup
     HLTPaths_.push_back("HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4"); // backup
     HLTPaths_.push_back("HLT_DoubleMu43NoFiltersNoVtx"); // backup (?)
     HLTPaths_.push_back("HLT_DoubleMu48NoFiltersNoVtx"); // backup (?)
     HLTPaths_.push_back("HLT_DoublePhoton33_CaloIdL");
     HLTPaths_.push_back("HLT_DoublePhoton70");
     HLTPaths_.push_back("HLT_DoublePhoton85");
     HLTPaths_.push_back("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90");
     HLTPaths_.push_back("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95");
     HLTPaths_.push_back("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto");
     HLTPaths_.push_back("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55");
     HLTPaths_.push_back("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55");
     HLTPaths_.push_back("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55");
     HLTPaths_.push_back("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55");
     HLTPaths_.push_back("HLT_Photon60_R9Id90_CaloIdL_IsoL");
     HLTPaths_.push_back("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL");
     HLTPaths_.push_back("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15");
     HLTPaths_.push_back("HLT_PFMET120_PFMHT120_IDTight"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET120_PFMHT120_IDTight_PFHT60"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_CaloMET350_HBHECleaned"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET200_HBHE_BeamHaloCleaned"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMET250_HBHECleaned"); // MET path for trigger efficiency
     HLTPaths_.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"); // MET path for trigger efficiency

   }




    ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Event_event", &Event_event, "Event_event/I");
    tree_out->Branch("Event_run", &Event_run, "Event_run/I");
    tree_out->Branch("Event_luminosityBlock", &Event_luminosityBlock, "Event_luminosityBlock/I");
    
    tree_out->Branch("nPU", &nPU, "nPU/I");
    tree_out->Branch("nPUTrue", &nPUTrue, "nPUTrue/I");
    tree_out->Branch("wPU", &wPU, "wPU/F");
    tree_out->Branch("genWeight", &genWeight, "genWeight/F");

    ///////////////////////////////// TRIGGER BRANCHES ///////////////////////////////

    for (unsigned int ihlt = 0; ihlt < HLTPaths_.size(); ihlt++) {
      tree_out->Branch(TString(HLTPaths_[ihlt]), &triggerPass[ihlt]);
    }

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
    
    if (_doCMSElectrons) {
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
  
    ///////////////////////////////// MUON BRANCHES /////////////////////////////////
    
    if (_doCMSMuons) {
       tree_out->Branch("nMuon", &nMuon, "nMuon/I");
       tree_out->Branch("MuonSel_pt", MuonSel_pt, "MuonSel_pt[nMuon]/F");
       tree_out->Branch("MuonSel_eta", MuonSel_eta, "MuonSel_eta[nMuon]/F");
       tree_out->Branch("MuonSel_phi", MuonSel_phi, "MuonSel_phi[nMuon]/F");
       tree_out->Branch("MuonSel_relIso", MuonSel_relIso, "MuonSel_relIso[nMuon]/F");
       tree_out->Branch("MuonSel_isGlobalMuon", MuonSel_isGlobalMuon, "MuonSel_isGlobalMuon[nMuon]/I");
       tree_out->Branch("MuonSel_isMediumMuon", MuonSel_isMediumMuon, "MuonSel_isMediumMuon[nMuon]/I");
       tree_out->Branch("MuonSel_isGoodMediumMuon", MuonSel_isGoodMediumMuon, "MuonSel_isGoodMediumMuon[nMuon]/I");
       tree_out->Branch("MuonSel_dB", MuonSel_dB, "MuonSel_dB[nMuon]/F");
       tree_out->Branch("MuonSel_edB", MuonSel_edB, "MuonSel_edB[nMuon]/F");
    }

    //////////////////////////// AOD MUON BRANCHES ///////////////////////////
    
    if (_DSAMode){

      tree_out->Branch("nDGM", &nDGM, "nDGM/I");
      tree_out->Branch("DGM_pt", DGM_pt, "DGM_pt[nDGM]/F");
      tree_out->Branch("DGM_ptError", DGM_ptError, "DGM_ptError[nDGM]/F");
      tree_out->Branch("DGM_eta", DGM_eta, "DGM_eta[nDGM]/F");
      tree_out->Branch("DGM_phi", DGM_phi, "DGM_phi[nDGM]/F");
      tree_out->Branch("DGM_dxy", DGM_dxy, "DGM_dxy[nDGM]/F");
      tree_out->Branch("DGM_dz", DGM_dz, "DGM_dz[nDGM]/F");
      tree_out->Branch("DGM_dxyError", DGM_dxyError, "DGM_dxyError[nDGM]/F");
      tree_out->Branch("DGM_dxy_PV", DGM_dxy_PV, "DGM_dxy_PV[nDGM]/F");
      tree_out->Branch("DGM_dxyError_PV", DGM_dxyError_PV, "DGM_dxyError_PV[nDGM]/F");
      tree_out->Branch("DGM_dxy_BS", DGM_dxy_BS, "DGM_dxy_BS[nDGM]/F");
      tree_out->Branch("DGM_dxyError_BS", DGM_dxyError_BS, "DGM_dxyError_BS[nDGM]/F");
      tree_out->Branch("DGM_relPFiso", DGM_relPFiso, "DGM_relPFiso[nDGM]/F");
      tree_out->Branch("DGM_relPFiso_noLT", DGM_relPFiso_noLT, "DGM_relPFiso_noLT[nDGM]/F");
      tree_out->Branch("DGM_normChi2", DGM_normChi2, "DGM_normChi2[nDGM]/F");
      tree_out->Branch("DGM_charge", DGM_charge, "DGM_charge[nDGM]/I");
      tree_out->Branch("DGM_numberOfValidHits", DGM_numberOfValidHits, "DGM_numberOfValidHits[nDGM]/I");
      tree_out->Branch("DGM_muonHits", DGM_muonHits, "DGM_muonHits[nDGM]/I");
      tree_out->Branch("DGM_outerTrackerHits", DGM_outerTrackerHits, "DGM_outerTrackerHits[nDGM]/I");
      tree_out->Branch("DGM_trackerHits", DGM_trackerHits, "DGM_trackerHits[nDGM]/I");
      tree_out->Branch("DGM_CSCHits", DGM_CSCHits, "DGM_CSCHits[nDGM]/I");
      tree_out->Branch("DGM_DTHits", DGM_DTHits, "DGM_DTHits[nDGM]/I");

      
    }


    //////////////////////////////// GENPARTICLE BRANCHES ///////////////////////////////

    tree_out->Branch("nGenLepton", &nGenLepton, "nGenLepton/I");
    tree_out->Branch("GenLeptonSel_pt", GenLeptonSel_pt, "GenLeptonSel_pt[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_E", GenLeptonSel_E, "GenLeptonSel_E[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_et", GenLeptonSel_et, "GenLeptonSel_et[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_eta", GenLeptonSel_eta, "GenLeptonSel_eta[nGenLepton]/F");
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

    
    tree_out->Branch("nHardProcessParticle", &nHardProcessParticle, "nHardProcessParticle/I");
    tree_out->Branch("HardProcessParticle_E", HardProcessParticle_E, "HardProcessParticle_E[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_pt", HardProcessParticle_pt, "HardProcessParticle_pt[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_eta", HardProcessParticle_eta, "HardProcessParticle_eta[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_phi", HardProcessParticle_phi, "HardProcessParticle_phi[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_vx", HardProcessParticle_vx, "HardProcessParticle_vx[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_vy", HardProcessParticle_vy, "HardProcessParticle_vy[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_vz", HardProcessParticle_vz, "HardProcessParticle_vz[nHardProcessParticle]/F");
    tree_out->Branch("HardProcessParticle_pdgId", HardProcessParticle_pdgId, "HardProcessParticle_pdgId[nHardProcessParticle]/I");

    //////////////////////////////////// MET BRANCHES ///////////////////////////////////

    tree_out->Branch("MET_pt", &MET_pt, "MET_pt/F");
    tree_out->Branch("MET_phi", &MET_phi, "MET_phi/F");

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
    tree_out->Branch("ElectronCandidate_pvAssociationQuality", ElectronCandidate_pvAssociationQuality, "ElectronCandidate_pvAssociationQuality[nElectronCandidate]/I");
    
    ///////////////////////////////// ACCEPTANCE CRITERIA //////////////////////////////

    tree_out->Branch("passAcceptanceCriteria", &passAcceptanceCriteria, "passAcceptanceCriteria/O");


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


    tree_out->Branch("nDMDM", &nDMDM, "nDMDM/I");
    if (_DSAMode) {
       tree_out->Branch("DMDM_idxA", DMDM_idxA, "DMDM_idxA[nDMDM]/I");
       tree_out->Branch("DMDM_idxB", DMDM_idxB, "DMDM_idxB[nDMDM]/I");
       tree_out->Branch("DMDM_Lxy_PV", DMDM_Lxy_PV, "DMDM_Lxy_PV[nDMDM]/F");
       tree_out->Branch("DMDM_Ixy_PV", DMDM_Ixy_PV, "DMDM_Ixy_PV[nDMDM]/F");
       tree_out->Branch("DMDM_Lxy_0", DMDM_Lxy_0, "DMDM_Lxy_0[nDMDM]/F");
       tree_out->Branch("DMDM_Ixy_0", DMDM_Ixy_0, "DMDM_Ixy_0[nDMDM]/F");
       tree_out->Branch("DMDM_Lxy_BS", DMDM_Lxy_BS, "DMDM_Lxy_BS[nDMDM]/F");
       tree_out->Branch("DMDM_Ixy_BS", DMDM_Ixy_BS, "DMDM_Ixy_BS[nDMDM]/F");
       tree_out->Branch("DMDM_trackDxy", DMDM_trackDxy, "DMDM_trackDxy[nDMDM]/F");
       tree_out->Branch("DMDM_trackIxy", DMDM_trackIxy, "DMDM_trackIxy[nDMDM]/F");
       tree_out->Branch("DMDM_trackDxy_PV", DMDM_trackDxy_PV, "DMDM_trackDxy_PV[nDMDM]/F");
       tree_out->Branch("DMDM_trackIxy_PV", DMDM_trackIxy_PV, "DMDM_trackIxy_PV[nDMDM]/F");
       tree_out->Branch("DMDM_trackDxy_BS", DMDM_trackDxy_BS, "DMDM_trackDxy_BS[nDMDM]/F");
       tree_out->Branch("DMDM_trackIxy_BS", DMDM_trackIxy_BS, "DMDM_trackIxy_BS[nDMDM]/F");
       tree_out->Branch("DMDM_vx", DMDM_vx, "DMDM_vx[nDMDM]/F");
       tree_out->Branch("DMDM_vy", DMDM_vy, "DMDM_vy[nDMDM]/F");
       tree_out->Branch("DMDM_mass", DMDM_mass, "DMDM_mass[nDMDM]/F");
       tree_out->Branch("DMDM_normalizedChi2", DMDM_normalizedChi2, "DMDM_normalizedChi2[nDMDM]/F");
       tree_out->Branch("DMDM_leadingPt", DMDM_leadingPt, "DMDM_leadingPt[nDMDM]/F");
       tree_out->Branch("DMDM_subleadingPt", DMDM_subleadingPt, "DMDM_subleadingPt[nDMDM]/F");
       tree_out->Branch("DMDM_cosAlpha", DMDM_cosAlpha, "DMDM_cosAlpha[nDMDM]/F");
       tree_out->Branch("DMDM_dPhi", DMDM_dPhi, "DMDM_dPhi[nDMDM]/F");
       tree_out->Branch("DMDM_dR", DMDM_dR, "DMDM_dR[nDMDM]/F");
       tree_out->Branch("DMDM_lldPhi", DMDM_lldPhi, "DMDM_lldPhi[nDMDM]/F");
       tree_out->Branch("DMDM_relisoA", DMDM_relisoA, "DMDM_relisoA[nDMDM]/F");
       tree_out->Branch("DMDM_relisoB", DMDM_relisoB, "DMDM_relisoB[nDMDM]/F");
    }

}
//=======================================================================================================================================================================================================================//

//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::endJob() 
{

    file_out->cd();
    tree_out->Write();
    counts->Write();
    sum2Weights->Write();
    file_out->Close();

}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//=======================================================================================================================================================================================================================//

std::string LongLivedAnalysis::getPathVersion(const edm:: TriggerNames &names, const std::string &rawPath) {

  std::string version;
  std::string finalPath;

  for (int s = 0; s < 20; s++) {

    version = rawPath + std::to_string(s);
    if (names.size() != names.triggerIndex(version)) {
      finalPath = version;
      break;
    }

  }

  return finalPath;

}



//=======================================================================================================================================================================================================================//

bool LongLivedAnalysis::passIsotrackSelection( const pat::IsolatedTrack &track) {

   // Quality cuts:
   const reco::HitPattern &hits = track.hitPattern();
   if (hits.numberOfValidTrackerHits() < 6) { return false; }
   //if (!track.isHighPurityTrack()) { return false;}

   // Isotrack must have packed candidate:
   const pat::PackedCandidateRef &pckCand = track.packedCandRef();
   if (!pckCand.isNonnull()) { return false; }
   if (!pckCand->hasTrackDetails()) { return false; }

   // Preselection cuts:
   if (track.pt() < 15) { return false; }
   if (fabs(track.eta()) > 2.4) { return false; }

   // To be noticed: Isolation cuts are applied later with the LLCandidate selection

   return true;
}


//=======================================================================================================================================================================================================================//

bool LongLivedAnalysis::passPhotonSelection( const pat::Photon &photon ) {

   // Quality cuts:
   if (photon.hadronicOverEm() > 0.05) { return false; } 
   if (photon.isEE() && photon.full5x5_sigmaIetaIeta() > 0.0425) { return false; }
   if (photon.isEB() && photon.full5x5_sigmaIetaIeta() > 0.0112) { return false; }

   // Preselection cuts:
   if (photon.et() < 15) {return false; }

   return true;

}

//=======================================================================================================================================================================================================================//


bool LongLivedAnalysis::passMuonSelection(const pat::Muon &muon) {

   if (muon.pt() < 20){ return false; }
   if (fabs(muon.eta()) > 2.4) { return false; }
   return true;
}


//=======================================================================================================================================================================================================================//

bool LongLivedAnalysis::passDGMSelection(const reco::Track &muon) {

   if (muon.pt() < 20){ return false; }
   if (fabs(muon.eta()) > 2.4) { return false; }
   if (muon.ptError()/muon.pt() > 0.3) {return false;}
   if (muon.hitPattern().numberOfValidMuonHits() < 13) {return false;}
   if (muon.hitPattern().numberOfValidStripHits() < 6) {return false;}
   //if (muon.numberOfValidHits() < 6) { return false; }
   return true;
}


//=======================================================================================================================================================================================================================//


float LongLivedAnalysis::computeDxy(const pat::IsolatedTrack & track, const reco::Vertex pv) {

   double vx = track.vx();
   double vy = track.vy();
   double phi = track.phi();
   double PVx = pv.x();
   double PVy = pv.y();

   double dxy = -(vx - PVx)*sin(phi) + (vy - PVy)*cos(phi);
   return dxy;
}


//=======================================================================================================================================================================================================================//


float LongLivedAnalysis::computeDxy(const reco::Track & track, const reco::Vertex pv) {

   double vx = track.vx();
   double vy = track.vy();
   double phi = track.phi();
   double PVx = pv.x();
   double PVy = pv.y();

   double dxy = -(vx - PVx)*sin(phi) + (vy - PVy)*cos(phi);
   return dxy;
}


//=======================================================================================================================================================================================================================//


float LongLivedAnalysis::computeDxyError(const pat::IsolatedTrack & track, const reco::Vertex pv) {

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


float LongLivedAnalysis::computeDxyError(const reco::Track & track, const reco::Vertex pv) {

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


float LongLivedAnalysis::computeRelIso(const reco::Track & track, edm::Handle<edm::View<pat::PackedCandidate> > pfs, edm::Handle<edm::View<pat::PackedCandidate> > lts, bool isPF, bool includeLostTracks) {


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

DEFINE_FWK_MODULE(LongLivedAnalysis);
