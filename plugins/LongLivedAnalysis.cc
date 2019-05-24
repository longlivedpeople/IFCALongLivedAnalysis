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
//                                                                                                                                                                                                                       //
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

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
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


//=======================================================================================================================================================================================================================//


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// FUNCTIONS ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

bool goodPhoton(const pat::Photon & photon)
{

    // Return true if the photon fulfills with the analysis requirements and false instead

    if (fabs(photon.eta()) > 1.4442 && fabs(photon.eta()) < 1.566) { return false; } // narrow EB region to be defined
    if (photon.hadronicOverEm() > 0.05) { return false; }
    if (photon.isEE() && photon.full5x5_sigmaIetaIeta() > 0.034) { return false; }
    if (photon.isEB() && photon.full5x5_sigmaIetaIeta() > 0.012) { return false; }

    return true;

}

bool goodTrack(const pat::IsolatedTrack & track)
{

    return true;

}


bool isGoodMuonTriggerObject( pat::TriggerObjectStandAlone obj)
{

    // Fill
    return true;

}


float dxy_value(const reco::GenParticle &p)
{

    float vx = p.vx();
    float vy = p.vy();
    float phi = p.phi();
  
    float dxy = -vx*sin(phi) + vy*cos(phi);
    return dxy;

}


bool isLongLivedLepton(const reco::GenParticle &p)
{

    if (!( abs(p.pdgId()) == 11 || abs(p.pdgId()) == 13)){ return false; }
    if (abs(p.mother()->pdgId()) != 1000022){ return false; }

    return true;

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
Float_t genWeight;

//-> TRIGGER TAGS
Bool_t Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6;
Bool_t Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8;

//-> PRIMARY VERTEX SELECTION
Int_t nPV;
Int_t nTruePV;
Float_t PV_vx;
Float_t PV_vy;
Float_t PV_vz;

// -> REFITTED PRIMARY VERTEX
Float_t RefittedPV_vx;
Float_t RefittedPV_vy;
Float_t RefittedPV_vz;
Int_t RefittedPV_nPFTrack;
Int_t RefittedPV_nLostTrack;


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
Float_t IsoTrackSel_phi[nIsoTrackMax];
Float_t IsoTrackSel_dxy[nIsoTrackMax];
Float_t IsoTrackSel_dxyError[nIsoTrackMax];
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
// Derived:
Float_t IsoTrackSel_dxySignificance[nIsoTrackMax];


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
Int_t PhotonSel_isGoodSC[nPhotonMax];

//-> ELECTRON SELECTION
const Int_t nElectronMax = 100;
Int_t nElectron;
Float_t ElectronSel_pt[nElectronMax];
Float_t ElectronSel_eta[nElectronMax];
Float_t ElectronSel_phi[nElectronMax];



//-> MUON TRIGGER OBJECT SELECTION
const Int_t nMuonTriggerObjectMax = 500;
Int_t nMuonTriggerObject;
Float_t MuonTriggerObjectSel_pt[nMuonTriggerObjectMax];
Float_t MuonTriggerObjectSel_eta[nMuonTriggerObjectMax];
Float_t MuonTriggerObjectSel_phi[nMuonTriggerObjectMax];

//-> GENLEPTON SELECTION
const Int_t nGenLeptonMax = 500;
Int_t nGenLepton;
Float_t GenLeptonSel_pt[nGenLeptonMax];
Float_t GenLeptonSel_et[nGenLeptonMax];
Float_t GenLeptonSel_eta[nGenLeptonMax];
Float_t GenLeptonSel_phi[nGenLeptonMax];
Int_t GenLeptonSel_pdgId[nGenLeptonMax];
Float_t GenLeptonSel_dxy[nGenLeptonMax];
Int_t GenLeptonSel_motherIdx[nGenLeptonMax];

//-> GENNEUTRALINO SELECTION
const Int_t nGenNeutralinoMax = 500;
Int_t nGenNeutralino;
Float_t GenNeutralinoSel_pt[nGenNeutralinoMax];
Float_t GenNeutralinoSel_eta[nGenNeutralinoMax];
Float_t GenNeutralinoSel_phi[nGenNeutralinoMax];
Int_t GenNeutralinoSel_pdgId[nGenNeutralinoMax];
Float_t GenNeutralinoSel_Lxy[nGenNeutralinoMax];

// -> GENERATION ACCEPTANCE CRITERIA
Bool_t passAcceptanceCriteria;


// -> MET 

Int_t MET_isPFMET;
Float_t MET_pt; // default Type 1 corrected MET
Float_t MET_phi;
Float_t MET_sumEt;
Float_t MET_genPt;
Float_t MET_genPhi;
Float_t MET_corPt;
Float_t MET_corPhi;
Float_t MET_uncorPt;
Float_t MET_uncorPhi;
Float_t MET_metSignificance;
Float_t MET_NeutralEMFraction;
Float_t MET_NeutralHadEtFraction;
Float_t MET_ChargedEMEtFraction;
Float_t MET_ChargedHadEtFraction;
Float_t MET_MuonEtFraction;
Float_t MET_Type6EtFraction;
Float_t MET_Type7EtFraction;

//-> ELECTRON CANDIDATE SELECTION
const Int_t nElectronCandidateMax = 100;
Int_t nElectronCandidate;
Float_t ElectronCandidate_pt[nElectronCandidateMax];
Float_t ElectronCandidate_et[nElectronCandidateMax];
Float_t ElectronCandidate_eta[nElectronCandidateMax];
Float_t ElectronCandidate_phi[nElectronCandidateMax];
Int_t ElectronCandidate_photonIdx[nElectronCandidateMax];
Int_t ElectronCandidate_isotrackIdx[nElectronCandidateMax];

//-> MUON CANDIDATE SELECTION
const Int_t nMuonCandidateMax = 100;
Int_t nMuonCandidate;
Float_t MuonCandidate_pt[nMuonCandidateMax];
Float_t MuonCandidate_eta[nMuonCandidateMax];
Float_t MuonCandidate_phi[nMuonCandidateMax];
Float_t MuonCandidate_triggerPt[nMuonCandidateMax];
Int_t MuonCandidate_muonTriggerObjectIdx[nMuonCandidateMax];
Int_t MuonCandidate_isotrackIdx[nMuonCandidateMax];

//-> LL Candidates                                                                                                                                                                                 
const Int_t nLLMax = 1000;

Int_t nLL;
Float_t LL_Lxy[nLLMax];
Float_t LL_Ixy[nLLMax];
Float_t LL_Mass[nLLMax];

Float_t LLSel_Lxy;
Float_t LLSel_Ixy;
Float_t LLSel_Mass;
Int_t LLSel_isMM;
Int_t LLSel_isEE;


Int_t nEE;
Float_t EE_Lxy[nLLMax];
Float_t EE_Ixy[nLLMax];
Float_t EE_Mass[nLLMax];

Int_t nMM;
Float_t MM_Lxy[nLLMax];
Float_t MM_Ixy[nLLMax];
Float_t MM_Mass[nLLMax];



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

      std::string output_filename;
      edm::ParameterSet parameters;

      edm::EDGetTokenT<edm::View<pat::Electron> > theElectronCollection;   
      edm::EDGetTokenT<edm::View<pat::Photon> > thePhotonCollection;
      edm::EDGetTokenT<edm::View<pat::IsolatedTrack> >  theIsoTrackCollection;
      edm::EDGetTokenT<edm::View<reco::Vertex> > thePrimaryVertexCollection;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > thePackedPFCandidateCollection;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > theLostTracksCollection;
      edm::EDGetTokenT<edm::View<pat::MET> > theMETCollection;

      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone> > triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales>  triggerPrescales_;

      edm::EDGetTokenT<reco::BeamSpot> theBeamSpot;

      // Gen collection
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  theGenParticleCollection;
      
      edm::EDGetTokenT<GenEventInfoProduct>  theGenEventInfoProduct;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  thePileUpSummary;

      


      //"Global" variables
      std::vector<int> iT; // track indexes
      edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;

      // Class functions
      bool buildLLcandidate(edm::Handle<edm::View<pat::IsolatedTrack> > const& isotracks, int idxA, int idxB, bool isEE);

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
   thePhotonCollection = consumes<edm::View<pat::Photon> > (parameters.getParameter<edm::InputTag>("PhotonCollection"));
   theIsoTrackCollection = consumes<edm::View<pat::IsolatedTrack> >  (parameters.getParameter<edm::InputTag>("IsoTrackCollection"));
   thePrimaryVertexCollection = consumes<edm::View<reco::Vertex> >  (parameters.getParameter<edm::InputTag>("PrimaryVertexCollection"));
   thePackedPFCandidateCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("PackedPFCandidateCollection"));
   theLostTracksCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("LostTracksCollection"));
   theMETCollection = consumes<edm::View<pat::MET> >  (parameters.getParameter<edm::InputTag>("METCollection"));


   triggerBits_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("bits"));
   triggerObjects_ = consumes<edm::View<pat::TriggerObjectStandAlone> > (parameters.getParameter<edm::InputTag>("objects"));
   triggerPrescales_ = consumes<pat::PackedTriggerPrescales > (parameters.getParameter<edm::InputTag>("prescales"));

   theBeamSpot = consumes<reco::BeamSpot>  (parameters.getParameter<edm::InputTag>("BeamSpot"));


   theGenParticleCollection = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genParticleCollection"));

   theGenEventInfoProduct = consumes<GenEventInfoProduct> (parameters.getParameter<edm::InputTag>("theGenEventInfoProduct"));

   thePileUpSummary = consumes<std::vector<PileupSummaryInfo> > (parameters.getParameter<edm::InputTag>("thePileUpSummary"));

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
   

   edm::Handle<edm::View<pat::Electron> > electrons;
   edm::Handle<edm::View<pat::Photon> > photons;
   edm::Handle<edm::View<pat::IsolatedTrack> > isotracks;
   edm::Handle<edm::View<reco::Vertex> > primaryvertices;
   edm::Handle<edm::View<pat::PackedCandidate> > packedPFCandidates;
   edm::Handle<edm::View<pat::PackedCandidate> > lostTracks;
   edm::Handle<edm::View<pat::MET> > METs;

   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<edm::View<pat::TriggerObjectStandAlone>  >triggerObjects;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

   edm::Handle<reco::BeamSpot> beamSpot;

   edm::Handle<edm::View<reco::GenParticle> > genParticles;

   edm::Handle<GenEventInfoProduct> genEvtInfo;

   edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;

   iEvent.getByToken(theElectronCollection, electrons);
   iEvent.getByToken(thePhotonCollection, photons);
   iEvent.getByToken(theIsoTrackCollection, isotracks);
   iEvent.getByToken(thePrimaryVertexCollection, primaryvertices);
   iEvent.getByToken(thePackedPFCandidateCollection, packedPFCandidates);
   iEvent.getByToken(theLostTracksCollection, lostTracks);
   iEvent.getByToken(theMETCollection, METs);

   iEvent.getByToken(triggerBits_, triggerBits);
   iEvent.getByToken(triggerObjects_, triggerObjects);
   //iEvent.getByToken(triggerPrescales_, triggerPrescales);

   iEvent.getByToken(theBeamSpot, beamSpot);

   iEvent.getByToken(theGenParticleCollection, genParticles);
    
   iEvent.getByToken(theGenEventInfoProduct, genEvtInfo);

   iEvent.getByToken(thePileUpSummary, puInfoH);



   /////////////////////////////////// MC OR DATA //////////////////////////////////////
   bool isMC = true;

   /////////////////////////////////// EVENT INFO //////////////////////////////////////

   
   Event_event = iEvent.id().event();
   Event_run = iEvent.id().run();
   Event_luminosityBlock = iEvent.id().luminosityBlock();
   
   genWeight = (float) genEvtInfo->weight();
   counts->Fill(0.5);
   sum2Weights->Fill(0.5, genWeight*genWeight);

   //////////////////////////////////// PILE UP ////////////////////////////////////////

   for(size_t i=0;i<puInfoH->size();++i) {
       if( puInfoH->at(i).getBunchCrossing() == 0) {
            nPU = puInfoH->at(i).getPU_NumInteractions();
            nPUTrue = puInfoH->at(i).getTrueNumInteractions();                                          
       }
   }
   //////////////////////////////////// BEAM SPOT //////////////////////////////////////

   reco::BeamSpot beamSpotObject = *beamSpot;
   BeamSpot_x0 = beamSpotObject.x0();
   BeamSpot_y0 = beamSpotObject.y0();
   BeamSpot_z0 = beamSpotObject.z0();
   BeamSpot_BeamWidthX = beamSpotObject.BeamWidthX();
   BeamSpot_BeamWidthY = beamSpotObject.BeamWidthY();


   /////////////////////////////// TRIGGER ACCEPTANCE //////////////////////////////////

   const std::string muonTriggerName = "HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6";
   //"HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v7"; // default 
   const std::string photonTriggerName = "HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8";
   //"HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v9";

   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6 = triggerBits->accept(names.triggerIndex(photonTriggerName));
   Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8 = triggerBits->accept(names.triggerIndex(muonTriggerName));


   ////////////////////////////// MUON TRIGGER OBJECTS /////////////////////////////////
 
   std::vector<int> iMT; // muon trigger object indexes 


   // Loop to get Muon Trigger objects:
   for (size_t i = 0; i < triggerObjects->size(); i++) 
   {
       

       pat::TriggerObjectStandAlone obj = (*triggerObjects)[i];


       obj.unpackPathNames(names);
       obj.unpackFilterLabels(iEvent, *triggerBits);    


       bool isMuonTriggerObject = obj.hasPathName( "HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6", true, true );

       if (!isMuonTriggerObject) { continue; }

       if (isGoodMuonTriggerObject(obj)) iMT.push_back(i);

   }
   


   // Sort the muon trigger objects by pt
   std::sort( std::begin(iMT), std::end(iMT), [&](int i1, int i2){ return triggerObjects->at(i1).pt() < triggerObjects->at(i2).pt(); });



   // Fill the muon trigger objects features
   nMuonTriggerObject = iMT.size();

   for (size_t i = 0; i < iMT.size(); i++){
   
       pat::TriggerObjectStandAlone obj = (*triggerObjects)[iMT.at(i)];

       obj.unpackPathNames(names);
       obj.unpackFilterLabels(iEvent, *triggerBits);

       MuonTriggerObjectSel_pt[i] = obj.pt();
       MuonTriggerObjectSel_eta[i] = obj.eta();
       MuonTriggerObjectSel_phi[i] = obj.phi();
   

   }


   /////  -> Code from twiki here: Print trigger info:
    
   //const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   /* 
   std::cout << "\n == TRIGGER PATHS= " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        std::cout << "Trigger " << names.triggerName(i) <<
                ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
                << std::endl;
    }

    */
    
   /* 

    std::cout << "\n TRIGGER OBJECTS " << std::endl;

    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(names);
        //obj.unpackFilterLabels(names);
        std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;

 
    // Print trigger object collection and type
    std::cout << "\t   Collection: " << obj.collection() << std::endl;
    std::cout << "\t   Type IDs:   ";
    for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
    std::cout << std::endl;
    // Print associated trigger filters
    obj.unpackFilterLabels(iEvent, *triggerBits); // added
    std::cout << "\t   Filters:    ";
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
    std::cout << std::endl;
    std::vector< std::string > pathNamesAll = obj.pathNames(false);
    std::vector< std::string > pathNamesLast = obj.pathNames(true);

    // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
    // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
    // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)

    std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
        bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
        bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
        bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
        bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
        std::cout << "   " << pathNamesAll[h];
        if (isBoth) std::cout << "(L,3)";
        if (isL3 && !isBoth) std::cout << "(*,3)";
        if (isLF && !isBoth) std::cout << "(L,*)";
        if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
    }
    std::cout << std::endl;
    }
    std::cout << std::endl;

   */


   ////////////////////////////// PRIMARY VERTEX FEATURES //////////////////////////////

   nTruePV = 0;
   nPV = primaryvertices->size();
   //const reco::Vertex &thevertex = (*primaryvertices)[0];
   for (size_t i = 0; i < primaryvertices->size(); i ++){

       const reco::Vertex &current_vertex = (*primaryvertices)[i];
       if(current_vertex.isValid()){ nTruePV++; }

   }


   ///////////////////////////////// ISOTRACK FEATURES /////////////////////////////////
   iT.clear();
   
   for (size_t i = 0; i < isotracks->size(); i++){

       const pat::IsolatedTrack & isotrack = (*isotracks)[i];
       if (goodTrack(isotrack)){ iT.push_back(i); }

   }


   nIsoTrack = iT.size(); // number of isotracks

   // Sort the isotracks by pt
   std::sort( std::begin(iT), std::end(iT), [&](int i1, int i2){ return isotracks->at(i1).pt() < isotracks->at(i2).pt(); });


   // Loop over the isotracks
   for (size_t i = 0; i < iT.size(); ++i){

       
       const pat::IsolatedTrack & isotrack = (*isotracks)[iT.at(i)];

       // Basic features:
       IsoTrackSel_pt[i] = isotrack.pt();
       IsoTrackSel_eta[i] = isotrack.eta();
       IsoTrackSel_phi[i] = isotrack.phi();
       
       // Isolation info:
       
       const pat::PFIsolation &pfiso = isotrack.pfIsolationDR03();
       const pat::PFIsolation &minipfiso = isotrack.miniPFIsolation();

       IsoTrackSel_pfIsolationDR03[i] = pfiso.chargedHadronIso() + pfiso.neutralHadronIso() + pfiso.photonIso() + pfiso.puChargedHadronIso();
       IsoTrackSel_miniPFIsolation[i] = minipfiso.chargedHadronIso() + minipfiso.neutralHadronIso() + minipfiso.photonIso() + minipfiso.puChargedHadronIso();
       IsoTrackSel_relPfIsolationDR03[i] = IsoTrackSel_pfIsolationDR03[i]/isotrack.pt();
       IsoTrackSel_relMiniPFIsolation[i] = IsoTrackSel_miniPFIsolation[i]/isotrack.pt();

       // Quality info:
       IsoTrackSel_isHighPurityTrack[i] = isotrack.isHighPurityTrack();

       // Impact parameter info:
       IsoTrackSel_dxy[i] = isotrack.dxy();
       IsoTrackSel_dxyError[i] = isotrack.dxyError();
       IsoTrackSel_dz[i] = isotrack.dz();
       IsoTrackSel_dzError[i] = isotrack.dzError();
       IsoTrackSel_dxySignificance[i] = fabs(isotrack.dxy())/isotrack.dxyError();


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
       // (PV(), vertex())
       IsoTrackSel_fromPV[i] = isotrack.fromPV(); 
       const pat::PackedCandidateRef &pckCand = isotrack.packedCandRef(); // access the packed candidate
      
 
       if (isotrack.fromPV() > -1){ // check it has a PV

           IsoTrackSel_vx[i] = (*pckCand).vx();
           IsoTrackSel_vy[i] = (*pckCand).vy();
           IsoTrackSel_vz[i] = (*pckCand).vz();

           const reco::VertexRef &PV = (*pckCand).vertexRef(); // access the PV of the candidate
           IsoTrackSel_PVx[i] = (*PV).x();
           IsoTrackSel_PVy[i] = (*PV).y();
           IsoTrackSel_PVz[i] = (*PV).z();



       }else{

           IsoTrackSel_vx[i] = -99;
           IsoTrackSel_vy[i] = -99;
           IsoTrackSel_vz[i] = -99;

           IsoTrackSel_PVx[i] = -99;
           IsoTrackSel_PVy[i] = -99;
           IsoTrackSel_PVz[i] = -99;

       }


   }





   ////////////////////////////////// PHOTON FEATURES //////////////////////////////////
   
   std::vector<int> iP; // photon indexes


   // Select good photons
   for (size_t i = 0; i < photons->size(); i++){
       
       // this is the place to put any preselection if required
       iP.push_back(i);

   }

   // Sort good lepton indexes by pt
   std::sort( std::begin(iP), std::end(iP), [&](int i1, int i2){ return photons->at(i1).et() < photons->at(i2).et(); });


   nPhoton = iP.size();
   // Loop over the photons
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

       PhotonSel_hcalIso[i] = photon.hcalIso();
       PhotonSel_ecalIso[i] = photon.ecalIso();
       PhotonSel_caloIso[i] = photon.caloIso();
       PhotonSel_relIso[i] = photon.caloIso()/photon.et();

       if (goodPhoton(photon)){PhotonSel_isGoodSC[i] = 1; }
       else {PhotonSel_isGoodSC[i] = 0; }

   }


   ///////////////////////////////// ELECTRON FEATURES /////////////////////////////////

   nElectron = electrons->size();

   for (size_t i = 0; i < electrons->size(); i++){

       const pat::Electron & electron = (* electrons)[i];

       ElectronSel_pt[i] = electron.pt();
       ElectronSel_eta[i] = electron.eta();
       ElectronSel_phi[i] = electron.phi();

   }


   //////////////////////////////// GENPARTICLE FEATURES ///////////////////////////////

   std::vector<int> iGL; // Generated lepton indexes
   std::vector<int> iGN; // Generated neutralino indexes

   for(size_t i = 0; i < genParticles->size(); i++) {


        const reco::GenParticle &genparticle = (*genParticles)[i];

        // Check if it is a lepton comming from a long lived neutralino:
        if (isLongLivedLepton(genparticle)){ iGL.push_back(i); continue; }

        // Check if it is a longlived neutralino (correct id and after all radiative emission)
        if (abs(genparticle.pdgId()) == 1000022 && (abs(genparticle.daughter(0)->pdgId()) != 1000022) && (abs(genparticle.daughter(0)->pdgId()) != 22)){
            iGN.push_back(i);
            continue;

        }
   } 


   // Number of generated leptons:
   nGenLepton = iGL.size();

   // Number of generated neutralinos:
   nGenNeutralino = iGN.size();

   // Sort the leptons and neutralinos by pt
   std::sort( std::begin(iGL), std::end(iGL), [&](int i1, int i2){ return genParticles->at(i1).pt() < genParticles->at(i2).pt(); });
   std::sort( std::begin(iGN), std::end(iGN), [&](int i1, int i2){ return genParticles->at(i1).pt() < genParticles->at(i2).pt(); });


   // Loop over the selected neutralinos
   for(size_t i = 0; i < iGN.size(); i++){

       const reco::GenParticle &genparticle = (*genParticles)[iGN.at(i)];

       GenNeutralinoSel_pt[i] = genparticle.pt();
       GenNeutralinoSel_eta[i] = genparticle.eta();
       GenNeutralinoSel_phi[i] = genparticle.phi();
       GenNeutralinoSel_pdgId[i] = genparticle.pdgId();

       // Generated transverse decay length of the neutralino

       const reco::Candidate *m = genparticle.mother();
   
       // To avoid radiative effects of the neutralino: 
       while((abs(m->pdgId()) == 100022)){ m = m->mother(); }

       GenNeutralinoSel_Lxy[i] = sqrt((m->vx()-genparticle.daughter(0)->vx())*(m->vx()-genparticle.daughter(0)->vx()) + (m->vy()-genparticle.daughter(0)->vy())*(m->vy()-genparticle.daughter(0)->vy()));

   }



   // Loop over the selected genleptons
   for(size_t i = 0; i < iGL.size(); i++){

       const reco::GenParticle &genparticle = (*genParticles)[iGL.at(i)];

       GenLeptonSel_pdgId[i] = genparticle.pdgId();
       GenLeptonSel_dxy[i] = dxy_value(genparticle);


       // Define the mother index
       GenLeptonSel_motherIdx[i] = -99;
       if(genparticle.mother()->pt() == GenNeutralinoSel_pt[0]){

           GenLeptonSel_motherIdx[i] = 0;

       } else if (genparticle.mother()->pt() == GenNeutralinoSel_pt[1]){

           GenLeptonSel_motherIdx[i] = 1;

       }
       

       // Get the last genlepton (to avoid radiative effects):
       if (genparticle.numberOfDaughters() > 0){

           const reco::Candidate *d = genparticle.daughter(0);
           while(d->numberOfDaughters()> 0 && d->daughter(0)->pdgId() == d->pdgId()){ d = d->daughter(0); }

           GenLeptonSel_pt[i] = d->pt();
           GenLeptonSel_et[i] = d->et();
           GenLeptonSel_eta[i] = d->eta();
           GenLeptonSel_phi[i] = d->phi();

       } else {

           GenLeptonSel_pt[i] = genparticle.pt();
           GenLeptonSel_et[i] = genparticle.et();
           GenLeptonSel_eta[i] = genparticle.eta();
           GenLeptonSel_phi[i] = genparticle.phi();

       }


   }



   //////////////////////////////// ACCEPTANCE CRITERIA ////////////////////////////////

   passAcceptanceCriteria = true; // true by default
   bool passLeadingElectron = false;

   for (size_t i = 0; i < iGL.size(); i++){

       if (abs(GenLeptonSel_pdgId[i]) == 11 && GenLeptonSel_et[i] > 40){ passLeadingElectron = true;}

       if (abs(GenLeptonSel_pdgId[i]) == 11 && GenLeptonSel_et[i] < 25){ passAcceptanceCriteria = false; break;}
       if (abs(GenLeptonSel_pdgId[i]) == 13 && GenLeptonSel_pt[i] < 26){ passAcceptanceCriteria = false; break;}
       if (fabs(GenLeptonSel_eta[i]) > 2){ passAcceptanceCriteria = false; break;}

   }

   if (passLeadingElectron == false) {passAcceptanceCriteria = false;}
   
   if (GenNeutralinoSel_Lxy[0] > 50 || GenNeutralinoSel_Lxy[1] > 50){ passAcceptanceCriteria = false;}


   //////////////////////////////////////// MET ////////////////////////////////////////

   const pat::MET &met = (*METs)[0]; // access the MET object

   MET_pt = met.pt();
   MET_phi = met.phi();
   MET_sumEt = met.sumEt();
   MET_corPt = met.corPt();
   MET_corPhi = met.corPhi();
   MET_uncorPt = met.uncorPt();
   MET_uncorPhi = met.uncorPhi();
   MET_metSignificance = met.metSignificance();

   if (isMC) {

        MET_genPt = met.genMET()->pt(); 
        MET_genPhi = met.genMET()->phi();

   } else { 

        MET_genPt = -99;
        MET_genPhi = -99;

   }


   // Specific PFMET stuff (Supposed to be filled always):
   if (met.isPFMET()){

       MET_NeutralEMFraction = met.NeutralEMFraction();
       MET_NeutralHadEtFraction = met.NeutralHadEtFraction();
       MET_ChargedEMEtFraction = met.ChargedEMEtFraction();
       MET_ChargedHadEtFraction = met.ChargedHadEtFraction();
       MET_MuonEtFraction = met.MuonEtFraction();
       MET_Type6EtFraction = met.Type6EtFraction();
       MET_Type7EtFraction = met.Type7EtFraction();

   } else {

       MET_NeutralEMFraction = -99;
       MET_NeutralHadEtFraction = -99;
       MET_ChargedEMEtFraction = -99;
       MET_ChargedHadEtFraction = -99;
       MET_MuonEtFraction = -99;
       MET_Type6EtFraction = -99;
       MET_Type7EtFraction = -99;

   }


   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////// LEPTON CANDIDATES /////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////


   // Variable initiallization:

   std::vector<int> matched_tracks, matched_SC, matched_triggerObjects; // std vectors with matched objects to avoid overlapping

   float dRMin = 99999; // dR to minimize as high as possible in the beginning
   float dRThreshold = 0.1; // Maximum dR to do the lepton matching
   float dR; // Computation of dR
   int matching_type = 0; // 0 if electron; 1 if muon
   int tmin, scmin, tomin, li; // minimum track, minimum SC, minimum trigger object, reconstructed lepton index


   // while loop that looks for the minimum matchings and stops when the dRThreshold is reached:

   while (1){

       dRMin = 99999; // redefinicion
       matching_type = 99;

       // Loop over the tracks
       for (size_t t = 0; t < iT.size(); t++){

           const pat::IsolatedTrack & isotrack = (*isotracks)[iT.at(t)];

           // pass if the track is associated already:
           if(std::find(matched_tracks.begin(), matched_tracks.end(), t) != matched_tracks.end()){ continue; }
           // pass if the track does not fulfil the prerequisites:
           if(!goodTrack(isotrack)){ continue; }


           // Loop over the superclusters:
           for (size_t sc = 0; sc < iP.size(); sc++){

               const pat::Photon & photon = (*photons)[iP.at(sc)];

               // pass if the SC is associated to other track:
               if(std::find(matched_SC.begin(), matched_SC.end(), sc) != matched_SC.end()){ continue; }

               // pass if the SC does not fulfil the prerequisites:
               if (!goodPhoton(photon)){ continue; }
              
               // ------------ SC matching -------------
               dR = getDeltaR(isotrack.phi(), isotrack.eta(), photon.phi(), photon.eta());
                      
               if (dR < dRMin){

                   dRMin = dR;
                   matching_type = 0;
                   scmin = sc;
                   tmin = t;

               }
           }


           // Loop over the muon trigger objecs
           for (size_t to = 0; to < iMT.size(); to++){

               pat::TriggerObjectStandAlone muon = (*triggerObjects)[iMT.at(to)];

               // pass if the trigger muon is associated to other track
               if(std::find(matched_triggerObjects.begin(), matched_triggerObjects.end(), to) != matched_triggerObjects.end()){ continue; }

               // ------------ trigger Object matching --------------
               dR = getDeltaR(isotrack.phi(), isotrack.eta(), muon.phi(), muon.eta());

               if (dR < dRMin){

                   /*
                   std::cout<< "new: " << "\t"<< isotrack.phi() << "\t" << isotrack.eta() << "\t" << muon.phi() << "\t" << muon.eta() << std::endl;
                   std::cout << dR << std::endl;
                   */

                   dRMin = dR;
                   matching_type = 1;
                   tomin = to; // careful cause it is not a muon trigger object of the ntuples
                   tmin = t;

               }
           }
       }

       

       if (dRMin > dRThreshold){ break; } // Here we go out the while loop

       if (matching_type == 0){

           li = matched_SC.size();
           ElectronCandidate_pt[li] = (*isotracks)[tmin].pt();
           ElectronCandidate_eta[li] = (*isotracks)[tmin].eta();
           ElectronCandidate_phi[li] = (*isotracks)[tmin].phi();
           ElectronCandidate_et[li] = (*photons)[scmin].et();
           ElectronCandidate_photonIdx[li] = scmin;
           ElectronCandidate_isotrackIdx[li] = tmin;

           matched_SC.push_back(scmin); matched_tracks.push_back(tmin);


       } else if (matching_type == 1){

           li = matched_triggerObjects.size();
           
           MuonCandidate_pt[li] = (*isotracks)[tmin].pt();
           MuonCandidate_eta[li] = (*isotracks)[tmin].eta();
           MuonCandidate_phi[li] = (*isotracks)[tmin].phi();
           MuonCandidate_triggerPt[li] = (*triggerObjects)[tomin].pt();
           MuonCandidate_muonTriggerObjectIdx[li] = tomin; // not well defined
           MuonCandidate_isotrackIdx[li] = tmin;

           matched_triggerObjects.push_back(tomin); matched_tracks.push_back(tmin);

       }
         
   }

   nElectronCandidate = matched_SC.size();
   nMuonCandidate = matched_triggerObjects.size();

 
   /////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////// VERTEX REFITTING /////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////

 
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);

   std::vector<reco::TransientTrack> refit_tracks; // tracks for refitting

   RefittedPV_nPFTrack = 0;
   RefittedPV_nLostTrack = 0;

   // Original PV values
   PV_vx = -99;
   PV_vy = -99;
   PV_vz = -99;

   // No suceed PV refitting default values
   RefittedPV_vx = -99;
   RefittedPV_vy = -99;
   RefittedPV_vz = -99;

 
   // Loop over packedPFCandidates
   for (size_t i = 0; i < packedPFCandidates->size(); i++ ){

       const pat::PackedCandidate &packedPFCandidate = (*packedPFCandidates)[i];

       if (!packedPFCandidate.hasTrackDetails()) continue;
       const reco::Track packedPFTrack = packedPFCandidate.pseudoTrack();

       // Selection criteria for the tracks
       if (packedPFCandidate.fromPV() != 3) continue;


       if(PV_vx == -99 && PV_vy == -99 && PV_vz == -99){ 

           const reco::VertexRef &PV = packedPFCandidate.vertexRef(); // access the PV of the candidate
           PV_vx = (*PV).x(); PV_vy = (*PV).y(); PV_vz = (*PV).z();      
 
       }

       RefittedPV_nPFTrack++;
       reco::TransientTrack  transientTrack = theTransientTrackBuilder->build(packedPFTrack);
       transientTrack.setBeamSpot(beamSpotObject);
       refit_tracks.push_back(transientTrack);


   }

   
   // Loop over lostTracks
   for (size_t i = 0; i < lostTracks->size(); i++){
       
       const pat::PackedCandidate &lostTrack = (*lostTracks)[i];

       if (!lostTrack.hasTrackDetails()) continue;
       const reco::Track packedLostTrack = lostTrack.pseudoTrack();

       if (lostTrack.fromPV() != 3) continue;

       RefittedPV_nLostTrack++;
       reco::TransientTrack  transientTrack = theTransientTrackBuilder->build(packedLostTrack);
       transientTrack.setBeamSpot(beamSpotObject);
       refit_tracks.push_back(transientTrack);



   }




   // Reffit the vertex

   if (refit_tracks.size() > 1){
       
       AdaptiveVertexFitter  theFitter(GeometricAnnealing(2.5));
       TransientVertex myVertex = theFitter.vertex(refit_tracks);

       if (myVertex.isValid()){

           RefittedPV_vx = myVertex.position().x();
           RefittedPV_vy = myVertex.position().y();
           RefittedPV_vz = myVertex.position().z();

       }

   }


   ////////////////////////////////// LL CANDIDATES //////////////////////////////////                                                    
   nLL = 0;
   nEE = 0;
   nMM = 0;

   LLSel_Lxy = -99;
   LLSel_Ixy = -99;
   LLSel_Mass = -99;
   LLSel_isEE = 0;
   LLSel_isMM = 0;
   for (int i = 1; i < nElectronCandidate; i++) {
     for (int j = 0; j < i; j++) {
       if (i != j) {
	 //if (i == 1 && j == 0) {//OOOOOOOOOOOOOOOOJOOOOOOOOOOOOOOOOOOO

	 //std::cout << "checking pair (" << i << ", " << j << ")" << std::endl;
	 bool goodpair = buildLLcandidate(isotracks, i, j, true);
         if (goodpair) {;}// do nothing, just avoid warning
	 //if (goodpair) std::cout << "valid candidate" << std::endl;
	 //else std::cout << "bad candidate" << std::endl;
       }
     }
   }

   for (int i = 1; i < nMuonCandidate; i++) {
     for (int j = 0; j < i; j++) {
       if (i != j) {
	 //if (i == 1 && j == 0) {//OOOOOOOOOOOOOOOOJOOOOOOOOOOOOOOOOOOO

	 //std::cout << "checking pair (" << i << ", " << j << ")" << std::endl;
	 bool goodpair = buildLLcandidate(isotracks, i, j, false);
         if (goodpair) {;}// do nothing, just avoid warning
	 //if (goodpair) std::cout << "valid candidate" << std::endl;
	 //else std::cout << "bad candidate" << std::endl;
       }
     }
   }


   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////// FILL THE TREE ///////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////
   tree_out->Fill();


}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::beginJob()
{

    // Output file definition
    output_filename = parameters.getParameter<std::string>("nameOfOutput");
    file_out = new TFile(output_filename.c_str(), "RECREATE");
    file_out->cd();

    
    // Output Tree definition
    tree_out = new TTree("Events", "Events");


    ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Event_event", &Event_event, "Event_event/I");
    tree_out->Branch("Event_run", &Event_run, "Event_run/I");
    tree_out->Branch("Event_luminosityBlock", &Event_luminosityBlock, "Event_luminosityBlock/I");
    
    tree_out->Branch("nPU", &nPU, "nPU/I");
    tree_out->Branch("nPUTrue", &nPUTrue, "nPUTrue/I");
    tree_out->Branch("genWeight", &genWeight, "genWeight/F");

    ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6", &Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6, "Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6/O");
    tree_out->Branch("Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8", &Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8, "Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8/O");

    ///////////////////////////////// BEAM SPOT BRANCHES ////////////////////////////////

    tree_out->Branch("BeamSpot_x0", &BeamSpot_x0, "BeamSpot_x0/F");
    tree_out->Branch("BeamSpot_y0", &BeamSpot_y0, "BeamSpot_y0/F");
    tree_out->Branch("BeamSpot_z0", &BeamSpot_z0, "BeamSpot_z0/F");
    tree_out->Branch("BeamSpot_BeamWidthX", &BeamSpot_BeamWidthX, "BeamSpot_BeamWidthX/F");
    tree_out->Branch("BeamSpot_BeamWidthY", &BeamSpot_BeamWidthY, "BeamSpot_BeamWidthY/F");


    ////////////////////////////// PRIMARY VERTEX BRANCHES //////////////////////////////

    tree_out->Branch("nPV", &nPV, "nPV/I");
    tree_out->Branch("nTruePV", &nTruePV, "nTruePV/I");
    tree_out->Branch("PV_vx", &PV_vx, "PV_vx/F");
    tree_out->Branch("PV_vy", &PV_vy, "PV_vy/F");
    tree_out->Branch("PV_vz", &PV_vz, "PV_vz/F");


    /////////////////////////// REFITTED PRIMARY VERTEX BRANCHES ////////////////////////

    tree_out->Branch("RefittedPV_vx", &RefittedPV_vx, "RefittedPV_vx/F");
    tree_out->Branch("RefittedPV_vy", &RefittedPV_vy, "RefittedPV_vy/F");
    tree_out->Branch("RefittedPV_vz", &RefittedPV_vz, "RefittedPV_vz/F");
    tree_out->Branch("RefittedPV_nPFTrack", &RefittedPV_nPFTrack, "RefittedPV_nPFTrack/I");
    tree_out->Branch("RefittedPV_nLostTrack", &RefittedPV_nLostTrack, "RefittedPV_nLostTrack/I");


    ///////////////////////////////// ISOTRACK BRANCHES /////////////////////////////////
    
    tree_out->Branch("nIsoTrack", &nIsoTrack, "nIsoTrack/I");
    tree_out->Branch("IsoTrackSel_pt", IsoTrackSel_pt, "IsoTrackSel_pt[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_eta", IsoTrackSel_eta, "IsoTrackSel_eta[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_phi", IsoTrackSel_phi, "IsoTrackSel_phi[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxy", IsoTrackSel_dxy, "IsoTrackSel_dxy[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxyError", IsoTrackSel_dxyError, "IsoTrackSel_dxyError[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dz", IsoTrackSel_dz, "IsoTrackSel_dz[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dzError", IsoTrackSel_dzError, "IsoTrackSel_dzError[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxySignificance", IsoTrackSel_dxySignificance, "IsoTrackSel_dxySignificance[nIsoTrack]/F");
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
    tree_out->Branch("IsoTrackSel_fromPV", IsoTrackSel_fromPV, "IsoTrackSel_fromPV[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_PVx", IsoTrackSel_PVx, "IsoTrackSel_PVx[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_PVy", IsoTrackSel_PVy, "IsoTrackSel_PVy[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_PVz", IsoTrackSel_PVz, "IsoTrackSel_PVz[nIsoTrack]/F");


    
    ////////////////////////////////// PHOTON BRANCHES //////////////////////////////////

    tree_out->Branch("nPhoton", &nPhoton, "nPhoton/I");
    tree_out->Branch("PhotonSel_et", PhotonSel_et, "PhotonSel_et[nPhoton]/F");
    tree_out->Branch("PhotonSel_eta", PhotonSel_eta, "PhotonSel_eta[nPhoton]/F");
    tree_out->Branch("PhotonSel_phi", PhotonSel_phi, "PhotonSel_phi[nPhoton]/F");
    tree_out->Branch("PhotonSel_hadronicOverEm", PhotonSel_hadronicOverEm, "PhotonSel_hadronicOverEm[nPhoton]/F");
    tree_out->Branch("PhotonSel_full5x5_sigmaIetaIeta", PhotonSel_full5x5_sigmaIetaIeta, "PhotonSel_full5x5_sigmaIetaIeta[nPhoton]/F");
    tree_out->Branch("PhotonSel_isEB", PhotonSel_isEB, "PhotonSel_isEB[nPhoton]/I");
    tree_out->Branch("PhotonSel_isEE", PhotonSel_isEE, "PhotonSel_isEE[nPhoton]/I");
    tree_out->Branch("PhotonSel_r9", PhotonSel_r9, "PhotonSel_r9[nPhoton]/F");
    tree_out->Branch("PhotonSel_ecalIso", PhotonSel_ecalIso, "PhotonSel_ecalIso[nPhoton]/F");
    tree_out->Branch("PhotonSel_hcalIso", PhotonSel_hcalIso, "PhotonSel_hcalIso[nPhoton]/F");
    tree_out->Branch("PhotonSel_caloIso", PhotonSel_caloIso, "PhotonSel_caloIso[nPhoton]/F");
    tree_out->Branch("PhotonSel_relIso", PhotonSel_relIso, "PhotonSel_relIso[nPhoton]/F");
    tree_out->Branch("PhotonSel_isGoodSC", PhotonSel_isGoodSC, "PhotonSel_isGoodSC[nPhoton]/I");


    ///////////////////////////////// ELECTRON BRANCHES /////////////////////////////////

    tree_out->Branch("nElectron", &nElectron, "nElectron/I");
    tree_out->Branch("ElectronSel_pt", ElectronSel_pt, "ElectronSel_pt[nElectron]/F");
    tree_out->Branch("ElectronSel_eta", ElectronSel_eta, "ElectronSel_eta[nElectron]/F");
    tree_out->Branch("ElectronSel_phi", ElectronSel_phi, "ElectronSel_phi[nElectron]/F");



    //////////////////////////// MUON TRIGGER OBJECT BRANCHES ///////////////////////////
    //
    tree_out->Branch("nMuonTriggerObject", &nMuonTriggerObject, "nMuonTriggerObject/I");
    tree_out->Branch("MuonTriggerObjectSel_pt", MuonTriggerObjectSel_pt, "MuonTriggerObjectSel_pt[nMuonTriggerObject]/F");
    tree_out->Branch("MuonTriggerObjectSel_eta", MuonTriggerObjectSel_eta, "MuonTriggerObjectSel_eta[nMuonTriggerObject]/F");
    tree_out->Branch("MuonTriggerObjectSel_phi", MuonTriggerObjectSel_phi, "MuonTriggerObjectSel_phi[nMuonTriggerObject]/F");


    //////////////////////////////// GENPARTICLE BRANCHES ///////////////////////////////

    tree_out->Branch("nGenLepton", &nGenLepton, "nGenLepton/I");
    tree_out->Branch("GenLeptonSel_pt", GenLeptonSel_pt, "GenLeptonSel_pt[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_et", GenLeptonSel_et, "GenLeptonSel_et[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_eta", GenLeptonSel_eta, "GenLeptonSel_eta[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_phi", GenLeptonSel_phi, "GenLeptonSel_phi[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_dxy", GenLeptonSel_dxy, "GenLeptonSel_dxy[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_pdgId", GenLeptonSel_pdgId, "GenLeptonSel_pdgId[nGenLepton]/I");
    tree_out->Branch("GenLeptonSel_motherIdx", GenLeptonSel_motherIdx, "GenLeptonSel_motherIdx[nGenLepton]/I");


    tree_out->Branch("nGenNeutralino", &nGenNeutralino, "nGenNeutralino/I");
    tree_out->Branch("GenNeutralinoSel_pt", GenNeutralinoSel_pt, "GenNeutralinoSel_pt[nGenNeutralino]/F");
    tree_out->Branch("GenNeutralinoSel_eta", GenNeutralinoSel_eta, "GenNeutralinoSel_eta[nGenNeutralino]/F");
    tree_out->Branch("GenNeutralinoSel_phi", GenNeutralinoSel_phi, "GenNeutralinoSel_phi[nGenNeutralino]/F");
    tree_out->Branch("GenNeutralinoSel_Lxy", GenNeutralinoSel_Lxy, "GenNeutralinoSel_Lxy[nGenNeutralino]/F");
    tree_out->Branch("GenNeutralinoSel_pdgId", GenNeutralinoSel_pdgId, "GenNeutralinoSel_pdgId[nGenNeutralino]/I");

    //////////////////////////////////// MET BRANCHES ///////////////////////////////////

    tree_out->Branch("MET_pt", &MET_pt, "MET_pt/F");
    tree_out->Branch("MET_phi", &MET_phi, "MET_phi/F");
    tree_out->Branch("MET_sumEt", &MET_sumEt, "MET_sumEt/F");
    tree_out->Branch("MET_genPt", &MET_genPt, "MET_genPt/F");
    tree_out->Branch("MET_genPhi", &MET_genPhi, "MET_genPhi/F");
    tree_out->Branch("MET_corPt", &MET_corPt, "MET_corPt/F");
    tree_out->Branch("MET_corPhi", &MET_corPhi, "MET_corPhi/F");
    tree_out->Branch("MET_uncorPt", &MET_uncorPt, "MET_uncorPt/F");
    tree_out->Branch("MET_uncorPhi", &MET_uncorPhi, "MET_uncorPhi/F");
    tree_out->Branch("MET_metSignificance", &MET_metSignificance, "MET_metSignificance/F");
    tree_out->Branch("MET_NeutralEMFraction", &MET_NeutralEMFraction, "MET_NeutralEMFraction/F");
    tree_out->Branch("MET_NeutralHadEtFraction", &MET_NeutralHadEtFraction, "MET_NeutralHadEtFraction/F");
    tree_out->Branch("MET_ChargedEMEtFraction", &MET_ChargedEMEtFraction, "MET_ChargedEMEtFraction/F");
    tree_out->Branch("MET_ChargedHadEtFraction", &MET_ChargedHadEtFraction, "MET_ChargedHadEtFraction/F");
    tree_out->Branch("MET_MuonEtFraction", &MET_MuonEtFraction, "MET_MuonEtFraction/F");
    tree_out->Branch("MET_Type6EtFraction", &MET_Type6EtFraction, "MET_Type6EtFraction/F");
    tree_out->Branch("MET_Type7EtFraction", &MET_Type7EtFraction, "MET_Type7EtFraction/F");


    //////////////////////////// ELECTRON CANDIDATE BRANCHES ////////////////////////////

    tree_out->Branch("nElectronCandidate", &nElectronCandidate, "nElectronCandidate/I");
    tree_out->Branch("ElectronCandidate_pt", ElectronCandidate_pt, "ElectronCandidate_pt[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_et", ElectronCandidate_et, "ElectronCandidate_et[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_eta", ElectronCandidate_eta, "ElectronCandidate_eta[nElectronCandidate]/F");    
    tree_out->Branch("ElectronCandidate_phi", ElectronCandidate_phi, "ElectronCandidate_phi[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_photonIdx", ElectronCandidate_photonIdx, "ElectronCandidate_photonIdx[nElectronCandidate]/I");
    tree_out->Branch("ElectronCandidate_isotrackIdx", ElectronCandidate_isotrackIdx, "ElectronCandidate_isotrackIdx[nElectronCandidate]/I");


    ///////////////////////////////// ACCEPTANCE CRITERIA //////////////////////////////

    tree_out->Branch("passAcceptanceCriteria", &passAcceptanceCriteria, "passAcceptanceCriteria/O");


    ////////////////////////////// MUON CANDIDATE BRANCHES /////////////////////////////

    tree_out->Branch("nMuonCandidate", &nMuonCandidate, "nMuonCandidate/I");
    tree_out->Branch("MuonCandidate_pt", MuonCandidate_pt, "MuonCandidate_pt[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_eta", MuonCandidate_eta, "MuonCandidate_eta[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_phi", MuonCandidate_phi, "MuonCandidate_phi[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_triggerPt", MuonCandidate_triggerPt, "MuonCandidate_triggerPt[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_muonTriggerObjectIdx", MuonCandidate_muonTriggerObjectIdx, "MuonCandidate_muonTriggerObjectIdx[nMuonCandidate]/I");
    tree_out->Branch("MuonCandidate_isotrackIdx", MuonCandidate_isotrackIdx, "MuonCandidate_isotrackIdx[nMuonCandidate]/I");

    ////////////////////////////// LL BRANCHES /////////////////////////////
    tree_out->Branch("nLL", &nLL, "nLL/I");
    tree_out->Branch("LL_Lxy", LL_Lxy, "LL_Lxy[nLL]/F");
    tree_out->Branch("LL_Ixy", LL_Ixy, "LL_Ixy[nLL]/F");
    tree_out->Branch("LL_Mass", LL_Mass, "LL_Mass[nLL]/F");
    tree_out->Branch("nEE", &nEE, "nEE/I");
    tree_out->Branch("EE_Lxy", EE_Lxy, "EE_Lxy[nEE]/F");
    tree_out->Branch("EE_Ixy", EE_Ixy, "EE_Ixy[nEE]/F");
    tree_out->Branch("EE_Mass", EE_Mass, "EE_Mass[nEE]/F");
    tree_out->Branch("nMM", &nMM, "nMM/I");
    tree_out->Branch("MM_Lxy", MM_Lxy, "MM_Lxy[nMM]/F");
    tree_out->Branch("MM_Ixy", MM_Ixy, "MM_Ixy[nMM]/F");
    tree_out->Branch("MM_Mass", MM_Mass, "MM_Mass[nMM]/F");
    tree_out->Branch("LLSel_Lxy", &LLSel_Lxy, "LLSel_Lxy/F");
    tree_out->Branch("LLSel_Ixy", &LLSel_Ixy, "LLSel_Ixy/F");
    tree_out->Branch("LLSel_Mass", &LLSel_Mass, "LLSel_Mass/F");
    tree_out->Branch("LLSel_isEE", &LLSel_isEE, "LLSel_isEE/I");
    tree_out->Branch("LLSel_isMM", &LLSel_isMM, "LLSel_isMM/I");
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


bool LongLivedAnalysis::buildLLcandidate(edm::Handle<edm::View<pat::IsolatedTrack> > const& isotracks, int idxA, int idxB, bool isEE) {

   std::vector<reco::TransientTrack> vec_refitTracks;
   std::vector<reco::Track> vec_refitRecoTracks;
  const pat::IsolatedTrack & it_A = (isEE)?(*isotracks)[iT.at(ElectronCandidate_isotrackIdx[idxA])]: (*isotracks)[iT.at(MuonCandidate_isotrackIdx[idxA])];
  const pat::IsolatedTrack & it_B = (isEE)?(*isotracks)[iT.at(ElectronCandidate_isotrackIdx[idxB])]: (*isotracks)[iT.at(MuonCandidate_isotrackIdx[idxB])];


  const pat::PackedCandidateRef &pckCandA = it_A.packedCandRef();
  const pat::PackedCandidateRef &pckCandB = it_B.packedCandRef();


  if (pckCandA.isNonnull() && pckCandA->hasTrackDetails() && pckCandB.isNonnull() && pckCandB->hasTrackDetails()) { 
    const reco::Track isorecotrkA = pckCandA->pseudoTrack();
    reco::TransientTrack isotransienttrackA = theTransientTrackBuilder->build(isorecotrkA);

    const reco::Track isorecotrkB = pckCandB->pseudoTrack();
    reco::TransientTrack isotransienttrackB = theTransientTrackBuilder->build(isorecotrkB);

    //std::cout << isorecotrkA.pt() << "<-|isotrack_A before SV isotrack_Bb|->" << isorecotrkB.pt() << std::endl;
    vec_refitTracks.push_back(isotransienttrackA);
    vec_refitTracks.push_back(isotransienttrackB);

    vec_refitRecoTracks.push_back(isorecotrkA);
    vec_refitRecoTracks.push_back(isorecotrkB);
  }

  if (vec_refitTracks.size() > 1) {
    AdaptiveVertexFitter  thefitterll(GeometricAnnealing(2.5));

    // Vertex refitting:                                                                                                                                                     
    TransientVertex myVertex = thefitterll.vertex(vec_refitTracks);
    if (myVertex.isValid()) {
      const reco::VertexRef &PV = (*pckCandA).vertexRef();
      const reco::Vertex* recopv = PV.get();
      const reco::Vertex pv = *recopv;

      GlobalVector axis(0,0,0); 
       
      const reco::Vertex secV = myVertex;

      axis = GlobalVector(secV.x(),secV.y(),secV.z());


      Measurement1D vMeas = reco::SecondaryVertex::computeDist2d(pv,secV,axis,true);



      reco::TrackKinematics secVkin(vec_refitRecoTracks);
      //std::cout << "secVkin.nTrack =" << secVkin.numberOfTracks() << std::endl;
      //std::cout << "secVkin.nTrack =" << secVkin.weightedVectorSum().M() << std::endl;

      LL_Lxy[nLL] = vMeas.value();
      LL_Ixy[nLL] = vMeas.significance();
      LL_Mass[nLL] = secVkin.weightedVectorSum().M(); 

      //we update the value if it has not been initialized or if it has been and the new LL is more displaced
      if ( (!LLSel_isEE && !LLSel_isMM) ||
	   ( (LLSel_isEE || LLSel_isMM) && (fabs(LL_Ixy[nLL]) < fabs(LLSel_Ixy)) ) )
	{
	  if (isEE) {
	    LLSel_isEE = true;
	    LLSel_isMM = false;
	  }
	  else {
	    LLSel_isMM = true;
	    LLSel_isEE = false;
	  }

	  LLSel_Lxy = vMeas.value();
	  LLSel_Ixy = vMeas.significance();
	  LLSel_Mass = secVkin.weightedVectorSum().M(); 
      }

      nLL++;

      if (isEE) {
	EE_Lxy[nEE] = vMeas.value();
	EE_Ixy[nEE] = vMeas.significance();
	EE_Mass[nEE] = secVkin.weightedVectorSum().M(); 
	nEE++;
      }
      else {
	MM_Lxy[nMM] = vMeas.value();
	MM_Ixy[nMM] = vMeas.significance();
	MM_Mass[nMM] = secVkin.weightedVectorSum().M(); 
	nMM++;
      }
    }
    else return false;
  } 
  else return false;

  return true;
  
}




DEFINE_FWK_MODULE(LongLivedAnalysis);
