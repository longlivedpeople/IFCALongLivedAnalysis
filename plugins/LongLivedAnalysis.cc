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


#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"


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



#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"


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


/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// DATA DEFINITION //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////// BRANCHES /////////////////////////////////////

//-> EVENT INFO
Int_t Event_event;
Int_t Event_luminosityBlock;
Int_t Event_run;


//-> PRIMARY VERTEX SELECTION
Int_t nPV;
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
Int_t IsoTrackSel_isHighPurityTrack[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidTrackerHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidPixelHits[nIsoTrackMax];
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

//-> GENPARTICLE SELECTION
const Int_t nGenParticleMax = 500;
Int_t nGenParticle;
Float_t GenParticleSel_pt[nGenParticleMax];
Float_t GenParticleSel_eta[nGenParticleMax];
Float_t GenParticleSel_phi[nGenParticleMax];
Int_t GenParticleSel_pdgId[nGenParticleMax];
Float_t GenParticleSel_dxy[nGenParticleMax];

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
Int_t MuonCandidate_muonTriggerObjectIdx[nMuonCandidateMax];
Int_t MuonCandidate_isotrackIdx[nMuonCandidateMax];




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

      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone> > triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales>  triggerPrescales_;

      edm::EDGetTokenT<reco::BeamSpot> theBeamSpot;

      // Gen collection
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  theGenParticleCollection;


};
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
LongLivedAnalysis::LongLivedAnalysis(const edm::ParameterSet& iConfig)
{
   usesResource("TFileService");
   
   parameters = iConfig;


   theElectronCollection = consumes<edm::View<pat::Electron> >  (parameters.getParameter<edm::InputTag>("ElectronCollection"));
   thePhotonCollection = consumes<edm::View<pat::Photon> > (parameters.getParameter<edm::InputTag>("PhotonCollection"));
   theIsoTrackCollection = consumes<edm::View<pat::IsolatedTrack> >  (parameters.getParameter<edm::InputTag>("IsoTrackCollection"));
   thePrimaryVertexCollection = consumes<edm::View<reco::Vertex> >  (parameters.getParameter<edm::InputTag>("PrimaryVertexCollection"));
   thePackedPFCandidateCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("PackedPFCandidateCollection"));
   theLostTracksCollection = consumes<edm::View<pat::PackedCandidate> >  (parameters.getParameter<edm::InputTag>("LostTracksCollection"));

   triggerBits_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("bits"));
   triggerObjects_ = consumes<edm::View<pat::TriggerObjectStandAlone> > (parameters.getParameter<edm::InputTag>("objects"));
   triggerPrescales_ = consumes<pat::PackedTriggerPrescales > (parameters.getParameter<edm::InputTag>("prescales"));

   theBeamSpot = consumes<reco::BeamSpot>  (parameters.getParameter<edm::InputTag>("BeamSpot"));


   theGenParticleCollection = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genParticleCollection"));


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


   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<edm::View<pat::TriggerObjectStandAlone>  >triggerObjects;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

   edm::Handle<reco::BeamSpot> beamSpot;

   edm::Handle<edm::View<reco::GenParticle> > genParticles;




   iEvent.getByToken(theElectronCollection, electrons);
   iEvent.getByToken(thePhotonCollection, photons);
   iEvent.getByToken(theIsoTrackCollection, isotracks);
   iEvent.getByToken(thePrimaryVertexCollection, primaryvertices);
   iEvent.getByToken(thePackedPFCandidateCollection, packedPFCandidates);
   iEvent.getByToken(theLostTracksCollection, lostTracks);


   iEvent.getByToken(triggerBits_, triggerBits);
   iEvent.getByToken(triggerObjects_, triggerObjects);
   iEvent.getByToken(triggerPrescales_, triggerPrescales);

   iEvent.getByToken(theBeamSpot, beamSpot);

   iEvent.getByToken(theGenParticleCollection, genParticles);


   /////////////////////////////////// EVENT INFO //////////////////////////////////////

   
   Event_event = iEvent.id().event();
   Event_run = iEvent.id().run();
   Event_luminosityBlock = iEvent.id().luminosityBlock();
   

   std::cout << iEvent.id().luminosityBlock() << std::endl;

   //////////////////////////////////// BEAM SPOT //////////////////////////////////////

   reco::BeamSpot beamSpotObject = *beamSpot;
   BeamSpot_x0 = beamSpotObject.x0();
   BeamSpot_y0 = beamSpotObject.y0();
   BeamSpot_z0 = beamSpotObject.z0();
   BeamSpot_BeamWidthX = beamSpotObject.BeamWidthX();
   BeamSpot_BeamWidthY = beamSpotObject.BeamWidthY();


   ////////////////////////////// MUON TRIGGER OBJECTS /////////////////////////////////
 
   std::vector<int> iMT; // muon trigger object indexes 
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   std::string muonTriggerName = "HLT_DoubleMu43NoFiltersNoVtx_v3"; // default 

   // Loop to get Muon Trigger objects:
   for (size_t i = 0; i < triggerObjects->size(); i++) 
   {
       

       pat::TriggerObjectStandAlone obj = (*triggerObjects)[i];


       obj.unpackPathNames(names);
       obj.unpackFilterLabels(iEvent, *triggerBits);    


       bool isMuonTriggerObject = obj.hasPathName( muonTriggerName, true, true );

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
   /* 
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   std::cout << "\n == TRIGGER PATHS= " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        std::cout << "Trigger " << names.triggerName(i) <<
                ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
                ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
                << std::endl;
    }


    

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


   nPV = primaryvertices->size();
   //const reco::Vertex &thevertex = (*primaryvertices)[0];



   ///////////////////////////////// ISOTRACK FEATURES /////////////////////////////////

   std::vector<int> iT; // track indexes

   for (size_t i = 0; i < isotracks->size(); i++){

       const pat::IsolatedTrack & isotrack = (*isotracks)[i];
       if (goodTrack(isotrack)){ iT.push_back(i); }

   }


   nIsoTrack = iT.size(); // number of isotracks


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


       // Info extracted form the packedCandidate of the isotrack
       // (PV(), vertex())
       IsoTrackSel_fromPV[i] = isotrack.fromPV(); 
       const pat::PackedCandidateRef &pckCand = isotrack.packedCandRef(); // access the packed candidate
      
 
       if (isotrack.fromPV() > 0){ // check it has a PV

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
   for (size_t i = 0; i < photons->size(); ++i){
       
       const pat::Photon & photon = (*photons)[i];

       if (goodPhoton(photon)) 
       { 
           iP.push_back(i); 
       }

   }

   // Sort good lepton indexes by pt
   std::sort( std::begin(iP), std::end(iP), [&](int i1, int i2){ return photons->at(i1).et() < photons->at(i2).et(); });


   nPhoton = iP.size();
   // Loop over the good photons
   for (size_t i = 0; i < iP.size(); ++i){

       const pat::Photon & photon = (*photons)[iP.at(i)];

       PhotonSel_et[i] = photon.et();
       PhotonSel_eta[i] = photon.eta();
       PhotonSel_phi[i] = photon.phi();
       PhotonSel_hadronicOverEm[i] = photon.hadronicOverEm();
       PhotonSel_full5x5_sigmaIetaIeta[i] = photon.full5x5_sigmaIetaIeta();
       PhotonSel_isEB[i] = photon.isEB();
       PhotonSel_isEE[i] = photon.isEE();


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

   std::vector<int> iGP;

   for(size_t i = 0; i < genParticles->size(); i++) {


        const reco::GenParticle &genparticle = (*genParticles)[i];

        if (isLongLivedLepton(genparticle)){ iGP.push_back(i); }


   } 



   nGenParticle = iGP.size();
   // Loop over the selected gen particles
   for(size_t i = 0; i < iGP.size(); i++){

       const reco::GenParticle &genparticle = (*genParticles)[iGP.at(i)];

       GenParticleSel_pdgId[i] = genparticle.pdgId();
       GenParticleSel_dxy[i] = dxy_value(genparticle);
       

       // Get the last genparticle (to avoid radiative effects):
       if (genparticle.numberOfDaughters() > 0){

           const reco::Candidate *d = genparticle.daughter(0);
           while(d->numberOfDaughters()> 0 && d->daughter(0)->pdgId() == d->pdgId()){ d = d->daughter(0); }

           GenParticleSel_pt[i] = d->pt();
           GenParticleSel_eta[i] = d->eta();
           GenParticleSel_phi[i] = d->phi();

       } else {

           GenParticleSel_pt[i] = genparticle.pt();
           GenParticleSel_eta[i] = genparticle.eta();
           GenParticleSel_phi[i] = genparticle.phi();

       }


   }



   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////// LEPTON CANDIDATES /////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////


   int e = 0; // index of the electron candidate
   int m = 0; // index of the muon candidate
   float deltaR_min; // variable ot found the minimum deltaR

   std::vector<int> matched_clusters; // free clusters to match

   std::vector<int> matched_triggerObjects; // free trigger objects to match

   int m_cluster;
   int m_triggerObject;


   // Loop over the isolated tracks to do a lepton matching
   for (size_t i = 0; i < iT.size(); ++i){

       const pat::IsolatedTrack & isotrack = (*isotracks)[iT.at(i)];
       
       // Matching variables initiallization:
       deltaR_min = 10.; 
       m_cluster = -99;
       m_triggerObject = -99;

       // Electron matching
       for (size_t jp = 0; jp < iP.size(); ++jp){


           const pat::Photon & photon = (*photons)[iP.at(jp)];

           float deltaPhi = fabs(photon.phi() - isotrack.phi());
           float deltaEta = fabs(photon.eta() - isotrack.eta());
           float deltaR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

           if (deltaR < 0.1 && deltaR < deltaR_min){

               m_cluster = jp;
               deltaR_min = deltaR;

           }           

       }


       // Muon matching
       for (size_t jm = 0; jm < iMT.size(); ++jm){


           pat::TriggerObjectStandAlone obj = (*triggerObjects)[iMT.at(jm)];

           obj.unpackPathNames(names);
           obj.unpackFilterLabels(iEvent, *triggerBits);

           float deltaPhi = fabs(obj.phi() - isotrack.phi());
           float deltaEta = fabs(obj.eta() - isotrack.eta());
           float deltaR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

           if (deltaR < 0.1 && deltaR < deltaR_min){

               m_triggerObject = jm; m_cluster = -99;
               deltaR_min = deltaR;

           }

       }


       // Final lepton matching
       if (m_cluster == -99 && m_triggerObject == -99){ continue; // no matching
       }
       else if(m_cluster == -99 && m_triggerObject != -99){ // muon candidate found
       
            pat::TriggerObjectStandAlone obj = (*triggerObjects)[iMT.at(m_triggerObject)];

            obj.unpackPathNames(names);
            obj.unpackFilterLabels(iEvent, *triggerBits);
  
            if(std::find(matched_triggerObjects.begin(), matched_triggerObjects.end(), m_triggerObject) != matched_triggerObjects.end()){ continue; }
 
            MuonCandidate_pt[m] = isotrack.pt();
            MuonCandidate_phi[m] = isotrack.phi();
            MuonCandidate_eta[m] = isotrack.eta();
            MuonCandidate_muonTriggerObjectIdx[m] = m_triggerObject;
            MuonCandidate_isotrackIdx[m] = i;
            matched_triggerObjects.push_back(m_triggerObject);

            m++; // Next muon candidate


       }
       else if(m_cluster != -99 && m_triggerObject == -99){ // electron candidate found
     
           if(std::find(matched_clusters.begin(), matched_clusters.end(), m_cluster) != matched_clusters.end()){ continue; }

           const pat::Photon & photon = (*photons)[iP.at(m_cluster)];

           ElectronCandidate_pt[e] = isotrack.pt();
           ElectronCandidate_et[e] = photon.et();
           ElectronCandidate_phi[e] = isotrack.phi();
           ElectronCandidate_eta[e] = isotrack.eta();
           ElectronCandidate_photonIdx[e] = m_cluster;
           ElectronCandidate_isotrackIdx[e] = i;
           matched_clusters.push_back(m_cluster);
               
           e++; // Next electron candidate

       }


   }
 
   nElectronCandidate = e; // number of electron candidates = last idx filled + 1
   nMuonCandidate = m; // number of muon candidates = last idx filled + 1


 
 
   /////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////// VERTEX REFITTING /////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////

 
   edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
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
       //transientTrack.setBeamSpot(beamSpotObject);
       refit_tracks.push_back(transientTrack);


   }

   std::cout << "PFTRack: " << RefittedPV_nPFTrack << "\t" << packedPFCandidates->size() << std::endl;

   
   // Loop over lostTracks
   for (size_t i = 0; i < lostTracks->size(); i++){
       
       const pat::PackedCandidate &lostTrack = (*lostTracks)[i];

       if (!lostTrack.hasTrackDetails()) continue;
       const reco::Track packedLostTrack = lostTrack.pseudoTrack();

       if (lostTrack.fromPV() !=  3) continue;

       RefittedPV_nLostTrack++;

       reco::TransientTrack  transientTrack = theTransientTrackBuilder->build(packedLostTrack);
       //transientTrack.setBeamSpot(beamSpotObject);
       refit_tracks.push_back(transientTrack);



   }


   // Reffit the vertex

   if (refit_tracks.size() > 1){
       AdaptiveVertexFitter  theFitter;
       TransientVertex myVertex = theFitter.vertex(refit_tracks);

       if (myVertex.isValid()){

           RefittedPV_vx = myVertex.position().x();
           RefittedPV_vy = myVertex.position().y();
           RefittedPV_vz = myVertex.position().z();

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

    std::cout << "the file is created" << std::endl;
    
    // Output Tree definition
    tree_out = new TTree("Events", "Events");

    std::cout << "The tree is created" << std::endl;

    ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Event_event", &Event_event, "Event_event/I");
    tree_out->Branch("Event_run", &Event_run, "Event_run/I");
    tree_out->Branch("Event_luminosityBlock", &Event_luminosityBlock, "Event_luminosityBlock/I");


    ///////////////////////////////// BEAM SPOT BRANCHES ////////////////////////////////

    tree_out->Branch("BeamSpot_x0", &BeamSpot_x0, "BeamSpot_x0/F");
    tree_out->Branch("BeamSpot_y0", &BeamSpot_y0, "BeamSpot_y0/F");
    tree_out->Branch("BeamSpot_z0", &BeamSpot_z0, "BeamSpot_z0/F");
    tree_out->Branch("BeamSpot_BeamWidthX", &BeamSpot_BeamWidthX, "BeamSpot_BeamWidthX/F");
    tree_out->Branch("BeamSpot_BeamWidthY", &BeamSpot_BeamWidthY, "BeamSpot_BeamWidthY/F");


    ////////////////////////////// PRIMARY VERTEX BRANCHES //////////////////////////////

    tree_out->Branch("nPV", &nPV, "nPV/I");
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
    tree_out->Branch("IsoTrackSel_isHighPurityTrack", IsoTrackSel_isHighPurityTrack, "IsoTrackSel_isHighPurityTrack[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidTrackerHits", IsoTrackSel_numberOfValidTrackerHits, "IsoTrackSel_numberOfValidTrackerHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidPixelHits", IsoTrackSel_numberOfValidPixelHits, "IsoTrackSel_numberOfValidPixelHits[nIsoTrack]/I");
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

    tree_out->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
    tree_out->Branch("GenParticleSel_pt", GenParticleSel_pt, "GenParticleSel_pt[nGenParticle]/F");
    tree_out->Branch("GenParticleSel_eta", GenParticleSel_eta, "GenParticleSel_eta[nGenParticle]/F");
    tree_out->Branch("GenParticleSel_phi", GenParticleSel_phi, "GenParticleSel_phi[nGenParticle]/F");
    tree_out->Branch("GenParticleSel_dxy", GenParticleSel_dxy, "GenParticleSel_dxy[nGenParticle]/F");
    tree_out->Branch("GenParticleSel_pdgId", GenParticleSel_pdgId, "GenParticleSel_pdgId[nGenParticle]/I");



    //////////////////////////// ELECTRON CANDIDATE BRANCHES ////////////////////////////

    tree_out->Branch("nElectronCandidate", &nElectronCandidate, "nElectronCandidate/I");
    tree_out->Branch("ElectronCandidate_pt", ElectronCandidate_pt, "ElectronCandidate_pt[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_et", ElectronCandidate_et, "ElectronCandidate_et[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_eta", ElectronCandidate_eta, "ElectronCandidate_eta[nElectronCandidate]/F");    
    tree_out->Branch("ElectronCandidate_phi", ElectronCandidate_phi, "ElectronCandidate_phi[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_photonIdx", ElectronCandidate_photonIdx, "ElectronCandidate_photonIdx[nElectronCandidate]/I");
    tree_out->Branch("ElectronCandidate_isotrackIdx", ElectronCandidate_isotrackIdx, "ElectronCandidate_isotrackIdx[nElectronCandidate]/I");



    ////////////////////////////// MUON CANDIDATE BRANCHES /////////////////////////////

    tree_out->Branch("nMuonCandidate", &nMuonCandidate, "nMuonCandidate/I");
    tree_out->Branch("MuonCandidate_pt", MuonCandidate_pt, "MuonCandidate_pt[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_eta", MuonCandidate_eta, "MuonCandidate_eta[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_phi", MuonCandidate_phi, "MuonCandidate_phi[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_muonTriggerObjectIdx", MuonCandidate_muonTriggerObjectIdx, "MuonCandidate_muonTriggerObjectIdx[nMuonCandidate]/I");
    tree_out->Branch("MuonCandidate_isotrackIdx", MuonCandidate_isotrackIdx, "MuonCandidate_isotrackIdx[nMuonCandidate]/I");



}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::endJob() 
{


    std::cout << "The event is writen" << std::endl;
    file_out->cd();
    tree_out->Write();
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





DEFINE_FWK_MODULE(LongLivedAnalysis);
