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
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"


//#include "DataFormats/MuonReco/interface/MuonFwd.h" 

#include <string>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"


//=======================================================================================================================================================================================================================//
/////// -----------------------> DATA DEFINITION
//
////////// BRANCHES


//-> Isotrack Selection
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
Float_t IsoTrackSel_pfIsolationDR03[nIsoTrackMax];
Float_t IsoTrackSel_miniPFIsolation[nIsoTrackMax];
Int_t IsoTrackSel_isHighPurityTrack[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidTrackerHits[nIsoTrackMax];
Int_t IsoTrackSel_numberOfValidPixelHits[nIsoTrackMax];
// Derived:
Float_t IsoTrackSel_dxySignificance[nIsoTrackMax];


//-> Muon Selection
const Int_t nMuonMax = 100;
Int_t nMuon;
Float_t MuonSel_pt[nMuonMax];
Float_t MuonSel_eta[nMuonMax];


////////// OUTPUT FILE AND OUTPUT TREE
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
      edm::EDGetTokenT<edm::View<pat::Muon> >  theMuonCollection;   
      edm::EDGetTokenT<edm::View<pat::Photon> > thePhotonCollection;
      edm::EDGetTokenT<edm::View<pat::IsolatedTrack> >  theIsoTrackCollection;

};
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
LongLivedAnalysis::LongLivedAnalysis(const edm::ParameterSet& iConfig)
{
   usesResource("TFileService");
   
   parameters = iConfig;
   theMuonCollection = consumes<edm::View<pat::Muon> >  (parameters.getParameter<edm::InputTag>("MuonCollection"));
   thePhotonCollection = consumes<edm::View<pat::Photon> > (parameters.getParameter<edm::InputTag>("PhotonCollection"));
   theIsoTrackCollection = consumes<edm::View<pat::IsolatedTrack> >  (parameters.getParameter<edm::InputTag>("IsoTrackCollection"));
 
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

   //////////////////////////////// GET THE COLLECTIONS ////////////////////////////////
   edm::Handle<edm::View<pat::Muon> > muons;
   edm::Handle<edm::View<pat::Photon> > photons;
   edm::Handle<edm::View<pat::IsolatedTrack> > isotracks;

   iEvent.getByToken(theMuonCollection, muons);
   iEvent.getByToken(thePhotonCollection, photons);
   iEvent.getByToken(theIsoTrackCollection, isotracks);



   ///////////////////////////////// ISOTRACK FEATURES /////////////////////////////////

   nIsoTrack = isotracks->size();

   // Loop over the isotracks
   for (size_t i = 0; i < isotracks->size(); ++i){

       const pat::IsolatedTrack & isotrack = (*isotracks)[i];

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
       //
       const reco::HitPattern &hits = isotrack.hitPattern();

       IsoTrackSel_numberOfValidTrackerHits[i] = hits.numberOfValidTrackerHits();
       IsoTrackSel_numberOfValidPixelHits[i] = hits.numberOfValidPixelHits();

   }


   /////////////////////////////////// MUON FEATURES ///////////////////////////////////

   nMuon = muons->size();

   // Loop over the muons
   for (size_t i = 0; i < muons->size(); ++i){

       const pat::Muon & muon = (*muons)[i];

       MuonSel_pt[i] = muon.pt();
       MuonSel_eta[i] = muon.eta(); 
   }



   /////////////////////////////////// FILL THE TREE ///////////////////////////////////
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


    // Isotrack branches definition
    tree_out->Branch("nIsoTrack", &nIsoTrack, "nIsoTrack/I");
    tree_out->Branch("IsoTrackSel_pt", IsoTrackSel_pt, "IsoTrackSel_pt[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_eta", IsoTrackSel_eta, "IsoTrackSel_eta[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_phi", IsoTrackSel_phi, "IsoTrackSel_phi[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxy", IsoTrackSel_dxy, "IsoTrackSel_dxy[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxyError", IsoTrackSel_dxyError, "IsoTrackSel_dxyError[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dz", IsoTrackSel_dz, "IsoTrackSel_dz[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dzError", IsoTrackSel_dzError, "IsoTrackSel_dzError[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_dxySignificance", IsoTrackSel_dxySignificance, "IsoTrackSel_dxySignificance[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_pfIsolationDR03", IsoTrackSel_pfIsolationDR03, "IsoTrackSel_pfIsolationDR03[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_miniPFIsolation", IsoTrackSel_miniPFIsolation, "IsoTrackSel_miniPFIsolation[nIsoTrack]/F");
    tree_out->Branch("IsoTrackSel_isHighPurityTrack", IsoTrackSel_isHighPurityTrack, "IsoTrackSel_isHighPurityTrack[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidTrackerHits", IsoTrackSel_numberOfValidTrackerHits, "IsoTrackSel_numberOfValidTrackerHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidPixelHits", IsoTrackSel_numberOfValidPixelHits, "IsoTrackSel_numberOfValidPixelHits[nIsoTrack]/I");

    // Muon branches definition
    tree_out->Branch("nMuon", &nMuon, "nMuon/I");
    tree_out->Branch("MuonSel_pt", MuonSel_pt, "MuonSel_pt[nMuon]/F");
    tree_out->Branch("MuonSel_eta", MuonSel_eta, "MuonSel_pt[nMuon]/F");


}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::endJob() 
{

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
