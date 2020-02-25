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


#include "MyAnalysis/IFCALongLivedAnalysis/interface/llCandidateDataFormat.h"


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

bool goodElectron(const pat::Electron & electron)
{

    float ecal_energy_inverse = 1.0/electron.ecalEnergy();
    float eSCoverP = electron.eSuperClusterOverP();

    // Barrel cuts:
    if (fabs(electron.superCluster()->eta()) <= 1.479)
    {

       // Combined isolation:
       //float comIso = (electron.dr03TkSumPt() + max(0., electron.dr03EcalRecHitSumEt() - 1.) + electron.dr03HcalTowerSumEt() ) / electron.pt()

       if (electron.full5x5_sigmaIetaIeta() > 0.011) { return false; }
       if (fabs(electron.deltaEtaSeedClusterTrackAtVtx()) > 0.00477) { return false; }
       if (fabs(electron.deltaPhiSuperClusterTrackAtVtx()) > 0.222) { return false; }
       if (electron.hadronicOverEm() > 0.298) { return false; }
       if (fabs(1.0 - eSCoverP)*ecal_energy_inverse > 0.241) {return false; }
       if (electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) > 1) { return false;}
       if (electron.passConversionVeto() == 0) {return false;}


    // Endcap cuts
    } else if (fabs(electron.superCluster()->eta()) > 1.479) {


       if (electron.full5x5_sigmaIetaIeta() > 0.0314) { return false; }
       if (fabs(electron.deltaEtaSeedClusterTrackAtVtx()) > 0.00868) { return false; }
       if (fabs(electron.deltaPhiSuperClusterTrackAtVtx()) > 0.213) { return false; }
       if (electron.hadronicOverEm() > 0.101) { return false; }
       if (fabs(1.0 - eSCoverP)*ecal_energy_inverse > 0.14) {return false; }
       if (electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) > 1) { return false;}
       if (electron.passConversionVeto() == 0) {return false;}

    }


    if (fabs(electron.eta()) > 2.4) {return false;}

    return true;

}


bool goodMediumMuon( pat::Muon muon)
{

    if (!muon.isLooseMuon()) { return false;}
    if (muon.innerTrack()->validFraction() < 0.8) { return false; }

    bool additional_requirements = false;

    if ( ( muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2() < 3 && muon.combinedQuality().chi2LocalPosition < 12 && muon.combinedQuality().trkKink < 20 && muon.segmentCompatibility() > 0.303) 
         || muon.segmentCompatibility() > 0.451)

    { additional_requirements = true; }

    if (!additional_requirements) {return false;}

    return true;
}


bool goodMuon( pat::Muon muon)
{

    //if (muon.pt() < 28) { return false; }
    if (fabs(muon.eta()) > 2) { return false;}

    return true;

}


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


float checkRelativeMatching(float ratio, float a, float b)
{

   // ratio: maximum ratio that is allowed for the matching to be true

   if (fabs(a - b)/a > ratio){ return false; }

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
Int_t nPU;
Int_t nPUTrue;
Float_t wPU;
Float_t genWeight;

//-> TRIGGER TAGS
Bool_t Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10;
Bool_t Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15;

//-> PRIMARY VERTEX SELECTION
Int_t nPV;
Int_t nTruePV;
Float_t PV_vx;
Float_t PV_vy;
Float_t PV_vz;
Float_t PV_xyFromBS;
Float_t PV_zFromBS;


// -> REFITTED PRIMARY VERTEX
Float_t RefittedPV_vx;
Float_t RefittedPV_vy;
Float_t RefittedPV_vz;
Int_t RefittedPV_nPFTrack;
Int_t RefittedPV_nLostTrack;
Int_t RefittedPV_nExcludedTrack;


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

//-> ELECTRON SELECTION
const Int_t nElectronMax = 100;
Int_t nElectron;
Float_t ElectronSel_pt[nElectronMax];
Float_t ElectronSel_et[nElectronMax];
Float_t ElectronSel_eta[nElectronMax];
Float_t ElectronSel_phi[nElectronMax];
Float_t ElectronSel_hadronicOverEm[nElectronMax];
Float_t ElectronSel_full5x5_sigmaIetaIeta[nElectronMax];
Int_t ElectronSel_isEB[nElectronMax];
Int_t ElectronSel_isEE[nElectronMax];
Float_t ElectronSel_r9[nElectronMax];
Float_t ElectronSel_trackIso[nElectronMax];
Float_t ElectronSel_ecalIso[nElectronMax];
Float_t ElectronSel_hcalIso[nElectronMax];
Float_t ElectronSel_caloIso[nElectronMax];
Float_t ElectronSel_relIso[nElectronMax];
Float_t ElectronSel_dxy[nElectronMax];
Float_t ElectronSel_dxyError[nElectronMax];
Float_t ElectronSel_dxySignificance[nElectronMax];
Float_t ElectronSel_dB[nElectronMax];
Float_t ElectronSel_edB[nElectronMax];
Int_t ElectronSel_isLoose[nElectronMax];


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
Float_t MuonSel_fractionOfValidTrackerHits[nMuonMax];
Float_t MuonSel_normGlobalTrackChi2[nMuonMax];

// -> GENHIGGS
Float_t GenHiggs_pt;
Float_t GenHiggs_eta;
Float_t GenHiggs_phi;

//-> MUON TRIGGER OBJECT SELECTION
const Int_t nMuonTriggerObjectMax = 500;
Int_t nMuonTriggerObject;
Float_t MuonTriggerObjectSel_pt[nMuonTriggerObjectMax];
Float_t MuonTriggerObjectSel_eta[nMuonTriggerObjectMax];
Float_t MuonTriggerObjectSel_phi[nMuonTriggerObjectMax];

//-> GENLEPTON SELECTION
const Int_t nGenLeptonMax = 10;
Int_t nGenLepton;
Float_t GenLeptonSel_pt[nGenLeptonMax];
Float_t GenLeptonSel_et[nGenLeptonMax];
Float_t GenLeptonSel_eta[nGenLeptonMax];
Float_t GenLeptonSel_phi[nGenLeptonMax];
Int_t GenLeptonSel_pdgId[nGenLeptonMax];
Float_t GenLeptonSel_dxy[nGenLeptonMax];
Int_t GenLeptonSel_motherIdx[nGenLeptonMax];
Int_t GenLeptonSel_objectMatch[nGenLeptonMax];
Int_t GenLeptonSel_trackMatch[nGenLeptonMax];
Float_t GenLeptonSel_objectdR[nGenLeptonMax];
Float_t GenLeptonSel_trackdR[nGenLeptonMax];
Int_t GenLeptonSel_hasValidPair[nGenLeptonMax];
Float_t GenLeptonSel_pairdR[nGenLeptonMax];
Int_t GenLeptonSel_trackDegeneration[nGenLeptonMax];
Int_t GenLeptonSel_objectDegeneration[nGenLeptonMax];
Float_t GenLeptonSel_vx[nGenLeptonMax];
Float_t GenLeptonSel_vy[nGenLeptonMax];
Float_t GenLeptonSel_vz[nGenLeptonMax];


//-> GENNEUTRALINO SELECTION
const Int_t nGenNeutralinoMax = 500;
Int_t nGenNeutralino;
Float_t GenNeutralinoSel_pt[nGenNeutralinoMax];
Float_t GenNeutralinoSel_eta[nGenNeutralinoMax];
Float_t GenNeutralinoSel_phi[nGenNeutralinoMax];
Int_t GenNeutralinoSel_pdgId[nGenNeutralinoMax];
Int_t GenNeutralinoSel_decaypdgId[nGenNeutralinoMax];
Bool_t GenNeutralinoSel_passAcceptance[nGenNeutralinoMax];
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
const Int_t nElectronCandidateMax = 1000;
Int_t nElectronCandidate;
Float_t ElectronCandidate_pt[nElectronCandidateMax];
Float_t ElectronCandidate_et[nElectronCandidateMax];
Float_t ElectronCandidate_eta[nElectronCandidateMax];
Float_t ElectronCandidate_phi[nElectronCandidateMax];
Float_t ElectronCandidate_dxy[nElectronCandidateMax];
Float_t ElectronCandidate_dxyError[nElectronCandidateMax];
Float_t ElectronCandidate_dxySignificance[nElectronCandidateMax];
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
Float_t MuonCandidate_dxySignificance[nMuonCandidateMax];
Float_t MuonCandidate_triggerPt[nMuonCandidateMax];
Int_t MuonCandidate_muonTriggerObjectIdx[nMuonCandidateMax];
Int_t MuonCandidate_isotrackIdx[nMuonCandidateMax];
Int_t MuonCandidate_pvAssociationQuality[nMuonCandidateMax];
Float_t MuonCandidate_ptDiff[nMuonCandidateMax];


// -> All EE candidates
Int_t nEE;
Int_t EE_idxA[20];
Int_t EE_idxB[20];
Float_t EE_Lxy[20];
Float_t EE_Ixy[20];
Float_t EE_trackDxy[20];
Float_t EE_trackIxy[20];
Float_t EE_mass[20];
Float_t EE_normalizedChi2[20];
Float_t EE_leadingPt[20];
Float_t EE_subleadingPt[20];
Float_t EE_leadingEt[20];
Float_t EE_subleadingEt[20];
Float_t EE_cosAlpha[20];
Float_t EE_dPhi[20];
Float_t EE_relisoA[20];
Float_t EE_relisoB[20];

// -> EE candidates that survive to baseline selection
Int_t nEEBase;
Int_t EEBase_maxIxy;
Int_t EEBase_idxA[20];
Int_t EEBase_idxB[20];
Float_t EEBase_Lxy[20];
Float_t EEBase_Ixy[20];
Float_t EEBase_trackDxy[20];
Float_t EEBase_trackIxy[20];
Float_t EEBase_vx[20];
Float_t EEBase_vy[20];
Float_t EEBase_mass[20];
Float_t EEBase_normalizedChi2[20];
Float_t EEBase_leadingPt[20];
Float_t EEBase_subleadingPt[20];
Float_t EEBase_leadingEt[20];
Float_t EEBase_subleadingEt[20];
Float_t EEBase_cosAlpha[20];
Float_t EEBase_dPhi[20];
Float_t EEBase_relisoA[20];
Float_t EEBase_relisoB[20];
Float_t EEBase_refittedDxy[20];
Float_t EEBase_refittedIxy[20];
Int_t EEBase_fromPVA[20];
Int_t EEBase_fromPVB[20];
Int_t EEBase_PVAssociation[20];


// -> All MM candidates
Int_t nMM;
Int_t MM_idxA[20];
Int_t MM_idxB[20];
Float_t MM_Lxy[20];
Float_t MM_Ixy[20];
Float_t MM_trackDxy[20];
Float_t MM_trackIxy[20];
Float_t MM_mass[20];
Float_t MM_normalizedChi2[20];
Float_t MM_leadingPt[20];
Float_t MM_subleadingPt[20];
Float_t MM_cosAlpha[20];
Float_t MM_dPhi[20];
Float_t MM_relisoA[20];
Float_t MM_relisoB[20];


// -> MM candidates that survive to baseline selection
Int_t nMMBase;
Int_t MMBase_maxIxy;
Int_t MMBase_idxA[20];
Int_t MMBase_idxB[20];
Float_t MMBase_Lxy[20];
Float_t MMBase_Ixy[20];
Float_t MMBase_trackDxy[20];
Float_t MMBase_trackIxy[20];
Float_t MMBase_vx[20];
Float_t MMBase_vy[20];
Float_t MMBase_mass[20];
Float_t MMBase_normalizedChi2[20];
Float_t MMBase_leadingPt[20];
Float_t MMBase_subleadingPt[20];
Float_t MMBase_cosAlpha[20];
Float_t MMBase_dPhi[20];
Float_t MMBase_relisoA[20];
Float_t MMBase_relisoB[20];
Float_t MMBase_refittedDxy[20];
Float_t MMBase_refittedIxy[20];
Int_t MMBase_fromPVA[20];
Int_t MMBase_fromPVB[20];
Int_t MMBase_PVAssociation[20];


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
      bool _isData;
      bool _BSMode;
      edm::ParameterSet parameters;

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

      // Gen collection
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  theGenParticleCollection;      
      edm::EDGetTokenT<GenEventInfoProduct>  theGenEventInfoProduct;

      // PU reweighting
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  thePileUpSummary;
      edm::LumiReWeighting lumi_weights = edm::LumiReWeighting("2016MCPileupHistogram.root", "2016DataPileupHistogram.root", "pileup", "pileup");
      //lumi_weights = edm::LumiReWeighting("2016MCPileupHistogram.root", "2016DataPileupHistogram.root", "pileup", "pileup");

      //"Global" variables
      std::vector<int> iT; // track indexes
      edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;

      // Class functions
      bool buildLLcandidate(edm::Handle<edm::View<pat::IsolatedTrack> > const& isotracks, int idxA, int idxB, bool isEE);
      bool isLooseElectron(const pat::Electron & electron); 
      bool isMediumElectron(const pat::Electron & electron); 
      bool passIsotrackSelection(const pat::IsolatedTrack &track);
      bool passPhotonSelection(const pat::Photon &photon);
      bool passL2MuonSelection( pat::TriggerObjectStandAlone obj); 
      bool passMuonSelection(const pat::Muon &muon);
      bool passBaselineSelection(llCandidate llc);
      float computeDxy(const pat::IsolatedTrack & track, const reco::Vertex pv);
      reco::Vertex getSVCandidate(const pat::PackedCandidateRef &pckCandA, const pat::PackedCandidateRef &pckCandB);
      float computeDxyError(const pat::IsolatedTrack & track, const reco::Vertex pv);

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
   //iEvent.getByToken(triggerPrescales_, triggerPrescales);
   iEvent.getByToken(theBeamSpot, beamSpot);

   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);

   if (!_isData){
      
       iEvent.getByToken(theGenParticleCollection, genParticles);    
       iEvent.getByToken(theGenEventInfoProduct, genEvtInfo);
       iEvent.getByToken(thePileUpSummary, puInfoH);


   }
   

   /////////////////////////////////// EVENT INFO //////////////////////////////////////

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

   //////////////////////////////////// PILE UP ////////////////////////////////////////

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

   //////////////////////////////////// BEAM SPOT //////////////////////////////////////

   reco::BeamSpot beamSpotObject = *beamSpot;
   BeamSpot_x0 = beamSpotObject.x0();
   BeamSpot_y0 = beamSpotObject.y0();
   BeamSpot_z0 = beamSpotObject.z0();
   BeamSpot_BeamWidthX = beamSpotObject.BeamWidthX();
   BeamSpot_BeamWidthY = beamSpotObject.BeamWidthY();


   /////////////////////////////// TRIGGER ACCEPTANCE //////////////////////////////////

   // Trigger names declaration
   std::string muonTriggerName;
   std::string photonTriggerName;
 
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);


   // Loop over trigger versions
   std::string version;

   // muons:
   for (int s = 0; s < 11; s++)
   {

      version = "HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v" + std::to_string(s);
      if (names.size() != names.triggerIndex(version)) 
      {

          muonTriggerName = version;
          break;         

      }
   }

   // electrons:
   for (int s = 0; s < 11; s++)
   {

      version = "HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v" + std::to_string(s);
      if (names.size() != names.triggerIndex(version)) 
      {

          photonTriggerName = version;
          break;         

      }
   }

   
   Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 = triggerBits->accept(names.triggerIndex(muonTriggerName));
   Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15 = triggerBits->accept(names.triggerIndex(photonTriggerName));



   ////////////////////////////// MUON TRIGGER OBJECTS /////////////////////////////////
 
   std::vector<int> iMT; // muon trigger object indexes 


   // Loop to get Muon Trigger objects:
   for (size_t i = 0; i < triggerObjects->size(); i++) 
   {
       

       pat::TriggerObjectStandAlone obj = (*triggerObjects)[i];


       obj.unpackPathNames(names);
       obj.unpackFilterLabels(iEvent, *triggerBits);    

       
       bool isMuonTriggerObject = obj.hasPathName(muonTriggerName, true, true );

       if (!isMuonTriggerObject) { continue; }

       if (passL2MuonSelection(obj)){ iMT.push_back(i); }

   }
   


   // Sort the muon trigger objects by pt:
   std::sort( std::begin(iMT), std::end(iMT), [&](int i1, int i2){ return triggerObjects->at(i1).pt() > triggerObjects->at(i2).pt(); });



   // Fill the muon trigger objects features:
   nMuonTriggerObject = iMT.size();

   for (size_t i = 0; i < iMT.size(); i++){
   
       pat::TriggerObjectStandAlone obj = (*triggerObjects)[iMT.at(i)];

       obj.unpackPathNames(names);
       obj.unpackFilterLabels(iEvent, *triggerBits);

       MuonTriggerObjectSel_pt[i] = obj.pt();
       MuonTriggerObjectSel_eta[i] = obj.eta();
       MuonTriggerObjectSel_phi[i] = obj.phi();
   

   }


    

   ////////////////////////////// PRIMARY VERTEX FEATURES //////////////////////////////

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

   // Distances with respect to the beam spot:
   PV_xyFromBS = sqrt((PV_vx - BeamSpot_x0)*(PV_vx - BeamSpot_x0) + (PV_vy - BeamSpot_y0)*(PV_vy - BeamSpot_y0));
   PV_zFromBS = PV_vz - BeamSpot_z0;


   ///////////////////////////////// ISOTRACK FEATURES /////////////////////////////////
   iT.clear();
   
   for (size_t i = 0; i < isotracks->size(); i++){

       const pat::IsolatedTrack & isotrack = (*isotracks)[i];

       if (!passIsotrackSelection(isotrack)){ continue; } // some quality cuts

       iT.push_back(i);

   }


   nIsoTrack = iT.size(); // number of isotracks

   // Sort the isotracks by pt
   std::sort( std::begin(iT), std::end(iT), [&](int i1, int i2){ return isotracks->at(i1).pt() > isotracks->at(i2).pt(); });


   // Loop over the isotracks
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
           
           // Impact parameter info:
           // Uncertainties are checked by means of a Transient Track Builder in method 'computeDxyError()' 
           IsoTrackSel_dxy[i] = fabs((*pckCand).dxy(thePrimaryVertex.position()));
           IsoTrackSel_dxyError[i] = computeDxyError(isotrack, thePrimaryVertex);
           IsoTrackSel_dz[i] = (*pckCand).dz(thePrimaryVertex.position());
           IsoTrackSel_dzError[i] = (*pckCand).dzError(); 
           IsoTrackSel_dxySignificance[i] = fabs((*pckCand).dxy(thePrimaryVertex.position()))/computeDxyError(isotrack, thePrimaryVertex);


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

       const pat::Photon & photon = (*photons)[i];
       
       if (!passPhotonSelection(photon)) { continue;}

       iP.push_back(i);

   }

   // Sort good lepton indexes by pt
   std::sort( std::begin(iP), std::end(iP), [&](int i1, int i2){ return photons->at(i1).et() > photons->at(i2).et(); });


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


   }


   ///////////////////////////////// ELECTRON FEATURES /////////////////////////////////

   /*
   std::vector<int> iE; // electron indexes


   // Select good photons
   for (size_t i = 0; i < electrons->size(); i++){

       const pat::Electron & electron = (*electrons)[i];

       // this is the place to put any preselection if required
       //if (goodElectron(electron)) { iE.push_back(i);}
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
       ElectronSel_hadronicOverEm[i] = electron.hadronicOverEm();
       ElectronSel_full5x5_sigmaIetaIeta[i] = electron.full5x5_sigmaIetaIeta();
       ElectronSel_isEB[i] = electron.isEB();
       ElectronSel_isEE[i] = electron.isEE();
       ElectronSel_r9[i] = electron.r9();

       ElectronSel_trackIso[i] = electron.trackIso();
       ElectronSel_hcalIso[i] = electron.hcalIso();
       ElectronSel_ecalIso[i] = electron.ecalIso();
       ElectronSel_caloIso[i] = electron.caloIso();
       ElectronSel_relIso[i] = electron.caloIso()/electron.pt();

       ElectronSel_dxyError[i] = electron.dxyError();
       ElectronSel_dxy[i] = electron.gsfTrack()->dxy();
       ElectronSel_dxySignificance[i] = fabs(ElectronSel_dxy[i])/electron.dxyError();

       ElectronSel_dB[i] = electron.dB();
       ElectronSel_edB[i] = electron.edB();

       ElectronSel_isLoose[i] = electron.electronID("cutBasedElectronID-Summer16-80X-V1-loose");


   }

   */

   ///////////////////////////////// MUON FEATURES /////////////////////////////////
   //PABLO: CHECK WHETHER WE WANT TO HAVE THIS ACTIVATED ALL THE TIME
   std::vector<int> iM; // muon indexes


   // Select good muons
   for (size_t i = 0; i < muons->size(); i++){

       const pat::Muon & muon = (*muons)[i];
       
       // this is the place to put any preselection if required
       if (goodMuon(muon)) { iM.push_back(i);}

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
       MuonSel_isGoodMediumMuon[i] = goodMediumMuon(muon);
       MuonSel_isPFMuon[i] = muon.isPFMuon();


       // Quality features:

       if (muon.isGlobalMuon() || muon.isTrackerMuon()) {MuonSel_fractionOfValidTrackerHits[i] = muon.innerTrack()->validFraction(); }

       if (muon.isGlobalMuon()){
 
          MuonSel_normGlobalTrackChi2[i] = muon.globalTrack()->normalizedChi2();

       } else {

          MuonSel_normGlobalTrackChi2[i] = -99;
       }

       MuonSel_dB[i] = muon.dB();
       MuonSel_edB[i] = muon.edB();
       MuonSel_dBSignificance[i] = muon.dB()/muon.edB();

   }


   //////////////////////////////// GENPARTICLE FEATURES ///////////////////////////////
   //PABLO: MAKE A POLICY FOR GENERATED PARTICLES
   std::vector<int> iGL; // Generated lepton indexes
   std::vector<int> iGN; // Generated neutralino indexes
   int iH = -1; // Higgs index

   if (!_isData)
   {

       for(size_t i = 0; i < genParticles->size(); i++) {


	   const reco::GenParticle &genparticle = (*genParticles)[i];

	   // Check if it is a lepton comming from a long lived neutralino:
	   if (isLongLivedLepton(genparticle)){ iGL.push_back(i); continue; }

	   // Check if it is a longlived neutralino (correct id and after all radiative emission)
	   if (abs(genparticle.pdgId()) == 1000022 && (abs(genparticle.daughter(0)->pdgId()) != 1000022) && (abs(genparticle.daughter(0)->pdgId()) != 22)){
	       iGN.push_back(i);
	       continue;
	   } else if (abs(genparticle.pdgId()) == 54 && (abs(genparticle.daughter(0)->pdgId()) != 54) && (abs(genparticle.daughter(0)->pdgId()) != 22)){
	       iGN.push_back(i);
	       continue;
	   }

           if (genparticle.pdgId() == 35 && genparticle.daughter(0)->pdgId() != 35 && genparticle.daughter(0)->pdgId() != 22)
           {
               iH = i;
           }

       } 

       // Higgs details
       GenHiggs_pt = -99;
       GenHiggs_eta = -99;
       GenHiggs_phi = -99;

       if (iH != -1)   
       {
           
           GenHiggs_pt =(*genParticles)[iH].pt();
           GenHiggs_eta = (*genParticles)[iH].eta();
           GenHiggs_phi = (*genParticles)[iH].phi();
       }

       // Number of generated leptons:
       nGenLepton = iGL.size();

       // Number of generated neutralinos:
       nGenNeutralino = iGN.size();

       // Sort the leptons and neutralinos by pt
       std::sort( std::begin(iGL), std::end(iGL), [&](int i1, int i2){ return genParticles->at(i1).pt() > genParticles->at(i2).pt(); });
       std::sort( std::begin(iGN), std::end(iGN), [&](int i1, int i2){ return genParticles->at(i1).pt() > genParticles->at(i2).pt(); });


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
	   while((abs(m->pdgId()) == 100022) or (abs(m->pdgId()) == 54)){ m = m->mother(); }

	   GenNeutralinoSel_Lxy[i] = sqrt((m->vx()-genparticle.daughter(0)->vx())*(m->vx()-genparticle.daughter(0)->vx()) + (m->vy()-genparticle.daughter(0)->vy())*(m->vy()-genparticle.daughter(0)->vy()));

       }


       // variables to avoid radiative effects on leptons
       int rad = 0;
       int rad_2 = 0;


       // Loop over the selected genleptons
       for(size_t i = 0; i < iGL.size(); i++){

	   const reco::GenParticle &genparticle = (*genParticles)[iGL.at(i)];

	   GenLeptonSel_pdgId[i] = genparticle.pdgId();
	   GenLeptonSel_dxy[i] = dxy_value(genparticle, thePrimaryVertex);
	   GenLeptonSel_vx[i] = genparticle.vx();
	   GenLeptonSel_vy[i] = genparticle.vy();
	   GenLeptonSel_vz[i] = genparticle.vz();




	   // Define the mother index
	   GenLeptonSel_motherIdx[i] = -99;
	   if(genparticle.mother()->pt() == GenNeutralinoSel_pt[0]){
	     if (abs(genparticle.pdgId()) == 11) GenNeutralinoSel_decaypdgId[0] = 11;
	     if (abs(genparticle.pdgId()) == 13) GenNeutralinoSel_decaypdgId[0] = 13;
	       GenLeptonSel_motherIdx[i] = 0;

	   } else if (genparticle.mother()->pt() == GenNeutralinoSel_pt[1]){
	     if (abs(genparticle.pdgId()) == 11) GenNeutralinoSel_decaypdgId[1] = 11;
	     if (abs(genparticle.pdgId()) == 13) GenNeutralinoSel_decaypdgId[1] = 13;
	       GenLeptonSel_motherIdx[i] = 1;

	   }


	   // Get the last genlepton (to avoid radiative effects):
	   if (genparticle.numberOfDaughters() > 0){

	       // look for the correct daughter
	       for (size_t j = 0; j < genparticle.numberOfDaughters(); j++)
	       {
		   if (genparticle.pdgId() == genparticle.daughter(j)->pdgId()) {rad = j; break; }
	       }

	       const reco::Candidate *d = genparticle.daughter(rad);

	       //while(d->numberOfDaughters()> 0 && d->daughter(0)->pdgId() == d->pdgId()){ d = d->daughter(0); }
	       while(d->numberOfDaughters()> 0)
	       {

		   for(size_t j = 0; j < d->numberOfDaughters(); j++)
		   {
		       if (d->pdgId() == d->daughter(j)->pdgId()) {rad_2 = j; break;}
		   }

		   d = d->daughter(rad_2);

	       }

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

   }

   //////////////////////////////// ACCEPTANCE CRITERIA ////////////////////////////////

   if (!_isData)
     {

       for (size_t iN = 0; iN < iGN.size(); iN++) {//loop to the 2 neutralinos

	 bool passLeadingElectron = false;
	 bool passAccLL = true;
	 for (size_t i = 0; i < iGL.size(); i++){

	   if (GenLeptonSel_motherIdx[i] != int(iN)) continue;//I only consider daughter leptons from the neutralino in the loop

	   if (abs(GenLeptonSel_pdgId[i]) == 11 && GenLeptonSel_et[i] > 40){ passLeadingElectron = true;}//isn't this a little silly if we have more than one LL candidate?

	   if (abs(GenLeptonSel_pdgId[i]) == 11 && GenLeptonSel_et[i] < 25){ passAccLL = false; /*break;*/}
	   if (abs(GenLeptonSel_pdgId[i]) == 13 && GenLeptonSel_pt[i] < 26){ passAccLL = false; /*break;*/}
	   if (fabs(GenLeptonSel_eta[i]) > 2){ passAccLL = false; /*break;*/}

	 }

	 if (passLeadingElectron == false && GenNeutralinoSel_decaypdgId[iN] == 11) {passAccLL = false;}
	 if (GenNeutralinoSel_Lxy[iN] > 50){ passAccLL = false;}

	 GenNeutralinoSel_passAcceptance[iN] = passAccLL;

       }


       passAcceptanceCriteria = true; // true by default
       for (size_t iN = 0; iN < iGN.size(); iN++) {//loop to the 2 neutralinos
	 if (!GenNeutralinoSel_passAcceptance[iN]) passAcceptanceCriteria = false;
       }
 

   }
   

   //////////////////////////////////////// MET ////////////////////////////////////////
   //CHECK: DO WE WANT TO USE PUPPIMET
   const pat::MET &met = (*METs)[0]; // access the MET object

   MET_pt = met.pt();
   MET_phi = met.phi();

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
           if(!passIsotrackSelection(isotrack)){ continue; }


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
                   matching_type = 0;
                   scmin = sc;
                   tmin = t;

               }
           }


           // Loop over the muon trigger objecs
           for (int to = 0; to < nMuonTriggerObject; to++){

               // pass if the trigger muon is associated to other track
               if(std::find(matched_triggerObjects.begin(), matched_triggerObjects.end(), to) != matched_triggerObjects.end()){ continue; }

               // ------------ trigger Object matching --------------
               dR = getDeltaR(isotrack.phi(), isotrack.eta(), MuonTriggerObjectSel_phi[to], MuonTriggerObjectSel_eta[to]);

               if (dR < dRMin){

                   dRMin = dR;
                   matching_type = 1;
                   tomin = to;
                   tmin = t;

               }
           }
       }

       

       if (dRMin > dRThreshold){ break; } // Here we go out the while loop

       if (matching_type == 0){

           li = matched_SC.size();
           ElectronCandidate_pt[li] = (*isotracks)[iT.at(tmin)].pt();
           ElectronCandidate_eta[li] = (*isotracks)[iT.at(tmin)].eta();
           ElectronCandidate_phi[li] = (*isotracks)[iT.at(tmin)].phi();
           ElectronCandidate_et[li] = (*photons)[iP.at(scmin)].et();
           ElectronCandidate_photonIdx[li] = scmin;
           ElectronCandidate_isotrackIdx[li] = tmin;
           ElectronCandidate_dxy[li] = IsoTrackSel_dxy[tmin];
           ElectronCandidate_dxyError[li] = IsoTrackSel_dxyError[tmin];
           ElectronCandidate_dxySignificance[li] = IsoTrackSel_dxySignificance[tmin];
 
           ElectronCandidate_pvAssociationQuality[li] = (*isotracks)[iT.at(tmin)].packedCandRef()->pvAssociationQuality();
           ElectronCandidate_ptDiff[li] = (*isotracks)[iT.at(tmin)].pt() - (*isotracks)[iT.at(tmin)].packedCandRef()->pseudoTrack().pt();

           matched_SC.push_back(scmin); matched_tracks.push_back(tmin);


       } else if (matching_type == 1){

           li = matched_triggerObjects.size();
           
           MuonCandidate_pt[li] = (*isotracks)[iT.at(tmin)].pt();
           MuonCandidate_eta[li] = (*isotracks)[iT.at(tmin)].eta();
           MuonCandidate_phi[li] = (*isotracks)[iT.at(tmin)].phi();
           MuonCandidate_triggerPt[li] = MuonTriggerObjectSel_pt[tomin];
           MuonCandidate_muonTriggerObjectIdx[li] = tomin; // Corrected
           MuonCandidate_isotrackIdx[li] = tmin;
           MuonCandidate_dxy[li] = IsoTrackSel_dxy[tmin];
           MuonCandidate_dxyError[li] = IsoTrackSel_dxyError[tmin];
           MuonCandidate_dxySignificance[li] = IsoTrackSel_dxySignificance[tmin];

           MuonCandidate_pvAssociationQuality[li] = (*isotracks)[iT.at(tmin)].packedCandRef()->pvAssociationQuality();
           MuonCandidate_ptDiff[li] = (*isotracks)[iT.at(tmin)].pt() - (*isotracks)[iT.at(tmin)].packedCandRef()->pseudoTrack().pt();
           matched_triggerObjects.push_back(tomin); matched_tracks.push_back(tmin);

       }
         
   }

   nElectronCandidate = matched_SC.size();
   nMuonCandidate = matched_triggerObjects.size();


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
   nEEBase = 0;
   EEBase_maxIxy = 0;
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
         //PABLO: CHECK THIS CLASS
         llCandidate testcandidate(thePrimaryVertex, theTransientTrackBuilder, it_i, it_j, true);

         if (!testcandidate.canFitVertex || !testcandidate.hasValidVertex) { continue ;} 

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

     llCandidate eeCandidate(thePrimaryVertex, theTransientTrackBuilder, it_A, it_B, true);
     // Additionally, for electrons we have Et:
     eeCandidate.leadingEt = (ElectronCandidate_et[min_i] > ElectronCandidate_et[min_j])? ElectronCandidate_et[min_i]: ElectronCandidate_et[min_j];
     eeCandidate.subleadingEt = (ElectronCandidate_et[min_i] < ElectronCandidate_et[min_j])? ElectronCandidate_et[min_i]: ElectronCandidate_et[min_j];

     if (!_BSMode){

        EE_idxA[nEE] = min_i;
        EE_idxB[nEE] = min_j;
        EE_Lxy[nEE] = eeCandidate.vertexLxy;
        EE_Ixy[nEE] = eeCandidate.vertexIxy;
        EE_trackDxy[nEE] = eeCandidate.trackDxy;
        EE_trackIxy[nEE] = eeCandidate.trackIxy;
        EE_normalizedChi2[nEE] = eeCandidate.normalizedChi2;
        EE_mass[nEE] = eeCandidate.mass;
        EE_leadingPt[nEE] = eeCandidate.leadingPt;
        EE_subleadingPt[nEE] = eeCandidate.subleadingPt;
        EE_cosAlpha[nEE] = eeCandidate.cosAlpha;
        EE_dPhi[nEE] = eeCandidate.dPhi;
        EE_relisoA[nEE] = eeCandidate.relisoA;
        EE_relisoB[nEE] = eeCandidate.relisoB;
        EE_leadingEt[nEE] = eeCandidate.leadingEt;
        EE_subleadingEt[nEE] = eeCandidate.subleadingEt;

     }
     nEE++;

     // -> Fill candidates that pass baseline selection:
     if ( passBaselineSelection(eeCandidate) ) {

        if (_BSMode) {

           leptonTracks.push_back(it_A); leptonTracks.push_back(it_B);

           EEBase_idxA[nEEBase] = min_i;
           EEBase_idxB[nEEBase] = min_j;
           EEBase_Lxy[nEEBase] = eeCandidate.vertexLxy;
           EEBase_Ixy[nEEBase] = eeCandidate.vertexIxy;
           EEBase_trackDxy[nEEBase] = eeCandidate.trackDxy;
           EEBase_trackIxy[nEEBase] = eeCandidate.trackIxy;
           EEBase_vx[nEEBase] = eeCandidate.vx;
           EEBase_vy[nEEBase] = eeCandidate.vy;
           EEBase_normalizedChi2[nEEBase] = eeCandidate.normalizedChi2;
           EEBase_mass[nEEBase] = eeCandidate.mass;
           EEBase_leadingPt[nEEBase] = eeCandidate.leadingPt;
           EEBase_subleadingPt[nEEBase] = eeCandidate.subleadingPt;
           EEBase_cosAlpha[nEEBase] = eeCandidate.cosAlpha;
           EEBase_dPhi[nEEBase] = eeCandidate.dPhi;
           EEBase_relisoA[nEEBase] = eeCandidate.relisoA;
           EEBase_relisoB[nEEBase] = eeCandidate.relisoB;
           EEBase_leadingEt[nEEBase] = eeCandidate.leadingEt;
           EEBase_subleadingEt[nEEBase] = eeCandidate.subleadingEt;
           EEBase_fromPVA[nEEBase] = eeCandidate.fromPVA;
           EEBase_fromPVB[nEEBase] = eeCandidate.fromPVB;
           EEBase_PVAssociation[nEEBase] = eeCandidate.PVAssociation;

           if ( fabs(EEBase_trackIxy[nEEBase]) > fabs(EEBase_trackIxy[EEBase_maxIxy]) ) { EEBase_maxIxy = nEEBase; }

        }
        nEEBase++;

     }

   } // end while 


   /////////////////////////////////////////////////////////
   // --------------------------------------------------- //
   // ----------- mmCandidates reconstruction ----------- //
   // --------------------------------------------------- //
   /////////////////////////////////////////////////////////
   nMM = 0;
   nMMBase = 0;
   MMBase_maxIxy = 0;
   std::vector<double> pairedM; // muons that are already paired


   while(2*nMM < nMuonCandidate - 1 ){

     // Init control variables:
     minChi2 = 10000;
     min_i = 99;
     min_j = 99;

     for (int i = 0; i < nMuonCandidate; i++) {
       for (int j = i+1; j < nMuonCandidate; j++) {
 
         if (i == j) { continue; }
         if ( std::find(pairedM.begin(), pairedM.end(), i) != pairedM.end() ) {continue;}
         if ( std::find(pairedM.begin(), pairedM.end(), j) != pairedM.end() ) {continue;}

         const pat::IsolatedTrack & it_i = (*isotracks)[iT.at(MuonCandidate_isotrackIdx[i])];
         const pat::IsolatedTrack & it_j = (*isotracks)[iT.at(MuonCandidate_isotrackIdx[j])];

         llCandidate testcandidate(thePrimaryVertex, theTransientTrackBuilder, it_i, it_j, false);

         if (!testcandidate.canFitVertex || !testcandidate.hasValidVertex) { continue ;} 

         // Check if the Chi2 is lower:
         if (testcandidate.normalizedChi2 < minChi2) {
           minChi2 = testcandidate.normalizedChi2;
           min_i = i;
           min_j = j;
         }

       } // end j muon loop
     } // end i muon loop

     if (min_i == 99 || min_j == 99) { break; }
     pairedM.push_back(min_i);
     pairedM.push_back(min_j);

     // -> Get LLP Candidate variables:
     const pat::IsolatedTrack &it_A = (*isotracks)[iT.at(MuonCandidate_isotrackIdx[min_i])];
     const pat::IsolatedTrack &it_B = (*isotracks)[iT.at(MuonCandidate_isotrackIdx[min_j])]; 

     llCandidate mmCandidate(thePrimaryVertex, theTransientTrackBuilder, it_A, it_B, false);

     if (!_BSMode){

        MM_idxA[nMM] = min_i;
        MM_idxB[nMM] = min_j;
        MM_Lxy[nMM] = mmCandidate.vertexLxy;
        MM_Ixy[nMM] = mmCandidate.vertexIxy;
        MM_trackDxy[nMM] = mmCandidate.trackDxy;
        MM_trackIxy[nMM] = mmCandidate.trackIxy;
        MM_normalizedChi2[nMM] = mmCandidate.normalizedChi2;
        MM_mass[nMM] = mmCandidate.mass;
        MM_leadingPt[nMM] = mmCandidate.leadingPt;
        MM_subleadingPt[nMM] = mmCandidate.subleadingPt;
        MM_cosAlpha[nMM] = mmCandidate.cosAlpha;
        MM_dPhi[nMM] = mmCandidate.dPhi;
        MM_relisoA[nMM] = mmCandidate.relisoA;
        MM_relisoB[nMM] = mmCandidate.relisoB;

     }
     nMM++;

     // -> Fill candidates that pass baseline selection:
     if ( passBaselineSelection(mmCandidate) ) {

        if (_BSMode){
           
           leptonTracks.push_back(it_A); leptonTracks.push_back(it_B);

           MMBase_idxA[nMMBase] = min_i;
           MMBase_idxB[nMMBase] = min_j;
           MMBase_Lxy[nMMBase] = mmCandidate.vertexLxy;
           MMBase_Ixy[nMMBase] = mmCandidate.vertexIxy;
           MMBase_trackDxy[nMMBase] = mmCandidate.trackDxy;
           MMBase_trackIxy[nMMBase] = mmCandidate.trackIxy;
           MMBase_vx[nMMBase] = mmCandidate.vx;
           MMBase_vy[nMMBase] = mmCandidate.vy;
           MMBase_normalizedChi2[nMMBase] = mmCandidate.normalizedChi2;
           MMBase_mass[nMMBase] = mmCandidate.mass;
           MMBase_leadingPt[nMMBase] = mmCandidate.leadingPt;
           MMBase_subleadingPt[nMMBase] = mmCandidate.subleadingPt;
           MMBase_cosAlpha[nMMBase] = mmCandidate.cosAlpha;
           MMBase_dPhi[nMMBase] = mmCandidate.dPhi;
           MMBase_relisoA[nMMBase] = mmCandidate.relisoA;
           MMBase_relisoB[nMMBase] = mmCandidate.relisoB;
           MMBase_fromPVA[nMMBase] = mmCandidate.fromPVA;
           MMBase_fromPVB[nMMBase] = mmCandidate.fromPVB;
           MMBase_PVAssociation[nMMBase] = mmCandidate.PVAssociation;

           if ( fabs(MMBase_trackIxy[nMMBase]) > fabs(MMBase_trackIxy[MMBase_maxIxy]) ) { MMBase_maxIxy = nMMBase; }

        }
        nMMBase++;

     }

   } // end while 


   ////////////////////////////////////////////////////////////////////////////////////////////
   //// ---------------------------------------------------------------------------------- ////
   //// -------------------------------- VERTEX REFITTING -------------------------------- ////                                                    
   //// ---------------------------------------------------------------------------------- ////
   ////////////////////////////////////////////////////////////////////////////////////////////

   // Auxiliary variables
   double ptDif = 0.1;
   double etaDif = 0.05;
   double phiDif = 0.1;

   std::vector<reco::TransientTrack> refit_tracks; // tracks for refitting
   bool excluded = false; // control bool variable to exclude tracks

   RefittedPV_nPFTrack = 0;
   RefittedPV_nLostTrack = 0;
   RefittedPV_nExcludedTrack = 0;


   // No suceed PV refitting default values:
   RefittedPV_vx = -99;
   RefittedPV_vy = -99;
   RefittedPV_vz = -99;

   for (int i = 0; i < nEEBase; i++) { EEBase_refittedDxy[i] = -99; EEBase_refittedIxy[i] = -99; }
   for (int i = 0; i < nMMBase; i++) { MMBase_refittedDxy[i] = -99; MMBase_refittedIxy[i] = -99; }


   // --> Only refit the vertex if running in BSMode:
   if (_BSMode) {

      // Loop over packedPFCandidates:
      for (size_t i = 0; i < packedPFCandidates->size(); i++ ){

         const pat::PackedCandidate & packedPFCandidate = (*packedPFCandidates)[i];

         if (!packedPFCandidate.hasTrackDetails()){ continue;}
         const reco::Track packedPFTrack = packedPFCandidate.pseudoTrack();

         if (abs(packedPFCandidate.pdgId()) == 11){
            if (packedPFCandidate.pvAssociationQuality() < 5 && packedPFCandidate.vertexRef().key() == 0){ continue; }
         } else {
            if (packedPFCandidate.pvAssociationQuality() < 6 && packedPFCandidate.vertexRef().key() == 0){ continue; }
         }
         //if (packedPFCandidate.fromPV(0) != 3){ continue;}

         excluded = false;
         for (size_t j = 0; j < leptonTracks.size(); j++) {
            pat::IsolatedTrack ltrack = leptonTracks.at(j);
            if (fabs(packedPFCandidate.pt() - ltrack.pt()) < ptDif && fabs(packedPFCandidate.eta() - ltrack.eta()) < etaDif && fabs(packedPFCandidate.phi() - ltrack.phi()) < phiDif) { excluded = true; break; }
         }

         if (excluded) {RefittedPV_nExcludedTrack++; continue;}

         RefittedPV_nPFTrack++;
         reco::TransientTrack transientTrack = theTransientTrackBuilder->build(packedPFTrack);
         transientTrack.setBeamSpot(beamSpotObject);
         refit_tracks.push_back(transientTrack);

      }

      // Loop over lostTracks:
      for (size_t i = 0; i < lostTracks->size(); i++){

         const pat::PackedCandidate &lostTrack = (*lostTracks)[i];

         if (!lostTrack.hasTrackDetails()) continue;
         const reco::Track packedLostTrack = lostTrack.pseudoTrack();

         if (abs(lostTrack.pdgId()) == 11){
            if (lostTrack.pvAssociationQuality() < 5 && lostTrack.vertexRef().key() == 0){ continue; }
         } else {
            if (lostTrack.pvAssociationQuality() < 6 && lostTrack.vertexRef().key() == 0){ continue; }
         }
         //if (lostTrack.fromPV(0) != 3){ continue;}

         excluded = false;
         for (size_t j = 0; j < leptonTracks.size(); j++) {
            pat::IsolatedTrack ltrack = leptonTracks.at(j);
            if ( fabs(lostTrack.pt() - ltrack.pt()) < ptDif && fabs(lostTrack.eta() - ltrack.eta()) < etaDif && fabs(lostTrack.phi() - ltrack.phi()) < phiDif) { excluded = true; break; }
         }

         if (excluded) {RefittedPV_nExcludedTrack++; continue;}

         RefittedPV_nLostTrack++;
         reco::TransientTrack  transientTrack = theTransientTrackBuilder->build(packedLostTrack);
         transientTrack.setBeamSpot(beamSpotObject);
         refit_tracks.push_back(transientTrack);

      }


      // -> Primary vertex refitting:

      if (refit_tracks.size() > 1){

         AdaptiveVertexFitter  theFitter(GeometricAnnealing(2.5));
         TransientVertex myVertex = theFitter.vertex(refit_tracks);

         if (myVertex.isValid()){
            RefittedPV_vx = myVertex.position().x();
            RefittedPV_vy = myVertex.position().y();
            RefittedPV_vz = myVertex.position().z();
         

            const reco::Vertex rePV = myVertex;


            // --> Recompute dxy with refitted vertex possition

            // Electrons:
            for (int i = 0; i < nEEBase; i++){

               if (RefittedPV_nExcludedTrack > 0) {

                  const pat::IsolatedTrack &it_A = (*isotracks)[iT.at(ElectronCandidate_isotrackIdx[EEBase_idxA[i]])];
                  const pat::IsolatedTrack &it_B = (*isotracks)[iT.at(ElectronCandidate_isotrackIdx[EEBase_idxB[i]])];
                  const pat::PackedCandidateRef &pckA = it_A.packedCandRef();
                  const pat::PackedCandidateRef &pckB = it_B.packedCandRef();

                  double redxyA = (*pckA).dxy(rePV.position());
                  double redxyB = (*pckB).dxy(rePV.position());
                  double reIxyA = fabs(redxyA/computeDxyError(it_A, rePV));
                  double reIxyB = fabs(redxyB/computeDxyError(it_B, rePV));

                  EEBase_refittedDxy[i] = (reIxyA < reIxyB) ? redxyA : redxyB;
                  EEBase_refittedIxy[i] = (reIxyA < reIxyB) ? reIxyA : reIxyB;


               } else {

                  EEBase_refittedDxy[i] = EEBase_trackDxy[i];
                  EEBase_refittedIxy[i] = EEBase_trackIxy[i];

               }

            } 
            // Muons:
            for (int i = 0; i < nMMBase; i++){

               if (RefittedPV_nExcludedTrack > 0) {

                  const pat::IsolatedTrack &it_A = (*isotracks)[iT.at(MuonCandidate_isotrackIdx[MMBase_idxA[i]])];
                  const pat::IsolatedTrack &it_B = (*isotracks)[iT.at(MuonCandidate_isotrackIdx[MMBase_idxB[i]])];
                  const pat::PackedCandidateRef &pckA = it_A.packedCandRef();
                  const pat::PackedCandidateRef &pckB = it_B.packedCandRef();

                  double redxyA = (*pckA).dxy(rePV.position());
                  double redxyB = (*pckB).dxy(rePV.position());
                  double reIxyA = fabs(redxyA/computeDxyError(it_A, rePV));
                  double reIxyB = fabs(redxyB/computeDxyError(it_B, rePV));

                  MMBase_refittedDxy[i] = (reIxyA < reIxyB) ? redxyA : redxyB;
                  MMBase_refittedIxy[i] = (reIxyA < reIxyB) ? reIxyA : reIxyB;

               } else {

                  MMBase_refittedDxy[i] = MMBase_trackDxy[i];
                  MMBase_refittedIxy[i] = MMBase_trackIxy[i];

               }

            }

         } // end myVertex.isValid() 

      } // end refit_tracks.size() > 1

   } // end _BSMode


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
    //file_out->cd();

    
    // Output Tree definition
    tree_out = new TTree("Events", "Events");

    // Analyzer parameters
    _isData = parameters.getParameter<bool>("isData");
    _BSMode = parameters.getParameter<bool>("BSMode");


    ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Event_event", &Event_event, "Event_event/I");
    tree_out->Branch("Event_run", &Event_run, "Event_run/I");
    tree_out->Branch("Event_luminosityBlock", &Event_luminosityBlock, "Event_luminosityBlock/I");
    
    tree_out->Branch("nPU", &nPU, "nPU/I");
    tree_out->Branch("nPUTrue", &nPUTrue, "nPUTrue/I");
    tree_out->Branch("wPU", &wPU, "wPU/F");
    tree_out->Branch("genWeight", &genWeight, "genWeight/F");

    ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10", &Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10, "Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10/O");
    tree_out->Branch("Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15", &Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15, "Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15/O");

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
    tree_out->Branch("PV_xyFromBS", &PV_xyFromBS, "PV_xyFromBS/F");
    tree_out->Branch("PV_zFromBS", &PV_zFromBS, "PV_zFromBS/F");


    /////////////////////////// REFITTED PRIMARY VERTEX BRANCHES ////////////////////////

    tree_out->Branch("RefittedPV_vx", &RefittedPV_vx, "RefittedPV_vx/F");
    tree_out->Branch("RefittedPV_vy", &RefittedPV_vy, "RefittedPV_vy/F");
    tree_out->Branch("RefittedPV_vz", &RefittedPV_vz, "RefittedPV_vz/F");
    tree_out->Branch("RefittedPV_nPFTrack", &RefittedPV_nPFTrack, "RefittedPV_nPFTrack/I");
    tree_out->Branch("RefittedPV_nLostTrack", &RefittedPV_nLostTrack, "RefittedPV_nLostTrack/I");
    tree_out->Branch("RefittedPV_nExcludedTrack", &RefittedPV_nExcludedTrack, "RefittedPV_nExcludedTrack/I");


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
    /*
    tree_out->Branch("IsoTrackSel_numberOfValidTrackerHits", IsoTrackSel_numberOfValidTrackerHits, "IsoTrackSel_numberOfValidTrackerHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidPixelHits", IsoTrackSel_numberOfValidPixelHits, "IsoTrackSel_numberOfValidPixelHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidPixelBarrelHits", IsoTrackSel_numberOfValidPixelBarrelHits, "IsoTrackSel_numberOfValidPixelBarrelHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidPixelEndcapHits", IsoTrackSel_numberOfValidPixelEndcapHits, "IsoTrackSel_numberOfValidPixelEndcapHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripHits", IsoTrackSel_numberOfValidStripHits, "IsoTrackSel_numberOfValidStripHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripTIBHits", IsoTrackSel_numberOfValidStripTIBHits, "IsoTrackSel_numberOfValidStripTIBHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripTIDHits", IsoTrackSel_numberOfValidStripTIDHits, "IsoTrackSel_numberOfValidStripTIDHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripTOBHits", IsoTrackSel_numberOfValidStripTOBHits, "IsoTrackSel_numberOfValidStripTOBHits[nIsoTrack]/I");
    tree_out->Branch("IsoTrackSel_numberOfValidStripTECHits", IsoTrackSel_numberOfValidStripTECHits, "IsoTrackSel_numberOfValidStripTECHits[nIsoTrack]/I");
    */

    
    ////////////////////////////////// PHOTON BRANCHES //////////////////////////////////
    
    tree_out->Branch("nPhoton", &nPhoton, "nPhoton/I");
    tree_out->Branch("PhotonSel_et", PhotonSel_et, "PhotonSel_et[nPhoton]/F");
    tree_out->Branch("PhotonSel_eta", PhotonSel_eta, "PhotonSel_eta[nPhoton]/F");
    tree_out->Branch("PhotonSel_phi", PhotonSel_phi, "PhotonSel_phi[nPhoton]/F");
    /*
    tree_out->Branch("PhotonSel_hadronicOverEm", PhotonSel_hadronicOverEm, "PhotonSel_hadronicOverEm[nPhoton]/F");
    tree_out->Branch("PhotonSel_full5x5_sigmaIetaIeta", PhotonSel_full5x5_sigmaIetaIeta, "PhotonSel_full5x5_sigmaIetaIeta[nPhoton]/F");
    tree_out->Branch("PhotonSel_isEB", PhotonSel_isEB, "PhotonSel_isEB[nPhoton]/I");
    tree_out->Branch("PhotonSel_isEE", PhotonSel_isEE, "PhotonSel_isEE[nPhoton]/I");
    tree_out->Branch("PhotonSel_r9", PhotonSel_r9, "PhotonSel_r9[nPhoton]/F");
    tree_out->Branch("PhotonSel_ecalIso", PhotonSel_ecalIso, "PhotonSel_ecalIso[nPhoton]/F");
    tree_out->Branch("PhotonSel_hcalIso", PhotonSel_hcalIso, "PhotonSel_hcalIso[nPhoton]/F");
    tree_out->Branch("PhotonSel_caloIso", PhotonSel_caloIso, "PhotonSel_caloIso[nPhoton]/F");
    tree_out->Branch("PhotonSel_relIso", PhotonSel_relIso, "PhotonSel_relIso[nPhoton]/F");
    */

    ///////////////////////////////// ELECTRON BRANCHES /////////////////////////////////

    /*
    tree_out->Branch("nElectron", &nElectron, "nElectron/I");
    tree_out->Branch("ElectronSel_pt", ElectronSel_pt, "ElectronSel_pt[nElectron]/F");
    tree_out->Branch("ElectronSel_et", ElectronSel_et, "ElectronSel_et[nElectron]/F");
    tree_out->Branch("ElectronSel_eta", ElectronSel_eta, "ElectronSel_eta[nElectron]/F");
    tree_out->Branch("ElectronSel_phi", ElectronSel_phi, "ElectronSel_phi[nElectron]/F");
    tree_out->Branch("ElectronSel_hadronicOverEm", ElectronSel_hadronicOverEm, "ElectronSel_hadronicOverEm[nElectron]/F");
    tree_out->Branch("ElectronSel_full5x5_sigmaIetaIeta", ElectronSel_full5x5_sigmaIetaIeta, "ElectronSel_full5x5_sigmaIetaIeta[nElectron]/F");
    tree_out->Branch("ElectronSel_isEB", ElectronSel_isEB, "ElectronSel_isEB[nElectron]/I");
    tree_out->Branch("ElectronSel_isEE", ElectronSel_isEE, "ElectronSel_isEE[nElectron]/I");
    tree_out->Branch("ElectronSel_r9", ElectronSel_r9, "ElectronSel_r9[nElectron]/F");
    tree_out->Branch("ElectronSel_trackIso", ElectronSel_trackIso, "ElectronSel_trackIso[nElectron]/F");
    tree_out->Branch("ElectronSel_ecalIso", ElectronSel_ecalIso, "ElectronSel_ecalIso[nElectron]/F");
    tree_out->Branch("ElectronSel_hcalIso", ElectronSel_hcalIso, "ElectronSel_hcalIso[nElectron]/F");
    tree_out->Branch("ElectronSel_caloIso", ElectronSel_caloIso, "ElectronSel_caloIso[nElectron]/F");
    tree_out->Branch("ElectronSel_relIso", ElectronSel_relIso, "ElectronSel_relIso[nElectron]/F");
    tree_out->Branch("ElectronSel_dxy", ElectronSel_dxy, "ElectronSel_dxy[nElectron]/F");
    tree_out->Branch("ElectronSel_dxyError", ElectronSel_dxyError, "ElectronSel_dxyError[nElectron]/F");
    tree_out->Branch("ElectronSel_dxySignificance", ElectronSel_dxySignificance, "ElectronSel_dxySignificance[nElectron]/F");
    tree_out->Branch("ElectronSel_dB", ElectronSel_dB, "ElectronSel_dB[nElectron]/F");
    tree_out->Branch("ElectronSel_edB", ElectronSel_edB, "ElectronSel_edB[nElectron]/F");
    tree_out->Branch("ElectronSel_isLoose", ElectronSel_isLoose, "ElectronSel_isLoose[nElectron]/I");
    */

    ///////////////////////////////// MUON BRANCHES /////////////////////////////////

    tree_out->Branch("nMuon", &nMuon, "nMuon/I");
    tree_out->Branch("MuonSel_pt", MuonSel_pt, "MuonSel_pt[nMuon]/F");
    tree_out->Branch("MuonSel_eta", MuonSel_eta, "MuonSel_eta[nMuon]/F");
    tree_out->Branch("MuonSel_phi", MuonSel_phi, "MuonSel_phi[nMuon]/F");
    tree_out->Branch("MuonSel_relIso", MuonSel_relIso, "MuonSel_relIso[nMuon]/F");
    tree_out->Branch("MuonSel_isMuon", MuonSel_isMuon, "MuonSel_isMuon[nMuon]/I");
    tree_out->Branch("MuonSel_isGlobalMuon", MuonSel_isGlobalMuon, "MuonSel_isGlobalMuon[nMuon]/I");
    tree_out->Branch("MuonSel_isTrackerMuon", MuonSel_isTrackerMuon, "MuonSel_isTrackerMuon[nMuon]/I");
    tree_out->Branch("MuonSel_isStandAloneMuon", MuonSel_isStandAloneMuon, "MuonSel_isStandAloneMuon[nMuon]/I");
    tree_out->Branch("MuonSel_isLooseMuon", MuonSel_isLooseMuon, "MuonSel_isLooseMuon[nMuon]/I");
    tree_out->Branch("MuonSel_isMediumMuon", MuonSel_isMediumMuon, "MuonSel_isMediumMuon[nMuon]/I");
    tree_out->Branch("MuonSel_isGoodMediumMuon", MuonSel_isGoodMediumMuon, "MuonSel_isGoodMediumMuon[nMuon]/I");

    tree_out->Branch("MuonSel_dB", MuonSel_dB, "MuonSel_dB[nMuon]/F");
    tree_out->Branch("MuonSel_edB", MuonSel_edB, "MuonSel_edB[nMuon]/F");
    tree_out->Branch("MuonSel_dBSignificance", MuonSel_dBSignificance, "MuonSel_dBSignificance[nMuon]/F");

    /*
    tree_out->Branch("MuonSel_isPFMuon", MuonSel_isPFMuon, "MuonSel_isPFMuon[nMuon]/I");
    tree_out->Branch("MuonSel_fractionOfValidTrackerHits", MuonSel_fractionOfValidTrackerHits, "MuonSel_fractionOfValidTrackerHits[nMuon]/F");
    tree_out->Branch("MuonSel_normGlobalTrackChi2", MuonSel_normGlobalTrackChi2, "MuonSel_normGlobalTrackChi2[nMuon]/F");
    */

    //////////////////////////// MUON TRIGGER OBJECT BRANCHES ///////////////////////////
    //
    tree_out->Branch("nMuonTriggerObject", &nMuonTriggerObject, "nMuonTriggerObject/I");
    tree_out->Branch("MuonTriggerObjectSel_pt", MuonTriggerObjectSel_pt, "MuonTriggerObjectSel_pt[nMuonTriggerObject]/F");
    tree_out->Branch("MuonTriggerObjectSel_eta", MuonTriggerObjectSel_eta, "MuonTriggerObjectSel_eta[nMuonTriggerObject]/F");
    tree_out->Branch("MuonTriggerObjectSel_phi", MuonTriggerObjectSel_phi, "MuonTriggerObjectSel_phi[nMuonTriggerObject]/F");


    //////////////////////////////// GENPARTICLE BRANCHES ///////////////////////////////

    tree_out->Branch("GenHiggs_pt", &GenHiggs_pt, "GenHiggs_pt/F");
    tree_out->Branch("GenHiggs_eta", &GenHiggs_eta, "GenHiggs_eta/F");
    tree_out->Branch("GenHiggs_phi", &GenHiggs_phi, "GenHiggs_phi/F");

    tree_out->Branch("nGenLepton", &nGenLepton, "nGenLepton/I");
    tree_out->Branch("GenLeptonSel_pt", GenLeptonSel_pt, "GenLeptonSel_pt[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_et", GenLeptonSel_et, "GenLeptonSel_et[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_eta", GenLeptonSel_eta, "GenLeptonSel_eta[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_phi", GenLeptonSel_phi, "GenLeptonSel_phi[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_dxy", GenLeptonSel_dxy, "GenLeptonSel_dxy[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_vx", GenLeptonSel_vx, "GenLeptonSel_vx[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_vy", GenLeptonSel_vy, "GenLeptonSel_vy[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_vz", GenLeptonSel_vz, "GenLeptonSel_vz[nGenLepton]/F");
    tree_out->Branch("GenLeptonSel_pdgId", GenLeptonSel_pdgId, "GenLeptonSel_pdgId[nGenLepton]/I");
    tree_out->Branch("GenLeptonSel_motherIdx", GenLeptonSel_motherIdx, "GenLeptonSel_motherIdx[nGenLepton]/I");

    
    tree_out->Branch("nGenNeutralino", &nGenNeutralino, "nGenNeutralino/I");  
    tree_out->Branch("GenNeutralinoSel_pt", GenNeutralinoSel_pt, "GenNeutralinoSel_pt[nGenNeutralino]/F");
    tree_out->Branch("GenNeutralinoSel_eta", GenNeutralinoSel_eta, "GenNeutralinoSel_eta[nGenNeutralino]/F");
    tree_out->Branch("GenNeutralinoSel_phi", GenNeutralinoSel_phi, "GenNeutralinoSel_phi[nGenNeutralino]/F"); 
    tree_out->Branch("GenNeutralinoSel_Lxy", GenNeutralinoSel_Lxy, "GenNeutralinoSel_Lxy[nGenNeutralino]/F");
    tree_out->Branch("GenNeutralinoSel_pdgId", GenNeutralinoSel_pdgId, "GenNeutralinoSel_pdgId[nGenNeutralino]/I");
    tree_out->Branch("GenNeutralinoSel_decaypdgId", GenNeutralinoSel_decaypdgId, "GenNeutralinoSel_decaypdgId[nGenNeutralino]/I");
    tree_out->Branch("GenNeutralinoSel_passAcceptance", GenNeutralinoSel_passAcceptance, "GenNeutralinoSel_passAcceptance[nGenNeutralino]/O");
    
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
    tree_out->Branch("ElectronCandidate_dxySignificance", ElectronCandidate_dxySignificance, "ElectronCandidate_dxySignificance[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_pvAssociationQuality", ElectronCandidate_pvAssociationQuality, "ElectronCandidate_pvAssociationQuality[nElectronCandidate]/I");
    tree_out->Branch("ElectronCandidate_ptDiff", ElectronCandidate_ptDiff, "ElectronCandidate_ptDiff[nElectronCandidate]/F");

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
    tree_out->Branch("MuonCandidate_dxy", MuonCandidate_dxy, "MuonCandidate_dxy[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_dxyError", MuonCandidate_dxyError, "MuonCandidate_dxyError[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_dxySignificance", MuonCandidate_dxySignificance, "MuonCandidate_dxySignificance[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_pvAssociationQuality", MuonCandidate_pvAssociationQuality, "MuonCandidate_pvAssociationQuality[nMuonCandidate]/I");
    tree_out->Branch("MuonCandidate_ptDiff", MuonCandidate_ptDiff, "MuonCandidate_ptDiff[nMuonCandidate]/F");

    ////////////////////////////// LL BRANCHES /////////////////////////////

    tree_out->Branch("nEE", &nEE, "nEE/I");
    if (!_BSMode) {
       tree_out->Branch("EE_idxA", EE_idxA, "EE_idxA[nEE]/I");
       tree_out->Branch("EE_idxB", EE_idxB, "EE_idxB[nEE]/I");
       tree_out->Branch("EE_Lxy", EE_Lxy, "EE_Lxy[nEE]/F");
       tree_out->Branch("EE_Ixy", EE_Ixy, "EE_Ixy[nEE]/F");
       tree_out->Branch("EE_trackDxy", EE_trackDxy, "EE_trackDxy[nEE]/F");
       tree_out->Branch("EE_trackIxy", EE_trackIxy, "EE_trackIxy[nEE]/F");
       tree_out->Branch("EE_mass", EE_mass, "EE_mass[nEE]/F");
       tree_out->Branch("EE_normalizedChi2", EE_normalizedChi2, "EE_normalizedChi2[nEE]/F");
       tree_out->Branch("EE_leadingPt", EE_leadingPt, "EE_leadingPt[nEE]/F");
       tree_out->Branch("EE_subleadingPt", EE_subleadingPt, "EE_subleadingPt[nEE]/F");
       tree_out->Branch("EE_leadingEt", EE_leadingEt, "EE_leadingEt[nEE]/F");
       tree_out->Branch("EE_subleadingEt", EE_subleadingEt, "EE_subleadingEt[nEE]/F");
       tree_out->Branch("EE_cosAlpha", EE_cosAlpha, "EE_cosAlpha[nEE]/F");
       tree_out->Branch("EE_dPhi", EE_dPhi, "EE_dPhi[nEE]/F");
       tree_out->Branch("EE_relisoA", EE_relisoA, "EE_relisoA[nEE]/F");
       tree_out->Branch("EE_relisoB", EE_relisoB, "EE_relisoB[nEE]/F");
    }

    tree_out->Branch("nEEBase", &nEEBase, "nEEBase/I");
    tree_out->Branch("EEBase_maxIxy", &EEBase_maxIxy, "EEBase_maxIxy/I");
    tree_out->Branch("EEBase_idxA", EEBase_idxA, "EEBase_idxA[nEEBase]/I");
    tree_out->Branch("EEBase_idxB", EEBase_idxB, "EEBase_idxB[nEEBase]/I");
    tree_out->Branch("EEBase_vx", EEBase_vx, "EEBase_vx[nEEBase]/F");
    tree_out->Branch("EEBase_vy", EEBase_vy, "EEBase_vy[nEEBase]/F");
    tree_out->Branch("EEBase_Lxy", EEBase_Lxy, "EEBase_Lxy[nEEBase]/F");
    tree_out->Branch("EEBase_Ixy", EEBase_Ixy, "EEBase_Ixy[nEEBase]/F");
    tree_out->Branch("EEBase_trackDxy", EEBase_trackDxy, "EEBase_trackDxy[nEEBase]/F");
    tree_out->Branch("EEBase_trackIxy", EEBase_trackIxy, "EEBase_trackIxy[nEEBase]/F");
    tree_out->Branch("EEBase_mass", EEBase_mass, "EEBase_mass[nEEBase]/F");
    tree_out->Branch("EEBase_normalizedChi2", EEBase_normalizedChi2, "EEBase_normalizedChi2[nEEBase]/F");
    tree_out->Branch("EEBase_leadingPt", EEBase_leadingPt, "EEBase_leadingPt[nEEBase]/F");
    tree_out->Branch("EEBase_subleadingPt", EEBase_subleadingPt, "EEBase_subleadingPt[nEEBase]/F");
    tree_out->Branch("EEBase_leadingEt", EEBase_leadingEt, "EEBase_leadingEt[nEEBase]/F");
    tree_out->Branch("EEBase_subleadingEt", EEBase_subleadingEt, "EEBase_subleadingEt[nEEBase]/F");
    tree_out->Branch("EEBase_cosAlpha", EEBase_cosAlpha, "EEBase_cosAlpha[nEEBase]/F");
    tree_out->Branch("EEBase_dPhi", EEBase_dPhi, "EEBase_dPhi[nEEBase]/F");
    tree_out->Branch("EEBase_relisoA", EEBase_relisoA, "EEBase_relisoA[nEEBase]/F");
    tree_out->Branch("EEBase_relisoB", EEBase_relisoB, "EEBase_relisoB[nEEBase]/F");
    tree_out->Branch("EEBase_refittedDxy", EEBase_refittedDxy, "EEBase_refittedDxy[nEEBase]/F");
    tree_out->Branch("EEBase_refittedIxy", EEBase_refittedIxy, "EEBase_refittedIxy[nEEBase]/F");
    tree_out->Branch("EEBase_fromPVA", EEBase_fromPVA, "EEBase_fromPVA[nEEBase]/I");
    tree_out->Branch("EEBase_fromPVB", EEBase_fromPVB, "EEBase_fromPVB[nEEBase]/I");
    tree_out->Branch("EEBase_PVAssociation", EEBase_PVAssociation, "EEBase_PVAssociation[nEEBase]/I");

    tree_out->Branch("nMM", &nMM, "nMM/I");
    if (!_BSMode) {
       tree_out->Branch("MM_idxA", MM_idxA, "MM_idxA[nMM]/I");
       tree_out->Branch("MM_idxB", MM_idxB, "MM_idxB[nMM]/I");
       tree_out->Branch("MM_Lxy", MM_Lxy, "MM_Lxy[nMM]/F");
       tree_out->Branch("MM_Ixy", MM_Ixy, "MM_Ixy[nMM]/F");
       tree_out->Branch("MM_trackDxy", MM_trackDxy, "MM_trackDxy[nMM]/F");
       tree_out->Branch("MM_trackIxy", MM_trackIxy, "MM_trackIxy[nMM]/F");
       tree_out->Branch("MM_mass", MM_mass, "MM_mass[nMM]/F");
       tree_out->Branch("MM_normalizedChi2", MM_normalizedChi2, "MM_normalizedChi2[nMM]/F");
       tree_out->Branch("MM_leadingPt", MM_leadingPt, "MM_leadingPt[nMM]/F");
       tree_out->Branch("MM_subleadingPt", MM_subleadingPt, "MM_subleadingPt[nMM]/F");
       tree_out->Branch("MM_cosAlpha", MM_cosAlpha, "MM_cosAlpha[nMM]/F");
       tree_out->Branch("MM_dPhi", MM_dPhi, "MM_dPhi[nMM]/F");
       tree_out->Branch("MM_relisoA", MM_relisoA, "MM_relisoA[nMM]/F");
       tree_out->Branch("MM_relisoB", MM_relisoB, "MM_relisoB[nMM]/F");
    }

    tree_out->Branch("nMMBase", &nMMBase, "nMMBase/I");
    tree_out->Branch("MMBase_maxIxy", &MMBase_maxIxy, "MMBase_maxIxy/I");
    tree_out->Branch("MMBase_idxA", MMBase_idxA, "MMBase_idxA[nMMBase]/I");
    tree_out->Branch("MMBase_idxB", MMBase_idxB, "MMBase_idxB[nMMBase]/I");
    tree_out->Branch("MMBase_vx", MMBase_vx, "MMBase_vx[nMMBase]/F");
    tree_out->Branch("MMBase_vy", MMBase_vy, "MMBase_vy[nMMBase]/F");
    tree_out->Branch("MMBase_Lxy", MMBase_Lxy, "MMBase_Lxy[nMMBase]/F");
    tree_out->Branch("MMBase_Ixy", MMBase_Ixy, "MMBase_Ixy[nMMBase]/F");
    tree_out->Branch("MMBase_trackDxy", MMBase_trackDxy, "MMBase_trackDxy[nMMBase]/F");
    tree_out->Branch("MMBase_trackIxy", MMBase_trackIxy, "MMBase_trackIxy[nMMBase]/F");
    tree_out->Branch("MMBase_mass", MMBase_mass, "MMBase_mass[nMMBase]/F");
    tree_out->Branch("MMBase_normalizedChi2", MMBase_normalizedChi2, "MMBase_normalizedChi2[nMMBase]/F");
    tree_out->Branch("MMBase_leadingPt", MMBase_leadingPt, "MMBase_leadingPt[nMMBase]/F");
    tree_out->Branch("MMBase_subleadingPt", MMBase_subleadingPt, "MMBase_subleadingPt[nMMBase]/F");
    tree_out->Branch("MMBase_cosAlpha", MMBase_cosAlpha, "MMBase_cosAlpha[nMMBase]/F");
    tree_out->Branch("MMBase_dPhi", MMBase_dPhi, "MMBase_dPhi[nMMBase]/F");
    tree_out->Branch("MMBase_relisoA", MMBase_relisoA, "MMBase_relisoA[nMMBase]/F");
    tree_out->Branch("MMBase_relisoB", MMBase_relisoB, "MMBase_relisoB[nMMBase]/F");
    tree_out->Branch("MMBase_refittedDxy", MMBase_refittedDxy, "MMBase_refittedDxy[nMMBase]/F");
    tree_out->Branch("MMBase_refittedIxy", MMBase_refittedIxy, "MMBase_refittedIxy[nMMBase]/F");
    tree_out->Branch("MMBase_fromPVA", MMBase_fromPVA, "MMBase_fromPVA[nMMBase]/I");
    tree_out->Branch("MMBase_fromPVB", MMBase_fromPVB, "MMBase_fromPVB[nMMBase]/I");
    tree_out->Branch("MMBase_PVAssociation", MMBase_PVAssociation, "MMBase_PVAssociation[nMMBase]/I");

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

bool isLooseElectron(const pat::Electron & electron) {

    float ecal_energy_inverse = 1.0/electron.ecalEnergy();
    float eSCoverP = electron.eSuperClusterOverP();

    // Barrel cuts:
    if (fabs(electron.superCluster()->eta()) <= 1.479) {

       // Combined isolation:
       //float comIso = (electron.dr03TkSumPt() + max(0., electron.dr03EcalRecHitSumEt() - 1.) + electron.dr03HcalTowerSumEt() ) / electron.pt()

       if (electron.full5x5_sigmaIetaIeta() > 0.011) { return false; }
       if (fabs(electron.deltaEtaSeedClusterTrackAtVtx()) > 0.00477) { return false; }
       if (fabs(electron.deltaPhiSuperClusterTrackAtVtx()) > 0.222) { return false; }
       if (electron.hadronicOverEm() > 0.298) { return false; }
       if (fabs(1.0 - eSCoverP)*ecal_energy_inverse > 0.241) {return false; }
       if (electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) > 1) { return false;}
       if (electron.passConversionVeto() == 0) {return false;}


    // Endcap cuts
    } else if (fabs(electron.superCluster()->eta()) > 1.479) {


       if (electron.full5x5_sigmaIetaIeta() > 0.0314) { return false; }
       if (fabs(electron.deltaEtaSeedClusterTrackAtVtx()) > 0.00868) { return false; }
       if (fabs(electron.deltaPhiSuperClusterTrackAtVtx()) > 0.213) { return false; }
       if (electron.hadronicOverEm() > 0.101) { return false; }
       if (fabs(1.0 - eSCoverP)*ecal_energy_inverse > 0.14) {return false; }
       if (electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) > 1) { return false;}
       if (electron.passConversionVeto() == 0) {return false;}

    }

    return true;

}

//=======================================================================================================================================================================================================================//

bool isMediumElectron(const pat::Electron & electron) {

    // Need to be filled correctly (Now they are the requirements for Loose Electrons)

    float ecal_energy_inverse = 1.0/electron.ecalEnergy();
    float eSCoverP = electron.eSuperClusterOverP();

    // Barrel cuts:
    if (fabs(electron.superCluster()->eta()) <= 1.479) {

       // Combined isolation:
       //float comIso = (electron.dr03TkSumPt() + max(0., electron.dr03EcalRecHitSumEt() - 1.) + electron.dr03HcalTowerSumEt() ) / electron.pt()

       if (electron.full5x5_sigmaIetaIeta() > 0.011) { return false; }
       if (fabs(electron.deltaEtaSeedClusterTrackAtVtx()) > 0.00477) { return false; }
       if (fabs(electron.deltaPhiSuperClusterTrackAtVtx()) > 0.222) { return false; }
       if (electron.hadronicOverEm() > 0.298) { return false; }
       if (fabs(1.0 - eSCoverP)*ecal_energy_inverse > 0.241) {return false; }
       if (electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) > 1) { return false;}
       if (electron.passConversionVeto() == 0) {return false;}


    // Endcap cuts
    } else if (fabs(electron.superCluster()->eta()) > 1.479) {


       if (electron.full5x5_sigmaIetaIeta() > 0.0314) { return false; }
       if (fabs(electron.deltaEtaSeedClusterTrackAtVtx()) > 0.00868) { return false; }
       if (fabs(electron.deltaPhiSuperClusterTrackAtVtx()) > 0.213) { return false; }
       if (electron.hadronicOverEm() > 0.101) { return false; }
       if (fabs(1.0 - eSCoverP)*ecal_energy_inverse > 0.14) {return false; }
       if (electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) > 1) { return false;}
       if (electron.passConversionVeto() == 0) {return false;}

    }

    return true;

}
//=======================================================================================================================================================================================================================//

bool LongLivedAnalysis::passIsotrackSelection( const pat::IsolatedTrack &track) {

   // Quality cuts:
   const reco::HitPattern &hits = track.hitPattern();
   if (hits.numberOfValidTrackerHits() < 6) { return false; }
   if (!track.isHighPurityTrack()) { return false;}

   // Isotrack must have packed candidate:
   const pat::PackedCandidateRef &pckCand = track.packedCandRef();
   if (!pckCand.isNonnull()) { return false; }

   // Preselection cuts:
   if (track.pt() < 28) { return false; }
   if (fabs(track.eta()) > 2) { return false; }

   // To be noticed: Isolation cuts are applied later with the LLCandidate selection

   return true;
}


//=======================================================================================================================================================================================================================//

bool LongLivedAnalysis::passPhotonSelection( const pat::Photon &photon ) {

   // Quality cuts:
   if (photon.hadronicOverEm() > 0.05) { return false; } 
   if (photon.isEE() && photon.full5x5_sigmaIetaIeta() > 0.034) { return false; }
   if (photon.isEB() && photon.full5x5_sigmaIetaIeta() > 0.012) { return false; }

   // Preselection cuts:
   if (fabs(photon.eta()) > 1.4442) { return false; }
   if (photon.et() < 28) {return false; }

   return true;

}

//=======================================================================================================================================================================================================================//

bool LongLivedAnalysis::passL2MuonSelection( pat::TriggerObjectStandAlone obj) {

   if (fabs(obj.eta()) > 2) { return false; }
   return true;
}

//=======================================================================================================================================================================================================================//

bool LongLivedAnalysis::passMuonSelection(const pat::Muon &muon) {

   if (muon.pt() < 31){ return false; }
   if (fabs(muon.eta()) > 2) { return false; }
   return true;
}

//=======================================================================================================================================================================================================================//

bool LongLivedAnalysis::passBaselineSelection(llCandidate llc) {

   // Electron selection:
   if (llc.type == 0) {

      if ( llc.leadingPt < 45 ) { return false; }
      if ( llc.subleadingPt < 28 ) { return false; }
      if ( llc.leadingEt < 45 ) { return false; }
      if ( llc.subleadingEt < 28 ) { return false; }
      if ( fabs(llc.etaA) > 1.442 || fabs(llc.etaB) > 1.442 ) { return false; }
      if ( fabs(llc.relisoA) > 0.2 || fabs(llc.relisoB) > 0.2 ) { return false; }
      if ( llc.normalizedChi2 > 10 ) { return false; }
      if ( llc.mass < 15 ) { return false; }

      return true;

   }

   // Muon selection:
   if (llc.type == 1) {

      if ( llc.leadingPt < 31 ) { return false; }
      if ( llc.subleadingPt < 31 ) { return false; }
      if ( fabs(llc.etaA) > 2 || fabs(llc.etaB) > 2 ) { return false; }
      if ( fabs(llc.relisoA) > 0.1 || fabs(llc.relisoB) > 0.1 ) { return false; }
      if ( llc.normalizedChi2 > 5 ) { return false; }
      if ( llc.mass < 15 ) { return false; }
      if ( llc.dR < 0.2 ) { return false; }
      if ( llc.cosAlpha < -0.79 ) { return false; }

      return true;

   }

   // Warning if not electron or muon.

   std::cout << "Warning: The selected llCandidate is not an electron or muon" << std::endl;
   return false;

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
//

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

      const reco::Track isorecotrkA = pckCandA->pseudoTrack();
      const reco::Track isorecotrkB = pckCandB->pseudoTrack();
      double trackDxy = (fabs(it_A.dxy()/it_A.dxyError()) < fabs(it_B.dxy())/it_B.dxyError())? it_A.dxy(): it_B.dxy();
      double trackIxy = (fabs(it_A.dxy()/it_A.dxyError()) < fabs(it_B.dxy())/it_B.dxyError())? it_A.dxy()/it_A.dxyError(): it_B.dxy()/it_B.dxyError();
      double leadingPt = (isorecotrkA.pt()>isorecotrkB.pt())? isorecotrkA.pt(): isorecotrkB.pt();
      double subleadingPt = (isorecotrkA.pt()<isorecotrkB.pt())? isorecotrkA.pt(): isorecotrkB.pt();



      // cosAlha and dPhi global computation:
      TVector3 vec3A(isorecotrkA.px(), isorecotrkA.py(), isorecotrkA.pz()); 
      TVector3 vec3B(isorecotrkB.px(), isorecotrkB.py(), isorecotrkB.pz()); 
      TVector3 divec3 = vec3A + vec3B;
      TVector3 vtxvec3(secV.x() - PV_vx, secV.y() - PV_vy, secV.z() - PV_vz);
      double cosAlpha = TMath::Cos(vec3A.Angle(vec3B));
      double dPhi = divec3.DeltaPhi(vtxvec3);

      // Isolation of both candidates
      const pat::PFIsolation &pfisoA = it_A.pfIsolationDR03();
      double relisoA = (fmax(0.0, pfisoA.photonIso() + pfisoA.neutralHadronIso() - 0.5*pfisoA.puChargedHadronIso()) + pfisoA.chargedHadronIso())/it_A.pt();
      const pat::PFIsolation &pfisoB = it_B.pfIsolationDR03();
      double relisoB = (fmax(0.0, pfisoB.photonIso() + pfisoB.neutralHadronIso() - 0.5*pfisoB.puChargedHadronIso()) + pfisoB.chargedHadronIso())/it_B.pt();

      // ########################## LLSel choice ################################

      // ####### Fill the dilepton candidates information

      if (isEE) {
        EE_idxA[nEE] = idxA;
        EE_idxB[nEE] = idxB;
	EE_Lxy[nEE] = vMeas.value();
	EE_Ixy[nEE] = vMeas.significance();
	EE_trackDxy[nEE] = trackDxy;
	EE_trackIxy[nEE] = trackIxy;
        EE_normalizedChi2[nEE] = myVertex.normalisedChiSquared();
	//EE_mass[nEE] = secVkin.weightedVectorSum().M(); 
        EE_leadingPt[nEE] = leadingPt;
        EE_subleadingPt[nEE] = subleadingPt;
        EE_leadingEt[nEE] = (ElectronCandidate_et[idxA] > ElectronCandidate_et[idxB])? ElectronCandidate_et[idxA]: ElectronCandidate_et[idxB];
        EE_subleadingEt[nEE] = (ElectronCandidate_et[idxA] < ElectronCandidate_et[idxB])? ElectronCandidate_et[idxA]: ElectronCandidate_et[idxB];
        EE_cosAlpha[nEE] = cosAlpha;
        EE_dPhi[nEE] = dPhi;
        EE_relisoA[nEE] = relisoA;
        EE_relisoB[nEE] = relisoB;

        // Invariant mass computed for electrons
        TLorentzVector TLA; 
        TLA.SetPtEtaPhiM(isorecotrkA.pt(), isorecotrkA.eta(), isorecotrkA.phi(), 0.510/1000.0);
        TLorentzVector TLB; 
        TLB.SetPtEtaPhiM(isorecotrkB.pt(), isorecotrkB.eta(), isorecotrkB.phi(), 0.510/1000.0);
        EE_mass[nEE] = (TLA + TLB).M();


        
	nEE++;
      }
      else {
        MM_idxA[nMM] = idxA;
        MM_idxB[nMM] = idxB;
	MM_Lxy[nMM] = vMeas.value();
	MM_Ixy[nMM] = vMeas.significance();
	MM_trackDxy[nMM] = trackDxy;
	MM_trackIxy[nMM] = trackIxy;
        MM_normalizedChi2[nMM] = myVertex.normalisedChiSquared();
	//MM_mass[nMM] = secVkin.weightedVectorSum().M(); 
        MM_leadingPt[nMM] = leadingPt;
        MM_subleadingPt[nMM] = subleadingPt;
        MM_cosAlpha[nMM] = cosAlpha;
        MM_dPhi[nMM] = dPhi;
        MM_relisoA[nMM] = relisoA;
        MM_relisoB[nMM] = relisoB;

        // Invariant mass computed for muons
        TLorentzVector TLA; 
        TLA.SetPtEtaPhiM(isorecotrkA.pt(), isorecotrkA.eta(), isorecotrkA.phi(), 105.658/1000.0);
        TLorentzVector TLB; 
        TLB.SetPtEtaPhiM(isorecotrkB.pt(), isorecotrkB.eta(), isorecotrkB.phi(), 105.658/1000.0);
        MM_mass[nMM] = (TLA + TLB).M();
 

	nMM++;
      }
    }
    else return false;
  } 
  else return false;

  return true;
  
}




DEFINE_FWK_MODULE(LongLivedAnalysis);
