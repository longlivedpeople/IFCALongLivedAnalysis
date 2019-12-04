#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <TROOT.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>


/////////////////////////////////////////////////////////
// --------------------------------------------------- //
// -------- llCandidate DataFormat definition -------- //
// --------------------------------------------------- //
/////////////////////////////////////////////////////////

struct llCandidate
{

   // ---- llCandidate DataFormat information

   int type = 99;  // 0:electrons 1:muons
   bool canFitVertex = false;
   bool hasValidVertex = false;
   double vertexLxy = 0;
   double vertexIxy = 0;
   double normalizedChi2 = 0;
   double trackDxy = 0;
   double trackIxy = 0;
   double etaA = 0;
   double etaB = 0;
   double leadingPt = 0;
   double subleadingPt = 0;
   double mass = 0;   
   double ptll = 0;
   double cosAlpha = 0;
   double dPhi = 0;
   double dR = 0;
   double relisoA = 0;
   double relisoB = 0;

   // Additionally for electrons ( type == 0 ) to be filled manually:
   double leadingEt = 0;
   double subleadingEt = 0;


   // ---- llCandidate constructors: 
   llCandidate(const reco::Vertex &pv, edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder, const pat::IsolatedTrack &it_A, const pat::IsolatedTrack &it_B, bool isEE)
   { 
      Init(theTransientTrackBuilder, pv, it_A, it_B, isEE); 
   };
 
   ~llCandidate(){};


   // ---- llCandidate DataFormat functions

   // -- Init llCandidate:
   void Init(edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder, const reco::Vertex &pv, const pat::IsolatedTrack &it_A, const pat::IsolatedTrack &it_B, bool isEE)
   {

      // Identify type:
      if (isEE) { type = 0; }
      else { type = 1; }

      const pat::PackedCandidateRef &pckCand_A = it_A.packedCandRef();
      const pat::PackedCandidateRef &pckCand_B = it_B.packedCandRef();

      if (pckCand_A.isNonnull() && pckCand_B.isNonnull() && pckCand_A->hasTrackDetails() && pckCand_B->hasTrackDetails()) {

         canFitVertex = true;

         // Get tracks:
         std::vector<reco::TransientTrack> vec_refitTracks;
         const reco::Track isorecotrkA = pckCand_A->pseudoTrack();
         const reco::Track isorecotrkB = pckCand_B->pseudoTrack();
         reco::TransientTrack isotransienttrackA = theTransientTrackBuilder->build(isorecotrkA);
         reco::TransientTrack isotransienttrackB = theTransientTrackBuilder->build(isorecotrkB);
         vec_refitTracks.push_back(isotransienttrackA); vec_refitTracks.push_back(isotransienttrackB);

         // Fit tracks:
         AdaptiveVertexFitter  thefitterll(GeometricAnnealing(2.5));
         TransientVertex myVertex = thefitterll.vertex(vec_refitTracks);
         const reco::Vertex secV = myVertex;

         // If the vertex is valid get the geometric information:
         if (secV.isValid()) {

             hasValidVertex = true;

             GlobalVector axis(0,0,0);
             axis = GlobalVector(secV.x(),secV.y(),secV.z());
             Measurement1D vMeas = reco::SecondaryVertex::computeDist2d(pv,secV,axis,true);
             
             // Values:
             vertexLxy = vMeas.value();
             vertexIxy = vMeas.significance();
             normalizedChi2 = myVertex.normalisedChiSquared();         

             leadingPt = (it_A.pt()>it_B.pt())? it_A.pt(): it_B.pt();
             subleadingPt = (it_A.pt()<it_B.pt())? it_A.pt(): it_B.pt();
             trackDxy = (fabs(pckCand_A->dxy(pv.position())/it_A.dxyError()) < fabs(pckCand_B->dxy(pv.position()))/it_B.dxyError())? pckCand_A->dxy(pv.position()): pckCand_B->dxy(pv.position());
             trackIxy = (fabs(pckCand_A->dxy(pv.position())/it_A.dxyError()) < fabs(pckCand_B->dxy(pv.position()))/it_B.dxyError())? pckCand_A->dxy(pv.position())/it_A.dxyError(): pckCand_B->dxy(pv.position())/it_B.dxyError();
             etaA = it_A.eta();
             etaB = it_B.eta();   

 
             // Vector angles:  
             TVector3 vec3A(it_A.px(), it_A.py(), it_A.pz());
             TVector3 vec3B(it_B.px(), it_B.py(), it_B.pz());
             TVector3 divec3 = vec3A + vec3B;
             TVector3 vtxvec3(secV.x() - pv.z(), secV.y() - pv.y(), secV.z() - pv.z());
             cosAlpha = TMath::Cos(vec3A.Angle(vec3B));
             dPhi = divec3.DeltaPhi(vtxvec3);
             dR = vec3A.DeltaR(vec3B);

             // Relative isolation:
             const pat::PFIsolation &pfisoA = it_A.pfIsolationDR03();
             relisoA = (fmax(0.0, pfisoA.photonIso() + pfisoA.neutralHadronIso() - 0.5*pfisoA.puChargedHadronIso()) + pfisoA.chargedHadronIso())/it_A.pt();
             const pat::PFIsolation &pfisoB = it_B.pfIsolationDR03();
             relisoB = (fmax(0.0, pfisoB.photonIso() + pfisoB.neutralHadronIso() - 0.5*pfisoB.puChargedHadronIso()) + pfisoB.chargedHadronIso())/it_B.pt();

             // Physical magnitudes:
             TLorentzVector la;
             TLorentzVector lb;

             if (isEE) {
                la.SetPtEtaPhiM(it_A.pt(), it_A.eta(), it_A.phi(), 0.510/1000.0);
                lb.SetPtEtaPhiM(it_B.pt(), it_B.eta(), it_B.phi(), 0.510/1000.0);
             } else {
                la.SetPtEtaPhiM(it_A.pt(), it_A.eta(), it_A.phi(), 105.658/1000.0);
                lb.SetPtEtaPhiM(it_B.pt(), it_B.eta(), it_B.phi(), 105.658/1000.0);
             }

             mass = (la + lb).M();
             ptll = (la + lb).Pt();



         } else {
            hasValidVertex = false;
         }


      } else {
         canFitVertex = false;
         hasValidVertex= false;
      }

   } // end GetTrackParemeters function 



}; // end struct llCandidate DataFormat definition











