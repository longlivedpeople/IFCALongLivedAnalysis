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
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


#include "TLorentzVector.h"
#include "TVector3.h"

#include <TROOT.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>


////////////////////////////////////////////////////
// ---------------------------------------------- //
// -------- trackPair DataFormat definition -------- //
// ---------------------------------------------- //
////////////////////////////////////////////////////

struct trackPair
{

   // ---- trackPair DataFormat information

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
   double vx = 0;                  // x coordinate of dilepton vertex
   double vy = 0;                  // y coordinate of dilepton vertex


   // PV association quality flags: 
   // 	{0: None of them has participated in the PV fit, 
   // 	 1: One of them has participated in the PV fit,
   // 	 2: Both of them have participated in the PV fit}
   int fromPVA = 0;
   int fromPVB = 0;
   int PVAssociation = 0;


   // ---- trackPair constructors: 
   trackPair(const reco::Vertex &pv, edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder, const reco::Track &tr_A, const reco::Track &tr_B, bool isEE)
   { 
      Init(theTransientTrackBuilder, pv, tr_A, tr_B, isEE); 
   };
 
   ~trackPair(){};


   // ---- trackPair DataFormat functions

   // -- Init trackPair:
   void Init(edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder, const reco::Vertex &pv, const reco::Track &tr_A, const reco::Track &tr_B, bool isEE)
   {

      if (isEE) { type = 0; }
      else { type = 1; }


      // Get tracks:
      std::vector<reco::TransientTrack> vec_refitTracks;
      reco::TransientTrack isotransienttrackA = theTransientTrackBuilder->build(tr_A);
      reco::TransientTrack isotransienttrackB = theTransientTrackBuilder->build(tr_B);
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
          vx = secV.x();
          vy = secV.y();

          leadingPt = (tr_A.pt()>tr_B.pt())? tr_A.pt(): tr_B.pt();
          subleadingPt = (tr_A.pt()<tr_B.pt())? tr_A.pt(): tr_B.pt();
          float trA_dxy = computeDxy(tr_A, pv);
          float trA_dxyError = computeDxyError(theTransientTrackBuilder, tr_A, pv); 
          float trB_dxy = computeDxy(tr_B, pv);
          float trB_dxyError = computeDxyError(theTransientTrackBuilder, tr_B, pv); 

          trackDxy = (fabs(trA_dxy/trA_dxyError) < fabs(trB_dxy/trB_dxyError))? trA_dxy: trB_dxy;
          trackIxy = (fabs(trA_dxy/trA_dxyError) < fabs(trB_dxy/trB_dxyError))? fabs(trA_dxy/trA_dxyError): fabs(trB_dxy/trB_dxyError);
          etaA = tr_A.eta();
          etaB = tr_B.eta();   

 
          // Vector angles:  
          TVector3 vec3A(tr_A.px(), tr_A.py(), tr_A.pz());
          TVector3 vec3B(tr_B.px(), tr_B.py(), tr_B.pz());
          TVector3 divec3 = vec3A + vec3B;
          TVector3 vtxvec3(secV.x() - pv.z(), secV.y() - pv.y(), secV.z() - pv.z());
          cosAlpha = TMath::Cos(vec3A.Angle(vec3B));
          dPhi = divec3.DeltaPhi(vtxvec3);
          dR = vec3A.DeltaR(vec3B);


          // Physical magnitudes:
          TLorentzVector la;
          TLorentzVector lb;

          if (isEE) {
             la.SetPtEtaPhiM(tr_A.pt(), tr_A.eta(), tr_A.phi(), 0.510/1000.0);
             lb.SetPtEtaPhiM(tr_B.pt(), tr_B.eta(), tr_B.phi(), 0.510/1000.0);
          } else {
             la.SetPtEtaPhiM(tr_A.pt(), tr_A.eta(), tr_A.phi(), 105.658/1000.0);
             lb.SetPtEtaPhiM(tr_B.pt(), tr_B.eta(), tr_B.phi(), 105.658/1000.0);
          }

          mass = (la + lb).M();
          ptll = (la + lb).Pt();

          /*
          // PV association quality flags:
          fromPVA = (*pckCand_A).fromPV(0);
          fromPVB = (*pckCand_B).fromPV(0);
          PVAssociation = 0;             

             if (fromPVA == 3 && fromPVB == 3) { PVAssociation = 2; }
             if ((fromPVA == 3 && fromPVB != 3) || (fromPVB == 3 && fromPVA != 3) ) { PVAssociation = 1; }

          */

      } else {
         hasValidVertex = false;
      }


   } // end Init function 



   float computeDxy(const reco::Track & track, const reco::Vertex pv) {

      double vx = track.vx();
      double vy = track.vy();
      double phi = track.phi();
      double PVx = pv.x();
      double PVy = pv.y();

      double dxy = -(vx - PVx)*sin(phi) + (vy - PVy)*cos(phi);
      return dxy;

   } // end computeDxy function


   float computeDxyError(edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder, const reco::Track & track, const reco::Vertex pv) {

      reco::TransientTrack isotk = theTransientTrackBuilder->build(track);
      GlobalPoint vert(pv.x(), pv.y(), pv.z());
      TrajectoryStateClosestToPoint  traj = isotk.trajectoryStateClosestToPoint(vert);

       float sigmaXY = traj.perigeeError().transverseImpactParameterError();

       return sigmaXY;

   } // end computeDxyError function



}; // end struct trackPair DataFormat definition











