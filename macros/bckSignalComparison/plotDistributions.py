import ROOT as r
from ROOT import gROOT, TCanvas, TFile, TGraphErrors, TLatex, TH1F, TLorentzVector, TVector3
#from plotTools import *
import os
import math


############################## Functions and constants definition #######################################

ELECTRON_MASS = 0.511*1E-3 # electron mass in GeV
MUON_MASS = 105.65*1E-3 # muon mass in GeV


def deltaPhi(phi1, phi2):

    deltaPhi = phi1 - phi2
    abs_deltaPhi = 2*math.pi - abs(deltaPhi)

    return abs_deltaPhi



############################## Open the files and create data structure #################################

## File structure:

files = {}
files['Signal'] = ['/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/IFCALongLivedAnalysis/output.root']
files['DY'] = ['/eos/cms/store/user/pablom/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v1.root']


## Tree structure:

trees = {}
for process in files.keys():
    
    trees[process] = r.TChain("Events", "")
    
    for name in files[process]:

        trees[process].Add(name)


############################### Define the LL-candidate variables #############################

nLL = 0
LL_leptons = []
LL_type = [] # 0 if EE or 1 if MM
LL_Ixy = []
LL_invariantMass = []
LL_leadingPt = []
LL_subleadingPt = []
LL_aperturePhi = []
LL_vertexPosition = []
LL_vertexChi2ndf = []

LLSel_leptons = 99
LLSel_type = 99 # 0 if EE or 1 if MM
LLSel_Ixy = 99
LLSel_invariantMass = 99
LLSel_leadingPt = 99
LLSel_subleadingPt = 99
LLSel_aperturePhi = 99
LLSel_vertexPosition = 99
LLSel_vertexChi2ndf = 99


############################### Histogram definition ################################

## histagram structure definition

histos = {}
for process in files.keys():

    histos[process] = {}

    histos[process]['nElectronCandidate'] = TH1F(process + '_' +'nElectronCandidate', '', 7, 0, 7)
    histos[process]['nMuonCandidate'] = TH1F(process + '_' +'nMuonCandidate', '', 7, 0, 7)
    histos[process]['ElectronCandidate_pt'] = TH1F(process + '_' +'ElectronCandidate_pt', '', 50, 0, 1200)
    histos[process]['MuonCandidate_pt'] = TH1F(process + '_' +'MuonCandidate_pt', '', 50, 0, 1200)
    histos[process]['nLL'] = TH1F(process + '_' +'nLL', '', 30, 0, 30)
    histos[process]['nEE'] = TH1F(process + '_' +'nEE', '', 30, 0, 30)
    histos[process]['nMM'] = TH1F(process + '_' +'nMM', '', 30, 0, 30)
    histos[process]['EE_dist'] = TH1F(process + '_' +'EE_dist', '', 20, 0, 10)
    histos[process]['EE_invMass'] = TH1F(process + '_' +'EE_invMass', '', 50, 0, 500)
    histos[process]['MM_invMass'] = TH1F(process + '_' +'MM_invMass', '', 50, 0, 500)
    histos[process]['EE_colinearity'] = TH1F(process + '_' +'EE_colinearity', '', 20, 0, math.pi)
    histos[process]['MM_colinearity'] = TH1F(process + '_' +'MM_colinearity', '', 20, 0, math.pi)
    histos[process]['EESel_colinearity'] = TH1F(process + '_' +'EESel_colinearity', '', 20, 0, math.pi)
    histos[process]['MMSel_colinearity'] = TH1F(process + '_' +'MMSel_colinearity', '', 20, 0, math.pi)


maxEntries = 18000

expected_nEE = 0
expected_nMM = 0

############################### Loop over the events ##############################


for process in files.keys():

    t = trees[process] # Select the tree process

    print('accesed: ' + process)

    for n in range(0, t.GetEntries()):

        print(n)

        if n > maxEntries: break

        t.GetEntry(n) # Access each event

        expected_nEE = expected_nEE + t.nElectronCandidate*(t.nElectronCandidate -1)/2.
        expected_nMM = expected_nMM + t.nMuonCandidate*(t.nMuonCandidate -1)/2.


        nLL = 0
        LL_leptons = []
        LL_type = [] # 0 if EE or 1 if MM
        LL_Ixy = []
        LL_invariantMass = []
        LL_pt = []
        LL_leadingPt = []
        LL_subleadingPt = []
        LL_aperturePhi = []
        LL_vertexPosition = []
        LL_vertexChi2ndf = []
        LL_vx = []
        LL_vy = []
        LL_vz = []
        LL_dist = [] # Distance of one track from the other
        LL_colinearity = []

        LLSel_leptons = 99
        LLSel_type = 99 # 0 if EE or 1 if MM
        LLSel_Ixy = 99
        LLSel_invariantMass = 99
        LLSel_leadingPt = 99
        LLSel_subleadingPt = 99
        LLSel_aperturePhi = 99
        LLSel_vertexPosition = 99
        LLSel_vertexChi2ndf = 99
        LLSel_colinearity = 99

        ## LL-candidates reconstruction

        # Electron channel
        for i in range(1, t.nElectronCandidate):
            for j in range(0, i):

                if i == j : continue

                if (t.IsoTrackSel_fromPV[t.ElectronCandidate_isotrackIdx[i]] < 0 or t.IsoTrackSel_fromPV[t.ElectronCandidate_isotrackIdx[j]] < 0): continue
                #if (t.IsoTrackSel_goodPackedCandidate[t.ElectronCandidate_isotrackIdx[i]] == 0 or t.IsoTrackSel_goodPackedCandidate[t.ElectronCandidate_isotrackIdx[j]] == 0 ): continue


                LL_leptons.append([i, j])
                LL_type.append(0)

                Ixy = t.ElectronCandidate_dxySignificance[i] if t.ElectronCandidate_dxySignificance[i] < t.ElectronCandidate_dxySignificance[j] else t.ElectronCandidate_dxySignificance[j]
                leadingPt = t.ElectronCandidate_pt[i] if t.ElectronCandidate_pt[i] > t.ElectronCandidate_pt[j] else t.ElectronCandidate_pt[j]
                subleadingPt = t.ElectronCandidate_pt[i] if t.ElectronCandidate_pt[i] < t.ElectronCandidate_pt[j] else t.ElectronCandidate_pt[j]
            
                LL_Ixy.append(Ixy)
                LL_leadingPt.append(leadingPt)
                LL_subleadingPt.append(subleadingPt)
                LL_aperturePhi.append(deltaPhi(t.ElectronCandidate_phi[i], t.ElectronCandidate_phi[j])) 

                # Invariant mass
                vi = TLorentzVector()
                vj = TLorentzVector()
                vi.SetPtEtaPhiM(t.ElectronCandidate_pt[i], t.ElectronCandidate_eta[i], t.ElectronCandidate_phi[i], ELECTRON_MASS)
                vj.SetPtEtaPhiM(t.ElectronCandidate_pt[j], t.ElectronCandidate_eta[j], t.ElectronCandidate_phi[j], ELECTRON_MASS)
                LL_invariantMass.append((vi+vj).M())
                LL_pt.append((vi+vj).Pt())

                # Fitted vertex
                dv_x = (t.IsoTrackSel_vx[t.ElectronCandidate_isotrackIdx[i]] + t.IsoTrackSel_vx[t.ElectronCandidate_isotrackIdx[j]])/2
                dv_y = (t.IsoTrackSel_vy[t.ElectronCandidate_isotrackIdx[i]] + t.IsoTrackSel_vy[t.ElectronCandidate_isotrackIdx[j]])/2
                dv_z = (t.IsoTrackSel_vz[t.ElectronCandidate_isotrackIdx[i]] + t.IsoTrackSel_vz[t.ElectronCandidate_isotrackIdx[j]])/2

                LL_vx.append(dv_x)
                LL_vy.append(dv_y)
                LL_vz.append(dv_z)

                LL_dist.append( math.sqrt( (t.IsoTrackSel_vx[t.ElectronCandidate_isotrackIdx[i]]- t.IsoTrackSel_vx[t.ElectronCandidate_isotrackIdx[j]])**2 + (t.IsoTrackSel_vy[t.ElectronCandidate_isotrackIdx[i]]- t.IsoTrackSel_vy[t.ElectronCandidate_isotrackIdx[j]])**2)  )

                # Colinearity
                vertexVector = TVector3()
                vertexVector.SetXYZ(dv_x - t.PV_vx, dv_y - t.PV_vy, dv_z - t.PV_vz)
                
                LLVector = TVector3()
                LLVector.SetPtEtaPhi((vi+vj).Pt(), (vi+vj).Eta(), (vi+vj).Phi())

                print('colinearity', abs(vertexVector.DeltaPhi(LLVector)))
                print(dv_x - t.PV_vx, dv_y - t.PV_vy)
                print((vi+vj).Phi())

                LL_colinearity.append(abs(vertexVector.DeltaPhi(LLVector)))


        # Muon channel
        for i in range(1, t.nMuonCandidate):
            for j in range(0, i):

                if i == j: continue

                if (t.IsoTrackSel_fromPV[t.MuonCandidate_isotrackIdx[i]] < 0 or t.IsoTrackSel_fromPV[t.MuonCandidate_isotrackIdx[j]] < 0): continue
                #if (t.IsoTrackSel_goodPackedCandidate[t.MuonCandidate_isotrackIdx[i]] == 0 or t.IsoTrackSel_goodPackedCandidate[t.MuonCandidate_isotrackIdx[j]] == 0 ): continue

                LL_leptons.append([i, j])
                LL_type.append(1)

                Ixy = t.MuonCandidate_dxySignificance[i] if t.MuonCandidate_dxySignificance[i] < t.MuonCandidate_dxySignificance[j] else t.MuonCandidate_dxySignificance[j]
                leadingPt = t.MuonCandidate_pt[i] if t.MuonCandidate_pt[i] > t.MuonCandidate_pt[j] else t.MuonCandidate_pt[j]
                subleadingPt = t.MuonCandidate_pt[i] if t.MuonCandidate_pt[i] < t.MuonCandidate_pt[j] else t.MuonCandidate_pt[j]
            
                LL_Ixy.append(Ixy)
                LL_leadingPt.append(leadingPt)
                LL_subleadingPt.append(subleadingPt)
                LL_aperturePhi.append(deltaPhi(t.MuonCandidate_phi[i], t.MuonCandidate_phi[j]))

                # Invariant mass
                vi = TLorentzVector()
                vj = TLorentzVector()
                vi.SetPtEtaPhiM(t.MuonCandidate_pt[i], t.MuonCandidate_eta[i], t.MuonCandidate_phi[i], MUON_MASS)
                vj.SetPtEtaPhiM(t.MuonCandidate_pt[j], t.MuonCandidate_eta[j], t.MuonCandidate_phi[j], MUON_MASS)
                LL_invariantMass.append((vi+vj).M())
                LL_pt.append((vi+vj).Pt())

                # Fitted vertex
                dv_x = (t.IsoTrackSel_vx[t.MuonCandidate_isotrackIdx[i]] + t.IsoTrackSel_vx[t.MuonCandidate_isotrackIdx[j]])/2
                dv_y = (t.IsoTrackSel_vy[t.MuonCandidate_isotrackIdx[i]] + t.IsoTrackSel_vy[t.MuonCandidate_isotrackIdx[j]])/2
                dv_z = (t.IsoTrackSel_vz[t.MuonCandidate_isotrackIdx[i]] + t.IsoTrackSel_vz[t.MuonCandidate_isotrackIdx[j]])/2

                LL_vx.append(dv_x)
                LL_vy.append(dv_y)
                LL_vz.append(dv_z)

                LL_dist.append( math.sqrt( (t.IsoTrackSel_vx[t.MuonCandidate_isotrackIdx[i]]- t.IsoTrackSel_vx[t.MuonCandidate_isotrackIdx[j]])**2 + (t.IsoTrackSel_vy[t.MuonCandidate_isotrackIdx[i]]- t.IsoTrackSel_vy[t.MuonCandidate_isotrackIdx[j]])**2)  )

                # Colinearity
                vertexVector = TVector3()
                vertexVector.SetXYZ(dv_x - t.PV_vx, dv_y - t.PV_vy, dv_z - t.PV_vz)
                
                LLVector = TVector3()
                LLVector.SetPtEtaPhi((vi+vj).Pt(), (vi+vj).Eta(), (vi+vj).Phi())

                LL_colinearity.append(abs(vertexVector.DeltaPhi(LLVector)))


        ## Counters:
        nLL = len(LL_type)
        nEE = 0
        nMM = 0
        for i in range(0,nLL): 
            if LL_type[i] == 0: nEE += 1
            if LL_type[i] == 1: nMM += 1


        ## Get the final LL-candidate (the one with the largest Ixy)
        maxIxy = 0
        for i in range(0, len(LL_type)):

            if LL_dist[i] > 1: continue

            ## Apply LL candidate selection

            if LL_Ixy[i] > maxIxy:

                maxIxy = LL_Ixy[i]

                LLSel_leptons = LL_leptons[i]
                LLSel_type = LL_type[i]
                LLSel_Ixy = LL_Ixy[i]
                LLSel_leadingPt = LL_leadingPt[i]
                LLSel_subleadingPt = LL_subleadingPt[i] 
                LLSel_aperturePhi = LL_aperturePhi[i]
                LLSel_invariantMass = LL_invariantMass[i]
                LLSel_colinearity = LL_colinearity[i]


        ################################ Fill histograms ###########################

        # Electron Candidates
        if t.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8: histos[process]['nElectronCandidate'].Fill(t.nElectronCandidate)
        for i in range(0, t.nElectronCandidate): histos[process]['ElectronCandidate_pt'].Fill(t.ElectronCandidate_pt[i])

        # Muon Candidates
        if t.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6:  histos[process]['nMuonCandidate'].Fill(t.nMuonCandidate)
        for i in range(0, t.nMuonCandidate): histos[process]['MuonCandidate_pt'].Fill(t.MuonCandidate_pt[i])

        histos[process]['nLL'].Fill(nLL)
        histos[process]['nEE'].Fill(nEE)
        histos[process]['nMM'].Fill(nMM)


        # Global LL-candidates

        for i in range(0, len(LL_type)):

            if LL_dist[i] > 1: continue

            if LL_type[i] == 0: # dielectron channel

                histos[process]['EE_dist'].Fill(LL_dist[i])
                histos[process]['EE_invMass'].Fill(LL_invariantMass[i])
                histos[process]['EE_colinearity'].Fill(LL_colinearity[i])
            
            if LL_type[i] == 1: # dimuon channel

                histos[process]['MM_invMass'].Fill(LL_invariantMass[i])
                histos[process]['MM_colinearity'].Fill(LL_colinearity[i])


        # Final LL-candidate

        if LLSel_type == 0: histos[process]['EESel_colinearity'].Fill(LLSel_colinearity)
        if LLSel_type == 1: histos[process]['MMSel_colinearity'].Fill(LLSel_colinearity)


#################################### Output definition ######################################

output_tree = TFile('output_histos.root', 'RECREATE')

for process in files.keys():

    for key in  histos[process].keys():

        histo = histos[process][key]
        #histo.SetName(process + '_' + histo.GetName())
        histo.Write()

output_tree.Close()


"""
### Compute total number of nLL
total_nLL = 0

for b in range(1, histos['Signal']['nLL'].GetNbinsX() +1):
    total_nLL = total_nLL + (b-1)*histos['Signal']['nLL'].GetBinContent(b)

print(total_nLL)

### Compute expected number of nLL
print(expected_nEE + expected_nMM)
"""
