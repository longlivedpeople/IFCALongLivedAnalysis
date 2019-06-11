import ROOT as r
from ROOT import gROOT, TCanvas, TFile, TGraphErrors, TLatex
from plotTools import *
import numpy as np
import os
import math


def getDeltaR(phi1, eta1, phi2, eta2):

    deltaEta = abs(eta1 - eta2)
    deltaPhi = abs(phi1 - phi2)

    if deltaPhi > math.pi: deltaPhi = 2*math.pi - deltaPhi

    return math.sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi)



############################### Open the files #################################

file_name = '/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/IFCALongLivedAnalysis/output.root'

File = TFile(file_name)
t = File.Get("Events")

print('Open the file: ' + file_name)
print("'Events' tree with " + str(t.GetEntries()) + ' entries \n')



############################### Output definition #################################

output_dir = '/eos/user/f/fernance/www/LLP/eff_and_fake_plots/'
if not os.path.exists(output_dir): os.mkdir('./'+output_dir)

############################### Histogram definition ##################################

dxy_bin = np.array([0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20., 30., 40., 50., 60., 70., 90., 110.0])



#### -> Efficiency

eff_dxy_genele = r.TEfficiency("eff_dxy_genele",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)
eff_dxy_genmu = r.TEfficiency("eff_dxy_genmu",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)


eff_dxy_ele = r.TEfficiency("eff_dxy_ele",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)
eff_dxy_mu = r.TEfficiency("eff_dxy_mu",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)


eff_dxy_cmsele = r.TEfficiency("eff_dxy_cmsele",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)
eff_dxy_cmsmu = r.TEfficiency("eff_dxy_cmsmu",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)

#### -> Fake rate
fake_dxy_ele = r.TEfficiency("fake_dxy_ele",";Generated |d_{xy}| (cm);Fake rate", len(dxy_bin) -1 , dxy_bin)
fake_dxy_mu = r.TEfficiency("fake_dxy_mu",";Generated |d_{xy}| (cm);Fake rate", len(dxy_bin) -1 , dxy_bin)


############################### Counters
total_ele = 0
total_muon = 0
reco_ele = 0
reco_muon = 0

############################### Loop over the tree entries ##################################

for n in range(0, t.GetEntries()):

    t.GetEntry(n)


    ##### Efficiency plots

    for g in range(0, t.nGenLepton):

        reconstructed = False # initially not reconstructed
        cms_reconstructed = False # initially not reconstructed by cmssw


        if abs(t.GenLeptonSel_pdgId[g]) == 11: total_ele+=1
        if abs(t.GenLeptonSel_pdgId[g]) == 13: total_muon+=1


        ### cmssw official reconstruction efficiency plots:

        # electron channel:

        if (abs(t.GenLeptonSel_pdgId[g]) == 11):

            for e in range(0, t.nElectron):

                dR = getDeltaR(t.GenLeptonSel_phi[g], t.GenLeptonSel_eta[g], t.ElectronSel_phi[e], t.ElectronSel_eta[e])
                if dR < 0.1:
                    cms_reconstructed = True
                    break

            eff_dxy_cmsele.Fill(cms_reconstructed, abs(t.GenLeptonSel_dxy[g]))

        # muon channel:

        if (abs(t.GenLeptonSel_pdgId[g]) == 13):

            for m in range(0, t.nMuon):

                if not t.MuonSel_isMediumMuon[m]: continue

                dR = getDeltaR(t.GenLeptonSel_phi[g], t.GenLeptonSel_eta[g], t.MuonSel_phi[m], t.MuonSel_eta[m])
                if dR < 0.1:
                    cms_reconstructed = True
                    break

            eff_dxy_cmsmu.Fill(cms_reconstructed, abs(t.GenLeptonSel_dxy[g]))


        # if the genlepton does not have a proper association to track or object the genlepton is not reconstructed directly:

        if (t.GenLeptonSel_trackMatch == 99 or t.GenLeptonSel_objectMatch == 99 or t.GenLeptonSel_trackdR[g] > 0.1): 

            if abs(t.GenLeptonSel_pdgId[g]) == 11: 
                eff_dxy_ele.Fill(reconstructed, abs(t.GenLeptonSel_dxy[g]))
                eff_dxy_genele.Fill(False, abs(t.GenLeptonSel_dxy[g]))
            if abs(t.GenLeptonSel_pdgId[g]) == 13: 
                eff_dxy_mu.Fill(reconstructed, abs(t.GenLeptonSel_dxy[g]))
                eff_dxy_genmu.Fill(False, abs(t.GenLeptonSel_dxy[g]))
            continue


        # if the genlepton has a track and an object the genele and genmu plots are filled:

        if abs(t.GenLeptonSel_pdgId[g]) == 11: 
            eff_dxy_genele.Fill(True, abs(t.GenLeptonSel_dxy[g]))
        if abs(t.GenLeptonSel_pdgId[g]) == 13: 
            eff_dxy_genmu.Fill(True, abs(t.GenLeptonSel_dxy[g]))



        # if it does, we access the track and object:

        gt = t.GenLeptonSel_trackMatch[g]
        go = t.GenLeptonSel_objectMatch[g]


        # electron channel:

        if (abs(t.GenLeptonSel_pdgId[g]) == 11):

            # look for a match in recoElectrons:

            for e in range(0, t.nElectronCandidate):

                rt = t.ElectronCandidate_isotrackIdx[e]
                ro = t.ElectronCandidate_photonIdx[e]

                if gt == rt and go == ro:

                    reconstructed = True
                    reco_ele += 1
                    break

            eff_dxy_ele.Fill(reconstructed, abs(t.GenLeptonSel_dxy[g]))


        # muon channel:

        if (abs(t.GenLeptonSel_pdgId[g]) == 13):

            # look for a match in recoMuons:

            for m in range(0, t.nMuonCandidate):

                rt = t.MuonCandidate_isotrackIdx[m]
                ro = t.MuonCandidate_muonTriggerObjectIdx[m]

                if gt == rt and go == ro:

                    reconstructed = True
                    reco_muon += 1
                    break


            eff_dxy_mu.Fill(reconstructed, abs(t.GenLeptonSel_dxy[g]))


    ##### Fake rate plots

    for e in range(t.nElectronCandidate):

        fake = True

        rt = t.ElectronCandidate_isotrackIdx[e]
        ro = t.ElectronCandidate_photonIdx[e]

        for g in range(0, t.nGenLepton):

            if (abs(t.GenLeptonSel_pdgId[g]) != 11): continue

            gt = t.GenLeptonSel_trackMatch[g]
            go = t.GenLeptonSel_objectMatch[g]

            if rt == gt and ro == go: fake = False 

        fake_dxy_ele.Fill(fake, abs(t.IsoTrackSel_dxy[rt]))



    for m in range(t.nMuonCandidate):

        fake = True

        rt = t.MuonCandidate_isotrackIdx[m]
        ro = t.MuonCandidate_muonTriggerObjectIdx[m]

        for g in range(0, t.nGenLepton):

            if (abs(t.GenLeptonSel_pdgId[g]) != 13): continue

            gt = t.GenLeptonSel_trackMatch[g]
            go = t.GenLeptonSel_objectMatch[g]

            if rt == gt and ro == go: fake = False 

        fake_dxy_mu.Fill(fake, abs(t.IsoTrackSel_dxy[rt]))


print('Electrons:', total_ele, reco_ele)
print('Muons:', total_muon, reco_muon)



############################### Plot in canvas ###############################################

gROOT.ProcessLine('.L include/tdrstyle.C')
gROOT.SetBatch(1)
r.setTDRStyle()
print('Using TDR style')


###### CMS custom:

latexCMS, latexCMSExtra = writeCMS(simulation = True)


###### Create initial canvas


c1 = TCanvas("c1", "", 800, 800)
c1.SetTicks(1,1)
r.gStyle.SetOptStat(0)


###### Plotting efficiencies:

efficiencies = [eff_dxy_ele,
                eff_dxy_mu,
                eff_dxy_cmsele,
                eff_dxy_cmsmu,
                eff_dxy_genele,
                eff_dxy_genmu]



for i,eff in enumerate(efficiencies):

    c1.Clear()
    c1.SetLogy(0)

    xlimit = tuneEfficiency(c1, eff)
    c1.Update()

    line_1 = r.TF1("line_1", "1.0", 0., xlimit)
    line_1.SetLineColor(r.kBlue)
    line_1.SetLineWidth(2)
    line_1.Draw("same")
    c1.Update()
    
    latexCMS.Draw()
    latexCMSExtra.Draw()


    c1.SaveAs(output_dir+eff.GetName()+'.png')



###### Plotting fake rates:

fakes = [fake_dxy_ele,
         fake_dxy_mu]



for i,fake in enumerate(fakes):

    c1.Clear()
    c1.SetLogy(0)

    xlimit = tuneEfficiency(c1, fake)
    c1.Update()

    line_1 = r.TF1("line_1", "1.0", 0., xlimit)
    line_1.SetLineColor(r.kBlue)
    line_1.SetLineWidth(2)
    line_1.Draw("same")
    c1.Update()
    
    latexCMS.Draw()
    latexCMSExtra.Draw()


    c1.SaveAs(output_dir+fake.GetName()+'.png')


