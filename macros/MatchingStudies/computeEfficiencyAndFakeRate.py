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

    dR = math.sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi) 
 
    #print('phi1: ', phi1, 'eta1 :', eta1)
    #print('phi2: ', phi2, 'eta2 :', eta2)
    #print('dR: ', dR)

    return dR



############################### Open the files #################################

#file_name = '/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/1500_494_20/output/merged_1500_494_20_output_N200000.root'
#file_name = '/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/1500_494_20/output/v4/merged.root'
file_name = '/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/IFCALongLivedAnalysis/output.root'


File = TFile(file_name)
t = File.Get("Events")

print('Open the file: ' + file_name)
print("'Events' tree with " + str(t.GetEntries()) + ' entries \n')



############################### Output definition #################################

output_dir = '/eos/user/f/fernance/www/LLP/eff_and_fake_plots/'
if not os.path.exists(output_dir): os.mkdir(output_dir)

############################### Histogram definition ##################################

dxy_bin = np.array([0.0,0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20., 30., 40., 50., 60., 70., 90., 110.0])
vxy_bin = np.array([0.0,0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20., 30., 40., 50., 60., 70., 90., 110.0])
v0_bin = np.array([0.0,0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20., 30., 40., 50., 60., 70., 90., 110.0])
#dxy_bin = np.arange(0.0, 1.0, 0.05)
#vxy_bin = np.arange(0.0, 1.0, 0.05)
#v0_bin = np.arange(0.0, 1.0, 0.05)
pt_bin = np.arange(0., 400., 10.)


#### -> Efficiency

## dxy

eff_dxy_genele = r.TEfficiency("eff_dxy_genele",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)
eff_dxy_genmu = r.TEfficiency("eff_dxy_genmu",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)


eff_v0_genele = r.TEfficiency("eff_v0_genele",";Generated L_{xy} (cm);Efficiency", len(v0_bin) -1 , v0_bin)
eff_v0_genmu = r.TEfficiency("eff_v0_genmu",";Generated L_{xy} (cm);Efficiency", len(v0_bin) -1 , v0_bin)

eff_vxy_genele = r.TEfficiency("eff_vxy_genele",";Generated |v_{xy}| (cm);Efficiency", len(vxy_bin) -1 , vxy_bin)
eff_vxy_genmu = r.TEfficiency("eff_vxy_genmu",";Generated |v_{xy}| (cm);Efficiency", len(vxy_bin) -1 , vxy_bin)

eff_dxy_ele = r.TEfficiency("eff_dxy_ele",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)
eff_dxy_mu = r.TEfficiency("eff_dxy_mu",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)

eff_v0_ele = r.TEfficiency("eff_v0_ele",";Generated L_{xy} (cm);Efficiency", len(v0_bin) -1 , v0_bin)
eff_v0_mu = r.TEfficiency("eff_v0_mu",";Generated L_{xy} (cm);Efficiency", len(v0_bin) -1 , v0_bin)

eff_vxy_ele = r.TEfficiency("eff_vxy_ele",";Generated |v_{xy}| (cm);Efficiency", len(vxy_bin) -1 , vxy_bin)
eff_vxy_mu = r.TEfficiency("eff_vxy_mu",";Generated |v_{xy}| (cm);Efficiency", len(vxy_bin) -1 , vxy_bin)

eff_dxy_cmsele = r.TEfficiency("eff_dxy_cmsele",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)
eff_dxy_cmsmu = r.TEfficiency("eff_dxy_cmsmu",";Generated |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)

eff_v0_cmsele = r.TEfficiency("eff_v0_cmsele",";Generated L_{xy} (cm);Efficiency", len(v0_bin) -1 , v0_bin)
eff_v0_cmsmu = r.TEfficiency("eff_v0_cmsmu",";Generated L_{xy} (cm);Efficiency", len(v0_bin) -1 , v0_bin)

eff_vxy_cmsele = r.TEfficiency("eff_vxy_cmsele",";Generated |v_{xy}| (cm);Efficiency", len(vxy_bin) -1 , vxy_bin)
eff_vxy_cmsmu = r.TEfficiency("eff_vxy_cmsmu",";Generated |v_{xy}| (cm);Efficiency", len(vxy_bin) -1 , vxy_bin)

## pt

eff_pt_genele = r.TEfficiency("eff_pt_genele",";Generated |p_{T}| (GeV);Efficiency", len(pt_bin) -1 , pt_bin)
eff_pt_genmu = r.TEfficiency("eff_pt_genmu",";Generated |p_{T}| (GeV);Efficiency", len(pt_bin) -1 , pt_bin)

eff_pt_ele = r.TEfficiency("eff_pt_ele",";Generated |p_{T}| (GeV);Efficiency", len(pt_bin) -1 , pt_bin)
eff_pt_mu = r.TEfficiency("eff_pt_mu",";Generated |p_{T}| (GeV);Efficiency", len(pt_bin) -1 , pt_bin)


eff_pt_cmsele = r.TEfficiency("eff_pt_cmsele",";Generated |p_{T}| (GeV);Efficiency", len(pt_bin) -1 , pt_bin)
eff_pt_cmsmu = r.TEfficiency("eff_pt_cmsmu",";Generated |p_{T}| (GeV);Efficiency", len(pt_bin) -1 , pt_bin)

#### -> Fake rate
fake_dxy_ele = r.TEfficiency("fake_dxy_ele",";Generated |d_{xy}| (cm);Fake rate", len(dxy_bin) -1 , dxy_bin)
fake_dxy_mu = r.TEfficiency("fake_dxy_mu",";Generated |d_{xy}| (cm);Fake rate", len(dxy_bin) -1 , dxy_bin)

#### -> Histograms

hist_deltaPt_ele = r.TH1F('hist_deltaPt_ele', '', 40, -2, 2)
hist_deltaPt_mu = r.TH1F('hist_deltaPt_mu', '', 40, -2, 2)

hist_deltaEta_ele = r.TH1F('hist_deltaEta_ele', '', 40, -1, 1)
hist_deltaEta_mu = r.TH1F('hist_deltaEta_mu', '', 40, -1, 1)

hist_deltaVx_ele = r.TH1F('hist_deltaVx_ele', '', 40, -5, 5)
hist_deltaVy_ele = r.TH1F('hist_deltaVy_ele', '', 40, -5, 5)
hist_deltaVz_ele = r.TH1F('hist_deltaVz_ele', '', 40, -5, 5)

############################### Counters
total_ele = 0
total_muon = 0
reco_ele = 0
reco_muon = 0

############################### Loop over the tree entries ##################################

for n in range(0, t.GetEntries()):
#for n in range(0, 20000):

    t.GetEntry(n)

    cmsele_matched = []
    cmsmu_matched = []

    #if (t.nGenLepton != 2): continue

    ##### Efficiency plots

    for g in range(0, t.nGenLepton):

        vxy = math.sqrt(t.GenLeptonSel_vx[g]**2 + t.GenLeptonSel_vy[g]**2)
        v0 = math.sqrt((t.GenLeptonSel_vx[g] - t.PV_vx)**2 + (t.GenLeptonSel_vy[g] - t.PV_vy)**2)
        #if(abs(t.GenLeptonSel_dxy[g]) > 1): continue
        #if(vxy > 1): continue

       
        if t.GenLeptonSel_pt[g] < 35: continue
 
        #if (abs(t.GenLeptonSel_pdgId[g]) == 11 and t.GenLeptonSel_pt[g] < 25): continue
        #if (abs(t.GenLeptonSel_pdgId[g]) == 13 and t.GenLeptonSel_pt[g] < 28): continue
        if (abs(t.GenLeptonSel_eta[g]) > 2): continue


        reconstructed = False # initially not reconstructed
        cms_reconstructed = False # initially not reconstructed by cmssw


        if abs(t.GenLeptonSel_pdgId[g]) == 11: total_ele+=1
        if abs(t.GenLeptonSel_pdgId[g]) == 13: total_muon+=1




        ### cmssw official reconstruction efficiency plots:

        # electron channel:

        if (abs(t.GenLeptonSel_pdgId[g]) == 11):
            
            for e in range(0, t.nElectron):

                if e in cmsele_matched: continue
                if not t.ElectronSel_isLoose[e]: continue

                dR = getDeltaR(t.GenLeptonSel_phi[g], t.GenLeptonSel_eta[g], t.ElectronSel_phi[e], t.ElectronSel_eta[e])
                if dR < 0.1:
                    cms_reconstructed = True
                    cmsele_matched.append(e)
                    break

            eff_dxy_cmsele.Fill(cms_reconstructed, abs(t.GenLeptonSel_dxy[g]))
            eff_vxy_cmsele.Fill(cms_reconstructed, vxy)
            eff_v0_cmsele.Fill(cms_reconstructed, v0)
            eff_pt_cmsele.Fill(cms_reconstructed, t.GenLeptonSel_pt[g])

        # muon channel:

        if (abs(t.GenLeptonSel_pdgId[g]) == 13):

            for m in range(0, t.nMuon):

                if m in cmsmu_matched: continue
                if not t.MuonSel_isMediumMuon[m]: continue

                dR = getDeltaR(t.GenLeptonSel_phi[g], t.GenLeptonSel_eta[g], t.MuonSel_phi[m], t.MuonSel_eta[m])
                if dR < 0.1:
                    cms_reconstructed = True
                    cmsmu_matched.append(m)
                    break

            eff_dxy_cmsmu.Fill(cms_reconstructed, abs(t.GenLeptonSel_dxy[g]))
            eff_pt_cmsmu.Fill(cms_reconstructed, t.GenLeptonSel_pt[g])
            eff_vxy_cmsmu.Fill(cms_reconstructed, vxy)
            eff_v0_cmsmu.Fill(cms_reconstructed, v0)

        # if the genlepton does not have a proper association to track or object the genlepton is not reconstructed directly:

        if (t.GenLeptonSel_trackMatch[g] == 99 or t.GenLeptonSel_objectMatch[g] == 99 or t.GenLeptonSel_trackdR[g] > 0.1 or t.GenLeptonSel_pairdR[g] > 0.1): 

            if abs(t.GenLeptonSel_pdgId[g]) == 11: 
                eff_dxy_ele.Fill(reconstructed, abs(t.GenLeptonSel_dxy[g]))
                eff_dxy_genele.Fill(False, abs(t.GenLeptonSel_dxy[g]))
                eff_pt_ele.Fill(reconstructed, t.GenLeptonSel_pt[g])
                eff_pt_genele.Fill(False, t.GenLeptonSel_pt[g])
                eff_vxy_ele.Fill(False, vxy)
                eff_v0_ele.Fill(False, v0)
                eff_vxy_genele.Fill(False, vxy)
                eff_v0_genele.Fill(False, v0)
            if abs(t.GenLeptonSel_pdgId[g]) == 13: 
                eff_dxy_mu.Fill(reconstructed, abs(t.GenLeptonSel_dxy[g]))
                eff_dxy_genmu.Fill(False, abs(t.GenLeptonSel_dxy[g]))
                eff_pt_mu.Fill(reconstructed, t.GenLeptonSel_pt[g])
                eff_pt_genmu.Fill(False, t.GenLeptonSel_pt[g])
                eff_vxy_mu.Fill(False, vxy)
                eff_v0_mu.Fill(False, v0)
                eff_vxy_genmu.Fill(False, vxy)
                eff_v0_genmu.Fill(False, v0)
            continue


        # if the genlepton has a track and an object the genele and genmu plots are filled:

        if abs(t.GenLeptonSel_pdgId[g]) == 11: 
            eff_dxy_genele.Fill(True, abs(t.GenLeptonSel_dxy[g]))
            eff_pt_genele.Fill(True, t.GenLeptonSel_pt[g])
            eff_vxy_genele.Fill(True, vxy)
            eff_v0_genele.Fill(True, v0)
        if abs(t.GenLeptonSel_pdgId[g]) == 13: 
            eff_dxy_genmu.Fill(True, abs(t.GenLeptonSel_dxy[g]))
            eff_pt_genmu.Fill(True, t.GenLeptonSel_pt[g])
            eff_vxy_genmu.Fill(True, vxy)
            eff_v0_genmu.Fill(True, v0)



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

                    hist_deltaPt_ele.Fill((t.GenLeptonSel_pt[g] - t.ElectronCandidate_pt[e])/t.GenLeptonSel_pt[g])
                    hist_deltaEta_ele.Fill((t.GenLeptonSel_eta[g] - t.ElectronCandidate_eta[e])/t.GenLeptonSel_eta[g])
                    hist_deltaVx_ele.Fill(t.GenLeptonSel_vx[g] - t.IsoTrackSel_vx[t.ElectronCandidate_isotrackIdx[e]])
                    hist_deltaVy_ele.Fill(t.GenLeptonSel_vy[g] - t.IsoTrackSel_vy[t.ElectronCandidate_isotrackIdx[e]])
                    hist_deltaVz_ele.Fill(t.GenLeptonSel_vz[g] - t.IsoTrackSel_vz[t.ElectronCandidate_isotrackIdx[e]])
                    reconstructed = True
                    reco_ele += 1
                    break

            eff_dxy_ele.Fill(reconstructed, abs(t.GenLeptonSel_dxy[g]))
            eff_pt_ele.Fill(reconstructed, t.GenLeptonSel_pt[g])
            eff_vxy_ele.Fill(reconstructed, vxy)
            eff_v0_ele.Fill(reconstructed, v0)


        # muon channel:

        if (abs(t.GenLeptonSel_pdgId[g]) == 13):

            # look for a match in recoMuons:

            for m in range(0, t.nMuonCandidate):

                rt = t.MuonCandidate_isotrackIdx[m]
                ro = t.MuonCandidate_muonTriggerObjectIdx[m]

                if gt == rt and go == ro:

                    reconstructed = True
                    hist_deltaPt_mu.Fill((t.GenLeptonSel_pt[g] - t.MuonCandidate_pt[m])/t.GenLeptonSel_pt[g])
                    hist_deltaEta_mu.Fill((t.GenLeptonSel_eta[g] - t.MuonCandidate_eta[m])/t.GenLeptonSel_eta[g])
                    reco_muon += 1
                    break


            eff_dxy_mu.Fill(reconstructed, abs(t.GenLeptonSel_dxy[g]))
            eff_pt_mu.Fill(reconstructed, t.GenLeptonSel_pt[g])
            eff_vxy_mu.Fill(reconstructed, vxy)
            eff_v0_mu.Fill(reconstructed, v0)

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


c1 = TCanvas("c1", "")
#c1.SetTicks(1,1)
#c1.SetGridx()
#c1.SetGridy()
#r.gStyle.SetOptStat(0)


###### Plotting efficiencies:

efficiencies = [eff_dxy_ele, # recoele
                eff_dxy_mu,
                eff_vxy_ele,
                eff_vxy_mu,
                eff_v0_ele,
                eff_v0_mu,
                eff_dxy_cmsele, # cmsele
                eff_dxy_cmsmu,
                eff_vxy_cmsele,
                eff_vxy_cmsmu,
                eff_v0_cmsele,
                eff_v0_cmsmu,
                eff_dxy_genele, # genele
                eff_dxy_genmu,
                eff_vxy_genele,
                eff_vxy_genmu,
                eff_v0_genele,
                eff_v0_genmu,
                eff_pt_ele,
                eff_pt_mu,
                eff_pt_cmsele,
                eff_pt_cmsmu,
                eff_pt_genele,
                eff_pt_genmu]



for i,eff in enumerate(efficiencies):

    c1.Clear()
    c1.SetLogy(0)

    xlimit = tuneEfficiency(c1, eff)
    c1.Update()

    #latexTotal, latexPassed = plotEffStat(eff)

    """
    line_1 = r.TF1("line_1", "1.0", 0., xlimit)
    line_1.SetLineColor(r.kBlue)
    line_1.SetLineWidth(2)
    line_1.Draw("same")
    c1.Update()
    """    

    latexCMS.Draw()
    latexCMSExtra.Draw()
    #latexTotal.Draw()
    #latexPassed.Draw()

    c1.SaveAs(output_dir+eff.GetName()+'.png')



###### Plotting fake rates:

fakes = [fake_dxy_ele,
         fake_dxy_mu]



for i,fake in enumerate(fakes):

    c1.Clear()
    c1.SetLogy(0)

    xlimit = tuneEfficiency(c1, fake)
    c1.Update()

    """
    line_1 = r.TF1("line_1", "1.0", 0., xlimit)
    line_1.SetLineColor(r.kBlue)
    line_1.SetLineWidth(2)
    line_1.Draw("same")
    c1.Update()
    """    

    latexCMS.Draw()
    latexCMSExtra.Draw()


    c1.SaveAs(output_dir+fake.GetName()+'.png')

 ###### Plot joint efficiencies

# -> Electron

ele_pairs = [[eff_dxy_cmsele, eff_dxy_ele],
             [eff_vxy_cmsele, eff_vxy_ele],
             [eff_v0_cmsele, eff_v0_ele]]

ele_tag_pairs = ['eff_dxy_ele_comparison',
                 'eff_vxy_ele_comparison',
                 'eff_v0_ele_comparison']

for p in range(0,len(ele_pairs)):

    c1.Clear()
    c1.SetLogy(0)
    xlimit = plotJointEfficiency(c1, ele_pairs[p])
    c1.Update()
    latexCMS.Draw()
    latexCMSExtra.Draw()
    leg = drawEffLegend(ele_pairs[p], ['Loose electrons', 'Electron matching'])
    leg.Draw()
    c1.SaveAs(output_dir+ele_tag_pairs[p]+'.png')


# -> Muon

mu_pairs = [[eff_dxy_cmsmu, eff_dxy_mu],
             [eff_vxy_cmsmu, eff_vxy_mu],
             [eff_v0_cmsmu, eff_v0_mu]]

mu_tag_pairs = ['eff_dxy_mu_comparison',
                 'eff_vxy_mu_comparison',
                 'eff_v0_mu_comparison']

for p in range(0,len(mu_pairs)):

    c1.Clear()
    c1.SetLogy(0)
    xlimit = plotJointEfficiency(c1, mu_pairs[p])
    c1.Update()
    latexCMS.Draw()
    latexCMSExtra.Draw()
    leg = drawEffLegend(mu_pairs[p], ['Medium muons', 'Muon matching'])
    leg.Draw()
    c1.SaveAs(output_dir+mu_tag_pairs[p]+'.png')


###### Plot pt difference histograms

tuneEmptyHisto(hist_deltaPt_ele, '(p_{T}^{gen} - p_{T}^{e})/p_{T}^{gen}', r.kRed, False)
tuneEmptyHisto(hist_deltaPt_mu, '(p_{T}^{gen} - p_{T}^{#mu})/p_{T}^{gen}', r.kRed, False)
tuneEmptyHisto(hist_deltaEta_ele, '(#eta^{gen} - #eta^{e})/#eta^{gen}', r.kRed, False)
tuneEmptyHisto(hist_deltaEta_mu, '(#eta^{gen} - #eta^{#mu})/#eta^{gen}', r.kRed, False)
tuneEmptyHisto(hist_deltaVz_ele, '(v_{z}^{gen} - v_{z}^{e})', r.kRed, False)
tuneEmptyHisto(hist_deltaVy_ele, '(v_{y}^{gen} - v_{y}^{e})', r.kRed, False)
tuneEmptyHisto(hist_deltaVx_ele, '(v_{x}^{gen} - v_{x}^{e})', r.kRed, False)

hist_deltaPt_list = [hist_deltaPt_ele,
                     hist_deltaPt_mu,
                     hist_deltaEta_ele,
                     hist_deltaEta_mu,
                     hist_deltaVz_ele,
                     hist_deltaVy_ele,
                     hist_deltaVx_ele]

c1.SetGridx(0)
c1.SetGridy(0)

for histo in hist_deltaPt_list:

    c1.Clear()
    c1.SetLogy(0)
    printHisto(histo)
    latexCMS.Draw()
    latexCMSExtra.Draw()

    c1.SaveAs(output_dir+histo.GetName()+'.png')

