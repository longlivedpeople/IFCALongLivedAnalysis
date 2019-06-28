""" Code to check the matching between the isoTrack collection and the genParticle collection"""

import ROOT
import math
import numpy as np



file_ = ROOT.TFile("/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/IFCALongLivedAnalysis/output.root")

tree = file_.Get("Events")


# Efficiency histograms
#dxy_bin = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20., 25.0, 30., 40., 50., 60., 70., 90., 110.0]
dxy_bin = [0.0, 10.0, 20., 30., 40., 50., 60., 70., 90., 110.0]
#dxy_bin = [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 20., 30., 40., 50., 60., 70., 90., 110.0]
eff_dxy_ele = ROOT.TEfficiency("eff",";Generated electron transverse impact parameter |d_{xy}| (cm);N_{matched} / N_{generated}", len(dxy_bin) -1 , np.array(dxy_bin))
eff_dxy_mu = ROOT.TEfficiency("eff",";Generated muon transverse impact parameter |d_{xy}| (cm);N_{matched} / N_{generated}", len(dxy_bin) -1 , np.array(dxy_bin))


total_ele = 0
total_mu = 0
matched_ele = 0
matched_mu = 0


for n,event in enumerate(tree):

    print(">>> Event number: " + str(n))


    ## Selection of the genleptons
    for i in range(0, event.nGenLepton):

        if event.GenLeptonSel_pt[i] < 20:
            continue
        if abs(event.GenLeptonSel_eta[i]) > 2.4:
            continue

        if abs(event.GenLeptonSel_pdgId[i]) == 11: # Electron channel
            
            total_ele += 1

            if event.nElectronCandidate == 0: matched == False

            for e in range(event.nElectronCandidate):

                matched = True    

                if(abs(event.GenLeptonSel_pt[i] - event.ElectronCandidate_pt[e])/event.GenLeptonSel_pt[i] > 0.4): matched = False
                if(abs(event.GenLeptonSel_pt[i] - event.ElectronCandidate_et[e])/event.GenLeptonSel_pt[i] > 0.4): matched = False

                deltaPhi = event.GenLeptonSel_phi[i] - event.ElectronCandidate_phi[e]
                if deltaPhi > 3.14: deltaPhi = deltaPhi - 2.0 * 3.14
                if deltaPhi < -3.14: deltaPhi = deltaPhi + 2.0 * 3.14

                deltaEta = abs(event.GenLeptonSel_eta[i] - event.ElectronCandidate_eta[e])
                deltaR = math.sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta)

                if deltaR > 0.4: matched = False

                if matched: 
                    matched_ele +=1
                    break


            eff_dxy_ele.Fill(matched, abs(event.GenLeptonSel_dxy[i]))    

        if abs(event.GenLeptonSel_pdgId[i]) == 13: #Muon channel

            total_mu += 1
            if event.nMuonCandidate == 0: matched == False

            for e in range(event.nMuonCandidate):

                matched = True

                if(abs(event.GenLeptonSel_pt[i] - event.MuonCandidate_pt[e])/event.GenLeptonSel_pt[i] > 0.4): matched = False
                if(abs(event.GenLeptonSel_pt[i] - event.MuonCandidate_triggerPt[e])/event.GenLeptonSel_pt[i] > 0.4): matched = False

                deltaPhi = event.GenLeptonSel_phi[i] - event.MuonCandidate_phi[e]
                if deltaPhi > 3.14: deltaPhi = deltaPhi - 2.0 * 3.14
                if deltaPhi < -3.14: deltaPhi = deltaPhi + 2.0 * 3.14

                deltaEta = abs(event.GenLeptonSel_eta[i] - event.MuonCandidate_eta[e])
                deltaR = math.sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta)

                if deltaR > 0.4: matched = False

                if matched: 
                    matched_mu +=1
                    break


            eff_dxy_mu.Fill(matched, abs(event.GenLeptonSel_dxy[i]))



print(matched_ele)
print(matched_mu)


c1 = ROOT.TCanvas("c1", "", 500, 500)
c1.SetTicks(1,1)
ROOT.gStyle.SetOptStat(0)


#######3# Efficiency histogram
c1.Clear()
c1.SetLogy(0)
# Histogram tunning
eff_dxy_ele.SetMarkerStyle(21)
eff_dxy_ele.SetMarkerColor(ROOT.kBlack)
eff_dxy_ele.SetLineWidth(2)
eff_dxy_ele.SetLineColor(ROOT.kBlack)
eff_dxy_ele.SetMarkerSize(1.)


eff_dxy_mu.SetMarkerStyle(21)
eff_dxy_mu.SetMarkerColor(ROOT.kBlack)
eff_dxy_mu.SetLineWidth(2)
eff_dxy_mu.SetLineColor(ROOT.kBlack)
eff_dxy_mu.SetMarkerSize(1.)

line_1 = ROOT.TF1("line_1", "1.0", 0., 110.)
line_1.SetLineColor(ROOT.kRed)
line_1.SetLineWidth(2)
#line_1.Draw()
#eff_dxy.GetListOfFunctions().Add(line_1)

eff_dxy_ele.Draw("AP")

c1.Update()
graph_ele = eff_dxy_ele.GetPaintedGraph(); 
graph_ele.SetMinimum(0.);
graph_ele.SetMaximum(1.1);
graph_ele.GetXaxis().SetLimits(0., 110.)
graph_ele.GetXaxis().SetTitleOffset(1.2)

c1.Update()
line_1.Draw("same")
c1.SaveAs("Ele_efficiency.png")

c1.Clear()
eff_dxy_mu.Draw("AP")
c1.Update()
graph_mu = eff_dxy_mu.GetPaintedGraph(); 
graph_mu.SetMinimum(0.);
graph_mu.SetMaximum(1.1);
graph_mu.GetXaxis().SetLimits(0., 110.)
graph_mu.GetXaxis().SetTitleOffset(1.1)

c1.Update()
line_1.Draw("same")
c1.SaveAs("Muon_efficiency.png")

# Text boxes
"""
t1 = ROOT.TText(50., 0.9, "Total generated leptons: " + str(particle))
t2 = ROOT.TText(50., 0.85, "Total matched leptons: "+ str(track))
t1.SetTextSize(0.03)
t2.SetTextSize(0.03)
t1.Draw()
t2.Draw()
"""






