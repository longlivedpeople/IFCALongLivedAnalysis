import ROOT as r
from ROOT import gROOT, TCanvas, TFile, TGraphErrors, TLatex, TVector3
from plotTools import *
import os
import numpy as np
import math

############################### Open the files #################################

#file_name = '/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/IFCALongLivedAnalysis/output.root'
file_name = '/afs/cern.ch/work/f/fernance/public/Samples/DisplacedSUSY/2016/1500-494-160/DisplacedSUSY_1500-494-160.root'

File = TFile(file_name)
t = File.Get("Events")

print('Open the file: ' + file_name)
print("'Events' tree with " + str(t.GetEntries()) + ' entries \n')



############################### Output definition #################################

output_dir = '/eos/user/f/fernance/www/LLP/DielectronEfficiency/'
if not os.path.exists(output_dir): os.mkdir(output_dir)


############################### TEfficiency declaration #################################

Lxy_bin = np.array([0.0,0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20., 30., 40., 50., 60., 70., 90., 110.0])

eff_electron_Lxy = r.TEfficiency("eff_electron_Lxy",";Generated L_{xy} (cm);Efficiency", len(Lxy_bin) -1 , Lxy_bin)

eff_dielectron_Lxy = r.TEfficiency("eff_dielectron_Lxy",";Generated L_{xy} (cm);Efficiency", 18, 0, 60)

fake_electron_Lxy = r.TEfficiency("fake_electron_Lxy",";Reconstructed L_{xy} (cm);Fake Rate", len(Lxy_bin) -1 , Lxy_bin)


############################### TProfile declaration #################################

prof_electron_deltaPt_vs_Lxy = r.TProfile("prof_electron_deltaPt_vxLxy", "", 20, 0, 60)
prof_electron_deltaPt_vs_Lxy.GetXaxis().SetTitle('Generated L_{xy} (cm)')
prof_electron_deltaPt_vs_Lxy.GetYaxis().SetTitle('(p_{T}^{gen} - p_{T}^{e})/p_{T}^{gen}')


################################# TH1F declaration ###################################

histo_deltaPt_Lxy_0_1 = r.TH1F("histo_deltaPt_Lxy_0_1", "", 40, -2.0, 2.0)
histo_deltaPt_Lxy_1_2 = r.TH1F("histo_deltaPt_Lxy_1_2", "", 40, -2.0, 2.0)
histo_deltaPt_Lxy_2_4 = r.TH1F("histo_deltaPt_Lxy_2_4", "", 40, -2.0, 2.0)
histo_deltaPt_Lxy_4_8 = r.TH1F("histo_deltaPt_Lxy_4_8", "", 40, -2.0, 2.0)
histo_deltaPt_Lxy_8_inf = r.TH1F("histo_deltaPt_Lxy_8_inf", "", 40, -2.0, 2.0)


############################### Loop over the tree entries ##################################

for n in range(0, t.GetEntries()):

    t.GetEntry(n)

    if not (t.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8): continue

    g_pair = []
    r_pair = []

    gen_e = []
    gen_ei = []

    candidate_status = [False for i in range(0, t.nElectronCandidate)]
    GenLepton_status = [False for i in range(0, t.nGenLepton)]    


    ### GenLepton Identification
    for g in range(0, t.nGenLepton):

        if t.GenLeptonSel_pt[g] < 35: continue
        if (abs(t.GenLeptonSel_eta[g]) > 2): continue

        if abs(t.GenLeptonSel_pdgId[g]) == 11:
            electron = TVector3()
            electron.SetPtEtaPhi(t.GenLeptonSel_pt[g], t.GenLeptonSel_eta[g], t.GenLeptonSel_phi[g])
            gen_e.append(electron)
            gen_ei.append(g)


    ### RecoElectron identification
    reco_e = []
    for i in range(0, t.nElectronCandidate):
        electron = TVector3()
        electron.SetPtEtaPhi(t.ElectronCandidate_pt[i], t.ElectronCandidate_eta[i], t.ElectronCandidate_phi[i])
        reco_e.append(electron)

    ### Single Electron matching efficiency plot
    for eg in range(0,len(gen_e)):

        genElectron = gen_e[eg]
        Lxy = math.sqrt((t.GenLeptonSel_vx[gen_ei[eg]] - t.PV_vx)**2 + (t.GenLeptonSel_vy[gen_ei[eg]] - t.PV_vy)**2)
        reconstructed = False
        reco_id = -99
        mindR = 99

        for er in range(0, len(reco_e)):

            recoElectron = reco_e[er]

            if recoElectron in r_pair: continue
 
            if recoElectron.DeltaR(genElectron) < 0.1 and recoElectron.DeltaR(genElectron) < mindR:

                # Update variables:
                mindR = recoElectron.DeltaR(genElectron)
                reco_id = er
                reconstructed = True

        if reconstructed:

            g_pair.append(genElectron)
            r_pair.append(reco_e[reco_id])
            candidate_status[reco_id] = True
            GenLepton_status[gen_ei[eg]] = True
            recoElectron = reco_e[reco_id]

            ##### Fill plots:
            deltaPt = (genElectron.Pt() - recoElectron.Pt())/genElectron.Pt()
            prof_electron_deltaPt_vs_Lxy.Fill(Lxy, deltaPt) 
            
            if Lxy < 1.0: histo_deltaPt_Lxy_0_1.Fill(deltaPt)
            elif Lxy > 1.0 and Lxy < 2.0: histo_deltaPt_Lxy_1_2.Fill(deltaPt)
            elif Lxy > 2.0 and Lxy < 4.0: histo_deltaPt_Lxy_2_4.Fill(deltaPt)
            elif Lxy > 4.0 and Lxy < 8.0: histo_deltaPt_Lxy_4_8.Fill(deltaPt)
            elif Lxy > 8.0: histo_deltaPt_Lxy_8_inf.Fill(deltaPt)



        eff_electron_Lxy.Fill(reconstructed, Lxy)


    ### Single Electron matching fake rate
    for er in range(0, t.nElectronCandidate):

        fake = not candidate_status[er]
        
        reco_Lxy = math.sqrt((t.IsoTrackSel_vx[t.ElectronCandidate_isotrackIdx[er]] - t.PV_vx)**2 + (t.IsoTrackSel_vy[t.ElectronCandidate_isotrackIdx[er]] - t.PV_vy)**2)
        fake_electron_Lxy.Fill(fake, reco_Lxy)


    ### Dilepton pairs efficiency 
    # Defined as the pairs of genElectrons that have been fully matched

    # pair variables definition
    pair0 = []
    pair1 = []

    for g in gen_ei:

        if t.GenLeptonSel_motherIdx[g] == 0:
            pair0.append(g)     
        if t.GenLeptonSel_motherIdx[g] == 1:
            pair1.append(g)


    if len(pair0) == 2: # First pair is a dielectron candidate

        Lxy = math.sqrt((t.GenLeptonSel_vx[pair0[0]] - t.PV_vx)**2 + (t.GenLeptonSel_vy[pair0[0]] - t.PV_vy)**2)

        if GenLepton_status[pair0[0]] and GenLepton_status[pair0[1]]:
            eff_dielectron_Lxy.Fill(True, Lxy)
        else:
            eff_dielectron_Lxy.Fill(False, Lxy)

    if len(pair1) == 2: # Second pair is a dielectron candidate

        Lxy = math.sqrt((t.GenLeptonSel_vx[pair1[0]] - t.PV_vx)**2 + (t.GenLeptonSel_vy[pair1[0]] - t.PV_vy)**2)

        if GenLepton_status[pair1[0]] and GenLepton_status[pair1[1]]:
            eff_dielectron_Lxy.Fill(True, Lxy)
        else:
            eff_dielectron_Lxy.Fill(False, Lxy)



##############################################################################################
############################### Save in root file ############################################
##############################################################################################

efficiency_plots = []
efficiency_plots.append(eff_electron_Lxy)
efficiency_plots.append(eff_dielectron_Lxy)

profile_plots = []
profile_plots.append(prof_electron_deltaPt_vs_Lxy)

histos_deltaPt_Lxy = []
histos_deltaPt_Lxy.append(histo_deltaPt_Lxy_0_1)
histos_deltaPt_Lxy.append(histo_deltaPt_Lxy_1_2)
histos_deltaPt_Lxy.append(histo_deltaPt_Lxy_2_4)
histos_deltaPt_Lxy.append(histo_deltaPt_Lxy_4_8)
histos_deltaPt_Lxy.append(histo_deltaPt_Lxy_8_inf)



output_file = TFile("dielectron_effPlots.root", "RECREATE")

for plot in efficiency_plots + profile_plots:
    plot.Write()

output_file.Close()

##############################################################################################
################################# Plot in canvas #############################################
##############################################################################################


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
#r.TGaxis.SetMaxDigits(4)

color_set = [r.kRed, r.kBlue, r.kGreen, r.kMagenta, r.kBlack]


###### Plot the distributions

for eff in efficiency_plots:

    c1.Clear()
    c1.SetLogy(0)
    xlimit = tuneEfficiency(c1, eff)
    c1.Update()
    latexCMS.Draw()
    latexCMSExtra.Draw()

    c1.SaveAs(output_dir+eff.GetName()+'.png')


for prof in profile_plots:

    c1.Clear()
    c1.SetLogy(0)
    prof.SetMarkerStyle(21)
    prof.SetMarkerColor(r.kBlack)
    prof.SetMarkerSize(1.)
    prof.Draw('P')
    c1.Update()
    
    latexCMS.Draw()
    latexCMSExtra.Draw()

    c1.SaveAs(output_dir+prof.GetName()+'.png')


###### Plot deltaPt Lxy histograms

c1.Clear()
c1.SetLogy(0)

tags = ['L_{xy}: [0, 1]',
         'L_{xy}: [1, 2]',
         'L_{xy}: [2, 4]',
         'L_{xy}: [4, 8]',
         'L_{xy}: [8, inf]']


for n in range(0, len(histos_deltaPt_Lxy)):

    tuneEmptyHisto(histos_deltaPt_Lxy[n], '(p_{T}^{gen} - p_{T}^{e})/p_{T}^{gen}', color_set[n], normed = True)

printJointHistos(histos_deltaPt_Lxy, log = False)

leg = drawLegend(histos_deltaPt_Lxy, tags)
leg.Draw()
latexCMS.Draw()
latexCMSExtra.Draw()

c1.Update()
c1.SaveAs(output_dir+'histos_deltaPt_Lxy.png')
