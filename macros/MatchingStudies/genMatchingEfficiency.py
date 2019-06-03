import ROOT as r
from ROOT import gROOT, TCanvas, TFile, TGraphErrors, TLatex
from plotTools import *
import os



############################### Open the files #################################

file_name = '/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/IFCALongLivedAnalysis/output.root'

File = TFile(file_name)
t = File.Get("Events")

print('Open the file: ' + file_name)
print("'Events' tree with " + str(t.GetEntries()) + ' entries \n')



############################### Output definition #################################

output_dir = '/eos/user/f/fernance/www/LLP/plots/'
if not os.path.exists(output_dir): os.mkdir('./'+output_dir)

############################### Histogram definition ##################################

#### -> Isolated histos
histo_nElectronCandidate = r.TH1F("histo_nElectronCandidate", "", 8, 0, 8)
histo_nMuonCandidate = r.TH1F("histo_nMuonCandidate", "", 8, 0, 8)

#### -> deltaR(gen-tr) for different dxy
histo_deltaRgt_dxy_0_1 = r.TH1F("histo_deltaRgt_dxy_0_1", "", 30, 0, 3)
histo_deltaRgt_dxy_1_2 = r.TH1F("histo_deltaRgt_dxy_1_2", "", 30, 0, 3)
histo_deltaRgt_dxy_2_4 = r.TH1F("histo_deltaRgt_dxy_2_4", "", 30, 0, 3)
histo_deltaRgt_dxy_4_8 = r.TH1F("histo_deltaRgt_dxy_4_8", "", 30, 0, 3)
histo_deltaRgt_dxy_8_inf = r.TH1F("histo_deltaRgt_dxy_8_inf", "", 30, 0, 3)

#### -> deltaR(tr-sc) for different dxy
histo_deltaRtsc_dxy_0_1 = r.TH1F("histo_deltaRtsc_dxy_0_1", "", 30, 0, 3)
histo_deltaRtsc_dxy_1_2 = r.TH1F("histo_deltaRtsc_dxy_1_2", "", 30, 0, 3)
histo_deltaRtsc_dxy_2_4 = r.TH1F("histo_deltaRtsc_dxy_2_4", "", 30, 0, 3)
histo_deltaRtsc_dxy_4_8 = r.TH1F("histo_deltaRtsc_dxy_4_8", "", 30, 0, 3)
histo_deltaRtsc_dxy_8_inf = r.TH1F("histo_deltaRtsc_dxy_8_inf", "", 30, 0, 3)

#### -> deltaR(tr-mu) for different dxy
histo_deltaRtmu_dxy_0_1 = r.TH1F("histo_deltaRtmu_dxy_0_1", "", 30, 0, 3)
histo_deltaRtmu_dxy_1_2 = r.TH1F("histo_deltaRtmu_dxy_1_2", "", 30, 0, 3)
histo_deltaRtmu_dxy_2_4 = r.TH1F("histo_deltaRtmu_dxy_2_4", "", 30, 0, 3)
histo_deltaRtmu_dxy_4_8 = r.TH1F("histo_deltaRtmu_dxy_4_8", "", 30, 0, 3)
histo_deltaRtmu_dxy_8_inf = r.TH1F("histo_deltaRtmu_dxy_8_inf", "", 30, 0, 3)


#### -> Profile
prof_deltaRtsc_vs_dxy = r.TProfile("prof_deltaRtsc_vs_dxy", "", 20, 0, 100)
prof_deltaRtmu_vs_dxy = r.TProfile("prof_deltaRtmu_vs_dxy", "", 20, 0, 100)



############################### Loop over the tree entries ##################################

for n in range(0, t.GetEntries()):

    t.GetEntry(n)


    histo_nElectronCandidate.Fill(t.nElectronCandidate)
    histo_nMuonCandidate.Fill(t.nMuonCandidate)


    for g in range(0, t.nGenLepton):

        if t.GenLeptonSel_trackMatch[g] == 99: continue

        dR = t.GenLeptonSel_trackdR[g]
        dxy = abs(t.GenLeptonSel_dxy[g])

        if (dxy <= 1): histo_deltaRgt_dxy_0_1.Fill(dR)
        elif (dxy > 1 and dxy <= 2): histo_deltaRgt_dxy_1_2.Fill(dR)
        elif (dxy > 2 and dxy <= 4): histo_deltaRgt_dxy_2_4.Fill(dR)
        elif (dxy > 4 and dxy <= 8): histo_deltaRgt_dxy_4_8.Fill(dR)
        else: histo_deltaRgt_dxy_8_inf.Fill(dR)


        if t.GenLeptonSel_objectMatch[g] == 99: continue

        if dR <= 0.1:

            if (abs(t.GenLeptonSel_pdgId[g]) == 11):
            
                if (dxy <= 1): histo_deltaRtsc_dxy_0_1.Fill(t.GenLeptonSel_pairdR[g])
                elif (dxy > 1 and dxy <= 2): histo_deltaRtsc_dxy_1_2.Fill(t.GenLeptonSel_pairdR[g])
                elif (dxy > 2 and dxy <= 4): histo_deltaRtsc_dxy_2_4.Fill(t.GenLeptonSel_pairdR[g])
                elif (dxy > 4 and dxy <= 8): histo_deltaRtsc_dxy_4_8.Fill(t.GenLeptonSel_pairdR[g])
                else: histo_deltaRtsc_dxy_8_inf.Fill(t.GenLeptonSel_pairdR[g])

            if (abs(t.GenLeptonSel_pdgId[g]) == 13):
            
                if (dxy <= 1): histo_deltaRtmu_dxy_0_1.Fill(t.GenLeptonSel_pairdR[g])
                elif (dxy > 1 and dxy <= 2): histo_deltaRtmu_dxy_1_2.Fill(t.GenLeptonSel_pairdR[g])
                elif (dxy > 2 and dxy <= 4): histo_deltaRtmu_dxy_2_4.Fill(t.GenLeptonSel_pairdR[g])
                elif (dxy > 4 and dxy <= 8): histo_deltaRtmu_dxy_4_8.Fill(t.GenLeptonSel_pairdR[g])
                else: histo_deltaRtmu_dxy_8_inf.Fill(t.GenLeptonSel_pairdR[g])


        if abs(t.GenLeptonSel_pdgId[g]) == 11: prof_deltaRtsc_vs_dxy.Fill(dxy, t.GenLeptonSel_objectdR[g])
        if abs(t.GenLeptonSel_pdgId[g]) == 13: prof_deltaRtmu_vs_dxy.Fill(dxy, t.GenLeptonSel_objectdR[g])


############################### Write the histos in a file ##################################

TList = r.TObjArray(0)
TList.Add(histo_deltaRgt_dxy_0_1)
TList.Add(histo_deltaRgt_dxy_1_2)
TList.Add(histo_deltaRgt_dxy_2_4)
TList.Add(histo_deltaRgt_dxy_4_8)
TList.Add(histo_deltaRgt_dxy_8_inf)

TList.Add(prof_deltaRtsc_vs_dxy)
TList.Add(prof_deltaRtmu_vs_dxy)


output_file = TFile("output_histos.root", "RECREATE")
TList.Write()


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

color_set = [r.kRed, r.kBlue, r.kGreen, r.kMagenta, r.kBlack]

###### Plot some distributions
c1.Clear()
c1.SetLogy(0)

tuneEmptyHisto(histo_nMuonCandidate, 'nMuonCandidate', r.kRed, normed = False)
tuneEmptyHisto(histo_nElectronCandidate, 'nElectronCandidate', r.kBlue, normed = False)

printJointHistos([histo_nMuonCandidate], ['nMuonCandidate'], log = False)

latexCMS.Draw()
latexCMSExtra.Draw()
c1.SaveAs(output_dir+'histo_nMuonCandidate.png')


c1.Clear()
c1.SetLogy(0)
printJointHistos([histo_nElectronCandidate], ['nElectronCandidate'], log = False)

latexCMS.Draw()
latexCMSExtra.Draw()
c1.SaveAs(output_dir+'histo_nElectronCandidate.png')

###### Plot deltaR(gen-track) distributions for dxy ranges:

c1.Clear()
c1.SetLogy(1)

tuneEmptyHisto(histo_deltaRgt_dxy_0_1, '#DeltaR (generation - track)  ', r.kRed, normed = True)
tuneEmptyHisto(histo_deltaRgt_dxy_1_2, '#DeltaR (generation - track)  ', r.kBlue, normed = True)
tuneEmptyHisto(histo_deltaRgt_dxy_2_4, '#DeltaR (generation - track)  ', r.kGreen, normed = True)
tuneEmptyHisto(histo_deltaRgt_dxy_4_8, '#DeltaR (generation - track)  ', r.kMagenta, normed = True)
tuneEmptyHisto(histo_deltaRgt_dxy_8_inf, '#DeltaR (generation - track)  ', r.kBlack, normed = True)

list_deltaRgt_dxy = [histo_deltaRgt_dxy_0_1,
                     histo_deltaRgt_dxy_1_2,
                     histo_deltaRgt_dxy_2_4,
                     histo_deltaRgt_dxy_4_8,
                     histo_deltaRgt_dxy_8_inf]

tags_deltaRgt_dxy = ['d_{xy}: [0, 1]',
                     'd_{xy}: [1, 2]',
                     'd_{xy}: [2, 4]',
                     'd_{xy}: [4, 8]',
                     'd_{xy}: [8, inf]']


printJointHistos(list_deltaRgt_dxy, tags_deltaRgt_dxy, log = True)

leg = drawLegend(list_deltaRgt_dxy, tags_deltaRgt_dxy)
leg.Draw()
latexCMS.Draw()
latexCMSExtra.Draw()

c1.Update()
c1.SaveAs(output_dir+'histo_deltaRgt_dxy.png')




###### Plot deltaR(track-sc) distributions for dxy ranges with deltaR(gen-track) < 0.1:

c1.Clear()
c1.SetLogy(1)


list_deltaRtsc_dxy = [histo_deltaRtsc_dxy_0_1,
                     histo_deltaRtsc_dxy_1_2,
                     histo_deltaRtsc_dxy_2_4,
                     histo_deltaRtsc_dxy_4_8,
                     histo_deltaRtsc_dxy_8_inf]

tags_deltaRtsc_dxy = ['d_{xy}: [0, 1]',
                     'd_{xy}: [1, 2]',
                     'd_{xy}: [2, 4]',
                     'd_{xy}: [4, 8]',
                     'd_{xy}: [8, inf]']

for n in range(0, len(list_deltaRtsc_dxy)):

    tuneEmptyHisto(list_deltaRtsc_dxy[n], '#DeltaR (track - supercluster)  ', color_set[n], normed = True)

printJointHistos(list_deltaRtsc_dxy, tags_deltaRtsc_dxy, log = True)

leg = drawLegend(list_deltaRtsc_dxy, tags_deltaRtsc_dxy)
leg.Draw()
latexCMS.Draw()
latexCMSExtra.Draw()

c1.Update()
c1.SaveAs(output_dir+'histo_deltaRtsc_dxy.png')



###### Plot deltaR(track-mu) distributions for dxy ranges with deltaR(gen-track) < 0.1:

c1.Clear()
c1.SetLogy(1)


list_deltaRtmu_dxy = [histo_deltaRtmu_dxy_0_1,
                     histo_deltaRtmu_dxy_1_2,
                     histo_deltaRtmu_dxy_2_4,
                     histo_deltaRtmu_dxy_4_8,
                     histo_deltaRtmu_dxy_8_inf]

tags_deltaRtmu_dxy = ['d_{xy}: [0, 1]',
                     'd_{xy}: [1, 2]',
                     'd_{xy}: [2, 4]',
                     'd_{xy}: [4, 8]',
                     'd_{xy}: [8, inf]']

for n in range(0, len(list_deltaRtmu_dxy)):

    tuneEmptyHisto(list_deltaRtmu_dxy[n], '#DeltaR (track - muon)  ', color_set[n], normed = True)

printJointHistos(list_deltaRtmu_dxy, tags_deltaRtmu_dxy, log = True)

leg = drawLegend(list_deltaRtmu_dxy, tags_deltaRtmu_dxy)
leg.Draw()
latexCMS.Draw()
latexCMSExtra.Draw()

c1.Update()
c1.SaveAs(output_dir+'histo_deltaRtmu_dxy.png')

