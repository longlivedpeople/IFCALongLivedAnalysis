import ROOT as r
from ROOT import gROOT, TCanvas, TFile, TGraphErrors, TLatex
from plotTools import *
import numpy as np



############################### Open the files #################################

file_name = '/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/IFCALongLivedAnalysis/output.root'

File = TFile(file_name)
t = File.Get("Events")

print('Open the file: ' + file_name)
print("'Events' tree with " + str(t.GetEntries()) + ' entries \n')



############################### Output definition #################################

output_dir = 'plots/'


############################### Histogram definition ##################################

dxy_bin = np.array([0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20., 30., 40., 50., 60., 70., 90., 110.0])



#### -> Efficiency

eff_dxy_ele = r.TEfficiency("eff_dxy_ele",";Generated transverse impact parameter |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)
eff_dxy_mu = r.TEfficiency("eff_dxy_mu",";Generated transverse impact parameter |d_{xy}| (cm);Efficiency", len(dxy_bin) -1 , dxy_bin)


#### -> Fake rate
fake_dxy_ele = r.TEfficiency("fake_dxy_ele",";Generated transverse impact parameter |d_{xy}| (cm);Fake rate", len(dxy_bin) -1 , dxy_bin)
fake_dxy_mu = r.TEfficiency("fake_dxy_mu",";Generated transverse impact parameter |d_{xy}| (cm);Fake rate", len(dxy_bin) -1 , dxy_bin)



############################### Loop over the tree entries ##################################

for n in range(0, t.GetEntries()):

    t.GetEntry(n)

    for g in range(0, t.nGenLepton):

        reconstructed = False

        if (t.GenLeptonSel_trackMatch == 99 or t.GenLeptonSel_objectMatch == 99): continue
        if (t.GenLeptonSel_trackdR[g] > 0.1): continue


        gt = t.GenLeptonSel_trackMatch[g]
        go = t.GenLeptonSel_objectMatch[g]

        if (abs(t.GenLeptonSel_pdgId[g]) == 11):


            for e in range(0, t.nElectronCandidate):
                rt = t.ElectronCandidate_isotrackIdx[e]
                ro = t.ElectronCandidate_photonIdx[e]

                if gt == rt and go == ro:

                    reconstructed = True
                    break

            eff_dxy_ele.Fill(reconstructed, t.GenLeptonSel_dxy[g])

        if (abs(t.GenLeptonSel_pdgId[g]) == 13):

            for m in range(0, t.nMuonCandidate):

                rt = t.MuonCandidate_isotrackIdx[m]
                ro = t.MuonCandidate_muonTriggerObjectIdx[m]

                if gt == rt and go == ro:

                    reconstructed = True
                    break


            eff_dxy_mu.Fill(reconstructed, t.GenLeptonSel_dxy[g])

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


###### Plotting:

efficiencies = [eff_dxy_ele,
                eff_dxy_mu]


eff_name = ['eff_dxy_ele',
            'eff_dxy_mu']



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


    c1.SaveAs(eff.GetName()+'.png')



