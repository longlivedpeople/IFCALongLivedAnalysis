import ROOT as r
from ROOT import TLatex



def tuneEmptyHisto(histo, xlabel, color, normed):

    histo.GetXaxis().SetTitleOffset(1.07)
    histo.GetYaxis().SetTitleOffset(1.3)
    
    histo.GetXaxis().SetTitleSize(0.04)
    histo.GetYaxis().SetTitleSize(0.04)

    # Assumed to be in format 'magnitude (units)'
    histo.GetXaxis().SetTitle(xlabel)

    width = (histo.GetXaxis().GetXmax() - histo.GetXaxis().GetXmin())/histo.GetNbinsX()

    units = xlabel.split()[-1]
    
    if '(' in units: units = units[1:-1]
    else: units = ''
    histo.GetYaxis().SetTitle("Events / %.1f"%width + " " + units )


    histo.SetLineWidth(1)
    histo.SetLineColor(color)


    if normed: histo.Scale(1/histo.GetEntries())



def printJointHistos(histo_list, histo_tags, log = False):


    ymax = 0

    for histo in histo_list:

        y = histo.GetMaximum()
        if y > ymax: ymax = y

    if log: histo_list[0].SetMaximum(10*ymax)
    else: histo_list[0].SetMaximum(1.3*ymax)

    histo_list[0].Draw('hist')


    for n in range(1, len(histo_list)):

        histo_list[n].Draw('same hist')



def drawLegend(histo_list, histo_tags):

    # Set legend dimensions
    yi = 0.88 - 0.04*len(histo_list)
    
    maxlen = 0
    for tag in histo_tags:
        if len(tag) > maxlen: maxlen = len(tag)

    xi = 0.85 - 0.01*maxlen


    leg = r.TLegend(xi, yi, 0.85, 0.88)
    leg.SetShadowColor(0)
    leg.SetFillColor(0)
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)

    for n in range(0, len(histo_list)):

        leg.AddEntry(histo_list[n], histo_tags[n], 'l')



    return leg


def writeCMS(simulation = True):

    latexCMS = TLatex(0.15, 0.93, "CMS")
    latexCMS.SetTextFont(61)
    latexCMS.SetTextAlign(11)
    latexCMS.SetTextSize(0.06)
    latexCMS.SetNDC(True)

    if simulation: 
        CMSExtra = 'Simulation'
        yextra = 0.81
    else:
        CMSExtra = 'Preliminary'
        yextra = 0.84

    latexCMSExtra = TLatex(0.28, 0.93, CMSExtra)
    latexCMSExtra.SetTextFont(52)
    latexCMSExtra.SetTextAlign(11)
    latexCMSExtra.SetTextSize(0.045)
    latexCMSExtra.SetNDC(True)

    return latexCMS, latexCMSExtra

