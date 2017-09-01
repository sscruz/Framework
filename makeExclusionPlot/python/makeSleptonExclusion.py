#####################################################################
######                                                              #
###### 88888888888         88                        888888888888   #  
###### 88                  88                                 ,88   #
###### 88                  88                               ,88"    #  
###### 88aaaaa     ,adPPYb,88  ,adPPYb,d8  ,adPPYba,      ,88"      #
###### 88"""""    a8"    `Y88 a8"    `Y88 a8P_____88    ,88"        #
###### 88         8b       88 8b       88 8PP"""""""  ,88"          #
###### 88         "8a,   ,d88 "8a,   ,d88 "8b,   ,aa 88"            #
###### 88888888888 `"8bbdP"Y8  `"YbbdP"Y8  `"Ybbd8"' 888888888888   #
######                       aa,    ,88                             #
######                         "Y8bbdP"                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraph, TGraphErrors, TGraphAsymmErrors, SetOwnership
import math, sys, optparse, array, copy
import gc, inspect
from array import array
import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample



class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'



def xsection(m, xs, mneu, re):
    for i in range(0, len(m)):
        if m[i] == mneu:
            return xs[i] * re
    return 0
    

if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting 1D plot...                          '
    print '#######################################################################' + bcolors.ENDC
    lumi = 35.9 ; maxrun = 999999; lumi_str = '35.9invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()
    m = []
    xs = []
    exs = []
    eys = []
    f = open("datacards/sleptonXsec.txt")
    xsecs = eval(f.read())

    for i in f.readlines():
        s = i.split()
        m.append(float(s[0]))
        xs.append(float(s[1]))
        exs.append(0)
        eys.append(float(s[2]))


    mval = []
    rexp = []
    rexp1m = []
    rexp2m = []
    rexp1p = []
    rexp2p = []
    robs = []
    rfile = r.TFile("limits_slepton2017.root")
    #rfile = r.TFile("limits_ChiZZ_Moriond2017.root")
    hExp = rfile.Get("ex_exp")
    hObs = rfile.Get("ex_obs")
    hExp1m = rfile.Get("ex_exp_m1s")
    hExp2m = rfile.Get("ex_exp_m2s")
    hExp1p = rfile.Get("ex_exp_p1s")
    hExp2p = rfile.Get("ex_exp_p2s")
    for binx in range(1, hExp.GetNbinsX()):
        mneu = hExp.GetXaxis().GetBinCenter(binx)
        mneu = hExp.GetXaxis().GetBinCenter(binx)
        bin = hExp.FindBin(mneu, 5)
        re = hExp.GetBinContent(bin)
        reobs = hObs.GetBinContent(bin)
        re1p = hExp1p.GetBinContent(bin)
        re2p = hExp2p.GetBinContent(bin)
        re1m = hExp1m.GetBinContent(bin)
        re2m = hExp2m.GetBinContent(bin)
        xsec    = xsecs[mneu][0] * re
        xsec1m  = xsecs[mneu][0] * re1m 
        xsec2m  = xsecs[mneu][0] * re2m 
        xsec1p  = xsecs[mneu][0] * re1p 
        xsec2p  = xsecs[mneu][0] * re2p 
        xsecobs = xsecs[mneu][0] * reobs
        if re != 0:
            print mneu, reobs, xsecobs 
            xs.append(xsecs[mneu][0])
            m .append(mneu)
            mval.append(mneu)
            rexp.append(xsec)
            robs.append(xsecobs)
            rexp1p.append(xsec1p)
            rexp2p.append(xsec2p)
            rexp1m.append(xsec1m)
            rexp2m.append(xsec2m)
            exs.append(0.)
            eys.append(xsecs[mneu][1])
        
    line1s = []
    eupline1s = []
    edownline1s = []
    line2s = []
    eupline2s = []
    edownline2s = []
    for i in range(0, len(rexp)):
        line1s.append(rexp[i])
        line2s.append(rexp[i])
        eupline1s.append(rexp1p[i] - rexp[i])    
        edownline1s.append(-rexp1m[i] + rexp[i])    
        eupline2s.append(rexp2p[i] - rexp[i])    
        edownline2s.append(-rexp2m[i] + rexp[i])    

    print "len(line1s)", len(line1s)
    print "array('d', mval)", array('d', mval)
    print "array('d', robs)", array('d', robs)
    graphobs = TGraph(len(line1s), array('d', mval), array('d', robs))
    graphexp = TGraphAsymmErrors(len(line1s), array('d', mval), array('d', line1s), array('d', exs), array('d', exs), array('d', edownline1s), array('d', eupline1s))
    graphexp1s = TGraphAsymmErrors(len(line1s), array('d', mval), array('d', line1s), array('d', exs), array('d', exs), array('d', edownline1s), array('d', eupline1s))
    graphexp2s = TGraphAsymmErrors(len(line2s), array('d', mval), array('d', line2s), array('d', exs), array('d', exs), array('d', edownline2s), array('d', eupline2s))
    graphexp1s.SetFillColor(r.kGreen)
    graphexp.SetLineColor(r.kBlack)
    graphexp.SetLineWidth(3)
    graphobs.SetLineColor(r.kBlack)
    graphobs.SetLineWidth(3)
    graphexp.SetLineStyle(2)
    
    label = "m(#tilde{#chi}^{0}_{1}) [GeV]"

    graphexp2s.SetFillColor(r.kYellow)
    graphexp1s.GetYaxis().SetTitle("#sigma [fb]")    
    graphexp2s.GetYaxis().SetTitle("#sigma [fb]")    
    graphexp.GetYaxis().SetTitle("#sigma [fb]")    
    graphobs.GetYaxis().SetTitle("#sigma [fb]")    
    graphexp1s.GetXaxis().SetTitle(label)    
    graphexp2s.GetXaxis().SetTitle(label)    
    graphexp.GetXaxis().SetTitle(label)    
    graphobs.GetXaxis().SetTitle(label)    
    
  

    graphtheory = TGraphErrors(len(xs), array('d', m), array('d', xs), array('d', exs), array('d', eys))
    graphtheory.SetFillColor(r.kBlue-3)
    graphtheory.GetYaxis().SetTitle("#sigma [fb]")    
    graphtheory.GetXaxis().SetTitle(label)    

    xmax = 450

    graphobs.GetXaxis().SetLimits(0, xmax)    
    graphexp.GetXaxis().SetLimits(0, xmax)    
    graphexp1s.GetXaxis().SetLimits(0, xmax)    
    graphexp2s.GetXaxis().SetLimits(0, xmax)    
    graphtheory.GetXaxis().SetLimits(0, xmax)    

    aux = r.TH1F("aux", "", 1, 100, xmax)
    aux.GetXaxis().SetTitle(label)
    aux.GetYaxis().SetTitle("#sigma [fb]")
    aux.GetYaxis().SetRangeUser(1, 300000)
    
    plot = Canvas.Canvas('TSlepton', 'png,pdf,root', 0.45, 0.5, 0.65, 0.7)
    plot.addHisto(aux, "H", "lab", "F", r.kBlack, 1, -1)
    plot.addGraph(graphexp2s, "Z3", "#pm 2 s.d. expected range", "F", r.kBlack, 1, 3)
    plot.addGraph(graphexp1s, "L3", "#pm 1 s.d. expected range", "F", r.kBlack, 1, 2)
    plot.addGraph(graphtheory, "Z3", "Theory", "F", r.kBlue, 1, 4)
    plot.addGraph(graphexp, "LX", "Expected 95% CL", "L", r.kBlack, 1, 1)
    plot.addGraph(graphobs, "LX", "Observed 95% CL", "L", r.kBlack, 1, 0)
    #plot.redrawAxis()
    plot.addLatex(0.2, 0.85, "pp #rightarrow #tilde{l}#tilde{l}")
    #plot.addLatex(0.2, 0.85, "pp #rightarrow #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1}")
    plot.addLatex(0.2, 0.78, "m(#tilde{#chi}^{0}_{1}) = 5 ")
    plot.addHisto(aux, "sameaxis", "lab", "F", r.kBlack, 1, -1)

    plot.save(1, 1, 1, lumi)

    outputfile = r.TFile("SUS_17_slepton.root", "recreate")
    #outputfile = r.TFile("SUS-PAS-16-0-21_Fig10.root", "recreate")
    outputfile.cd()
    graphexp2s.SetName("Two-sigma-expected-range")
    graphexp2s.SetTitle("Two-sigma-expected-range")
    graphexp1s.SetName("One-sigma-expected-range")
    graphexp1s.SetTitle("One-sigma-expected-range")
    graphexp.SetName("Expected")
    graphexp.SetTitle("Expected")
    graphobs.SetName("Observed")
    graphobs.SetTitle("Observed")
    graphtheory.SetName("Theory")
    graphtheory.SetTitle("Theory")
    graphexp2s.Write()
    graphexp1s.Write()
    graphexp.Write()
    graphobs.Write()
    graphtheory.Write()
    outputfile.Close()



