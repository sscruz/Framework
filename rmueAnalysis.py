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
import math as math
import Rounder as rounder


from optparse import OptionParser
from ROOT import gROOT, TCanvas, TFile
from Sample import Sample, Block, Tree
from CutManager import CutManager
from Canvas import Canvas
from ROOT import TGraphErrors
from array import array

def rmue(Nmm, Nee, Emm, Eee):

    val = [0, 0]
    if(Nee != 0):
        val[0] = math.sqrt(Nmm/Nee)
        val[1] = math.sqrt(0.25 * Emm * Emm / (Nmm * Nee) + 0.25 * Eee * Eee * Nmm / (Nee * Nee * Nee))

    return val


def make_rmue(histo_mm, histo_ee):

    ratio = histo_mm.Clone("rmue_" + histo_mm.GetName())
    ratio.GetYaxis().SetTitle("r_{#mu e}")
  
    for i in range(0, histo_mm.GetNbinsX()+1):
        Nmm = histo_mm.GetBinContent(i)
        Nee = histo_ee.GetBinContent(i)
        Emm = histo_mm.GetBinError(i)
        Eee = histo_ee.GetBinError(i)
        if(Nee != 0):
            val = rmue(Nmm, Nee, Emm, Eee)
            ratio.SetBinContent(i, val[0])
            ratio.SetBinError(i, val[1])
    return ratio


def make_rmue_signal(ee_l, mm_l, ee_Z, mm_Z, ee_h, mm_h):

    rmue_l = rmue(mm_l[0], ee_l[0], mm_l[1], ee_l[1])
    rmue_Z = rmue(mm_Z[0], ee_Z[0], mm_Z[1], ee_Z[1])
    rmue_h = rmue(mm_h[0], ee_h[0], mm_h[1], ee_h[1])

    rmue_x = array("d", [(20+70)/2.0, (101.0+81.0)/2.0, (110+300)/2.0])
    rmue_ex = array("d", [25, 5, 100])
    rmue_y = array("d", [rmue_l[0], rmue_Z[0], rmue_h[0]])
    rmue_ey = array("d", [rmue_l[1], rmue_Z[1], rmue_h[1]])
    
    print "r_mue in low mass signal region " + str(rmue_l[0]) + " +/- " + str(rmue_l[1])
    print "r_mue in onZ mass signal region " + str(rmue_Z[0]) + " +/- " + str(rmue_Z[1])
    print "r_mue in high mass signal region " + str(rmue_h[0]) + " +/- " + str(rmue_h[1])

    gr = TGraphErrors(3, rmue_x, rmue_y, rmue_ex, rmue_ey)
    gr.GetYaxis().SetTitle("r_{#mu e}")
    gr.GetXaxis().SetTitle("m_{ll} [GeV]")

    return [gr, rmue_l, rmue_Z, rmue_h]

def make_rmue_meas(ee, mm, syst):

    rmue_m = rmue(mm[0], ee[0], mm[1], ee[1])

    rmue_x = array("d", [(300+20)/2.0])
    rmue_ex = array("d", [140])
    rmue_x2 = array("d", [(100)/2.0])
    rmue_ex2 = array("d", [50])
    rmue_y = array("d", [rmue_m[0]])
    error = math.sqrt(rmue_m[1]*rmue_m[1] + (syst*rmue_m[0])*(syst*rmue_m[0]))
    rmue_ey = array("d", [error])
    
    print "r_mue in DY region (mean) " + str(rmue_m[0]) + " +/- " + str(rmue_m[1]) + " +/- " + str(syst*rmue_m[0]) 

    
    gr = TGraphErrors(1, rmue_x, rmue_y, rmue_ex, rmue_ey)
    gr.GetYaxis().SetTitle("r_{#mu e}")
    gr.GetXaxis().SetTitle("m_{ll} [GeV]")
    gr2 = TGraphErrors(1, rmue_x2, rmue_y, rmue_ex2, rmue_ey)
    gr2.GetYaxis().SetTitle("r_{#mu e}")
    gr2.GetXaxis().SetTitle("MET [GeV]")
    
    return [gr, gr2, rmue_m]

if __name__ == "__main__":

    parser = OptionParser(usage="usage: %prog [options] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-m", "--mode", action="store", dest="mode", default="rmue", help="Operation mode")
    (options, args) = parser.parse_args()

    if len(args) != 1:
      parser.error("wrong number of arguments")

    inputFileName = args[0]
    tree = Tree(inputFileName, "MC", 0)
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager()

    bins = [20,30,40, 50, 60, 70, 81, 101, 120, 150, 180, 220, 260, 300]
    #####Attention: setting the first bin in 50 because we only have DY mll>50 at this point
    #mll_ee_central = tree.getTH1F(4, "mll_ee_central", "t.lepsMll_Edge", 25, 50, 300, cuts.Add(cuts.DYControlNoMassLeptonee(), cuts.Central()), "", "m_{ll} [GeV]")
    #mll_mm_central = tree.getTH1F(4, "mll_mm_central", "t.lepsMll_Edge", 25, 50, 300, cuts.Add(cuts.DYControlNoMassLeptonmm(), cuts.Central()), "", "m_{ll} [GeV]")
    mll_ee_central = tree.getTH1F(4, "mll_ee_central", "t.lepsMll_Edge", bins, 1, 1, cuts.Add(cuts.DYControlNoMassLeptonee(), cuts.Central()), "", "m_{ll} [GeV]")
    mll_mm_central = tree.getTH1F(4, "mll_mm_central", "t.lepsMll_Edge", bins, 1, 1,cuts.Add(cuts.DYControlNoMassLeptonmm(), cuts.Central()), "", "m_{ll} [GeV]")
    met_ee_central = tree.getTH1F(4, "met_ee_central", "met_pt", 10, 0, 100, cuts.Add(cuts.DYControlNoMetLeptonee(), cuts.Central()), "", "m_{ll} [GeV]")
    met_mm_central = tree.getTH1F(4, "met_mm_central", "met_pt", 10, 0, 100, cuts.Add(cuts.DYControlNoMetLeptonmm(), cuts.Central()), "", "m_{ll} [GeV]")
    mll_ee_central_lowmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalLowMassee(), cuts.Central()))
    mll_mm_central_lowmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalLowMassmm(), cuts.Central()))
    mll_ee_central_highmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalHighMassee(), cuts.Central()))
    mll_mm_central_highmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalHighMassmm(), cuts.Central()))
    mll_ee_central_Zmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalZMassee(), cuts.Central()))
    mll_mm_central_Zmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalZMassmm(), cuts.Central()))
    mll_ee_central_DYmeas = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.AddList([cuts.DYControlNoMassLeptonee(), cuts.Central(), cuts.DYMass()]))
    mll_mm_central_DYmeas = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.AddList([cuts.DYControlNoMassLeptonmm(), cuts.Central(), cuts.DYMass()]))
    
    #mll_ee_forward = tree.getTH1F(4, "mll_ee_forward", "t.lepsMll_Edge", 25, 50, 300, cuts.Add(cuts.DYControlNoMassLeptonee(), cuts.Forward()), "", "m_{ll} [GeV]")
    #mll_mm_forward = tree.getTH1F(4, "mll_mm_forward", "t.lepsMll_Edge", 25, 50, 300, cuts.Add(cuts.DYControlNoMassLeptonmm(), cuts.Forward()), "", "m_{ll} [GeV]")
    mll_ee_forward = tree.getTH1F(4, "mll_ee_forward", "t.lepsMll_Edge", bins, 1, 1, cuts.Add(cuts.DYControlNoMassLeptonee(), cuts.Forward()), "", "m_{ll} [GeV]")
    mll_mm_forward = tree.getTH1F(4, "mll_mm_forward", "t.lepsMll_Edge", bins, 1, 1, cuts.Add(cuts.DYControlNoMassLeptonmm(), cuts.Forward()), "", "m_{ll} [GeV]")
    met_ee_forward = tree.getTH1F(4, "met_ee_forward", "met_pt", 10, 0, 100, cuts.Add(cuts.DYControlNoMetLeptonee(), cuts.Forward()), "", "m_{ll} [GeV]")
    met_mm_forward = tree.getTH1F(4, "met_mm_forward", "met_pt", 10, 0, 100, cuts.Add(cuts.DYControlNoMetLeptonmm(), cuts.Forward()), "", "m_{ll} [GeV]")
    mll_ee_forward_lowmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalLowMassee(), cuts.Forward()))
    mll_mm_forward_lowmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalLowMassmm(), cuts.Forward()))
    mll_ee_forward_highmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalHighMassee(), cuts.Forward()))
    mll_mm_forward_highmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalHighMassmm(), cuts.Forward()))
    mll_ee_forward_Zmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalZMassee(), cuts.Forward()))
    mll_mm_forward_Zmass = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.Add(cuts.SignalZMassmm(), cuts.Forward()))
    mll_ee_forward_DYmeas = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.AddList([cuts.DYControlNoMassLeptonee(), cuts.Forward(), cuts.DYMass()]))
    mll_mm_forward_DYmeas = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, cuts.AddList([cuts.DYControlNoMassLeptonmm(), cuts.Forward(), cuts.DYMass()]))


    rmue_mll_central = make_rmue(mll_mm_central, mll_ee_central)
    rmue_mll_central.GetYaxis().SetRangeUser(0, 2)
    plot_rmue_mll_central = Canvas("plot_rmue_mll_central", "png", 0.6, 0.6, 0.8, 0.8)
    plot_rmue_mll_central.addHisto(rmue_mll_central, "E1,SAME", "DY", "L", r.kBlack, 1, 0)
    plot_rmue_mll_central.save(0, 0, 0, 4.0)
    
    rmue_met_central = make_rmue(met_mm_central, met_ee_central)
    rmue_met_central.GetYaxis().SetRangeUser(0, 2)
    plot_rmue_met_central = Canvas("plot_rmue_met_central", "png", 0.6, 0.6, 0.8, 0.8)
    plot_rmue_met_central.addHisto(rmue_met_central, "E1,SAME", "DY", "L", r.kBlack, 1, 0)
    plot_rmue_met_central.save(0, 0, 0, 4.0)
    
    rmue_mll_forward = make_rmue(mll_mm_forward, mll_ee_forward)
    rmue_mll_forward.GetYaxis().SetRangeUser(0, 2)
    plot_rmue_mll_forward = Canvas("plot_rmue_mll_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plot_rmue_mll_forward.addHisto(rmue_mll_forward, "E1,SAME", "DY", "L", r.kBlack, 1, 0)
    plot_rmue_mll_forward.save(0, 0, 0, 4.0)
    
    rmue_met_forward = make_rmue(met_mm_forward, met_ee_forward)
    rmue_met_forward.GetYaxis().SetRangeUser(0, 2)
    plot_rmue_met_forward = Canvas("plot_rmue_met_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plot_rmue_met_forward.addHisto(rmue_met_forward, "E1,SAME", "DY", "L", r.kBlack, 1, 0)
    plot_rmue_met_forward.save(0, 0, 0, 4.0)
    
    [rmue_central_signal, rmue_central_l, rmue_central_Z, rmue_central_h] = make_rmue_signal(mll_ee_central_lowmass, mll_mm_central_lowmass, mll_ee_central_Zmass, mll_mm_central_Zmass, mll_ee_central_highmass, mll_mm_central_highmass)
    rmue_central_signal.GetYaxis().SetRangeUser(0, 2)
    rmue_central_signal.SetMarkerSize(0)
    [rmue_forward_signal, rmue_forward_l, rmue_forward_Z, rmue_forward_h] = make_rmue_signal(mll_ee_forward_lowmass, mll_mm_forward_lowmass, mll_ee_forward_Zmass, mll_mm_forward_Zmass, mll_ee_forward_highmass, mll_mm_forward_highmass)
    rmue_forward_signal.GetYaxis().SetRangeUser(0, 2)
    rmue_forward_signal.SetMarkerSize(0)
    
    [rmue_central_meas, rmue_central_met, rmue_central_DY_meas] = make_rmue_meas(mll_ee_central_DYmeas, mll_mm_central_DYmeas, 0.1)
    rmue_central_meas.GetYaxis().SetRangeUser(0, 2)
    rmue_central_meas.GetXaxis().SetRangeUser(20, 300)
    rmue_central_meas.SetMarkerSize(0)
    rmue_central_meas.SetFillColor(r.kBlue-9)
    rmue_central_met.GetYaxis().SetRangeUser(0, 2)
    rmue_central_met.GetXaxis().SetRangeUser(0, 100)
    rmue_central_met.SetMarkerSize(0)
    rmue_central_met.SetFillColor(r.kBlue-9)
    [rmue_forward_meas, rmue_forward_met, rmue_forward_DY_meas] = make_rmue_meas(mll_ee_forward_DYmeas, mll_mm_forward_DYmeas, 0.2)
    rmue_forward_meas.GetYaxis().SetRangeUser(0, 2)
    rmue_forward_meas.GetXaxis().SetRangeUser(20, 300)
    rmue_forward_meas.SetMarkerSize(0)
    rmue_forward_meas.SetFillColor(r.kBlue-9)
    rmue_forward_met.GetYaxis().SetRangeUser(0, 2)
    rmue_forward_met.GetXaxis().SetRangeUser(0, 100)
    rmue_forward_met.SetMarkerSize(0)
    rmue_forward_met.SetFillColor(r.kBlue-9)

    
    rmue_mll_central.SetMarkerSize(1.2)
    rmue_mll_forward.SetMarkerSize(1.2)

    finalplot_rmue_mll_central = Canvas("finalplot_rmue_mll_central", "png", 0.6, 0.65, 0.8, 0.85)
    finalplot_rmue_mll_central.addGraph(rmue_central_meas, "AP2", "<DY Region>", "L", r.kBlue-9, 1, 0)
    finalplot_rmue_mll_central.addLine(20, rmue_central_DY_meas[0], 300, rmue_central_DY_meas[0], r.kBlue-4)
    finalplot_rmue_mll_central.addGraph(rmue_central_signal, "P", "Signal Region", "L", r.kRed, 1, 2)
    finalplot_rmue_mll_central.addHisto(rmue_mll_central, "E1,SAME", "DY Region", "L", r.kBlack, 1, 1)
    finalplot_rmue_mll_central.save(1, 0, 0, 4.0)


    finalplot_rmue_mll_forward = Canvas("finalplot_rmue_mll_forward", "png", 0.6, 0.15, 0.8, 0.35)
    finalplot_rmue_mll_forward.addGraph(rmue_forward_meas, "AP2", "<DY Region>", "L", r.kBlue-9, 1, 0)
    finalplot_rmue_mll_forward.addLine(20, rmue_forward_DY_meas[0], 300, rmue_forward_DY_meas[0], r.kBlue-4)
    finalplot_rmue_mll_forward.addGraph(rmue_forward_signal, "P", "Signal Region", "L", r.kRed, 1, 2)
    finalplot_rmue_mll_forward.addHisto(rmue_mll_forward, "E1,SAME", "DY Region", "L", r.kBlack, 1, 1)
    finalplot_rmue_mll_forward.save(1, 0, 0, 4.0)


    finalplot_rmue_met_central = Canvas("finalplot_rmue_met_central", "png", 0.6, 0.65, 0.8, 0.85)
    finalplot_rmue_met_central.addGraph(rmue_central_met, "AP2", "<DY Region>", "L", r.kBlue-9, 1, 0)
    finalplot_rmue_met_central.addLine(0, rmue_central_DY_meas[0], 100, rmue_central_DY_meas[0], r.kBlue-4)
    finalplot_rmue_met_central.addHisto(rmue_met_central, "E1,SAME", "DY Region", "L", r.kBlack, 1, 1)
    finalplot_rmue_met_central.save(1, 0, 0, 4.0)


    finalplot_rmue_met_forward = Canvas("finalplot_rmue_met_forward", "png", 0.6, 0.15, 0.8, 0.35)
    finalplot_rmue_met_forward.addGraph(rmue_forward_met, "AP2", "<DY Region>", "L", r.kBlue-9, 1, 0)
    finalplot_rmue_met_forward.addLine(0, rmue_forward_DY_meas[0], 100, rmue_forward_DY_meas[0], r.kBlue-4)
    finalplot_rmue_met_forward.addHisto(rmue_met_forward, "E1,SAME", "DY Region", "L", r.kBlack, 1, 1)
    finalplot_rmue_met_forward.save(1, 0, 0, 4.0)


 
    #####Printing results
    a = rounder.Rounder()
    print "Measurement in DY region 60 GeV < mll < 120 GeV central: " + a.toStringB(rmue_central_DY_meas[0], rmue_central_DY_meas[1]) + " +/- " + a.toString(rmue_central_DY_meas[0]*0.1)
    print "Measurement in Signal region 50 GeV < mll < 70 GeV central: " + a.toStringB(rmue_central_l[0], rmue_central_l[1])
    print "Measurement in Signal region 81 GeV < mll < 101 GeV central: " + a.toStringB(rmue_central_Z[0], rmue_central_Z[1])
    print "Measurement in Signal region mll > 120 GeV central: " + a.toStringB(rmue_central_h[0], rmue_central_h[1])
    print "Measurement in DY region 60 GeV < mll < 120 GeV forward: " + a.toStringB(rmue_forward_DY_meas[0], rmue_forward_DY_meas[1]) + " +/- " + a.toString(rmue_forward_DY_meas[0]*0.2)
    print "Measurement in Signal region 50 GeV < mll < 70 GeV forward: " + a.toStringB(rmue_forward_l[0], rmue_forward_l[1])
    print "Measurement in Signal region 81 GeV < mll < 101 GeV forward: " + a.toStringB(rmue_forward_Z[0], rmue_forward_Z[1])
    print "Measurement in Signal region mll > 120 GeV forward: " + a.toStringB(rmue_forward_h[0], rmue_forward_h[1])




