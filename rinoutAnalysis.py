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


def r_inout(lowmass, highmass, zmass)


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
    mll_SF_central = tree.getTH1F(4, "mll_SF_central", "t.lepsMll_Edge", bins, 1, 1, cuts.Add(cuts.DYControlNoMassLeptonSF(), cuts.Central()), "", "m_{ll} [GeV]")
    mll_OF_central = tree.getTH1F(4, "mll_OF_central", "t.lepsMll_Edge", bins, 1, 1, cuts.Add(cuts.DYControlNoMassLeptonOF(), cuts.Central()), "", "m_{ll} [GeV]")
    mll_SF_forward = tree.getTH1F(4, "mll_SF_forward", "t.lepsMll_Edge", bins, 1, 1, cuts.Add(cuts.DYControlNoMassLeptonSF(), cuts.Forward()), "", "m_{ll} [GeV]")
    mll_OF_forward = tree.getTH1F(4, "mll_OF_forward", "t.lepsMll_Edge", bins, 1, 1, cuts.Add(cuts.DYControlNoMassLeptonOF(), cuts.Forward()), "", "m_{ll} [GeV]")

    plot_mll_central = Canvas("plot_mll_central", "png", 0.6, 0.6, 0.8, 0.8)
    plot_mll_central.addHisto(mll_SF_central, "HIST", "SF", "L", r.kRed, 1, 0)
    plot_mll_central.addHisto(mll_OF_central, "HIST,SAME", "OF", "L", r.kBlue, 1, 0)
    plot_mll_central.save(1, 0, 1, 4.0)
    
    plot_mll_forward = Canvas("plot_mll_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plot_mll_forward.addHisto(mll_SF_forward, "HIST", "SF", "L", r.kRed, 1, 0)
    plot_mll_forward.addHisto(mll_OF_forward, "HIST,SAME", "OF", "L", r.kBlue, 1, 0)
    plot_mll_forward.save(1, 0, 1, 4.0)


    #Get yields 
    rinZ_SF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlZSF(), cuts.Central()))
    rinZ_OF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlZOF(), cuts.Central()))
    rinlow_SF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlLowSF(), cuts.Central()))
    rinlow_OF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlLowOF(), cuts.Central()))
    rinhigh_SF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlHighSF(), cuts.Central()))
    rinhigh_OF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlHighOF(), cuts.Central()))

    rinout_low _central = 


 
    rinZ_SF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlZSF(), cuts.Forward()))
    rinZ_OF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlZOF(), cuts.Forward()))
    rinlow_SF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlLowSF(), cuts.Forward()))
    rinlow_OF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlLowOF(), cuts.Forward()))
    rinhigh_SF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlHighSF(), cuts.Forward()))
    rinhigh_OF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000,  cuts.Add(cuts.DYControlHighOF(), cuts.Forward()))
    
    met_rinZ_SF_central = tree.getTH1F(4, "met_rinZ_SF_central", "met_pt", 10, 0, 100,  cuts.Add(cuts.DYControlNoMetZSF(), cuts.Central()), "", "MET [GeV]")
    met_rinZ_OF_central = tree.getTH1F(4, "met_rinZ_OF_central", "met_pt", 10, 0, 100,  cuts.Add(cuts.DYControlNoMetZOF(), cuts.Central()), "", "MET [GeV]")
    met_rinlow_SF_central = tree.getTH1F(4, "met_rinlow_SF_central", "met_pt", 10, 0, 100,  cuts.Add(cuts.DYControlNoMetLowSF(), cuts.Central()), "", "MET [GeV]")
    met_rinlow_OF_central = tree.getTH1F(4, "met_rinlow_OF_central", "met_pt", 10, 0,100,  cuts.Add(cuts.DYControlNoMetLowOF(), cuts.Central()), "", "MET [GeV]")
    met_rinhigh_SF_central = tree.getTH1F(4, "met_rinhigh_SF_central", "met_pt", 10, 0, 100, cuts.Add(cuts.DYControlNoMetHighSF(), cuts.Central()), "", "MET [GeV]")
    met_rinhigh_OF_central = tree.getTH1F(4, "met_rinhigh_OF_central", "met_pt", 10, 0, 100, cuts.Add(cuts.DYControlNoMetHighOF(), cuts.Central()), "", "MET [GeV]")

    met_rinZ_central = met_rinZ_SF_central.Clone("met_rinZ_central")
    met_rinZ_central.Add(met_rinZ_OF_central, -1.0)

    met_rinlow_central = met_rinlow_SF_central.Clone("met_rinlow_central")
    met_rinlow_central.Add(met_rinlow_OF_central, -1.0)

    met_rinhigh_central = met_rinhigh_SF_central.Clone("met_rinhigh_central")
    met_rinhigh_central.Add(met_rinhigh_OF_central, -1.0)

    met_rinhigh_central.Divide(met_rinZ_central)
    met_rinlow_central.Divide(met_rinZ_central)

    plot_metrinhigh_central = Canvas("plot_metrinhigh_central", "png",0.6, 0.6, 0.8, 0.8)
    plot_metrinhigh_central.addHisto(met_rinhigh_central, "E1", "data", "P", r.kRed, 1, 0)
    plot_metrinhigh_central.save(1, 0, 0, 4.0)



