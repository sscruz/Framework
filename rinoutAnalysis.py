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


def r_inout(massSF, massOF, refmassSF, refmassOF):

    themass = massSF[0]-massOF[0]
    themassE = math.sqrt(massSF[1]*massSF[1] + massOF[1]*massOF[1])
    therefmass = refmassSF[0]-refmassOF[0]
    therefmassE = math.sqrt(refmassSF[1]*refmassSF[1] + refmassOF[1]*refmassOF[1])

    sol = [0, 0]
    sol[0] = themass/therefmass
    sol[1] = math.sqrt((themassE/therefmass)*(themassE/therefmass) + (themass*therefmassE/(therefmass*therefmass))*(themass*therefmassE/(therefmass*therefmass)))

    return sol

def prepareDefault(x, y, ex, ey, label):
 
    xarr = array("d", [x])
    xarr_e = array("d", [ex])
    yarr = array("d", [y])
    yarr_e = array("d", [ey])

    graph = TGraphErrors(1, xarr, yarr, xarr_e, yarr_e)
    graph.GetYaxis().SetTitle("r_{#mu e}")
    graph.GetXaxis().SetTitle(label)
    graph.SetFillColor(r.kBlue-9)
    graph.SetMarkerSize(0)
    graph.GetXaxis().SetRangeUser(x-ex, x+ex)
    graph.GetYaxis().SetRangeUser(0, 2*y)

    return graph


if __name__ == "__main__":

    parser = OptionParser(usage="usage: %prog [options] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-m", "--mode", action="store", dest="mode", default="rmue", help="Operation mode")
    (options, args) = parser.parse_args()

    if len(args) != 2:
      parser.error("wrong number of arguments")

    if args[1] != "MC" and args[1] != "DATA":
      parser.error("Second argument must be MC or DATA")

    inputFileName = args[0]
    typeOfSample = args[1]
    isData = 1 if typeOfSample == 'DATA' else 0

    tree = Tree(inputFileName, typeOfSample, isData)
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 

    ####Cuts needed by rinout
    cuts = CutManager()

    regionZ_CentralSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYControlRegion, cuts.Zmass, cuts.Central()])
    regionZ_ForwardSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYControlRegion, cuts.Zmass, cuts.Forward()])

    regionHigh_CentralSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYControlRegion, cuts.highmass, cuts.Central()])
    regionHigh_ForwardSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYControlRegion, cuts.highmass, cuts.Forward()])

    regionLow_CentralSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYControlRegion, cuts.lowmass, cuts.Central()])
    regionLow_ForwardSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYControlRegion, cuts.lowmass, cuts.Forward()])

    region_Central_nomassSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYControlRegion, cuts.Central()])
    region_Forward_nomassSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYControlRegion, cuts.Forward()])

    regionZ_Central_nometSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.nj2, cuts.Zmass, cuts.Central()])
    regionZ_Forward_nometSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.nj2, cuts.Zmass, cuts.Forward()])

    regionHigh_Central_nometSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.nj2, cuts.highmass, cuts.Central()])
    regionHigh_Forward_nometSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.nj2, cuts.highmass, cuts.Forward()])

    regionLow_Central_nometSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.nj2, cuts.lowmass, cuts.Central()])
    regionLow_Forward_nometSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.nj2, cuts.lowmass, cuts.Forward()])

    regionZ_Central_nojetsSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYmet, cuts.Zmass, cuts.Central()])
    regionZ_Forward_nojetsSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYmet, cuts.Zmass, cuts.Forward()])

    regionHigh_Central_nojetsSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYmet, cuts.highmass, cuts.Central()])
    regionHigh_Forward_nojetsSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYmet, cuts.highmass, cuts.Forward()])

    regionLow_Central_nojetsSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYmet, cuts.lowmass, cuts.Central()])
    regionLow_Forward_nojetsSF = cuts.AddList([cuts.GoodLeptonSF(), cuts.DYmet, cuts.lowmass, cuts.Forward()])

    regionZ_CentralOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYControlRegion, cuts.Zmass, cuts.Central()])
    regionZ_ForwardOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYControlRegion, cuts.Zmass, cuts.Forward()])

    regionHigh_CentralOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYControlRegion, cuts.highmass, cuts.Central()])
    regionHigh_ForwardOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYControlRegion, cuts.highmass, cuts.Forward()])

    regionLow_CentralOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYControlRegion, cuts.lowmass, cuts.Central()])
    regionLow_ForwardOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYControlRegion, cuts.lowmass, cuts.Forward()])

    region_Central_nomassOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYControlRegion, cuts.Central()])
    region_Forward_nomassOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYControlRegion, cuts.Forward()])

    regionZ_Central_nometOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.nj2, cuts.Zmass, cuts.Central()])
    regionZ_Forward_nometOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.nj2, cuts.Zmass, cuts.Forward()])

    regionHigh_Central_nometOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.nj2, cuts.highmass, cuts.Central()])
    regionHigh_Forward_nometOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.nj2, cuts.highmass, cuts.Forward()])

    regionLow_Central_nometOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.nj2, cuts.lowmass, cuts.Central()])
    regionLow_Forward_nometOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.nj2, cuts.lowmass, cuts.Forward()])

    regionZ_Central_nojetsOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYmet, cuts.Zmass, cuts.Central()])
    regionZ_Forward_nojetsOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYmet, cuts.Zmass, cuts.Forward()])

    regionHigh_Central_nojetsOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYmet, cuts.highmass, cuts.Central()])
    regionHigh_Forward_nojetsOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYmet, cuts.highmass, cuts.Forward()])

    regionLow_Central_nojetsOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYmet, cuts.lowmass, cuts.Central()])
    regionLow_Forward_nojetsOF = cuts.AddList([cuts.GoodLeptonOF(), cuts.DYmet, cuts.lowmass, cuts.Forward()])


    ############# Mass plots to show the invariant mass distribution for SF and OF events #################################
    mll_SF_central = tree.getTH1F(4, "mll_SF_central", "t.lepsMll_Edge", 28, 20, 300, region_Central_nomassSF, "", "m_{ll} [GeV]")
    mll_OF_central = tree.getTH1F(4, "mll_OF_central", "t.lepsMll_Edge", 28, 20, 300, region_Central_nomassOF, "", "m_{ll} [GeV]")
    mll_SF_forward = tree.getTH1F(4, "mll_SF_forward", "t.lepsMll_Edge", 28, 20, 300, region_Forward_nomassSF, "", "m_{ll} [GeV]")
    mll_OF_forward = tree.getTH1F(4, "mll_OF_forward", "t.lepsMll_Edge", 28, 20, 300, region_Forward_nomassOF, "", "m_{ll} [GeV]")

    plot_mll_central = Canvas("plot_mll_central", "png", 0.6, 0.6, 0.8, 0.8)
    plot_mll_central.addHisto(mll_SF_central, "HIST", "SF", "L", r.kRed, 1, 0)
    plot_mll_central.addHisto(mll_OF_central, "HIST,SAME", "OF", "L", r.kBlue, 1, 0)
    plot_mll_central.save(1, 0, 1, 4.0)
    
    plot_mll_forward = Canvas("plot_mll_forward", "png", 0.6, 0.6, 0.8, 0.8)
    plot_mll_forward.addHisto(mll_SF_forward, "HIST", "SF", "L", r.kRed, 1, 0)
    plot_mll_forward.addHisto(mll_OF_forward, "HIST,SAME", "OF", "L", r.kBlue, 1, 0)
    plot_mll_forward.save(1, 0, 1, 4.0)



    ########### Get yields ############################################################################################### 
    rinZ_SF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionZ_CentralSF)
    rinZ_OF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionZ_CentralOF)
    rinlow_SF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionLow_CentralSF)
    rinlow_OF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionLow_CentralOF)
    rinhigh_SF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionHigh_CentralSF)
    rinhigh_OF_central = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionHigh_CentralOF)

    rinZ_SF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionZ_ForwardSF)
    rinZ_OF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionZ_ForwardOF)
    rinlow_SF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionLow_ForwardSF)
    rinlow_OF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionLow_ForwardOF)
    rinhigh_SF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionHigh_ForwardSF)
    rinhigh_OF_forward = tree.getYields(4, "t.lepsMll_Edge", 20, 1000, regionHigh_ForwardOF)

    rinout_low_central = r_inout(rinlow_SF_central, rinlow_OF_central, rinZ_SF_central, rinZ_OF_central)
    rinout_high_central = r_inout(rinhigh_SF_central, rinhigh_OF_central, rinZ_SF_central, rinZ_OF_central)
    rinout_low_forward = r_inout(rinlow_SF_forward, rinlow_OF_forward, rinZ_SF_forward, rinZ_OF_forward)
    rinout_high_forward = r_inout(rinhigh_SF_forward, rinhigh_OF_forward, rinZ_SF_forward, rinZ_OF_forward)

    a = rounder.Rounder()
    print "r_in/out central, low: " + a.toStringB(rinout_low_central[0], rinout_low_central[1]) + " +/- " + a.toString(rinout_low_central[0]*0.2)
    print "r_in/out central, high: " + a.toStringB(rinout_high_central[0], rinout_high_central[1]) + " +/- " + a.toString(rinout_high_central[0]*0.2)
    print "r_in/out forward, low: " + a.toStringB(rinout_low_forward[0], rinout_low_forward[1]) + " +/- " + a.toString(rinout_low_forward[0]*0.2)
    print "r_in/out forward, high: " + a.toStringB(rinout_high_forward[0], rinout_high_forward[1]) +  " +/- " + a.toString(rinout_high_forward[0]*0.2)


    

    ########### Default histograms ############################################################################################
    met_default_low_central = prepareDefault(50, rinout_low_central[0], 50, rinout_low_central[1], "MET [GeV]")
    met_default_high_central = prepareDefault(50, rinout_high_central[0], 50, rinout_high_central[1], "MET [GeV]")
    met_default_low_forward = prepareDefault(50, rinout_low_forward[0], 50, rinout_low_forward[1], "MET [GeV]")
    met_default_high_forward = prepareDefault(50, rinout_high_forward[0], 50, rinout_high_forward[1], "MET [GeV]")
    
    jets_default_low_central = prepareDefault(5, rinout_low_central[0], 5, rinout_low_central[1], "MET [GeV]")
    jets_default_high_central = prepareDefault(5, rinout_high_central[0], 5, rinout_high_central[1], "MET [GeV]")
    jets_default_low_forward = prepareDefault(5, rinout_low_forward[0], 5, rinout_low_forward[1], "MET [GeV]")
    jets_default_high_forward = prepareDefault(5, rinout_high_forward[0], 5, rinout_high_forward[1], "MET [GeV]")



    ########### Get Systematics ############################################################################################### 
    met_rinZ_SF_central = tree.getTH1F(4, "met_rinZ_SF_central", "met_pt", 10, 0, 100, regionZ_Central_nometSF, "", "MET [GeV]")
    met_rinZ_OF_central = tree.getTH1F(4, "met_rinZ_OF_central", "met_pt", 10, 0, 100, regionZ_Central_nometOF, "", "MET [GeV]")
    met_rinlow_SF_central = tree.getTH1F(4, "met_rinlow_SF_central", "met_pt", 10, 0, 100, regionLow_Central_nometSF, "", "MET [GeV]")
    met_rinlow_OF_central = tree.getTH1F(4, "met_rinlow_OF_central", "met_pt", 10, 0, 100, regionLow_Central_nometOF, "", "MET [GeV]")
    met_rinhigh_SF_central = tree.getTH1F(4, "met_rinhigh_SF_central", "met_pt", 10, 0, 100, regionHigh_Central_nometSF , "", "MET [GeV]")
    met_rinhigh_OF_central = tree.getTH1F(4, "met_rinhigh_OF_central", "met_pt", 10, 0, 100, regionHigh_Central_nometOF , "", "MET [GeV]")

    met_rinZ_central = met_rinZ_SF_central.Clone("met_rinZ_central")
    met_rinZ_central.Add(met_rinZ_OF_central, -1.0)
    met_rinlow_central = met_rinlow_SF_central.Clone("met_rinlow_central")
    met_rinlow_central.Add(met_rinlow_OF_central, -1.0)
    met_rinhigh_central = met_rinhigh_SF_central.Clone("met_rinhigh_central")
    met_rinhigh_central.Add(met_rinhigh_OF_central, -1.0)
    met_rinhigh_central.Divide(met_rinZ_central)
    met_rinlow_central.Divide(met_rinZ_central)

    plot_metrinhigh_central = Canvas("plot_metrinhigh_central", "png",0.6, 0.6, 0.8, 0.8)
    plot_metrinhigh_central.addGraph(met_default_high_central, "AP2", "Measured", "F", r.kBlue-9, 1, 0)
    plot_metrinhigh_central.addLine(0, rinout_high_central[0], 100, rinout_high_central[0], r.kBlue-4)
    plot_metrinhigh_central.addHisto(met_rinhigh_central, "E1,SAME", "data", "P", r.kRed, 1, 0)
    plot_metrinhigh_central.save(1, 0, 0, 4.0)



 
