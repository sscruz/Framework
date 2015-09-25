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
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors
import math, sys, optparse, array
import Rounder as rounder


import include.helper     as helper
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample

def calc_rmue(Nmm, Nee, Emm, Eee):

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
            if(Nee == 0 or Nmm == 0):
              continue
            val = calc_rmue(Nmm, Nee, Emm, Eee)
            ratio.SetBinContent(i, val[0])
            ratio.SetBinError(i, val[1])
    return ratio


if __name__ == "__main__":


    parser = optparse.OptionParser(usage="usage: %prog [options] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-m", "--mode", action="store", dest="mode", default="rmue", help="Operation mode")
    (options, args) = parser.parse_args()

    inputFileName = args[0]
    doEta = args[1]

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTJets', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    #tree = treeMC
    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ####Cuts needed by rmue
    cuts = CutManager.CutManager()

    lumi = 0.020


    #for eta in ['central']:#, 'forward']:
    #for eta in ['forward']:#, 'forward']:
    for eta in [doEta]:#, 'forward']:
        regions = []
        dy_nomass = helper.rmueRegion('DY_nomass', eta, 
                           [cuts.DYControlRegion],
                           [20, 30, 40, 50, 60, 70, 81, 101, 120, 150, 180, 220, 260, 300],
                           False, 'mll', True)
        regions.append(dy_nomass)
        dy_onZ    = helper.rmueRegion('DY_onZ'   , eta, 
                           [cuts.DYControlRegion, cuts.DYmass],
                           [60, 120],
                           True, 'mll', True)
        regions.append(dy_onZ)
        dy_nomet  = helper.rmueRegion('DY_nomet' , eta, 
                           [cuts.nj2, cuts.DYmass],
                           range(0,110,10),
                           True, 'met', True)
        regions.append(dy_nomet)
        sig_lm    = helper.rmueRegion('Signal_lowmass' , eta, 
                           [cuts.METJetsSignalRegion, cuts.lowmass],
                           [20,  70],
                           True, 'mll', False)
        regions.append(sig_lm)
        sig_onZ   = helper.rmueRegion('Signal_onZ' , eta, 
                           [cuts.METJetsSignalRegion, cuts.Zmass],
                           [81, 101],
                           True, 'mll', False)
        regions.append(sig_onZ)
        sig_hm    = helper.rmueRegion('Signal_highmass' , eta, 
                           [cuts.METJetsSignalRegion, cuts.highmass],
                           [120, 300],
                           True, 'mll', False)
        regions.append(sig_hm)


        for region in regions:
            print 'i am at region', region.name
            cuts_ee = cuts.AddList([cuts.GoodLeptonee()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+region.cuts)
            cuts_mm = cuts.AddList([cuts.GoodLeptonmm()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+region.cuts)


            for tree in ([treeMC, treeDA] if region.doData else [treeMC]):

                if region.doMll:
                    region.mll_ee = tree.getTH1F(lumi, "mll_ee_"+region.cenFwd, "t.lepsMll_Edge", region.bins, 1, 1, cuts_ee, "", "m_{ll} (GeV)")
                    region.mll_mm = tree.getTH1F(lumi, "mll_mm_"+region.cenFwd, "t.lepsMll_Edge", region.bins, 1, 1, cuts_mm, "", "m_{ll} (GeV)")

                    region.set_rmue_mll(make_rmue(region.mll_mm, region.mll_ee), False if tree == treeMC else True)

                if region.doMet:
                    region.met_ee = tree.getTH1F(lumi, "met_ee_"+region.cenFwd, "met_pt", region.bins, 1, 1, cuts_ee, "", "ME_{T} (GeV)")
                    region.met_mm = tree.getTH1F(lumi, "met_mm_"+region.cenFwd, "met_pt", region.bins, 1, 1, cuts_mm, "", "ME_{T} (GeV)")

                    region.set_rmue_met(make_rmue(region.met_mm, region.met_ee), False if tree == treeMC else True)
    


        ## =================
        ## MAKE THE Mll PLOT
        ## =================
        meas_rmue_mc   = dy_onZ.rmue_mll.GetBinContent(1)
        meas_rmue_mc_e = math.sqrt(dy_onZ.rmue_mll.GetBinError(1)**2 + 0.1**2)
        meas_rmue_da   = dy_onZ.rmue_mll_data.GetBinContent(1)
        meas_rmue_da_e = math.sqrt(dy_onZ.rmue_mll_data.GetBinError(1)**2 + 0.1**2)

        plot_rmue_mll = Canvas.Canvas("rmue/plot_rmue_mll_"+dy_nomass.cenFwd, "png,pdf", 0.6, 0.15, 0.8, 0.35)
        plot_rmue_mll.addHisto(dy_nomass.rmue_mll     , "E,SAME", "DY"       , "PL", r.kRed+1 , 1, 0)
        plot_rmue_mll.addHisto(dy_nomass.rmue_mll_data, "E,SAME", "DY - data", "PL", r.kBlack , 1, 5)
        plot_rmue_mll.addGraph(sig_lm .rmue_mll_gr, "PZ", "<SR low>"  , "PL", r.kBlue-7, 1, 1)
        plot_rmue_mll.addGraph(sig_onZ.rmue_mll_gr, "PZ", "<SR onZ>"  , "PL", r.kBlue-8, 1, 2)
        plot_rmue_mll.addGraph(sig_hm .rmue_mll_gr, "PZ", "<SR high>" , "PL", r.kBlue-9, 1, 3)
        plot_rmue_mll.addGraph(dy_onZ .rmue_mll_gr, "PZ", "<DY onZ>"  , "PL", r.kCyan+1, 1, 4)
        plot_rmue_mll.addLine (dy_nomass.rmue_mll.GetXaxis().GetXmin(), 1.                   , dy_nomass.rmue_mll.GetXaxis().GetXmax(), 1.                   , r.kGreen)
        plot_rmue_mll.addBand (dy_nomass.rmue_mll.GetXaxis().GetXmin(), meas_rmue_da-meas_rmue_da_e, dy_nomass.rmue_mll.GetXaxis().GetXmax(), meas_rmue_da+meas_rmue_da_e, r.kGray+1, 0.2)
        plot_rmue_mll.addLine (dy_nomass.rmue_mll.GetXaxis().GetXmin(), meas_rmue_da               , dy_nomass.rmue_mll.GetXaxis().GetXmax(), meas_rmue_da               , r.kGray+1)
        plot_rmue_mll.addLatex(0.2, 0.2, dy_nomass.cenFwd)
        plot_rmue_mll.save(1, 0, 0, lumi)

        ## =================
        ## MAKE THE MET PLOT
        ## =================
        plot_rmue_met = Canvas.Canvas("rmue/plot_rmue_met_"+region.cenFwd, "png,pdf", 0.6, 0.15, 0.8, 0.35)
        plot_rmue_met.addHisto(dy_nomet.rmue_met     , "E,SAME", "DY"       , "PL", r.kRed+1 , 1, 0)
        plot_rmue_met.addHisto(dy_nomet.rmue_met_data, "E,SAME", "DY - data", "PL", r.kBlack , 1, 1)
        plot_rmue_met.addLine (dy_nomet.rmue_met.GetXaxis().GetXmin(), 1.           , dy_nomet.rmue_met.GetXaxis().GetXmax(), 1.           , r.kGreen)
        plot_rmue_met.addBand (dy_nomet.rmue_met.GetXaxis().GetXmin(), meas_rmue_da-meas_rmue_da_e, dy_nomet.rmue_met.GetXaxis().GetXmax(), meas_rmue_da+meas_rmue_da_e, r.kGray+1, 0.2)
        plot_rmue_met.addLine (dy_nomet.rmue_met.GetXaxis().GetXmin(), meas_rmue_da               , dy_nomet.rmue_met.GetXaxis().GetXmax(), meas_rmue_da               , r.kGray+1)
        plot_rmue_met.addLatex(0.2, 0.2, dy_nomass.cenFwd)
        plot_rmue_met.save(1, 0, 0, lumi)

        dy_onZ    .printValues()
        dy_nomass .printValues()
        sig_onZ   .printValues()
        sig_lm    .printValues()
