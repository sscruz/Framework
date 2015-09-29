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
import include.Region     as Region
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

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTJets', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C',
                  'DoubleMuon_Run2015D', 'DoubleEG_Run2015D', 'MuonEG_Run2015D']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    #tree = treeMC
    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ####Cuts needed by rmue
    cuts = CutManager.CutManager()

    lumi = 0.150


    regions = []
    dy_nomass = Region.region('DY_nomass',
                       [cuts.DYControlRegion],
                       ['mll'],
                       [[20, 30, 40, 50, 60, 70, 81, 101, 120, 150, 180, 220, 260, 300]],
                       True)
    regions.append(dy_nomass)
    dy_onZ    = Region.region('DY_onZ',
                       [cuts.DYControlRegion, cuts.DYmass],
                       ['mll'],
                       [[60, 120]],
                       True)
    regions.append(dy_onZ)
    dy_nomet  = Region.region('DY_nomet',
                       [cuts.nj2, cuts.DYmass],
                       ['met'],
                       [range(0,110,10)],
                       True)
    regions.append(dy_nomet)
    sig_lm    = Region.region('Signal_lowmass',
                       [cuts.METJetsSignalRegion, cuts.lowmass],
                       ['mll'],
                       [[20,  70]],
                       False)
    regions.append(sig_lm)
    sig_onZ   = Region.region('Signal_onZ',
                       [cuts.METJetsSignalRegion, cuts.Zmass],
                       ['mll'],
                       [[81, 101]],
                       False)
    regions.append(sig_onZ)
    sig_hm    = Region.region('Signal_highmass',
                       [cuts.METJetsSignalRegion, cuts.highmass],
                       ['mll'],
                       [[120, 300]],
                       False)
    regions.append(sig_hm)
    
    
    for reg in regions:
        print 'i am at region', reg.name
        for eta in ['central', 'forward']:
            print '... in %s' %(eta)

            cuts_ee = cuts.AddList([cuts.GoodLeptonee()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
            cuts_mm = cuts.AddList([cuts.GoodLeptonmm()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
    
    
            for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):
    
                dataMC = 'DATA' if tree == treeDA else 'MC'

                if 'mll' in reg.rvars:
                    reg.mll_ee = tree.getTH1F(lumi, "mll_ee_"+eta, "t.lepsMll_Edge", reg.bins[reg.rvars.index('mll')], 1, 1, cuts_ee, "", "m_{ll} (GeV)")
                    reg.mll_mm = tree.getTH1F(lumi, "mll_mm_"+eta, "t.lepsMll_Edge", reg.bins[reg.rvars.index('mll')], 1, 1, cuts_mm, "", "m_{ll} (GeV)")
    
                    reg.mll.setHisto(make_rmue(reg.mll_mm, reg.mll_ee), dataMC, eta)
    
                if 'met' in reg.rvars:
                    reg.met_ee = tree.getTH1F(lumi, "met_ee_"+eta, "met_pt", reg.bins[reg.rvars.index('met')], 1, 1, cuts_ee, "", "ME_{T} (GeV)")
                    reg.met_mm = tree.getTH1F(lumi, "met_mm_"+eta, "met_pt", reg.bins[reg.rvars.index('met')], 1, 1, cuts_mm, "", "ME_{T} (GeV)")
    
                    reg.met.setHisto(make_rmue(reg.met_mm, reg.met_ee), dataMC, eta)
    
    
    
    for eta in ['central', 'forward']:
        ## =================
        ## MAKE THE Mll PLOT
        ## =================
        meas_rmue_mc   = dy_onZ.mll.getHisto('MC'  , eta).GetBinContent(1)
        meas_rmue_da   = dy_onZ.mll.getHisto('DATA', eta).GetBinContent(1)
        meas_rmue_mc_e = math.sqrt(dy_onZ.mll.getHisto('MC'  , eta).GetBinError(1)**2 + 0.1**2)
        meas_rmue_da_e = math.sqrt(dy_onZ.mll.getHisto('DATA', eta).GetBinError(1)**2 + 0.1**2)
        upEdge = meas_rmue_da+meas_rmue_da_e
        dnEdge = meas_rmue_da-meas_rmue_da_e
        middle = meas_rmue_da
        
        plot_rmue_mll = Canvas.Canvas("rmue/plot_rmue_mll_"+eta, "png,pdf", 0.6, 0.15, 0.8, 0.35)
        plot_rmue_mll.addHisto(dy_nomass.mll.getHisto('MC'  , eta), "E,SAME", "DY"       , "PL", r.kRed+1 , 1, 0)
        plot_rmue_mll.addHisto(dy_nomass.mll.getHisto('DATA', eta), "E,SAME", "DY - data", "PL", r.kBlack , 1, 5)
        plot_rmue_mll.addGraph(sig_lm   .mll.getGraph('MC'  ,eta), "PZ", "SR low - MC"  , "PL", r.kBlue-7, 1, 1)
        plot_rmue_mll.addGraph(sig_onZ  .mll.getGraph('MC'  ,eta), "PZ", "SR onZ - MC"  , "PL", r.kBlue-8, 1, 2)
        plot_rmue_mll.addGraph(sig_hm   .mll.getGraph('MC'  ,eta), "PZ", "SR high - MC" , "PL", r.kBlue-9, 1, 3)
        plot_rmue_mll.addGraph(dy_onZ   .mll.getGraph('MC'  ,eta), "PZ", "DY onZ - MC"  , "PL", r.kCyan+1, 1, 4)
        plot_rmue_mll.addLine (dy_nomass.mll.getHisto('MC', eta).GetXaxis().GetXmin(),     1., dy_nomass.mll.getHisto('MC', eta).GetXaxis().GetXmax(),     1., r.kGreen)
        plot_rmue_mll.addBand (dy_nomass.mll.getHisto('MC', eta).GetXaxis().GetXmin(), dnEdge, dy_nomass.mll.getHisto('MC', eta).GetXaxis().GetXmax(), upEdge, r.kGray+1, 0.2)
        plot_rmue_mll.addLine (dy_nomass.mll.getHisto('MC', eta).GetXaxis().GetXmin(), middle, dy_nomass.mll.getHisto('MC', eta).GetXaxis().GetXmax(), middle, r.kGray+1)
        plot_rmue_mll.addLatex(0.2, 0.2, eta)
        plot_rmue_mll.save(1, 0, 0, lumi)
        
        ## =================
        ## MAKE THE MET PLOT
        ## =================
        plot_rmue_met = Canvas.Canvas("rmue/plot_rmue_met_"+eta, "png,pdf", 0.6, 0.15, 0.8, 0.35)
        plot_rmue_met.addHisto(dy_nomet.met.getHisto('MC'  , eta), "E,SAME", "DY"       , "PL", r.kRed+1 , 1, 0)
        plot_rmue_met.addHisto(dy_nomet.met.getHisto('DATA', eta), "E,SAME", "DY - data", "PL", r.kBlack , 1, 1)
        plot_rmue_met.addLine (dy_nomet.met.getHisto('MC'  , eta).GetXaxis().GetXmin(),     1., dy_nomet.met.getHisto('MC', eta).GetXaxis().GetXmax(),     1., r.kGreen)
        plot_rmue_met.addBand (dy_nomet.met.getHisto('MC'  , eta).GetXaxis().GetXmin(), dnEdge, dy_nomet.met.getHisto('MC', eta).GetXaxis().GetXmax(), upEdge, r.kGray+1, 0.2)
        plot_rmue_met.addLine (dy_nomet.met.getHisto('MC'  , eta).GetXaxis().GetXmin(), middle, dy_nomet.met.getHisto('MC', eta).GetXaxis().GetXmax(), middle, r.kGray+1)
        plot_rmue_met.addLatex(0.2, 0.2, eta)
        plot_rmue_met.save(1, 0, 0, lumi)
        
    ## =================
    ## PRINT AND SAVE ==
    ## =================
    dy_onZ    .mll.printValues()
    dy_nomass .mll.printValues()
    sig_onZ   .mll.printValues()
    sig_lm    .mll.printValues()
    
    dy_onZ .mll.saveInFile(['rmue', 'dycr_dym'], 0.1)
    sig_lm .mll.saveInFile(['rmue', 'sr_lm'   ], 0.1)
    sig_onZ.mll.saveInFile(['rmue', 'sr_onZ'  ], 0.1)
    sig_hm .mll.saveInFile(['rmue', 'sr_hm'   ], 0.1)

