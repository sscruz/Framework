#####################################################################
######                                                              #
###### 8=========D         OO                        8==========D   #  
###### OO                  ||                                ,88   #
###### ||                  ||                              ,88'    #  
###### ||O---\     ,adPPYb,|| ,adPPYb,d8  ,adPPYba,      ,88'      #
###### ||O---/    a8'    `Y||a8'    `Y88 a8P_____88    ,88'        #
###### ||         8b       ||8b       88 8PP'''''''  ,88'          #
###### \/         '8a,   ,d||'8a,   ,d88 '8b,   ,aa 88'            #
###### 8=========D `'8bbdP'\/  `'YbbdP'Y8  `'Ybbd8'' 8==========D   #
######                       aa,    ,88                             #
######                         'Y8bbdP'                             #
######                                                              #
#####################################################################

import ROOT as r
import math as math
import sys
import Canvas, CutManager

from optparse import OptionParser
from ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
from Sample import Sample, Block, Tree, selectSamples


def make_rsfof(histo_sf, histo_of, isData):

    ratio = histo_sf.Clone('rsfof_' + histo_sf.GetName())
    ratio.Divide(histo_of)
    ratio.GetYaxis().SetTitle('r_{SFOF}')

    doFit = 0
    if doFit:
        fit = TF1('myfit','pol0', ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
        fit.SetLineColor(r.kBlack if isData else r.kRed+1)
        fit.SetLineStyle(2)
        
        ratio.Fit('myfit','0')

    ratio.GetYaxis().SetRangeUser(0.5,1.5)

    f = open('txts/'+ratio.GetName()+'_values.txt', 'w')
    for i in range(1, ratio.GetNbinsX()+1):
        min, max = ratio.GetBinLowEdge(i), ratio.GetBinLowEdge(i)+ratio.GetBinWidth(i)
        print    '%10s : R_SFOF in [%.2f, %.2f] GeV:\t%.3f +- %.3f'    %('DATA' if isData else 'MC', min, max, ratio.GetBinContent(i), ratio.GetBinError(i) )
        f.write( '%10s : R_SFOF in [%.2f, %.2f] GeV:\t%.3f +- %.3f \n' %('DATA' if isData else 'MC', min, max, ratio.GetBinContent(i), ratio.GetBinError(i) ) )
    f.close()


    return ratio


class rsfofRegion:
    def __init__(self, name, cenFwd, cuts, bins, isGraph, var, doData):
        self.name      = name
        self.cuts      = cuts
        self.bins      = bins
        self.isGraph   = isGraph
        self.cenFwd    = cenFwd
        self.isCentral = (cenFwd == 'central')
        self.doData    = doData

    def set_rsfof(self, rsfof, data=0):
        if not data:
            self.rsfof         = rsfof
            self.rsfof_gr      = r.TGraphErrors(rsfof)
            self.rsfof.GetYaxis().SetRangeUser(0, 2)
        else:
            self.rsfof_data    = rsfof
            self.rsfof_data_gr = r.TGraphErrors(rsfof)
            self.rsfof_data.GetYaxis().SetRangeUser(0, 2)

    def printValues(self):
        print 'REGION', self.name, self.cenFwd
        for _bin in range(1,self.rsfof.GetNbinsX()+1):
            print 'r_sfof in [%.0f, %.0f] in %s: %.3f +- %.3f' %(
                  self.rsfof.GetXaxis().GetBinLowEdge(_bin),
                  self.rsfof.GetXaxis().GetBinUpEdge (_bin),
                  self.cenFwd,
                  self.rsfof.GetBinContent(_bin),
                  self.rsfof.GetBinError(_bin))


if __name__ == '__main__':

    print 'Starting r_SFOF analysis...'
    parser = OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    parser.add_option('-t', '--trigger', action='store', type='int', dest='triggerFlag', default='1', help='Trigger cut. Set to 0 if you want to run without trigger')
    (options, args) = parser.parse_args()


    if len(args) != 2:
      parser.error('wrong number of arguments')

    inputFileName = args[0]
    centralForward = args[1]

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTJets']#, 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C']
    treeMC = Tree(selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Tree(selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    print 'Trees successfully loaded...'


    lumi = 0.020
    print 'Running with an integrated luminosity of', lumi,'fb-1'
   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    for eta in [centralForward]:#, 'forward']:
        regions = []
        ttjets_meas    = rsfofRegion('ttjets_meas', eta, 
                                     [cuts.METJetsControlRegion],
                                     [20, 70, 81, 101, 120, 300],
                                     False, 'mll', True)
        regions.append(ttjets_meas)
        ttjets_sig_lm  = rsfofRegion('ttjets_sig_lm', eta, 
                                     [cuts.METJetsSignalRegion, cuts.lowmass],
                                     [20, 70],
                                     True, 'mll', False)
        regions.append(ttjets_sig_lm)
        ttjets_sig_onZ = rsfofRegion('ttjets_sig_onZ', eta, 
                                     [cuts.METJetsSignalRegion, cuts.Zmass],
                                     [81, 101],
                                     True, 'mll', False)
        regions.append(ttjets_sig_onZ)
        ttjets_sig_hm  = rsfofRegion('ttjets_sig_hm', eta, 
                                     [cuts.METJetsSignalRegion, cuts.highmass],
                                     [120, 300],
                                     True, 'mll', False)
        regions.append(ttjets_sig_hm)
        ## ttjets_inc     = rsfofRegion('ttjets_inc', eta, 
        ##                              [cuts.nj2],
        ##                              [20, 70, 81, 101, 120, 300],
        ##                              False, 'mll', True)
        ## regions.append(ttjets_inc)

        print '...in', eta

        for region in regions:

            print 'i am at region', region.name
            cuts_sf = cuts.AddList([cuts.GoodLeptonSF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+region.cuts)
            cuts_of = cuts.AddList([cuts.GoodLeptonOF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+region.cuts)

            for tree in ([treeMC, treeDA] if region.doData else [treeMC]):

                region.mll_sf = tree.getTH1F(lumi, 'mll_sf_'+region.cenFwd, 't.lepsMll_Edge', region.bins, 1, 1, cuts_sf, '', "m_{ll} (GeV)")
                region.mll_of = tree.getTH1F(lumi, 'mll_of_'+region.cenFwd, 't.lepsMll_Edge', region.bins, 1, 1, cuts_of, '', "m_{ll} (GeV)")

                isData = False if tree == treeMC else True
                region.set_rsfof(make_rsfof(region.mll_sf, region.mll_of, isData), isData)
    

        ## =================
        ## MAKE THE Mll PLOT
        ## =================
        #meas_rsfof = 
        plot_rsfof = Canvas.Canvas('rsfof/plot_rsfof_mll_'+ttjets_meas.cenFwd, 'png,pdf', 0.5, 0.2, 0.75, 0.4)
        plot_rsfof.addHisto(ttjets_meas   .rsfof     , 'PE'     , 'ttjets region - MC'  , 'PL', r.kRed+1 , 1, 0)
        plot_rsfof.addHisto(ttjets_meas   .rsfof_data, 'PE,SAME', 'ttjets region - DATA', 'PL', r.kBlack , 1, 1)
        plot_rsfof.addGraph(ttjets_sig_lm .rsfof     , 'PZ,SAME', 'signal lowmass - MC' , 'PL', r.kBlue-8, 1, 2)
        plot_rsfof.addGraph(ttjets_sig_onZ.rsfof     , 'PZ,SAME', 'signal onZ - MC'     , 'PL', r.kBlue-7, 1, 3)
        plot_rsfof.addGraph(ttjets_sig_hm .rsfof     , 'PZ,SAME', 'signal highmass - MC', 'PL', r.kBlue-9, 1, 4)
        #plot_rsfof.addHisto(ttjets_inc .rsfof     , 'E,SAME', 'incl.  region - MC'  , 'PL', r.kBlack , 1, 1)
        plot_rsfof.addBand (ttjets_meas.rsfof.GetXaxis().GetXmin(), 0.9, ttjets_meas.rsfof.GetXaxis().GetXmax(), 1.1, r.kGreen, 0.2)
        plot_rsfof.addLine (ttjets_meas.rsfof.GetXaxis().GetXmin(), 1.0, ttjets_meas.rsfof.GetXaxis().GetXmax(), 1.0, r.kGreen)
        plot_rsfof.addLatex(0.2, 0.2, ttjets_meas.cenFwd)
        plot_rsfof.save(1, 1, 0, lumi)
        
