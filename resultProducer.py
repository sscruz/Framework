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
import math, sys, optparse
import Canvas, CutManager, Sample

from ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats

class valErrs:
    def __init__(self, val, sys, stat):
        self.val  = val
        self.sys  = sys
        self.stat = stat
    def err(self):
        return math.sqrt(sys*sys + stat*stat)

class ingredients:
    def readValues(self):
        f = open(self.infile, 'r')
        lines = f.read().splitlines()
        for line in lines:
            if '#' in line: continue

    def __init__(self, infile, isData):
        self.infile = infile
        self.isData = isData
        ## rmue
        self.rmue_sr_lm       = valErrs(-1., -1., -1.)
        self.rmue_sr_onZ      = valErrs(-1., -1., -1.)
        self.rmue_sr_hm       = valErrs(-1., -1., -1.)
        self.rmue_dycr_dym    = valErrs(-1., -1., -1.)

        #rsfof
        self.rsfof_sr_lm      = valErrs(-1., -1., -1.)
        self.rsfof_sr_onZ     = valErrs(-1., -1., -1.)
        self.rsfof_sr_hm      = valErrs(-1., -1., -1.)
        self.rsfof_ttcr_lm    = valErrs(-1., -1., -1.)
        self.rsfof_ttcr_onZ   = valErrs(-1., -1., -1.)
        self.rsfof_ttcr_hm    = valErrs(-1., -1., -1.)

        # RT
        self.rt_region        = valErrs(-1., -1., -1.)

if __name__ == '__main__':

    print 'Starting r_SFOF analysis...'
    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples FilenameWithIngredients', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()


    if len(args) != 2:
      parser.error('wrong number of arguments')

    inputFileName = args[0]
    centralForward = args[1]

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTJets']#, 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C']

    treeMC = Sample.Tree(Sample.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(Sample.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'


    lumi = 0.020
    print 'Running with an integrated luminosity of', lumi,'fb-1'


   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()


    doClosure = True
    isBlinded = True

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
        
