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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample


def makePrediction(of_histo, ing, eta):
    sf_pred = copy.deepcopy(of_histo)
    for _bin in range(1,of_histo.GetNbinsX()+1):
        tmp_cont  = of_histo.GetBinContent(_bin)
        tmp_err   = of_histo.GetBinError  (_bin)
        tmp_pred  = 1.2*tmp_cont
        tmp_prede = 1.2*tmp_err
        sf_pred.SetBinContent(_bin, tmp_pred )
        sf_pred.SetBinError  (_bin, tmp_prede)
    return sf_pred



if __name__ == '__main__':

    print 'Starting r_SFOF analysis...'
    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples FilenameWithIngredients', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()


    if len(args) != 2:
      parser.error('wrong number of arguments')

    sampleFile     = args[0]
    ingredientFile = args[1]

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTJets']#, 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C']

    treeMC = Sample.Tree(helper.selectSamples(sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'


    lumi = 0.150
    print 'Running with an integrated luminosity of %.2f fb-1' %(lumi)


   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()


    doClosure = True
    isBlinded = True

    ingMC = helper.ingredients(ingredientFile, 'MC'  )
    ingDA = helper.ingredients(ingredientFile, 'DATA')


    regions = []
    signalRegion = Region.region('signalRegion', 
                                 [cuts.METJetsSignalRegion],
                                 ['mll', 'met'],
                                 [range(0,310,10), range(0,210,10)],
                                 True if not isBlinded else False)
    regions.append(signalRegion)


    for reg in regions:
        print 'i am at region %s' %(reg.name)
        for eta in ['central', 'forward']:
            print '... in %s' %(eta)

            cuts_sf = cuts.AddList([cuts.GoodLeptonSF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
            cuts_of = cuts.AddList([cuts.GoodLeptonOF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)

            for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):

                if 'mll' in reg.rvars:
                    reg.mll_sf = tree.getTH1F(lumi, 'mll_sf_'+eta, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_sf, '', "m_{ll} (GeV)")
                    reg.mll_of = tree.getTH1F(lumi, 'mll_of_'+eta, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_of, '', "m_{ll} (GeV)")

                    dataMC = 'DATA' if tree == treeDA else 'MC'
                    reg.mll     .setHisto(reg.mll_sf, dataMC, eta)
                    reg.mll_pred.setHisto(makePrediction(reg.mll_of, ingMC, eta), dataMC, eta)

                    del reg.mll_sf, reg.mll_of



    for eta in ['central', 'forward']:
        ## =================
        ## MAKE THE Mll PLOT
        ## =================
        signalRegion.mll_pred.getHisto('MC', eta).GetYaxis().SetRangeUser(0., 15. if eta == 'central' else 10.)
        plot_result = Canvas.Canvas('results/plot_results_mll_'+eta, 'png,pdf', 0.6, 0.5, 0.80, 0.7)
        plot_result.addHisto(signalRegion.mll_pred.getHisto('MC', eta), 'PE'     , 'SR - predicted', 'PL', r.kRed+1 , 1, 0)
        plot_result.addHisto(signalRegion.mll     .getHisto('MC', eta), 'PE,SAME', 'SR - observed' , 'PL', r.kBlack , 1, 1)
        plot_result.addLatex(0.2, 0.2, eta)
        plot_result.saveRatio(1, 1, 0, lumi, signalRegion.mll_pred.getHisto('MC', eta), signalRegion.mll.getHisto('MC', eta))
