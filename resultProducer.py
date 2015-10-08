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
    central = (eta == 'central')
    sf_pred = copy.deepcopy(of_histo)
    for _bin in range(1,of_histo.GetNbinsX()+1):
        tmp_cont  = of_histo.GetBinContent(_bin)
        tmp_err   = of_histo.GetBinError  (_bin)
        tmp_mass  = of_histo.GetXaxis().GetBinCenter(_bin)
        if   20 <= tmp_mass <=  70.:
            scale = ing.rsfof_ttcr_lm.cen_val  if central else ing.rsfof_ttcr_lm.fwd_val
            s_err = ing.rsfof_ttcr_lm.cen_err  if central else ing.rsfof_ttcr_lm.fwd_err
        elif 70 <  tmp_mass <= 120.:
            scale = ing.rsfof_ttcr_onZ.cen_val if central else ing.rsfof_ttcr_onZ.fwd_val
            s_err = ing.rsfof_ttcr_onZ.cen_err if central else ing.rsfof_ttcr_onZ.fwd_err
        elif 120 < tmp_mass:
            scale = ing.rsfof_ttcr_hm.cen_val  if central else ing.rsfof_ttcr_hm.fwd_val
            s_err = ing.rsfof_ttcr_hm.cen_err  if central else ing.rsfof_ttcr_hm.fwd_err
        tmp_pred  = scale*tmp_cont
        tmp_prede = math.sqrt(tmp_err*tmp_err + s_err*s_err*tmp_cont*tmp_cont)
        sf_pred.SetBinContent(_bin, tmp_pred )
        sf_pred.SetBinError  (_bin, tmp_prede)
    return sf_pred



if __name__ == '__main__':

    print 'Starting to produce some good ol\' results...'
    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples FilenameWithIngredients', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()


    if len(args) != 2:
      parser.error('wrong number of arguments')

    sampleFile     = args[0]
    ingredientFile = args[1]

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTJets']#, 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C',
                  'DoubleMuon_Run2015D', 'DoubleEG_Run2015D', 'MuonEG_Run2015D']

    treeMC = Sample.Tree(helper.selectSamples(sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'


    lumi = 0.225
    print 'Running with an integrated luminosity of %.2f fb-1' %(lumi)

    isBlinded = False


   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()


    ingMC = helper.ingredients(ingredientFile, 'MC'  )
    ingDA = helper.ingredients(ingredientFile, 'DATA')


    regions = []
    signalRegion = Region.region('signalRegion', 
                                 [cuts.METJetsSignalRegion],
                                 ['mll', 'met'],
                                 [range(20,310, 5), range(0,210,10)],
                                 True)
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
                    isData = (dataMC == 'DATA')
                    reg.mll     .setHisto(reg.mll_sf, dataMC, eta)
                    reg.mll_pred.setHisto(makePrediction(reg.mll_of, ingDA if isData else ingMC, eta), dataMC, eta)

                    del reg.mll_sf, reg.mll_of



    for eta in ['central', 'forward']:
        ## =================
        ## MAKE THE Mll PLOT for MC only
        ## =================
        signalRegion.mll_pred.getHisto('MC', eta).GetYaxis().SetRangeUser(0.,  8. if eta == 'central' else  5.)
        signalRegion.mll_pred.getGraph('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.2)
        signalRegion.mll_pred.getGraph('MC', eta).SetFillStyle(3001)
        plot_closure = Canvas.Canvas('results/plot_results_mll_MCClosure_'+eta, 'png,pdf', 0.6, 0.5, 0.80, 0.7)
        plot_closure.addHisto(signalRegion.mll_pred.getHisto('MC', eta), 'hist'   , 'SR - predicted', 'L', r.kBlue+1 , 1, 0)
        plot_closure.addGraph(signalRegion.mll_pred.getGraph('MC', eta), '2,same' , 'SR - predicted', 'L', r.kBlue+1  , 1, -1)
        plot_closure.addHisto(signalRegion.mll     .getHisto('MC', eta), 'PE,SAME', 'SR - observed' , 'PL', r.kBlack, 1, 1)
        plot_closure.addLatex(0.7, 0.8, eta)
        plot_closure.saveRatio(1, 0, 0, lumi, signalRegion.mll.getHisto('MC', eta), signalRegion.mll_pred.getHisto('MC', eta), 0.5, 1.5)

        ## these two only when not blinded!!!
        if not isBlinded:
            ## =================
            ## MAKE THE Mll PLOT for DATA MC comparison
            ## =================
            signalRegion.mll_pred.getHisto('MC', eta).GetYaxis().SetRangeUser(0.,  8. if eta == 'central' else  5.)
            signalRegion.mll_pred.getGraph('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.2)
            signalRegion.mll_pred.getGraph('MC', eta).SetFillStyle(3001)
            plot_dataMC = Canvas.Canvas('results/plot_results_mll_dataMC_'+eta, 'png,pdf', 0.6, 0.5, 0.80, 0.7)
            plot_dataMC.addHisto(signalRegion.mll_pred.getHisto('MC', eta)  , 'hist'   , 'SR - predicted (MC)'  , 'L', r.kBlue+1 , 1, 0)
            plot_dataMC.addGraph(signalRegion.mll_pred.getGraph('MC', eta)  , '2,same' , 'SR - predicted (MC)'  , 'L', r.kBlue+1  , 1, -1)
            plot_dataMC.addHisto(signalRegion.mll_pred.getHisto('DATA', eta), 'PE,SAME', 'SR - predicted (DATA)', 'PL', r.kBlack, 1, 1)
            plot_dataMC.addLatex(0.7, 0.8, eta)
            plot_dataMC.saveRatio(1, 0, 0, lumi, signalRegion.mll.getHisto('DATA', eta), signalRegion.mll_pred.getHisto('DATA', eta), 0.5, 1.5)


            ## =================
            ## MAKE THE Mll PLOT for DATA only
            ## =================
            signalRegion.mll_pred.getHisto('DATA', eta).GetYaxis().SetRangeUser(0.,  8. if eta == 'central' else  5.)
            signalRegion.mll_pred.getGraph('DATA', eta).SetFillColorAlpha(r.kBlue+1, 0.2)
            signalRegion.mll_pred.getGraph('DATA', eta).SetFillStyle(3001)
            plot_result = Canvas.Canvas('results/plot_results_mll_data_'+eta, 'png,pdf', 0.6, 0.5, 0.80, 0.7)
            plot_result.addHisto(signalRegion.mll_pred.getHisto('DATA', eta), 'hist'   , 'SR - predicted (DATA)', 'L', r.kBlue+1 , 1, 0)
            plot_result.addGraph(signalRegion.mll_pred.getGraph('DATA', eta), '2,same' , 'SR - predicted (DATA)', 'L', r.kBlue+1  , 1, -1)
            plot_result.addHisto(signalRegion.mll     .getHisto('DATA', eta), 'PE,SAME', 'SR - observed  (DATA)' , 'PL', r.kBlack, 1, 1)
            plot_result.addLatex(0.7, 0.8, eta)
            plot_result.saveRatio(1, 0, 0, lumi, signalRegion.mll.getHisto('DATA', eta), signalRegion.mll_pred.getHisto('DATA', eta), 0.5, 1.5)
        else:
            print 'you sneaky bastard. stop this!!!'
            sys.exit('...exiting!')

