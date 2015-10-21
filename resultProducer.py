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

def makeTable(myregion, dataMC, eta): ## for this to make sense the region should be properly binned!!
    line0 = ' %30s &           & '  %('')
    line1 = ' %30s & MC pred.  & ' %('\multirow{3}{*}{%s}' %(myregion.name))
    line2 = ' %30s & MC obs.   & ' %('')
    line3 = ' %30s &  ratio    & ' %('')
    tmp_histo_obs  = myregion.mll     .getHisto(dataMC, eta)
    tmp_histo_pred = myregion.mll_pred.getHisto(dataMC, eta)
    my_range = range(1,tmp_histo_obs.GetNbinsX()+1)
    for i in my_range:
        tmp_ratio   = tmp_histo_obs.GetBinContent(i) / tmp_histo_pred.GetBinContent(i)
        tmp_ratio_e = math.sqrt( (tmp_histo_obs.GetBinError(i)/tmp_histo_obs.GetBinContent(i))**2  +  (tmp_histo_pred.GetBinError(i)/tmp_histo_obs.GetBinContent(i))**2)  * tmp_ratio
        line0 += '%.0f $<$ \\mll $<$ %.0f ' %(tmp_histo_pred.GetXaxis().GetBinLowEdge(i), tmp_histo_pred.GetXaxis().GetBinUpEdge(i))
        line1 += '  %.2f $\\pm$ %.2f      ' %(tmp_histo_pred.GetBinContent(i), tmp_histo_pred.GetBinError(i))
        line2 += '  %.2f $\\pm$ %.2f      ' %(tmp_histo_obs.GetBinContent(i), tmp_histo_obs.GetBinError(i))
        line3 += '  %.2f $\\pm$ %.2f      ' %(tmp_ratio, tmp_ratio_e)
        if i != max(my_range):
            line0+=' & '
            line1+=' & '
            line2+=' & '
            line3+=' & '
        else:
            line0+=' \\\\ '
            line1+=' \\\\ '
            line2+=' \\\\ '
            line3+=' \\\\ '
    return line0, line1, line2, line3

def makePrediction(of_histo, ing, eta):
    central = (eta == 'central')
    sf_pred = copy.deepcopy(of_histo)
    for _bin in range(1,of_histo.GetNbinsX()+1):
        tmp_cont  = of_histo.GetBinContent(_bin)
        tmp_err   = of_histo.GetBinError  (_bin)
        tmp_mass  = of_histo.GetXaxis().GetBinCenter(_bin)
        if   20 <= tmp_mass <=  70.:
            scale = ing.rsfof_ttcr_lm.cen_val   if central else ing.rsfof_ttcr_lm.fwd_val
            s_err = ing.rsfof_ttcr_lm.cen_stat  if central else ing.rsfof_ttcr_lm.fwd_stat
        elif 70 <  tmp_mass <= 120.:
            scale = ing.rsfof_ttcr_onZ.cen_val  if central else ing.rsfof_ttcr_onZ.fwd_val
            s_err = ing.rsfof_ttcr_onZ.cen_stat if central else ing.rsfof_ttcr_onZ.fwd_stat
        elif 120 < tmp_mass:
            scale = ing.rsfof_ttcr_hm.cen_val   if central else ing.rsfof_ttcr_hm.fwd_val
            s_err = ing.rsfof_ttcr_hm.cen_stat  if central else ing.rsfof_ttcr_hm.fwd_stat
        print 'using rsfof: %.4f' %(scale)
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
    mcDatasets = ['TTJets_LO']#, 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
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
                                 [range(20,310,10), range(0,210,10)],
                                 True)
    binnedSR     = Region.region('SR8TeV', 
                                 [cuts.METJetsSignalRegion],
                                 ['mll'],
                                 [ [20., 70., 81., 101., 120., 500.] ],
                                 True)
    regions.append(signalRegion)
    regions.append(binnedSR)


    for reg in regions:
        print 'i am at region %s' %(reg.name)
        for eta in ['central', 'forward']:
            print '... in %s' %(eta)

            cuts_sf = cuts.AddList([cuts.GoodLeptonSF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts+[cuts.trigger])
            cuts_of = cuts.AddList([cuts.GoodLeptonOF()]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts+[cuts.trigger])

            for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):

                if 'mll' in reg.rvars:
                    reg.mll_sf = tree.getTH1F(lumi, 'mll_sf_'+eta, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_sf, '', "m_{ll} (GeV)")
                    reg.mll_of = tree.getTH1F(lumi, 'mll_of_'+eta, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_of, '', "m_{ll} (GeV)")

                    dataMC = 'DATA' if tree == treeDA else 'MC'
                    isData = (dataMC == 'DATA')
                    reg.mll     .setHisto(reg.mll_sf, dataMC, eta)
                    reg.mll_pred.setHisto(makePrediction(reg.mll_of, ingDA if isData else ingMC, eta), dataMC, eta)

                    del reg.mll_sf, reg.mll_of


    tables = []

    for eta in ['central', 'forward']:
        ## =================
        ## MAKE THE Mll PLOT for MC only
        ## =================
        signalRegion.mll_pred.getHisto('MC', eta).GetYaxis().SetRangeUser(0.,  24. if eta == 'central' else  15.)
        signalRegion.mll_pred.getGraph('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.2)
        signalRegion.mll_pred.getGraph('MC', eta).SetFillStyle(3001)
        plot_closure = Canvas.Canvas('results/plot_results_mll_MCClosure_'+eta, 'png,pdf', 0.6, 0.5, 0.80, 0.7)
        plot_closure.addHisto(signalRegion.mll_pred.getHisto('MC', eta), 'hist'   , 'SR - predicted (MC)', 'L', r.kBlue+1 , 1, 0)
        plot_closure.addGraph(signalRegion.mll_pred.getGraph('MC', eta), '2,same' , 'SR - predicted (MC)', 'L', r.kBlue+1  , 1, -1)
        plot_closure.addHisto(signalRegion.mll     .getHisto('MC', eta), 'PE,SAME', 'SR - observed  (MC)', 'PL', r.kBlack, 1, 1)
        plot_closure.addLatex(0.7, 0.8, eta)
        plot_closure.saveRatio(1, 0, 0, lumi, signalRegion.mll.getHisto('MC', eta), signalRegion.mll_pred.getHisto('MC', eta), 0.5, 1.5)

        ## =================
        ## MAKE THE Mll PLOT for DATA MC comparison
        ## =================
        signalRegion.mll_pred.getHisto('MC', eta).GetYaxis().SetRangeUser(0.,  24. if eta == 'central' else  15.)
        signalRegion.mll_pred.getGraph('MC', eta).SetFillColorAlpha(r.kBlue+1, 0.2)
        signalRegion.mll_pred.getGraph('MC', eta).SetFillStyle(3001)
        plot_dataMC = Canvas.Canvas('results/plot_results_mll_dataMC_'+eta, 'png,pdf', 0.6, 0.5, 0.80, 0.7)
        plot_dataMC.addHisto(signalRegion.mll_pred.getHisto('MC', eta)  , 'hist'   , 'SR - predicted (MC)'  , 'L', r.kBlue+1 , 1, 0)
        plot_dataMC.addGraph(signalRegion.mll_pred.getGraph('MC', eta)  , '2,same' , 'SR - predicted (MC)'  , 'L', r.kBlue+1  , 1, -1)
        plot_dataMC.addHisto(signalRegion.mll_pred.getHisto('DATA', eta), 'PE,SAME', 'SR - predicted (DATA)', 'PL', r.kBlack, 1, 1)
        plot_dataMC.addLatex(0.7, 0.8, eta)
        plot_dataMC.saveRatio(1, 0, 0, lumi, signalRegion.mll.getHisto('DATA', eta), signalRegion.mll_pred.getHisto('DATA', eta), 0.5, 1.5)

        tables.append(makeTable(binnedSR, 'MC', eta))

        ## these two only when not blinded!!!
        isBlinded = True
        if not isBlinded:
            ## =================
            ## MAKE THE Mll PLOT for DATA only
            ## =================
            signalRegion.mll_pred.getHisto('DATA', eta).GetYaxis().SetRangeUser(0.,  24. if eta == 'central' else  15.)
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

