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
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy
import gc, inspect

import include.LeptonSF
import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample



class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts, labelx, logx, protection):

    MC   = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx)
    MCS  = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx)
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx)

    maxValmc = MC.GetMaximum() #GetBinContent(MC.GetMaximumBin())
    maxValdata = DATA.GetMaximum() #GetBinContent(DATA.GetMaximumBin())
    maxVal = max(maxValmc, maxValdata)
    if not logx:
        MC.GetYaxis().SetRangeUser(0.1, 1.3*maxVal)
        MCS.SetMaximum(1.3*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 1.3*maxVal)
    else:
        MC.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
        MCS.SetMaximum(2.0*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)

    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('DataMC/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.6, 0.85, 0.9)
    plot.addStack(MCS, "HIST", 1, 1)
    plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA, MC)


    

if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting r_SFOF analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    DYDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    ttDatasets = ['TTJets_DiLepton']
    mcDatasets = ['ZZTo4L', 'GGHZZ4L',  'WZTo3LNu', 'WWW', 'WWZ','ZZZ', 'tZq_ll','WWTo2L2Nu', 'ZZTo2L2Nu', 'WZTo2L2Q','TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TTTT',  'TTWToQQ', 'TTZToLLNuNu' ,'TTWToLNu', 'WJetsToLNu_LO']
    mcDatasets += ttDatasets
    mcDatasets += DYDatasets

    daDatasetsB = ['DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                   'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                   'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                   'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                   'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376']

    daDatasetsC = ['DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016C_23Sep2016_v1_runs_271036_284044',
                   'MuonEG_Run2016C_23Sep2016_v1_runs_271036_284044']

    daDatasetsD = ['DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part2',
                   'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part1',
                   'MuonEG_Run2016D_23Sep2016_v1_runs_271036_284044']

    daDatasetsE = ['DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part2',
                   'MuonEG_Run2016E_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part1']

    daDatasetsF = ['DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044',
                   'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044',
                   'MuonEG_Run2016F_23Sep2016_v1_runs_271036_284044']

    daDatasetsG = ['DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part2',
                   'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                   'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part3',
                   'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part2',
                   'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                   'MuonEG_Run2016G_23Sep2016_v1_runs_271036_284044']

    daDatasetsH = ['DoubleEG_Run2016H-PromptReco-v2_runs_281207_284035_part1',
                   'DoubleEG_Run2016H-PromptReco-v2_runs_281207_284035_part2',
                   'DoubleEG_Run2016H-PromptReco-v3_runs_284036_284044',
                   'DoubleMuon_Run2016H-PromptReco-v2_runs_281207_284035_part1',
                   'DoubleMuon_Run2016H-PromptReco-v3_runs_284036_284044',
                   'DoubleMuon_Run2016H-PromptReco-v2_runs_281207_284035_part2',
                   'MuonEG_Run2016H-PromptReco-v3_runs_284036_284044',
                   'DoubleMuon_Run2016H-PromptReco-v2_runs_281207_284035_part3',
                   'MuonEG_Run2016H-PromptReco-v2_runs_281207_284035']                      



    daDatasets = daDatasetsB + daDatasetsC + daDatasetsD +daDatasetsE + daDatasetsF + daDatasetsG + daDatasetsH        

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isScan = 0)
    #treeDAB = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsB, 'DA'), 'DATA', 1, isScan = 0)
    #treeDAC = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsC, 'DA'), 'DATA', 1, isScan = 0)
    #treeDAD = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsD, 'DA'), 'DATA', 1, isScan = 0)
    #treeDAE = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsE, 'DA'), 'DATA', 1, isScan = 0)
    #treeDAF = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsF, 'DA'), 'DATA', 1, isScan = 0)
    #treeDAG = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsG, 'DA'), 'DATA', 1, isScan = 0)
    #treeDAH = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsH, 'DA'), 'DATA', 1, isScan = 0)

    lumiB = 5.89
    lumiC = 2.65                                                                                                     
    lumiD = 4.35
    lumiE = 4.05                                                                                                     
    lumiF = 3.16
    lumiG = 7.55                                                                                           
    lumiH = 8.76

    lumi_strB = '5.89invfb'
    lumi_strC = '2.65invfb'
    lumi_strD = '4.35invfb'
    lumi_strE = '4.05invfb'
    lumi_strF = '3.16invfb'
    lumi_strG = '7.55invfb'
    lumi_strH = '8.76invfb'


    lumi = 36.4 ; maxrun = 99276811; lumi_str = '36.4invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelmt2 = 'MT2 [GeV]'
    labelpt2 = 'p_{T} (2nd) [GeV]'
    labelnjet = "N. Jets"



    ###### DY region #######

    #makePlot(lumiB, lumi_strB, treeDAB, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiB, lumi_strB, treeDAB, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiB, lumi_strB, treeDAB, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiC, lumi_strC, treeDAC, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiC, lumi_strC, treeDAC, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiC, lumi_strC, treeDAC, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiC, lumi_strC, treeDAC, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiD, lumi_strD, treeDAD, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiD, lumi_strD, treeDAD, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiD, lumi_strD, treeDAD, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiD, lumi_strD, treeDAD, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiE, lumi_strE, treeDAE, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiE, lumi_strE, treeDAE, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiE, lumi_strE, treeDAE, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiE, lumi_strE, treeDAE, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiF, lumi_strF, treeDAF, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiF, lumi_strF, treeDAF, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiF, lumi_strF, treeDAF, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiF, lumi_strF, treeDAF, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiG, lumi_strG, treeDAG, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiG, lumi_strG, treeDAG, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiG, lumi_strG, treeDAG, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiG, lumi_strG, treeDAG, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiH, lumi_strH, treeDAH, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiH, lumi_strH, treeDAH, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiH, lumi_strH, treeDAH, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    #makePlot(lumiH, lumi_strH, treeDAH, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)




    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
    
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton]), cuts, labelmll,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton]), cuts, labelmll,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton]), cuts, labelmll,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_mm", 20, 20, 300, cuts.AddList([cuts.mm]), cuts, labelmll,  1, 0)

    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "mll_pt2_SF", 8, 20, 200, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.DYControlRegion]), cuts, labelpt2,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "mll_pt2_OF", 8, 20, 200, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.DYControlRegion]), cuts, labelpt2,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "mll_pt2_ee", 8, 20, 200, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.DYControlRegion]), cuts, labelpt2,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "mll_pt2_mm", 8, 20, 200, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.DYControlRegion]), cuts, labelpt2,  1, 0)

    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_SF", 25, 0, 300, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1, 1)
    makePlot(lumi,  lumi_str,  treeDA, treeMC, "met_Edge", "met_ttbar_OF", 25, 0, 300, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_ee", 25, 0, 300, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_mm", 25, 0, 300, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1, 1)

    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.SignalRegion]), cuts, labelmll,  1, 1)
    makePlot(lumi,  lumi_str,  treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.SignalRegion]), cuts, labelmll,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.SignalRegion]), cuts, labelmll,  1, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.SignalRegion]), cuts, labelmll,  1, 1)

    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_signal_SF", 20, 20, 200, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.SignalRegionNoMT2]), cuts, labelmt2,  1, 1)
    makePlot(lumi,  lumi_str,  treeDA, treeMC, "mt2_Edge", "mt2_signal_OF", 20, 20, 200, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.SignalRegionNoMT2]), cuts, labelmt2,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_signal_ee", 20, 20, 200, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.SignalRegionNoMT2]), cuts, labelmt2,  1, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_signal_mm", 20, 20, 200, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.SignalRegionNoMT2]), cuts, labelmt2,  1, 1)

    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_SF", 7, 2, 9, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2]), cuts, labelnjet, 1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_OF", 7, 2, 9, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2]), cuts, labelnjet, 1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_ee", 7, 2, 9, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj2]), cuts, labelnjet, 1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_mm", 7, 2, 9, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj2]), cuts, labelnjet, 1, 0)                                              



 










    # =================================================================

    # r0b1b_of_cuts = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 0', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # r0b1b_of_den_cuts_e1b = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # #makePlot(lumi, lumi_str, treeDA, treeMC, "nBJetMedium25_Edge", "r0b1b_of_nb_btagscaledPOWHEG", 4, -0.5, 3.5, r0b1b_of_cuts, cuts, 'n_{b-jets}', 0, True, True)
    # #makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "r0b1b_of_den_e1b", 14, 60, 130, r0b1b_of_den_cuts_e1b, cuts, labelmll, 0, True)

    # fmll_of_1b_cuts = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # fmll_of_0b_cuts = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 0', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "fmll_of_1b", 14, 60, 200, fmll_of_1b_cuts, cuts, labelmll, 0, True)
    # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "fmll_of_0b", 14, 60, 200, fmll_of_0b_cuts, cuts, labelmll, 0, True)

    #ewino_of_eq1b = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    #a =makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "ewino_of_1b_eq1b", 14, 60, 130, ewino_of_eq1b, cuts, labelmll, 0, True)
    # ewino_sf = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', cuts.Zveto, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # c = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "ewino_sf_1b_zveto", 14, 60, 130, ewino_sf, cuts, labelmll, 0, True)
    # ewino_of = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # b = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "ewino_of_1b", 14, 60, 130, ewino_of, cuts, labelmll, 0, True)

   # nll = '(nll_Edge > 21)'
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF_nll", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.SignalRegionBaseLine, nll]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF_nll", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.SignalRegionBaseLine, nll]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OFSF_nll", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLine]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.SignalRegionBaseLine]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.SignalRegionBaseLine]), cuts, labelmll, 0)
   # tightCharge = '(Lep1_tightCharge_Edge > 0 && Lep2_tightCharge_Edge > 0)'
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF_nll_tightCharge_ttpow", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.SignalRegionBaseLine, nll, tightCharge]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF_nll_tightCharge_ttpow", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.SignalRegionBaseLine, nll, tightCharge]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee_nll_tightCharge_ttpow", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.SignalRegionBaseLine, nll, tightCharge]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm_nll_tightCharge_ttpow", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.SignalRegionBaseLine, nll, tightCharge]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2]), cuts, labelmet, 1)



