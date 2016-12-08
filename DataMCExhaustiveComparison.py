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


def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts, labelx, logx):

    theCutDATA = cuts.AddList([theCut, cuts.trigger ])
    MC   = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx)
    MCS  = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx)
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCutDATA, '', labelx)

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
    #dyDatasets = ['DYJetsToLL_M10to50','DYJetsToLL_M50_HT100to200', 'DYJetsToLL_M50_HT200to400', 'DYJetsToLL_M50_HT400to600', 'DYJetsToLL_M50_HT600toInf']
    #ttDatasets = ['TTJets_DiLepton', 'TTJets_DiLepton_ext']
    ttDatasets = ['TTJets_DiLepton']
    #mcDatasets = ['ZZTo4L', 'WZTo3LNu', 'WWW', 'WWZ','WZZ', 'ZZZ',  'TTZToLLNuNu' ,'TTWToLNu', 'T_tWch', 'TBar_tWch' , 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TToLeptons_sch', 'TToLeptons_tch_powheg', 'TBarToLeptons_tch_powheg', 'TTHnobb_pow', 'VHToNonbb', 'WJetsToLNu_LO']
    mcDatasets = ['ZZTo4L', 'WZTo3LNu', 'WWW', 'WWZ','WZZ', 'ZZZ',  'TTZToLLNuNu' ,'TTWToLNu', 'T_tWch', 'TBar_tWch' , 'TToLeptons_sch', 'TToLeptons_tch_powheg', 'TBarToLeptons_tch_powheg', 'TTHnobb_pow', 'VHToNonbb', 'WJetsToLNu_LO']
    mcDatasets += ttDatasets
    mcDatasets += DYDatasets
    daDatasets = ['DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part5',
                  'MuonEG_Run2016F_23Sep2016_v1_runs_271036_284044',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part3',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part4',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part5',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part10',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part11',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part3',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part4',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part5',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part7',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part1',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part4',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part5',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part6',
                  'DoubleEG_Run2016H-PromptReco-v3_runs_284036_284044',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part1',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part10',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part3',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part4',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part5',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part7',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part8',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part9',
                  'DoubleMuon_Run2016H-PromptReco-v3_runs_284036_284044',
                   #'MuonEG_Run2016H-PromptReco-v2_runs_281613_284035',
                   #'MuonEG_Run2016H-PromptReco-v3_runs_284036_284044',
                  'MuonEG_Run2016H-PromptReco-v2_runs_281613_284035',
                  'MuonEG_Run2016H-PromptReco-v3_runs_284036_284044',
                  'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                  'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                  'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part3',
                  'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part4',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part8',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part9',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part2',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part6',
                  'DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044_part4',
                  #'DoubleMuon_Run2016C_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016C_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016C_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part5',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part6',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part7',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part5',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part6',
                  #'MuonEG_Run2016C_23Sep2016_v1_runs_271036_284044',
                  #'MuonEG_Run2016D_23Sep2016_v1_runs_271036_284044',
                  #'MuonEG_Run2016E_23Sep2016_v1_runs_271036_284044',
                  'MuonEG_Run2016C_23Sep2016_v1_runs_271036_284044',
                  'MuonEG_Run2016D_23Sep2016_v1_runs_271036_284044',
                  'MuonEG_Run2016E_23Sep2016_v1_runs_271036_284044',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part2',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part3',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part6',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part5',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part6',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part5',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part6',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part7',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part8',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part9',
                  #'MuonEG_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                  #'MuonEG_Run2016G_23Sep2016_v1_runs_271036_284044_part2']
                  'MuonEG_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                  'MuonEG_Run2016G_23Sep2016_v1_runs_271036_284044_part2']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isScan = 0)

    lumi = 36.4 ; maxrun = 276811; lumi_str = '36.4invfb'
    lumir = 17.0 ; maxrun = 276811; lumi_strr = '17.0invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelnjet = "N. Jets"



    ###### DY region #######

    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  1)

    makePlot(lumir, lumi_strr, treeDA, treeMC, "met_Edge", "met_ttbar_SF", 25, 0, 200, cuts.AddList([cuts.SF, cuts.RSFOFDirectControlRegionNoMET, cuts.protection]), cuts, labelmet, 1)
    makePlot(lumi,  lumi_str,  treeDA, treeMC, "met_Edge", "met_ttbar_OF", 25, 0, 200, cuts.AddList([cuts.OF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet,  1)
    makePlot(lumir, lumi_strr, treeDA, treeMC, "met_Edge", "met_ttbar_ee", 25, 0, 200, cuts.AddList([cuts.ee, cuts.RSFOFDirectControlRegionNoMET, cuts.protection]), cuts, labelmet,  1)
    makePlot(lumir, lumi_strr, treeDA, treeMC, "met_Edge", "met_ttbar_mm", 25, 0, 200, cuts.AddList([cuts.mm, cuts.RSFOFDirectControlRegionNoMET, cuts.protection]), cuts, labelmet,  1)

    makePlot(lumir, lumi_strr, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.SignalRegion]), cuts, labelmll,  1)
    makePlot(lumi,  lumi_str,  treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.SignalRegion]), cuts, labelmll,  1)
    makePlot(lumir, lumi_strr, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.SignalRegion]), cuts, labelmll,  1)
    makePlot(lumir, lumi_strr, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.SignalRegion]), cuts, labelmll,  1)

    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_SF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.SF]), cuts, labelnjet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_OF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.OF]), cuts, labelnjet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_ee", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.ee]), cuts, labelnjet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_mm", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.mm]), cuts, labelnjet, 1)                                              



 










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



