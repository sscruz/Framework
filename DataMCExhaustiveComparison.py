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
from itertools import product


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

    print 'getting stack'
    MCS  = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx,'1',0)
    print' getting mc'
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx,'1',0)
    MC = MCS.GetStack().Last()


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
    print "DATA integral ", DATA.Integral()
    print "MC   integral ", MC.Integral()
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = kee
    DATA.SetMarkerStyle(r.kFullCircle)
    plot = Canvas.Canvas('DataMC_trigger/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.6, 0.85, 0.9)
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

    #DYDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    #ttDatasets = ['TTJets_DiLepton']
    #mcDatasets = ['TTJets_SingleLeptonFromT',  'T_tch_powheg', 'TBar_tch_powheg', 'WWTo2L2Nu',  'WZTo3LNu', 'ZZTo4L', 'ZZTo2L2Nu', 'WWW', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll', 'TTZToLLNuNu' , 'TTLLJets_m1to10', 'TTWToLNu', 'TTTT', 'TTHnobb_pow', 'GGHZZ4L',  'WJetsToLNu_LO']
    #mcDatasets += ttDatasets
    #mcDatasets += DYDatasets

    dyDatasets = ['DYJetsToLL_M50_ext_part1+DYJetsToLL_M50_ext_part2+DYJetsToLL_M50_ext_part3','DYJetsToLL_M10to50']
    ttDatasets = ['TTTo2L2Nu_part1+TTTo2L2Nu_part2','TTToSemiLeptonic'] # tt1l missing
    stDatasets = ['TW','TbarW'] # t and s (lol) channel missing 
    ttzDatasets = ['TTZ_LO_ext1','TTW_LO']#,'TTWZ','TTGJets_newpmx']
    zz2lDatasets = ['ZZTo2L2Q','ZZTo2L2Nu']
    zz4lDatasets = ['ZZTo4L_ext1_part1+ZZTo4L_ext1_part2',
                    #'GluGluToContinToZZTo2e2mu+GluGluToContinToZZTo2e2mu_ext1',
                    #'GluGluToContinToZZTo2e2nu+GluGluToContinToZZTo2e2nu_ext1',
                    #'GluGluToContinToZZTo2mu2nu+GluGluToContinToZZTo2mu2nu_ext1',
                    #'GluGluToContinToZZTo4e+GluGluToContinToZZTo4e_ext1',
                    #'GluGluToContinToZZTo4mu+GluGluToContinToZZTo4mu_ext1'
    ]
    wwDatasets = ['WWTo2L2Nu']
    wzDatasets = ['WZTo3LNu','WZTo2L2Q']
    raDatasets = [ 'TTH_amc']
    mcDatasets = zz4lDatasets + zz2lDatasets + ttzDatasets + raDatasets + wwDatasets +wzDatasets + stDatasets+  ttDatasets + dyDatasets
    
    daDatasets = [ '%s_Run2017%s'%(x,y) for x,y in product('DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(','), 'B,C,D,E,F'.split(','))]
    daDatasetsB = [ '%s_Run2017B'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]
    daDatasetsC = [ '%s_Run2017C'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]
    daDatasetsD = [ '%s_Run2017D'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]
    daDatasetsE = [ '%s_Run2017E'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]
    daDatasetsF = [ '%s_Run2017F'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isScan = 0)
    treeDAB = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsB, 'DA'), 'DATA', 1, isScan = 0)
    treeDAC = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsC, 'DA'), 'DATA', 1, isScan = 0)
    treeDAD = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsD, 'DA'), 'DATA', 1, isScan = 0)
    treeDAE = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsE, 'DA'), 'DATA', 1, isScan = 0)
    treeDAF = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsF, 'DA'), 'DATA', 1, isScan = 0)
    #treeDAG = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsG, 'DA'), 'DATA', 1, isScan = 0)
    #treeDAH = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsH, 'DA'), 'DATA', 1, isScan = 0)



    lumi = 41.4 ; maxrun = 999999999999999999999; lumi_str = '41.4invfb'
    
    lumiB = 5.89; lumi_strB = '5.89invfb'
    lumiC = 2.65; lumi_strC = '2.65invfb'                                                                                                     
    lumiD = 4.35; lumi_strD = '4.35invfb'
    lumiE = 4.05; lumi_strE = '4.05invfb'                                                                                                     
    lumiF = 3.16; lumi_strF = '3.16invfb'
    
    
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelmt2 = 'MT2 [GeV]'
    labelpt2 = 'p_{T} (2nd) [GeV]'
    labelpt1 = 'p_{T} (1nd) [GeV]'
    labelnjet = "N. Jets"
    labelnbjet = "N. B Jets"
    labelnbjetnjet = "(N_{jets},N_{b-jets})"
    labelfatjet = "N. AK8 jets" 



    makePlot(lumi, lumi_str, treeDA, treeMC, "MET_pt_Edge", "MET_DYJets_OF", 20, 0, 500, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, labelmet,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "MET_pt_Edge", "MET_DYJets_ee", 20, 0, 500, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, labelmet,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "MET_pt_Edge", "MET_DYJets_mm", 20, 0, 500, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, labelmet,  1, 0)

    ###### DY region #######
    makePlot(lumi, lumi_str, treeDA, treeMC, "njnb(nJetSel_Edge, nBJetMedium25_Edge)", "nbnj_DYJets_OF", 7, -0.5,6.5, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, labelnbjetnjet,  1, 0)

    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, labelmll,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, labelmll,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, labelmll,  1, 0)

    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep1_pt_Edge", "lep1_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, labelpt1,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep1_pt_Edge", "lep1_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, labelpt1,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "Lep1_pt_Edge", "lep1_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, labelpt1,  1, 0)

    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_DYJets_OF", 10,-0.4,9.5, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, labelnjet,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_DYJets_ee", 10,-0.4,9.5, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, labelnjet,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_DYJets_mm", 10,-0.4,9.5, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, labelnjet,  1, 0)


    makePlot(lumi, lumi_str, treeDA, treeMC, "njnb(nJetSel_Edge, nBJetMedium25_Edge)", "nbnj_DYJets_ee", 7, -0.5,6.5, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, labelnbjetnjet,  1, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "njnb(nJetSel_Edge, nBJetMedium25_Edge)", "nbnj_DYJets_mm", 7, -0.5,6.5, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, labelnbjetnjet,  1, 0)
    
    makePlot(lumi, lumi_str, treeDA, treeMC, "nFatJetSel_Edge", "nFatJet_OF", 5,-0.5,4.5, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, labelfatjet,  1, 0)

    

    # makePlot(lumi, lumi_str, treeDA, treeMC, "1", "tt_xsec_OF", 7, -0.5,6.5, cuts.AddList([cuts.OF, cuts.goodLepton17, 'nJetSel_Edge > 1', 'nBJetMedium25_Edge > 0']), cuts, labelnbjetnjet,  1, 0)
    # makePlot(lumi, lumi_str, treeDA, treeMC, "1", "tt_xsec_ee", 7, -0.5,6.5, cuts.AddList([cuts.ee, cuts.goodLepton17, 'nJetSel_Edge > 1', 'nBJetMedium25_Edge > 0','abs(lepsMll_Edge - 91) > 15', 'MET_pt_Edge > 50']), cuts, labelnbjetnjet,  1, 0)
    # makePlot(lumi, lumi_str, treeDA, treeMC, "1", "tt_xsec_mm", 7, -0.5,6.5, cuts.AddList([cuts.mm, cuts.goodLepton17, 'nJetSel_Edge > 1', 'nBJetMedium25_Edge > 0','abs(lepsMll_Edge - 91) > 15', 'MET_pt_Edge > 50']), cuts, labelnbjetnjet,  1, 0)





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



    #bins = [20, 60, 86, 96, 150, 200, 300, 400]
    #makePlot(lumi, lumi_str, treeMC, treeMC, "lepsMll_Edge", "rares_SF", bins, 1, 1, cuts.AddList([cuts.goodLepton, cuts.JETMETBaseline,  cuts.METg150,  cuts.SF]), cuts, labelmll,  1, 0)                           
    #makePlot(lumi, lumi_str, treeMC, treeMC, "lepsMll_Edge", "rares_SF_fromZ", bins, 1, 1, cuts.AddList([cuts.goodLepton, cuts.JETMETBaseline, cuts.lepsFromZ,  cuts.METg150,  cuts.SF]), cuts, labelmll,  1, 0)                           
    #makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_rt_SF", 20, 0, 150, cuts.AddList([cuts.goodLeptonNoTrigger, cuts.donot(cuts.JETMETBaseline),cuts.HT,  cuts.donot(cuts.RSFOFDirectControlRegion),  cuts.trigger, cuts.SF]), cuts, labelmll,  1, 0)                           
    #makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_rt_OF", 20, 0, 150, cuts.AddList([cuts.goodLeptonNoTrigger, cuts.donot(cuts.JETMETBaseline),cuts.HT,  cuts.donot(cuts.RSFOFDirectControlRegion),  cuts.trigger, cuts.OF]), cuts, labelmll,  1, 0)                           
    #makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_rt_ee", 20, 0, 150, cuts.AddList([cuts.goodLeptonNoTrigger, cuts.donot(cuts.JETMETBaseline),cuts.HT,  cuts.donot(cuts.RSFOFDirectControlRegion),  cuts.trigger, cuts.ee]), cuts, labelmll,  1, 0)                           
    #makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_rt_mm", 20, 0, 150, cuts.AddList([cuts.goodLeptonNoTrigger, cuts.donot(cuts.JETMETBaseline),cuts.HT,  cuts.donot(cuts.RSFOFDirectControlRegion),  cuts.trigger, cuts.mm]), cuts, labelmll,  1, 0)                           
    #makePlot(lumi, lumi_str, treeDA, treeMC, "Lep1_pt_Edge", "pt1_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)                           
    #makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "pt2_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)                           
    #makePlot(lumi, lumi_str, treeDA, treeMC, "Lep1_eta_Edge", "eta1_DYJets_mm", 20, -2.4, 2.4, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)                       
    #makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_eta_Edge", "eta2_DYJets_mm", 20, -2.4, 2.4, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)




















   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.DYControlRegionNoMll]), cuts, labelmll,  1, 0)
   # 
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsZPt_Edge", "ZPt_incl_SF", 20, 0, 200, cuts.AddList([cuts.SF, cuts.goodLepton]), cuts, "Z p_{T}")
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsZPt_Edge", "ZPt_incl_OF", 20, 0, 200, cuts.AddList([cuts.OF, cuts.goodLepton]), cuts, "Z p_{T}")
    #makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton]), cuts, labelmll,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton]), cuts, labelmll,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton]), cuts, labelmll,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_mm", 20, 20, 300, cuts.AddList([cuts.mm]), cuts, labelmll,  1, 0)

   # makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "mll_pt2_SF", 8, 20, 200, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.DYControlRegion]), cuts, labelpt2,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "mll_pt2_OF", 8, 20, 200, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.DYControlRegion]), cuts, labelpt2,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "mll_pt2_ee", 8, 20, 200, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.DYControlRegion]), cuts, labelpt2,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "Lep2_pt_Edge", "mll_pt2_mm", 8, 20, 200, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.DYControlRegion]), cuts, labelpt2,  1, 0)

    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsZPt_Edge", "ZPt_ttbar_SF", 20, 0, 200, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.RSFOFDirectControlRegion]), cuts, "Z p_{T}")
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsZPt_Edge", "ZPt_ttbar_OF", 20, 0, 200, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.RSFOFDirectControlRegion]), cuts, "Z p_{T}")
    #makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_SF", 25, 0, 300, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1, 1)
   # makePlot(lumi,  lumi_str,  treeDA, treeMC, "met_Edge", "met_ttbar_OF", 25, 0, 300, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_ee", 25, 0, 300, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1, 1)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_mm", 25, 0, 300, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1, 1)

    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsZPt_Edge", "ZPt_signal_SF", 20, 0, 200, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.SignalRegion]), cuts, "Z p_{T}")
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsZPt_Edge", "ZPt_signal_OF", 20, 0, 200, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.SignalRegion]), cuts, "Z p_{T}")
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.SignalRegion]), cuts, labelmll,  1, 1)
   # makePlot(lumi,  lumi_str,  treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.SignalRegion]), cuts, labelmll,  1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.SignalRegion]), cuts, labelmll,  1, 1)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.SignalRegion]), cuts, labelmll,  1, 1)

    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_signal_SF", 20, 20, 200, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.SignalRegionNoMT2]), cuts, labelmt2)
    makePlot(lumi,  lumi_str,  treeDA, treeMC, "mt2_Edge", "mt2_signal_OF", 20, 20, 200, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.SignalRegionNoMT2]), cuts, labelmt2)
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_signal_ee", 20, 20, 200, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.SignalRegionNoMT2]), cuts, labelmt2)
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_signal_mm", 20, 20, 200, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.SignalRegionNoMT2]), cuts, labelmt2)   
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_ttbar_SF", 20, 20, 200, cuts.AddList([cuts.SF, cuts.goodLepton, cuts.RSFOFDirectControlRegion]), cuts, labelmt2)
    makePlot(lumi,  lumi_str,  treeDA, treeMC, "mt2_Edge", "mt2_ttbar_OF", 20, 20, 200, cuts.AddList([cuts.OF, cuts.goodLepton, cuts.RSFOFDirectControlRegion]), cuts, labelmt2)
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_ttbar_ee", 20, 20, 200, cuts.AddList([cuts.ee, cuts.goodLepton, cuts.RSFOFDirectControlRegion]), cuts, labelmt2)
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_ttbar_mm", 20, 20, 200, cuts.AddList([cuts.mm, cuts.goodLepton, cuts.RSFOFDirectControlRegion]), cuts, labelmt2)   
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_incl_SF", 20, 20, 200, cuts.AddList([cuts.SF, cuts.goodLepton]), cuts, labelmt2)
    makePlot(lumi,  lumi_str,  treeDA, treeMC, "mt2_Edge", "mt2_incl_OF", 20, 20, 200, cuts.AddList([cuts.OF, cuts.goodLepton]), cuts, labelmt2)
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_incl_ee", 20, 20, 200, cuts.AddList([cuts.ee, cuts.goodLepton]), cuts, labelmt2)
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_incl_mm", 20, 20, 200, cuts.AddList([cuts.mm, cuts.goodLepton]), cuts, labelmt2)   

   # makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_SF", 7, 2, 9, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2]), cuts, labelnjet, 1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_OF", 7, 2, 9, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2]), cuts, labelnjet, 1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_ee", 7, 2, 9, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj2]), cuts, labelnjet, 1, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_mm", 7, 2, 9, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj2]), cuts, labelnjet, 1, 0)                                              



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



