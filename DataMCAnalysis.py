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
import math, sys, optparse, array
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


def dumpObjects():
    gc.collect()
    oo = gc.get_objects()
    for o in oo:
        if getattr(o, "__class__", None):
            name = o.__class__.__name__
            if name not in exclude:
                filename = inspect.getabsfile(o.__class__)            
                print "Object of class:", name, "...",
                print "defined in file:", filename                

def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx):

    theCutDATA = cuts.AddList([theCut, cuts.trigger])
    MC = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx)
    MCS = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx)
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCutDATA, '', labelx)

    maxValmc = MC.GetBinContent(MC.GetMaximumBin())
    maxValdata = DATA.GetBinContent(DATA.GetMaximumBin())
    maxVal = max(maxValmc, maxValdata)
    if logx == 0:
        MC.GetYaxis().SetRangeUser(0, 1.3*maxVal)
        MCS.GetYaxis().SetRangeUser(0, 1.3*maxVal)
        DATA.GetYaxis().SetRangeUser(0, 1.3*maxVal)
    else:
        MC.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
        MCS.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
   
    print name
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('DataMC/%s/plot_%s'%(lumi_str,name), 'png,pdf', 0.7, 0.7, 0.95, 0.9)
    plot.addStack(MCS, "HIST", 1, 1)
    plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA, MC)
    
    del plot
    return



if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting r_SFOF analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    mcDatasets = ['TTJets_DiLepton', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50', 'WW', 'WZ', 'ZZ']
    daDatasets = ['DoubleMuon_Run2016B_PromptReco_v2_runs_273150_273730', 'DoubleEG_Run2016B_PromptReco_v2_runs_273150_273730', 'MuonEG_Run2016B_PromptReco_v2_runs_273150_273730']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 0.589 ; maxrun = 999999
    lumi_str = '0.6invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'MET [GeV]'
    labelnjet = "N. Jets"

    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_SF", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, labelmll, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_OF", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, labelmll, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_ee", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, labelmll, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_mm", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, labelmll, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.JetMETBaseline]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.JetMETBaseline]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.JetMETBaseline]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.JetMETBaseline]), cuts, labelmll, 0)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj2]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj2]), cuts, labelmet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_SF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.SF]), cuts, labelnjet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_OF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.OF]), cuts, labelnjet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_ee", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.ee]), cuts, labelnjet, 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_mm", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.mm]), cuts, labelnjet, 1)

 

      
