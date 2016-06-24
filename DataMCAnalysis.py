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

def makeSummaryTable3l4l(plot_3l, plot_4l):
    wz3l, wz3l_e, zz4l, zz4l_e = 0., 0., 0., 0.,
    bkg3l, bkg3l_e, bkg4l, bkg4l_e = 0., 0., 0., 0.,
    obs3l, obs4l = 0, 0
    rat3l, rat3l_e, rat4l, rat4l_e = 0., 0., 0., 0.,
    
    for i, h_3l in enumerate(plot_3l.histos):
        if not i: continue
        tmp_double = r.Double()
        if 'WZ' in h_3l.GetName():
            wz3l = h_3l.IntegralAndError(1, h_3l.GetNbinsX()+1, tmp_double)
            wz3l_e = tmp_double
        elif not 'DATA' in h_3l.GetName(): 
            tmp_yield = h_3l.IntegralAndError(1, h_3l.GetNbinsX()+1, tmp_double)
            bkg3l  += tmp_yield
            bkg3l_e = math.sqrt(bkg3l_e**2 + tmp_double**2)
        elif 'DATA' in h_3l.GetName():
            obs3l = h_3l.Integral()
        
    for i, h_4l in enumerate(plot_4l.histos):
        if not i: continue
        tmp_double = r.Double()
        if 'ZZ' in h_4l.GetName():
            zz4l = h_4l.IntegralAndError(1, h_4l.GetNbinsX()+1, tmp_double)
            zz4l_e = tmp_double
        elif not 'DATA' in h_4l.GetName(): 
            tmp_yield = h_4l.IntegralAndError(1, h_4l.GetNbinsX()+1, tmp_double)
            bkg4l  += tmp_yield
            bkg4l_e = math.sqrt(bkg4l_e**2 + tmp_double**2)
        elif 'DATA' in h_4l.GetName():
            obs4l = h_4l.Integral()
        
    obsSub3l, obsSub3l_e, obsSub4l, obsSub4l_e = obs3l-bkg3l, math.sqrt(bkg3l_e**2 + obs3l), obs4l-bkg4l, math.sqrt(bkg4l_e**2 + obs4l)
    rat3l, rat3l_e = helper.ratioError(wz3l, wz3l_e, obsSub3l, obsSub3l_e)
    rat4l, rat4l_e = helper.ratioError(zz4l, zz4l_e, obsSub4l, obsSub4l_e)

    table3l4lComparisonString = '''\\documentclass[12pt,a4paper]{{article}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{3-lepton and 4-lepton control regions. Signal MC is WZ-$>$3lnu in the 3-lepton region and ZZ-$>$4l in the 4-lepton region.}} 
\\label{{tab:tab3l4lcontrol}} 
\\begin{{tabular}}{{l| c c }} 
              & 3-lepton region & 4-lepton region \\\\ \\hline \\hline
signal MC        & {wz3l:.2f}     $\\pm$  {wz3l_e:.2f}       & {zz4l:.2f}     $\\pm$  {zz4l_e:.2f}       \\\\ \\hline
bkg. MC          & {bkg3l:.2f}  $\\pm$  {bkg3l_e:.2f}        & {bkg4l:.2f}  $\\pm$  {bkg4l_e:.2f}    \\\\ \\hline \\hline
\\textbf{{data}}             & \\textbf{{{obs3l}}}                               & \\textbf{{{obs4l}}}      \\\\ \\hline
data-bkg.        & {obsSub3l:.2f}   $\\pm$  {obsSub3l_e:.2f} & {obsSub4l:.2f}   $\\pm$  {obsSub4l_e:.2f}     \\\\ \\hline \\hline
sig./(data-bkg.) & {rat3l:.2f}   $\\pm$  {rat3l_e:.2f}       & {rat4l:.2f}   $\\pm$  {rat4l_e:.2f}     \\\\ \\hline

\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(wz3l=wz3l, wz3l_e=wz3l_e, zz4l=zz4l, zz4l_e=zz4l_e,
                            bkg3l=bkg3l, bkg3l_e=bkg3l_e, bkg4l=bkg4l,bkg4l_e=bkg4l_e, 
                            obs3l=int(obs3l),obs4l=int(obs4l),
                            obsSub3l=obsSub3l, obsSub3l_e=obsSub3l_e, obsSub4l=obsSub4l, obsSub4l_e=obsSub4l_e,
                            rat3l=rat3l, rat3l_e=rat3l_e, rat4l=rat4l, rat4l_e=rat4l_e)
    compTableFile = open('plots/DataMC/tables/controlRegion3l4l_comparison.tex','w')
    compTableFile.write(table3l4lComparisonString)
    compTableFile.close()

def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, returnplot = False):

    theCutDATA = cuts.AddList([theCut, cuts.trigger])
    MC = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx)
    MCS = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx)
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCutDATA, '', labelx)

    maxValmc = MC.GetMaximum() #GetBinContent(MC.GetMaximumBin())
    maxValdata = DATA.GetMaximum() #GetBinContent(DATA.GetMaximumBin())
    maxVal = max(maxValmc, maxValdata)
    if not logx:
        MC.GetYaxis().SetRangeUser(0, 1.3*maxVal)
        MCS.SetMaximum(1.3*maxVal)
        DATA.GetYaxis().SetRangeUser(0, 1.3*maxVal)
    else:
        MC.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
        MCS.SetMaximum(2.0*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
   
    print name
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('DataMC/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.7, 0.7, 0.95, 0.9)
    plot.addStack(MCS, "HIST", 1, 1)
    plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA, MC)
    
    if returnplot:
        return plot
    else:
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

    #dyDatasets = ['DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    dyDatasets = ['DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext', 'DYJetsToLL_M50_HT400to600_ext', 'DYJetsToLL_M50_HT600toInf_ext']
    mcDatasets = ['ZZ', 'TTJets_DiLepton_total', 'WWTo2L2Nu', 'WZTo2L2Q', 'WZTo3LNu', 'TTZToLLNuNu', 'TTWToLNu']
    mcDatasets += dyDatasets
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2', 'DoubleEG_Run2016B-PromptReco-v2', 'MuonEG_Run2016B-PromptReco-v2']



    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 2.6 ; maxrun = 999999
    lumi_str = '2.6invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'MET [GeV]'
    labelnjet = "N. Jets"
    
    ## ewino_SR = makePlot(0.8, '0.8invfb', treeDA, treeMC, "met_Edge", "met_ewino_SR", 12, 0, 300, cuts.AddList([cuts.ewinoSR, 'run_Edge <= 274240']), cuts, labelmet, 0, True)
    
    ## #3l and 4l plots
    ## plot_3l = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion_AF", 14, 60, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l]), cuts, labelmet, 0, True)
    ## plot_4l = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_4lregion_AF", 10,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l]), cuts, labelmll, 0, True)
    ## makeSummaryTable3l4l(plot_3l, plot_4l)

    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_met50_2jets_OF", 20,  0,  200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2, 'met_Edge > 50']), cuts, 'M_{T2}', 1)
    makePlot(lumi, lumi_str, treeDA, treeMC, "mt2_Edge", "mt2_met50_2jets_SF", 20,  0,  200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2, 'met_Edge > 50']), cuts, 'M_{T2}', 1)

    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_OF", 40,  0,  40, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, 'n_{vertices}', 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_SF", 40,  0,  40, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, 'n_{vertices}', 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_EE", 40,  0,  40, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, 'n_{vertices}', 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_MM", 40,  0,  40, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, 'n_{vertices}', 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_SF", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, labelmll, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_OF", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, labelmll, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_ee", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, labelmll, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_mm", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, labelmll, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.DYControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.JetMETBaseline]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.JetMETBaseline]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.JetMETBaseline]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.JetMETBaseline]), cuts, labelmll, 0)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj2]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj2]), cuts, labelmet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_SF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.SF]), cuts, labelnjet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_OF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.OF]), cuts, labelnjet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_ee", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.ee]), cuts, labelnjet, 1)
    ## makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_mm", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.mm]), cuts, labelnjet, 1)
