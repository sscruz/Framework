#
###################################################################
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
import include.LeptonSF
import include.FastSimSF


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

def compareMMEEPlots(plotmm, plotee):
    for _i,h in enumerate(plotmm.histos):
        if 'hDATA' in h.GetName():
            ind = _i
    h_mm = plotmm.histos[_i]
    h_ee = plotee.histos[_i]
    name = plotmm.name.split('/')[-1].replace('plot_','').replace('mm_','')
    lumistr = plotmm.name.split('/')[-2]

    plot = Canvas.Canvas('DataMC/%s/plot_mmeeComparison_%s'%(lumistr,name), 'png,pdf,root', 0.65, 0.7, 0.85, 0.9)
    plot.addHisto(h_ee, "PE1,SAME", "elel", "P", r.kYellow+1, 1, 1)
    plot.addHisto(h_mm, "PE1,SAME", "mumu", "P", r.kBlack   , 1, 0)
    plot.saveRatio(1, 1, 0, 4, h_ee, h_mm)
    

def makeSummaryTable3l4l(plot):
    sig, sig_e = 0., 0.
    bkg, bkg_e = 0., 0.
    obs= 0
    rat, rat_e = 0., 0.         
    
        
    for i, h in enumerate(plot.histos):
        if not i: continue
        tmp_double = r.Double()
        print "h.GetName() ", h.GetName()
        if 'WZTo3l' in h.GetName():
            sig = h.IntegralAndError(1, h.GetNbinsX()+1, tmp_double)
            sig_e = tmp_double
        elif not 'DATA' in h.GetName(): 
            tmp_yield = h.IntegralAndError(1, h.GetNbinsX()+1, tmp_double)
            bkg  += tmp_yield
            bkg_e = math.sqrt(bkg_e**2 + tmp_double**2)
        elif 'DATA' in h.GetName():
            obs = h.Integral()                                                      
   
    
    obsSub,     obsSub_e = obs-bkg, math.sqrt(obs + bkg_e**2)
    rat, rat_e = helper.ratioError(obsSub, obsSub_e, sig, sig_e)

    table3l4lComparisonString = '''\\documentclass[12pt,a4paper]{{article}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{Scale factor derived in a 3 lepton control region.}} 
\\label{{tab:tab3lcontrol}} 
\\begin{{tabular}}{{l| c }} 
              & 3-lepton region \\\\ \\hline \\hline
signal MC        & {sig:.2f}     $\\pm$  {sig_e:.2f} \\\\ \\hline
bkg. MC          & {bkg:.2f}  $\\pm$  {bkg_e:.2f}\\\\ \\hline \\hline
\\textbf{{data}}       & \\textbf{{{obs}}}  \\\\ \\hline
data-bkg.        &  {obsSub:.2f}   $\\pm$  {obsSub_e:.2f} \\\\ \\hline \\hline
(data-bkg.)/sig. & {rat:.2f}   $\\pm$  {rat_e:.2f}\\\\ \\hline

\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(sig=sig, sig_e=sig_e,
                            bkg=bkg,bkg_e=bkg_e,
                            obs=int(obs),
                            obsSub=obsSub, obsSub_e=obsSub_e,
                            rat=rat, rat_e=rat_e)                                                          
    
    helper.ensurePath('plots/CRs/35.9invfb/tables/')
    compTableFile = open('plots/CRs/35.9invfb/tables/controlRegion3l_comparison.txt','w')
    compTableFile.write(table3l4lComparisonString)
    compTableFile.close()

def makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, treeMCForPlots, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, returnplot = False, scaleToData = False, normalized = False, cumulative = False, onlyMC = False):

    print 'returnplot', returnplot
    print 'scaletodata', scaleToData
    print 'normalized', normalized
    print 'cumulative', cumulative
    print 'onlyMC', onlyMC

    theCutDATA = cuts.AddList([theCut, cuts.trigger, cuts.nj25eq0])
    theCut = cuts.AddList([theCut, cuts.nj25eq0])
    MC4l        = tree4l.getTH1F(lumi, "hMC4l_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx, "1", 'ZZpt')
    MCBKG       = treeBKG.getTH1F(lumi, "hMCBKG_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx, "1", 'ZZpt')
    MC          = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx, "1", 'ZZpt')
    MCS         = treeMCForPlots.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx)
    DATA        = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCutDATA, '', labelx, "1", 'ZZpt')
    
    print "MCBKG ", MCBKG.Integral() 
    print "MC ", MC.Integral() 
    for _i,_h in enumerate(MCS.GetHists()):
        errorInte =  r.Double()
        inte = _h.IntegralAndError(1, _h.GetNbinsX()+1, errorInte)
        print "sample ", _h.GetName(), " has integral ", _h.Integral(), " and error ", errorInte
    zz4l_e = r.Double()
    bkg4l_e = r.Double()
    zz4l = MC4l.IntegralAndError(1, MC4l.GetNbinsX()+1, zz4l_e)
    bkg4l = MCBKG.IntegralAndError(1, MCBKG.GetNbinsX()+1, bkg4l_e)
    #zz4l_e = tmp_double
    obs4l = DATA.Integral()                                                     
    obsSub4l, obsSub4l_e = obs4l-bkg4l, math.sqrt(bkg4l_e**2 + obs4l)
    print "wz3l %.2f +- %.2f "%(zz4l, zz4l_e) 
    print "da %.2f "%(obs4l) 
    print "bkg %.2f "%(bkg4l) 
    print "da with bkg sub%.2f "%(obsSub4l) 
    rat4l, rat4l_e = helper.ratioError(obsSub4l, obsSub4l_e,zz4l, zz4l_e)
    print "rat %.2f +- %.2f "%(rat4l, rat4l_e)

    maxValmc = MC.GetMaximum() #GetBinContent(MC.GetMaximumBin())
    maxValdata = DATA.GetMaximum() #GetBinContent(DATA.GetMaximumBin())
    maxVal = max(maxValmc, maxValdata)
    if not logx:
        MC.GetYaxis().SetRangeUser(0.1, 1.3*maxVal)
        MCS.SetMaximum(1.3*maxVal)
        if scaleToData:
            newStack.SetMaximum(1.3*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 1.3*maxVal)
    else:
        MC.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
        MCS.SetMaximum(2.0*maxVal)
        if scaleToData:
            newStack.SetMaximum(2.0*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)

    if normalized:
        MC.GetYaxis().SetRangeUser(0.01, 1.2)
   
    print name
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('CRs/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.58, 0.87, 0.9)
    plot.addStack(MCS if not scaleToData else newStack, "HIST", 1, 1)
    plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA, MC)                                                                                                         

    

    makeSummaryTable3l4l(plot)
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

    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    threelDatasets = ['WZTo3LNu']
   
    #mcDatasets = [ 'GGHZZ4L', 'QQHZZ4L', 'ZZTo4L',  'ZGTo2LG', 'WGToLNuG','WZTo3LNu', 'WWG',  'WWW', 'WWZ', 'WWTo2L2Nu',  'GGWWTo2L2Nu',  'GGHWWTo2L2Nu', 'WWDouble','WpWpJJ',   'WZZ', 'ZZZ', 'WGToLNuG', 'TTHnobb_pow','DYJetsToLL_M50_LO', 'DYJetsToLL_M10to50_LO',  'VHToNonbb', 'TWZ', 'tZq_ll', 'TBar_tch_powheg',  'T_tch_powheg',  'TTTT', 'TTLLJets_m1to10', 'TTZToLLNuNu_ext2', 'TTWToLNu_ext2', 'GluGluToContinToZZTo2e2tau',  'GluGluToContinToZZTo2e2mu', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau', 'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo2e2nu']
    mcDatasets = [ 'GGHZZ4L', 'QQHZZ4L', 'ZZTo4L',  'ZGTo2LG', 'WGToLNuG','WZTo3LNu', 'WZTo2L2Q' , 'WWG',  'WWW', 'WWZ', 'WWTo2L2Nu',  'GGWWTo2L2Nu',  'GGHWWTo2L2Nu', 'WWDouble','WpWpJJ',   'WZZ', 'ZZZ', 'WGToLNuG', 'TTHnobb_pow','TTTo2L2Nu', 'TTJets_SingleLeptonFromT', 'TTJets_SingleLeptonFromTbar',  'DYJetsToLL_M50_LO', 'DYJetsToLL_M10to50_LO',  'VHToNonbb', 'TWZ', 'tZq_ll', 'TBar_tch_powheg',  'T_tch_powheg',  'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTZToLLNuNu_ext2', 'TTWToLNu_ext2', 'GluGluToContinToZZTo2e2tau',  'GluGluToContinToZZTo2e2mu', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau', 'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo2e2nu']
    #bkgDatasets = [ 'GGHZZ4L', 'QQHZZ4L', 'ZZTo4L',  'ZGTo2LG', 'WGToLNuG', 'WWG',  'WWW', 'WWZ', 'WWTo2L2Nu',  'GGWWTo2L2Nu',  'GGHWWTo2L2Nu', 'WWDouble','WpWpJJ',   'WZZ', 'ZZZ', 'WGToLNuG', 'TTHnobb_pow','DYJetsToLL_M50_LO', 'DYJetsToLL_M10to50_LO',  'VHToNonbb', 'TWZ', 'tZq_ll', 'TBar_tch_powheg',  'T_tch_powheg',  'TTTT', 'TTLLJets_m1to10', 'TTZToLLNuNu_ext2', 'TTWToLNu_ext2', 'GluGluToContinToZZTo2e2tau',  'GluGluToContinToZZTo2e2mu', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau', 'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo2e2nu']
    bkgDatasets = [ 'GGHZZ4L', 'QQHZZ4L', 'ZZTo4L',  'ZGTo2LG', 'WGToLNuG','WZTo2L2Q' , 'WWG',  'WWW', 'WWZ', 'WWTo2L2Nu',  'GGWWTo2L2Nu',  'GGHWWTo2L2Nu', 'WWDouble','WpWpJJ',   'WZZ', 'ZZZ', 'WGToLNuG', 'TTHnobb_pow','TTTo2L2Nu', 'TTJets_SingleLeptonFromT', 'TTJets_SingleLeptonFromTbar',  'DYJetsToLL_M50_LO', 'DYJetsToLL_M10to50_LO',  'VHToNonbb', 'TWZ', 'tZq_ll', 'TBar_tch_powheg',  'T_tch_powheg',  'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTZToLLNuNu_ext2', 'TTWToLNu_ext2', 'GluGluToContinToZZTo2e2tau',  'GluGluToContinToZZTo2e2mu', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau', 'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo2e2nu']
    
                                                                              
    daDatasetsB = ['DoubleEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'DoubleMuon_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'MuonEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376']    
                                                                              
    daDatasetsC = ['DoubleEG_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016C_03Feb2017_v1_runs_271036_284044']    
    
    daDatasetsD = ['DoubleEG_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016D_03Feb2017_v1_runs_271036_284044']    
                                                                              
    daDatasetsE = ['DoubleEG_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016E_03Feb2017_v1_runs_271036_284044']    
                                                                              
    daDatasetsF = ['DoubleEG_Run2016F_03Feb2017_v1_runs_271036_284044',
                  'DoubleMuon_Run2016F_03Feb2017_v1_runs_271036_284044',
                  'MuonEG_Run2016F_03Feb2017_v1_runs_271036_284044']  
                                                                              
    daDatasetsG = ['DoubleEG_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016G_03Feb2017_v1_runs_271036_284044']    
                                                                              
    daDatasetsH = ['DoubleEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'DoubleMuon_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleMuon_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'MuonEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035', 
                   'MuonEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044']    

                                                                                                                   
                                                                                                                   
    daDatasets = daDatasetsB + daDatasetsC + daDatasetsD +daDatasetsE + daDatasetsF + daDatasetsG + daDatasetsH    


    tree3l  = Sample.Tree(helper.selectSamples("samplesSleptonFor3l.dat", threelDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 0)
    treeMC  = Sample.Tree(helper.selectSamples("samplesSleptonFor3l.dat", mcDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 0)
    treeMCForPlots  = Sample.Tree(helper.selectSamples("samplesSleptonFor3lForPlots.dat", mcDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 0)
    treeBKG = Sample.Tree(helper.selectSamples("samplesSleptonFor3l.dat", bkgDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 0)
    treeDA  = Sample.Tree(helper.selectSamples("samplesSleptonFor3lForPlots.dat", daDatasets, 'DA'), 'DATA', 1, isScan = 0, isOnEOS = 0)
    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 35.9 ; maxrun = 999999; lumi_str = '35.9invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    cuts2l = CutManager.CutManager()
    cuts4l = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelnjet = "N. Jets"
    labelht = "H_{T} [GeV]"
    ## =================================================================
    

    ## 3l and 4l plots
    ## =================================================================
    #plot_3l = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion", 13, 70, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.region3lSlepton]), cuts, labelmet, 0, True)
    #plot_3l_leadingLepsMll_ZmassExcl = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "lepsMll_Edge", "leadingLepsMll_3l_ZmassExcl", 100,  0, 300, cuts.AddList([cuts.goodLepton3l, cuts.region3lSlepton, "(lepsMll_Edge > 101 || (lepsMll_Edge < 81 && lepsMll_Edge > 0))"]), cuts, "m_{ll} leading leptons [GeV]", 0, True)
    #plot_3l_mll = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "mZ1_Edge", "mll_3l", 100,  0, 300, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton]), cuts, "m_{ll} of best Z candidate [GeV]", 0, True)
   # plot_3l_WmT_Zmass_met65 = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "WmT_Edge", "WmT_3l_Zmass_met65", 50,  0, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton, "mZ1_Edge < 101     && mZ1_Edge > 81 && met_Edge > 65"]), cuts, "M_{T} (w/ W lep) [GeV]", 0, True)
    #plot_3l_met_ZmassCut_WmT = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "met_Edge", "met_3l_ZmassCut_WmT40", 30,  50, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton, "mZ1_Edge < 101 && mZ1_Edge > 81 && WmT_Edge > 50"]), cuts, "p_{T}^{miss} [GeV]", 0, True)
    #plot_3l_mll_ZmassCut = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "mZ1_Edge", "mll_3l_ZmassCut", 50,  80, 102, cuts.AddList([cuts.goodLepton3l, cuts.region3lSlepton, "mZ1_Edge <  101 && mZ1_Edge > 81"]), cuts, "m_{ll} of best Z candidate [GeV]", 0, True)
    plot_3l_met_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, treeMCForPlots, "met_Edge", "met_3l_Zmass", 13,  70, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton, "mZ1_Edge < 106     && mZ1_Edge > 76"]), cuts, "p_{T}^{miss} [GeV]", 0, True)
    #makeSummaryTable3l4l(plot_3l_met_Zmass)
    #plot_3l_WmT_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, treeMCForPlots, "WmT_Edge", "WmT_3l_Zmass", 15,  50, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton, "mZ1_Edge < 106   && mZ1_Edge > 76"]), cuts, "M_{T} (w/ W lep) [GeV]", 0, True)
    #plot_3l_WMT2_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, treeMCForPlots, "WZ_ll_MT2_Edge", "WZMT2_3l_Zmass", 20,  0, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton, "mZ1_Edge < 106     && mZ1_Edge > 76"]), cuts, "M_{T2} (w/ W lep) [GeV]", 0, True)
    #plot_3l_ZMT2_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, treeMCForPlots, "Z_ll_MT2_Edge", "ZMT2_3l_Zmass", 19,  10, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton, "mZ1_Edge < 106     && mZ1_Edge > 76"]), cuts, "M_{T2} (w/ Z lep) [GeV]", 0, True)
    #plot_3l_mll = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG,treeMCForPlots,  "mZ1_Edge", "mll_3l", 25,  70, 120, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton]), cuts, "m_{ll} of best Z candidate [GeV]", 0, True)
    #makeSummaryTable3l4l(plot_3l_mll_ZmassCut)
    #plot_3l_WmT = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "WmT_Edge", "WmT_3l", 30,  50, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton, "mZ1_Edge < 101     && mZ1_Edge > 81"]), cuts, "M_{T} (w/ W lep) [GeV]", 0, True)
    #plot_3l_WmT = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "WmT_Edge", "WmT_3l", 50,  0, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton]), cuts, "M_{T} (w/ W lep) [GeV]", 0, True)
    #plot_3l_minMT = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "minMT_Edge", "minMT_3l", 50,  0, 300, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton]), cuts, "min M_{T} [GeV]", 0, True)
    #plot_3l_met_ZmassCut = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "met_Edge", "met_3l_ZmassCut", 30,  65, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton, "mZ1_Edge < 101 && mZ1_Edge > 81"]), cuts, "p_{T}^{miss} [GeV]", 0, True)
    #plot_3l_met = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "met_Edge", "met_3l", 30,  50, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton]), cuts, "p_{T}^{miss} [GeV]", 0, True)
    #plot_3l_mll = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "mZ1_Edge", "mll_3l", 50,  0, 200, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton]), cuts, "m_{ll} of best Z candidate [GeV]", 0, True)
    #plot_3l_mll_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, tree3l, treeBKG, "mZ1_Edge", "mll_3l_ZmassCut", 80,  22, 102, cuts.AddList([cuts.goodLepton, cuts.region3lSlepton]), cuts, "m_{ll} of best Z candidate [GeV]", 0, True)
    #makeSummaryTable3l4l(plot_3l_mll_ZmassCut)
#    plot_3l_Zveto = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion", 13, 70, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.region3lSlepton, "mZ1_Edge > 96 && mZ1_Edge < 86"]), cuts, labelmet, 0, True)
#    plot_3l_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion_Zmass", 13, 70, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.region3lSlepton, "mZ1_Edge < 96 && mZ1_Edge > 86"]), cuts, labelmet, 0, True)
#    ##plot_3l_njets = makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge"    , "met_3lregion_njets",9, 0.5, 9.5, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l, cuts.Zmass]), cuts, labelmet, 0, True)
    #makeSummaryTable3l4l(plot_3l_met_ZmassCut)
