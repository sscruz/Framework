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
    

def makeSummaryTable3l4l(plot, plot_sig, up, dn, plot_up, plot_dn):

    sig, sig_e = 0., 0.
    sigUp, sigDn = 0., 0.
    bkg, bkg_e = 0., 0.
    bkgUp, bkgDn = 0., 0.
    obs= 0
    rat, rat_e = 0., 0.         


        
    for i, h in enumerate(plot.histos):
        if not i: continue
        tmp_double = r.Double()
        print "h.GetName() ", h.GetName()
        if 'ZZTo4l' in h.GetName():
            old = h.IntegralAndError(1, h.GetNbinsX()+1, tmp_double)
            old_e = tmp_double
        elif not 'DATA' in h.GetName(): 
            tmp_yield = h.IntegralAndError(1, h.GetNbinsX()+1, tmp_double)
            bkg  += tmp_yield
            bkg_e = math.sqrt(bkg_e**2 + tmp_double**2)
        elif 'DATA' in h.GetName():
            obs = h.Integral()                                                      
    
    for i, h in enumerate(plot_sig.histos):
        sig = h.IntegralAndError(1, 4, tmp_double)
        print "sig int", sig
        sig_e = tmp_double                            
    for i, h in enumerate(plot_up.histos):
        sigUp = h.IntegralAndError(1, 4, tmp_double)
        print "sigUp int", sigUp
        sigUp_e = tmp_double                            
    for i, h in enumerate(plot_dn.histos):
        sigDn = h.IntegralAndError(1, 4, tmp_double)
        print "sigDn int", sigDn
        sigDn_e = tmp_double                            
    
    
    for i, h in enumerate(up.histos):                           
        if not i: continue
        print "jec up h.GetName() ", h.GetName()
        if 'ZZTo4l' in h.GetName():
            sigUpOld = h.Integral()
        elif not 'DATA' in h.GetName(): 
            tmp_yieldUp = h.Integral()
            bkgUp  += tmp_yieldUp                         

    for i, h in enumerate(dn.histos):                      
        if not i: continue
        print "jec down h.GetName() ", h.GetName()
        if 'ZZTo4l' in h.GetName():
            sigDnOld = h.Integral()
        elif not 'DATA' in h.GetName(): 
            tmp_yieldDn = h.Integral()
            bkgDn  += tmp_yieldDn                               

    print "bkgUp-bkg", abs(bkgUp-bkg)
    print "bkg-bkgDn", abs(bkg-bkgDn)
    print "sigUp-sig", abs(sigUp-sig)
    print "sig-sigDn", abs(sig-sigDn)
    print "bkg_e ", bkg_e**2              

    if abs(bkgUp-bkg) > abs(bkg-bkgDn): bkgSyst_e = abs(bkgUp-bkg) 
    else:bkgSyst_e = bkg-bkgDn
    if abs(sigUp-sig) > abs(sig-sigDn): sigSyst_e = abs(sigUp-sig) 
    else:sigSyst_e = abs(sig-sigDn)
    obsSub,     obsSub_e = obs-bkg, math.sqrt(obs + bkg_e**2)
    obsSubSyst, obsSubSyst_e = obs-bkg, math.sqrt(bkgSyst_e**2)
    rat, rat_e = helper.ratioError(obsSub, obsSub_e, sig, sig_e)
    ratSyst, ratSyst_e = helper.ratioError(obsSub, obsSubSyst_e, sig, sigSyst_e)



#    obsSub, obsSub_e = obs-bkg4l, math.sqrt(bkg4l_e**2 + obs4l)
#    rat4l, rat4l_e = helper.ratioError(obsSub4l, obsSub4l_e,zz4l, zz4l_e)

    table3l4lComparisonString = '''\\documentclass[12pt,a4paper]{{article}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{Scale factor derived in a 4 lepton control region.}} 
\\label{{tab:tab4lcontrol}} 
\\begin{{tabular}}{{l| c }} 
              & 4-lepton region \\\\ \\hline \\hline
signal MC        & {sig:.2f}     $\\pm$  {sig_e:.2f}   $\\pm$  {sigSyst_e:.2f}  \\\\ \\hline
bkg. MC          & {bkg:.2f}  $\\pm$  {bkg_e:.2f}  $\\pm$  {bkgSyst_e:.2f}   \\\\ \\hline \\hline
\\textbf{{data}}       & \\textbf{{{obs}}}  \\\\ \\hline
data-bkg.        &  {obsSub:.2f}   $\\pm$  {obsSub_e:.2f} $\\pm$  {obsSubSyst_e:.2f}   \\\\ \\hline \\hline
(data-bkg.)/sig. & {rat:.2f}   $\\pm$  {rat_e:.2f}  $\\pm$  {ratSyst_e:.2f}   \\\\ \\hline

\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(sig=sig, sig_e=sig_e,sigSyst_e=sigSyst_e,
                            bkg=bkg,bkg_e=bkg_e,bkgSyst_e=bkgSyst_e,
                            obs=int(obs),
                            obsSub=obsSub, obsSub_e=obsSub_e, obsSubSyst_e=obsSubSyst_e,
                            rat=rat, rat_e=rat_e,  ratSyst_e=ratSyst_e)                                                          

    helper.ensurePath('plots/CRs/35.9invfb/tables/')
    compTableFile = open('plots/CRs/35.9invfb/tables/controlRegion4l.txt','w')
    compTableFile.write(table3l4lComparisonString)
    compTableFile.close()

def makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, treeMCForPlots, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, doKFactor):

    returnplot =True
    normalized =False
    theCutDATA = cuts.AddList([theCut, cuts.trigger ])
    theCutZ = cuts.AddList([theCut, cuts.Zmass ])
    MC4l   = tree4l.getTH1F(lumi, "hMC4l_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx, "1",  'ZZpt')
    MCBKG   = treeBKG.getTH1F(lumi, "hMCBKG_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx,  "1",  'ZZpt')
    MC   = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx, "1", 'ZZpt')
    MC_jecUp   = treeMC.getTH1F(lumi, "hMC_jecUp%s"%(name), "met_jecUp_Edge", nbin, xmin, xmax, theCut, '', labelx, "1", 'ZZpt')
    MC_jecDn   = treeMC.getTH1F(lumi, "hMC_jecDn%s"%(name), "met_jecDn_Edge", nbin, xmin, xmax, theCut, '', labelx, "1", 'ZZpt')
    MCS  = treeMCForPlots.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx)
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCutDATA, '', labelx, "1", 'ZZpt')
    MCS_jecUp   = treeMCForPlots.getStack(lumi, "hMCS_jecUp%s"%(name), "met_jecUp_Edge", nbin, xmin, xmax, theCut, "", labelx)
    MCS_jecDn   = treeMCForPlots.getStack(lumi, "hMCS_jecDn%s"%(name), "met_jecDn_Edge", nbin, xmin, xmax, theCut, "", labelx)


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
    print "zz4l %.2f +- %.2f "%(zz4l, zz4l_e) 
    print "da %.2f "%(obs4l) 
    print "da with bkg sub%.2f "%(obsSub4l) 
    rat4l, rat4l_e = helper.ratioError(obsSub4l, obsSub4l_e,zz4l, zz4l_e)
    print "rat %.2f +- %.2f "%(rat4l, rat4l_e)

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

    if normalized:
        MC.GetYaxis().SetRangeUser(0.01, 1.2)
   
    print name
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('CRs/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.65, 0.85, 0.9)
    if not normalized:
        plot.addStack(MCS , "HIST", 1, 1)
        plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    else:
        for h in hists:
            plot.addHisto(h, 'hist,same' if not h.GetName() == 'data' else 'p,same', h.GetName(), 'PL', h.GetLineColor(), 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA, MC)

    plot_jecUpKF = Canvas.Canvas('CRs/%s/plotupKF_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.58, 0.87, 0.9)
    plot_jecUpKF.addHisto(MC_jecUp , "HIST", "4lUp", "P", r.kBlack, 1, 0)                                                         
    plot_jecDnKF = Canvas.Canvas('CRs/%s/plotdnKF_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.58, 0.87, 0.9)
    plot_jecDnKF.addHisto(MC_jecDn , "HIST", "4lDn", "P", r.kBlack, 1, 0)                                                         
    plot_jecUp = Canvas.Canvas('CRs/%s/plotup_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.58, 0.87, 0.9)
    plot_jecUp.addStack(MCS_jecUp , "HIST", 1, 1)                                                         
    plot_jecDn = Canvas.Canvas('CRs/%s/plotdn_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.58, 0.87, 0.9)
    plot_jecDn.addStack(MCS_jecDn , "HIST", 1, 1)                                                         
    plot_sig = Canvas.Canvas('CRs/%s/plotsig_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.58, 0.87, 0.9)
    plot_sig.addHisto(MC4l , "HIST", "4l", "P", r.kBlack, 1, 0)                                                         
    
    #makeSummaryTable3l4l(plot, plot_sig, plot_jecUp, plot_jecDn, plot_jecUpKF, plot_jecDnKF )
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

    #mcDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    #ttDatasets = ['TTTo2L2Nu']
    #fourlDatasets = ['ZZTo4L']
    fourlDatasets = ['ZZTo4L','GluGluToContinToZZTo2e2mu', 'GluGluToContinToZZTo2e2tau', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau']
    #fourlDatasets = ['GGHZZ4L', 'ZZTo4L','GluGluToContinToZZTo2e2mu',  'GluGluToContinToZZTo2e2tau', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau', 'WWG', 'WWW', 'WWZ', 'WZZ', 'ZZZ', 'VHToNonbb']
   
    mcDatasets = [ 'GGHZZ4L', 'QQHZZ4L', 'ZZTo4L',  'ZGTo2LG', 'WGToLNuG','WZTo3LNu', 'WZTo2L2Q' , 'WWG',  'WWW', 'WWZ', 'WZZ', 'ZZZ', 'WGToLNuG', 'TTHnobb_pow','VHToNonbb', 'TWZ', 'tZq_ll', 'TBar_tch_powheg',  'T_tch_powheg',  'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTZToLLNuNu_ext2', 'TTWToLNu_ext2', 'GluGluToContinToZZTo2e2tau',  'GluGluToContinToZZTo2e2mu', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau']
    bkgDatasets = [ 'GGHZZ4L', 'QQHZZ4L',  'ZGTo2LG', 'WGToLNuG','WZTo3LNu', 'WZTo2L2Q' , 'WWG',  'WWW', 'WWZ', 'WZZ', 'ZZZ', 'WGToLNuG', 'TTHnobb_pow','VHToNonbb', 'TWZ', 'tZq_ll', 'TBar_tch_powheg',  'T_tch_powheg',  'TTTT',  'TTZToQQ', 'TTWToQQ','TTLLJets_m1to10', 'TTZToLLNuNu_ext2', 'TTWToLNu_ext2']
    #mcDatasets = [ 'GGHZZ4L', 'QQHZZ4L', 'ZZTo4L',  'ZGTo2LG', 'WGToLNuG','WZTo3LNu', 'WZTo2L2Q' , 'WWG',  'WWW', 'WWZ', 'WWTo2L2Nu',  'GGWWTo2L2Nu',  'GGHWWTo2L2Nu', 'WWDouble','WpWpJJ',   'WZZ', 'ZZZ', 'WGToLNuG', 'TTHnobb_pow','TTTo2L2Nu','TTJets_SingleLeptonFromT','TTJets_SingleLeptonFromTbar',  'DYJetsToLL_M50_LO', 'DYJetsToLL_M10to50_LO',  'VHToNonbb', 'TWZ', 'tZq_ll', 'WJetsToLNu_LO',  'TBar_tch_powheg',  'T_tch_powheg',  'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTZToLLNuNu_ext2', 'TTWToLNu_ext2', 'GluGluToContinToZZTo2e2tau',  'GluGluToContinToZZTo2e2mu', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau', 'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo2e2nu']
    #mcDatasets = ['GGHZZ4L', 'QQHZZ4L', 'ZZTo4L', 'WZTo3LNu','WWG',  'WWW', 'WWZ',  'WZZ', 'ZZZ', 'WGToLNuG', 'TTHnobb_pow', 'VHToNonbb', 'TWZ',   'T_tch_powheg',  'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTWToLNu_ext2', 'GluGluToContinToZZTo2e2tau',  'GluGluToContinToZZTo2e2mu', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau']
    #bkgDatasets = ['WZTo3LNu','WWG',  'WWW', 'WWZ',  'WZZ', 'ZZZ',  'WGToLNuG', 'TTHnobb_pow', 'VHToNonbb', 'TWZ',  'T_tch_powheg', 'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTWToLNu_ext2']
    
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


    tree4l  = Sample.Tree(helper.selectSamples("samplesSleptonFor4l.dat", fourlDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 0)
    treeMC  = Sample.Tree(helper.selectSamples("samplesSleptonFor4l.dat", mcDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 0)
    treeMCForPlots  = Sample.Tree(helper.selectSamples("samplesSleptonFor4lForPlots.dat", mcDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 0)
    treeBKG = Sample.Tree(helper.selectSamples("samplesSleptonFor4l.dat", bkgDatasets, 'MC'), 'MC'  , 0, isScan = 0, isOnEOS = 0)
    treeDA  = Sample.Tree(helper.selectSamples("samplesSleptonFor4l.dat", daDatasets, 'DA'), 'DATA', 1, isScan = 0, isOnEOS = 0)
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
    #makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OFSF_nll", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLine]), cuts, labelmll, 0)
    
    ## EWino signal region and table
    ## =================================================================
    #ewino_SR = makePlot(12.9, '12.9invfb', treeDA, treeMC, "met_Edge", "met_ewino_SR", 12, 0, 300, cuts.AddList([cuts.ewinoSR, 'run_Edge <= 999999']), cuts, labelmet, 0, True)
    #makeSummaryEWino(ewino_SR)
    

    #makePlot(lumi, lumi_str, treeSI, treeMC, "bestMjj_Edge", "mjj_ewinoSR", 20,  0,  300, cuts.AddList([cuts.ewinoSR]), cuts, 'm_{jj}', 0, True)
    #makePlot(lumi, lumi_str, treeSI, treeMC, "bestMjj_Edge", "mjj", 20,  0,  300, cuts.AddList([cuts.goodLepton]), cuts, 'm_{jj}', 0, True)
    ## 3l and 4l plots
    ## =================================================================
#    plot_3l = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion", 13, 70, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.region3lSlepton]), cuts, labelmet, 0, True)
#    plot_3l_Zveto = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion", 13, 70, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.region3lSlepton, "mZ1_Edge > 96 && mZ1_Edge < 86"]), cuts, labelmet, 0, True)
#    plot_3l_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion_Zmass", 13, 70, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.region3lSlepton, "mZ1_Edge < 96 && mZ1_Edge > 86"]), cuts, labelmet, 0, True)
#    ##plot_3l_njets = makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge"    , "met_3lregion_njets",9, 0.5, 9.5, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l, cuts.Zmass]), cuts, labelmet, 0, True)
    #plot_2l2nu_DPhiLep2Met = makePlot(lumi, lumi_str, treeMC, treeMC, "metl2DPhi_Edge", "dPhiLep2Met_2l2nu", 40,  0, 10, cuts.AddList([cuts.goodLepton, cuts.SF,  cuts.region2lnuSlepton]), cuts, "dPhi(l2, MET)", 0, True)
    #plot_2l2nu_DPhiLep1Met = makePlot(lumi, lumi_str, treeMC, treeMC, "metl1DPhi_Edge", "dPhiLep1Met_2l2nu", 40,  0, 10, cuts.AddList([cuts.goodLepton, cuts.SF,  cuts.region2lnuSlepton]), cuts, "dPhi(l1, MET)", 0, True)
#    plot_2l2nu_Mll = makePlot(lumi, lumi_str, treeMC, treeMC, "mZ1_Edge", "bestMll_2l2nu", 40,  0, 400, cuts.AddList([cuts.goodLepton, cuts.SF,  cuts.region2lnuSlepton]), cuts, "mll of best Z candidate", 0, True)
#    plot_2l2nu_met = makePlot(lumi, lumi_str, treeMC, treeMC, "met_Edge", "met_2l2nu", 40,  0, 400, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.region2lnuSlepton]), cuts, "met", 0, True)
#    plot_2l2nu_lepsDPhi = makePlot(lumi, lumi_str, treeMC, treeMC, "lepsDPhi_Edge", "lepsDPhi_2l2nu", 40,  0, 10, cuts.AddList([cuts.goodLepton, cuts.SF,  cuts.region2lnuSlepton]), cuts, "dPhi(l1, l2)", 0, True)
#    plot_4l_newMet = makePlot(lumi, lumi_str, tree4l, tree4l, "newMet_Edge", "newMet_4l", 40,  0, 400, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, cuts.mZ1inZ]), cuts, "new met", 0, True)
#    plot_4l_newMetPhi = makePlot(lumi, lumi_str, tree4l, tree4l, "newMetPhi_Edge", "newMetPhi_4l", 40,  -5, 5, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, cuts.mZ1inZ]), cuts, "new met phi", 0, True)
#    plot_2l2nu_met  = makePlot(lumi, lumi_str, tree2l, tree2l, "met_Edge", "met_2l2nu", 40,  0, 400, cuts.AddList([cuts.goodLepton, cuts.region2l2nuSlepton]), cuts, "met 2l2nu", 0, True)
#    plot_2l2nu_metPhi  = makePlot(lumi, lumi_str, tree2l, tree2l, "met_phi_Edge", "metPhi_2l2nu", 40,  -5, 5, cuts.AddList([cuts.goodLepton, cuts.region2l2nuSlepton]), cuts, "met phi 2l2nu", 0, True)
    #plot_4l_met  = makePlot(lumi, lumi_str, treeDA, tree4l, "newMet_Edge", "newMet_4l", 15,  0, 300, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, " mZ1_Edge < 106 && mZ1_Edge > 76 && mZ2_Edge < 106 && mZ2_Edge > 76"]), cuts, "new E_{T}^{miss} [GeV]", 0, True)
    #plot_4l_njets  = makePlot(lumi, lumi_str, treeDA, tree4l, "nJet25_Edge", "njets25_4l", 7,  0, 7, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton]), cuts, "njets (pt > 25)", 0, True)
    #plot_4l_njets35  = makePlot(lumi, lumi_str, treeDA, tree4l, "nJet35_Edge", "njets35_4l", 7,  0, 7, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton]), cuts, "njets (pt > 35)", 0, True)
    #plot_4l_njets  = makePlot(lumi, lumi_str, treeDA, tree4l, "nJet25_Edge", "njets25_4l", 7,  0, 7, cuts.AddList([cuts.goodLepton4l, cuts.region4lSlepton, "mZ2_Edge < 101 && mZ2_Edge > 81 && nPairLep_Edge > 0"]), cuts, "njets (pt > 25)", 0, True)
    #plot_4l_njets35  = makePlot(lumi, lumi_str, treeDA, tree4l, "nJet35_Edge", "njets35_4l", 7,  0, 7, cuts.AddList([cuts.goodLepton4l, cuts.region4lSlepton, "mZ2_Edge < 101 && mZ2_Edge > 81 && nPairLep_Edge > 0"]), cuts, "njets (pt > 35)", 0, True)
    #plot_4l_jet1Pt = makePlot(lumi, lumi_str, treeDA, tree4l, "JetSel_Edge_pt[0]", "jet1Pt_4l", 100,  0, 200, cuts.AddList([cuts.goodLepton4l, cuts.region4lSlepton, " mZ1_Edge < 101 && mZ1_Edge > 81 && mZ2_Edge < 101 && mZ2_Edge > 81"]), cuts, "Leading Jet p_{T} [GeV]", 0, True)
    #plot_4l_jet2Pt = makePlot(lumi, lumi_str, treeDA, tree4l, "JetSel_Edge_pt[1]", "jet2Pt_4l", 100,  0, 200, cuts.AddList([cuts.goodLepton4l, cuts.region4lSlepton, " mZ1_Edge < 101 && mZ1_Edge > 81 && mZ2_Edge < 101 && mZ2_Edge > 81"]), cuts, "Subleading Jet p_{T} [GeV]", 0, True)
    #plot_4l_jet3Pt = makePlot(lumi, lumi_str, treeDA, tree4l, "JetSel_Edge_pt[2]", "jet3Pt_4l", 100,  0, 200, cuts.AddList([cuts.goodLepton4l, cuts.region4lSlepton, " mZ1_Edge < 101 && mZ1_Edge > 81 && mZ2_Edge < 101 && mZ2_Edge > 81"]), cuts, "Third Jet Jet p_{T} [GeV]", 0, True)
    plot_4l_bestMll  = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG,treeMCForPlots, "mZ1_Edge", "bestMll_4l", 100,  70, 110, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, "mZ1_Edge < 106 && mZ1_Edge > 76 &&mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "m_{ll} of best Z candidate [GeV]", 0, True)
# start here
#    plot_4l_newMET = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, "newMet_Edge", "newMET_4l_Kfactor", 16,  0, 400, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, " mZ1_Edge < 106 && mZ1_Edge > 76 && mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "New p_{T}^{miss} [GeV]", 0, True, True)
#    plot_4l_newMT2 = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, "mt2BestZ_Edge", "newMT2_4l_Kfactor", 16,  0, 400, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, " mZ1_Edge < 106 && mZ1_Edge > 76 && mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "New M_{T2} [GeV]", 0, True, True)
    plot_4l_otherMll = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, treeMCForPlots,"mZ2_Edge", "otherMll_4l_Kfactor", 39,  50, 128, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, " mZ1_Edge < 106 && mZ1_Edge > 76 && mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "m_{ll} of other Z candidate [GeV]", 0, True)
    #plot_4l_met_noKF  = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, treeMCForPlots, "met_Edge", "met_4l_noKfactor", 16,  0, 400, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, "mZ1_Edge < 106 && mZ1_Edge > 76 &&  mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "m_{ll} of best Z candidate [GeV]", 0, 'noKFactor')
    #plot_4l_bestMll_noKF  = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, treeMCForPlots, "mZ1_Edge", "bestMll_4l_noKfactor", 35,  75, 110, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, "mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "m_{ll} of best Z candidate [GeV]", 0, 'noKFactor')
    #makeSummaryTable3l4l(plot_4l_bestMll_noKF, 'noKFactor')
#    plot_4l_bestMll_ZZpt  = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, "mZ1_Edge", "bestMll_4l_ptKfactor", 35,  75, 110, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, "mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "m_{ll} of best Z candidate [GeV]", 0, 'ZZpt')
#    makeSummaryTable3l4l(plot_4l_bestMll_ZZpt, 'ZZpt')
#    plot_4l_bestMll_ZZmass  = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, "mZ1_Edge", "bestMll_4l_massKfactor", 35,  75, 110, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, "mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "m_{ll} of best Z candidate [GeV]", 0, 'ZZmass')
#    makeSummaryTable3l4l(plot_4l_bestMll_ZZmass, 'ZZmass')
#    plot_4l_bestMll_ZZdPhi  = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, "mZ1_Edge", "bestMll_4l_dPhiKfactor", 35,  75, 110, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, "mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "m_{ll} of best Z candidate [GeV]", 0, 'ZZdPhi')
#    makeSummaryTable3l4l(plot_4l_bestMll_ZZdPhi, 'ZZdPhi')
    #plot_4l_bestMT2 = makePlot(lumi, lumi_str, treeDA, tree4l, "mt2BestZ_Edge", "bestMT2_4l", 15,  0, 300, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton , " mZ1_Edge < 106 && mZ1_Edge > 76 && mZ2_Edge < 106 && mZ2_Edge > 76"]), cuts, "new MT2", 0, True)
    #plot_4l_otherMll = makePlot(lumi, lumi_str, treeDA, tree4l, "mZ2_Edge", "otherMll_4l", 100,  0, 200, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton]), cuts, "m_{ll} of other Z candidate [GeV]", 0, True)
    #plot_4l_otherMll = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG, "mZ2_Edge", "otherMll_4l", 35,  60, 120, cuts.AddList([cuts.goodLepton4l, cuts.region4lSlepton, " mZ1_Edge < 106 && mZ1_Edge > 76 && mZ2_Edge < 120 && mZ2_Edge > 60"]), cuts, "m_{ll} of other Z candidate [GeV]", 0, True)
    plot_4l_lep1_pt = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG,treeMCForPlots, "Lep1_pt_Edge", "Lep1_pt_4l", 35,  0, 200, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, " mZ1_Edge < 106 && mZ1_Edge > 76 && mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "Leading lepton p_{T} [GeV]", 0, True)
    plot_4l_lep2_pt = makePlot(lumi, lumi_str, treeDA, treeMC, tree4l, treeBKG,treeMCForPlots, "Lep2_pt_Edge", "Lep2_pt_4l", 35,  0, 200, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, " mZ1_Edge < 106 && mZ1_Edge > 76 && mZ2_Edge < 130 && mZ2_Edge > 50"]), cuts, "Subleading lepton p_{T} [GeV]", 0, True)
    #plot_4l_jet1Pt = makePlot(lumi, lumi_str, treeDA, tree4l, "JetSel_Edge_pt[0]", "jet1Pt_4l", 100,  0, 200, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, " mZ1_Edge < 101 && mZ1_Edge > 81 && mZ2_Edge < 101 && mZ2_Edge > 81"]), cuts, "Leading Jet p_{T} [GeV]", 0, True)
    #plot_4l_jet2Pt = makePlot(lumi, lumi_str, treeDA, tree4l, "JetSel_Edge_pt[1]", "jet2Pt_4l", 100,  0, 200, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, " mZ1_Edge < 101 && mZ1_Edge > 81 && mZ2_Edge < 101 && mZ2_Edge > 81"]), cuts, "Subleading Jet p_{T} [GeV]", 0, True)
#    plot_2l2nu_met  = makePlot(lumi, lumi_str, treeDA, tree4l, "met_Edge", "met_2l2nu", 40,  0, 400, cuts.AddList([cuts.goodLepton, cuts.region2l2nuSlepton]), cuts, "met 2l2nu", 0, True)
#    plot_4l_met_bestMll   = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_4l_bestMll", 40,  0, 400, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, cuts.mZ1inZ]), cuts, "met ", 0, True)
#    plot_4l_met_otherMll   = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_4l_otherMll", 40,  0, 400, cuts.AddList([cuts.goodLepton, cuts.region4lSlepton, cuts.mZ2inZ]), cuts, "met ", 0, True)
