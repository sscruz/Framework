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
    

def makeSummaryTable3l4l(plot_3l, plot_4l, plot_ttZ):
    wz3l, wz3l_e, zz4l, zz4l_e, ttZ, ttZ_e = 0., 0., 0., 0.,0., 0.
    bkg3l, bkg3l_e, bkg4l, bkg4l_e, bkgttZ, bkgttZ_e = 0., 0., 0., 0.,0., 0.
    obs3l, obs4l, obsttZ = 0, 0, 0
    rat3l, rat3l_e, rat4l, rat4l_e, ratttZ, ratttZ_e = 0., 0., 0., 0.,0., 0.
    
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
        print "h_4l.GetName() ", h_4l.GetName()
        if 'ZZ' in h_4l.GetName():
            zz4l = h_4l.IntegralAndError(1, h_4l.GetNbinsX()+1, tmp_double)
            zz4l_e = tmp_double
        elif not 'DATA' in h_4l.GetName(): 
            tmp_yield = h_4l.IntegralAndError(1, h_4l.GetNbinsX()+1, tmp_double)
            bkg4l  += tmp_yield
            bkg4l_e = math.sqrt(bkg4l_e**2 + tmp_double**2)
        elif 'DATA' in h_4l.GetName():
            obs4l = h_4l.Integral()                                                      
    
    for i, h_ttZ in enumerate(plot_ttZ.histos):
        if not i: continue
        tmp_double = r.Double()
        print "h_ttZ.GetName() ", h_ttZ.GetName()
        if 'TTZ' in h_ttZ.GetName():
            ttZ = h_ttZ.IntegralAndError(1, h_ttZ.GetNbinsX()+1, tmp_double)
            ttZ_e = tmp_double
        elif not 'DATA' in h_ttZ.GetName(): 
            print "now here!!! "
            tmp_yield = h_ttZ.IntegralAndError(1, h_ttZ.GetNbinsX()+1, tmp_double)
            bkgttZ  += tmp_yield
            bkgttZ_e = math.sqrt(bkgttZ_e**2 + tmp_double**2)
        elif 'DATA' in h_ttZ.GetName():
            print "here obs!!!"
            obsttZ = h_ttZ.Integral()     


    obsSub3l, obsSub3l_e, obsSub4l, obsSub4l_e = obs3l-bkg3l, math.sqrt(bkg3l_e**2 + obs3l), obs4l-bkg4l, math.sqrt(bkg4l_e**2 + obs4l)
    obsSubttZ, obsSubttZ_e = obsttZ-bkgttZ, math.sqrt(bkgttZ_e**2 + obsttZ)
    rat3l, rat3l_e = helper.ratioError(obsSub3l, obsSub3l_e, wz3l, wz3l_e)
    rat4l, rat4l_e = helper.ratioError(obsSub4l, obsSub4l_e,zz4l, zz4l_e)
    ratttZ, ratttZ_e = helper.ratioError(obsSubttZ, obsSubttZ_e,ttZ, ttZ_e)

    table3l4lComparisonString = '''\\documentclass[12pt,a4paper]{{article}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{3-lepton, 4-lepton and ttZ control regions. Signal MC is WZ-$>$3lnu in the 3-lepton region, ZZ-$>$4l in the 4-lepton region and ttZ-$>$2l2nu and ttZ-$>$qq in the ttZ region}} 
\\label{{tab:tab3l4lcontrol}} 
\\begin{{tabular}}{{l| c c c}} 
              & 3-lepton region & 4-lepton region &  ttZ region\\\\ \\hline \\hline
signal MC        & {wz3l:.2f}     $\\pm$  {wz3l_e:.2f}       & {zz4l:.2f}     $\\pm$  {zz4l_e:.2f}    & {ttZ:.2f}     $\\pm$  {ttZ_e:.2f}     \\\\ \\hline
bkg. MC          & {bkg3l:.2f}  $\\pm$  {bkg3l_e:.2f}        & {bkg4l:.2f}  $\\pm$  {bkg4l_e:.2f}    &   {bkgttZ:.2f}  $\\pm$  {bkgttZ_e:.2f}\\\\ \\hline \\hline
\\textbf{{data}}             & \\textbf{{{obs3l}}}                               & \\textbf{{{obs4l}}}    & \\textbf{{{obsttZ}}}   \\\\ \\hline
data-bkg.        & {obsSub3l:.2f}   $\\pm$  {obsSub3l_e:.2f} & {obsSub4l:.2f}   $\\pm$  {obsSub4l_e:.2f}   & {obsSubttZ:.2f}   $\\pm$  {obsSubttZ_e:.2f}    \\\\ \\hline \\hline
(data-bkg.)/sig. & {rat3l:.2f}   $\\pm$  {rat3l_e:.2f}       & {rat4l:.2f}   $\\pm$  {rat4l_e:.2f}     & {ratttZ:.2f}   $\\pm$  {ratttZ_e:.2f}    \\\\ \\hline

\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(wz3l=wz3l, wz3l_e=wz3l_e, zz4l=zz4l, zz4l_e=zz4l_e,ttZ=ttZ, ttZ_e=ttZ_e,
                            bkg3l=bkg3l, bkg3l_e=bkg3l_e, bkg4l=bkg4l,bkg4l_e=bkg4l_e, bkgttZ=bkgttZ,bkgttZ_e=bkgttZ_e,
                            obs3l=int(obs3l),obs4l=int(obs4l),obsttZ=int(obsttZ),
                            obsSub3l=obsSub3l, obsSub3l_e=obsSub3l_e, obsSub4l=obsSub4l, obsSub4l_e=obsSub4l_e,obsSubttZ=obsSubttZ, obsSubttZ_e=obsSubttZ_e,
                            rat3l=rat3l, rat3l_e=rat3l_e, rat4l=rat4l, rat4l_e=rat4l_e, ratttZ=ratttZ, ratttZ_e=ratttZ_e)
    
    helper.ensurePath('plots/CRs/36.4invfb/tables/')
    compTableFile = open('plots/CRs/36.4invfb/tables/controlRegion3l4l_comparison.txt','w')
    compTableFile.write(table3l4lComparisonString)
    compTableFile.close()

def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, returnplot = False, scaleToData = False, normalized = False, cumulative = False, onlyMC = False):

    print 'returnplot', returnplot
    print 'scaletodata', scaleToData
    print 'normalized', normalized
    print 'cumulative', cumulative
    print 'onlyMC', onlyMC

    theCutDATA = cuts.AddList([theCut, cuts.triggerForCR ])
    theCutZ = cuts.AddList([theCut, cuts.Zmass ])
    MC   = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx, "1", kf)
    MCS  = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx)
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCutDATA, '', labelx, "1", kf)

    if scaleToData:
        print "you're scaling mc to data!! are you sure you wanna do that??"
        newStack = r.THStack()
        mcInt = MC  .Integral()
        daInt = DATA.Integral()
        for _i,_h in enumerate(MCS.GetHists()):
            _hnew = _h.Clone(_h.GetName()+'_new')
            _tmpInt = _h.Integral()
            _hnew.Scale(daInt/mcInt)
            newStack.Add(_hnew)
        MC.Scale(daInt/mcInt)

    for _i,_h in enumerate(MCS.GetHists()):
        errorInte =  r.Double()
        inte = _h.IntegralAndError(1, _h.GetNbinsX()+1, errorInte)
        print "sample ", _h.GetName(), " has integral ", _h.Integral(), " and error ", errorInte
        
    if normalized:
        hists = []
        for _i,_h in enumerate(MCS.GetHists()):
            yie_e = r.Double()
            yie = _h.IntegralAndError(1, _h.GetNbinsX()+1, yie_e)
            _h.SetName(_h.GetName().split('_')[-2]+' %.1f +- %.1f'%(yie, yie_e))
            _h.Scale(1./_h.Integral())
            _h.SetLineColor(_h.GetFillColor())
            _h.SetFillColor(0)
            _h.SetLineWidth(2)
            hists.append(_h)
        yie_e = r.Double()
        yie = MC.IntegralAndError(1, MC.GetNbinsX()+1, yie_e)
        MC.SetName('Full MC')
        MC  .Scale(1./MC.Integral())
        MC  .SetLineColor(r.kRed)
        MC  .SetLineWidth(2)
        DATA.SetName('Data')
        DATA.Scale(1./DATA.Integral())
        DATA.SetMarkerSize(0.8)
        DATA.SetMarkerColor(r.kBlack)
        DATA.SetLineColor(r.kBlack)
        hists.append(MC)
        hists.append(DATA)                                                                   

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
    plot = Canvas.Canvas('CRs/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.5, 0.85, 0.9)
    if not normalized:
        plot.addStack(MCS if not scaleToData else newStack, "HIST", 1, 1)
        plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
    else:
        for h in hists:
            plot.addHisto(h, 'hist,same' if not h.GetName() == 'data' else 'p,same', h.GetName(), 'PL', h.GetLineColor(), 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA, MC)

    if cumulative:
        cum_plot = Canvas.Canvas('DataMC/%s/plot_cumulative_%s'%(lumi_str,name), 'png,pdf,root', 0.60, 0.15, 0.80, 0.45)
        if not normalized:
            cum_plot.addStack(MCS if not scaleToData else newStack, "HIST", 1, 1)
            cum_plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
        else:
            for h in hists:
                if onlyMC and not h.GetName() == 'data':
                    cum_plot.addHisto(h.GetCumulative(), 'hist,same' if not h.GetName() == 'w dPhi cut' else 'p,same', h.GetName(), 'PL' if h.GetName() == 'w dPhi cut' else 'L', h.GetLineColor(), 1, 0)
        cum_plot.save(1, 1, logx, lumi)
    
    if returnplot:
        return plot
    else:
        del plot
        return

def rebinEWino(histo):
    #a= [0,50,100,150,225,300,325]
    a= [50,100,150,225,300,325]
    _arr = array.array('d', a)
    h_ret = r.TH1F(histo.GetName(),histo.GetTitle(), len(a)-1, _arr)
    for ii,ib in enumerate(a):
        if ib == a[-1]: continue
        tmp_err = r.Double()
        if not ib == 300:
            tmp_cont = histo.IntegralAndError(histo.FindBin(ib), histo.FindBin(a[ii+1]-1.), tmp_err )
        else:
            tmp_cont = histo.IntegralAndError(histo.FindBin(ib), histo.GetNbinsX()+1, tmp_err )
        h_ret.SetBinContent(ii+1, tmp_cont)
        h_ret.SetBinError  (ii+1, tmp_err )
    return copy.deepcopy(h_ret)

    
    

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

    mcDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    ttDatasets = ['TTTo2L2Nu_1']
    rareDatasets = ['GGHZZ4L',  'ZZTo4L',  'ZZTo2L2Nu', 'TTZToLLNuNu_ext1','TTLLJets_m1to10', 'TTZToQQ', 'WZTo3LNu' , 'WWTo2L2Nu', 'WZTo2L2Q', 'TWZ', 'WWW', 'WWZ', 'WZZ', 'ZZZ', 'tZq_ll','TTTT','TTWToQQ', 'TTWToLNu_ext2',  'TTHnobb_pow', 'VHToNonbb',  'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'T_tch_powheg', 'TBar_tch_powheg',  'WJetsToLNu_LO']
    mcDatasets += ttDatasets
    mcDatasets += rareDatasets
    kf = "noKFactor"   
    
                                                                              
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



    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isScan = 0)
    #treeSI = Sample.Tree(helper.selectSamples(opts.sampleFile, siDatasets, 'SI'), 'SI', 0, isScan = 1)
    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 36.4 ; maxrun = 999999; lumi_str = '36.4invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
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
    #plot_ttZ_met = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttZregion", 20,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ, cuts.Zmass]), cuts, labelmet, 0, True)
    #plot_ttZ_ee = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttZregion_ee", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ, cuts.ee]), cuts, labelmll, 0, True)
    #plot_ttZ_mm = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttZregion_mm", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ, cuts.mm]), cuts, labelmll, 0, True)
    #plot_ttZ = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttZregion", 3,  61, 121, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZLesya]), cuts, labelmll, 0, True)
    #plot_ttZ_pt = makePlot(lumi, lumi_str, treeDA, treeMC, "Lep1_pt_Edge", "pt_ttZregion", 20,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZLesya]), cuts, "pt", 0, True)
    #plot_ttZ_met = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttZregion", 20,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZLesya]), cuts, labelmet, 0, True)
    #plot_ttZ_njets = makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njets_ttZregion", 9,  0.5, 9.5, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj0, cuts.regionttZ]), cuts, labelnjet, 0, True)
    plot_ttZ_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "mZ1_Edge", "mll_ttZregion_Zmass", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ, cuts.mZ1inZOLD]), cuts, labelmll, 0, True)
    plot_ttZ = makePlot(lumi, lumi_str, treeDA, treeMC, "mZ1_Edge", "mll_ttZregion", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ]), cuts, labelmll, 0, True)
#    plot_ttZ_OldZmass = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttZregion_ZmassOld", 6,  0, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ, cuts.Zmass]), cuts, "pt", 0, True)
#    #plot_ttZ_ZmassOld_met = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttZregion_ZmassOld", 20,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ, cuts.Zmass]), cuts, labelmet, 0, True)
    plot_3l = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion", 14, 60, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l]), cuts, labelmet, 0, True)
    plot_3l_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion_Zmass", 14, 60, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l, cuts.mZ1inZOLD]), cuts, labelmet, 0, True)
#    ##plot_3l_njets = makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge"    , "met_3lregion_njets",9, 0.5, 9.5, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l, cuts.Zmass]), cuts, labelmet, 0, True)
    plot_4l = makePlot(lumi, lumi_str, treeDA, treeMC, "mZ1_Edge", "mll_4lregion", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l]), cuts, labelmll, 0, True)
    plot_4l_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "mZ1_Edge", "mll_4lregion_Zmass", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l, cuts.mZ1inZOLD]), cuts, labelmll, 0, True)
#    #plot_4l_njets = makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "mll_4lregion_njets",9,  0.5, 9.5, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l, cuts.Zmass]), cuts, labelmll, 0, True)
# ## ## ## plot_4l = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_4lregion_met0to30_AF", 10,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l, 'met_Edge < 30']), cuts, labelmll, 0, True)
    makeSummaryTable3l4l(plot_3l_Zmass, plot_4l_Zmass, plot_ttZ_Zmass)
    #print "made the table"
    #ttbar_region = cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2, 'run_Edge <= %d'%maxrun, 'met_Edge >150'])
    #makePlot(lumi, lumi_str, treeDA, treeMC, "nBJetMedium35_Edge", "ttbar_region_met150_2jets_of_dataScaled", 4, 0, 4, ttbar_region, cuts, 'n_{b-jets}', 0, False, True)

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
