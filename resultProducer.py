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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats, TStyle
import math, sys, optparse, copy, re, array


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables


def makeResultsTable(binnedSR, dyShapes, dataMC, eta, nbs):
    line0 = ' %30s &           & '  %('')
    line1 = ' %30s & OF pred.  & ' %('\multirow{3}{*}{%s}' %(binnedSR.name))
    line2 = ' %30s & DY pred.  & ' %('')
    line3 = ' %30s & total     & ' %('')
    line4 = ' %30s & obs.      & ' %('')
    tmp_histo_obs  = binnedSR.mll     .getHisto(dataMC, eta)
    tmp_histo_pred = binnedSR.mll_pred.getHisto(dataMC, eta)
    tmp_histo_dy   = dyShapes['%db_%s_%s_binned'%(nbs, dataMC[:2].lower(), eta)]

    my_range = range(1,tmp_histo_obs.GetNbinsX()+1)
    for i in my_range:
        tmp_dy   = tmp_histo_dy.GetBinContent(i)  ; tmp_dy_e   = tmp_histo_dy.GetBinError(i)
        tmp_of   = tmp_histo_pred.GetBinContent(i); tmp_of_e   = tmp_histo_pred.GetBinError(i)
        tmp_full = tmp_dy + tmp_of                ; tmp_full_e = math.sqrt(tmp_dy_e**2 + tmp_of_e**2)
        tmp_obs  = tmp_histo_obs.GetBinContent(i) ; tmp_obs_e  = tmp_histo_obs.GetBinError(i)
        mll_low , mll_high = tmp_histo_pred.GetXaxis().GetBinLowEdge(i), tmp_histo_pred.GetXaxis().GetBinUpEdge(i)

        line0 += '%.0f $<$ \\mll $<$ %.0f %s' %(mll_low, mll_high   , ' & ' if i != max(my_range) else '\\\\')
        line1 += '  %.2f $\\pm$ %.2f      %s' %(tmp_of  , tmp_of_e  , ' & ' if i != max(my_range) else '\\\\')
        line2 += '  %.2f $\\pm$ %.2f      %s' %(tmp_dy  , tmp_dy_e  , ' & ' if i != max(my_range) else '\\\\')
        line3 += '  %.2f $\\pm$ %.2f      %s' %(tmp_full, tmp_full_e, ' & ' if i != max(my_range) else '\\\\')
        line4 += '  %.2f $\\pm$ %.2f      %s' %(tmp_obs , tmp_obs_e , ' & ' if i != max(my_range) else '\\\\')
    line0 += '\\hline'; line2 += '\\hline'; line3 += '\\hline \\hline'

    return line0, line1, line2, line3, line4

    

def makeRatioTable(myregion, dataMC, eta, nbs): ## for this to make sense the region should be properly binned!!
    header= 'THIS IS THE TABLE FOR %s in %s for %s b-tags'%(dataMC, eta, str(nbs))
    line0 = ' %30s &           & '  %('')
    line1 = ' %30s & %s pred.  & ' %('\multirow{3}{*}{%s}' %(myregion.name), dataMC)
    line2 = ' %30s & %s obs.   & ' %('', dataMC)
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
    return header, line0, line1, line2, line3


def scaleByRSFOF(histo, rsfof, rsfof_err):
    h_rsfof = copy.deepcopy(histo)
    h_rsfof.SetName('h_rsfof')
    for i in range(1, h_rsfof.GetNbinsX()+1):
        h_rsfof.SetBinContent(i,rsfof)
        h_rsfof.SetBinError  (i,rsfof_e)
    histo.Multiply(h_rsfof)
    return histo

def makeResultTable(resultPlotLoNLL, resultPlotHiNLL, lint, lint_str):
    for i,h in enumerate(resultPlotLoNLL.histos):
        if 'da_SF' in h.GetName(): ob_ind = i
        if 'da_OF' in h.GetName(): pr_ind = i
    h_obsLoNLL = resultPlotLoNLL.histos[ob_ind]
    h_preLoNLL = resultPlotLoNLL.histos[pr_ind]
    h_obsHiNLL = resultPlotHiNLL.histos[ob_ind]
    h_preHiNLL = resultPlotHiNLL.histos[pr_ind]
    bin01 = 1
    bin90 = h_obsLoNLL.FindBin(90.)
    binZZ = h_obsLoNLL.GetNbinsX()+1

    ## dyTotal = 11.8; dyTotal_e = 3.1
    #dyTotal = 102.9; dyTotal_e = 14.2
    ## dyTotal = 112.2; dyTotal_e = 34.0
    dyTotal = 188.8; dyTotal_e = 44.2 #12.9 fb-1
    #rinoutLoM = 0.109; rinoutLoM_e = math.sqrt(0.001**2 + 0.027**2)
    #rinoutHiM = 0.063; rinoutHiM_e = math.sqrt(0.001**2 + 0.016**2)
    rinoutLoM = 0.110; rinoutLoM_e = math.sqrt(0.001**2 + 0.027**2)
    rinoutHiM = 0.062; rinoutHiM_e = math.sqrt(0.001**2 + 0.016**2)
    dyEffLoNLL = 0.65
    dyEffHiNLL = 1.-dyEffLoNLL
    
    prFSloNloM_e, prFSloNhiM_e = r.Double(),  r.Double()
    prFShiNloM_e, prFShiNhiM_e = r.Double(),  r.Double()

    #print fapppasdfsf
    prFSloNloM = h_preLoNLL.IntegralAndError(bin01, bin90, prFSloNloM_e)
    prDYloNloM = dyTotal*rinoutLoM*dyEffLoNLL; prDYloNloM_e = dyTotal_e*rinoutLoM*dyEffLoNLL
    prTOloNloM = prFSloNloM+prDYloNloM; prTOloNloM_e = math.sqrt(prFSloNloM_e**2 + prDYloNloM_e**2)
    obTOloNloM = h_obsLoNLL.Integral        (bin01, bin90              )
    prFSloNhiM = h_preLoNLL.IntegralAndError(bin90, binZZ, prFSloNhiM_e)
    prDYloNhiM = dyTotal*rinoutHiM*dyEffLoNLL; prDYloNhiM_e = math.sqrt(prDYloNhiM)
    prTOloNhiM = prFSloNhiM+prDYloNhiM; prTOloNhiM_e = math.sqrt(prFSloNhiM_e**2 + prDYloNhiM_e**2)
    obTOloNhiM = h_obsLoNLL.Integral        (bin90, binZZ              )

    prFShiNloM = h_preHiNLL.IntegralAndError(bin01, bin90, prFShiNloM_e)
    prDYhiNloM = dyTotal*rinoutLoM*dyEffHiNLL; prDYhiNloM_e = math.sqrt(prDYhiNloM)
    prTOhiNloM = prFShiNloM+prDYhiNloM; prTOhiNloM_e = math.sqrt(prFShiNloM_e**2 + prDYhiNloM_e**2)
    obTOhiNloM = h_obsHiNLL.Integral        (bin01, bin90              )
    prFShiNhiM = h_preHiNLL.IntegralAndError(bin90, binZZ, prFShiNhiM_e)
    prDYhiNhiM = dyTotal*rinoutHiM*dyEffHiNLL; prDYhiNhiM_e = math.sqrt(prDYhiNhiM)
    prTOhiNhiM = prFShiNhiM+prDYhiNhiM; prTOhiNhiM_e = math.sqrt(prFShiNhiM_e**2 + prDYhiNhiM_e**2)
    obTOhiNhiM = h_obsHiNLL.Integral        (bin90, binZZ              )
    
    resultTable = '''\\documentclass[12pt,a4paper]{{article}}
\\usepackage{{multirow}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{Predicted and observed yields for {lint} fb$^{{-1}}$ of 2016 data.}} 
\\label{{tab:resultTableData}} 
\\begin{{tabular}}{{r l c c }} 
              &                 & ttbar-like  & non-ttbar-like\\\\ \\hline
\\multirow{{4}}{{*}}{{mll $<$ 81 GeV}}       & pred. FS        & {prFSloNloM:.1f}  $\\pm$  {prFSloNloM_e:.1f}    & {prFShiNloM:.1f}     $\\pm$  {prFShiNloM_e:.1f}  \\\\
                                             & pred. DY        & {prDYloNloM:.1f}  $\\pm$  {prDYloNloM_e:.1f}    & {prDYhiNloM:.1f}     $\\pm$  {prDYhiNloM_e:.1f}  \\\\
                                             & pred. total     & {prTOloNloM:.1f}  $\\pm$  {prTOloNloM_e:.1f}    & {prTOhiNloM:.1f}     $\\pm$  {prTOhiNloM_e:.1f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNloM}}}                        & \\textbf{{{obTOhiNloM}}}                         \\\\ \\hline
\\multirow{{4}}{{*}}{{mll $>$ 101 GeV}}      & pred. FS     & {prFSloNhiM:.1f}  $\\pm$  {prFSloNhiM_e:.1f}       & {prFShiNhiM:.1f}     $\\pm$  {prFShiNhiM_e:.1f}  \\\\
                                             & pred. DY     & {prDYloNhiM:.1f}  $\\pm$  {prDYloNhiM_e:.1f}       & {prDYhiNhiM:.1f}     $\\pm$  {prDYhiNhiM_e:.1f}  \\\\
                                             & pred. total  & {prTOloNhiM:.1f}  $\\pm$  {prTOloNhiM_e:.1f}       & {prTOhiNhiM:.1f}     $\\pm$  {prTOhiNhiM_e:.1f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNhiM}}}                        & \\textbf{{{obTOhiNhiM}}}                         \\\\
\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(
    prFSloNloM_e=prFSloNloM_e, prDYloNloM_e=prDYloNloM_e, prTOloNloM_e=prTOloNloM_e, prFSloNhiM_e=prFSloNhiM_e, prDYloNhiM_e=prDYloNhiM_e, prTOloNhiM_e=prTOloNhiM_e,
    prFShiNloM_e=prFShiNloM_e, prDYhiNloM_e=prDYhiNloM_e, prTOhiNloM_e=prTOhiNloM_e, prFShiNhiM_e=prFShiNhiM_e, prDYhiNhiM_e=prDYhiNhiM_e, prTOhiNhiM_e=prTOhiNhiM_e,

    prFSloNloM = prFSloNloM,      prFShiNloM = prFShiNloM,
    prDYloNloM = prDYloNloM,      prDYhiNloM = prDYhiNloM,
    prTOloNloM = prTOloNloM,      prTOhiNloM = prTOhiNloM,
    obTOloNloM = int(obTOloNloM), obTOhiNloM = int(obTOhiNloM),
    prFSloNhiM = prFSloNhiM,      prFShiNhiM = prFShiNhiM,
    prDYloNhiM = prDYloNhiM,      prDYhiNhiM = prDYhiNhiM,
    prTOloNhiM = prTOloNhiM,      prTOhiNhiM = prTOhiNhiM,
    obTOloNhiM = int(obTOloNhiM), obTOhiNhiM = int(obTOhiNhiM), lint=lint)

    helper.ensureDirectory('plots/results/%s/'%lint_str); 
    helper.ensureDirectory('plots/results/%s/tables/'%lint_str)
    tablename = (resultPlotLoNLL.name.split('_')[-1]).replace('nllABove21','').replace('nllBelow21','')
    compTableFile = open('plots/results/%s/tables/resultTable_%s%s.tex'%(lint_str, str(lint).replace('.','p'), tablename),'w')
    compTableFile.write(resultTable)
    compTableFile.close()

def makeClosureTests(var, specialcut = '', scutstring = '', doCumulative = False):
    #rsfof_mc = 1.060; rsfof_mc_err = 0.046

    if var == 'mll':
        treevar = 'lepsMll_Edge'
        nbins, xmin, xmax = 23, 20, 250
        xlabel = 'm_{ll} (GeV)'
    elif var == 'nll':
        treevar = 'nll_Edge'
        nbins, xmin, xmax = 26, 10, 36
        xlabel = 'NLL'
    elif var == 'nllMC':
        treevar = 'nll_mc_Edge'
        nbins, xmin, xmax = 26, 10, 36
        xlabel = 'NLL'
    elif var == 'nllMCSF':
        treevar = 'nll_mc_sf_Edge'
        nbins, xmin, xmax = 26, 10, 36
        xlabel = 'NLL'
        

    ## ## mll ditributions
    mc_OF = treeMC.getTH1F(lumi, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.OF, cuts.Zveto]), '', xlabel)
    mc_SF = treeMC.getTH1F(lumi, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.SF, cuts.Zveto]), '', xlabel)
    da_OF = treeDA.getTH1F(lumi, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.OF, cuts.Zveto]), '', xlabel)
    da_SF = treeDA.getTH1F(lumi, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.SF, cuts.Zveto]), '', xlabel)
    dy_SF = treeDY.getTH1F(lumi, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.SF, cuts.Zveto]), '', xlabel)
    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_err.SetFillStyle(3004); mc_OF_err.SetMarkerSize(0.)
    dy_SF.SetFillColorAlpha(r.kGreen+2,0.5)

    ## da_mll_nllInc_OF_err = copy.deepcopy(da_mll_nllInc_OF)
    ## da_mll_nllInc_OF_err.SetFillColorAlpha(r.kBlack, 0.8)
    ## da_mll_nllInc_OF_err.SetFillStyle(3004); da_mll_nllInc_OF_err.SetMarkerSize(0.)

    mc_SF.GetYaxis().SetRangeUser(0., 1.3*mc_SF.GetMaximum())
    print helper.bcolors.HEADER + '[MC only closure test not scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure_noRSFOF = Canvas.Canvas('closure/%s/plot_closure_%s_mcPredmcObs%s_noRSFOF'%(lumi_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_closure_noRSFOF.addHisto(mc_SF    , 'PE'       , 'MC - SF', 'PL', r.kRed+1  , 1,  0)
    plot_closure_noRSFOF.addHisto(mc_OF_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure_noRSFOF.addHisto(mc_OF    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_closure_noRSFOF.addHisto(dy_SF    , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
    plot_closure_noRSFOF.addLatex (0.2, 0.8, 'no R_{SFOF} scaling')
    plot_closure_noRSFOF.saveRatio(1, 0, 1, lumi, mc_SF, mc_OF, 0.2, 1.8)

    mc_OF_rsfofScaled = copy.deepcopy(mc_OF)
    mc_OF_rsfofScaled = scaleByRSFOF(mc_OF_rsfofScaled, rsfof_mc, rsfof_mc_e)

    mc_OF_rsfofScaled_err = copy.deepcopy(mc_OF_rsfofScaled)
    mc_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_rsfofScaled_err.SetFillStyle(3004); mc_OF_rsfofScaled_err.SetMarkerSize(0.)

    print helper.bcolors.HEADER + '[MC only closure test scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure = Canvas.Canvas('closure/%s/plot_closure_%s_mcPredmcObs%s'%(lumi_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_closure.addHisto(mc_SF                , 'PE'       , 'MC - SF', 'PL', r.kRed+1  , 1,  0)
    plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_closure.addHisto(dy_SF                , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
    plot_closure.addLatex (0.2, 0.8, 'R_{SFOF} scaled')
    plot_closure.saveRatio(1, 0, 0, lumi, mc_SF, mc_OF_rsfofScaled, 0.2, 1.8)

    if True:
        mc_OF_rsfofScaled    .Scale(0.8/10.)
        mc_OF_rsfofScaled_err.Scale(0.8/10.)
        dy_SF                .Scale(0.8/10.)
        plot_closure = Canvas.Canvas('closure/%s/plot_closure_%s_mcPreddaObs%s'%(lumi_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
        plot_closure.addHisto(da_SF                , 'PE'       , 'data-SF', 'PL', r.kRed+1  , 1,  0)
        plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
        plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
        plot_closure.addHisto(dy_SF                , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
        plot_closure.addLatex (0.2, 0.8, 'R_{SFOF} scaled')
        plot_closure.saveRatio(1, 0, 0, 0.8 , mc_SF, mc_OF_rsfofScaled, 0.2, 1.8)

    ## make cumulative distributions
    if doCumulative:
        mc_SF_cum                 = copy.deepcopy(mc_SF                ).GetCumulative()
        mc_OF_rsfofScaled_err_cum = copy.deepcopy(mc_OF_rsfofScaled_err).GetCumulative()
        mc_OF_rsfofScaled_cum     = copy.deepcopy(mc_OF_rsfofScaled    ).GetCumulative()
        dy_SF_cum                 = copy.deepcopy(dy_SF                ).GetCumulative()
        mc_SF_cum            .Scale(1./mc_SF            .Integral()); mc_SF_cum            .SetLineWidth(2)
        mc_OF_rsfofScaled_cum.Scale(1./mc_OF_rsfofScaled.Integral()); mc_OF_rsfofScaled_cum.SetLineWidth(2)
        dy_SF_cum            .Scale(1./dy_SF            .Integral()); dy_SF_cum            .SetLineWidth(2); dy_SF_cum.SetFillColor(r.kWhite)
        mc_SF_cum.GetYaxis() .SetRangeUser(0.9, 1.01)
        print helper.bcolors.HEADER + '[MC only cumulative distribution, scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
        plot_cumulative = Canvas.Canvas('closure/%s/plot_cumulative_%s_mcPredmcObs%s'%(lumi_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.8, 0.2, 0.95, 0.4)
        plot_cumulative.addHisto(mc_SF_cum                 , 'hist,same', 'MC - SF', 'L' , r.kRed+1  , 1,  0)
        plot_cumulative.addHisto(mc_OF_rsfofScaled_cum     , 'hist,SAME', 'MC - OF', 'L', r.kBlue+1 , 1,  1)
        #plot_cumulative.addHisto(dy_SF_cum                 , 'hist,SAME', 'DY - SF', 'L', r.kGreen+2, 1,  2)
        plot_cumulative.addLatex (0.2, 0.8, 'R_{SFOF} scaled')
        plot_cumulative.saveRatio(1, 0, 0, lumi, mc_SF_cum, mc_OF_rsfofScaled_cum, 0.2, 1.8)
        return mc_SF_cum

def makeResultData(var, maxrun = 274240, lint = 0.864, specialcut = '', scutstring = '', _options = ''):
    returnplot, addRares, splitFlavor, makeTable, printIntegral = False, False, False, False, False
    if 'returnplot'    in _options: print 'found option %s'%'returnplot'    ;returnplot    = True
    if 'addRares'      in _options: print 'found option %s'%'addRares'      ;addRares      = True
    if 'splitFlavor'   in _options: print 'found option %s'%'splitFlavor'   ;splitFlavor   = True
    if 'makeTable'     in _options: print 'found option %s'%'makeTable'     ;makeTable     = True
    if 'printIntegral' in _options: print 'found option %s'%'printIntegral' ;printIntegral = True

    # if   var == 'mll'      : treevar = 'lepsMll_Edge'        ; nbins, xmin, xmax = 23, 20 , 250 ; xlabel = 'm_{ll} (GeV)'
    if   var == 'mll'      : treevar = 'lepsMll_Edge'        ; nbins, xmin, xmax = 20, 20 , 320 ; xlabel = 'm_{ll} (GeV)'
    elif var == 'nll'      : treevar = 'nll_Edge'            ; nbins, xmin, xmax = 15, 12 , 27  ; xlabel = 'NLL'
    elif var == 'nb'       : treevar = 'nBJetMedium35_Edge'  ; nbins, xmin, xmax =  3,  0 ,  3  ; xlabel = 'n_{b-jets,35}'
    elif var == 'nj'       : treevar = 'nJetSel_Edge'        ; nbins, xmin, xmax =  6, 0.5, 6.5 ; xlabel = 'n_{jets}'
    elif var == 'zpt'      : treevar = 'lepsZPt_Edge'        ; nbins, xmin, xmax = 10,  0 ,1000 ; xlabel = 'p_{T}^{ll}'
    elif var == 'mlb'      : treevar = 'sum_mlb_Edge'        ; nbins, xmin, xmax = 15,  0 ,1500 ; xlabel = '#Sigma m_{lb}'
    elif var == 'met'      : treevar = 'met_Edge'            ; nbins, xmin, xmax = 10,100 ,1000 ; xlabel = 'E_{T}^{miss.}'
    elif var == 'metraw'   : treevar = 'met_raw_Edge'        ; nbins, xmin, xmax = 10,100 ,1000 ; xlabel = 'E_{T}^{miss.} raw'
    elif var == 'ldp'      : treevar = 'abs(lepsDPhi_Edge)'  ; nbins, xmin, xmax = 10,  0 , 3.15; xlabel = '#Delta #phi ll'
    elif var == 'pt1'      : treevar = 'Lep1_pt_Edge'        ; nbins, xmin, xmax = 20,  0 , 500 ; xlabel = 'p_{T} leading'
    elif var == 'pt2'      : treevar = 'Lep2_pt_Edge'        ; nbins, xmin, xmax = 10, 20 , 200 ; xlabel = 'p_{T} trailing'
    elif var == 'ldr'      : treevar = 'lepsDR_Edge'         ; nbins, xmin, xmax = 10,  0 , 6.  ; xlabel = '#Delta R (ll)'
    elif var == 'iso1'     : treevar = 'Lep1_miniRelIso_Edge'; nbins, xmin, xmax = 10,  0 , 0.05; xlabel = 'mini iso l1'
    elif var == 'iso2'     : treevar = 'Lep2_miniRelIso_Edge'; nbins, xmin, xmax = 10,  0 , 0.05; xlabel = 'mini iso l2'
    elif var == 'l1reliso03'     : treevar = 'Lep1_relIso03_Edge'; nbins, xmin, xmax = 10,  0 , 0.20; xlabel = 'rel.Iso 03 leading'
    elif var == 'l2reliso03'     : treevar = 'Lep2_relIso03_Edge'; nbins, xmin, xmax = 10,  0 , 0.20; xlabel = 'rel.Iso 03 trailing'
    elif var == 'l1reliso04'     : treevar = 'Lep1_relIso04_Edge'; nbins, xmin, xmax = 10,  0 , 0.20; xlabel = 'rel.Iso 04 leading'
    elif var == 'l2reliso04'     : treevar = 'Lep2_relIso04_Edge'; nbins, xmin, xmax = 10,  0 , 0.20; xlabel = 'rel.Iso 04 trailing'
    elif var == 'eta1'     : treevar = 'Lep1_eta_Edge'       ; nbins, xmin, xmax = 10,-2.5, 2.5 ; xlabel = '#eta leading'
    elif var == 'eta2'     : treevar = 'Lep2_eta_Edge'       ; nbins, xmin, xmax = 10,-2.5, 2.5 ; xlabel = '#eta trailing'
    elif var == 'l1metdphi': treevar = 'abs(metl1DPhi_Edge)' ; nbins, xmin, xmax = 10,  0., 3.15; xlabel = '#Delta #phi_{MET,lead.}'
    elif var == 'l2metdphi': treevar = 'abs(metl2DPhi_Edge)' ; nbins, xmin, xmax = 10,  0., 3.15; xlabel = '#Delta #phi_{MET,trail.}'
    elif var == 'mt2'      : treevar = 'mt2_Edge'            ; nbins, xmin, xmax = 10,  0 , 100.; xlabel = 'M_{T2}^{ll}'
    elif var == 'l1dxy'    : treevar = 'abs(Lep1_dxy_Edge)'  ; nbins, xmin, xmax = 10,  0., 0.03; xlabel = 'd_{xy} lead.'
    elif var == 'l2dxy'    : treevar = 'abs(Lep2_dxy_Edge)'  ; nbins, xmin, xmax = 10,  0., 0.03; xlabel = 'd_{xy} trail.'
    elif var == 'l1dz'     : treevar = 'abs(Lep1_dz_Edge)'   ; nbins, xmin, xmax = 10,  0., 0.03; xlabel = 'd_{z} lead.'
    elif var == 'l2dz'     : treevar = 'abs(Lep2_dz_Edge)'   ; nbins, xmin, xmax = 10,  0., 0.03; xlabel = 'd_{z} trail.'
    elif var == 'l1sip3d'  : treevar = 'Lep1_sip3d_Edge'     ; nbins, xmin, xmax = 10,  0., 5.  ; xlabel = 'SIP-3D lead.'
    elif var == 'l2sip3d'  : treevar = 'Lep2_sip3d_Edge'     ; nbins, xmin, xmax = 10,  0., 5.  ; xlabel = 'SIP-3D trail.'

    elif var == 'nll_noMET': treevar = '-1.*TMath::Log(lh_ana_ldp_data_Edge*lh_ana_zpt_data_Edge*lh_ana_mlb_data_Edge)'; nbins, xmin, xmax = 20, 10, 30   ; xlabel = 'NLL - no MET pdf'
    elif var == 'nll_noMLB': treevar = '-1.*TMath::Log(lh_ana_ldp_data_Edge*lh_ana_zpt_data_Edge*lh_ana_met_data_Edge)'; nbins, xmin, xmax = 20, 10, 30   ; xlabel = 'NLL - no MLB pdf'
    elif var == 'nll_noZPT': treevar = '-1.*TMath::Log(lh_ana_ldp_data_Edge*lh_ana_mlb_data_Edge*lh_ana_met_data_Edge)'; nbins, xmin, xmax = 20, 10, 30   ; xlabel = 'NLL - no ZPT pdf'
    elif var == 'nll_noLDP': treevar = '-1.*TMath::Log(lh_ana_zpt_data_Edge*lh_ana_mlb_data_Edge*lh_ana_met_data_Edge)'; nbins, xmin, xmax = 20, 10, 30   ; xlabel = 'NLL - no ZPT pdf'
        

    newLumiString = str(lint)+'invfb'
    if not specialcut:
        specialcut = 'run_Edge <= {run}'.format(run=maxrun)
    else:
        specialcut += ' && run_Edge <= {run}'.format(run=maxrun)
    ## ## mll ditributions
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.SF, cuts.Zveto]), '', xlabel)
    da_mm = treeDA.getTH1F(lint, var+"da_mm"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.mm, cuts.Zveto]), '', xlabel)
    da_ee = treeDA.getTH1F(lint, var+"da_ee"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.ee, cuts.Zveto]), '', xlabel)
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.OF, cuts.Zveto]), '', xlabel)
    if addRares:
        ra_OF = treeRA.getTH1F(lint, var+"ra_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.OF, cuts.Zveto]), '', xlabel)
        ra_SF = treeRA.getTH1F(lint, var+"ra_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.SF, cuts.Zveto]), '', xlabel)
        ra_SF.SetFillColorAlpha(r.kRed+1, 0.8)
        ra_SF.SetFillStyle(3017); ra_SF.SetMarkerSize(0.)
        ra_SF.Add(ra_OF, -1.)
    da_OF_err = copy.deepcopy(da_OF)
    da_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    da_OF_err.SetFillStyle(3017); da_OF_err.SetMarkerSize(0.)

    da_OF_rsfofScaled = copy.deepcopy(da_OF)
    da_OF_rsfofScaled.SetName(da_OF.GetName()+'_scaledRSFOF')
    da_OF_rsfofScaled = scaleByRSFOF(da_OF_rsfofScaled, rsfof, rsfof_e)
    da_OF_rsfofScaled.SetLineWidth(2)

    da_OF_rsfofScaled_err = copy.deepcopy(da_OF_rsfofScaled)
    da_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    da_OF_rsfofScaled_err.SetFillStyle(3017); da_OF_rsfofScaled_err.SetMarkerSize(0.)

    maxrat = 0.5
    for ib in range(1,da_SF.GetNbinsX()+1):
        tmp_rat = da_SF.GetBinContent(ib)/( da_OF_rsfofScaled.GetBinContent(ib) if da_OF_rsfofScaled.GetBinContent(ib) > 0 else 1. )
        if tmp_rat > maxrat:
            maxrat = tmp_rat
    
    maxCont = max(da_OF_rsfofScaled_err.GetMaximum(), da_SF.GetMaximum())
    da_OF_rsfofScaled_err.GetYaxis().SetRangeUser(0.1, 1.30*maxCont)
    da_SF                .GetYaxis().SetRangeUser(0.1, 1.30*maxCont)

    print helper.bcolors.HEADER + '[result scaled by RSFOF for DATA] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    sstring = '' if not addRares else 'withRaresFromMC'
    plot_result = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s'%(newLumiString+sstring, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.60, 0.65, 0.80, 0.85)
    plot_result.addHisto(da_OF_rsfofScaled_err, 'e2,same'  , ''         , 'PL', r.kBlue+1 , 1, -1)
    plot_result.addHisto(da_OF_rsfofScaled    , 'hist,SAME', 'predicted', 'L' , r.kBlue+1 , 1,  1)
    if addRares:
        plot_result.addHisto(ra_SF, 'hist,SAME', 'rares - SF', 'L' , r.kRed+1 , 1,  2)
    if splitFlavor:
        plot_result.addHisto(da_ee, 'pe,same', 'data- elel'  , 'PL', r.kYellow+1  , 1,  2)
        plot_result.addHisto(da_mm, 'pe,same', 'data- #mu#mu', 'PL', r.kGreen+1   , 1,  3)
    plot_result.addHisto(da_SF                , 'PE,same'    , 'observed data', 'PL', r.kBlack  , 1,  0)
    #plot_result.addLatex (0.2, 0.80, 'R_{SFOF} scaled')
    if maxrun < 999999: plot_result.addLatex (0.2, 0.85, 'max. run {run}'.format(run=maxrun))
    if printIntegral  : plot_result.addLatex (0.2, 0.65, 'SF: {SF:.1f} (#mu#mu:{mm:.1f},ee:{ee:.1f})'.format(SF=da_SF.Integral(),mm=da_mm.Integral(),ee=da_ee.Integral() ) )
    if printIntegral  : plot_result.addLatex (0.2, 0.60, 'OF: {OF:.2f}'.format(OF=da_OF_rsfofScaled.Integral() ) )
    plot_result.saveRatio(1, 1, 0, lint, da_SF, da_OF_rsfofScaled, 0. , int(maxrat+1.0) )

    if makeTable:
        makeSimpleTable(plot_result, addRares)
    if returnplot:
        return plot_result

def makeSimpleTable(plot, addRares=False):
    for i,h in enumerate(plot.histos):
        if 'da_OF' in h.GetName(): fs_ind = i
        if 'da_SF' in h.GetName(): ob_ind = i
        if 'da_ee' in h.GetName(): ee_ind = i
        if 'da_mm' in h.GetName(): mm_ind = i
        if 'ra_SF' in h.GetName(): ra_ind = i

    h_fs= plot.histos[fs_ind]
    h_ee= plot.histos[ee_ind]
    h_mm= plot.histos[mm_ind]
    h_ob= plot.histos[ob_ind]

    ob_e, fs_e, ra_e, mm_e, ee_e = r.Double(), r.Double(), r.Double(), r.Double(), r.Double()

    ob = h_ob.IntegralAndError(1, h_ob.GetNbinsX()+1, ob_e)
    mm = h_mm.IntegralAndError(1, h_mm.GetNbinsX()+1, ob_e)
    ee = h_ee.IntegralAndError(1, h_ee.GetNbinsX()+1, ob_e)
    fs = h_fs.IntegralAndError(1, h_ob.GetNbinsX()+1, fs_e)
    if addRares: 
        h_ra= plot.histos[ra_ind]
        ra = h_ra.IntegralAndError(1, h_ra.GetNbinsX()+1, ra_e) 
    else: 
        ra, ra_e = 0., 0.

    pr = fs+ra
    pr_e = math.sqrt(fs_e*fs_e + ra_e*ra_e)

    sign = (ob-pr)/math.sqrt(ob)

    regionDesc = plot.name.split('daObs_')[-1]
    tableName='plots/'+plot.name.replace('plot','table')+'.tex'

    resultTable = '''\\documentclass[12pt,a4paper]{{article}}
\\usepackage{{multirow}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{Predicted and observed yields for {lint} fb$^{{-1}}$ of 2016 data. In region: {regionDesc}. Significance is just sig/sqrt(sig+bkg)}} 
\\label{{tab:resultTableData}} 
\\begin{{tabular}}{{l c }} 
    pred. FS         & {fs:.2f}     $\\pm$  {fs_e:.2f}  \\\\
    pred. rare       & {ra:.2f}     $\\pm$  {ra_e:.2f}  \\\\
    pred. total      & {pr:.2f}     $\\pm$  {pr_e:.2f}  \\\\
    \\textbf{{ob}}  & \\textbf{{{ob}}}          \\\\ \\hline
    \\textbf{{$\\mu\\mu$}}  & \\textbf{{{mm}}}          \\\\
    \\textbf{{$ee$      }}  & \\textbf{{{ee}}}          \\\\ \\hline
    sign.                   & {sig:.2f}                     \\\\
\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format( regionDesc = regionDesc.replace('_', ' '), ob = ob, pr=pr, ra=ra, fs=fs, mm=mm, ee=ee,
                                                 pr_e=pr_e, ra_e=ra_e, fs_e=fs_e, lint='4fb-1', sig=sign)

    #helper.ensureDirectory('plots/results/%s/'%lint_str); 
    #helper.ensureDirectory('plots/results/%s/tables/'%lint_str)
    print 'printing table', tableName
    compTableFile = open(tableName,'w')
    compTableFile.write(resultTable)
    compTableFile.close()

    

def makePlotsCombinedSR(srlist):
    for i,sR in enumerate(srlist):
        if not i: newSR = copy.deepcopy(sR)
        else:
            newSR = newSR.add(sR)
    print 'this is newSR.mll', newSR.mll
    print 'this is newSR.mll_pred', newSR.mll_pred
    print 'this is newSR.mll_pred.mc', newSR.mll_pred.mc
    makePlots([newSR])

if __name__ == '__main__':

    print 'Starting to produce some good ol\' results...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples FilenameWithIngredients', version='%prog 1.0')
    parser.add_option('-t', action='store_true', dest='onlyTTbar', default=False, help='just do OF closure test')
    parser.add_option('-c', action='store_true', dest='onlyClosure', default=False, help='just do the closure test. don\'t bother with data')
    parser.add_option('-b', '--nbs', action='store', type=int, dest='nbs', default=0, help='do this for different numbers of b\'s')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-l', '--loadShapes', action='store_true', dest='loadShapes', default=False, help='reload dy shapes. default is off since this takes a while')
    parser.add_option('-M', '--maxRun', action='store', type=int, dest='maxRun', default=999999, help='max run to use for analysis (run is included)')
    parser.add_option('-m', '--minRun', action='store', type=int, dest='minRun', default=-1    , help='min run to use for analysis (run not included)')

    ## make the options globa.. also the lumi
    global opts, lumi, lumi_str, dy_shapes, nbstring
    global ingMC, ingDA, onZ, treeDA, treeMC, treeDY, treeTT
    (opts, args) = parser.parse_args()

    print 'running with these options \n'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'

    ## redp this for 2016 ingMC = helper.ingredients(opts.ingredientFile, 'MC'  )
    ## redp this for 2016 ingDA = helper.ingredients(opts.ingredientFile, 'DATA')

    ## rsfofTable = Tables.makeRSFOFTable(ingDA, ingMC)

    ##print asdf
    print 'Going to load DATA and MC trees...'
    dyDatasets = ['DYJetsToLL_M10to50', 'DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext', 'DYJetsToLL_M50_HT400to600_ext', 'DYJetsToLL_M50_HT600toInf_ext', 'WZTo3LNu', 'TTWToLNu', 'TTZToLLNuNu']
    #ttDatasets = ['TTJets_DiLepton_total']
    ttDatasets = ['TT_pow_ext34']
    mcDatasets = ttDatasets + ([] if opts.onlyTTbar else dyDatasets)
    raDatasets = ['WZTo3LNu', 'TTWToLNu', 'TTZToLLNuNu', 'ZZ', 'WZTo2L2Q', 'T_tWch', 'TBar_tWch']
    #daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_271036_276097', 'DoubleEG_Run2016B-PromptReco-v2_runs_271036_276097', 'MuonEG_Run2016B-PromptReco-v2_runs_271036_276097']
    #daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_271036_276097', 'DoubleEG_Run2016B-PromptReco-v2_runs_271036_276097', 'MuonEG_Run2016B-PromptReco-v2_runs_271036_276097',
    #              'DoubleMuon_Run2016C-PromptReco-v2_runs_275420_276098', 'DoubleEG_Run2016C-PromptReco-v2_runs_275420_276098', 'MuonEG_Run2016C-PromptReco-v2_runs_275420_276098']
    #daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_271036_276097', 'DoubleEG_Run2016B-PromptReco-v2_runs_271036_276097', 'MuonEG_Run2016B-PromptReco-v2_runs_271036_276097',
    #              'DoubleMuon_Run2016C-PromptReco-v2_runs_271036_276811', 'DoubleEG_Run2016C-PromptReco-v2_runs_271036_276811', 'MuonEG_Run2016C-PromptReco-v2_runs_271036_276811',
    #              'DoubleMuon_Run2016D-PromptReco-v2_runs_271036_276811', 'DoubleEG_Run2016D-PromptReco-v2_runs_271036_276811', 'MuonEG_Run2016D-PromptReco-v2_runs_271036_276811']
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_273150_275376', 'DoubleEG_Run2016B-PromptReco-v2_runs_273150_275376', 'MuonEG_Run2016B-PromptReco-v2_runs_273150_275376',
                  'DoubleMuon_Run2016C-PromptReco-v2_runs_275420_276283', 'DoubleEG_Run2016C-PromptReco-v2_runs_275420_276283', 'MuonEG_Run2016C-PromptReco-v2_runs_275420_276283',
                  'DoubleMuon_Run2016D-PromptReco-v2_runs_276315_276811', 'DoubleEG_Run2016D-PromptReco-v2_runs_276315_276811', 'MuonEG_Run2016D-PromptReco-v2_runs_276315_276811']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0)
    treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT'  , 0)
    treeRA = Sample.Tree(helper.selectSamples(opts.sampleFile, raDatasets, 'RA'), 'RA'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'

    lumi     = 10.
    lumi_str = str(lumi).replace('.', 'p')+'invfb'
    print 'Running with an integrated luminosity of %.2f fb-1' %(lumi)

    isBlinded = True

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()
    #lint = 0.8 ; maxrun = 274240; lint_str = '0.8invfb'
    #lint = 4.0 ; maxrun = 999999; lint_str = '4.0invfb'
    #lint = 7.65 ; maxrun = 999999; lint_str = '7.65invfb'
    lint = 12.9  ; maxrun = 999999; lint_str = '12.9invfb'

    ## ============================================================
    ## ========== set RSFOF globally ==============================
    ## ============================================================
    global rsfof, rsfof_e, rsfof_mc, rsfof_mc_e
    #rsfof   = 1.102
    #rsfof_e = 0.076
    rsfof   = 1.085
    rsfof_e = 0.060

    rsfof_mc   = 1.070
    rsfof_mc_e = 0.050

    ## result plots in different variables:
    ## ===========================================
    #for v in ['mll']:#'l2sip3d']:#'l1metdphi', 'l2metdphi']:#'eta1', 'eta2']: #'iso1', 'iso2', 'mll', 'nll', 'nb', 'nj', 'zpt', 'mlb', 'met', 'ldp', 'pt1', 'pt2']:
    #    resultPlot      = makeResultData(v , maxrun , lint , specialcut = ''                                      , scutstring = ''                   , _options='returnplot,splitFlavor,printIntegral,makeTable')
    #    resultPlotLoNLL = makeResultData(v , maxrun , lint , specialcut = 'nll_Edge > 21.'                        , scutstring = 'nllAbove21'         , _options='returnplot,splitFlavor,printIntegral,makeTable')
    #    resultPlotLoNLL = makeResultData(v , maxrun , lint , specialcut = 'nll_Edge < 21.'                        , scutstring = 'nllBelow21'         , _options='returnplot,splitFlavor,printIntegral,makeTable')
    #    ## e = makeResultData(v , maxrun , lint , specialcut = 'nll_Edge > 21. && lepsMll_Edge > 101.' , scutstring = 'highMassnllAbove21' , _options='returnplot,splitFlavor,printIntegral,makeTable')
    #    ## a = makeResultData(v , maxrun , lint , specialcut = '                  lepsMll_Edge > 101.' , scutstring = 'highMassnllIncl'    , _options='returnplot,splitFlavor,printIntegral,makeTable')
    #    ## b = makeResultData(v , maxrun , lint , specialcut = '                  lepsMll_Edge <  81.' , scutstring = 'lowMassnllIncl'     , _options='returnplot,splitFlavor,printIntegral,makeTable')
    #    ## f = makeResultData(v , maxrun , lint , specialcut = 'nll_Edge > 21. && lepsMll_Edge < 81.'  , scutstring = 'lowMassnllAbove21'  , _options='returnplot,splitFlavor,printIntegral,makeTable')
    ## resultPlotHiNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge > 21.' , scutstring = 'nllAbove21'         , _options='returnplot,splitFlavor,printIntegral,makeTable')
    ## resultPlotLoNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge < 21.' , scutstring = 'nllBelow21'         , _options='returnplot,splitFlavor,printIntegral,makeTable')
    ## makeResultTable(resultPlotLoNLL, resultPlotHiNLL, lint, lint_str)
    # resultPlotHiNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge > 21. && run_Edge <= 275125' , scutstring = 'nllAbove21_first4p3invfb', _options='returnplot,splitFlavor,printIntegral,makeTable')
    # resultPlotLoNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge < 21. && run_Edge <= 275125' , scutstring = 'nllBelow21_first4p3invfb', _options='returnplot,splitFlavor,printIntegral,makeTable')
    # makeResultTable(resultPlotLoNLL, resultPlotHiNLL, lint, lint_str)
    # resultPlotHiNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge > 21. && run_Edge > 275125' , scutstring = 'nllAbove21_second3p4invfb', _options='returnplot,splitFlavor,printIntegral,makeTable')
    # resultPlotLoNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge < 21. && run_Edge > 275125' , scutstring = 'nllBelow21_second3p4invfb', _options='returnplot,splitFlavor,printIntegral,makeTable')
    # makeResultTable(resultPlotLoNLL, resultPlotHiNLL, lint, lint_str)
    ## resultPlotHiNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge > 21. && run_Edge > 276097' , scutstring = 'nllAbove21_last5p25invfb', _options='returnplot,splitFlavor,printIntegral,makeTable')
    ## resultPlotLoNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge < 21. && run_Edge > 276097' , scutstring = 'nllBelow21_last5p25invfb', _options='returnplot,splitFlavor,printIntegral,makeTable')
    ## makeResultTable(resultPlotLoNLL, resultPlotHiNLL, lint, lint_str)
    ##   resultPlotHiNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge > 21.' , scutstring = 'nllAbove21', _options='returnplot,splitFlavor,printIntegral,makeTable')
    ##   resultPlotLoNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge < 21.' , scutstring = 'nllBelow21', _options='returnplot,splitFlavor,printIntegral,makeTable')
    ##   makeResultTable(resultPlotLoNLL, resultPlotHiNLL, lint, lint_str)
    resultPlotLoNLL = makeResultData('nll' , maxrun , lint , specialcut = ''               , scutstring = 'mllInc'    , _options='returnplot')
    resultPlotHiNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge > 21.' , scutstring = 'nllAbove21', _options='returnplot')
    resultPlotLoNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge < 21.' , scutstring = 'nllBelow21', _options='returnplot')
    print ads
    makeResultTable(resultPlotLoNLL, resultPlotHiNLL, lint, lint_str)

    ## for v in ['nj']:
    ##     makeResultData(v , maxrun , lint , specialcut = 'met_Edge > 250 && nBJetMedium35_Edge ==0 && lepsMll_Edge > 101', scutstring = 'ASRmet250_nb0_highMass', returnplot = True , addRares = True, splitFlavor = True, makeTable = True)
    ##     makeResultData(v , maxrun , lint , specialcut = 'met_Edge > 250 && nBJetMedium35_Edge ==1 && lepsMll_Edge > 101', scutstring = 'ASRmet250_nb1_highMass', returnplot = True , addRares = True, splitFlavor = True, makeTable = True)
    ##     makeResultData(v , maxrun , lint , specialcut = 'met_Edge > 250 && nBJetMedium35_Edge >=2 && lepsMll_Edge > 101', scutstring = 'ASRmet250_nb2_highMass', returnplot = True , addRares = True, splitFlavor = True, makeTable = True)

    ## for v in ['mll']:#'nll_noMET', 'nll_noMLB', 'nll_noZPT', 'nll_noLDP']: #'iso1', 'iso2', 'mll', 'nll', 'nb', 'nj', 'zpt', 'mlb', 'met', 'ldp', 'pt1', 'pt2']:
    ##     for _i,(_min,_max) in enumerate([(0, 274250), (274250,274388), (274388,274971), (274971,275125)]):
    ##         _rstr = 'run_Edge > {_minrun} && run_Edge <= {_maxrun}'.format(_minrun=_min,_maxrun=_max)
    ##         chrono = '_fb{i}'.format(i=_i+1)
    ##         makeResultData(v , _max , lint , specialcut = _rstr+''                                        , scutstring = ''+chrono                   , _options='returnplot,addRares,splitFlavor,printIntegral')
    ##         makeResultData(v , _max , lint , specialcut = _rstr+'&&nll_Edge > 21.'                        , scutstring = 'nllAbove21'+chrono         , _options='returnplot,addRares,splitFlavor,printIntegral')
    ##         makeResultData(v , _max , lint , specialcut = _rstr+'&&nll_Edge > 21. && lepsMll_Edge > 101.' , scutstring = 'highMassnllAbove21'+chrono , _options='returnplot,addRares,splitFlavor,printIntegral')
    ##         makeResultData(v , _max , lint , specialcut = _rstr+'&&nll_Edge > 21. && lepsMll_Edge < 81.'  , scutstring = 'lowMassnllAbove21'+chrono  , _options='returnplot,addRares,splitFlavor,printIntegral')
        
    ##makeClosureTests('nll'    ,'lepsMll_Edge > 101.', 'highMass', True)
    ##makeClosureTests('nllMC'  ,'lepsMll_Edge > 101.', 'highMass', True)
    ##makeClosureTests('nllMCSF','lepsMll_Edge > 101.', 'highMass', True)
    print adsfas

    ## makeResultData('mll', maxrun, lint, specialcut = '' , scutstring = '', returnplot = True)

    ## resultPlot = makeResultData('mll', maxrun, lint, specialcut = '' , scutstring = '', returnplot = True)
    ## ## makeResultData('nll', maxrun = 275125, lint = 4.0)
    ## resultPlotLoNLL = makeResultData('mll',    maxrun,        lint, 'nll_Edge < 21.', 'nllBelow21', returnplot = True)
    ## resultPlotHiNLL = makeResultData('mll',    maxrun,        lint, 'nll_Edge > 21.', 'nllAbove21', returnplot = True)
    ## print addsf
    ## make for region with fixed trigger:
    #### random shit resultPlot = makeResultData('mll', maxrun, lint, specialcut = 'run_Edge >= 274094' , scutstring = 'maxrun274094', returnplot = True)
    #### random shit ##makeResultData('nll', maxrun = 274240, lint = 0.864)
    #### random shit resultPlotLoNLL = makeResultData('mll',    maxrun,        lint, 'nll_Edge < 21. && run_Edge >= 274094', 'nllBelow21_maxrun274094', returnplot = True)
    #### random shit resultPlotHiNLL = makeResultData('mll',    maxrun,        lint, 'nll_Edge > 21. && run_Edge >= 274094', 'nllAbove21_maxrun274094', returnplot = True)
    #####print asdfsadf
    ## makeClosureTests('mll','run_Edge <= 999999', 'withWZ')
    ##makeClosureTests('mll','nll_Edge > 21. && run_Edge <= 999999', 'nllAbove21')
    ##makeClosureTests('mll','nll_Edge > 21. && run_Edge <= 999999 && nBJetMedium25_Edge > 1', 'nllAbove21_nb2')
    #makeClosureTests('nll','lepsMll_Edge > 101. && run_Edge <= 274240', 'highMass_0p8fb-1')

    ##a = makeClosureTests('nll', '', '', True)
    ##makeClosureTests('mll')
    ##makeClosureTests('mll','nll_Edge > 21.', 'nllAbove21')
    ##makeClosureTests('mll','nll_Edge < 21.', 'nllBelow21')
    ## makeClosureTests('nll', '', '', True)
    ##makeClosureTests('nll','lepsMll_Edge < 81.' , 'lowMass')
    ##makeClosureTests('nll','lepsMll_Edge > 101.', 'highMass')
    print asdfasdfsdf

