#####################################################################
######                                                              #
###### 8=========D         OO                        8==========D   #  
###### OO                  ||                                ,88    #
###### ||                  ||                              ,88'     #  
###### ||O---\     ,adPPYb,|| ,adPPYb,d8  ,adPPYba,      ,88'       #
###### ||O---/    a8'    `Y||a8'    `Y88 a8P_____88    ,88'         #
###### ||         8b       ||8b       88 8PP'''''''  ,88'           #
###### \/         '8a,   ,d||'8a,   ,d88 '8b,   ,aa 88'             #
###### 8=========D `'8bbdP'\/  `'YbbdP'Y8  `'Ybbd8'' 8==========D   #
######                       aa,    ,88                             #
######                         'Y8bbdP'                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats, TStyle
import math, sys, optparse, copy, re, array, subprocess


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables


def saveInFile(theFile, measuredValueMC, measuredValueUncMC,  measuredValueData, measuredValueUncData):

    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rsfof") != -1 and  line.find("final") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rsfof        final           DATA        %.4f      %0.4f      0.0000\n'%(measuredValueData, measuredValueUncData))
            else:
                foutput.write('rsfof        final           MC          %.4f      %0.4f      0.0000\n'%(measuredValueMC, measuredValueUncMC))           
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)                                                                                      

                                                                                                                      
def getFinalError(stat, syst):
    err = math.sqrt(stat**2 + syst**2)
    return err                           

def getFactor(x,  xstat, xsyst, y, ystat, ysyst):
    z = x*y
    err = z* math.sqrt((xstat/x)**2 + (xsyst/x)**2+ (ystat/y)**2 + (ysyst/y)**2)
    return z, err                           

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

def makeFactorsTable(): ## for this to make sense the region should be properly binned!!
    rmue_da = helper.readFromFileRmue("ingredients.dat", "DATA") 
    rmue_mc = helper.readFromFileRmue("ingredients.dat", "MC") 
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rsfof_mc = helper.readFromFileRsfofD("ingredients.dat", "MC")   
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rt_mc = helper.readFromFileRT("ingredients.dat", "MC")  
    rsfof_fac_da, rsfof_fac_da_e = getFactor(rmue_da[0], rmue_da[1],rmue_da[2], rt_da[0], rt_da[1], rt_da[2]) 
    rsfof_fac_mc, rsfof_fac_mc_e = getFactor(rmue_mc[0], rmue_mc[1],rmue_mc[2], rt_da[0], rt_da[1], rt_da[2]) 
    line0 = ' & Data & MC '
    line05 = '\\hline '
    line1 = 'rmue  &  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f \\\  ' %(rmue_da[0],getFinalError(rmue_da[1], rmue_da[2]), rmue_mc[0],getFinalError(rmue_mc[1], rmue_mc[2]))
    line2 = 'RT    &  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f  \\\ ' %(rt_da[0],getFinalError(rt_da[1], rt_da[2]), rt_mc[0],getFinalError(rt_mc[1], rt_mc[2]))
    line25 = '\\hline '
    line3 = 'RSFOF factorization &  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f  \\\ ' %(rsfof_fac_da, rsfof_fac_da_e, rsfof_fac_mc, rsfof_fac_mc_e)
    line4 = 'RSFOF direct        &  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f  \\\ ' %(rsfof_da[0],getFinalError(rsfof_da[1], rsfof_da[2]), rsfof_mc[0],getFinalError(rsfof_mc[1], rsfof_mc[2]))
    line5 = 'RSFOF weighted average&  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f  \\\ ' %((rsfof_da[0]+rsfof_fac_da)/2, getFinalError(rsfof_fac_da_e,  getFinalError(rsfof_da[1], rsfof_da[2])), (rsfof_mc[0]+rsfof_fac_mc)/2, getFinalError( rsfof_fac_mc_e, getFinalError(rsfof_mc[1], rsfof_mc[2])))
    line55 = '\\hline'
    print line0                                                                                                                                                      
    print line05                                                                                                                                                      
    print line1                                                                                                                                                      
    print line2                                                                                                                                                      
    print line25                                                                                                                                                      
    print line3                                                                                                                                                      
    print line4                                                                                                                                                      
    print line5                                                                                                                                                      
    print line55 
    saveInFile("ingredients.dat", rsfof_fac_mc, rsfof_fac_mc_e, rsfof_fac_da, rsfof_fac_da_e)
    return (rsfof_da[0]+rsfof_fac_da)/2, getFinalError(rsfof_fac_da_e,  getFinalError(rsfof_da[1], rsfof_da[2])), (rsfof_mc[0]+rsfof_fac_mc)/2, getFinalError(rsfof_fac_mc_e,  getFinalError(rsfof_mc[1], rsfof_mc[2]))



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
    bin86 = h_obsLoNLL.FindBin(86.)
    print "bin 86", bin86, "  ", h_obsLoNLL.GetBinContent(bin86)
    bin96 = h_obsLoNLL.FindBin(96.)
    print "bin 96", bin96, "  ", h_obsLoNLL.GetBinContent(bin96)
    bin150 = h_obsLoNLL.FindBin(150.)
    print "bin 150", bin150, "  ", h_obsLoNLL.GetBinContent(bin150)
    bin200 = h_obsLoNLL.FindBin(200.)                              
    print "bin 200", bin200, "  ", h_obsLoNLL.GetBinContent(bin200)
    bin300 = h_obsLoNLL.FindBin(300.)
    print "bin 300", bin300, "  ", h_obsLoNLL.GetBinContent(bin300)
    bin400 = h_obsLoNLL.FindBin(400.)
    print "bin 400", bin400, "  ", h_obsLoNLL.GetBinContent(bin400)
    bin500 = h_obsLoNLL.GetNbinsX()+1
    print "bin 500", bin500, "  ", h_obsLoNLL.GetBinContent(bin500)

    dyTotal = 98.9; dyTotal_e = 25.1 #12.9 fb-1
    global rinoutbin1, rinoutbin1_e,  rinoutbin2, rinoutbin2_e,  rinoutbin3,rinoutbin3_e,  rinoutbin4, rinoutbin4_e,  rinoutbin5, rinoutbin5_e,  rinoutbin6, rinoutbin6_e 
    rinoutbin1 = readFromFileRinout(ingredientsFile, "DATA", "dy_m1")[0]
    rinoutbin1_stat = readFromFileRinout(ingredientsFile, "DATA", "dy_m1")[1]
    rinoutbin1_syst = readFromFileRinout(ingredientsFile, "DATA", "dy_m1")[2]
    rinoutbin1_e = math.sqrt(rinoutbin1_stat**2 + rinoutbin1_syst**2)           
    rinoutbin2 = readFromFileRinout(ingredientsFile, "DATA", "dy_m2")[0]
    rinoutbin2_stat = readFromFileRinout(ingredientsFile, "DATA", "dy_m2")[1]
    rinoutbin2_syst = readFromFileRinout(ingredientsFile, "DATA", "dy_m2")[2]
    rinoutbin2_e = math.sqrt(rinoutbin2_stat**2 + rinoutbin2_syst**2)           
    rinoutbin3 = readFromFileRinout(ingredientsFile, "DATA", "dy_m3")[0]
    rinoutbin3_stat = readFromFileRinout(ingredientsFile, "DATA", "dy_m3")[1]
    rinoutbin3_syst = readFromFileRinout(ingredientsFile, "DATA", "dy_m3")[2]
    rinoutbin3_e = math.sqrt(rinoutbin3_stat**2 + rinoutbin3_syst**2)           
    rinoutbin4 = readFromFileRinout(ingredientsFile, "DATA", "dy_m4")[0]
    rinoutbin4_stat = readFromFileRinout(ingredientsFile, "DATA", "dy_m4")[1]
    rinoutbin4_syst = readFromFileRinout(ingredientsFile, "DATA", "dy_m4")[2]
    rinoutbin4_e = math.sqrt(rinoutbin4_stat**2 + rinoutbin4_syst**2)           
    rinoutbin5 = readFromFileRinout(ingredientsFile, "DATA", "dy_m5")[0]
    rinoutbin5_stat = readFromFileRinout(ingredientsFile, "DATA", "dy_m5")[1]
    rinoutbin5_syst = readFromFileRinout(ingredientsFile, "DATA", "dy_m5")[2]
    rinoutbin5_e = math.sqrt(rinoutbin5_stat**2 + rinoutbin5_syst**2)           
    rinoutbin6 = readFromFileRinout(ingredientsFile, "DATA", "dy_m6")[0]
    rinoutbin6_stat = readFromFileRinout(ingredientsFile, "DATA", "dy_m6")[1]
    rinoutbin6_syst = readFromFileRinout(ingredientsFile, "DATA", "dy_m6")[2]
    rinoutbin6_e = math.sqrt(rinoutbin6_stat**2 + rinoutbin6_syst**2)           
    dyEffLoNLL = 0.65
    dyEffHiNLL = 1.-dyEffLoNLL
    
    prFSloNbin1_e, prFSloNbin1_e , prFSloNbin2_e, prFSloNbin2_e ,prFSloNbin3_e, prFSloNbin3_e ,prFSloNbin4_e, prFSloNbin4_e, prFSloNbin5_e, prFSloNbin5_e ,prFSloNbin6_e, prFSloNbin6_e = r.Double(),  r.Double(),  r.Double(),      r.Double(), r.Double(),      r.Double(), r.Double(),      r.Double(), r.Double(),      r.Double(), r.Double(),      r.Double()
    prFShiNbin1_e, prFShiNbin1_e , prFShiNbin2_e, prFShiNbin2_e ,prFShiNbin3_e, prFShiNbin3_e ,prFShiNbin4_e, prFShiNbin4_e, prFShiNbin5_e, prFShiNbin5_e ,prFShiNbin6_e, prFShiNbin6_e = r.Double(),  r.Double(),  r.Double(),      r.Double(), r.Double(),      r.Double(), r.Double(),      r.Double(), r.Double(),      r.Double(), r.Double(),      r.Double()

    #print fapppasdfsf
    prFSloNbin1 = h_preLoNLL.IntegralAndError(bin01, bin86, prFSloNbin1_e)
    prDYloNbin1 = dyTotal*rinoutbin1*dyEffLoNLL; prDYloNbin1_e = dyTotal_e*rinoutbin1*dyEffLoNLL
    prTOloNbin1 = prFSloNbin1+prDYloNbin1; prTOloNbin1_e = math.sqrt(prFSloNbin1_e**2 + prDYloNbin1_e**2)
    obTOloNbin1 = h_obsLoNLL.Integral        (bin01, bin86              )
    prFSloNbin2 = h_preLoNLL.IntegralAndError(bin96, bin150, prFSloNbin2_e)
    prDYloNbin2 = dyTotal*rinoutbin2*dyEffLoNLL; prDYloNbin2_e = math.sqrt(prDYloNbin2)
    prTOloNbin2 = prFSloNbin2+prDYloNbin2; prTOloNbin2_e = math.sqrt(prFSloNbin2_e**2 + prDYloNbin2_e**2)
    obTOloNbin2 = h_obsLoNLL.Integral        (bin96, bin150              )                                     
    prFSloNbin3 = h_preLoNLL.IntegralAndError(bin150, bin200, prFSloNbin3_e)
    prDYloNbin3 = dyTotal*rinoutbin3*dyEffLoNLL; prDYloNbin3_e = math.sqrt(prDYloNbin3)
    prTOloNbin3 = prFSloNbin3+prDYloNbin3; prTOloNbin3_e = math.sqrt(prFSloNbin3_e**2 + prDYloNbin3_e**2)
    obTOloNbin3 = h_obsLoNLL.Integral        (bin150, bin200              )                                 
    prFSloNbin4 = h_preLoNLL.IntegralAndError(bin200, bin300, prFSloNbin4_e)
    prDYloNbin4 = dyTotal*rinoutbin4*dyEffLoNLL; prDYloNbin4_e = math.sqrt(prDYloNbin4)
    prTOloNbin4 = prFSloNbin4+prDYloNbin4; prTOloNbin4_e = math.sqrt(prFSloNbin4_e**2 + prDYloNbin4_e**2)
    obTOloNbin4 = h_obsLoNLL.Integral        (bin200, bin300              )                                 
    prFSloNbin5 = h_preLoNLL.IntegralAndError(bin300, bin400, prFSloNbin5_e)
    prDYloNbin5 = dyTotal*rinoutbin5*dyEffLoNLL; prDYloNbin5_e = math.sqrt(prDYloNbin5)
    prTOloNbin5 = prFSloNbin5+prDYloNbin5; prTOloNbin5_e = math.sqrt(prFSloNbin5_e**2 + prDYloNbin5_e**2)
    obTOloNbin5 = h_obsLoNLL.Integral        (bin300, bin400              )                                 
    prFSloNbin6 = h_preLoNLL.IntegralAndError(bin400, bin500, prFSloNbin6_e)                                  
    prDYloNbin6 = dyTotal*rinoutbin6*dyEffLoNLL; prDYloNbin6_e = math.sqrt(prDYloNbin6)
    prTOloNbin6 = prFSloNbin6+prDYloNbin6; prTOloNbin6_e = math.sqrt(prFSloNbin6_e**2 + prDYloNbin6_e**2)
    obTOloNbin6 = h_obsLoNLL.Integral        (bin400, bin500              )                                 

    prFShiNbin1 = h_preHiNLL.IntegralAndError(bin01, bin86, prFShiNbin1_e)
    prDYhiNbin1 = dyTotal*rinoutbin1*dyEffHiNLL; prDYhiNbin1_e = math.sqrt(prDYhiNbin1)
    prTOhiNbin1 = prFShiNbin1+prDYhiNbin1; prTOhiNbin1_e = math.sqrt(prFShiNbin1_e**2 + prDYhiNbin1_e**2)
    obTOhiNbin1 = h_obsHiNLL.Integral        (bin01, bin86              )
    prFShiNbin2 = h_preHiNLL.IntegralAndError(bin96, bin150, prFShiNbin2_e)
    prDYhiNbin2 = dyTotal*rinoutbin2*dyEffHiNLL; prDYhiNbin2_e = math.sqrt(prDYhiNbin2)
    prTOhiNbin2 = prFShiNbin2+prDYhiNbin2; prTOhiNbin2_e = math.sqrt(prFShiNbin2_e**2 + prDYhiNbin2_e**2)
    obTOhiNbin2 = h_obsHiNLL.Integral        (bin96, bin150              )                                 
    prFShiNbin3 = h_preHiNLL.IntegralAndError(bin150, bin200, prFShiNbin3_e)                                  
    prDYhiNbin3 = dyTotal*rinoutbin3*dyEffHiNLL; prDYhiNbin3_e = math.sqrt(prDYhiNbin3)
    prTOhiNbin3 = prFShiNbin3+prDYhiNbin3; prTOhiNbin3_e = math.sqrt(prFShiNbin3_e**2 + prDYhiNbin3_e**2)
    obTOhiNbin3 = h_obsHiNLL.Integral        (bin150, bin200              )                                 
    prFShiNbin4 = h_preHiNLL.IntegralAndError(bin200, bin300, prFShiNbin4_e)
    prDYhiNbin4 = dyTotal*rinoutbin4*dyEffHiNLL; prDYhiNbin4_e = math.sqrt(prDYhiNbin4)
    prTOhiNbin4 = prFShiNbin4+prDYhiNbin4; prTOhiNbin4_e = math.sqrt(prFShiNbin4_e**2 + prDYhiNbin4_e**2)
    obTOhiNbin4 = h_obsHiNLL.Integral        (bin200, bin300              )                                 
    prFShiNbin5 = h_preHiNLL.IntegralAndError(bin300, bin400, prFShiNbin5_e)
    prDYhiNbin5 = dyTotal*rinoutbin5*dyEffHiNLL; prDYhiNbin5_e = math.sqrt(prDYhiNbin5)
    prTOhiNbin5 = prFShiNbin5+prDYhiNbin5; prTOhiNbin5_e = math.sqrt(prFShiNbin5_e**2 + prDYhiNbin5_e**2)
    obTOhiNbin5 = h_obsHiNLL.Integral        (bin300, bin400              )                                 
    prFShiNbin6 = h_preHiNLL.IntegralAndError(bin400, bin500, prFShiNbin6_e)
    prDYhiNbin6 = dyTotal*rinoutbin6*dyEffHiNLL; prDYhiNbin6_e = math.sqrt(prDYhiNbin6)
    prTOhiNbin6 = prFShiNbin6+prDYhiNbin6; prTOhiNbin6_e = math.sqrt(prFShiNbin6_e**2 + prDYhiNbin6_e**2)
    obTOhiNbin6 = h_obsHiNLL.Integral        (bin400, bin500              )                                 

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
\\multirow{{4}}{{*}}{{ 20 $<$ mll $<$ 86 GeV}}       & pred. FS        & {prFSloNbin1:.1f}  $\\pm$  {prFSloNbin1_e:.1f}    & {prFShiNbin1:.1f}     $\\pm$  {prFShiNbin1_e:.1f}  \\\\
                                             & pred. DY        & {prDYloNbin1:.1f}  $\\pm$  {prDYloNbin1_e:.1f}    & {prDYhiNbin1:.1f}     $\\pm$  {prDYhiNbin1_e:.1f}  \\\\
                                             & pred. total     & {prTOloNbin1:.1f}  $\\pm$  {prTOloNbin1_e:.1f}    & {prTOhiNbin1:.1f}     $\\pm$  {prTOhiNbin1_e:.1f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNbin1}}}                        & \\textbf{{{obTOhiNbin1}}}                         \\\\ \\hline
\\multirow{{4}}{{*}}{{ 96 $<$ mll $<$ 150 GeV}}      & pred. FS     & {prFSloNbin2:.1f}  $\\pm$  {prFSloNbin2_e:.1f}       & {prFShiNbin2:.1f}     $\\pm$  {prFShiNbin2_e:.1f}  \\\\
                                             & pred. DY     & {prDYloNbin2:.1f}  $\\pm$  {prDYloNbin2_e:.1f}       & {prDYhiNbin2:.1f}     $\\pm$  {prDYhiNbin2_e:.1f}  \\\\
                                             & pred. total  & {prTOloNbin2:.1f}  $\\pm$  {prTOloNbin2_e:.1f}       & {prTOhiNbin2:.1f}     $\\pm$  {prTOhiNbin2_e:.1f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNbin2}}}                        & \\textbf{{{obTOhiNbin2}}}                         \\\\ \\hline                

\\multirow{{4}}{{*}}{{ 150 $<$ mll $<$ 200 GeV}}      & pred. FS     & {prFSloNbin3:.1f}  $\\pm$  {prFSloNbin3_e:.1f}       & {prFShiNbin3:.1f}     $\\pm$  {prFShiNbin3_e:.1f}  \\\\
                                             & pred. DY     & {prDYloNbin3:.1f}  $\\pm$  {prDYloNbin3_e:.1f}       & {prDYhiNbin3:.1f}     $\\pm$  {prDYhiNbin3_e:.1f}  \\\\
                                             & pred. total  & {prTOloNbin3:.1f}  $\\pm$  {prTOloNbin3_e:.1f}       & {prTOhiNbin3:.1f}     $\\pm$  {prTOhiNbin3_e:.1f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNbin3}}}                        & \\textbf{{{obTOhiNbin3}}}                         \\\\ \\hline 
\\multirow{{4}}{{*}}{{ 200 $<$ mll $<$ 300 GeV}}      & pred. FS     & {prFSloNbin4:.1f}  $\\pm$  {prFSloNbin4_e:.1f}       & {prFShiNbin4:.1f}     $\\pm$  {prFShiNbin4_e:.1f}  \\\\
                                             & pred. DY     & {prDYloNbin4:.1f}  $\\pm$  {prDYloNbin4_e:.1f}       & {prDYhiNbin4:.1f}     $\\pm$  {prDYhiNbin4_e:.1f}  \\\\
                                             & pred. total  & {prTOloNbin4:.1f}  $\\pm$  {prTOloNbin4_e:.1f}       & {prTOhiNbin4:.1f}     $\\pm$  {prTOhiNbin4_e:.1f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNbin4}}}                        & \\textbf{{{obTOhiNbin4}}}                         \\\\ \\hline 
\\multirow{{4}}{{*}}{{ 300 $<$ mll $<$ 400 GeV}}      & pred. FS     & {prFSloNbin5:.1f}  $\\pm$  {prFSloNbin5_e:.1f}       & {prFShiNbin5:.1f}     $\\pm$  {prFShiNbin5_e:.1f}  \\\\
                                             & pred. DY     & {prDYloNbin5:.1f}  $\\pm$  {prDYloNbin5_e:.1f}       & {prDYhiNbin5:.1f}     $\\pm$  {prDYhiNbin5_e:.1f}  \\\\
                                             & pred. total  & {prTOloNbin5:.1f}  $\\pm$  {prTOloNbin5_e:.1f}       & {prTOhiNbin5:.1f}     $\\pm$  {prTOhiNbin5_e:.1f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNbin5}}}                        & \\textbf{{{obTOhiNbin5}}}                         \\\\ \\hline 
\\multirow{{4}}{{*}}{{ mll $>$ 400 GeV}}      & pred. FS     & {prFSloNbin6:.1f}  $\\pm$  {prFSloNbin6_e:.1f}       & {prFShiNbin6:.1f}     $\\pm$  {prFShiNbin6_e:.1f}  \\\\
                                             & pred. DY     & {prDYloNbin6:.1f}  $\\pm$  {prDYloNbin6_e:.1f}       & {prDYhiNbin6:.1f}     $\\pm$  {prDYhiNbin6_e:.1f}  \\\\
                                             & pred. total  & {prTOloNbin6:.1f}  $\\pm$  {prTOloNbin6_e:.1f}       & {prTOhiNbin6:.1f}     $\\pm$  {prTOhiNbin6_e:.1f}  \\\\
                                             & \\textbf{{obs}} & \\textbf{{{obTOloNbin6}}}                        & \\textbf{{{obTOhiNbin6}}}                         \\\\ 

\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(
    prFSloNbin1_e=prFSloNbin1_e, prDYloNbin1_e=prDYloNbin1_e,   prTOloNbin1_e=prTOloNbin1_e, prFSloNbin2_e=prFSloNbin2_e, prDYloNbin2_e=prDYloNbin2_e, prTOloNbin2_e=prTOloNbin2_e, prFSloNbin3_e=prFSloNbin3_e, prDYloNbin3_e=prDYloNbin3_e, prTOloNbin3_e=prTOloNbin3_e,prFSloNbin4_e=prFSloNbin4_e, prDYloNbin4_e=prDYloNbin4_e, prTOloNbin4_e=prTOloNbin4_e,prFSloNbin5_e=prFSloNbin5_e, prDYloNbin5_e=prDYloNbin5_e, prTOloNbin5_e=prTOloNbin5_e,prFSloNbin6_e=prFSloNbin6_e, prDYloNbin6_e=prDYloNbin6_e, prTOloNbin6_e=prTOloNbin6_e,
    
    prFShiNbin1_e=prFShiNbin1_e, prDYhiNbin1_e=prDYhiNbin1_e,   prTOhiNbin1_e=prTOhiNbin1_e, prFShiNbin2_e=prFShiNbin2_e, prDYhiNbin2_e=prDYhiNbin2_e, prTOhiNbin2_e=prTOhiNbin2_e, prFShiNbin3_e=prFShiNbin3_e, prDYhiNbin3_e=prDYhiNbin3_e, prTOhiNbin3_e=prTOhiNbin3_e,prFShiNbin4_e=prFShiNbin4_e, prDYhiNbin4_e=prDYhiNbin4_e, prTOhiNbin4_e=prTOhiNbin4_e,prFShiNbin5_e=prFShiNbin5_e, prDYhiNbin5_e=prDYhiNbin5_e, prTOhiNbin5_e=prTOhiNbin5_e,prFShiNbin6_e=prFShiNbin6_e, prDYhiNbin6_e=prDYhiNbin6_e, prTOhiNbin6_e=prTOhiNbin6_e,
    
    prFSloNbin1 = prFSloNbin1,      prFShiNbin1 = prFShiNbin1,
    prDYloNbin1 = prDYloNbin1,      prDYhiNbin1 = prDYhiNbin1,
    prTOloNbin1 = prTOloNbin1,      prTOhiNbin1 = prTOhiNbin1,
    obTOloNbin1 = int(obTOloNbin1), obTOhiNbin1 = int(obTOhiNbin1),
    prFSloNbin2 = prFSloNbin2,      prFShiNbin2 = prFShiNbin2,
    prDYloNbin2 = prDYloNbin2,      prDYhiNbin2 = prDYhiNbin2,
    prTOloNbin2 = prTOloNbin2,      prTOhiNbin2 = prTOhiNbin2,
    obTOloNbin2 = int(obTOloNbin2), obTOhiNbin2 = int(obTOhiNbin2),
    prFSloNbin3 = prFSloNbin3,      prFShiNbin3 = prFShiNbin3,
    prDYloNbin3 = prDYloNbin3,      prDYhiNbin3 = prDYhiNbin3,
    prTOloNbin3 = prTOloNbin3,      prTOhiNbin3 = prTOhiNbin3,
    obTOloNbin3 = int(obTOloNbin3), obTOhiNbin3 = int(obTOhiNbin3),
    prFSloNbin4 = prFSloNbin4,      prFShiNbin4 = prFShiNbin4,
    prDYloNbin4 = prDYloNbin4,      prDYhiNbin4 = prDYhiNbin4,
    prTOloNbin4 = prTOloNbin4,      prTOhiNbin4 = prTOhiNbin4,
    obTOloNbin4 = int(obTOloNbin4), obTOhiNbin4 = int(obTOhiNbin4),
    prFSloNbin5 = prFSloNbin5,      prFShiNbin5 = prFShiNbin5,
    prDYloNbin5 = prDYloNbin5,      prDYhiNbin5 = prDYhiNbin5,
    prTOloNbin5 = prTOloNbin5,      prTOhiNbin5 = prTOhiNbin5,
    obTOloNbin5 = int(obTOloNbin5), obTOhiNbin5 = int(obTOhiNbin5),
    prFSloNbin6 = prFSloNbin6,      prFShiNbin6 = prFShiNbin6,
    prDYloNbin6 = prDYloNbin6,      prDYhiNbin6 = prDYhiNbin6,
    prTOloNbin6 = prTOloNbin6,      prTOhiNbin6 = prTOhiNbin6,
    obTOloNbin6 = int(obTOloNbin6), obTOhiNbin6 = int(obTOhiNbin6), lint = lint)

    print "made the table!!!, its here plots/results/%s/tables/"%(lint_str)
    helper.ensureDirectory('plots/results/%s/'%lint_str); 
    helper.ensureDirectory('plots/results/%s/tables/'%lint_str)
    tablename = (resultPlotLoNLL.name.split('_')[-1]).replace('nllABove21','').replace('nllBelow21','')
    compTableFile = open('plots/results/%s/tables/resultTable_%s%s.tex'%(lint_str, str(lint).replace('.','p'), tablename),'w')
    compTableFileWWW = open('/afs/cern.ch/user/m/mvesterb/www/edgeZ/resultsTable/resultTable_%s%s.tex'%(str(lint).replace('.','p'), tablename),'w')
    compTableFile.write(resultTable)
    compTableFileWWW.write(resultTable)
    compTableFile.close()
    compTableFileWWW.close()

def makeClosureTests(var, specialcut = '', scutstring = '', doCumulative = False):

    if var == 'mll':
        treevar = 'lepsMll_Edge'
        nbins, xmin, xmax = 28, 20, 300
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
    mc_OF = treeMC.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.Zveto, cuts.OF]), '', xlabel)
    mc_SF = treeMC.getTH1F(lint, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.Zveto, cuts.SF]), '', xlabel)
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.Zveto, cuts.OF]), '', xlabel)
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.Zveto, cuts.SF]), '', xlabel)
    #mc_ee = treeMC.getTH1F(lint, var+"mc_ee"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.Zveto, cuts.ee]), '', xlabel)
    #mc_mm = treeMC.getTH1F(lint, var+"mc_mm"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.Zveto, cuts.mm]), '', xlabel)
    dy_SF = treeDY.getTH1F(lint, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.Zveto, cuts.SF]), '', xlabel)
    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_err.SetFillStyle(3004); mc_OF_err.SetMarkerSize(0.)
    dy_SF.SetFillColorAlpha(r.kGreen+2,0.5)

    mc_SF.GetYaxis().SetRangeUser(0., 1.3*mc_SF.GetMaximum())
    print helper.bcolors.HEADER + '[MC only closure test not scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure_noRSFOF = Canvas.Canvas('closure/%s/plot_closure_%s_mcPredmcObs%s_noRSFOF'%(lint_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_closure_noRSFOF.addHisto(mc_SF    , 'PE'       , 'MC - SF', 'PL', r.kRed+1  , 1,  0)
    plot_closure_noRSFOF.addHisto(mc_OF_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure_noRSFOF.addHisto(mc_OF    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_closure_noRSFOF.addHisto(dy_SF    , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
    plot_closure_noRSFOF.addLatex (0.61, 0.82, 'no R_{SFOF} scaling')
    plot_closure_noRSFOF.saveRatio(1, 0, 1, lint, mc_SF, mc_OF, 0.2, 1.8)

    mc_OF_rsfofScaled = copy.deepcopy(mc_OF)
    mc_OF_rsfofScaled = scaleByRSFOF(mc_OF_rsfofScaled, rsfof_mc, rsfof_mc_e)
    mc_OF_rsfofScaled_err = copy.deepcopy(mc_OF_rsfofScaled)
    mc_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_rsfofScaled_err.SetFillStyle(3004); mc_OF_rsfofScaled_err.SetMarkerSize(0.)

    print helper.bcolors.HEADER + '[MC only closure test scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure = Canvas.Canvas('closure/%s/plot_closure_%s_mcPredmcObs%s'%(lint_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_closure.addHisto(mc_SF                , 'PE'       , 'MC - SF', 'PL', r.kRed+1  , 1,  0)
    plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_closure.addHisto(dy_SF                , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
    plot_closure.addLatex (0.61, 0.82, 'R_{SFOF} scaled')
    plot_closure.saveRatio(1, 0, 0, lint, mc_SF, mc_OF_rsfofScaled, 0.2, 1.8)

    if False:
        plot_closure = Canvas.Canvas('closure/%s/plot_closure_%s_mcPreddaObs%s'%(lint_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
        plot_closure.addHisto(da_SF                , 'PE'       , 'data-SF', 'PL', r.kRed+1  , 1,  0)
        plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
        plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
        plot_closure.addHisto(dy_SF                , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
        plot_closure.addLatex(0.61, 0.82, 'R_{SFOF} scaled')
        plot_closure.saveRatio(1, 0, 0, lint , mc_SF, mc_OF_rsfofScaled, 0.2, 1.8)

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
        plot_cumulative = Canvas.Canvas('closure/%s/plot_cumulative_%s_mcPredmcObs%s'%(lint_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.8, 0.2, 0.95, 0.4)
        plot_cumulative.addHisto(mc_SF_cum                 , 'hist,same', 'MC - SF', 'L' , r.kRed+1  , 1,  0)
        plot_cumulative.addHisto(mc_OF_rsfofScaled_cum     , 'hist,SAME', 'MC - OF', 'L', r.kBlue+1 , 1,  1)
        #plot_cumulative.addHisto(dy_SF_cum                 , 'hist,SAME', 'DY - SF', 'L', r.kGreen+2, 1,  2)
        plot_cumulative.addLatex (0.61, 0.82, 'R_{SFOF} scaled')
        plot_cumulative.saveRatio(1, 0, 0, lint, mc_SF_cum, mc_OF_rsfofScaled_cum, 0.2, 1.8)
        return mc_SF_cum

def makeResultData(var, maxrun = 999999, lint = 12.9, specialcut = '', scutstring = '', _options = ''):
    returnplot, addRares, splitFlavor, makeTable, printIntegral = False, True, False, False, False
    if 'returnplot'    in _options: print 'found option %s'%'returnplot'    ;returnplot    = True
    addRares      = True
    if 'splitFlavor'   in _options: print 'found option %s'%'splitFlavor'   ;splitFlavor   = True
    if 'makeTable'     in _options: print 'found option %s'%'makeTable'     ;makeTable     = True
    if 'printIntegral' in _options: print 'found option %s'%'printIntegral' ;printIntegral = True

    # if   var == 'mll'      : treevar = 'lepsMll_Edge'        ; nbins, xmin, xmax = 23, 20 , 250 ; xlabel = 'm_{ll} (GeV)'
    if   var == 'mll'      : treevar = 'lepsMll_Edge'        ; nbins, xmin, xmax = 30, 20 , 420 ; xlabel = 'm_{ll} [GeV]'
    elif var == 'nll'      : treevar = 'nll_Edge'            ; nbins, xmin, xmax = 15, 12 , 27  ; xlabel = 'NLL'
    elif var == 'nb'       : treevar = 'nBJetMedium35_Edge'  ; nbins, xmin, xmax =  3,  0 ,  3  ; xlabel = 'n_{b-jets,35}'
    elif var == 'nj'       : treevar = 'nJetSel_Edge'        ; nbins, xmin, xmax =  6, 0.5, 6.5 ; xlabel = 'n_{jets}'
    elif var == 'zpt'      : treevar = 'lepsZPt_Edge'        ; nbins, xmin, xmax = 10,  0 ,1000 ; xlabel = 'p_{T}^{ll}'
    elif var == 'mlb'      : treevar = 'sum_mlb_Edge'        ; nbins, xmin, xmax = 15,  0 ,1500 ; xlabel = '#Sigma m_{lb}'
    elif var == 'met'      : treevar = 'met_Edge'            ; nbins, xmin, xmax = 10,100 ,1000 ; xlabel = 'E_{T}^{miss.} [GeV]'
    elif var == 'metraw'   : treevar = 'met_raw_Edge'        ; nbins, xmin, xmax = 10,100 ,1000 ; xlabel = 'E_{T}^{miss.} raw'
    elif var == 'ldp'      : treevar = 'abs(lepsDPhi_Edge)'  ; nbins, xmin, xmax = 10,  0 , 3.15; xlabel = '#Delta #phi ll'
    elif var == 'pt1'      : treevar = 'Lep1_pt_Edge'        ; nbins, xmin, xmax = 20,  0 , 500 ; xlabel = 'p_{T} leading'
    elif var == 'pt2'      : treevar = 'Lep2_pt_Edge'        ; nbins, xmin, xmax = 10, 20 , 200 ; xlabel = 'p_{T} trailing'
    elif var == 'ldr'      : treevar = 'lepsDR_Edge'         ; nbins, xmin, xmax = 10,  0 , 6.  ; xlabel = '#Delta R (ll)'
    elif var == 'iso1'     : treevar = 'Lep1_miniRelIso_Edge'; nbins, xmin, xmax = 10,  0 , 0.05; xlabel = 'mini iso l1'
        

    newLumiString = str(lint)+'invfb'
    if not specialcut:
        specialcut = 'run_Edge <= {run}'.format(run=maxrun)
    else:
        specialcut += ' && run_Edge <= {run}'.format(run=maxrun)
    ## ## mll distributions
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.SF, cuts.Zveto]), '', xlabel)
    da_mm = treeDA.getTH1F(lint, var+"da_mm"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.mm, cuts.Zveto]), '', xlabel)
    da_ee = treeDA.getTH1F(lint, var+"da_ee"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.ee, cuts.Zveto]), '', xlabel)
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.OF, cuts.Zveto]), '', xlabel)
    if addRares:
        print "Adding the rares"
        ra_OF = treeRA.getTH1F(lint, var+"ra_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.OF, cuts.Zveto]), '', xlabel)
        ra_SF = treeRA.getTH1F(lint, var+"ra_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.SF, cuts.Zveto]), '', xlabel)
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
    sstring = '' if not addRares else ''
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
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-l', '--loadShapes', action='store_true', dest='loadShapes', default=False, help='reload dy shapes. default is off since this takes a while')
    parser.add_option('-M', '--maxRun', action='store', type=int, dest='maxRun', default=999999, help='max run to use for analysis (run is included)')
    parser.add_option('-m', '--minRun', action='store', type=int, dest='minRun', default=-1    , help='min run to use for analysis (run not included)')
    ## make the options globa.. also the lumi
    global opts, lumi, lumi_str, dy_shapes, nbstring
    global ingMC, ingDA, onZ, treeDA, treeMC, treeDY, treeTT
    (opts, args) = parser.parse_args()

    ingredientsFile = opts.ingredientsFile
    print 'running with these options \n'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'

    ## rsfofTable = Tables.makeRSFOFTable(ingDA, ingMC)

    ##print asdf
    print 'Going to load DATA and MC trees...'
    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_HT100to200', 'DYJetsToLL_M50_HT200to400', 'DYJetsToLL_M50_HT400to600', 'DYJetsToLL_M50_HT600toInf', 'DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext', 'DYJetsToLL_M50_HT400to600_ext', 'DYJetsToLL_M50_HT600toInf_ext']
    raDatasets = ['ZZTo4L', 'WZTo3LNu', 'WWW', 'WWZ','WZZ', 'ZZZ',  'TTZToLLNuNu' ,'TTWToLNu', 'T_tWch', 'TBar_tWch' , 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TToLeptons_sch', 'TToLeptons_tch_powheg', 'TBarToLeptons_tch_powheg', 'TTHnobb_pow', 'VHToNonbb', 'WJetsToLNu_LO']
    ttDatasets = ['TTJets_DiLepton', 'TTJets_DiLepton_ext']
    mcDatasets = ttDatasets + ([] if opts.onlyTTbar else dyDatasets)
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_273150_275376', 'DoubleEG_Run2016B-PromptReco-v2_runs_273150_275376', 'MuonEG_Run2016B-PromptReco-v2_runs_273150_275376',
                  'DoubleMuon_Run2016C-PromptReco-v2_runs_275420_276283', 'DoubleEG_Run2016C-PromptReco-v2_runs_275420_276283', 'MuonEG_Run2016C-PromptReco-v2_runs_275420_276283',
                  'DoubleMuon_Run2016D-PromptReco-v2_runs_276315_276811', 'DoubleEG_Run2016D-PromptReco-v2_runs_276315_276811', 'MuonEG_Run2016D-PromptReco-v2_runs_276315_276811']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0)
    treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT'  , 0)
    treeRA = Sample.Tree(helper.selectSamples(opts.sampleFile, raDatasets, 'RA'), 'RA'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'

    isBlinded = True

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()
    lint = 12.9  ; maxrun = 999999; lint_str = '12.9invfb'
    print 'Running with an integrated luminosity of %.2f fb-1' %(lint)

    ## ============================================================
    ## ========== set RSFOF globally ==============================
    ## ============================================================
    global rsfof, rsfof_e, rsfof_mc, rsfof_mc_e
    rsfof, rsfof_e, rsfof_mc, rsfof_mc_e =  makeFactorsTable()
    print rsfof, rsfof_e, rsfof_mc, rsfof_mc_e
    ## result plots in different variables:
 
 #detta har du gjort leo!!
    #resultPlotLoNLL = makeResultData('nll' , maxrun , lint , specialcut = ''               , scutstring = 'mllInc'    , _options='returnplot')
    #resultPlotHiNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge > 21.' , scutstring = 'nllAbove21', _options='returnplot')
    #resultPlotLoNLL = makeResultData('mll' , maxrun , lint , specialcut = 'nll_Edge < 21.' , scutstring = 'nllBelow21', _options='returnplot')
    #makeResultTable(resultPlotLoNLL, resultPlotHiNLL, lint, lint_str)

 #   #for v in ['mll']:#'nll_noMET', 'nll_noMLB', 'nll_noZPT', 'nll_noLDP']: #'iso1', 'iso2', 'mll', 'nll', 'nb', 'nj', 'zpt', 'mlb', 'met', 'ldp', 'pt1', 'pt2']:
 #   #    makeResultData(v , maxrun , lint , specialcut =''                                     , scutstring = ''                 , _options='returnplot,splitFlavor,printIntegral')
 #   #    makeResultData(v , maxrun , lint , specialcut ='nll_Edge > 21.'                        , scutstring = 'nllAbove21'        , _options='returnplot,splitFlavor,printIntegral')
 #   #    makeResultData(v , maxrun , lint , specialcut ='nll_Edge <= 21.'                        , scutstring = 'nllBelow21'        , _options='returnplot,splitFlavor,printIntegral')
 #       
 #   for v in ['met']:#'nll_noMET', 'nll_noMLB', 'nll_noZPT', 'nll_noLDP']: #'iso1', 'iso2', 'mll', 'nll', 'nb', 'nj', 'zpt', 'mlb', 'met', 'ldp', 'pt1', 'pt2']:
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge > 21. && lepsMll_Edge <= 86. && lepsMll_Edge > 20' , scutstring = 'mll20-86_nllAbove21', _options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge <= 21. && lepsMll_Edge <= 86. && lepsMll_Edge > 20' , scutstring = 'mll20-86_nllBelow21', _options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge > 21. && lepsMll_Edge <= 150. && lepsMll_Edge > 96' , scutstring = 'mll96-150_nllAbove21', _options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge <= 21. && lepsMll_Edge <= 150. && lepsMll_Edge > 96' , scutstring = 'mll96-150_nllBelow21', _options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge > 21. && lepsMll_Edge <= 200. && lepsMll_Edge > 150' , scutstring = 'mll150-200_nllAbove21',_options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge <= 21. && lepsMll_Edge<= 200. && lepsMll_Edge > 150' , scutstring = 'mll150-200_nllBelow21',_options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge > 21. && lepsMll_Edge <= 300. && lepsMll_Edge > 200' , scutstring ='mll200-300_nllAbove21',_options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge <= 21. && lepsMll_Edge <= 300.&& lepsMll_Edge > 200' , scutstring ='mll200-300_nllBelow21', _options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge > 21. && lepsMll_Edge <= 400. && lepsMll_Edge > 300' , scutstring ='mll300-400_nllAbove21',_options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge <= 21. && lepsMll_Edge <= 400. && lepsMll_Edge > 300' , scutstring ='mll300-400_nllBelow21',_options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge > 21. && lepsMll_Edge > 400' , scutstring = 'mll400_nllAbove21', _options='returnplot,splitFlavor')
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge <= 21. && lepsMll_Edge > 400' , scutstring = 'mll400_nllBelow21', _options='returnplot,splitFlavor')


    makeClosureTests('nll','', 'inclusive', True)
    makeClosureTests('nll','lepsMll_Edge <= 86. && lepsMll_Edge > 20', 'mll20-86', True)
    makeClosureTests('nll','lepsMll_Edge <= 150. && lepsMll_Edge > 96', 'mll96-150', True)
    makeClosureTests('nll','lepsMll_Edge <= 200. && lepsMll_Edge > 150', 'mll150-200', True)
    makeClosureTests('nll','lepsMll_Edge <= 300. && lepsMll_Edge > 200', 'mll200-300', True)
    makeClosureTests('nll','lepsMll_Edge <= 400. && lepsMll_Edge > 300', 'mll300-400', True)
    makeClosureTests('nll','lepsMll_Edge > 400', 'mll400', True)
    makeClosureTests('mll','', 'inclusive', True)
    makeClosureTests('mll','nll_Edge > 21.' , 'nllAbove21', True)
    makeClosureTests('mll','nll_Edge <= 21.' , 'nllBelow21', True)
    
    #makeResultData('mll', maxrun, lint, specialcut = '' , scutstring = '', returnplot = True)

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

