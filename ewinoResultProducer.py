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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats, TStyle, TH1,   SetOwnership
import math, sys, optparse, copy, re, array, subprocess
from array import array

import include.LeptonSF
import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables
import include.Scans      as Scans


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

def saveInFileFactors(theFile, measuredValueMC, measuredValueUncMC,  measuredValueData, measuredValueUncData, var, reg):

    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find(var) != -1 and  line.find(reg) != -1 and  line.find(reg) != -1:
            if line.find("DATA") != -1:
                foutput.write('%s       %s    DATA        %.4f      %0.4f      0.0000\n'%(var, reg, measuredValueData, measuredValueUncData))
            else:
                foutput.write('%s       %s    MC          %.4f      %0.4f      0.0000\n'%(var, reg, measuredValueMC, measuredValueUncMC))           
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)


def makeClosureTable(obs_tt, pred_tt, pred_da):
    obs_tt_lowMet_e = r.Double(); pred_tt_lowMet_e = r.Double();pred_da_lowMet_e = r.Double();obs_tt_highMet_e = r.Double(); pred_tt_highMet_e = r.Double();pred_da_highMet_e = r.Double();
    
    obs_tt_lowMet = obs_tt.IntegralAndError(obs_tt.FindBin(100.), obs_tt.FindBin(300.), obs_tt_lowMet_e)
    pred_tt_lowMet = pred_tt.IntegralAndError(pred_tt.FindBin(100.), pred_tt.FindBin(300.), pred_tt_lowMet_e)
    pred_da_lowMet = pred_da.IntegralAndError(pred_da.FindBin(100.), pred_da.FindBin(300.), pred_da_lowMet_e)
    obs_tt_highMet = obs_tt.IntegralAndError(obs_tt.FindBin(150.), obs_tt.FindBin(300.), obs_tt_highMet_e)
    pred_tt_highMet = pred_tt.IntegralAndError(pred_tt.FindBin(150.), pred_tt.FindBin(300.), pred_tt_highMet_e)
    pred_da_highMet = pred_da.IntegralAndError(pred_da.FindBin(150.), pred_da.FindBin(300.), pred_da_highMet_e)
    
    line0 = '          &  MET 100 GeV & MET 150 GeV \\\\\hline'
    line1 = ' obs(ttbar)  &  %.3f $\\pm$ %.3f   &   %.3f $\\pm$ %.3f \\\\ ' %(obs_tt_lowMet  , obs_tt_lowMet_e , obs_tt_highMet , obs_tt_highMet_e  )
    line2 = ' pred(ttbar) &  %.3f $\\pm$ %.3f   &   %.3f $\\pm$ %.3f \\\\ ' %(pred_tt_lowMet  , pred_tt_lowMet_e , pred_tt_highMet  , pred_tt_highMet_e  )
    line3 = ' pred(data)  &  %.3f $\\pm$ %.3f   &   %.3f $\\pm$ %.3f \\\\\hline ' %(pred_da_lowMet, pred_da_lowMet_e  , pred_da_highMet, pred_da_highMet_e)
    line4 = ' obs./pred (ttbar)  &  %.3f $\\pm$ %.3f   &   %.3f $\\pm$ %.3f \\\\ ' %(getFraction(obs_tt_lowMet, obs_tt_lowMet_e, pred_tt_lowMet, pred_tt_lowMet_e)[0], getFraction(obs_tt_lowMet, obs_tt_lowMet_e, pred_tt_lowMet, pred_tt_lowMet_e)[1], getFraction(obs_tt_highMet, obs_tt_highMet_e, pred_tt_highMet, pred_tt_highMet_e)[0], getFraction(obs_tt_highMet, obs_tt_highMet_e, pred_tt_highMet, pred_tt_highMet_e)[1])
    line5 = ' obs./pred (data)  &  %.3f $\\pm$ %.3f   &   %.3f $\\pm$ %.3f \\\\ ' %(getFraction(obs_tt_lowMet, obs_tt_lowMet_e, pred_da_lowMet, pred_da_lowMet_e)[0],getFraction(obs_tt_lowMet, obs_tt_lowMet_e, pred_da_lowMet, pred_da_lowMet_e)[1], getFraction(obs_tt_highMet, obs_tt_highMet_e, pred_da_highMet, pred_da_highMet_e)[0], getFraction(obs_tt_highMet, obs_tt_highMet_e, pred_da_highMet, pred_da_highMet_e)[0])
    line0 += '\\hline';
    print line0
    print line1
    print line2
    print line3                                                                                                             
    print line4                                                                                                             
    print line5                                                                                                             


def makeEWKFactorsTable(h1_mc, h1_mc_e, h2_mc, h2_mc_e, h3_mc, h3_mc_e ,h1_da, h1_da_e, h2_da, h2_da_e, h3_da, h3_da_e , var ): 
    if var == 'fmll':
        line0 = ' region         &  MC & Data \\\\\hline'
        line1 = ' OF $\geq$ 1 b  &  %.3f $\\pm$ %.3f   &   %.3f $\\pm$ %.3f \\\\ ' %(h1_mc  , h1_mc_e , h1_da  , h1_da_e  )
        line2 = ' OF bveto       &  %.3f $\\pm$ %.3f   &   %.3f $\\pm$ %.3f \\\\ ' %(h2_mc  , h2_mc_e , h2_da  , h2_da_e  )
        line3 = ' SF bveto       &  %.3f $\\pm$ %.3f   & blinded   \\\\ ' %(h3_mc  , h3_mc_e )
        line0 += '\\hline';
        print "fmll"
        print line0
        print line1
        print line2
        print line3                                                                                                             
    else:
        line0 = ' region         &  MC & Data \\\\\hline'
        line1 = ' OF onZ         &  %.3f $\\pm$ %.3f   &   %.3f $\\pm$ %.3f \\\\ ' %(h1_mc  , h1_mc_e , h1_da  , h1_da_e  )
        line2 = ' OF mll 61-121  &  %.3f $\\pm$ %.3f   &   %.3f $\\pm$ %.3f \\\\ ' %(h2_mc  , h2_mc_e , h2_da  , h2_da_e  )
        line3 = ' SF mll 61-121  &  %.3f $\\pm$ %.3f   & blinded  \\\\  ' %(h3_mc  , h3_mc_e )
        line0 += '\\hline'; 
        print "r0b1b"
        print line0
        print line1
        print line2
        print line3                                                                                                                                                                                   


def getFinalError(stat, syst):
    err = math.sqrt(stat**2 + syst**2)
    return err                           

def getFactor(x,  xstat, xsyst, y, ystat, ysyst):
    z = x*y
    err = z* math.sqrt((xstat/x)**2 + (xsyst/x)**2 + (ystat/y)**2 + (ysyst/y)**2)
    return z, err


def getRMueError(norm, up, dn):
    final = copy.deepcopy(norm)
    for i in range(1,norm.GetNbinsX()+1):
        stat = norm.GetBinError(i)
        syst = max( abs(norm.GetBinContent(i) - up.GetBinContent(i)), abs(norm.GetBinContent(i) - dn.GetBinContent(i)))
        final.SetBinError(i, math.sqrt(stat**2 + syst**2))
    return final                                                                                                                 





def weightedAverage( factor, direct, unwght):
    
    rsfof_factor = copy.deepcopy(factor)
    rsfof_direct = copy.deepcopy(factor)
    rsfof_final  = copy.deepcopy(factor)
    nof_final    = copy.deepcopy(unwght)
    
    for i in range(1, unwght.GetNbinsX()+1):
        if not unwght.GetBinContent(i): 
            rsfof_factor.SetBinContent( i, 0.)             
            rsfof_direct.SetBinContent( i, 0.) 
            rsfof_final .SetBinContent( i, 0.) 
            nof_final   .SetBinContent( i, 0.)    
            continue 
        rsfof_factor.SetBinContent(i,factor.GetBinContent(i) / unwght.GetBinContent(i))
        rsfof_factor.SetBinError  (i, math.sqrt(abs(factor.GetBinError(i)**2 - (rsfof_factor.GetBinContent(i)*unwght.GetBinError(i))**2)) / unwght.GetBinContent(i) )
        rsfof_direct.SetBinContent(i, direct.GetBinContent(i) / unwght.GetBinContent(i))
        rsfof_direct.SetBinError  (i, math.sqrt(abs(direct.GetBinError(i)**2 - (rsfof_direct.GetBinContent(i)*unwght.GetBinError(i))**2)) / unwght.GetBinContent(i) )
        
        rsfof_final .SetBinContent(i, (rsfof_factor.GetBinContent(i) / rsfof_factor.GetBinError(i)**2 + rsfof_direct.GetBinContent(i) / rsfof_factor.GetBinError(i)**2) / ( 1. / rsfof_factor.GetBinError(i)**2 + 1. / rsfof_factor.GetBinError(i)**2 ))
        rsfof_final .SetBinError  (i, 1 / math.sqrt(1/rsfof_factor.GetBinError(i)**2 + 1/rsfof_direct.GetBinError(i)**2))
        
        nof_final.SetBinContent(i, unwght.GetBinContent(i)*rsfof_final.GetBinContent(i))
        nof_final.SetBinError  (i, math.sqrt( (unwght.GetBinError(i)*rsfof_final.GetBinContent(i))**2 + (unwght.GetBinContent(i)*rsfof_final.GetBinError(i))**2))

    return [rsfof_factor, rsfof_direct, rsfof_final, nof_final]                                                                              


def makeResultsTable(da, fs, dy, zz, wz, ttz, others, vvv, region):
    line00 = '\\begin{tabular}{r l c c c c}  '
    if region == "TChiWZ":
        line0 = 'TChiWZ  & 50-100 & 100-150 & 150-250 & 250-350 & 350+\\\\ \\hline  '
    if region == "TChiZH":
        line0 = 'TChiZH  & 50-100 & 100-150 & 150-250 & 250+\\\\ \\hline  '
    
    line1 = '  Templates  & '
    line2 = '  FS  & ' 
    line3 = '  ZZ  & ' 
    line4 = '  WZ  & ' 
    line5 = '  ttZ  & ' 
    line6 = '  VVV  & ' 
    line7 = '  others  & ' 
    line8 = '  total     & ' 
    line9 = '  obs.      & ' 
    line10 = '  \end{tabular}  ' 

    my_range = range(1,da.GetNbinsX()+1)
    for i in my_range:
        tmp_dy   = dy.GetBinContent(i)  ; tmp_dy_e   = dy.GetBinError(i)
        tmp_fs   = fs.GetBinContent(i)  ; tmp_fs_e   = fs.GetBinError(i)
        tmp_zz   = zz.GetBinContent(i)  ; tmp_zz_e   = zz.GetBinError(i)
        tmp_wz   = wz.GetBinContent(i)  ; tmp_wz_e   = wz.GetBinError(i)
        tmp_ttz  = ttz.GetBinContent(i) ; tmp_ttz_e = ttz.GetBinError(i)
        tmp_vvv  = vvv.GetBinContent(i) ; tmp_vvv_e = vvv.GetBinError(i)
        tmp_others  = others.GetBinContent(i) ; tmp_others_e = others.GetBinError(i)
        tmp_da   = da.GetBinContent(i)  ; tmp_da_e   = da.GetBinError(i)
        tmp_full = tmp_dy + tmp_fs + tmp_zz + tmp_wz + tmp_ttz + tmp_others  ; tmp_full_e = math.sqrt(tmp_dy_e**2 + tmp_fs_e**2+ tmp_zz_e**2+ tmp_wz_e**2+ tmp_ttz_e**2+ tmp_others_e**2)
        #tmp_full = tmp_dy + tmp_fs ; tmp_full_e = math.sqrt(tmp_dy_e**2 + tmp_fs_e**2)
        line1 += '  %.2f $\\pm$ %.2f      %s' %(tmp_dy  , tmp_dy_e  , ' & ' if i != max(my_range) else '\\\\')
        line2 += '  %.2f $\\pm$ %.2f      %s' %(tmp_fs  , tmp_fs_e  , ' & ' if i != max(my_range) else '\\\\')
        line3 += '  %.2f $\\pm$ %.2f      %s' %(tmp_zz  , tmp_zz_e  , ' & ' if i != max(my_range) else '\\\\')
        line4 += '  %.2f $\\pm$ %.2f      %s' %(tmp_wz  , tmp_wz_e  , ' & ' if i != max(my_range) else '\\\\')
        line5 += '  %.2f $\\pm$ %.2f      %s' %(tmp_ttz , tmp_ttz_e , ' & ' if i != max(my_range) else '\\\\')
        line6 += '  %.2f $\\pm$ %.2f      %s' %(tmp_vvv , tmp_vvv_e , ' & ' if i != max(my_range) else '\\\\')
        line7 += '  %.2f $\\pm$ %.2f      %s' %(tmp_others , tmp_others_e , ' & ' if i != max(my_range) else '\\\\')
        line8 += '  %.2f $\\pm$ %.2f      %s' %(tmp_full, tmp_full_e, ' & ' if i != max(my_range) else '\\\\')
        line9 += '  %.2f       %s' %(tmp_da , ' & ' if i != max(my_range) else '\\\\')
    line7 += '\\hline'; line9 += '\\hline'
    #line8 += '\\hline'

    print line00
    print line0
    print line1
    print line2
    print line3
    print line4
    print line5
    print line6
    print line7
    print line8
    print line9
    print line10
    

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


def scaleByRSFOF(histo, rsfof, rsfof_err):
    h_rsfof = copy.deepcopy(histo)
    h_rsfof.SetName('h_rsfof')
    for i in range(1, h_rsfof.GetNbinsX()+1):
        h_rsfof.SetBinContent(i,rsfof)
        h_rsfof.SetBinError  (i,rsfof_e)
    histo.Multiply(h_rsfof)
    return histo                                   
def scaleByEWKFactors(histo, kappa):
    h_scaled = copy.deepcopy(histo)
    h_scaled.SetName('h_scaled')
    for i in range(1, h_scaled.GetNbinsX()+1):
        h_scaled.SetBinContent(i,kappa)
        h_scaled.SetBinError  (i,kappa*0.3)
    histo.Multiply(h_scaled)
    return histo                                  

def makeTheFactors():
    lint = 35.9  ; maxrun = 999999; lint_str = '35.9invfb'
    
    mll = 'm_{ll} [GeV]' 
    bins = [100, 150, 250]

    mll_da=treeDA.getTH1F(lint,"mll_da",'met_Edge', bins, 1, 1,  cuts.AddList([cuts.goodLepton,cuts.baseline,cuts.Zmass, cuts.ThirdLeptonVeto, cuts.OF]), '',mll)
    mll_mc=treeMC.getTH1F(lint,"mll_mc",'met_Edge', bins, 1, 1,  cuts.AddList([cuts.goodLepton,cuts.baseline,cuts.Zmass, cuts.ThirdLeptonVeto,cuts.OF]), '',mll)
    mll_ext_da=treeDA.getTH1F(lint,"mll_ext_da",'met_Edge', bins, 1, 1,  cuts.AddList([cuts.goodLepton,cuts.baseline, cuts.ThirdLeptonVeto,cuts.OF]), '',mll)
    mll_ext_mc=treeMC.getTH1F(lint,"mll_ext_mc",'met_Edge', bins, 1, 1,  cuts.AddList([cuts.goodLepton,cuts.baseline, cuts.ThirdLeptonVeto,cuts.OF]), '',mll)     
    mll_ZH_da=treeDA.getTH1F(lint,"mll_ZH_da",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.ewinoZH, cuts.OF,cuts.METg100]), '',mll)
    mll_ZH_mc=treeMC.getTH1F(lint,"mll_ZH_mc",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.ewinoZH, cuts.OF,cuts.METg100]), '',mll)     
    mll_ext_ZH_da=treeDA.getTH1F(lint,"mll_ZH_ext_da",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.ewinoZHExtMll, cuts.OF, cuts.METg100]), '',mll)
    mll_ext_ZH_mc=treeMC.getTH1F(lint,"mll_ZH_ext_mc",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.ewinoZHExtMll, cuts.OF, cuts.METg100]), '',mll)                    
    mll_WZ_da=treeDA.getTH1F(lint,"mll_WZ_da",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.ewinoWZ, cuts.OF,cuts.METg100]), '',mll)
    mll_WZ_mc=treeMC.getTH1F(lint,"mll_WZ_mc",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.ewinoWZ, cuts.OF,cuts.METg100]), '',mll)     
    mll_ext_WZ_da=treeDA.getTH1F(lint,"mll_WZ_ext_da",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.ewinoWZExtMll, cuts.OF, cuts.METg100]), '',mll)
    mll_ext_WZ_mc=treeMC.getTH1F(lint,"mll_WZ_ext_mc",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.ewinoWZExtMll, cuts.OF, cuts.METg100]), '',mll)                   
    mll_onZ_b_da=treeDA.getTH1F(lint,"mll_onZ_b_da",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.strongOnZWithB, cuts.METg100, cuts.Zmass,  cuts.OF]), '',mll)
    mll_onZ_b_mc=treeMC.getTH1F(lint,"mll_onZ_b_mc",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.strongOnZWithB, cuts.METg100, cuts.Zmass,  cuts.OF]), '',mll)     
    mll_ext_onZ_b_da=treeDA.getTH1F(lint,"mll_onZ_b_ext_da",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.strongOnZWithB, cuts.METg100 , cuts.OF]), '',mll)
    mll_ext_onZ_b_mc=treeMC.getTH1F(lint,"mll_onZ_b_ext_mc",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.strongOnZWithB, cuts.METg100 , cuts.OF]), '',mll)     
    mll_onZ_bveto_da=treeDA.getTH1F(lint,"mll_onZ_bveto_da",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.strongOnZBVeto, cuts.Zmass,cuts.OF]), '',mll)
    mll_onZ_bveto_mc=treeMC.getTH1F(lint,"mll_onZ_bveto_mc",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.strongOnZBVeto, cuts.Zmass,cuts.OF]), '',mll)     
    mll_ext_onZ_bveto_da=treeDA.getTH1F(lint,"mll_onZ_bveto_ext_da",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.strongOnZBVeto, cuts.METg100 , cuts.OF]), '',mll)
    mll_ext_onZ_bveto_mc=treeMC.getTH1F(lint,"mll_onZ_bveto_ext_mc",'lepsMll_Edge', 1, 0, 1000,  cuts.AddList([cuts.goodLepton,cuts.strongOnZBVeto, cuts.METg100 , cuts.OF]), '',mll)   

    k_base_mc = mll_mc.Clone("mll_mc");k_base_mc.Divide(mll_mc, mll_ext_mc, 1, 1, "B") ; 
    k_base_da = mll_da.Clone("mll_da");k_base_da.Divide(mll_da, mll_ext_da, 1, 1, "B") ; 
    k_met100_150_mc = k_base_mc.GetBinContent(1); k_met100_150_mc_e = k_base_mc.GetBinError(1);
    k_met100_150_da = k_base_da.GetBinContent(1); k_met100_150_da_e = k_base_da.GetBinError(1);
    k_met150_250_mc = k_base_mc.GetBinContent(2); k_met150_250_mc_e = k_base_mc.GetBinError(2);
    k_met150_250_da = k_base_da.GetBinContent(2); k_met150_250_da_e = k_base_da.GetBinError(2);
    k_met250_mc = k_base_mc.GetBinContent(3); k_met250_mc_e = k_base_mc.GetBinError(3);             
    k_met250_da = k_base_da.GetBinContent(3); k_met250_da_e = k_base_da.GetBinError(3);             
    ZH_mc = mll_ZH_mc.Clone("mll_ZH_mc");ZH_mc.Divide(mll_ZH_mc, mll_ext_ZH_mc, 1, 1, "B");k_ZH_mc = ZH_mc.GetBinContent(1); k_ZH_mc_e = ZH_mc.GetBinError(1); 
    ZH_da = mll_ZH_da.Clone("mll_ZH_da");ZH_da.Divide(mll_ZH_da, mll_ext_ZH_da, 1, 1, "B");k_ZH_da = ZH_da.GetBinContent(1); k_ZH_da_e = ZH_da.GetBinError(1);   
    WZ_mc = mll_WZ_mc.Clone("mll_WZ_mc");WZ_mc.Divide(mll_WZ_mc, mll_ext_WZ_mc, 1, 1, "B");k_WZ_mc = WZ_mc.GetBinContent(1); k_WZ_mc_e = WZ_mc.GetBinError(1); 
    WZ_da = mll_WZ_da.Clone("mll_WZ_da");WZ_da.Divide(mll_WZ_da, mll_ext_WZ_da, 1, 1, "B");k_WZ_da = WZ_da.GetBinContent(1); k_WZ_da_e = WZ_da.GetBinError(1);   
    onZ_b_mc = mll_onZ_b_mc.Clone("mll_onZ_b_mc");onZ_b_mc.Divide(mll_onZ_b_mc, mll_ext_onZ_b_mc, 1, 1, "B");k_onZ_b_mc = onZ_b_mc.GetBinContent(1); k_onZ_b_mc_e = onZ_b_mc.GetBinError(1); 
    onZ_b_da = mll_onZ_b_da.Clone("mll_onZ_b_da");onZ_b_da.Divide(mll_onZ_b_da, mll_ext_onZ_b_da, 1, 1, "B");k_onZ_b_da = onZ_b_da.GetBinContent(1); k_onZ_b_da_e = onZ_b_da.GetBinError(1); 
    onZ_bveto_mc = mll_onZ_bveto_mc.Clone("mll_onZ_bveto_mc");onZ_bveto_mc.Divide(mll_onZ_bveto_mc, mll_ext_onZ_bveto_mc, 1, 1, "B");k_onZ_bveto_mc = onZ_bveto_mc.GetBinContent(1); k_onZ_bveto_mc_e = onZ_bveto_mc.GetBinError(1); 
    onZ_bveto_da = mll_onZ_bveto_da.Clone("mll_onZ_bveto_da");onZ_bveto_da.Divide(mll_onZ_bveto_da, mll_ext_onZ_bveto_da, 1, 1, "B");k_onZ_bveto_da = onZ_bveto_da.GetBinContent(1); k_onZ_bveto_da_e = onZ_bveto_da.GetBinError(1); 

    c1 = r.TCanvas()
    l1 = r.TLegend(0.75,0.8,0.9,0.9)
    fmll = r.TH1F('k_mc','f_{mll}',7,0,7)
    fmll.SetBinContent(1, k_met100_150_mc);fmll.SetBinError(1,k_met100_150_mc_e);fmll.GetXaxis().SetBinLabel(1, 'Baseline E_{T}^{miss} 100-150 GeV');
    fmll.SetBinContent(2, k_met150_250_mc);fmll.SetBinError(2,k_met150_250_mc_e);fmll.GetXaxis().SetBinLabel(2, 'Baseline E_{T}^{miss} 150-250 GeV');
    fmll.SetBinContent(3, k_met250_mc);fmll.SetBinError(3,k_met250_mc_e);fmll.GetXaxis().SetBinLabel(3, 'Baseline E_{T}^{miss} > 250 GeV');
    fmll.SetBinContent(4, k_onZ_bveto_mc);fmll.SetBinError(4,k_onZ_bveto_mc_e);fmll.GetXaxis().SetBinLabel(4, 'OnZ bveto ');
    fmll.SetBinContent(5, k_onZ_b_mc);fmll.SetBinError(5,k_onZ_b_mc_e);fmll.GetXaxis().SetBinLabel(5, 'OnZ with b ');
    fmll.SetBinContent(6, k_ZH_mc);fmll.SetBinError(6,k_ZH_mc_e);fmll.GetXaxis().SetBinLabel(6, 'TChiZH ');
    fmll.SetBinContent(7, k_WZ_mc);fmll.SetBinError(7,k_WZ_mc_e);fmll.GetXaxis().SetBinLabel(7, 'TChiWZ ');
    fmll.SetMarkerColor(r.kTeal-5)
    fmll.SetLineColor  (r.kTeal-5)
    fmll.SetLineWidth  (2)
    r.gStyle.SetPaintTextFormat('.3f')
    fmll.Draw('E1, X0')
    fmll.Draw('text00 same')
    fmll_da = fmll.Clone('data fracs')
    l1.AddEntry(fmll, 'MC', 'pl')
    fmll_da.SetBinContent(1, k_met100_150_da);fmll_da.SetBinError(1,k_met100_150_da_e);fmll_da.GetXaxis().SetBinLabel(1, 'Baseline E_{T}^{miss} 100-150 GeV');
    fmll_da.SetBinContent(2, k_met150_250_da);fmll_da.SetBinError(2,k_met150_250_da_e);fmll_da.GetXaxis().SetBinLabel(2, 'Baseline E_{T}^{miss} 150-250 GeV');
    fmll_da.SetBinContent(3, k_met250_da);fmll_da.SetBinError(3,k_met250_da_e);fmll_da.GetXaxis().SetBinLabel(3, 'Baseline E_{T}^{miss} >250 GeV');
    fmll_da.SetBinContent(4, k_onZ_bveto_da);fmll_da.SetBinError(4,k_onZ_bveto_da_e);fmll_da.GetXaxis().SetBinLabel(4, 'OnZ bveto');
    fmll_da.SetBinContent(5, k_onZ_b_da);fmll_da.SetBinError(5,k_onZ_b_da_e);fmll_da.GetXaxis().SetBinLabel(5, 'OnZ with b');
    fmll_da.SetBinContent(6, k_ZH_da);fmll_da.SetBinError(6,k_ZH_da_e);fmll_da.GetXaxis().SetBinLabel(6, 'TChiZH');
    fmll_da.SetBinContent(7, k_WZ_da);fmll_da.SetBinError(7,k_WZ_da_e);fmll_da.GetXaxis().SetBinLabel(7, 'TChiWZ');
    fmll_da.GetYaxis().SetRangeUser(0.0,0.2);
    fmll_da.GetYaxis().SetTitle('f_{mll}')
    fmll_da.SetMarkerStyle(22)
    fmll_da.SetMarkerColor(r.kBlack)
    fmll_da.SetLineColor  (r.kBlack)
    fmll_da.SetLineWidth  (2)
    fmll_da.Draw('E1, X0')
    fmll_da.Draw('text00 same ')
    fmll.Draw('E1, X0 same')
    fmll.Draw('text00 same')
    fmll_da.Draw('text00 same ')
    fmll_da.Draw('E1, X0 same ')
    line1 =  r.TLine(fmll_da.GetXaxis().GetXmin(), 0.065, fmll_da.GetXaxis().GetXmax (), 0.065);line1.SetLineColor (r.kBlack);  line1.Draw();
    line2 =  r.TLine(fmll_da.GetXaxis().GetXmin(), 0.1, fmll_da.GetXaxis().GetXmax (), 0.1);line2.SetLineColor (r.kBlack); line2.SetLineStyle(2); line2.Draw();
    line3 =  r.TLine(fmll_da.GetXaxis().GetXmin(), 0.03, fmll_da.GetXaxis().GetXmax (), 0.03);line3.SetLineColor (r.kBlack); line3.SetLineStyle(2); line3.Draw();
    l1.AddEntry(fmll_da, 'data', 'pl')
    l1.Draw('same')
    lat1 = r.TLatex(); lat1.SetNDC(); lat1.SetTextFont(fmll.GetYaxis().GetTitleFont()); lat1.SetTextSize(fmll.GetYaxis().GetTitleSize())
    lat1.DrawLatex(0.76, 0.93, '%.2f fb^{-1}'%(lint))
    path = 'plots/ewino/fmll.png'
    helper.ensurePath(path)
    c1.SaveAs(path)                                                                                                                                  
    
    
    return k_met150_250_da, k_met150_250_da_e, k_met150_250_mc, k_met150_250_mc_e

def makeDYMETShape(var, specialcut = '', scutstring = '', doCumulative = False, region = ''):
    if var == 'met':
        treevar = 'met_Edge'
        xlabel = 'E_{T}^{miss} [GeV]'                    
        nbins, xmin, xmax = 15, 0, 300
    if region == "TChiWZ":
        bin1  =795.0 ; bin1_e = 32.6;bin2  = 28.1; bin2_e = 5.5;bin3  = 4.1; bin3_e = 2.8;bin4  = 0.8; bin4_e = 1.0;bin5  = 0.00; bin5_e = 0.5;
        metbins = [50.0, 100.0, 150.0, 250.0, 350.0, 450.0]
        metbins_ = array('d', [50.0, 100.0, 150.0, 250.0, 350.0, 450.0])
        nbin = 5
    if region == "TChiZH":
        bin1  = 77.8; bin1_e = 9.5;bin2  = 2.9; bin2_e =2.1;bin3  = 0.2; bin3_e = 0.3;bin4  = 0.1; bin4_e = 0.1
        metbins = [50.0, 100.0, 150.0, 250.0, 350.0]
        metbins_ = array('d', [50.0, 100.0, 150.0, 250.0, 350.0])
        nbin = 4
    lint = 35.9  ; maxrun = 999999; lint_str = '35.9invfb'
    dy = r.TH1F('dy'+region,'dy'+region, nbin, metbins_)
    dy.SetBinContent(1, bin1);dy.SetBinError(1,bin1_e);
    dy.SetBinContent(2, bin2);dy.SetBinError(2,bin2_e);
    dy.SetBinContent(3, bin3);dy.SetBinError(3,bin3_e);
    dy.SetBinContent(4, bin4);dy.SetBinError(4,bin4_e);
    if region == "TChiWZ":
        dy.SetBinContent(5, bin5);dy.SetBinError(5,bin5_e);

    plot = Canvas.Canvas('ewino/%s/templates_%s'%(lint_str, region), 'png,pdf,root', 0.65, 0.6, 0.85, 0.9)
    plot.addHisto(dy, "HIST, SAME", "templates", "L", r.kBlack,  1, 0)
    plot.save(1, 1, 1, lint)                                                                              
    return dy                             
    del dy, plot
    
    
def makeClosureTests(var, specialcut = '', scutstring = '', doCumulative = False, region = ''):

    if var == 'mll':
        treevar = 'lepsMll_Edge'
        xlabel = 'm_{ll} (GeV)'              
    elif var == 'met':
        treevar = 'met_Edge'
        xlabel = 'E_{T}^{miss} [GeV]'                    
    lint = 35.9  ; maxrun = 999999 ; lint_str = '35.9invfb'
    if region == "TChiWZ":
        regioncut =  cuts.ewinoWZExtMll
        bins = [50.0, 100.0, 150.0, 250.0, 350.0]
    if region == "TChiZH":
        regioncut =  cuts.ewinoZHExtMll
        bins = [50.0, 100.0, 150.0, 250.0]
    kappa_da, kappa_da_e, kappa_mc, kappa_mc_e = makeTheFactors()
    kappa_da, kappa_da_e, kappa_mc, kappa_mc_e = 0.065, 0.01, 0.065, 0.01
    ## mll distributions
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","A","DATA")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","B","DATA")
    rsfof_mc = helper.readFromFileRsfofD("ingredients.dat", "MC") 
    rt_mc = helper.readFromFileRT("ingredients.dat", "MC")
    rmue_a_mc = helper.readFromFileRmueCoeff("ingredients.dat","A","MC")
    rmue_b_mc = helper.readFromFileRmueCoeff("ingredients.dat","B","MC")

    da_OF =        treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF]), '', xlabel)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    print 'trivial check'
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue
                                                                                                 
    da_OF_direct = copy.deepcopy(da_OF)
    da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], rsfof_da[1])
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    result = weightedAverage( da_OF_factor, da_OF_direct, da_OF)                                 
    prediction = copy.deepcopy( da_OF )
    da_OF.SetBinErrorOption( TH1.kPoisson)                                                                                                                                  
    for i in range(1, prediction.GetNbinsX()+1):
        prediction.SetBinError(i,0.)
    prediction.Multiply(result[2])                                                                                                                                                                                                                                                                                                      
   
    mc_OF =        treeFS.getTH1F(lint, var+"mc_OF"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF]), '', xlabel)
    mc_OF_factor = treeFS.getTH1F(lint, var+"mc_OF_factor"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))
    mc_OF_factorUp = treeFS.getTH1F(lint, var+"mc_OF_factorUp"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))
    mc_OF_factorDn = treeFS.getTH1F(lint, var+"mc_OF_factorDn"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))
    print 'trivial check'
    for i in range(1, mc_OF.GetNbinsX()+1):
        if not  mc_OF.GetBinContent(i): continue
                                                                                                 
    mc_OF_direct = copy.deepcopy(mc_OF)
    mc_OF_direct = scaleByRSFOF(mc_OF_direct, rsfof_mc[0], rsfof_mc[1])
    mc_OF_factor = getRMueError(mc_OF_factor, mc_OF_factorUp, mc_OF_factorDn)
    mc_OF_factor = scaleByRSFOF(mc_OF_factor, rt_mc[0], getFinalError(rt_mc[1],rt_mc[2]))
    mc_result = weightedAverage( mc_OF_factor, mc_OF_direct, mc_OF)                                 
    mc_prediction = copy.deepcopy( mc_OF )
    mc_OF.SetBinErrorOption( TH1.kPoisson)                                                                                                                                  
    for i in range(1, mc_prediction.GetNbinsX()+1):
        mc_prediction.SetBinError(i,0.)
    mc_prediction.Multiply(mc_result[2])                                                                                                                                                          

    mc_SF = treeFS.getTH1F(lint, var+"mc_SF"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut ,cuts.Zmass, cuts.SF]), '', xlabel)
    mc_SF.GetYaxis().SetRangeUser(0., 1.5*mc_SF.GetMaximum())
    print helper.bcolors.HEADER + '[MC only closure test not scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC

    mc_OF_fmllScaled = copy.deepcopy(mc_prediction)
    mc_OF_fmllScaled = scaleByEWKFactors(mc_OF_fmllScaled, kappa_mc)
    mc_OF_fmllScaled_err = copy.deepcopy(mc_OF_fmllScaled)
    mc_OF_fmllScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_fmllScaled_err.SetFillStyle(3004); mc_OF_fmllScaled_err.SetMarkerSize(0.)                                                                                    
    da_OF_fmllScaled = copy.deepcopy(prediction)
    da_OF_fmllScaled = scaleByEWKFactors(da_OF_fmllScaled, kappa_da)

    for i in range(1, mc_OF.GetNbinsX()+1):
        if not  mc_OF.GetBinContent(i): continue


    print helper.bcolors.HEADER + '[MC only closure test scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure = Canvas.Canvas('ewino/%s/plot_closure_%s_%s_%s'%(lint_str, var, '' if not scutstring else '_'+scutstring, region), 'png,pdf', 0.6, 0.55, 0.75, 0.8)
    plot_closure.addHisto(mc_SF                , 'PE'       , 'MC SF', 'PL', r.kRed+1  , 1,  0)
    plot_closure.addHisto(mc_OF_fmllScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure.addHisto(mc_OF_fmllScaled    , 'hist,SAME', 'MC (OF*R_{SFOF}*#kappa)', 'L' , r.kBlue+1 , 1,  1)
    plot_closure.saveRatio(1, 0, 0, lint, mc_SF, mc_OF_fmllScaled, 0.2, 1.8, "MC SF/MC OF")
    return da_OF_fmllScaled                                                                                                                              


def makeResultData(analysis, var, maxrun = 999999, lint = 35.9, specialcut = '', scutstring = '', region = '', _options = ''):
    scan = Scans.Scan(analysis)
    lint = 35.9
    print "Doing region: ", region
    dy_shape =  makeDYMETShape('met','', '', True, region)
    fs_shape  = makeClosureTests('met','','', True, region)
    kappa_da, kappa_da_e, kappa_mc, kappa_mc_e = 0.065, 0.01, 0.065, 0.01
    #kappa_da, kappa_da_e, kappa_mc, kappa_mc_e = makeTheFactors()
    returnplot, addRares, splitFlavor, makeTable, printIntegral = False, True, False, False, False
    if 'returnplot'    in _options: print 'found option %s'%'returnplot'    ;returnplot    = True
    if 'splitFlavor'   in _options: print 'found option %s'%'splitFlavor'   ;splitFlavor   = True
    if 'makeTable'     in _options: print 'found option %s'%'makeTable'     ;makeTable     = True
    if 'printIntegral' in _options: print 'found option %s'%'printIntegral' ;printIntegral = True

    if var == 'met'      : treevar = 'met_Edge'     ; xlabel = 'E_{T}^{miss.} [GeV]'
    if region == "TChiWZ":
        bins = [50.0, 100.0, 150.0, 250.0, 350.0]
        regioncut = cuts.ewinoWZExtMll
    if region == "TChiZH":
        bins = [50.0, 100.0, 150.0, 250.0]
        regioncut = cuts.ewinoZHExtMll
    mc_stack = r.THStack() 
    newLumiString = '35.9invfb'

    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","A","DATA")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","B","DATA")
                                                                                                      
    da_OF =        treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF, cuts.trigger]), '', xlabel)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF, cuts.trigger]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF, cuts.trigger]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, regioncut, cuts.OF, cuts.trigger]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    print 'trivial check'
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue
                                                                                                 
    da_OF_direct = copy.deepcopy(da_OF)
    da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], rsfof_da[1])
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    result = weightedAverage( da_OF_factor, da_OF_direct, da_OF)                                 
    prediction = copy.deepcopy( da_OF )
    da_OF.SetBinErrorOption( TH1.kPoisson)                                                                                                                                  
    for i in range(1, prediction.GetNbinsX()+1):
        prediction.SetBinError(i,0.)
    prediction.Multiply(result[2])                                                                                                                                                                                                                                                                                                       
    for bin, label in scan.SRLabels.items():
        bin = prediction.FindBin(bin)
        # to get asymmetric errors (i dont know why it doesnt work another way)
        dummyHisto = r.TH1F()
        dummyHisto.SetBinErrorOption(TH1.kPoisson)
        for i in range(1, int(da_OF.GetBinContent(bin))+1):
            dummyHisto.Fill(0.5)
        dummyHisto.GetBinContent(1), '+/-', dummyHisto.GetBinErrorUp(1), dummyHisto.GetBinErrorLow(1)
        syst = '  {value:4.1f}^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(value = prediction.GetBinContent(bin),
                                                                             errUp = prediction.GetBinErrorUp(bin),
                                                                             errDn = prediction.GetBinErrorLow(bin))
        stat = '^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(errUp = dummyHisto.GetBinErrorUp(1) * result[2].GetBinContent(bin),
                                                               errDn = dummyHisto.GetBinErrorLow(1) * result[2].GetBinContent(bin))
        del dummyHisto                                                                                                                           

    da_OF_fmllScaled = copy.deepcopy(prediction)
    da_OF_fmllScaled = scaleByEWKFactors(da_OF_fmllScaled, kappa_da)
    da_OF_fmllScaled_err = copy.deepcopy(da_OF_fmllScaled)
    da_OF_fmllScaled_err.SetFillColorAlpha(r.kBlack, 0.8)
    da_OF_fmllScaled_err.SetFillStyle(3004); da_OF_fmllScaled_err.SetMarkerSize(0.)       
    
    
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar,    bins, 1,1, cuts.AddList([specialcut, cuts.Zmass, cuts.goodLepton, regioncut, cuts.SF]), '', xlabel)
    ra_SF = treeRA.getTH1F(lint, var+"ra_SF"+scutstring, treevar,    bins, 1,1, cuts.AddList([specialcut, cuts.lepsFromZ, cuts.Zmass, cuts.goodLepton, regioncut, cuts.SF]), '', xlabel)
    ttz_SF = treeTTZ.getTH1F(lint, var+"ttz_SF"+scutstring, treevar, bins, 1,1, cuts.AddList([specialcut, cuts.lepsFromZ, cuts.Zmass, cuts.goodLepton, regioncut, cuts.SF]), '', xlabel)
    vvv_SF = treeVVV.getTH1F(lint, var+"vvv_SF"+scutstring, treevar, bins, 1,1, cuts.AddList([specialcut, cuts.lepsFromZ, cuts.Zmass, cuts.goodLepton, regioncut, cuts.SF]), '', xlabel)
    zz_SF = treeZZ.getTH1F(lint, var+"zz_SF"+scutstring, treevar,    bins, 1,1, cuts.AddList([specialcut, cuts.lepsFromZ, cuts.Zmass, cuts.goodLepton, regioncut, cuts.SF]), '', xlabel)
    wz_SF = treeWZ.getTH1F(lint, var+"wz_SF"+scutstring, treevar,    bins, 1,1, cuts.AddList([specialcut, cuts.lepsFromZ, cuts.Zmass, cuts.goodLepton, regioncut, cuts.SF]), '', xlabel)
    da_SF.SetTitle("data SF")
    vvv_SF.SetFillColorAlpha(r.kGreen-3, 0.8);vvv_SF.SetTitle("VVV");vvv_SF.SetLineColor(r.kBlack)
    ra_SF.SetFillColorAlpha(r.kGreen-5, 0.8);ra_SF.SetTitle("rares");ra_SF.SetLineColor(r.kBlack)
    ttz_SF.SetFillColorAlpha(r.kBlue-7 , 0.8);ttz_SF.SetTitle("ttZ");ttz_SF.SetLineColor(r.kBlack)                                                             
    zz_SF.SetFillColorAlpha(r.kCyan+2 , 0.8);zz_SF.SetTitle("ZZ");zz_SF.SetLineColor(r.kBlack)                                                             
    wz_SF.SetFillColorAlpha(r.kGreen-8 , 0.8);wz_SF.SetTitle("WZ");wz_SF.SetLineColor(r.kBlack)                                                             
    dy_shape.SetFillColorAlpha(r.kYellow-9, 0.8);dy_shape.SetTitle("E_{T}^{miss} templates")
    da_OF_fmllScaled.SetFillColorAlpha(r.kRed-9, 0.8);da_OF_fmllScaled.SetTitle("FS")                                  
    #sorry these are the multilepton control region scale factor and they are hard coded for now, it's uglyy
    ttz_SF.Scale(1.36);wz_SF.Scale(1.06);zz_SF.Scale(1.71) 
    for bin, label in scan.SRLabels.items():
        if ttz_SF.GetBinContent(bin) > 0:
            ttz_SF.SetBinError(bin, ttz_SF.GetBinContent(bin)*math.sqrt(((ttz_SF.GetBinContent(bin)*0.3)/ttz_SF.GetBinContent(bin))**2 + (ttz_SF.GetBinError(bin)/ttz_SF.GetBinContent(bin)**2)) )
        else: ttz_SF.SetBinError(bin , ttz_SF.GetBinError(bin))                                                                                                                                    
        if zz_SF.GetBinContent(bin) > 0:
            zz_SF.SetBinError(bin, zz_SF.GetBinContent(bin)*math.sqrt(((zz_SF.GetBinContent(bin)*0.5)/zz_SF.GetBinContent(bin))**2 + (zz_SF.GetBinError(bin)/zz_SF.GetBinContent(bin)**2)) )
        else: zz_SF.SetBinError(bin , zz_SF.GetBinError(bin))                                                                                                                                    
        if wz_SF.GetBinContent(bin) > 0:
            wz_SF.SetBinError(bin, wz_SF.GetBinContent(bin)*math.sqrt(((wz_SF.GetBinContent(bin)*0.3)/wz_SF.GetBinContent(bin))**2 + (wz_SF.GetBinError(bin)/wz_SF.GetBinContent(bin)**2)) )
        else: wz_SF.SetBinError(bin , wz_SF.GetBinError(bin)) 

        if ra_SR.GetBinContent(bin) > 0:
            ra_SR.SetBinError(bin, ra_SR.GetBinContent(bin)*math.sqrt(((ra_SR.GetBinContent(bin)*0.3)/ra_SR.GetBinContent(bin))**2 + (ra_SR.GetBinError(bin)/ra_SR.GetBinContent(bin)**2)) )
        else: ra_SR.SetBinError(bin , ra_SR.GetBinError(bin)) 
        
        if vvv_SF.GetBinContent(bin) > 0:
            vvv_SF.SetBinError(bin, vvv_SF.GetBinContent(bin)*math.sqrt(((vvv_SF.GetBinContent(bin)*0.3)/vvv_SF.GetBinContent(bin))**2 + (vvv_SF.GetBinError(bin)/vvv_SF.GetBinContent(bin)**2)) )
        else: vvv_SF.SetBinError(bin , vvv_SF.GetBinError(bin)) 


    mc_full = copy.deepcopy(zz_SF)
    other_datacard = copy.deepcopy(zz_SF)

    mc_stack.Add(vvv_SF); mc_full.Add(vvv_SF, 1.); other_datacard.Add(vvv_SF,1.)
    mc_stack.Add(ttz_SF); mc_full.Add(ttz_SF, 1.); other_datacard.Add(ttz_SF,1.)
    mc_stack.Add(ra_SF); mc_full.Add(ra_SF, 1.)  ; other_datacard.Add(ra_SF ,1.)
    mc_stack.Add(wz_SF); mc_full.Add(wz_SF, 1.)  ; other_datacard.Add(wz_SF ,1.)
    mc_stack.Add(da_OF_fmllScaled); mc_full.Add(da_OF_fmllScaled, 1.) 
    mc_stack.Add(zz_SF);                           
    mc_stack.Add(dy_shape); mc_full.Add(dy_shape); 
    mc_stack.Draw()
    mc_stack.GetXaxis().SetTitle(xlabel)
    mc_full_e = copy.deepcopy(mc_full)
    mc_full_e.SetFillColorAlpha(r.kBlue+1, 0.8);mc_full_e.SetFillStyle(3017); mc_full_e.SetMarkerSize(0.)
    
    maxCont = max(da_OF_fmllScaled.GetMaximum(), da_SF.GetMaximum())
    da_SF.GetYaxis().SetRangeUser(0.1, 1.30*maxCont)
    mc_stack.SetMaximum(1.3*maxCont)
    SetOwnership(mc_stack, 0 );SetOwnership(da_SF, 0 )                                                                                                                       

    print helper.bcolors.HEADER + '[result scaled by RSFOF for DATA] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_result = Canvas.Canvas('ewino/%s/plot_result_%s_%s_%s'%(newLumiString, var, '' if not scutstring else '_'+scutstring, region), 'png,pdf', 0.64, 0.59, 0.91, 0.85)
    plot_result.addStack(mc_stack, "HIST" , 1, 1)
    plot_result.addHisto(mc_full_e, 'e2,same'  , '', 'PL', r.kBlack , 1, -1)
    plot_result.addHisto(da_SF, 'E1, SAME', 'data SF' , 'P', r.kBlack , 1, 0)
    plot_result.saveRatio(1, 1, 0, lint, da_SF, mc_full) 
    makeResultsTable(da_SF, da_OF_fmllScaled, dy_shape, zz_SF, wz_SF, ttz_SF, ra_SF, vvv_SF, region )

    #### for the datacards
    sergiosLabel = 'CharNeu_Moriond2017' if region == 'TChiWZ' else 'NeuNeu_Moriond2017'
    forDataCards = TFile.Open('datacards/forDatacards_%s.root'%sergiosLabel,'recreate')
    tf_CR_SR = scaleByEWKFactors(result[2], kappa_da)
    da_SF   .SetName('da_SF'   );     da_SF    .Write();      
    da_OF   .SetName('da_OF'   );     da_OF    .Write();  
    tf_CR_SR.SetName('tf_CR_SR');     tf_CR_SR .Write(); 
    dy_shape.SetName('dy_shape');     dy_shape .Write(); 
    other_datacard.SetName('mc_full' );     other_datacard.Write(); 
    forDataCards.Close()


    if makeTable:
        makeSimpleTable(plot_result, addRares)
    if returnplot:
        return plot_result
    del da_SF, ra_SF, ttz_SF, zz_SF, wz_SF, dy_shape, vvv_SF

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
    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    fsDatasets = [ 'TTJets_DiLepton_ext1', 'WWTo2L2Nu', 'WWW', 'TTWToLNu_ext2',  'TTWToQQ', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT']
    #fsDatasets = ['TTTT', 'TTHnobb_pow', 'VHToNonbb',  'TTJets_DiLepton_ext', 'TBar_tch_powheg', 'T_tch_powheg', 'WWTo2L2Nu', 'WWW', 'TTZToQQ', 'TTWToLNu',  'TTWToQQ', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT',   'WJetsToLNu_LO']
    zzDatasets = ['ZZTo2L2Nu']
    wzDatasets = ['WZTo3LNu']
    ttzDatasets = ['TTZToLLNuNu_ext1']
    raDatasets = ['TWZ', 'tZq_ll']
    vvvDatasets = ['WWZ','WZZ', 'ZZZ']
    mcDatasets = fsDatasets+dyDatasets + raDatasets + zzDatasets + wzDatasets + ttzDatasets + vvvDatasets
    
                                                                                    
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

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0)
    treeRA = Sample.Tree(helper.selectSamples(opts.sampleFile, raDatasets, 'RA'), 'RA'  , 0)
    treeVVV = Sample.Tree(helper.selectSamples(opts.sampleFile, vvvDatasets, 'VVV'), 'VVV'  , 0)
    treeFS = Sample.Tree(helper.selectSamples(opts.sampleFile, fsDatasets, 'FS'), 'FS'  , 0)
    treeWZ = Sample.Tree(helper.selectSamples(opts.sampleFile, wzDatasets, 'WZ'), 'WZ'  , 0)
    treeZZ = Sample.Tree(helper.selectSamples(opts.sampleFile, zzDatasets, 'ZZ'), 'ZZ'  , 0)
    treeTTZ = Sample.Tree(helper.selectSamples(opts.sampleFile, ttzDatasets, 'TTZ'), 'TTZ'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'

    isBlinded = True

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()
    lint = 35.9  ; maxrun = 999999; lint_str = '35.9invfb'
    print 'Running with an integrated luminosity of %.2f fb-1' %(lint)

    ## ============================================================
    ## ========== set RSFOF globally ==============================
    ## ============================================================
    global rsfof, rsfof_e, rsfof_mc, rsfof_mc_e
    rsfof, rsfof_e, rsfof_mc, rsfof_mc_e =  1.0449 ,     0.0641 , 1.0444 ,     0.0640
    ## result plots in different variables:
 
    #makeClosureTests('met','', 'inclusive', True)
    #makeClosureTests('met','', 'inclusive', True)
    #makeDYMETShape('met','', 'inclusive', True)
    
    makeResultData('Edge_Moriond2017', 'met', maxrun, lint, specialcut = "" , scutstring = '', region = 'TChiWZ')
    makeResultData('Edge_Moriond2017', 'met', maxrun, lint, specialcut = "" , scutstring = '', region = 'TChiZH')

