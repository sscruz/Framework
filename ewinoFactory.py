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

def getFraction(x,  xstat, y, ystat):
    if y > 0:
        z = x/y
        err = z* math.sqrt((xstat/x)**2 + (ystat/y)**2)
        return z, err  
    else: 
        print "denominator < 0!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        return x, xstat

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


def scaleByRSFOF(histo, rsfof, rsfof_err):
    h_rsfof = copy.deepcopy(histo)
    h_rsfof.SetName('h_rsfof')
    for i in range(1, h_rsfof.GetNbinsX()+1):
        h_rsfof.SetBinContent(i,rsfof)
        h_rsfof.SetBinError  (i,rsfof_e)
    histo.Multiply(h_rsfof)
    return histo                                   
def scaleByEWKFactors(histo, rsfof, rsfof_err, fmll, fmll_e, r0b1b, r0b1b_e):
    h_rsfof = copy.deepcopy(histo)
    h_rsfof.SetName('h_rsfof')
    factor = rsfof*fmll*r0b1b
    factor_e = rsfof*fmll*r0b1b*math.sqrt((rsfof_e/rsfof)**2 + (fmll_e/fmll)**2+ (r0b1b_e/r0b1b)**2 )
    for i in range(1, h_rsfof.GetNbinsX()+1):
        h_rsfof.SetBinContent(i,factor)
        h_rsfof.SetBinError  (i,factor_e)
    histo.Multiply(h_rsfof)
    return histo                                  



def makeTheFactors():
    lint = 12.9  ; maxrun = 999999; lint_str = '12.9invfb'
    
    mll = 'm_{ll} [GeV]' 
    mll_bveto_OF_mc=treeTT.getTH1F(lint,"mll_mc_bveto_OF",'lepsMll_Edge', 1, 0, 300, cuts.AddList([cuts.goodLepton, cuts.ewinoCR,cuts.Zmass,cuts.bVeto, cuts.OF]),'',mll)
    mll_bveto_SF_mc=treeTT.getTH1F(lint,"mll_mc_bveto_SF",'lepsMll_Edge',1, 0, 300, cuts.AddList([cuts.goodLepton, cuts.ewinoCR,cuts.Zmass,cuts.bVeto, cuts.SF]),'',mll)
    mll_binv_OF_mc =treeTT.getTH1F(lint,"mll_mc_binv_OF", 'lepsMll_Edge', 1, 0, 300, cuts.AddList([cuts.goodLepton, cuts.ewinoCR,cuts.Zmass,cuts.binv, cuts.OF]),'', mll)
    mll_ext_bveto_OF_mc=treeTT.getTH1F(lint,"mll_x_mc_bveto_OF",'lepsMll_Edge',1,0,300,cuts.AddList([cuts.goodLepton,cuts.ewinoCR,cuts.ZmassExtended,cuts.bVeto,cuts.OF]), '',mll)
    mll_ext_bveto_SF_mc=treeTT.getTH1F(lint,"mll_x_mc_bveto_SF",'lepsMll_Edge',1,0,300,cuts.AddList([cuts.goodLepton,cuts.ewinoCR,cuts.ZmassExtended,cuts.bVeto,cuts.SF]), '',mll)
    mll_ext_binv_OF_mc =treeTT.getTH1F(lint,"mll_x_mc_binv_OF", 'lepsMll_Edge',1,0,300,cuts.AddList([cuts.goodLepton,cuts.ewinoCR,cuts.ZmassExtended,cuts.binv,cuts.OF]), '', mll)
    mll_ext_binv_SF_mc =treeTT.getTH1F(lint,"mll_x_mc_binv_SF", 'lepsMll_Edge',1,0,300,cuts.AddList([cuts.goodLepton,cuts.ewinoCR,cuts.ZmassExtended,cuts.binv,cuts.SF]), '', mll)
    mll_bveto_OF_da=treeDA.getTH1F(lint,"mll_da_bveto_OF",'lepsMll_Edge', 1, 0, 300, cuts.AddList([cuts.goodLepton,cuts.ewinoCR,cuts.Zmass,cuts.bVeto, cuts.OF]), '',mll)
    mll_binv_OF_da =treeDA.getTH1F(lint,"mll_da_binv_OF", 'lepsMll_Edge', 1, 0, 300, cuts.AddList([cuts.goodLepton,cuts.ewinoCR,cuts.Zmass,cuts.binv, cuts.OF]), '', mll)
    mll_ext_bveto_OF_da=treeDA.getTH1F(lint,"mll_ext_da_bveto_OF",'lepsMll_Edge', 1,0, 300, cuts.AddList([cuts.goodLepton,cuts.ewinoCR,cuts.ZmassExtended,cuts.bVeto, cuts.OF]), '',mll)
    mll_ext_binv_OF_da =treeDA.getTH1F(lint,"mll_ext_da_binv_OF", 'lepsMll_Edge', 1,0, 300, cuts.AddList([cuts.goodLepton,cuts.ewinoCR,cuts.ZmassExtended,cuts.binv, cuts.OF]), '', mll)
    mll_ext_binv_SF_da =treeDA.getTH1F(lint,"mll_x_da_binv_SF", 'lepsMll_Edge',1,0,300,cuts.AddList([cuts.goodLepton,cuts.ewinoCR,cuts.ZmassExtended,cuts.binv,cuts.SF]), '', mll)

    mll_bveto_OF_mc_int = mll_bveto_OF_mc.GetBinContent(1); mll_bveto_OF_mc_e = mll_bveto_OF_mc.GetBinError(1);
    mll_bveto_SF_mc_int = mll_bveto_SF_mc.GetBinContent(1); mll_bveto_SF_mc_e = mll_bveto_SF_mc.GetBinError(1);
    mll_binv_OF_mc_int = mll_binv_OF_mc.GetBinContent(1); mll_binv_OF_mc_e = mll_binv_OF_mc.GetBinError(1);     
    mll_ext_bveto_OF_mc_int = mll_ext_bveto_OF_mc.GetBinContent(1); mll_ext_bveto_OF_mc_e = mll_ext_bveto_OF_mc.GetBinError(1);
    mll_ext_bveto_SF_mc_int = mll_ext_bveto_SF_mc.GetBinContent(1); mll_ext_bveto_SF_mc_e = mll_ext_bveto_SF_mc.GetBinError(1);
    mll_ext_binv_OF_mc_int =  mll_ext_binv_OF_mc.GetBinContent(1); mll_ext_binv_OF_mc_e = mll_ext_binv_OF_mc.GetBinError(1);     
    mll_ext_binv_SF_mc_int =  mll_ext_binv_SF_mc.GetBinContent(1); mll_ext_binv_SF_mc_e = mll_ext_binv_SF_mc.GetBinError(1);     
    mll_bveto_OF_da_int = mll_bveto_OF_da.GetBinContent(1); mll_bveto_OF_da_e = mll_bveto_OF_mc.GetBinError(1);
    mll_binv_OF_da_int = mll_binv_OF_da.GetBinContent(1); mll_binv_OF_da_e = mll_binv_OF_mc.GetBinError(1);     
    mll_ext_bveto_OF_da_int = mll_ext_bveto_OF_da.GetBinContent(1); mll_ext_bveto_OF_da_e = mll_ext_bveto_OF_da.GetBinError(1);
    mll_ext_binv_OF_da_int =  mll_ext_binv_OF_da.GetBinContent(1); mll_ext_binv_OF_da_e = mll_ext_binv_OF_da.GetBinError(1);     
    
    fmll_bveto_OF_mc,  fmll_bveto_OF_mc_e = getFraction(mll_bveto_OF_mc_int, mll_bveto_OF_mc_e, mll_ext_bveto_OF_mc_int, mll_ext_bveto_OF_mc_e)
    fmll_bveto_SF_mc,  fmll_bveto_SF_mc_e = getFraction(mll_bveto_SF_mc_int, mll_bveto_SF_mc_e, mll_ext_bveto_SF_mc_int, mll_ext_bveto_SF_mc_e)
    fmll_binv_OF_mc ,  fmll_binv_OF_mc_e  = getFraction(mll_binv_OF_mc_int,  mll_binv_OF_mc_e, mll_ext_binv_OF_mc_int, mll_ext_binv_OF_mc_e)     
    fmll_bveto_OF_da,  fmll_bveto_OF_da_e = getFraction(mll_bveto_OF_da_int, mll_bveto_OF_da_e, mll_ext_bveto_OF_da_int, mll_ext_bveto_OF_da_e)
    fmll_binv_OF_da ,  fmll_binv_OF_da_e  = getFraction(mll_binv_OF_da_int,  mll_binv_OF_da_e, mll_ext_binv_OF_da_int, mll_ext_binv_OF_da_e)     

    print "fmll_bveto_OF_mc ",  fmll_bveto_OF_mc
    print "fmll_bveto_SF_mc ",  fmll_bveto_SF_mc
    print "fmll_binv_OF_mc  ",  fmll_binv_OF_mc 
    print "fmll_bveto_OF_da ",  fmll_bveto_OF_da
    print "fmll_binv_OF_da  ",  fmll_binv_OF_da  
    c1 = r.TCanvas()
    l1 = r.TLegend(0.4,0.8,0.6,0.89)
    
    fmll = r.TH1F('fmll_mc','f_{mll}',3,0,3)
    fmll.SetBinContent(1, fmll_binv_OF_mc);fmll.SetBinError(1,fmll_binv_OF_mc_e);fmll.GetXaxis().SetBinLabel(1, 'OF #geq 1b');
    fmll.SetBinContent(2, fmll_bveto_OF_mc);fmll.SetBinError(2,fmll_bveto_OF_mc_e);fmll.GetXaxis().SetBinLabel(2, 'OF = 0b');
    fmll.SetBinContent(3, fmll_bveto_SF_mc);fmll.SetBinError(3,fmll_bveto_SF_mc_e);fmll.GetXaxis().SetBinLabel(3, 'SF = 0b'); 
    fmll.SetMarkerColor(r.kBlue+2)
    fmll.SetLineColor  (r.kBlue-7)
    fmll.SetLineWidth  (2)
    r.gStyle.SetPaintTextFormat('.3f')
    fmll.Draw('e1')
    fmll.Draw('text00 same')
    fmll_da = fmll.Clone('data fracs')
    l1.AddEntry(fmll, 'MC', 'pl')
    fmll.SetBinContent(1, fmll_binv_OF_da);fmll.SetBinError(1,fmll_binv_OF_da_e);fmll.GetXaxis().SetBinLabel(1, 'OF #geq 1b');
    fmll.SetBinContent(2, fmll_bveto_OF_da);fmll.SetBinError(2,fmll_bveto_OF_da_e);fmll.GetXaxis().SetBinLabel(2, 'OF = 0b');
    #fmll.GetYaxis().SetRangeUser(0.10, 0.60)
    fmll.SetMinimum(0.10)
    fmll.SetMaximum(0.60)
    fmll.GetYaxis().SetTitle('f_{mll}')
    #fmll.GetYaxis().SetTitleOffset(1.2*fmll.GetYaxis().GetTitleOffset())
    fmll_da.SetMarkerStyle(22)
    fmll_da.SetMarkerColor(r.kBlack)
    fmll_da.SetLineColor  (r.kGray+1)
    fmll_da.SetLineWidth  (2)
    fmll_da.Draw('e1 ')
    fmll_da.Draw('text00 same ')
    fmll.Draw('e1 same')
    fmll.Draw('text00 same')
    l1.AddEntry(fmll_da, 'data', 'pl')
    l1.Draw('same')
    lat1 = r.TLatex(); lat1.SetNDC(); lat1.SetTextFont(fmll.GetYaxis().GetTitleFont()); lat1.SetTextSize(fmll.GetYaxis().GetTitleSize())
    lat1.DrawLatex(0.76, 0.93, '%.2f fb^{-1}'%(lint))
    lat1.DrawLatex (0.1, 0.9, 'fmll')
    c1.SaveAs('plots/ewino/%s/factors/fmll.png' %(lint_str))                                                                                                                                  

    r0b1b_mll_ext_OF_mc,  r0b1b_mll_ext_OF_mc_e = getFraction(mll_ext_bveto_OF_mc_int, mll_ext_bveto_OF_mc_e, mll_ext_binv_OF_mc_int, mll_ext_binv_OF_mc_e)
    r0b1b_mll_OF_mc,      r0b1b_mll_OF_mc_e =     getFraction(mll_bveto_OF_mc_int, mll_bveto_OF_mc_e, mll_binv_OF_mc_int, mll_binv_OF_mc_e)
    r0b1b_mll_ext_SF_mc,  r0b1b_mll_ext_SF_mc_e = getFraction(mll_ext_bveto_SF_mc_int, mll_ext_bveto_SF_mc_e, mll_ext_binv_SF_mc_int, mll_ext_binv_SF_mc_e)
    r0b1b_mll_ext_OF_da,  r0b1b_mll_ext_OF_da_e = getFraction(mll_ext_bveto_OF_da_int, mll_ext_bveto_OF_da_e, mll_ext_binv_OF_da_int, mll_ext_binv_OF_da_e)
    r0b1b_mll_OF_da,      r0b1b_mll_OF_da_e =     getFraction(mll_bveto_OF_da_int, mll_bveto_OF_da_e, mll_binv_OF_da_int, mll_binv_OF_da_e)

    c2 = r.TCanvas()
    l2 = r.TLegend(0.4,0.8,0.6,0.89)
    
    r0b1b = r.TH1F('r0b1b_mc','f_{mll}',3,0,3)
    r0b1b.SetBinContent(1, r0b1b_mll_ext_OF_mc);r0b1b.SetBinError(1,r0b1b_mll_ext_OF_mc_e);r0b1b.GetXaxis().SetBinLabel(1, 'OF mll 61-121');
    r0b1b.SetBinContent(2, r0b1b_mll_OF_mc);r0b1b.SetBinError(2,r0b1b_mll_OF_mc_e);r0b1b.GetXaxis().SetBinLabel(2, 'OF on-Z');
    r0b1b.SetBinContent(3, r0b1b_mll_ext_SF_mc);r0b1b.SetBinError(3,r0b1b_mll_ext_SF_mc_e);r0b1b.GetXaxis().SetBinLabel(3, 'SF mll 61-121'); 
    r0b1b.SetMarkerColor(r.kRed+2)
    r0b1b.SetLineColor  (r.kRed-7)
    r0b1b.SetLineWidth  (2)
    r.gStyle.SetPaintTextFormat('.3f')
    r0b1b.Draw('e1')
    r0b1b.Draw('text00 same')
    r0b1b_da = r0b1b.Clone('data fracs')
    l2.AddEntry(r0b1b, 'MC', 'pl')
    r0b1b.SetBinContent(1, r0b1b_mll_ext_OF_da);r0b1b.SetBinError(1,r0b1b_mll_ext_OF_da_e);r0b1b.GetXaxis().SetBinLabel(1, 'OF mll 61-121');
    r0b1b.SetBinContent(2, r0b1b_mll_OF_da);r0b1b.SetBinError(2,r0b1b_mll_OF_da_e);r0b1b.GetXaxis().SetBinLabel(2, 'OF on-Z');
    r0b1b.GetYaxis().SetRangeUser(0.10, 0.60)
    r0b1b.SetMinimum(0.10)
    r0b1b.SetMaximum(0.60)
    r0b1b.GetYaxis().SetTitle('r_{0b1b}')
    #r0b1b.GetYaxis().SetTitleOffset(1.2*r0b1b.GetYaxis().GetTitleOffset())
    r0b1b_da.SetMarkerStyle(22)
    r0b1b_da.SetMarkerColor(r.kBlack)
    r0b1b_da.SetLineColor  (r.kGray+1)
    r0b1b_da.SetLineWidth  (2)
    r0b1b_da.Draw('e1 ')
    r0b1b_da.Draw('text00 same ')
    r0b1b.Draw('e1 same')
    r0b1b.Draw('text00 same')
    l2.AddEntry(r0b1b_da, 'data', 'pl')
    l2.Draw('same')
    lat2 = r.TLatex(); lat2.SetNDC(); lat2.SetTextFont(r0b1b.GetYaxis().GetTitleFont()); lat2.SetTextSize(r0b1b.GetYaxis().GetTitleSize())
    lat2.DrawLatex(0.76, 0.93, '%.2f fb^{-1}'%(lint))
    lat2.DrawLatex (0.1, 0.9, 'r_{0b1b}')
    c2.SaveAs('plots/ewino/%s/factors/r0b1b%s.png' %(lint_str, ''))         
    saveInFileFactors("ingredients.dat", r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e,  r0b1b_mll_ext_OF_da, r0b1b_mll_ext_OF_da_e, "r0b1b", "rb_mll_ext_OF")
    saveInFileFactors("ingredients.dat", r0b1b_mll_OF_mc, r0b1b_mll_OF_mc_e,  r0b1b_mll_OF_da, r0b1b_mll_OF_da_e, "r0b1b", "rb_mll_reg_OF")
    saveInFileFactors("ingredients.dat", r0b1b_mll_ext_SF_mc, r0b1b_mll_ext_SF_mc_e,  0.0000, 0.0000,  "r0b1b","rb_mll_ext_SF")
    saveInFileFactors("ingredients.dat", fmll_bveto_OF_mc, fmll_bveto_OF_mc_e,  fmll_bveto_OF_da, fmll_bveto_OF_da_e,"fmll_",  "fmll_bveto_OF")
    saveInFileFactors("ingredients.dat", fmll_bveto_SF_mc, fmll_bveto_SF_mc_e,  0.0000, 0.0000,"fmll_",  "fmll_bveto_SF")
    saveInFileFactors("ingredients.dat", fmll_binv_OF_mc, fmll_binv_OF_mc_e,  fmll_binv_OF_da, fmll_binv_OF_da_e, "fmll_", "fmll_binv__OF")
    makeEWKFactorsTable(fmll_bveto_OF_mc, fmll_bveto_OF_mc_e,  fmll_binv_OF_mc, fmll_binv_OF_mc_e , fmll_bveto_SF_mc, fmll_bveto_SF_mc_e, fmll_bveto_OF_da, fmll_bveto_OF_da_e, fmll_binv_OF_da, fmll_binv_OF_da_e , 0.0000, 0.0000, 'fmll' )    
    makeEWKFactorsTable(r0b1b_mll_OF_mc, r0b1b_mll_OF_mc_e,  r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e , r0b1b_mll_ext_SF_mc, r0b1b_mll_ext_SF_mc_e, r0b1b_mll_OF_da, r0b1b_mll_OF_da_e,  r0b1b_mll_ext_OF_da, r0b1b_mll_ext_OF_da_e , 0.0000, 0.0000,'r0b1b' )    
    return fmll_binv_OF_mc, fmll_binv_OF_mc_e, r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e

def makeClosureTests(var, specialcut = '', scutstring = '', doCumulative = False):

    if var == 'mll':
        treevar = 'lepsMll_Edge'
        nbins, xmin, xmax = 28, 20, 300
        xlabel = 'm_{ll} (GeV)'              
    elif var == 'met':
        treevar = 'met_Edge'
        nbins, xmin, xmax = 15, 0, 300
        xlabel = 'E_{T}^{miss} [GeV]'                    
    lint = 12.9  ; maxrun = 999999; lint_str = '12.9invfb'
    makeFactors = False 
    if makeFactors:
        fmll_binv_OF_mc, fmll_binv_OF_mc, r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e = makeTheFactors()
    else:
        fmll_binv_OF_mc, fmll_binv_OF_mc_e         = helper.readFromFileEWKFactors("ingredients.dat", "MC", "fmll_", "fmll_binv__OF")
        r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e = helper.readFromFileEWKFactors("ingredients.dat", "MC", "r0b1b", "rb_mll_ext_OF") 
    ## mll distributions
    mc_OF = treeTT.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.ewinoClosureExt,  cuts.OF]), '', xlabel)
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.ewinoClosureExt,  cuts.OF]), '', xlabel)
    print "mc_OF", mc_OF.Integral()
    print "da_OF", da_OF.Integral()
    print "mc_OF_corr", mc_OF.Integral()*fmll_binv_OF_mc*r0b1b_mll_ext_OF_mc
    print "da_OF_corr", da_OF.Integral()*fmll_binv_OF_mc*r0b1b_mll_ext_OF_mc
    mc_SF = treeTT.getTH1F(lint, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.ewinoClosureReg,  cuts.SF]), '', xlabel)
    print "mc_SF", mc_SF.Integral()
    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_err.SetFillStyle(3004); mc_OF_err.SetMarkerSize(0.)
    da_OF_err = copy.deepcopy(da_OF)
    da_OF_err.SetFillColorAlpha(r.kBlack+1, 0.8)
    da_OF_err.SetFillStyle(3004); da_OF_err.SetMarkerSize(0.)

    mc_SF.GetYaxis().SetRangeUser(0., 2*mc_SF.GetMaximum())
    print helper.bcolors.HEADER + '[MC only closure test not scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC

    mc_OF_rsfofScaled = copy.deepcopy(mc_OF)
    mc_OF_rsfofScaled = scaleByEWKFactors(mc_OF_rsfofScaled, rsfof_mc, rsfof_mc_e, fmll_binv_OF_mc, fmll_binv_OF_mc, r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e)
    mc_OF_rsfofScaled_err = copy.deepcopy(mc_OF_rsfofScaled)
    mc_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_rsfofScaled_err.SetFillStyle(3004); mc_OF_rsfofScaled_err.SetMarkerSize(0.)                                                                                    
    da_OF_rsfofScaled = copy.deepcopy(da_OF)
    da_OF_rsfofScaled = scaleByEWKFactors(da_OF_rsfofScaled, rsfof_mc, rsfof_mc_e, fmll_binv_OF_mc, fmll_binv_OF_mc, r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e)
    da_OF_rsfofScaled_err = copy.deepcopy(da_OF_rsfofScaled)
    da_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlack, 0.8)
    da_OF_rsfofScaled_err.SetFillStyle(3004); da_OF_rsfofScaled_err.SetMarkerSize(0.)                                                                                    

    print helper.bcolors.HEADER + '[MC only closure test scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure = Canvas.Canvas('ewino/closure/%s/plot_closure_%s_mcPredmcObs%s'%(lint_str, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.55, 0.75, 0.8)
    plot_closure.addHisto(mc_SF                , 'PE'       , 'MC - SF', 'PL', r.kRed+1  , 1,  0)
    plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_closure.addHisto(da_OF_rsfofScaled    , 'hist,SAME', 'Data - OF', 'L' , r.kBlack , 1,  1)
    plot_closure.addLatex (0.61, 0.82, 'R_{SFOF} scaled')
    plot_closure.saveRatio(1, 0, 0, lint, mc_SF, mc_OF_rsfofScaled, 0.2, 1.8)
    makeClosureTable(mc_SF, mc_OF_rsfofScaled, da_OF_rsfofScaled)


def makeResultData(var, maxrun = 999999, lint = 12.9, specialcut = '', scutstring = '', _options = ''):
    
    makeFactors = False 
    if makeFactors:
        fmll_binv_OF_mc, fmll_binv_OF_mc, r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e = makeTheFactors()
    else:
        fmll_binv_OF_mc, fmll_binv_OF_mc_e         = helper.readFromFileEWKFactors("ingredients.dat", "MC", "fmll_", "fmll_binv__OF")
        r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e = helper.readFromFileEWKFactors("ingredients.dat", "MC", "r0b1b", "rb_mll_ext_OF") 
    returnplot, addRares, splitFlavor, makeTable, printIntegral = False, True, False, False, False
    if 'returnplot'    in _options: print 'found option %s'%'returnplot'    ;returnplot    = True
    if 'splitFlavor'   in _options: print 'found option %s'%'splitFlavor'   ;splitFlavor   = True
    if 'makeTable'     in _options: print 'found option %s'%'makeTable'     ;makeTable     = True
    if 'printIntegral' in _options: print 'found option %s'%'printIntegral' ;printIntegral = True

    # if   var == 'mll'      : treevar = 'lepsMll_Edge'        ; nbins, xmin, xmax = 23, 20 , 250 ; xlabel = 'm_{ll} (GeV)'
    if   var == 'mll'      : treevar = 'lepsMll_Edge'        ; nbins, xmin, xmax = 30, 20 , 420 ; xlabel = 'm_{ll} [GeV]'
    elif var == 'zpt'      : treevar = 'lepsZPt_Edge'        ; nbins, xmin, xmax = 10,  0 ,1000 ; xlabel = 'p_{T}^{ll}'
    elif var == 'mlb'      : treevar = 'sum_mlb_Edge'        ; nbins, xmin, xmax = 15,  0 ,1500 ; xlabel = '#Sigma m_{lb}'
    elif var == 'met'      : treevar = 'met_Edge'            ; nbins, xmin, xmax = 18, 100, 1000 ; xlabel = 'E_{T}^{miss.} [GeV]'
        
    newLumiString = str(lint)+'invfb'
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.ewinoClosureReg,cuts.MET150,  cuts.SF]), '', xlabel)
    print "da_SF ", da_SF.Integral()
    da_mm = treeDA.getTH1F(lint, var+"da_mm"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.ewinoClosureReg,cuts.MET150,  cuts.mm]), '', xlabel)
    da_ee = treeDA.getTH1F(lint, var+"da_ee"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.ewinoClosureReg,cuts.MET150,  cuts.ee]), '', xlabel)
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.ewinoClosureExt,cuts.MET150,  cuts.OF]), '', xlabel)
    print "da_OF ", da_OF.Integral()
    print "fmll_binv_OF_mc*r0b1b_mll_ext_OF_mc", fmll_binv_OF_mc*r0b1b_mll_ext_OF_mc
    print "da_OF*fmll_binv_OF_mc*r0b1b_mll_ext_OF_mc ", da_OF.Integral()*fmll_binv_OF_mc*r0b1b_mll_ext_OF_mc
    ra_OF = treeRA.getTH1F(lint, var+"ra_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.ewinoClosureReg, cuts.MET150, cuts.OF]), '', xlabel)
    ra_SF = treeRA.getTH1F(lint, var+"ra_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.ewinoClosureReg, cuts.MET150, cuts.SF]), '', xlabel)
    ra_SF.SetFillColorAlpha(r.kCyan+1, 0.8)
    print "rare_OF ", ra_OF.Integral()
    print "rare_SF ", ra_SF.Integral()
    ra_SF.SetFillStyle(3017); ra_SF.SetMarkerSize(0.)
    ra_SF.Add(ra_OF, -1.)                                                                                                                                                                      
    print "rare_OF_corr ", ra_SF.Integral()
    
    da_OF_err = copy.deepcopy(da_OF)
    da_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    da_OF_err.SetFillStyle(3017); da_OF_err.SetMarkerSize(0.)

    da_OF_rsfofScaled = copy.deepcopy(da_OF)
    da_OF_rsfofScaled.SetName(da_OF.GetName()+'_scaledRSFOF')
    da_OF_rsfofScaled = scaleByEWKFactors(da_OF_rsfofScaled, rsfof, rsfof_e, fmll_binv_OF_mc, fmll_binv_OF_mc, r0b1b_mll_ext_OF_mc, r0b1b_mll_ext_OF_mc_e)
    da_OF_rsfofScaled.SetLineWidth(2)
    print "da_OF_corr ", da_OF_rsfofScaled.Integral()

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
    plot_result = Canvas.Canvas('ewino/results/%s/plot_result_%s_daPreddaObs%s'%(newLumiString, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.60, 0.65, 0.80, 0.85)
    plot_result.addHisto(da_OF_rsfofScaled_err, 'e2,same'  , ''         , 'PL', r.kBlue+1 , 1, -1)
    plot_result.addHisto(da_OF_rsfofScaled    , 'hist,SAME', 'predicted', 'L' , r.kBlue+1 , 1,  1)
    plot_result.addHisto(da_ee, 'pe,same', 'data- elel'  , 'PL', r.kYellow+1  , 1,  2)
    plot_result.addHisto(da_mm, 'pe,same', 'data- #mu#mu', 'PL', r.kGreen+1   , 1,  3)
    plot_result.addHisto(da_SF                , 'PE,same'    , 'observed data', 'PL', r.kBlack  , 1,  0)
    plot_result.addHisto(ra_SF                , 'hist, SAME'    , 'rares SF', 'PL', r.kCyan+1  , 1,  0)
    #plot_result.addLatex (0.2, 0.80, 'R_{SFOF} scaled')
    if printIntegral  : plot_result.addLatex (0.2, 0.65, 'SF: {SF:.1f} (#mu#mu:{mm:.1f},ee:{ee:.1f})'.format(SF=da_SF.Integral(),mm=da_mm.Integral(),ee=da_ee.Integral() ) )
    if printIntegral  : plot_result.addLatex (0.2, 0.60, 'OF: {OF:.2f}'.format(OF=da_OF_rsfofScaled.Integral() ) )
    plot_result.saveRatio(1, 1, 0, lint, da_SF, da_OF_rsfofScaled, 0. , int(maxrat+1.0) )

    if makeTable:
        makeSimpleTable(plot_result, addRares)
    if returnplot:
        return plot_result

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
    ttDatasets = ['TTJets_DiLepton', 'TTJets_DiLepton_ext']
    raDatasets = ['ZZTo4L', 'WZTo3LNu', 'WWW', 'WWZ','WZZ', 'ZZZ',  'TTZToLLNuNu' ,'TTWToLNu', 'T_tWch', 'TBar_tWch' , 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TToLeptons_sch', 'TToLeptons_tch_powheg', 'TBarToLeptons_tch_powheg', 'TTHnobb_pow', 'VHToNonbb', 'WJetsToLNu_LO']
    mcDatasets = ttDatasets+ ([] if opts.onlyTTbar else dyDatasets + raDatasets)
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_273150_275376', 'DoubleEG_Run2016B-PromptReco-v2_runs_273150_275376', 'MuonEG_Run2016B-PromptReco-v2_runs_273150_275376',
                  'DoubleMuon_Run2016C-PromptReco-v2_runs_275420_276283', 'DoubleEG_Run2016C-PromptReco-v2_runs_275420_276283', 'MuonEG_Run2016C-PromptReco-v2_runs_275420_276283',
                  'DoubleMuon_Run2016D-PromptReco-v2_runs_276315_276811', 'DoubleEG_Run2016D-PromptReco-v2_runs_276315_276811', 'MuonEG_Run2016D-PromptReco-v2_runs_276315_276811']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0)
    treeRA = Sample.Tree(helper.selectSamples(opts.sampleFile, raDatasets, 'RA'), 'RA'  , 0)
    treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT'  , 0)
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
 

 #   #for v in ['mll']:#'nll_noMET', 'nll_noMLB', 'nll_noZPT', 'nll_noLDP']: #'iso1', 'iso2', 'mll', 'nll', 'nb', 'nj', 'zpt', 'mlb', 'met', 'ldp', 'pt1', 'pt2']:
 #   #    makeResultData(v , maxrun , lint , specialcut =''                                     , scutstring = ''                 , _options='returnplot,splitFlavor,printIntegral')
 #   #    makeResultData(v , maxrun , lint , specialcut ='nll_Edge > 21.'                        , scutstring = 'nllAbove21'        , _options='returnplot,splitFlavor,printIntegral')
 #   #    makeResultData(v , maxrun , lint , specialcut ='nll_Edge <= 21.'                        , scutstring = 'nllBelow21'        , _options='returnplot,splitFlavor,printIntegral')
 #       
 #       makeResultData(v,maxrun,lint,specialcut='nll_Edge <= 21. && lepsMll_Edge > 400' , scutstring = 'mll400_nllBelow21', _options='returnplot,splitFlavor,printIntegral')


    makeClosureTests('met','', 'inclusive', True)
    
    #makeResultData('met', maxrun, lint, specialcut = '' , scutstring = '')

