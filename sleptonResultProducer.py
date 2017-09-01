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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats, TStyle, TH1, SetOwnership
import math, sys, optparse, copy, re, array, subprocess

import include.nll      
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


#def makeFactorsTable(): ## for this to make sense the region should be properly binned!!
#    rmue_da = helper.readFromFileRmue("ingredients.dat", "DATA") 
#    rmue_mc = helper.readFromFileRmue("ingredients.dat", "MC") 
#    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA","DATA")
#    rmue_a_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffA","MC")
#    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB","DATA")
#    rmue_b_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffB","MC")
#    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
#    rsfof_mc = helper.readFromFileRsfofD("ingredients.dat", "MC")   
#    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
#    rt_mc = helper.readFromFileRT("ingredients.dat", "MC")  
#    rsfof_fac_da, rsfof_fac_da_e = getFactor(rmue_da[0], rmue_da[1],rmue_da[2], rt_da[0], rt_da[1], rt_da[2]) 
#    rsfof_fac_mc, rsfof_fac_mc_e = getFactor(rmue_mc[0], rmue_mc[1],rmue_mc[2], rt_da[0], rt_da[1], rt_da[2]) 
#    line0 = ' & Data & MC '
#    line05 = '\\hline '
#    line1 = 'rmue  &  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f \\\  ' %(rmue_da[0],getFinalError(rmue_da[1], rmue_da[2]), rmue_mc[0],getFinalError(rmue_mc[1], rmue_mc[2]))
#    line2 = 'RT    &  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f  \\\ ' %(rt_da[0],getFinalError(rt_da[1], rt_da[2]), rt_mc[0],getFinalError(rt_mc[1], rt_mc[2]))
#    line25 = '\\hline '
#    line3 = 'RSFOF factorization &  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f  \\\ ' %(rsfof_fac_da, rsfof_fac_da_e, rsfof_fac_mc, rsfof_fac_mc_e)
#    line4 = 'RSFOF direct        &  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f  \\\ ' %(rsfof_da[0],getFinalError(rsfof_da[1], rsfof_da[2]), rsfof_mc[0],getFinalError(rsfof_mc[1], rsfof_mc[2]))
#    line5 = 'RSFOF weighted average&  %.3f $\\pm$ %.3f  & %.3f $\\pm$ %.3f  \\\ ' %((rsfof_da[0]+rsfof_fac_da)/2, getFinalError(rsfof_fac_da_e,  getFinalError(rsfof_da[1], rsfof_da[2])), (rsfof_mc[0]+rsfof_fac_mc)/2, getFinalError( rsfof_fac_mc_e, getFinalError(rsfof_mc[1], rsfof_mc[2])))
#    line55 = '\\hline'
#    print line0                                                                                                                                                      
#    print line05                                                                                                                                                      
#    print line1                                                                                                                                                      
#    print line2                                                                                                                                                      
#    print line25                                                                                                                                                      
#    print line3                                                                                                                                                      
#    print line4                                                                                                                                                      
#    print line5                                                                                                                                                      
#    print line55 
#    saveInFile("ingredients.dat", rsfof_fac_mc, rsfof_fac_mc_e, rsfof_fac_da, rsfof_fac_da_e)
#    return (rsfof_da[0]+rsfof_fac_da)/2, getFinalError(rsfof_fac_da_e,  getFinalError(rsfof_da[1], rsfof_da[2])), (rsfof_mc[0]+rsfof_fac_mc)/2, getFinalError(rsfof_fac_mc_e,  getFinalError(rsfof_mc[1], rsfof_mc[2])) , rmue_a_da, rmue_a_mc, rmue_b_da, rmue_b_mc                                                                                   
#


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
        h_rsfof.SetBinError  (i,rsfof_err)
    histo.Multiply(h_rsfof)
    return histo

def makeResultsTableSig(da, fs, zz, wz, rare, mc_full, sig1, sig2, sig3,signames):
    
    line0 = '\\begin{tabular}{c c c c}  '                                                               
    line1 = '\\\\ \hline\hline   '                                                               
    line2 = ' p_{T}^{miss} & 100-175& 175-250 & 250+  \\\\ \hline  '                                                               
    metLabels = ['100-175 ', '175-250 ','+250']
    lines = []
    print line0                                                                                                                                                      
    print line1                                                                                                                                                      
    print line2                                                                                                                                                      
    print 'FS &   %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  \\\  ' %(fs.GetBinContent(1),fs.GetBinError(1),fs.GetBinContent(2),fs.GetBinError(2),fs.GetBinContent(3),fs.GetBinError(3))
    print 'ZZ &   %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  \\\  ' %(zz.GetBinContent(1),zz.GetBinError(1),zz.GetBinContent(2),zz.GetBinError(2),zz.GetBinContent(3),zz.GetBinError(3))
    print 'WZ &   %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  \\\  ' %(wz.GetBinContent(1),wz.GetBinError(1),wz.GetBinContent(2),wz.GetBinError(2),wz.GetBinContent(3),wz.GetBinError(3))
    print 'rare & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  \\\\ \hline ' %(rare.GetBinContent(1),rare.GetBinError(1),rare.GetBinContent(2),rare.GetBinError(2),rare.GetBinContent(3),rare.GetBinError(3))
    print 'Total Pred.& %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  \\\\ \hline  ' %(mc_full.GetBinContent(1),mc_full.GetBinError(1),mc_full.GetBinContent(2),mc_full.GetBinError(2),mc_full.GetBinContent(3),mc_full.GetBinError(3))
    print 'm_{l}: %s, m_{{\Chi_{1}^{0}}} : %s &  %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  \\\  ' %( signames[0][0], signames[0][1],   sig1.GetBinContent(1),sig1.GetBinError(1),sig1.GetBinContent(2),sig1.GetBinError(2),sig1.GetBinContent(3),sig1.GetBinError(3))
    print 'm_{l}: %s, m_{{\Chi_{1}^{0}}} : %s &  %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  \\\  ' %( signames[1][0], signames[1][1],   sig2.GetBinContent(1),sig2.GetBinError(1),sig2.GetBinContent(2),sig2.GetBinError(2),sig2.GetBinContent(3),sig2.GetBinError(3))
    print 'm_{l}: %s, m_{{\Chi_{1}^{0}}} : %s &  %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  \\\  ' %( signames[2][0], signames[2][1],   sig3.GetBinContent(1),sig3.GetBinError(1),sig3.GetBinContent(2),sig3.GetBinError(2),sig3.GetBinContent(3),sig3.GetBinError(3))
    print '\\hline \hline'
    print '\end{tabular}'                                                                                                                                                         

def makeResultsTable(da, fs, zz, wz ,rare, mc):
    
    line0 = '\\begin{tabular}{c c c c c c c}  '                                                               
    line2 = ' p_{T}^{miss} & FS & ZZ & WZ & Rares & Total & Data \\\\ \hline  '                                                               
    line1 = '\\\\ \hline\hline   '                                                               
    metLabels = ['100-175 ', '175-250 ','+250']
    lines = []
    print line0                                                                                                                                                      
    print line1                                                                                                                                                      
    print line2                                                                                                                                                      
    for bin,metLabel in enumerate(metLabels):
        if metLabel == 'onz': continue
        print '%s &   %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f  & %.1f $\\pm$ %.1f   & %.1f $\\pm$ %.1f   & %.1f $\\pm$ %.1f   & blinded!  \\\  ' %(metLabel,
                                                                                                                      fs.GetBinContent(bin+1),
                                                                                                                      fs.GetBinError(bin+1),
                                                                                                                      zz.GetBinContent(bin+1),
                                                                                                                      zz.GetBinError(bin+1),
                                                                                                                      wz.GetBinContent(bin+1),
                                                                                                                      wz.GetBinError(bin+1),
                                                                                                                      rare.GetBinContent(bin+1),
                                                                                                                      rare.GetBinError(bin+1),
                                                                                                                      mc.GetBinContent(bin+1),
                                                                                                                      mc.GetBinError(bin+1) ,
                                                                                                                      #da.GetBinContent(bin+1)
                                                                                                                      )

    print '\\hline \hline'
    print '\end{tabular}'                                                                                                                                                         



def makeClosureTests(var, specialcut = '', scutstring = '', doCumulative = False, nbins=0, xmin=0, xmax=0, save=True):

    if var == 'met':
        treevar = 'met_Edge'
        nbins, xmin, xmax = 29, 16, 306
        xlabel = 'm_{ll} (GeV)'
    else: 
        treevar = var
        xlabel = ''
    
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rsfof_mc = helper.readFromFileRsfofD("ingredients.dat", "MC")   
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rt_mc = helper.readFromFileRT("ingredients.dat", "MC")  
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "DATA")
    rmue_a_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "MC")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "DATA")
    rmue_b_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "MC")

                                                                                                                                
    mc_OF = treeFS.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.OF]), '', xlabel,"1", kf)
    mc_SF = treeFS.getTH1F(lint, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel,"1", kf)
    dy_SF = treeDY.getTH1F(lint, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel,"1", kf)
    mc_OF_factor = treeFS.getTH1F(lint, var+"mc_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel, '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]),kf)
    mc_OF_factorUp = treeFS.getTH1F(lint, var+"mc_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]), kf)
    mc_OF_factorDn = treeFS.getTH1F(lint, var+"mc_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]), kf)

    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_direct = copy.deepcopy(mc_OF) 
    mc_OF_direct = scaleByRSFOF(mc_OF_direct, rsfof_mc[0], rsfof_mc[1])
    mc_OF_factor = getRMueError(mc_OF_factor, mc_OF_factorUp, mc_OF_factorDn)
    mc_OF_factor = scaleByRSFOF(mc_OF_factor, rt_mc[0], getFinalError(rt_mc[1],rt_mc[2]))
    result = weightedAverage( mc_OF_factor, mc_OF_direct, mc_OF)
    return result


def makeClosureTestPlots(analysis, var, specialcut = '', scutstring = '', doCumulative = False, nbins=0, xmin=0, xmax=0, save=True):
    treevar = 'met_Edge'        ; nbins = [100, 175, 250]; xmin =1 ; xmax = 1; xlabel = 'p_{T}^{miss} [GeV]'
    scan = Scans.Scan(analysis) 
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rsfof_mc = helper.readFromFileRsfofD("ingredients.dat", "MC")   
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rt_mc = helper.readFromFileRT("ingredients.dat", "MC")  
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "DATA")
    rmue_a_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "MC")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "DATA")
    rmue_b_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "MC")

    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline,  cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline,  cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    print '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0])
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue
    da_OF_direct = copy.deepcopy(da_OF)
    da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], rsfof_da[1])
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    da_result = weightedAverage( da_OF_factor, da_OF_direct, da_OF)                                   
    ### 
    da_prediction = copy.deepcopy( da_OF )
    da_OF.SetBinErrorOption( TH1.kPoisson)
    for i in range(1, da_prediction.GetNbinsX()+1):
        da_prediction.SetBinError(i,da_result[3].GetBinError(i))
    da_prediction.Multiply(da_result[2])                       
                                                                                                                                
    for bin, label in scan.SRLabels.items():
        bin = da_prediction.FindBin(bin)
        # to get asymmetric errors (i dont know why it doesnt work another way)
        dummyHisto = r.TH1F()
        dummyHisto.SetBinErrorOption(TH1.kPoisson)
        for i in range(1, int(da_OF.GetBinContent(bin))+1):
            dummyHisto.Fill(0.5)
        dummyHisto.GetBinContent(1), '+/-', dummyHisto.GetBinErrorUp(1), dummyHisto.GetBinErrorLow(1)
        syst = '  {value:4.1f}^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(value = da_prediction.GetBinContent(bin),
                                                                             errUp = da_prediction.GetBinErrorUp(bin),
                                                                             errDn = da_prediction.GetBinErrorLow(bin))
        stat = '^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(errUp = dummyHisto.GetBinErrorUp(1) * da_result[2].GetBinContent(bin),
                                                               errDn = dummyHisto.GetBinErrorLow(1) * da_result[2].GetBinContent(bin))
        del dummyHisto                                                                                                                 

    ## ## mll distributions
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.SF]), '', xlabel)
    mc_OF = treeFS.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel)
    mc_SF = treeFS.getTH1F(lint, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.SF]), '', xlabel)
    mc_ee = treeFS.getTH1F(lint, var+"mc_ee"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.ee]), '', xlabel)
    mc_mm = treeFS.getTH1F(lint, var+"mc_mm"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.mm]), '', xlabel)
    dy_SF = treeDY.getTH1F(lint, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.SF]), '', xlabel)
    mc_OF_factor = treeFS.getTH1F(lint, var+"mc_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))
    mc_OF_factorUp = treeFS.getTH1F(lint, var+"mc_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))
    mc_OF_factorDn = treeFS.getTH1F(lint, var+"mc_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))

    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_err.SetFillStyle(3004); mc_OF_err.SetMarkerSize(0.)
    dy_SF.SetFillColorAlpha(r.kGreen+2,0.5)

    mc_OF_direct = copy.deepcopy(mc_OF) 
    mc_OF_direct = scaleByRSFOF(mc_OF_direct, rsfof_mc[0], getFinalError(rsfof_mc[1], rsfof_mc[2]))
    mc_OF_factor = getRMueError(mc_OF_factor, mc_OF_factorUp, mc_OF_factorDn)
    mc_OF_factor = scaleByRSFOF(mc_OF_factor, rt_mc[0], getFinalError(rt_mc[1],rt_mc[2]))
    result = weightedAverage( mc_OF_factor, mc_OF_direct, mc_OF)
    mc_OF_rsfofScaled = result[3] 
    mc_OF_rsfofScaled_err = copy.deepcopy( mc_OF_rsfofScaled ) 
    mc_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_rsfofScaled_err.SetFillStyle(3004); mc_OF_rsfofScaled_err.SetMarkerSize(0.)          

    maxCont = max(da_prediction.GetMaximum(), mc_OF_rsfofScaled.GetMaximum())
    mc_SF.GetYaxis().SetRangeUser(0.1, 1.50*maxCont)
    mc_OF_rsfofScaled.GetYaxis().SetRangeUser(0.1, 1.50*maxCont)

    maxrat = 0.5
    for ib in range(1,da_SF.GetNbinsX()+1):
        tmp_rat = da_SF.GetBinContent(ib)/( mc_OF_rsfofScaled.GetBinContent(ib) if mc_OF_rsfofScaled.GetBinContent(ib) > 0 else 1. )
        if tmp_rat > maxrat:
            maxrat = tmp_rat                                                                                          
    mc_SF.GetYaxis().SetRangeUser(0.01, 1.5*mc_SF.GetMaximum())
    print helper.bcolors.HEADER + '[MC only closure test not scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC

    print helper.bcolors.HEADER + '[MC only closure test scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure = Canvas.Canvas('closure/%s/plot_closure_met_mcPredmcObs%s'%(lint_str,  '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_closure.addHisto(mc_SF                , 'PE'       , 'MC - SF', 'PL', r.kGreen+1  , 1,  0)
    plot_closure.addHisto(mc_ee                , 'PE, same'       , 'MC - ee', 'PL', r.kRed+1  , 1,  0)
    plot_closure.addHisto(mc_mm                , 'PE, same'       , 'MC - mm', 'PL', r.kCyan+1  , 1,  0)
    plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_closure.addLatex (0.61, 0.82, 'R_{SFOF} scaled')
    plot_closure.saveRatio(1, 0, 0, lint, mc_SF, mc_OF_rsfofScaled,  0. , 2, "SF/OF")                                                                                              
                                                                                                                                                                             
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
    prediction_final    = copy.deepcopy(unwght)
    
    for i in range(1, unwght.GetNbinsX()+1):
        if not unwght.GetBinContent(i): 
            rsfof_factor.SetBinContent( i, 0.)             
            rsfof_direct.SetBinContent( i, 0.) 
            rsfof_final .SetBinContent( i, 0.) 
            prediction_final   .SetBinContent( i, 0.)    
            continue 
        rsfof_factor.SetBinContent(i,factor.GetBinContent(i) / unwght.GetBinContent(i))
        rsfof_factor.SetBinError  (i, math.sqrt(abs(factor.GetBinError(i)**2 - (rsfof_factor.GetBinContent(i)*unwght.GetBinError(i))**2)) / unwght.GetBinContent(i) )
        rsfof_direct.SetBinContent(i, direct.GetBinContent(i) / unwght.GetBinContent(i))
        rsfof_direct.SetBinError  (i, math.sqrt(abs(direct.GetBinError(i)**2 - (rsfof_direct.GetBinContent(i)*unwght.GetBinError(i))**2)) / unwght.GetBinContent(i) )
        
        rsfof_final .SetBinContent(i, (rsfof_factor.GetBinContent(i) / rsfof_factor.GetBinError(i)**2 + rsfof_direct.GetBinContent(i) / rsfof_factor.GetBinError(i)**2) / ( 1. / rsfof_factor.GetBinError(i)**2 + 1. / rsfof_factor.GetBinError(i)**2 ))
        rsfof_final .SetBinError  (i, 1 / math.sqrt(1/rsfof_factor.GetBinError(i)**2 + 1/rsfof_direct.GetBinError(i)**2))
        
        prediction_final.SetBinContent(i, unwght.GetBinContent(i)*rsfof_final.GetBinContent(i))
        prediction_final.SetBinError  (i, math.sqrt( (unwght.GetBinError(i)*rsfof_final.GetBinContent(i))**2 + (unwght.GetBinContent(i)*rsfof_final.GetBinError(i))**2))

    return [rsfof_factor, rsfof_direct, rsfof_final, prediction_final]                                                                                                    


def makePred(scan, specialcut = '', scutstring = '', doCumulative = False, nbins=0, xmin=0, xmax=0):
    kf = "noKFactor"
    var = scan.srID
    if var == 'mll':
        treevar = 'lepsMll_Edge'
        xlabel = 'm_{ll} (GeV)'
    elif var == 'nll':
        treevar = 'nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge)'
        xlabel = 'NLL'
    elif var == 'nllMC':
        treevar = 'nll_mc_Edge'
        nbins, xmin, xmax = 26, 10, 36
        xlabel = 'NLL'
    elif var == 'nllMCSF':
        treevar = 'nll_mc_sf_Edge'
        nbins, xmin, xmax = 26, 10, 36
        xlabel = 'NLL'
    else: 
        treevar = var
        xlabel = ''
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA","DATA")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB","DATA")

    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel, "1", kf)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,'(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    print '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0])
    print 'trivial check'
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue

    da_OF_direct = copy.deepcopy(da_OF)
    da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], getFinalError(rsfof_da[1], rsfof_da[2]))
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    result = weightedAverage( da_OF_factor, da_OF_direct, da_OF)                                   
    ### 

    prediction = copy.deepcopy( da_OF )
    da_OF.SetBinErrorOption( TH1.kPoisson)
    for i in range(1, prediction.GetNbinsX()+1):
        prediction.SetBinError(i,result[3].GetBinError(i))
    prediction.Multiply(result[2])                       

    for bin, label in scan.SRLabels.items():
        bin = prediction.FindBin(bin)
        # to get asymmetric errors (i dont know why it doesnt work another way)
        dummyHisto = r.TH1F()
        dummyHisto.SetBinErrorOption(TH1.kPoisson)
        for i in range(1, int(da_OF.GetBinContent(bin))+1):
            dummyHisto.Fill(0.5)
        dummyHisto.GetBinContent(1), '+/-', dummyHisto.GetBinErrorUp(1), dummyHisto.GetBinErrorLow(1)
        print da_OF.GetBinContent(bin), '+', da_OF.GetBinErrorUp(bin), '-', da_OF.GetBinErrorLow(bin)
        syst = '  {value:4.1f}^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(value = prediction.GetBinContent(bin),
                                                                             errUp = prediction.GetBinErrorUp(bin),
                                                                             errDn = prediction.GetBinErrorLow(bin))
        stat = '^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(errUp = dummyHisto.GetBinErrorUp(1) * result[2].GetBinContent(bin),
                                                               errDn = dummyHisto.GetBinErrorLow(1) * result[2].GetBinContent(bin))
        del dummyHisto
        print label, syst, stat                                                                                                           
    
    return result                                                                                                                                                             

def makeRSOFTable(analysis):
    scan = Scans.Scan(analysis)
    rsfof_factor_mc, rsfof_direct_mc, rsfof_final_mc, nof_final_mc = makeClosureTests(scan.srID, specialcut = '', scutstring = '', doCumulative = False, nbins=scan.srIDMax+1, xmin=-0.5, xmax=scan.srIDMax+0.5,save=False)
    rsfof_factor_da, rsfof_direct_da, rsfof_final_da, prediction_final_da = makePred(scan, specialcut = '', scutstring = '', doCumulative = False, nbins=scan.srIDMax+1, xmin=-0.5, xmax=scan.srIDMax+0.5)
    print '                &  Data       & & & MC  & \\\ \\hline'
    print '     m_{ll}  &    $\mathrm{r_{SF/OF}^{fact}}$  & $\mathrm{r_{SF/OF}^{dict}}$  &  $\mathrm{r_{SF/OF}}$ &  $\mathrm{r_{SF/OF}^{fact}}$  & $\mathrm{r_{SF/OF}^{dict}}$  &  $\mathrm{r_{SF/OF}}$ \\\ \\hline'
    for bin, label in scan.SRLabels.items():
        print ' %s &      %4.3f $\pm$ %4.3f   & %4.3f $\pm$ %4.3f & %4.3f $\pm$ %4.3f & %4.3f $\pm$ %4.3f & %4.3f $\pm$ %4.3f & %4.3f $\pm$ %4.3f \\\\'%(label, 
                rsfof_factor_da.GetBinContent(rsfof_factor_da.FindBin(bin)), rsfof_factor_da.GetBinError(rsfof_factor_da.FindBin(bin)), 
                rsfof_direct_da.GetBinContent(rsfof_direct_da.FindBin(bin)), rsfof_direct_da.GetBinError(rsfof_direct_da.FindBin(bin)), 
                rsfof_final_da.GetBinContent(  rsfof_final_da.FindBin(bin)), rsfof_final_da.GetBinError(rsfof_final_da.FindBin(bin)), 
                rsfof_factor_mc.GetBinContent(rsfof_factor_mc.FindBin(bin)), rsfof_factor_mc.GetBinError(rsfof_factor_mc.FindBin(bin)), 
                rsfof_direct_mc.GetBinContent(rsfof_direct_mc.FindBin(bin)), rsfof_direct_mc.GetBinError(rsfof_direct_mc.FindBin(bin)), 
                rsfof_final_mc.GetBinContent(  rsfof_final_mc.FindBin(bin)), rsfof_final_mc.GetBinError(rsfof_final_mc.FindBin(bin)))


def makeResultData(analysis, var, signames, maxrun = 999999, lint = 35.9, specialcut = '', scutstring = '', _options = ''):
    returnplot, addRares, splitFlavor, makeTable, printIntegral = True, True, False, False, False
    if   var == 'met'      : treevar = 'met_Edge'        ; nbins = [100, 175, 250]; xmin =1 ; xmax = 1; xlabel = 'p_{T}^{miss} [GeV]'
    if not specialcut:
        specialcut = specialcut 
    else:
        specialcut = specialcut 
    scan = Scans.Scan(analysis)
    mc_stack = r.THStack()
    mc_stack_perGeV = r.THStack() 
    newLumiString = '35.9invfb'
    ##get the ingredients
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA","DATA")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB","DATA")
    #get the fs prediction
    sig1 = treeSIG1.getTH1F(lint, var+"sig1"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLeptonSignal, cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    sig2 = treeSIG2.getTH1F(lint, var+"sig2"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLeptonSignal, cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    sig3 = treeSIG3.getTH1F(lint, var+"sig3"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLeptonSignal, cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    sig1.SetLineStyle(2);sig2.SetLineStyle(2);sig3.SetLineStyle(2);
    sig1.SetLineWidth(2);sig2.SetLineWidth(2);sig3.SetLineWidth(2);
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slep0jet, cuts.OF]), '', xlabel, "1",kf)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slep0jet, cuts.OF]), '', xlabel,'(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton,cuts.slep0jet,  cuts.OF]), '', xlabel, '(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton,cuts.slep0jet,  cuts.OF]), '', xlabel, '(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    print '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0])
    print "da_OF ", da_OF.Integral()
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue
    da_OF_direct = copy.deepcopy(da_OF)
    da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], getFinalError(rsfof_da[1], rsfof_da[2]))
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    
    result = weightedAverage( da_OF_factor, da_OF_direct, da_OF) 




    ### 
    prediction = copy.deepcopy( da_OF )
    da_OF.SetBinErrorOption( TH1.kPoisson)
    for i in range(1, prediction.GetNbinsX()+1):
        prediction.SetBinError(i,result[3].GetBinError(i))
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
        print label, syst, stat                                                                                                           
    others = treeOTHERS.getTH1F(lint, var+"others"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton,  cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    zz_SF = treeZZ.getTH1F(lint, var+"zz_SF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1",kf)
    zz_SF_qcdUp = treeZZ.getTH1F(lint, var+"zz_SF_qcdUp"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "LHEweight_wgt_Edge*(LHEweight_id_Edge == 1005)", 'ZZpt')
    zz_SF_qcdDn = treeZZ.getTH1F(lint, var+"zz_SF_qcdDn"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "LHEweight_wgt_Edge*(LHEweight_id_Edge == 1009)", 'ZZpt')
    wz_SF = treeWZ.getTH1F(lint, var+"wz_SF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel,"1", 'ZZpt')
    wz_SF_qcdUp = treeWZ.getTH1F(lint, var+"wz_SF_qcdUp"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF,cuts.lepsFromZ]),'', xlabel, "LHEweight_wgt_Edge*(LHEweight_id_Edge == 1005)", 'ZZpt')
    wz_SF_qcdUp = treeWZ.getTH1F(lint, var+"wz_SF_qcdDn"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF,cuts.lepsFromZ]),'', xlabel, "LHEweight_wgt_Edge*(LHEweight_id_Edge == 1009)", 'ZZpt')
    #scale wz and zz according to SF 
    wz_SF.Scale(1.05);zz_SF.Scale(1.08)


    for bin, label in scan.SRLabels.items():
        if zz_SF.GetBinContent(bin) > 0:
            zz_SF.SetBinError(bin, zz_SF.GetBinContent(bin)*math.sqrt(((zz_SF.GetBinContent(bin)*0.1)/zz_SF.GetBinContent(bin))**2 + zz_SF.GetBinError(bin)/zz_SF.GetBinContent(bin)**2 + abs(zz_SF_qcdUp.GetBinContent(bin) - zz_SF.GetBinContent(bin)) + abs(zz_SF_qcdDn.GetBinContent(bin) - zz_SF.GetBinContent(bin))))
        else: zz_SF.SetBinError(bin , zz_SF.GetBinError(bin))                                                                                                                                    
        if wz_SF.GetBinContent(bin) > 0:
            wz_SF.SetBinError(bin, wz_SF.GetBinContent(bin)*math.sqrt(((wz_SF.GetBinContent(bin)*0.06)/wz_SF.GetBinContent(bin))**2 + wz_SF.GetBinError(bin)/wz_SF.GetBinContent(bin)**2 + abs(wz_SF_qcdUp.GetBinContent(bin) - wz_SF.GetBinContent(bin)) + abs(wz_SF_qcdDn.GetBinContent(bin) - wz_SF.GetBinContent(bin))  ))
        else: wz_SF.SetBinError(bin , wz_SF.GetBinError(bin))                                                                                                                                         
        if others.GetBinContent(bin) > 0:
            others.SetBinError(bin, others.GetBinContent(bin)*math.sqrt(((others.GetBinContent(bin)*0.5)/others.GetBinContent(bin))**2 + (others.GetBinError(bin)/others.GetBinContent(bin)**2)) )
        else: others.SetBinError(bin , others.GetBinError(bin))                                                                                                                                  
    
    # aesthetics
    
    rare = copy.deepcopy(others)
    rare.SetFillColorAlpha(r.kCyan-5, 1);rare.SetTitle("rares");rare.SetLineColor(r.kBlack)
    zz_SF.SetFillColorAlpha(r.kCyan+2, 1);zz_SF.SetTitle("ZZ");zz_SF.SetLineColor(r.kBlack)
    wz_SF.SetFillColorAlpha(r.kGreen-8, 1);wz_SF.SetTitle("WZ");wz_SF.SetLineColor(r.kBlack)
    prediction.SetFillColorAlpha(r.kRed-9, 0.7);prediction.SetTitle("FS");
    da_SF.SetTitle("data SF")       
    da_perGeV = copy.deepcopy(da_SF);fs_perGeV = copy.deepcopy(prediction);rare_perGeV = copy.deepcopy(rare) ;wz_SF_perGeV = copy.deepcopy(wz_SF);zz_SF_perGeV = copy.deepcopy(zz_SF)
    for ib in range(1,da_SF.GetNbinsX()+1):
        da_perGeV.SetBinContent(ib,  da_SF.GetBinContent(ib)/(da_SF.GetBinLowEdge(ib+1)-da_SF.GetBinLowEdge(ib)))
        fs_perGeV.SetBinContent(ib,  prediction.GetBinContent(ib)/(prediction.GetBinLowEdge(ib+1)-prediction.GetBinLowEdge(ib)))
        rare_perGeV.SetBinContent(ib, rare.GetBinContent(ib)/(rare.GetBinLowEdge(ib+1)-rare.GetBinLowEdge(ib)))                            
        wz_SF_perGeV.SetBinContent(ib, wz_SF.GetBinContent(ib)/(wz_SF.GetBinLowEdge(ib+1)-wz_SF.GetBinLowEdge(ib)))                            
        zz_SF_perGeV.SetBinContent(ib, zz_SF.GetBinContent(ib)/(zz_SF.GetBinLowEdge(ib+1)-zz_SF.GetBinLowEdge(ib)))                            
        #da_perGeV.SetBinError(ib,  da_SF.GetBinError(ib)/(da_SF.GetBinLowEdge(ib+1)-da_SF.GetBinLowEdge(ib)))
        #fs_perGeV.SetBinError(ib,  prediction.GetBinError(ib)/(prediction.GetBinLowEdge(ib+1)-prediction.GetBinLowEdge(ib)))
        #rare_perGeV.SetBinError(ib,  rare.GetBinError(ib)/(rare.GetBinLowEdge(ib+1)-rare.GetBinLowEdge(ib)))

    mc_stack.Add(rare      );                      mc_stack_perGeV.Add(rare_perGeV);
    mc_stack.Add(wz_SF      );                     mc_stack_perGeV.Add(wz_SF_perGeV);
    mc_stack.Add(zz_SF      );                     mc_stack_perGeV.Add(zz_SF_perGeV);
    mc_stack.Add(prediction);                      mc_stack_perGeV.Add(fs_perGeV); 
    mc_stack.Draw();                               mc_stack_perGeV.Draw("SAME");
    mc_stack.GetXaxis().SetTitle('p_{T}^{miss} [GeV]');  mc_stack_perGeV.GetXaxis().SetTitle('p_{T}^{miss} [GeV]');
    mc_full = copy.deepcopy(prediction);             mc_full_perGeV =  copy.deepcopy(fs_perGeV);
    mc_full.Add(rare, 1.);                         mc_full_perGeV.Add(rare_perGeV, 1.); 
    mc_full.Add(wz_SF, 1.);                        mc_full_perGeV.Add(wz_SF_perGeV, 1.); 
    mc_full.Add(zz_SF, 1.);                        mc_full_perGeV.Add(zz_SF_perGeV, 1.); 
    mc_full_e = copy.deepcopy(mc_full);            mc_full_perGeV_e = copy.deepcopy(mc_full_perGeV);
    mc_full_e.SetFillColorAlpha(r.kBlue+1, 0.8);mc_full_e.SetFillStyle(3017); mc_full_e.SetMarkerSize(0.);
    mc_full_perGeV_e.SetFillColorAlpha(r.kBlue+1, 0.8);mc_full_perGeV_e.SetFillStyle(3017); mc_full_perGeV_e.SetMarkerSize(0.);
    maxrat = 0.5
    for ib in range(1,da_SF.GetNbinsX()+1):
        tmp_rat = da_SF.GetBinContent(ib)/( mc_full.GetBinContent(ib) if mc_full.GetBinContent(ib) > 0 else 1. )
        if tmp_rat > maxrat:
            maxrat = tmp_rat                                                                            
        
    maxCont = max(da_SF.GetMaximum(), mc_full.GetMaximum())
    maxCont_perGeV = max(da_perGeV.GetMaximum(), mc_full_perGeV.GetMaximum())
    da_SF.GetYaxis().SetRangeUser(0.5, 1.50*maxCont)
    da_perGeV.GetYaxis().SetRangeUser(0.5, 1.50*maxCont_perGeV)
    mc_stack.SetMaximum(1.5*maxCont)
    mc_stack.SetMinimum(0.5)
    mc_stack_perGeV.SetMaximum(1.5*maxCont_perGeV)
    SetOwnership(mc_stack, 0);SetOwnership(da_SF, 0);SetOwnership(mc_stack_perGeV, 0);SetOwnership(da_perGeV, 0);SetOwnership(mc_full, 0); SetOwnership(mc_full_perGeV, 0);              



    print helper.bcolors.HEADER + '[result scaled by RSFOF for DATA] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_result = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s'%(newLumiString, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.67, 0.59, 0.90, 0.85)
    plot_result.addStack(mc_stack, "HIST" , 1, 1)
    plot_result.addHisto(mc_full_e, 'E2,SAME'  , '' , 'PL', r.kBlack , 1, -1)
    #plot_result.addHisto(da_SF    , 'E1,SAME'  , 'observed data', 'PL', r.kBlack  , 1,  0)
    plot_result.saveRatio(1, 1, 0, lint, mc_full, mc_full, 0. , int(maxrat+1.0)) 
    #plot_result.saveRatio(1, 1, 0, lint, da_SF, mc_full, 0. , int(maxrat+1.0)) 
    #THIS IS TO UNBLIND!
    #plot_result.saveRatio(1, 1, 0, lint, da_SF, mc_full, 0. , int(maxrat+1.0))                                                                                                       
    makeResultsTable(da_SF, prediction, zz_SF, wz_SF, rare, mc_full)
     
    plot_resultSigs = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s_withSignal'%(newLumiString, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.67, 0.59, 0.90, 0.85)
    plot_resultSigs.addStack(mc_stack, "HIST" , 1, 1)
    plot_resultSigs.addHisto(mc_full_e, 'E2,SAME'  , '' , 'PL', r.kBlack , 1, -1)
    plot_resultSigs.addHisto(sig1, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(signames[0][0], signames[0][1]), "L", r.kGreen-9, 1, 0)
    plot_resultSigs.addHisto(sig2, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(signames[1][0], signames[1][1]), "L", r.kCyan-7, 1, 0)
    plot_resultSigs.addHisto(sig3, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(signames[2][0], signames[2][1]), "L", r.kPink-4, 1, 0)
    plot_resultSigs.saveRatio(1, 1, 0, lint, mc_full, mc_full, 0. , int(maxrat+1.0)) 
    makeResultsTableSig(da_SF, prediction, zz_SF, wz_SF, rare, mc_full, sig1, sig2, sig3,signames)


#    plot_result = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s_perGeV'%(newLumiString, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.67, 0.59, 0.90, 0.85)
#    plot_result.addStack(mc_stack_perGeV, "HIST" , 1, 1)
#    plot_result.addHisto(mc_full_perGeV_e, 'E2,SAME'  , ''  , 'PL', r.kBlack , 1, -1)
#    #plot_result.addHisto(da_perGeV       , 'E1,SAME'  , 'observed data', 'PL', r.kBlack  , 1,  0)
#    plot_result.saveRatio(1, 1, 0, lint, mc_full_perGeV, mc_full_perGeV, 0. , int(maxrat+1.0)) 
#    #plot_result.saveRatio(1, 1, 0, lint, da_perGeV, mc_full_perGeV, 0. , int(maxrat+1.0)) 
    
    tfile = r.TFile('datacards/forDatacards_slepton_2017_%s.root'%scutstring,'recreate')
    tfile.cd()
    tfile.WriteTObject(rare      , 'rare')
    tfile.WriteTObject(wz_SF     , 'wz_SF')
    tfile.WriteTObject(zz_SF     , 'zz_SF')
    tfile.WriteTObject(da_SF     , 'da_SF')
    tfile.WriteTObject(da_OF     , 'da_OF')
    tfile.WriteTObject(result[2] , 'tf_CR_SR')
    tfile.Close()

    if returnplot:
        return plot_result                                                                                                                                                                         
        del mc_stack;da_SF;mc_stack_perGeV;da_perGeV;mc_full; mc_full_perGeV;                                                                                                                                                                                                                                                                                   


if __name__ == '__main__':

    print 'Starting to produce some good ol\' results...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples FilenameWithIngredients', version='%prog 1.0')
    parser.add_option('-t', action='store_true', dest='onlyTTbar', default=False, help='just do OF closure test')
    parser.add_option('-c', action='store_true', dest='onlyClosure', default=False, help='just do the closure test. don\'t bother with data')
    parser.add_option('-b', '--nbs', action='store', type=int, dest='nbs', default=0, help='do this for different numbers of b\'s')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samplesSlepton.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-l', '--loadShapes', action='store_true', dest='loadShapes', default=False, help='reload dy shapes. default is off since this takes a while')
    parser.add_option('-M', '--maxRun', action='store', type=int, dest='maxRun', default=999999, help='max run to use for analysis (run is included)')
    parser.add_option('-m', '--minRun', action='store', type=int, dest='minRun', default=-1    , help='min run to use for analysis (run not included)')
    ## make the options globa.. also the lumi
    global opts, lumi, lumi_str, dy_shapes, nbstring
    global ingMC, ingDA, onZ, treeDA, treeMC, treeDY, treeTT, treeRA
    (opts, args) = parser.parse_args()

    ingredientsFile = opts.ingredientsFile
    print 'running with these options \n'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'
    kf = "noKFactor"
    print 'Going to load DATA and MC trees...'
    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    zzDatasets = ['ZZTo2L2Nu',  'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo2e2nu',]
    wzDatasets = ['WZTo3LNu']
    othersDatasets = ['WWZ', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll', 'TTZToLLNuNu_ext2', 'TTZToQQ', 'TTLLJets_m1to10', 'TTHnobb_pow', 'VHToNonbb']
    fsDatasets = ['TTTT',  'TTTo2L2Nu', 'TBar_tch_powheg', 'T_tch_powheg', 'WWTo2L2Nu', 'WWW', 'WWDouble', 'WpWpJJ', 'TTWToLNu_ext2',  'TTWToQQ', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT']                       
    
    mcDatasets = fsDatasets+dyDatasets + othersDatasets + zzDatasets + wzDatasets
 
                                                                                     
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
    doSleptons = True
    if doSleptons: opts.sampleFile = "samplesSlepton.dat"
    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isOnEOS = 0)
    treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0, isOnEOS = 0)
    treeFS = Sample.Tree(helper.selectSamples(opts.sampleFile, fsDatasets, 'FS'), 'FS'  , 0, isOnEOS = 0)
    treeOTHERS = Sample.Tree(helper.selectSamples(opts.sampleFile, othersDatasets, 'OTHERS'), 'OTHERS'  , 0, isOnEOS = 0)
    treeWZ = Sample.Tree(helper.selectSamples(opts.sampleFile, wzDatasets, 'WZ'), 'WZ'  , 0, isOnEOS = 0)
    treeZZ = Sample.Tree(helper.selectSamples(opts.sampleFile, zzDatasets, 'ZZ'), 'ZZ'  , 0, isOnEOS = 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isOnEOS = 0)
    
    signals = ['SMS_450_50', 'SMS_350_150', 'SMS_250_180']
    signames = [[450, 50], [350, 150], [250, 180]]
    xsec = [2*(0.77+ 0.3),2*(2.33+0.89), 2*(9.210+3.470)]
    treeSIG1 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[0]], 'SI'), 'SI', 0, isScan = xsec[0])
    treeSIG2 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[1]], 'SI'), 'SI', 0, isScan = xsec[1])
    treeSIG3 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[2]], 'SI'), 'SI', 0, isScan = xsec[2])
    sigs = [treeSIG1, treeSIG2, treeSIG3]                                                               

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
    #rsfof, rsfof_e, rsfof_mc, rsfof_mc_e, rmue_a_da, rmue_a_mc, rmue_b_da, rmue_b_mc =  makeFactorsTable()
    ## result plots in different variables:
    for v in ['met']:#'nll_noMET', 'nll_noMLB', 'nll_noZPT', 'nll_noLDP']: #'iso1', 'iso2', 'mll', 'nll', 'nb', 'nj', 'zpt', 'mlb', 'met', 'ldp', 'pt1', 'pt2']:
        makeResultData('slepton2017', v,signames, maxrun,lint,specialcut='' , scutstring = '',    _options='returnplot,splitFlavor')
        #makeClosureTestPlots('slepton2017','met0jet','nJet25_Edge  == 0 && mt2_Edge > 90' , 'ht0', True)
        #makeClosureTestPlots('slepton2017',v,'', 'inclusive', True)
        #resultPlotLoNll = makeResultData('Edge_Moriond2017', v,maxrun,lint,specialcut='nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) < 21.' , scutstring = 'nllBelow21',    _options='returnplot,splitFlavor')
        #resultPlotHiNll = makeResultData('Edge_Moriond2017', v,maxrun,lint,specialcut='nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) >= 21.' , scutstring = 'nllAbove21',    _options='returnplot,splitFlavor')

   # makeRSOFTable('Edge_Moriond2017')
    
