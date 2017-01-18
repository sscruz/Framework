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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats, TStyle, SetOwnership
import math, sys, optparse, copy, re, array, subprocess

import include.nll      
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

def makeDYMllShape(var, specialcut = '', scutstring = ''):
    isBlinded  = True
    if var == 'mll':
        treevar = 'lepsMll_Edge'
        xlabel = 'm_{ll} [GeV]'                    
        nbins = [20, 60, 86, 96, 150, 200, 300, 400]
    if isBlinded:
        lint = 18.1  ; maxrun = 999999; lint_str = '18.1invfb'
        pred  = 46.5; pred_e =14.92 
    else:
        lint = 36.4  ; maxrun = 999999; lint_str = '36.4invfb'
        pred = 96.5; pred_e =38.24
    
    rinout1 = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m20_60__")[0]
    rinout1_stat = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m20_60__")[1]
    rinout1_syst = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m20_60__")[2]
    rinout1_e = math.sqrt(rinout1_stat**2 + rinout1_syst**2)                          
    rinout2 = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m60_86__")[0]
    rinout2_stat = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m60_86__")[1]
    rinout2_syst = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m60_86__")[2]
    rinout2_e = math.sqrt(rinout2_stat**2 + rinout2_syst**2)               
    rinout3 = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m96_150_")[0]
    rinout3_stat = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m96_150_")[1]
    rinout3_syst = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m96_150_")[2]
    rinout3_e = math.sqrt(rinout3_stat**2 + rinout3_syst**2)                         
    rinout4 = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m150_200")[0]
    rinout4_stat = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m150_200")[1]
    rinout4_syst = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m150_200")[2]
    rinout4_e = math.sqrt(rinout4_stat**2 + rinout4_syst**2)                         
    rinout5 = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m200_300")[0]
    rinout5_stat = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m200_300")[1]
    rinout5_syst = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m200_300")[2]
    rinout5_e = math.sqrt(rinout5_stat**2 + rinout5_syst**2)                         
    rinout6 = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m300_400")[0]
    rinout6_stat = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m300_400")[1]
    rinout6_syst = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m300_400")[2]
    rinout6_e = math.sqrt(rinout6_stat**2 + rinout6_syst**2)                             
    rinout7 = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m400___")[0]
    rinout7_stat = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m400___")[1]
    rinout7_syst = helper.readFromFileRinout(ingredientsFile, "DATA", "dy_m400___")[2]
    rinout7_e = math.sqrt(rinout7_stat**2 + rinout7_syst**2)                            


    dy_shape=treeDY.getTH1F(lint,"mll_shape",'lepsMll_Edge', nbins, 1, 1,  cuts.AddList([cuts.goodLepton, cuts.SignalRegion, cuts.Zmass, cuts.SF]), '',xlabel)
    dy_loNll=treeDY.getTH1F(lint,"mll_loNll",'lepsMll_Edge', nbins, 1, 1,  cuts.AddList(["nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge)  <= 21.", cuts.goodLepton, cuts.SignalRegion, cuts.Zmass, cuts.SF]), '',xlabel)
    dy_hiNll=treeDY.getTH1F(lint,"mll_hiNll",'lepsMll_Edge', nbins, 1, 1,  cuts.AddList(["nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) >  21.", cuts.goodLepton, cuts.SignalRegion, cuts.Zmass, cuts.SF]), '',xlabel)
    if scutstring == "nllAbove21":
        pred = pred * (dy_hiNll.Integral()/dy_shape.Integral()) 
    
    if scutstring == "nllBelow21":
        pred = pred * (dy_loNll.Integral()/dy_shape.Integral())
    dy_shape.SetBinContent(1, pred*rinout1); dy_shape.SetBinError(1, pred*rinout1*math.sqrt((pred_e/pred)**2 + (rinout1_e/rinout1)**2))
    dy_shape.SetBinContent(2, pred*rinout2); dy_shape.SetBinError(2, pred*rinout2*math.sqrt((pred_e/pred)**2 + (rinout2_e/rinout1)**2))
    dy_shape.SetBinContent(3, 0); dy_shape.SetBinError(3, 0)
    dy_shape.SetBinContent(4, pred*rinout3); dy_shape.SetBinError(4, pred*rinout3*math.sqrt((pred_e/pred)**2 + (rinout3_e/rinout1)**2))
    dy_shape.SetBinContent(5, pred*rinout4); dy_shape.SetBinError(5, pred*rinout4*math.sqrt((pred_e/pred)**2 + (rinout4_e/rinout1)**2))
    dy_shape.SetBinContent(6, pred*rinout5); dy_shape.SetBinError(6, pred*rinout5*math.sqrt((pred_e/pred)**2 + (rinout5_e/rinout1)**2))
    dy_shape.SetBinContent(7, pred*rinout6); dy_shape.SetBinError(7, pred*rinout6*math.sqrt((pred_e/pred)**2 + (rinout6_e/rinout1)**2))
    dy_shape.SetBinContent(8, pred*rinout7); dy_shape.SetBinError(8, pred*rinout7*math.sqrt((pred_e/pred)**2 + (rinout7_e/rinout1)**2))
    plot_dy_shape = Canvas.Canvas('results/18.1invfb/plot_Templates', 'png,pdf', 0.60, 0.65, 0.80, 0.85)
    plot_dy_shape.addHisto(dy_shape  , 'HIST'       , 'Templates', 'PL', r.kBlack  , 1,  0)
    plot_dy_shape.save(1, 0, 0, lint)                            
    
    
    return  dy_shape                                                                                                                                  


def makeResultsTable(da, fs, dy, ra, mc, nll = ''):
    
    line0 = '\\begin{tabular}{c c c c c c }  '                                                               
    line2 = ' M_{ll} & Flavour symmetric & Drell-Yan & Rares & Total & Data \\\\ \hline  '                                                               
    print "making table ", nll
    if nll == "nllAbove21":
        line1 =  '&& ttbar-like  &&& \\\\ \hline'
    if nll == "nllBelow21":
        line1 =  '&& non ttbar-like  &&& \\\\ \hline'
    if nll == "":
        line1 =  '&& Total (inclusive in nll)  &&& \\\\ \hline'
    line3 = '20-60  &  %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f   & %.f  \\\  ' %(fs.GetBinContent(1),fs.GetBinError(1),dy.GetBinContent(1),dy.GetBinError(1),ra.GetBinContent(1),ra.GetBinError(1), mc.GetBinContent(1),mc.GetBinError(1) , da.GetBinContent(1))
    line4 = '60-86  &  %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f   & %.f  \\\  ' %(fs.GetBinContent(2),fs.GetBinError(2),dy.GetBinContent(2),dy.GetBinError(2),ra.GetBinContent(2),ra.GetBinError(2), mc.GetBinContent(2),mc.GetBinError(2) , da.GetBinContent(2))
    line5 = '96-150  &  %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f   & %.f  \\\  ' %(fs.GetBinContent(4),fs.GetBinError(4),dy.GetBinContent(4),dy.GetBinError(4),ra.GetBinContent(4),ra.GetBinError(4), mc.GetBinContent(4),mc.GetBinError(4) , da.GetBinContent(4))
    line6 = '150-200  &  %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f   & %.f  \\\  ' %(fs.GetBinContent(5),fs.GetBinError(5),dy.GetBinContent(5),dy.GetBinError(5),ra.GetBinContent(5),ra.GetBinError(5), mc.GetBinContent(5),mc.GetBinError(5) , da.GetBinContent(5))
    line7 = '200-300  &  %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f   & %.f  \\\  ' %(fs.GetBinContent(6),fs.GetBinError(6),dy.GetBinContent(6),dy.GetBinError(6),ra.GetBinContent(6),ra.GetBinError(6), mc.GetBinContent(6),mc.GetBinError(6) , da.GetBinContent(6))
    line8 = '300-400  &  %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f   & %.f  \\\  ' %(fs.GetBinContent(7),fs.GetBinError(7),dy.GetBinContent(7),dy.GetBinError(7),ra.GetBinContent(7),ra.GetBinError(7), mc.GetBinContent(7),mc.GetBinError(7) , da.GetBinContent(7))
    line9 = ' + 400  &  %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f   & %.f  \\\  ' %(fs.GetBinContent(8),fs.GetBinError(8),dy.GetBinContent(8),dy.GetBinError(8),ra.GetBinContent(8),ra.GetBinError(8), mc.GetBinContent(8),mc.GetBinError(8) , da.GetBinContent(8))
    line10 = '\\hline'
    line11 = '\end{tabular}'
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
    print line11                                                                                                                                                      


def makeClosureTests(var, specialcut = '', scutstring = '', doCumulative = False):

    if var == 'mll':
        treevar = 'lepsMll_Edge'
        nbins, xmin, xmax = 28, 20, 300
        xlabel = 'm_{ll} (GeV)'
    elif var == 'nll':
        treevar = 'nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge)'
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
    mc_OF = treeFS.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.OF]), '', xlabel)
    mc_SF = treeFS.getTH1F(lint, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel)
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.trigger,  cuts.Zveto, cuts.OF]), '', xlabel)
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.trigger,  cuts.Zveto, cuts.SF]), '', xlabel)
    #mc_ee = treeMC.getTH1F(lint, var+"mc_ee"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.Zveto, cuts.ee]), '', xlabel)
    #mc_mm = treeMC.getTH1F(lint, var+"mc_mm"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.Zveto, cuts.mm]), '', xlabel)
    dy_SF = treeDY.getTH1F(lint, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel)
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

def makeResultData(var, maxrun = 999999, lint = 36.4, specialcut = '', scutstring = '', _options = ''):
    returnplot, addRares, splitFlavor, makeTable, printIntegral = False, True, False, False, False
    if 'returnplot'    in _options: print 'found option %s'%'returnplot'    ;returnplot    = True;addRares      = True
    if 'splitFlavor'   in _options: print 'found option %s'%'splitFlavor'   ;splitFlavor   = True
    if 'makeTable'     in _options: print 'found option %s'%'makeTable'     ;makeTable     = True
    if 'printIntegral' in _options: print 'found option %s'%'printIntegral' ;printIntegral = True

    # if   var == 'mll'      : treevar = 'lepsMll_Edge'        ; nbins, xmin, xmax = 23, 20 , 250 ; xlabel = 'm_{ll} (GeV)'
    if   var == 'mll'      : treevar = 'lepsMll_Edge'        ; nbins, xmin, xmax = 30, 20 , 420 ; xlabel = 'm_{ll} [GeV]'
    elif var == 'nll'      : treevar = 'nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge)'     ; nbins, xmin, xmax = 15, 12 , 27  ; xlabel = 'NLL'
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
        
    nbins = [20, 60, 86, 96, 150, 200, 300, 400]
    newLumiString = '18.1invfb'
    if not specialcut:
        print "blind stuff!!!!"
        specialcut = specialcut + '((run_Edge <=276811) ||  (278820<=run_Edge && run_Edge<=279931))'
    else:
        specialcut = specialcut + '&&((run_Edge <=276811) ||  (278820<=run_Edge && run_Edge<=279931))'
    ## ## mll distributions
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.trigger, cuts.SF, cuts.Zveto]), '', xlabel)
    da_mm = treeDA.getTH1F(lint, var+"da_mm"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.trigger,cuts.mm, cuts.Zveto]), '', xlabel)
    da_ee = treeDA.getTH1F(lint, var+"da_ee"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.trigger,cuts.ee, cuts.Zveto]), '', xlabel)
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.trigger,cuts.OF, cuts.Zveto]), '', xlabel)
    ra_OF = treeRA.getTH1F(lint, var+"ra_OF"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.OF, cuts.Zveto]), '', xlabel)
    ra_SF = treeRA.getTH1F(lint, var+"ra_SF"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.SF, cuts.Zveto]), '', xlabel)
    ra_SF.Add(ra_OF, -1.)                                                                                                                                                                      
    da_OF_err = copy.deepcopy(da_OF)
    da_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    da_OF_err.SetFillStyle(3017); da_OF_err.SetMarkerSize(0.)

    da_OF_rsfofScaled = copy.deepcopy(da_OF)
    da_OF_rsfofScaled = scaleByRSFOF(da_OF_rsfofScaled, rsfof, rsfof_e)

    da_OF_rsfofScaled_err = copy.deepcopy(da_OF_rsfofScaled)
    da_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    da_OF_rsfofScaled_err.SetFillStyle(3017); da_OF_rsfofScaled_err.SetMarkerSize(0.)
    
    dy_shape = makeDYMllShape('mll',specialcut,scutstring )
    
    mc_stack = r.THStack() 
   
    da_SF.SetTitle("data SF")                                                                           
    ra_SF.SetFillColorAlpha(r.kGreen-5, 0.8);ra_SF.SetTitle("rares");ra_SF.SetLineColor(r.kBlack)
    dy_shape.SetFillColorAlpha(r.kYellow-9, 0.8);dy_shape.SetTitle("E_{T}^{miss} templates")
    da_OF_rsfofScaled.SetFillColorAlpha(r.kRed-9, 0.8);da_OF_rsfofScaled.SetTitle("FS")                                  
    mc_full = copy.deepcopy(ra_SF)
    mc_stack.Add(ra_SF); 
    mc_stack.Add(dy_shape); mc_full.Add(dy_shape, 1.) 
    mc_stack.Add(da_OF_rsfofScaled); mc_full.Add(da_OF_rsfofScaled, 1.) 
    mc_full_e = copy.deepcopy(mc_full)
    mc_full_e.SetFillColorAlpha(r.kBlue+1, 0.8);mc_full_e.SetFillStyle(3017); mc_full_e.SetMarkerSize(0.)
    maxrat = 0.5
    for ib in range(1,da_SF.GetNbinsX()+1):
        tmp_rat = da_SF.GetBinContent(ib)/( da_OF_rsfofScaled.GetBinContent(ib) if da_OF_rsfofScaled.GetBinContent(ib) > 0 else 1. )
        if tmp_rat > maxrat:
            maxrat = tmp_rat
    
    maxCont = max(da_SF.GetMaximum(), mc_full.GetMaximum())
    da_OF_rsfofScaled_err.GetYaxis().SetRangeUser(0.1, 1.50*maxCont)
    mc_stack.SetMaximum(1.5*maxCont)
    SetOwnership(mc_stack, 0 );SetOwnership(da_SF, 0 )                                                                                                                       

    print helper.bcolors.HEADER + '[result scaled by RSFOF for DATA] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    sstring = '' if not addRares else ''
    plot_result = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s'%(newLumiString+sstring, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.60, 0.65, 0.80, 0.85)
    plot_result.addStack(mc_stack, "HIST" , 1, 1)
    plot_result.addHisto(mc_full_e, 'e2,same'  , ''         , 'PL', r.kBlue+1 , 1, -1)
    #plot_result.addHisto(da_OF_rsfofScaled    , 'hist,SAME', 'predicted', 'L' , r.kBlue+1 , 1,  1)
    plot_result.addHisto(da_SF                , 'E1,SAME'    , 'observed data', 'PL', r.kBlack  , 1,  0)
    if maxrun < 999999: plot_result.addLatex (0.2, 0.85, 'max. run %s'%("run <=276811) ||  (278820 < run && run < 279931)"))
    if printIntegral  : plot_result.addLatex (0.2, 0.65, 'SF: {SF:.1f} (#mu#mu:{mm:.1f},ee:{ee:.1f})'.format(SF=da_SF.Integral(),mm=da_mm.Integral(),ee=da_ee.Integral() ) )
    if printIntegral  : plot_result.addLatex (0.2, 0.60, 'OF: {OF:.2f}'.format(OF=da_OF_rsfofScaled.Integral() ) )
    plot_result.saveRatio(1, 1, 0, lint, da_SF, mc_full) 
    print "making table ", scutstring
    makeResultsTable(da_SF, da_OF_rsfofScaled, dy_shape, ra_SF, mc_full, scutstring)
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
    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    #dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext', 'DYJetsToLL_M50_HT400to600_ext', 'DYJetsToLL_M50_HT600toInf_ext']
    raDatasets = ['GGHZZ4L', 'ZZTo2L2Q', 'ZZTo2L2Nu', 'ZZTo4L', 'TTTT', 'tZq_ll', 'TWZ', 'WZTo3LNu',  'WWZ','WZZ', 'ZZZ',  'TTZToLLNuNu' ,'TTWToLNu', 'TTHnobb_pow', 'VHToNonbb']
    fsDatasets = ['TTJets_DiLepton_ext', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'WWTo2L2Nu', 'WWW', 'TTWToLNu','TTWToQQ', 'T_tWch', 'TBar_tWch' ,'TToLeptons_sch', 'TToLeptons_tch_powheg', 'TBarToLeptons_tch_powheg',  'WJetsToLNu_LO']
    mcDatasets = fsDatasets + dyDatasets + raDatasets

    daDatasets = ['DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleMuon_Run2016F_23Sep2016_v1_runs_271036_284044_part5',
                  'MuonEG_Run2016F_23Sep2016_v1_runs_271036_284044',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part3',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part4',
                  'DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376_part5',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part10',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part11',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part3',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part4',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part5',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part7',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part1',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part4',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part5',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part6',
                  'DoubleEG_Run2016H-PromptReco-v3_runs_284036_284044',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part1',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part10',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part3',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part4',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part5',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part7',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part8',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part9',
                  'DoubleMuon_Run2016H-PromptReco-v3_runs_284036_284044',
                  'MuonEG_Run2016H-PromptReco-v2_runs_281613_284035',
                  'MuonEG_Run2016H-PromptReco-v3_runs_284036_284044',
                  'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part1',
                  'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part2',
                  'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part3',
                  'MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part4',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part8',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part9',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part2',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part6',
                  'DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleMuon_Run2016C_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016C_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part5',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part6',
                  'DoubleMuon_Run2016D_23Sep2016_v1_runs_271036_284044_part7',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part5',
                  'DoubleMuon_Run2016E_23Sep2016_v1_runs_271036_284044_part6',
                  'MuonEG_Run2016C_23Sep2016_v1_runs_271036_284044',
                  'MuonEG_Run2016D_23Sep2016_v1_runs_271036_284044',
                  'MuonEG_Run2016E_23Sep2016_v1_runs_271036_284044',
                  'DoubleMuon_Run2016H-PromptReco-v2_runs_281613_284035_part2',
                  'DoubleEG_Run2016H-PromptReco-v2_runs_281613_284035_part3',
                  'DoubleMuon_Run2016B_23Sep2016_v3_runs_273150_275376_part6',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part5',
                  'DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044_part6',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part2',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part3',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part4',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part5',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part6',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part7',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part8',
                  'DoubleMuon_Run2016G_23Sep2016_v1_runs_271036_284044_part9',
                  'MuonEG_Run2016G_23Sep2016_v1_runs_271036_284044_part1',
                  'MuonEG_Run2016G_23Sep2016_v1_runs_271036_284044_part2']                



    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0)
    treeFS = Sample.Tree(helper.selectSamples(opts.sampleFile, fsDatasets, 'FS'), 'FS'  , 0)
    treeRA = Sample.Tree(helper.selectSamples(opts.sampleFile, raDatasets, 'RA'), 'RA'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'

    isBlinded = True

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()
    lint = 18.1  ; maxrun = 999999; lint_str = '18.1invfb'
    #lint = 36.4  ; maxrun = 999999; lint_str = '36.4invfb'
    print 'Running with an integrated luminosity of %.2f fb-1' %(lint)

    ## ============================================================
    ## ========== set RSFOF globally ==============================
    ## ============================================================
    global rsfof, rsfof_e, rsfof_mc, rsfof_mc_e
    rsfof, rsfof_e, rsfof_mc, rsfof_mc_e =  makeFactorsTable()
    print rsfof, rsfof_e, rsfof_mc, rsfof_mc_e
    ## result plots in different variables:
        
    for v in ['mll']:#'nll_noMET', 'nll_noMLB', 'nll_noZPT', 'nll_noLDP']: #'iso1', 'iso2', 'mll', 'nll', 'nb', 'nj', 'zpt', 'mlb', 'met', 'ldp', 'pt1', 'pt2']:
        makeResultData(v,maxrun,lint,specialcut='' , scutstring = '',    _options='returnplot,splitFlavor')
        resultPlotHiNll = makeResultData(v,maxrun,lint,specialcut='nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21.' , scutstring = 'nllAbove21',    _options='returnplot,splitFlavor')
        resultPlotLoNll = makeResultData(v,maxrun,lint,specialcut='nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) <= 21.' , scutstring = 'nllBelow21',    _options='returnplot,splitFlavor')


    makeClosureTests('nll','', 'inclusive', True)
    makeClosureTests('nll','lepsMll_Edge <= 60. && lepsMll_Edge > 20', 'mll20-60', True)
    makeClosureTests('nll','lepsMll_Edge <= 86. && lepsMll_Edge > 60', 'mll60-86', True)
    makeClosureTests('nll','lepsMll_Edge <= 150. && lepsMll_Edge > 96', 'mll96-150', True)
    makeClosureTests('nll','lepsMll_Edge <= 200. && lepsMll_Edge > 150', 'mll150-200', True)
    makeClosureTests('nll','lepsMll_Edge <= 300. && lepsMll_Edge > 200', 'mll200-300', True)
    makeClosureTests('nll','lepsMll_Edge <= 400. && lepsMll_Edge > 300', 'mll300-400', True)
    makeClosureTests('nll','lepsMll_Edge > 400', 'mll400', True)
    makeClosureTests('mll','', 'inclusive', True)
    makeClosureTests('mll','nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21.' , 'nllAbove21', True)
    makeClosureTests('mll','nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) <= 21.' , 'nllBelow21', True)
    
    #makeResultData('mll', maxrun, lint, specialcut = '' , scutstring = '', returnplot = True)

    ## resultPlot = makeResultData('mll', maxrun, lint, specialcut = '' , scutstring = '', returnplot = True)
    ## ## makeResultData('nll', maxrun = 275125, lint = 4.0)
    ## resultPlotLoNLL = makeResultData('mll',    maxrun,        lint, 'nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) < 21.', 'nllBelow21', returnplot = True)
    ## resultPlotHiNLL = makeResultData('mll',    maxrun,        lint, 'nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21.', 'nllAbove21', returnplot = True)
    ## print addsf
    ## make for region with fixed trigger:
    #### random shit resultPlot = makeResultData('mll', maxrun, lint, specialcut = 'run_Edge >= 274094' , scutstring = 'maxrun274094', returnplot = True)
    #### random shit ##makeResultData('nll', maxrun = 274240, lint = 0.864)
    #### random shit resultPlotLoNLL = makeResultData('mll',    maxrun,        lint, 'nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) < 21. && run_Edge >= 274094', 'nllBelow21_maxrun274094', returnplot = True)
    #### random shit resultPlotHiNLL = makeResultData('mll',    maxrun,        lint, 'nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21. && run_Edge >= 274094', 'nllAbove21_maxrun274094', returnplot = True)
    #####print asdfsadf
    ## makeClosureTests('mll','run_Edge <= 999999', 'withWZ')
    ##makeClosureTests('mll','nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21. && run_Edge <= 999999', 'nllAbove21')
    ##makeClosureTests('mll','nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21. && run_Edge <= 999999 && nBJetMedium25_Edge > 1', 'nllAbove21_nb2')
    #makeClosureTests('nll','lepsMll_Edge > 101. && run_Edge <= 274240', 'highMass_0p8fb-1')

