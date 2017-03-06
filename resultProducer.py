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


def makeFactorsTable(): ## for this to make sense the region should be properly binned!!
    rmue_da = helper.readFromFileRmue("ingredients.dat", "DATA") 
    rmue_mc = helper.readFromFileRmue("ingredients.dat", "MC") 
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA","DATA")
    rmue_a_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffA","MC")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB","DATA")
    rmue_b_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffB","MC")
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
    return (rsfof_da[0]+rsfof_fac_da)/2, getFinalError(rsfof_fac_da_e,  getFinalError(rsfof_da[1], rsfof_da[2])), (rsfof_mc[0]+rsfof_fac_mc)/2, getFinalError(rsfof_fac_mc_e,  getFinalError(rsfof_mc[1], rsfof_mc[2])) , rmue_a_da, rmue_a_mc, rmue_b_da, rmue_b_mc                                                                                   



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


def makeDYMllShape(var, specialcut = '', scutstring = ''):
    isBlinded  = True
    if var == 'mll':
        treevar = 'lepsMll_Edge'
        xlabel = 'm_{ll} [GeV]'                    
        nbins = [20, 60, 86, 96, 150, 200, 300, 400]
    lint = 36.8  ; maxrun = 999999; lint_str = '36.8invfb'
    pred = 21.6; pred_e =13.1
    
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


    dy_shape=treeDY.getTH1F(lint,"mll_shape",'lepsMll_Edge', nbins, 1, 1,  cuts.AddList([cuts.goodLepton, cuts.SignalRegionNoDPhi, cuts.Zmass, cuts.SF]), '',xlabel)
    dy_loNll=treeDY.getTH1F(lint,"mll_loNll",'lepsMll_Edge', nbins, 1, 1,  cuts.AddList(["nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) < 21.", cuts.goodLepton, cuts.SignalRegionNoDPhi, cuts.Zmass, cuts.SF]), '',xlabel)
    dy_hiNll=treeDY.getTH1F(lint,"mll_hiNll",'lepsMll_Edge', nbins, 1, 1,  cuts.AddList(["nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) >=  21.", cuts.goodLepton, cuts.SignalRegionNoDPhi, cuts.Zmass, cuts.SF]), '',xlabel)
    if scutstring == "nllAbove21":
        print "nll above 21 ", (dy_hiNll.Integral()/dy_shape.Integral()) 
        pred = pred * (dy_hiNll.Integral()/dy_shape.Integral()) 
    if scutstring == "nllBelow21":
        print "nll below 21 ", (dy_loNll.Integral()/dy_shape.Integral())
        pred = pred * (dy_loNll.Integral()/dy_shape.Integral())
    dy_shape.SetBinContent(1, pred*rinout1); dy_shape.SetBinError(1, pred*rinout1*math.sqrt((pred_e/pred)**2 + (rinout1_e/rinout1)**2))
    dy_shape.SetBinContent(2, pred*rinout2); dy_shape.SetBinError(2, pred*rinout2*math.sqrt((pred_e/pred)**2 + (rinout2_e/rinout2)**2))
    dy_shape.SetBinContent(3, 0); dy_shape.SetBinError(3, 0)
    dy_shape.SetBinContent(4, pred*rinout3); dy_shape.SetBinError(4, pred*rinout3*math.sqrt((pred_e/pred)**2 + (rinout3_e/rinout3)**2))
    dy_shape.SetBinContent(5, pred*rinout4); dy_shape.SetBinError(5, pred*rinout4*math.sqrt((pred_e/pred)**2 + (rinout4_e/rinout4)**2))
    dy_shape.SetBinContent(6, pred*rinout5); dy_shape.SetBinError(6, pred*rinout5*math.sqrt((pred_e/pred)**2 + (rinout5_e/rinout5)**2))
    dy_shape.SetBinContent(7, pred*rinout6); dy_shape.SetBinError(7, pred*rinout6*math.sqrt((pred_e/pred)**2 + (rinout6_e/rinout6)**2))
    dy_shape.SetBinContent(8, pred*rinout7); dy_shape.SetBinError(8, pred*rinout7*math.sqrt((pred_e/pred)**2 + (rinout7_e/rinout7)**2))
    plot_dy_shape = Canvas.Canvas('results/36.8invfb/plot_Templates', 'png,pdf', 0.60, 0.65, 0.80, 0.85)
    plot_dy_shape.addHisto(dy_shape  , 'HIST'       , 'Templates', 'PL', r.kBlack  , 1,  0)
    plot_dy_shape.save(1, 0, 0, lint)                            
    return  dy_shape                                                                                                                                  


def makeResultsTable(da, fs, dy,rare, mc, nll = ''):
    
    line0 = '\\begin{tabular}{c c c c c c}  '                                                               
    line2 = ' M_{ll} & Flavour symmetric & Drell-Yan & Rares & Total & Data \\\\ \hline  '                                                               
    print "making table ", nll
    if nll == "nllAbove21":
        line1 =  '&& ttbar-like  &&& \\\\ \hline'
    if nll == "nllBelow21":
        line1 =  '&& non ttbar-like  &&& \\\\ \hline'
    if nll == "":
        line1 =  '&& Total (inclusive in nll)  &&& \\\\ \hline'
    mllLabels = ['20-60 ', '60-86 ','onz', '96-150 ', '150-200','200-300','300-400', '+400']
    lines = []
    print line0                                                                                                                                                      
    print line1                                                                                                                                                      
    print line2                                                                                                                                                      
    for bin,mllLabel in enumerate(mllLabels):
        if mllLabel == 'onz': continue
        print '%s &   %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f  & %.2f $\\pm$ %.2f   & %.2f $\\pm$ %.2f   & %.f  \\\  ' %(mllLabel,
                                                                                                                      fs.GetBinContent(bin+1),
                                                                                                                      fs.GetBinError(bin+1),
                                                                                                                      dy.GetBinContent(bin+1),
                                                                                                                      dy.GetBinError(bin+1),
                                                                                                                      rare.GetBinContent(bin+1),
                                                                                                                      rare.GetBinError(bin+1),
                                                                                                                      mc.GetBinContent(bin+1),
                                                                                                                      mc.GetBinError(bin+1) ,
                                                                                                                      da.GetBinContent(bin+1))

    print '\\hline'
    print '\end{tabular}'



def makeClosureTests(var, specialcut = '', scutstring = '', doCumulative = False, nbins=0, xmin=0, xmax=0, save=True):

    if var == 'mll':
        treevar = 'lepsMll_Edge'
        nbins, xmin, xmax = 29, 16, 306
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

                                                                                                                                
    mc_OF = treeFS.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.OF]), '', xlabel)
    mc_SF = treeFS.getTH1F(lint, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel)
    dy_SF = treeDY.getTH1F(lint, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel)
    mc_OF_factor = treeFS.getTH1F(lint, var+"mc_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))
    mc_OF_factorUp = treeFS.getTH1F(lint, var+"mc_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))
    mc_OF_factorDn = treeFS.getTH1F(lint, var+"mc_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))

    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_direct = copy.deepcopy(mc_OF) 
    mc_OF_direct = scaleByRSFOF(mc_OF_direct, rsfof_mc[0], rsfof_mc[1])
    mc_OF_factor = getRMueError(mc_OF_factor, mc_OF_factorUp, mc_OF_factorDn)
    mc_OF_factor = scaleByRSFOF(mc_OF_factor, rt_mc[0], getFinalError(rt_mc[1],rt_mc[2]))
    result = weightedAverage( mc_OF_factor, mc_OF_direct, mc_OF)
    return result


def makeClosureTestPlots(analysis, var, specialcut = '', scutstring = '', doCumulative = False, nbins=0, xmin=0, xmax=0, save=True):

    if var == 'mll':
        treevar = 'lepsMll_Edge'
        xlabel = 'm_{ll} (GeV)'
        nbins, xmin, xmax = 30, 20, 300
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
    else: 
        treevar = var
        xlabel = ''
    scan = Scans.Scan(analysis) 
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rsfof_mc = helper.readFromFileRsfofD("ingredients.dat", "MC")   
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rt_mc = helper.readFromFileRT("ingredients.dat", "MC")  
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "DATA")
    rmue_a_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "MC")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "DATA")
    rmue_b_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "MC")

    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.OF, cuts.Zveto]), '', xlabel)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.OF, cuts.Zveto]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
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
        da_prediction.SetBinError(i,0.)
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
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel)
    mc_OF = treeFS.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.OF]), '', xlabel)
    mc_SF = treeFS.getTH1F(lint, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel)
    dy_SF = treeDY.getTH1F(lint, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel)
    mc_OF_factor = treeFS.getTH1F(lint, var+"mc_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))
    mc_OF_factorUp = treeFS.getTH1F(lint, var+"mc_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))
    mc_OF_factorDn = treeFS.getTH1F(lint, var+"mc_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]))

    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_err.SetFillStyle(3004); mc_OF_err.SetMarkerSize(0.)
    dy_SF.SetFillColorAlpha(r.kGreen+2,0.5)

    mc_OF_direct = copy.deepcopy(mc_OF) 
    mc_OF_direct = scaleByRSFOF(mc_OF_direct, rsfof_mc[0], rsfof_mc[1])
    mc_OF_factor = getRMueError(mc_OF_factor, mc_OF_factorUp, mc_OF_factorDn)
    mc_OF_factor = scaleByRSFOF(mc_OF_factor, rt_mc[0], getFinalError(rt_mc[1],rt_mc[2]))
    result = weightedAverage( mc_OF_factor, mc_OF_direct, mc_OF)
    mc_OF_rsfofScaled = result[3] 
    mc_OF_rsfofScaled_err = copy.deepcopy( mc_OF_rsfofScaled ) 
    mc_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_rsfofScaled_err.SetFillStyle(3004); mc_OF_rsfofScaled_err.SetMarkerSize(0.)          

    maxCont = max(da_prediction.GetMaximum(), mc_OF_rsfofScaled.GetMaximum())
    mc_SF.GetYaxis().SetRangeUser(0.1, 1.30*maxCont)
    mc_OF_rsfofScaled.GetYaxis().SetRangeUser(0.1, 1.30*maxCont)

    maxrat = 0.5
    for ib in range(1,da_SF.GetNbinsX()+1):
        tmp_rat = da_SF.GetBinContent(ib)/( mc_OF_rsfofScaled.GetBinContent(ib) if mc_OF_rsfofScaled.GetBinContent(ib) > 0 else 1. )
        if tmp_rat > maxrat:
            maxrat = tmp_rat                                                                                          
    mc_SF.GetYaxis().SetRangeUser(0., 1.3*mc_SF.GetMaximum())
    print helper.bcolors.HEADER + '[MC only closure test not scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
   # plot_closure_noRSFOF = Canvas.Canvas('closure/%s/plot_closure_mll_mcPredmcObs%s_noRSFOF'%(lint_str,  '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
   # plot_closure_noRSFOF.addHisto(mc_SF    , 'PE'       , 'MC - SF', 'PL', r.kRed+1  , 1,  0)
   # plot_closure_noRSFOF.addHisto(mc_OF_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
   # plot_closure_noRSFOF.addHisto(mc_OF    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
   # #plot_closure_noRSFOF.addHisto(dy_SF    , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
   # plot_closure_noRSFOF.addLatex (0.61, 0.82, 'no R_{SFOF} scaling')
   # plot_closure_noRSFOF.saveRatio(1, 0, 1, lint, mc_SF, mc_OF, 0.2, 1.8)

    print helper.bcolors.HEADER + '[MC only closure test scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure = Canvas.Canvas('closure/%s/plot_closure_mll_mcPredmcObs%s'%(lint_str,  '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_closure.addHisto(mc_SF                , 'PE'       , 'MC - SF', 'PL', r.kRed+1  , 1,  0)
    plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure.addHisto(da_prediction, 'hist,SAME'  , 'data - OF'       , 'L', r.kBlack , 1, 1)
    plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    #plot_closure.addHisto(dy_SF                , 'hist,SAME', 'DY - SF', 'FL', r.kGreen+2, 1,  2)
    plot_closure.addLatex (0.61, 0.82, 'R_{SFOF} scaled')
    plot_closure.saveRatio(1, 0, 0, lint, da_prediction, mc_OF_rsfofScaled,  0. , int(maxrat+1.0))                                                                                              
                                                                                                                                                                             
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
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","A","DATA")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","B","DATA")

    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    print '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0])
    print 'trivial check'
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue

    da_OF_direct = copy.deepcopy(da_OF)
    da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], rsfof_da[1])
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    result = weightedAverage( da_OF_factor, da_OF_direct, da_OF)                                   
    ### 

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


def makeResultData(analysis, var, maxrun = 999999, lint = 36.8, specialcut = '', scutstring = '', _options = ''):
    returnplot, addRares, splitFlavor, makeTable, printIntegral = True, True, False, False, False
    if   var == 'mll'      : treevar = 'lepsMll_Edge'        ; nbins = [20, 60, 86, 96, 150, 200, 300, 400]; xmin =1 ; xmax = 1; xlabel = 'm_{ll} [GeV]'
    if not specialcut:
        specialcut = specialcut 
    else:
        specialcut = specialcut 
    scan = Scans.Scan(analysis)
    mc_stack = r.THStack()
    mc_stack_perGeV = r.THStack() 
    newLumiString = '36.8invfb'
    ##get the ingredients
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","A","DATA")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","B","DATA")
    #get the fs prediction
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.OF, cuts.Zveto]), '', xlabel)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.OF, cuts.Zveto]), '', xlabel,extraWeight='(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,extraWeight='(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]))
    print '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0])
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue
    da_OF_direct = copy.deepcopy(da_OF)
    da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], rsfof_da[1])
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    result = weightedAverage( da_OF_factor, da_OF_direct, da_OF)                                   
    ### 
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
    
    # get the dy shape, data and rares
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.SF, cuts.Zveto]), '', xlabel)
    dy_shape = makeDYMllShape('mll',specialcut,scutstring )
    others = treeOTHERS.getTH1F(lint, var+"others"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.SF, cuts.Zveto, cuts.lepsFromZ]), '', xlabel)
    ttz_SF = treeTTZ.getTH1F(lint, var+"ttz_SF"+scutstring, treevar, nbins, 1,1, cuts.AddList([specialcut, cuts.SignalRegion, cuts.goodLepton, cuts.SF, cuts.Zveto,cuts.lepsFromZ]), '', xlabel)
    zz_SF = treeZZ.getTH1F(lint, var+"zz_SF"+scutstring, treevar,    nbins, 1,1, cuts.AddList([specialcut, cuts.SignalRegion, cuts.goodLepton, cuts.SF, cuts.Zveto,cuts.lepsFromZ]), '', xlabel)
    wz_SF = treeWZ.getTH1F(lint, var+"wz_SF"+scutstring, treevar,    nbins, 1,1, cuts.AddList([specialcut, cuts.SignalRegion, cuts.goodLepton, cuts.SF, cuts.Zveto,cuts.lepsFromZ]), '', xlabel)
    ttz_SF.Scale(1.36);wz_SF.Scale(1.0);zz_SF.Scale(1.57) 
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
    # aesthetics
    
    rare = copy.deepcopy(others)
    rare.Add(ttz_SF, 1.)
    rare.Add(zz_SF, 1.)
    rare.Add(wz_SF, 1.)
    rare.SetFillColorAlpha(r.kCyan+2, 0.7);rare.SetTitle("rares");rare.SetLineColor(r.kBlack)
    dy_shape.SetFillColorAlpha(r.kYellow-9, 0.7);dy_shape.SetTitle("E_{T}^{miss} templates");
    prediction.SetFillColorAlpha(r.kRed-9, 0.7);prediction.SetTitle("FS");
    da_SF.SetTitle("data SF")       
    da_perGeV = copy.deepcopy(da_SF); dy_perGeV = copy.deepcopy(dy_shape);fs_perGeV = copy.deepcopy(prediction);rare_perGeV = copy.deepcopy(rare)
    for ib in range(1,da_SF.GetNbinsX()+1):
        da_perGeV.SetBinContent(ib,  da_SF.GetBinContent(ib)/(da_SF.GetBinLowEdge(ib+1)-da_SF.GetBinLowEdge(ib)))
        dy_perGeV.SetBinContent(ib,  dy_shape.GetBinContent(ib)/(dy_shape.GetBinLowEdge(ib+1)-dy_shape.GetBinLowEdge(ib)))
        fs_perGeV.SetBinContent(ib,  prediction.GetBinContent(ib)/(prediction.GetBinLowEdge(ib+1)-prediction.GetBinLowEdge(ib)))
        rare_perGeV.SetBinContent(ib, rare.GetBinContent(ib)/(rare.GetBinLowEdge(ib+1)-rare.GetBinLowEdge(ib)))                            
        da_perGeV.SetBinError(ib,  da_SF.GetBinError(ib)/(da_SF.GetBinLowEdge(ib+1)-da_SF.GetBinLowEdge(ib)))
        dy_perGeV.SetBinError(ib,  dy_shape.GetBinError(ib)/(dy_shape.GetBinLowEdge(ib+1)-dy_shape.GetBinLowEdge(ib)))
        fs_perGeV.SetBinError(ib,  prediction.GetBinError(ib)/(prediction.GetBinLowEdge(ib+1)-prediction.GetBinLowEdge(ib)))
        rare_perGeV.SetBinError(ib,  rare.GetBinError(ib)/(rare.GetBinLowEdge(ib+1)-rare.GetBinLowEdge(ib)))


    mc_stack.Add(dy_shape  );                      mc_stack_perGeV.Add(dy_perGeV); 
    mc_stack.Add(rare      );                      mc_stack_perGeV.Add(rare_perGeV);
    mc_stack.Add(prediction);                      mc_stack_perGeV.Add(fs_perGeV); 
    mc_stack.Draw();                               mc_stack_perGeV.Draw("SAME");
    mc_stack.GetXaxis().SetTitle('m_{ll} [GeV]');  mc_stack_perGeV.GetXaxis().SetTitle('m_{ll} [GeV]');
    mc_full = copy.deepcopy(dy_shape);             mc_full_perGeV =  copy.deepcopy(dy_perGeV);
    mc_full.Add(prediction, 1.);                   mc_full_perGeV.Add(fs_perGeV, 1.); 
    mc_full.Add(rare, 1.);                         mc_full_perGeV.Add(rare_perGeV, 1.); 
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
    da_SF.GetYaxis().SetRangeUser(0.1, 1.30*maxCont)
    da_perGeV.GetYaxis().SetRangeUser(0.1, 1.30*maxCont_perGeV)
    mc_stack.SetMaximum(1.3*maxCont)
    mc_stack_perGeV.SetMaximum(1.3*maxCont_perGeV)
    SetOwnership(mc_stack, 0);SetOwnership(da_SF, 0);SetOwnership(mc_stack_perGeV, 0);SetOwnership(da_perGeV, 0);SetOwnership(mc_full, 0); SetOwnership(mc_full_perGeV, 0);              

    print helper.bcolors.HEADER + '[result scaled by RSFOF for DATA] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_result = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s'%(newLumiString, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.67, 0.59, 0.90, 0.85)
    plot_result.addStack(mc_stack, "HIST" , 1, 1)
    plot_result.addHisto(mc_full_e, 'E2,SAME'  , '' , 'PL', r.kBlack , 1, -1)
    plot_result.addHisto(da_SF    , 'E1,SAME'  , 'observed data', 'PL', r.kBlack  , 1,  0)
    plot_result.saveRatio(1, 1, 0, lint, da_SF, mc_full, 0. , int(maxrat+1.0)) 
    makeResultsTable(da_SF, prediction, dy_shape, rare, mc_full, scutstring)
    plot_result = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s_perGeV'%(newLumiString, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.67, 0.59, 0.90, 0.85)
    plot_result.addStack(mc_stack_perGeV, "HIST" , 1, 1)
    plot_result.addHisto(mc_full_perGeV_e, 'E2,SAME'  , ''  , 'PL', r.kBlack , 1, -1)
    plot_result.addHisto(da_perGeV       , 'E1,SAME'  , 'observed data', 'PL', r.kBlack  , 1,  0)
    plot_result.saveRatio(1, 1, 0, lint, da_perGeV, mc_full_perGeV, 0. , int(maxrat+1.0)) 
    
    tfile = r.TFile('datacards/forDatacards_Edge_Moriond2017_%s.root'%scutstring,'recreate')
    tfile.cd()
    tfile.WriteTObject(dy_shape  , 'dy_shape')
    tfile.WriteTObject(rare      , 'mc_full')
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
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
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

    print 'Going to load DATA and MC trees...'
    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    zzDatasets = ['ZZTo2L2Nu']
    wzDatasets = ['WZTo3LNu']
    ttzDatasets = ['TTZToLLNuNu']
    othersDatasets = ['WWZ', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll']
    fsDatasets = ['TTTT', 'TTHnobb_pow', 'VHToNonbb',  'TTJets_DiLepton', 'TTLLJets_m1to10', 'TBar_tch_powheg', 'T_tch_powheg', 'WWTo2L2Nu', 'WWW', 'TTWToLNu',  'TTWToQQ', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT',   'WJetsToLNu_LO']
    
    mcDatasets = fsDatasets+dyDatasets + othersDatasets + zzDatasets + wzDatasets + ttzDatasets
 
 
                                                                                     
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
    treeFS = Sample.Tree(helper.selectSamples(opts.sampleFile, fsDatasets, 'FS'), 'FS'  , 0)
    treeOTHERS = Sample.Tree(helper.selectSamples(opts.sampleFile, othersDatasets, 'OTHERS'), 'OTHERS'  , 0)
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
    lint = 36.8  ; maxrun = 999999; lint_str = '36.8invfb'
    print 'Running with an integrated luminosity of %.2f fb-1' %(lint)

    ## ============================================================
    ## ========== set RSFOF globally ==============================
    ## ============================================================
    rsfof, rsfof_e, rsfof_mc, rsfof_mc_e, rmue_a_da, rmue_a_mc, rmue_b_da, rmue_b_mc =  makeFactorsTable()
    ## result plots in different variables:
    for v in ['mll']:#'nll_noMET', 'nll_noMLB', 'nll_noZPT', 'nll_noLDP']: #'iso1', 'iso2', 'mll', 'nll', 'nb', 'nj', 'zpt', 'mlb', 'met', 'ldp', 'pt1', 'pt2']:
        makeResultData('Edge_Moriond2017', v,maxrun,lint,specialcut='' , scutstring = '',    _options='returnplot,splitFlavor')
        resultPlotLoNll = makeResultData('Edge_Moriond2017', v,maxrun,lint,specialcut='nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) < 21.' , scutstring = 'nllBelow21',    _options='returnplot,splitFlavor')
        resultPlotHiNll = makeResultData('Edge_Moriond2017', v,maxrun,lint,specialcut='nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) >= 21.' , scutstring = 'nllAbove21',    _options='returnplot,splitFlavor')

   # makeRSOFTable('Edge_Moriond2017')
   # makeClosureTestPlots('Edge_Moriond2017','nll','', 'inclusive', True)
   # makeClosureTestPlots('Edge_Moriond2017','nll','lepsMll_Edge <= 60. && lepsMll_Edge > 20', 'mll20-60', True)
   # makeClosureTestPlots('Edge_Moriond2017','nll','lepsMll_Edge <= 86. && lepsMll_Edge > 60', 'mll60-86', True)
   # makeClosureTestPlots('Edge_Moriond2017','nll','lepsMll_Edge <= 150. && lepsMll_Edge > 96', 'mll96-150', True)
   # makeClosureTestPlots('Edge_Moriond2017','nll','lepsMll_Edge <= 200. && lepsMll_Edge > 150', 'mll150-200', True)
   # makeClosureTestPlots('Edge_Moriond2017','nll','lepsMll_Edge <= 300. && lepsMll_Edge > 200', 'mll200-300', True)
   # makeClosureTestPlots('Edge_Moriond2017','nll','lepsMll_Edge <= 400. && lepsMll_Edge > 300', 'mll300-400', True)
   # makeClosureTestPlots('Edge_Moriond2017','nll','lepsMll_Edge > 400', 'mll400', True)
   # makeClosureTestPlots('Edge_Moriond2017','mll','', 'inclusive', True)
   # makeClosureTestPlots('Edge_Moriond2017','mll','nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21.' , 'nllAbove21', True)
   # makeClosureTestPlots('Edge_Moriond2017','mll','nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) <= 21.' , 'nllBelow21', True)
    
