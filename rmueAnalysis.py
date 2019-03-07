#####################################################################
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
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors
import math, sys, optparse, array, copy, subprocess

import include.LeptonSF
import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Rounder    as Rounder

from itertools import product

class bcolors: 
    HEADER = '\033[95m' 
    OKBLUE = '\033[94m' 
    OKGREEN = '\033[92m' 
    WARNING = '\033[93m' 
    FAIL = '\033[91m' 
    ENDC = '\033[0m' 
    BOLD = '\033[1m' 
    UNDERLINE = '\033[4m' 


def makeTable( mc_a, mc_a_e, data_a, data_a_e,  mc_b, mc_b_e, data_b, data_b_e):
    line0 = '  \hline'
    line1 = '  Data &' 
    line2 = '  MC   &' 
    line0 += ' & A & B & '
    line1 += ' %.4f $\\pm$ %.4f& %.4f $\\pm$ %.4f   %s' %(data_a, data_a_e, data_b, data_b_e , '\\\\')
    line2 += ' %.4f $\\pm$ %.4f& %.4f $\\pm$ %.4f   %s' %(mc_a, mc_a_e, mc_b, mc_b_e ,  '\\\\')
    line0 += '\\hline'; line2 += '\\hline';

    helper.ensureDirectory('plots/rmue/%s/'%lumi_str)
    helper.ensureDirectory('plots/rmue/%s/tables/'%lumi_str)
    compTableFile = open('plots/rmue/%s/tables/resultTable_%s%s.txt'%(lumi_str, str(lumi).replace('.','p'), "rmue"),'w')
    compTableFile.write(line0+'\n')
    compTableFile.write(line1+'\n')
    compTableFile.write(line2+'\n')                                                                                             
    print line0
    print line1
    print line2                                                                                                                                                                      


def saveInFile(theFile, measuredValueMC, measuredValueUncMC, measuredValueUncSystMC, measuredValueData, measuredValueUncData, measuredValueUncSystData, tag):

    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rmue") != -1 and line.find(tag) != -1 and  line.find("alone") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rmue        alone           DATA        %.4f      %0.4f       %.4f\n'%(measuredValueData, measuredValueUncData, measuredValueUncSystData))
            else:
                foutput.write('rmue        alone           MC          %.4f      %0.4f       %.4f\n'%(measuredValueMC, measuredValueUncMC, measuredValueUncSystMC))           
        elif line.find("rmue") != -1 and line.find(tag) != -1 and  line.find("factor") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rmue        factor           DATA        %.4f      %0.4f       \n'%(getFactor(measuredValueData, measuredValueUncData, measuredValueUncSystData)))
            else:
                foutput.write('rmue        factor           MC          %.4f      %0.4f       \n'%(getFactor(measuredValueMC, measuredValueUncMC, measuredValueUncSystMC)))           
        
        elif line.find("rmue") != -1 and line.find(tag) != -1 and  line.find("coeffA") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rmue        coeffA           DATA        %.4f      %0.4f       \n'%(measuredValueData, measuredValueUncSystData))
            else:
                foutput.write('rmue        coeffA           MC        %.4f      %0.4f       \n'%(measuredValueMC, measuredValueUncSystMC))

        elif line.find("rmue") != -1 and line.find(tag) != -1 and  line.find("coeffB") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rmue        coeffB           DATA        %.4f      %0.4f       \n'%(measuredValueData, measuredValueUncSystData))
            else:
                foutput.write('rmue        coeffB           MC        %.4f      %0.4f       \n'%(measuredValueMC, measuredValueUncSystMC))
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)
                                                                                                                                                                                               


def calc_rmue(Nmm, Nee, Emm, Eee):

    val = [0, 0]
    print Nmm, Nee
    if(Nmm >= 0 and Nee > 0):
        val[0] = math.sqrt(Nmm/Nee)
        val[1] = math.sqrt(0.25 * Emm * Emm / (Nmm * Nee) + 0.25 * Eee * Eee * Nmm / (Nee * Nee * Nee))
    else:
        print bcolors.HEADER + "[rmueAnalysis] " + bcolors.FAIL + " The yields in Nee and Nmm are <=0 " + bcolors.ENDC
    return val


def make_rmue(histo_mm, histo_ee):

    ratio = histo_mm.Clone("rmue_" + histo_mm.GetName())
    ratio.GetYaxis().SetTitle("r_{#mu e}")
    ratio.GetXaxis().SetTitle(histo_mm.GetXaxis().GetTitle())
  
    for i in range(0, histo_mm.GetNbinsX()+1):
        Nmm = histo_mm.GetBinContent(i)
        Nee = histo_ee.GetBinContent(i)
        Emm = histo_mm.GetBinError(i)
        Eee = histo_ee.GetBinError(i)
        if(Nee != 0):
            if(Nee == 0 or Nmm == 0):
              continue
            val = calc_rmue(Nmm, Nee, Emm, Eee)
            ratio.SetBinContent(i, val[0])
            ratio.SetBinError(i, val[1])
    return ratio



def convertToFactor(histo, getGraph=False):
    tmp_histo = copy.deepcopy(histo)
    for i in range(tmp_histo.GetNbinsX()+1):
        rmue     = tmp_histo.GetBinContent(i)
        rmue_err = tmp_histo.GetBinError  (i)
        if rmue:
            fac = 0.5*(rmue + 1./rmue)
            err = 0.5*((1. - 1./(rmue**2))*rmue_err)
        else:
            fac = 0.
            err = 0.
        tmp_histo.SetBinContent(i, fac)
        tmp_histo.SetBinError  (i, err)
        tmp_histo.GetYaxis().SetTitle("0.5(r_{#mu e} + 1/r_{#mu e})")
    if getGraph:
        return TGraphErrors(tmp_histo)
    else:
        return tmp_histo                                                         


def getFactor(rmue, rmue_err, rmue_err_syst):
    fac = 0.5*(rmue + 1./rmue)
    err = 0.5*((1. - 1./(rmue**2))*math.sqrt(rmue_err**2 + rmue_err_syst**2) )
    return fac, err                                                     


def makeAnalysis(treeDA, treeMC, cuts, specialcut, tag, save, ingredientsFile):
    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Producing histograms...' + bcolors.ENDC
    MCSDYControlMllee =   treeMC.getStack(lumi,"MCDYControlMllee", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut,cuts.DYControlRegionNoMll,cuts.goodLepton17,cuts.ee]), '', labelx, "1", kf)
    MCSDYControlMllmm =   treeMC.getStack(lumi,"MCDYControlMllmm", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut,cuts.DYControlRegionNoMll,cuts.goodLepton17,cuts.mm]), '', labelx, "1", kf)
    DADYControlMllee =  treeDA.getTH1F(lumi,"DADYControlMllee", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut,cuts.DYControlRegionNoMll,cuts.goodLepton17,cuts.ee]),'', labelx, "1", kf)
    DADYControlMllmm =  treeDA.getTH1F(lumi,"DADYControlMllmm", "lepsMll_Edge", 20, 20, 200, cuts.AddList([specialcut,cuts.DYControlRegionNoMll,cuts.goodLepton17,cuts.mm]),'', labelx, "1", kf)

    MCSDYControlMlleevalue = treeMC.getStack(lumi,"MCDYControlMlleevalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut,cuts.DYControlRegion,cuts.goodLepton17,cuts.ee]), '', labelx, "1", kf)
    MCSDYControlMllmmvalue = treeMC.getStack(lumi,"MCDYControlMllmmvalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut,cuts.DYControlRegion,cuts.goodLepton17,cuts.mm]), '', labelx, "1", kf)
    DADYControlMlleevalue =treeDA.getTH1F(lumi,"DADYCMlleevalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.DYControlRegion,cuts.goodLepton17,cuts.ee]),'', labelx, "1", kf)
    DADYControlMllmmvalue =treeDA.getTH1F(lumi,"DADYCMllmmvalue", "lepsMll_Edge", 1, 20, 300, cuts.AddList([specialcut, cuts.DYControlRegion,cuts.goodLepton17,cuts.mm]),'', labelx, "1", kf)
    ##
    MCSDYControlMETee =   treeMC.getStack(lumi,"MCDYControlMETee", "MET_pt_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.DYControlRegionNoMET,cuts.goodLepton17,cuts.ee]), '', labelmet, "1", kf)
    MCSDYControlMETmm =   treeMC.getStack(lumi,"MCDYControlMETmm", "MET_pt_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.DYControlRegionNoMET,cuts.goodLepton17,cuts.mm]), '', labelmet, "1", kf)
    DADYControlMETee =  treeDA.getTH1F(lumi,"DADYControlMETee", "MET_pt_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.DYControlRegionNoMET,cuts.goodLepton17,cuts.ee]), '', labelmet, "1", kf)
    DADYControlMETmm =  treeDA.getTH1F(lumi,"DADYControlMETmm", "MET_pt_Edge", 5, 20, 200, cuts.AddList([specialcut, cuts.DYControlRegionNoMET,cuts.goodLepton17,cuts.mm]), '', labelmet, "1", kf)
    ##
    MCSDYControlmt2ee =   treeMC.getStack(lumi,"MCDYControlmt2ee", "mt2_Edge", 8, 0, 160, cuts.AddList([specialcut, cuts.DYControlRegion,cuts.goodLepton17,cuts.ee]), '', labelmt2, "1", kf)
    MCSDYControlmt2mm =   treeMC.getStack(lumi,"MCDYControlmt2mm", "mt2_Edge", 8, 0, 160, cuts.AddList([specialcut, cuts.DYControlRegion,cuts.goodLepton17,cuts.mm]), '', labelmt2, "1", kf)
    DADYControlmt2ee =  treeDA.getTH1F(lumi,"DADYControlmt2ee", "mt2_Edge", 8, 0, 160, cuts.AddList([specialcut, cuts.DYControlRegion,cuts.goodLepton17,cuts.ee]), '', labelmt2, "1", kf)
    DADYControlmt2mm =  treeDA.getTH1F(lumi,"DADYControlmt2mm", "mt2_Edge", 8, 0, 160, cuts.AddList([specialcut, cuts.DYControlRegion,cuts.goodLepton17,cuts.mm]), '', labelmt2, "1", kf)
    ##
    MCSDYControlPT2ee =   treeMC.getStack(lumi, "MCDYControlPT2ee", "Lep2_pt_Edge", 8, 20, 200, cuts.AddList([specialcut,  cuts.DYControlRegion,cuts.goodLepton17,cuts.ee]), '', labelpt2, "1", kf)
    MCSDYControlPT2mm =   treeMC.getStack(lumi, "MCDYControlPT2mm", "Lep2_pt_Edge", 8, 20, 200, cuts.AddList([specialcut,  cuts.DYControlRegion,cuts.goodLepton17,cuts.mm]), '', labelpt2, "1", kf)
    DADYControlPT2ee =  treeDA.getTH1F(lumi, "DADYControlPT2ee", "Lep2_pt_Edge", 8, 20, 200, cuts.AddList([specialcut, cuts.DYControlRegion,cuts.goodLepton17,cuts.ee]), '', labelpt2, "1", kf)
    DADYControlPT2mm =  treeDA.getTH1F(lumi, "DADYControlPT2mm", "Lep2_pt_Edge", 8, 20, 200, cuts.AddList([specialcut, cuts.DYControlRegion,cuts.goodLepton17,cuts.mm]), '', labelpt2, "1", kf)

    DYControlMllee = 0;DYControlMllmm = 0; DYControlMlleevalue = 0; DYControlMllmmvalue = 0;DYControlMETee = 0; DYControlMETmm = 0; DYControlmt2ee = 0; DYControlmt2mm = 0;DYControlmt2ee = 0;DYControlmt2mm =0 ; DYControlPT2ee = 0; DYControlPT2mm = 0;
    for _i,_h in enumerate(MCSDYControlMllee.GetHists()):
        if not DYControlMllee: MCDYControlMllee = copy.deepcopy(_h)
        else:MCDYControlMllee.Add(_h, 1.)               
        DYControlMllee = 1                                             
    for _i,_h in enumerate(MCSDYControlMllmm.GetHists()):
        if not DYControlMllmm: MCDYControlMllmm = copy.deepcopy(_h)
        else:MCDYControlMllmm.Add(_h, 1.)               
        DYControlMllmm = 1                                             
    for _i,_h in enumerate(MCSDYControlMlleevalue.GetHists()):
        if not DYControlMlleevalue: MCDYControlMlleevalue = copy.deepcopy(_h)
        else:MCDYControlMlleevalue.Add(_h, 1.)               
        DYControlMlleevalue = 1                                                 
    for _i,_h in enumerate(MCSDYControlMllmmvalue.GetHists()):
        if not DYControlMllmmvalue: MCDYControlMllmmvalue = copy.deepcopy(_h)
        else:MCDYControlMllmmvalue.Add(_h, 1.)               
        DYControlMllmmvalue = 1                                                
    for _i,_h in enumerate(MCSDYControlMETee.GetHists()):
        if not DYControlMETee: MCDYControlMETee = copy.deepcopy(_h)
        else:MCDYControlMETee.Add(_h, 1.)               
        DYControlMETee = 1                                                                           
    for _i,_h in enumerate(MCSDYControlMETmm.GetHists()):
        if not DYControlMETmm: MCDYControlMETmm = copy.deepcopy(_h)
        else:MCDYControlMETmm.Add(_h, 1.)               
        DYControlMETmm = 1                                               
    for _i,_h in enumerate(MCSDYControlmt2ee.GetHists()):
        if not DYControlmt2ee: MCDYControlmt2ee = copy.deepcopy(_h)
        else:MCDYControlmt2ee.Add(_h, 1.)               
        DYControlmt2ee = 1                                               
    for _i,_h in enumerate(MCSDYControlmt2mm.GetHists()):
        if not DYControlmt2mm: MCDYControlmt2mm = copy.deepcopy(_h)
        else:MCDYControlmt2mm.Add(_h, 1.)               
        DYControlmt2mm = 1                                               
    for _i,_h in enumerate(MCSDYControlPT2ee.GetHists()):
        if not DYControlPT2ee: MCDYControlPT2ee = copy.deepcopy(_h)
        else:MCDYControlPT2ee.Add(_h, 1.)               
        DYControlPT2ee = 1                                           
    for _i,_h in enumerate(MCSDYControlPT2mm.GetHists()):
        if not DYControlPT2mm: MCDYControlPT2mm = copy.deepcopy(_h)
        else:MCDYControlPT2mm.Add(_h, 1.)               
        DYControlPT2mm = 1                                           
    print "got all TH1Fs" 
    MCDYControlMll =           make_rmue(MCDYControlMllmm, MCDYControlMllee)
    DADYControlMll =         make_rmue(DADYControlMllmm, DADYControlMllee)
    MCDYControlMllvalue =      make_rmue(MCDYControlMllmmvalue, MCDYControlMlleevalue)
    DADYControlMllvalue =    make_rmue(DADYControlMllmmvalue, DADYControlMlleevalue)
    MCDYControlMET =           make_rmue(MCDYControlMETmm, MCDYControlMETee)
    DADYControlMET =         make_rmue(DADYControlMETmm, DADYControlMETee)
    #MCDYControlJet =           make_rmue(MCDYControlJetmm, MCDYControlJetee)
    #DADYControlJet =         make_rmue(DADYControlJetmm, DADYControlJetee)
    MCDYControlmt2 =           make_rmue(MCDYControlmt2mm, MCDYControlmt2ee)
    DADYControlmt2 =         make_rmue(DADYControlmt2mm, DADYControlmt2ee)
    MCDYControlPT2 =           make_rmue(MCDYControlPT2mm, MCDYControlPT2ee)
    DADYControlPT2   =       make_rmue(DADYControlPT2mm, DADYControlPT2ee)


    print 'fits to pt2'
    fun_da = r.TF1('fun_da','[0]+[1]/x', 20, 400)
    fun_mc = r.TF1('fun_mc','[0]+[1]/x', 20, 400)
    fun_da.SetLineWidth(2); fun_mc.SetLineWidth(2)
    DADYControlPT2.Fit('fun_da','rQ'); 
    data_a = fun_da.GetParameter(0); data_a_e = fun_da.GetParError(0)
    data_b = fun_da.GetParameter(1); data_b_e = fun_da.GetParError(1)
    MCDYControlPT2.Fit('fun_mc','rQ')
    mc_a = fun_mc.GetParameter(0); mc_a_e = fun_mc.GetParError(0)
    mc_b = fun_mc.GetParameter(1); mc_b_e = fun_mc.GetParError(1)
    da_string = 'Fit on data: (%4.2f #pm %4.2f) + (%4.2f #pm %4.2f) p_{T}^{-1}'%(data_a,data_a_e, data_b, data_b_e)
    mc_string = 'Fit on MC:   (%4.2f #pm %4.2f) + (%4.2f #pm %4.2f) p_{T}^{-1}'%(mc_a,mc_a_e, mc_b, mc_b_e)

    plot_rmue_pt2 = Canvas.Canvas('rmue/%s/plot_rmue_pT2_%s'%(lumi_str, tag), 'png,pdf', 0.6, 0.2, 0.75, 0.35)
    plot_rmue_pt2.addHisto(MCDYControlPT2, 'PE', 'MC', 'PL', r.kRed+1, 1, 0)
    plot_rmue_pt2.addHisto(DADYControlPT2, 'PE,SAME', 'DATA', 'PL', r.kBlue , 1, 1)
    plot_rmue_pt2.addHisto(fun_da, 'L,SAME', 'Fit - Data', 'L', r.kBlack, 1,2)
    plot_rmue_pt2.addHisto(fun_mc, 'L,SAME', 'Fit - MC', 'L', r.kRed+1, 1,3)
    plot_rmue_pt2.addLatex(0.2,0.2, da_string, 42, 0.02)
    plot_rmue_pt2.addLatex(0.2,0.25, mc_string, 42, 0.02)
    #plot_rmue_pt2.save(1, 1, 0, lumi)
    plot_rmue_pt2.save(1, 1, 0, lumi, "",  0.2, 1.8)




    print "made all the rmues" 
    factorMCDYControlMll =          convertToFactor(MCDYControlMll)
    factorDADYControlMll =        convertToFactor(DADYControlMll)
    factorMCDYControlMllvalue =     convertToFactor(MCDYControlMllvalue)
    factorDADYControlMllvalue =   convertToFactor(DADYControlMllvalue)
    factorMCDYControlMET =          convertToFactor(MCDYControlMET)
    factorDADYControlMET =        convertToFactor(DADYControlMET)
    #factorMCDYControlJet =          convertToFactor(MCDYControlJet)
    #factorDADYControlJet =        convertToFactor(DADYControlJet)
    factorMCDYControlmt2 =          convertToFactor(MCDYControlmt2)
    factorDADYControlmt2 =        convertToFactor(DADYControlmt2)

    print "convert all factors" 
    systematicForrmue = 0.1 
    MCrmuemeasured =               MCDYControlMllvalue.GetBinContent(1)
    MCrmuemeasuredUnc =            MCDYControlMllvalue.GetBinError(1)
    MCrmuemeasuredUncSyst =        MCrmuemeasured * systematicForrmue
    DArmuemeasured =             DADYControlMllvalue.GetBinContent(1)
    DArmuemeasuredUnc =          DADYControlMllvalue.GetBinError(1)
    DArmuemeasuredUncSyst =      DArmuemeasured * systematicForrmue
    DArmuemeasuredUncTot  =      math.sqrt(DArmuemeasuredUnc**2 + DArmuemeasuredUncSyst**2)
    
    factorMCrmuemeasured =             factorMCDYControlMllvalue.GetBinContent(1)
    factorMCrmuemeasuredUnc =          factorMCDYControlMllvalue.GetBinError(1)
    factorMCrmuemeasuredUncSyst =      (1.0-1.0/MCrmuemeasured**2) * MCrmuemeasuredUncSyst if MCrmuemeasured != 0 else 1000.00
    factorMCrmuemeasuredUncTot  =      math.sqrt(MCrmuemeasuredUnc**2 + MCrmuemeasuredUncSyst**2)
    factorDArmuemeasured =             factorDADYControlMllvalue.GetBinContent(1)
    factorDArmuemeasuredUnc =          factorDADYControlMllvalue.GetBinError(1)
    factorDArmuemeasuredUncSyst =      (1.0-1.0/DArmuemeasured**2) * DArmuemeasuredUncSyst if DArmuemeasured != 0 else 1000.00
    factorDArmuemeasuredUncTot  =      math.sqrt(DArmuemeasuredUnc**2 + DArmuemeasuredUncSyst**2)                                    
    
    
    print "ready to save"
    if save==True:
        #saveInFile(ingredientsFile, MCrmuemeasured, MCrmuemeasuredUnc, MCrmuemeasuredUncSyst, DA16rmuemeasured, DA16rmuemeasuredUnc, DA16rmuemeasuredUncSyst, 'alone')
        #saveInFile(ingredientsFile, factorMCrmuemeasured, factorMCrmuemeasuredUnc, factorMCrmuemeasuredUncSyst, factorDA16rmuemeasured, factorDA16rmuemeasuredUnc, factorDA16rmuemeasuredUncSyst, 'factor')
        saveInFile(ingredientsFile, mc_a, 0., mc_a_e, data_a, 0., data_a_e,  'coeffA')
        saveInFile(ingredientsFile, mc_b, 0., mc_b_e, data_b, 0., data_b_e,  'coeffB')
        makeTable( mc_a, mc_a_e, data_a, data_a_e,  mc_b,mc_b_e, data_b, data_b_e)
    MCDYControlPT2.Fit('fun_mc','rQ')
    mc_a = fun_mc.GetParameter(0); mc_a_e = fun_mc.GetParError(0)
    mc_b = fun_mc.GetParameter(1); mc_b_e = fun_mc.GetParError(1)
      
    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Producing plots...' + bcolors.ENDC
    plot_rmue_mll = Canvas.Canvas('rmue/%s/plot_rmue_mll_%s'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_rmue_mll.addHisto(MCDYControlMll, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    plot_rmue_mll.addHisto(DADYControlMll, 'PE,SAME', 'DATA', 'PL', r.kBlue , 1, 0)
    #plot_rmue_mll.addBand(MCDYControlMll.GetXaxis().GetXmin(), DA16rmuemeasured-DA16rmuemeasuredUncTot, MCDYControlMll.GetXaxis().GetXmax(), DA16rmuemeasured+DA16rmuemeasuredUncTot, r.kGreen, 0.2)
    #plot_rmue_mll.addLine(MCDYControlMll.GetXaxis().GetXmin(), DA16rmuemeasured, MCDYControlMll.GetXaxis().GetXmax(), DA16rmuemeasured,r.kGreen)
    plot_rmue_mll.save(1, 1, 0, lumi, "",0.2, 1.8)
    
    plot_rmue_met = Canvas.Canvas('rmue/%s/plot_rmue_met_%s'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_rmue_met.addHisto(MCDYControlMET, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    #plot_rmue_met.addHisto(DA16DYControlMET, 'PE,SAME', 'DATA 2016', 'PL', r.kBlack , 1, 0)
    plot_rmue_met.addHisto(DADYControlMET, 'PE,SAME', 'DATA', 'PL', r.kBlue , 1, 0)
    #plot_rmue_met.addBand(MCDYControlMET.GetXaxis().GetXmin(), DA16rmuemeasured-DA16rmuemeasuredUncTot, MCDYControlMET.GetXaxis().GetXmax(), DA16rmuemeasured+DA16rmuemeasuredUncTot, r.kGreen, 0.2)
    #plot_rmue_met.addLine(MCDYControlMET.GetXaxis().GetXmin(), DA16rmuemeasured, MCDYControlMET.GetXaxis().GetXmax(), DA16rmuemeasured,r.kGreen)
    plot_rmue_met.save(1, 1, 0, lumi, "",0.2, 1.8)
  
#    plot_rmue_jet = Canvas.Canvas('rmue/%s_%s/plot_rmue_jet'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
#    plot_rmue_jet.addHisto(MCDYControlJet, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
#    plot_rmue_jet.addHisto(DA16DYControlJet, 'PE,SAME', 'DATA 2016', 'PL', r.kBlack , 1, 0)
#    plot_rmue_jet.addHisto(DADYControlJet, 'PE,SAME', 'DATA 2017', 'PL', r.kBlue , 1, 0)
#    plot_rmue_jet.addBand(MCDYControlJet.GetXaxis().GetXmin(), DA16rmuemeasured-DA16rmuemeasuredUncTot, MCDYControlJet.GetXaxis().GetXmax(), DA16rmuemeasured+DA16rmuemeasuredUncTot, r.kGreen, 0.2)
#    plot_rmue_jet.addLine(MCDYControlJet.GetXaxis().GetXmin(), DA16rmuemeasured, MCDYControlJet.GetXaxis().GetXmax(), DA16rmuemeasured,r.kGreen)
#    plot_rmue_jet.save(1, 1, 0, lumi, 0.2, 1.8)
    
    plot_factor_mll = Canvas.Canvas('rmue/%s/plot_factor_mll_%s'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_factor_mll.addHisto(factorMCDYControlMll, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    #plot_factor_mll.addHisto(factorDA16DYControlMll, 'PE,SAME', 'DATA 2016', 'PL', r.kBlack , 1, 0)
    plot_factor_mll.addHisto(factorDADYControlMll, 'PE,SAME', 'DATA', 'PL', r.kBlue , 1, 0)
    #plot_factor_mll.addBand(factorMCDYControlMll.GetXaxis().GetXmin(), factorDA16rmuemeasured-factorDA16rmuemeasuredUncTot, factorMCDYControlMll.GetXaxis().GetXmax(), factorDA16rmuemeasured+factorDA16rmuemeasuredUncTot, r.kGreen, 0.2)
    #lot_factor_mll.addLine(factorMCDYControlMll.GetXaxis().GetXmin(), factorDA16rmuemeasured, factorMCDYControlMll.GetXaxis().GetXmax(), factorDA16rmuemeasured,r.kGreen)
    plot_factor_mll.save(1, 1, 0, lumi, "",0.2, 1.8)

    plot_factor_met = Canvas.Canvas('rmue/%s/plot_factor_met_%s'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_factor_met.addHisto(factorMCDYControlMET, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    #plot_factor_met.addHisto(factorDA16DYControlMET, 'PE,SAME', 'DATA 2016', 'PL', r.kBlack , 1, 0)
    plot_factor_met.addHisto(factorDADYControlMET, 'PE,SAME', 'DATA', 'PL', r.kBlue , 1, 0)
    #plot_factor_met.addBand(factorMCDYControlMET.GetXaxis().GetXmin(), factorDA16rmuemeasured-factorDA16rmuemeasuredUncTot, factorMCDYControlMET.GetXaxis().GetXmax(), factorDA16rmuemeasured+factorDA16rmuemeasuredUncTot, r.kGreen, 0.2)
    plot_factor_met.addLine(factorMCDYControlMET.GetXaxis().GetXmin(), factorDA16rmuemeasured, factorMCDYControlMET.GetXaxis().GetXmax(), factorDA16rmuemeasured,r.kGreen)
    plot_factor_met.save(1, 1, 0, lumi, "",0.2, 1.8)

#    plot_factor_jet = Canvas.Canvas('rmue/%s_%s/plot_factor_jet'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
#    plot_factor_jet.addHisto(factorMCDYControlJet, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
#    plot_factor_jet.addHisto(factorDA16DYControlJet, 'PE,SAME', 'DATA 2016', 'PL', r.kBlack , 1, 0)
#    plot_factor_jet.addHisto(factorDADYControlJet, 'PE,SAME', 'DATA 2017', 'PL', r.kBlue , 1, 0)
#    plot_factor_jet.addBand(factorMCDYControlJet.GetXaxis().GetXmin(), factorDA16rmuemeasured-factorDA16rmuemeasuredUncTot, factorMCDYControlJet.GetXaxis().GetXmax(), factorDA16rmuemeasured+factorDA16rmuemeasuredUncTot, r.kGreen, 0.2)
#    plot_factor_jet.addLine(factorMCDYControlJet.GetXaxis().GetXmin(), factorDA16rmuemeasured, factorMCDYControlJet.GetXaxis().GetXmax(), factorDA16rmuemeasured,r.kGreen)
#    plot_factor_jet.save(1, 1, 0, lumi, "Data/Prediction",0.2, 1.8)                                                                                                                                                                                                 
    plot_factor_mt2 = Canvas.Canvas('rmue/%s/plot_factor_mt2_%s'%(lumi_str, tag), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
    plot_factor_mt2.addHisto(factorMCDYControlmt2, 'PE', 'MC', 'PL', r.kRed+1 , 1, 0)
    #plot_factor_mt2.addHisto(factorDA16DYControlmt2, 'PE,SAME', 'DATA 2016', 'PL', r.kBlack , 1, 0)
    plot_factor_mt2.addHisto(factorDADYControlmt2, 'PE,SAME', 'DATA', 'PL', r.kBlue , 1, 0)
    #plot_factor_mt2.addBand(factorMCDYControlmt2.GetXaxis().GetXmin(), factorDA16rmuemeasured-factorDA16rmuemeasuredUncTot, factorMCDYControlmt2.GetXaxis().GetXmax(), factorDA16rmuemeasured+factorDA16rmuemeasuredUncTot, r.kGreen, 0.2)
    #plot_factor_mt2.addLine(factorMCDYControlmt2.GetXaxis().GetXmin(), factorDA16rmuemeasured, factorMCDYControlmt2.GetXaxis().GetXmax(), factorDA16rmuemeasured,r.kGreen)
    plot_factor_mt2.save(1, 1, 0, lumi, "",0.2, 1.8)                   

if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting r_mue analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    dyDatasets = ['DYJetsToLL_M50_ext_part1+DYJetsToLL_M50_ext_part2+DYJetsToLL_M50_ext_part3']
    ttDatasets = ['TTTo2L2Nu_part1+TTTo2L2Nu_part2','TTToSemiLeptonic'] # tt1l missing
    stDatasets = ['TW','TbarW'] # t and s (lol) channel missing 
    ttzDatasets = ['TTZ_LO_ext1','TTW_LO']#,'TTWZ','TTGJets_newpmx']
    zz2lDatasets = ['ZZTo2L2Q','ZZTo2L2Nu']
    zz4lDatasets = ['ZZTo4L_ext1_part1+ZZTo4L_ext1_part2',
                    #'GluGluToContinToZZTo2e2mu+GluGluToContinToZZTo2e2mu_ext1',
                    #'GluGluToContinToZZTo2e2nu+GluGluToContinToZZTo2e2nu_ext1',
                    #'GluGluToContinToZZTo2mu2nu+GluGluToContinToZZTo2mu2nu_ext1',
                    #'GluGluToContinToZZTo4e+GluGluToContinToZZTo4e_ext1',
                    #'GluGluToContinToZZTo4mu+GluGluToContinToZZTo4mu_ext1'
    ]
    wwDatasets = ['WWTo2L2Nu']
    wzDatasets = ['WZTo3LNu','WZTo2L2Q']
    raDatasets = [ 'TTH_amc']
    mcDatasets = zz4lDatasets + zz2lDatasets + ttzDatasets + raDatasets + wwDatasets +wzDatasets + stDatasets+  ttDatasets + dyDatasets

    daDatasets = [ '%s_Run2017%s'%(x,y) for x,y in product('DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(','), 'B,C,D,E,F'.split(','))]
    daDatasetsB = [ '%s_Run2017B'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]
    daDatasetsC = [ '%s_Run2017C'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]
    daDatasetsD = [ '%s_Run2017D'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]
    daDatasetsE = [ '%s_Run2017E'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]
    daDatasetsF = [ '%s_Run2017F'%(x) for x in 'DoubleEG,DoubleMuon,SingleElectron,SingleMuon,MuonEG'.split(',')]


    treeMC = Sample.Tree(helper.selectSamples("samples.dat", mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples("samples.dat", daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[rmueAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 41.9 ; maxrun = 276811; lumi_str = '41.9invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()
    cuts = CutManager.CutManager()

    kf = "noKFactor"
    labelx = "m_{ll} [GeV]"
    labelmet = "p_{T}^{miss} [GeV]"
    labelnjet = "N. Jets"
    labelmt2 = "M_{T2} [GeV]"
    labelpt2 = "p_{T}^{lep2}"
    #Cuts needed by rmue
    cuts = CutManager.CutManager()

    makeAnalysis(treeDA, treeMC,  cuts, '', '', True, opts.ingredientsFile)


    

